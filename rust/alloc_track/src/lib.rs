//! Dev-only allocation tracking global allocator.
//!
//! Opt-in via the `enabled` feature. When off, `TrackingAllocator` is a
//! zero-cost passthrough and `snap!` / `report_snapshot_diff` compile to
//! nothing. See `.plans/2026-04-18-alloc-track.md` for design notes.

use std::alloc::{
    GlobalAlloc,
    Layout,
    System,
};

pub struct TrackingAllocator<A: GlobalAlloc = System> {
    inner: A,
}

impl TrackingAllocator<System> {
    pub const fn new() -> Self {
        Self { inner: System }
    }
}

impl<A: GlobalAlloc> TrackingAllocator<A> {
    pub const fn with_inner(inner: A) -> Self {
        Self { inner }
    }
}

/// Snapshot of allocator activity since process start.
///
/// All fields are cumulative (monotonic) except `live_bytes`.
#[derive(Clone, Copy, Debug)]
pub struct Counters {
    /// Total bytes passed to `alloc` calls (never decreases).
    pub total_bytes: u64,
    /// Total `alloc` calls (never decreases). Deallocs NOT counted here.
    pub total_count: u64,
    /// Currently-live bytes (alloc bytes minus dealloc bytes).
    pub live_bytes: u64,
    /// High-water mark of `live_bytes` ever observed.
    ///
    /// Updated with Relaxed ordering after `live_bytes`, so a reader may
    /// briefly observe `peak_bytes < live_bytes` between the two stores.
    /// Self-corrects on the next alloc; acceptable for a dev-only tool.
    /// Note: if dealloc underflow drift has occurred (see
    /// `report_snapshot_diff`'s WARNING line), this value is a lower bound
    /// relative to reality and does not self-correct from dealloc-only workloads.
    pub peak_bytes: u64,
    /// Per-bucket alloc count histogram. See `BUCKET_LABELS` for bin edges.
    pub hist: [u64; 16],
}

impl Counters {
    pub const ZERO: Self = Self {
        total_bytes: 0,
        total_count: 0,
        live_bytes: 0,
        peak_bytes: 0,
        hist: [0; 16],
    };
}

impl Default for Counters {
    fn default() -> Self {
        Self::ZERO
    }
}

// --- enabled impl ---------------------------------------------------------

#[cfg(feature = "enabled")]
mod imp {
    use super::*;
    use std::sync::Mutex;
    use std::sync::atomic::{
        AtomicU64,
        Ordering,
    };

    // All counters use Relaxed ordering: tracker is self-contained (no
    // cross-thread synchronization is published through these values), and
    // snapshots tolerate transient skew between counters. Atomicity alone
    // suffices for eventual consistency.
    static ALLOCATED: AtomicU64 = AtomicU64::new(0);
    static COUNT: AtomicU64 = AtomicU64::new(0);
    static LIVE: AtomicU64 = AtomicU64::new(0);
    static PEAK: AtomicU64 = AtomicU64::new(0);
    /// Number of `on_dealloc` calls whose size exceeded observed `LIVE` —
    /// i.e. accounting drift (layout mismatch, custom allocator bug, etc.).
    /// Reported by [`report_snapshot_diff`]; nonzero means numbers are suspect.
    static DEALLOC_UNDERFLOWS: AtomicU64 = AtomicU64::new(0);

    const ZERO_ATOMIC: AtomicU64 = AtomicU64::new(0);
    static HIST: [AtomicU64; 16] = [ZERO_ATOMIC; 16];

    static LAST: Mutex<Counters> = Mutex::new(Counters::ZERO);

    #[inline]
    fn bucket(size: usize) -> usize {
        if size == 0 {
            return 0;
        }
        let bits = usize::BITS - size.leading_zeros();
        (bits >> 2).min(15) as usize
    }

    #[inline]
    pub(crate) fn on_alloc(size: usize) {
        let s = size as u64;
        ALLOCATED.fetch_add(s, Ordering::Relaxed);
        COUNT.fetch_add(1, Ordering::Relaxed);
        let new_live = LIVE.fetch_add(s, Ordering::Relaxed) + s;
        PEAK.fetch_max(new_live, Ordering::Relaxed);
        HIST[bucket(size)].fetch_add(1, Ordering::Relaxed);
    }

    #[inline]
    pub(crate) fn on_dealloc(size: usize) {
        let s = size as u64;
        // Atomic saturating sub. Drift (curr < s) means the inner allocator
        // reported mismatched layouts — accounting is already wrong, so we
        // clamp LIVE to 0, bump DEALLOC_UNDERFLOWS, and let report_snapshot_diff
        // surface it. Cannot panic/log from the hook: both may allocate and
        // recurse into this function. The debug_assert catches drift in tests;
        // on CAS retries under contention the closure runs multiple times and
        // the assert re-checks each time, which is intentional.
        let _ = LIVE.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |curr| {
            debug_assert!(
                curr >= s,
                "alloc_track: dealloc size {} exceeds live {}",
                s,
                curr
            );
            if curr < s {
                DEALLOC_UNDERFLOWS.fetch_add(1, Ordering::Relaxed);
            }
            Some(curr.saturating_sub(s))
        });
    }

    pub fn current() -> Counters {
        let mut hist = [0u64; 16];
        for (i, h) in HIST.iter().enumerate() {
            hist[i] = h.load(Ordering::Relaxed);
        }
        Counters {
            total_bytes: ALLOCATED.load(Ordering::Relaxed),
            total_count: COUNT.load(Ordering::Relaxed),
            live_bytes: LIVE.load(Ordering::Relaxed),
            peak_bytes: PEAK.load(Ordering::Relaxed),
            hist,
        }
    }

    pub fn report_snapshot_diff(name: &'static str) {
        let curr = current();
        let prev = {
            let mut last = LAST.lock().unwrap_or_else(|e| e.into_inner());
            std::mem::replace(&mut *last, curr)
        };
        let d_bytes = curr.total_bytes.wrapping_sub(prev.total_bytes);
        let d_count = curr.total_count.wrapping_sub(prev.total_count);
        let d_live = curr.live_bytes as i64 - prev.live_bytes as i64;
        // churn = bytes allocated per byte retained; high = throwaway, ~1 = all retained.
        let churn = if d_live > 0 {
            format!("{:.1}x", d_bytes as f64 / d_live as f64)
        } else {
            "inf".to_string()
        };
        let mut d_hist = [0u64; 16];
        for i in 0..16 {
            d_hist[i] = curr.hist[i].wrapping_sub(prev.hist[i]);
        }
        eprintln!(
            "[alloc] {name:<32} d_bytes=+{} d_live={} churn={} d_n=+{} peak={} hist=[{}]",
            fmt_bytes(d_bytes),
            fmt_signed_bytes(d_live),
            churn,
            fmt_count(d_count),
            fmt_bytes(curr.peak_bytes),
            fmt_hist(&d_hist),
        );
        let drift = DEALLOC_UNDERFLOWS.load(Ordering::Relaxed);
        if drift > 0 {
            eprintln!(
                "[alloc] WARNING: {drift} dealloc underflow(s) observed — live/peak numbers above are suspect."
            );
        }
    }

    const BUCKET_LABELS: [&str; 16] = [
        "<=8B", "<=128B", "<=2K", "<=32K", "<=512K", "<=8M", "<=128M", "<=2G", "<=32G", "<=512G",
        "<=8T", "<=128T", "<=2P", "<=32P", "<=512P", "huge",
    ];

    fn fmt_hist(hist: &[u64; 16]) -> String {
        let mut parts = Vec::new();
        for (i, &n) in hist.iter().enumerate() {
            if n > 0 {
                parts.push(format!("{}:{}", BUCKET_LABELS[i], fmt_count(n)));
            }
        }
        parts.join(" ")
    }

    fn fmt_signed_bytes(n: i64) -> String {
        if n >= 0 {
            format!("+{}", fmt_bytes(n as u64))
        } else {
            format!("-{}", fmt_bytes(n.unsigned_abs()))
        }
    }

    fn fmt_bytes(n: u64) -> String {
        const K: u64 = 1024;
        const M: u64 = K * 1024;
        const G: u64 = M * 1024;
        const T: u64 = G * 1024;
        if n >= T {
            format!("{:.2} TB", n as f64 / T as f64)
        } else if n >= G {
            format!("{:.2} GB", n as f64 / G as f64)
        } else if n >= M {
            format!("{:.2} MB", n as f64 / M as f64)
        } else if n >= K {
            format!("{:.2} KB", n as f64 / K as f64)
        } else {
            format!("{} B", n)
        }
    }

    fn fmt_count(n: u64) -> String {
        let s = n.to_string();
        let bytes = s.as_bytes();
        let mut out = String::with_capacity(s.len() + s.len() / 3);
        for (i, b) in bytes.iter().enumerate() {
            if i > 0 && (bytes.len() - i) % 3 == 0 {
                out.push(',');
            }
            out.push(*b as char);
        }
        out
    }
}

#[cfg(feature = "enabled")]
pub use imp::{
    current,
    report_snapshot_diff,
};

#[cfg(feature = "enabled")]
unsafe impl<A: GlobalAlloc> GlobalAlloc for TrackingAllocator<A> {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = unsafe { self.inner.alloc(layout) };
        if !p.is_null() {
            imp::on_alloc(layout.size());
        }
        p
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        imp::on_dealloc(layout.size());
        unsafe { self.inner.dealloc(ptr, layout) }
    }
}

// --- disabled (no-op) impl -----------------------------------------------

#[cfg(not(feature = "enabled"))]
#[inline(always)]
pub fn current() -> Counters {
    Counters::ZERO
}

#[cfg(not(feature = "enabled"))]
#[inline(always)]
pub fn report_snapshot_diff(_name: &'static str) {}

#[cfg(not(feature = "enabled"))]
unsafe impl<A: GlobalAlloc> GlobalAlloc for TrackingAllocator<A> {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        unsafe { self.inner.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { self.inner.dealloc(ptr, layout) }
    }
}

// --- macro ---------------------------------------------------------------
// `#[cfg]` inside a macro_rules! body is evaluated in the caller's crate,
// so we gate the whole macro definition on our own feature instead.

#[cfg(feature = "enabled")]
#[macro_export]
macro_rules! snap {
    ($name:expr) => {
        $crate::report_snapshot_diff($name)
    };
}

#[cfg(not(feature = "enabled"))]
#[macro_export]
macro_rules! snap {
    ($name:expr) => {
        ()
    };
}
