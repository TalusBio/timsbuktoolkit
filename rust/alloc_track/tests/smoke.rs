//! End-to-end smoke test. Runs only with `--features enabled` because it
//! installs `TrackingAllocator` as the global allocator for the test binary.

#![cfg(feature = "enabled")]

use alloc_track::{
    Counters,
    TrackingAllocator,
    current,
};

#[global_allocator]
static GLOBAL: TrackingAllocator = TrackingAllocator::new();

#[test]
fn counts_bytes_for_vec() {
    let before = current();

    let v: Vec<u8> = vec![0u8; 1_000_000];
    let mid = current();
    assert!(
        mid.total_bytes >= before.total_bytes + 1_000_000,
        "expected total_bytes to grow by at least 1_000_000, got {} -> {}",
        before.total_bytes,
        mid.total_bytes,
    );
    assert!(mid.total_count > before.total_count);
    assert!(mid.peak_bytes >= mid.live_bytes);

    drop(v);
    let after = current();
    assert!(
        after.live_bytes <= mid.live_bytes,
        "live_bytes should shrink after drop: {} -> {}",
        mid.live_bytes,
        after.live_bytes,
    );
    assert!(
        after.peak_bytes >= mid.peak_bytes,
        "peak is monotone: {} -> {}",
        mid.peak_bytes,
        after.peak_bytes,
    );
}

#[test]
fn counters_zero_const() {
    let z = Counters::ZERO;
    assert_eq!(z.total_bytes, 0);
    assert_eq!(z.hist, [0u64; 16]);
}
