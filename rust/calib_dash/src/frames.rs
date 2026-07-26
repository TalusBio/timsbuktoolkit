//! Per-batch calibrant snapshots, stored in one preallocated slab.
//!
//! The batch axis is bounded by BYTES, not by a frame count: frame size scales
//! with `n_calibrants`, so a count that is comfortable at 2000 calibrants is
//! not at 200k. The stride is derived once, up front, from quantities all known
//! before the first chunk.
//!
//! Decimation is by stride, not by eviction. A ring buffer would retain the
//! most recent batches and discard the earliest, which is backwards — the early
//! batches are where the fit is volatile and where "when did this stabilize"
//! gets answered.

use std::ops::Range;

/// One heap entry, flattened. `speclib_index` is carried because churn diffing
/// needs a stable identity for a calibrant: RT coordinates are not unique and
/// cannot distinguish "same peptide, re-scored" from "different peptide, same
/// RT".
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibrantPoint {
    pub library_rt: f64,
    pub observed_rt: f64,
    pub score: f64,
    pub speclib_index: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FrameIndex {
    pub chunk: usize,
    /// Speclib indices this chunk covered.
    pub range: Range<usize>,
    offset: usize,
    pub len: usize,
}

pub struct FrameStore {
    slab: Vec<CalibrantPoint>,
    index: Vec<FrameIndex>,
    n_calibrants: usize,
    keep_every: usize,
    /// Number of on-stride spans the slab has room for. This, not
    /// `index.capacity()`, is the bound `record` gates on: `Vec::capacity()`
    /// is only guaranteed to be `>=` the requested amount and may
    /// over-allocate, which would let an on-stride write drift into the
    /// reserved final span or past the end of the slab.
    retained: usize,
    /// Span reserved for the always-kept final frame.
    last_span: usize,
    last_index: Option<FrameIndex>,
    seen: usize,
}

impl FrameStore {
    pub fn new(n_frames: usize, n_calibrants: usize, budget_bytes: usize) -> Self {
        // Defensive floor, not a semantic default: callers are expected to
        // pass a real calibrant count. This only keeps a degenerate zero
        // from producing a zero-sized frame_bytes / divide-by-zero below.
        let n_calibrants = n_calibrants.max(1);
        let frame_bytes = n_calibrants * std::mem::size_of::<CalibrantPoint>();
        let max_frames = (budget_bytes / frame_bytes.max(1)).max(1);
        let keep_every = n_frames.max(1).div_ceil(max_frames).max(1);
        let retained = n_frames.max(1).div_ceil(keep_every);
        // +1 span for the always-kept last frame, which need not land on the
        // stride.
        let slab = vec![
            CalibrantPoint {
                library_rt: 0.0,
                observed_rt: 0.0,
                score: 0.0,
                speclib_index: 0
            };
            (retained + 1) * n_calibrants
        ];
        Self {
            slab,
            index: Vec::with_capacity(retained),
            n_calibrants,
            keep_every,
            retained,
            last_span: retained * n_calibrants,
            last_index: None,
            seen: 0,
        }
    }

    pub fn keep_every(&self) -> usize {
        self.keep_every
    }

    pub fn len(&self) -> usize {
        self.index.len()
    }

    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    pub fn dropped(&self) -> usize {
        self.seen.saturating_sub(self.index.len())
    }

    #[cfg(test)]
    pub fn slab_capacity(&self) -> usize {
        self.slab.capacity()
    }

    #[cfg(test)]
    pub fn index_capacity(&self) -> usize {
        self.index.capacity()
    }

    /// Copy this chunk's heap contents into the slab. Returns whether the frame
    /// was retained on the stride. The final frame is handled by `finish`.
    pub fn record(
        &mut self,
        chunk: usize,
        range: Range<usize>,
        points: impl Iterator<Item = CalibrantPoint>,
    ) -> bool {
        self.seen += 1;
        let on_stride = chunk.is_multiple_of(self.keep_every) && self.index.len() < self.retained;
        let offset = if on_stride {
            self.index.len() * self.n_calibrants
        } else {
            self.last_span
        };
        let mut len = 0;
        for p in points.take(self.n_calibrants) {
            self.slab[offset + len] = p;
            len += 1;
        }
        let idx = FrameIndex {
            chunk,
            range,
            offset,
            len,
        };
        if on_stride {
            self.index.push(idx);
        } else {
            self.last_index = Some(idx);
        }
        on_stride
    }

    /// Promote the reserved span into the index if the final chunk was not
    /// already retained on the stride.
    ///
    /// `last_chunk` must be the chunk passed to the most recent call to
    /// `record`. The reserved slot is only ever promoted when it actually
    /// holds that chunk — never blindly, because a stale non-stride entry
    /// left over from an earlier chunk would otherwise get appended after
    /// chunks already in the index, breaking the monotonic chunk order a
    /// history scrubber depends on.
    pub fn finish(&mut self, last_chunk: usize) {
        if self.index.last().map(|f| f.chunk) == Some(last_chunk) {
            return;
        }
        let Some(idx) = self.last_index.take() else {
            return;
        };
        if idx.chunk == last_chunk {
            self.index.push(idx);
        } else {
            // Only reachable if a caller passes a `last_chunk` that matches
            // neither the index tail nor the reserved slot — i.e. a value
            // that was never actually the last chunk recorded. That is a
            // caller bug: flag it loudly in debug builds, and in release
            // builds drop the stale entry instead of corrupting order.
            debug_assert!(
                false,
                "FrameStore::finish({last_chunk}) called but the reserved slot holds \
                 chunk {}; dropping it rather than corrupting index order",
                idx.chunk
            );
            tracing::warn!(
                last_chunk,
                reserved_chunk = idx.chunk,
                "FrameStore::finish called with a last_chunk matching neither the index \
                 tail nor the reserved slot; dropping the stale reserved entry"
            );
        }
    }

    pub fn frame(&self, i: usize) -> Option<(&FrameIndex, &[CalibrantPoint])> {
        let idx = self.index.get(i)?;
        Some((idx, &self.slab[idx.offset..idx.offset + idx.len]))
    }

    pub fn last(&self) -> Option<(&FrameIndex, &[CalibrantPoint])> {
        self.frame(self.index.len().checked_sub(1)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(i: usize) -> CalibrantPoint {
        CalibrantPoint {
            library_rt: i as f64,
            observed_rt: i as f64 * 2.0,
            score: 1.0,
            speclib_index: i,
        }
    }

    fn budget_for(frames: usize, n_cal: usize) -> usize {
        frames * n_cal * std::mem::size_of::<CalibrantPoint>()
    }

    #[test]
    fn a_generous_budget_keeps_every_frame() {
        let store = FrameStore::new(10, 4, budget_for(100, 4));
        assert_eq!(store.keep_every(), 1);
    }

    #[test]
    fn a_tight_budget_strides() {
        // Room for 3 frames, 12 chunks to record.
        let store = FrameStore::new(12, 4, budget_for(3, 4));
        assert_eq!(store.keep_every(), 4);
    }

    #[test]
    fn retained_frames_never_exceed_the_budget() {
        let mut store = FrameStore::new(12, 4, budget_for(3, 4));
        for chunk in 0..12 {
            store.record(chunk, chunk..chunk + 1, (0..4).map(pt));
        }
        store.finish(11);
        assert!(
            store.len() <= 4,
            "3 strided + 1 reserved last, got {}",
            store.len()
        );
        // Hand-computed for this fixture: chunks 0, 4, 8 land on the stride
        // (keep_every = 4), chunk 11 is promoted as the reserved last frame,
        // and the other 8 of the 12 recorded chunks are dropped.
        assert_eq!(store.dropped(), 8);
    }

    #[test]
    fn the_last_frame_is_always_retained() {
        let mut store = FrameStore::new(10, 4, budget_for(3, 4));
        for chunk in 0..10 {
            store.record(chunk, chunk..chunk + 1, (0..4).map(pt));
        }
        store.finish(9);
        let (idx, _) = store.last().expect("a last frame must exist");
        assert_eq!(idx.chunk, 9);
    }

    #[test]
    fn recorded_points_round_trip() {
        let mut store = FrameStore::new(1, 4, budget_for(1, 4));
        store.record(0, 0..1, (0..4).map(pt));
        let (idx, pts) = store.frame(0).unwrap();
        assert_eq!(idx.len, 4);
        assert_eq!(pts, &[pt(0), pt(1), pt(2), pt(3)]);
    }

    #[test]
    fn a_short_frame_records_its_real_length() {
        // Early chunks have not filled the heap yet.
        let mut store = FrameStore::new(1, 4, budget_for(1, 4));
        store.record(0, 0..1, (0..2).map(pt));
        let (idx, pts) = store.frame(0).unwrap();
        assert_eq!(idx.len, 2);
        assert_eq!(pts.len(), 2);
    }

    #[test]
    fn an_overlong_frame_is_truncated_not_overrun() {
        let mut store = FrameStore::new(1, 4, budget_for(1, 4));
        store.record(0, 0..1, (0..9).map(pt));
        let (_, pts) = store.frame(0).unwrap();
        assert_eq!(pts.len(), 4, "clamped to n_calibrants");
    }

    #[test]
    fn a_budget_smaller_than_one_frame_keeps_exactly_one() {
        let store = FrameStore::new(10, 100, 8);
        assert_eq!(
            store.keep_every(),
            10,
            "one retained frame across ten chunks"
        );
    }

    #[test]
    fn the_slab_is_allocated_once() {
        let mut store = FrameStore::new(12, 4, budget_for(3, 4));
        let slab_cap = store.slab_capacity();
        let index_cap = store.index_capacity();
        for chunk in 0..12 {
            store.record(chunk, chunk..chunk + 1, (0..4).map(pt));
        }
        assert_eq!(store.slab_capacity(), slab_cap, "slab must not reallocate");
        // `record`'s safety bound is `index.len() < self.retained`, gated on
        // an explicit field rather than `index.capacity()` — but if `index`
        // reallocated anyway, that would defeat the "allocate once" promise
        // this store is built on.
        assert_eq!(
            store.index_capacity(),
            index_cap,
            "index must not reallocate"
        );
    }

    #[test]
    fn finish_drops_a_mismatched_stale_entry_without_corrupting_order() {
        // keep_every = 4, retained = 3: chunks 0, 4, 8 land on stride and
        // fill the index outright, while chunks 5..7 leave a stale
        // non-stride entry (chunk 7) sitting in the reserved slot that is
        // never promoted, since chunk 8 was already retained on the stride.
        let mut store = FrameStore::new(12, 4, budget_for(3, 4));
        for chunk in 0..9 {
            store.record(chunk, chunk..chunk + 1, (0..4).map(pt));
        }
        // 20 matches neither the index tail (chunk 8) nor the reserved slot
        // (chunk 7): a caller bug. The old code would push the stale chunk-7
        // entry after chunk 8 regardless, producing out-of-order chunks.
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            store.finish(20);
        }));
        assert!(
            result.is_err(),
            "a last_chunk matching neither slot must be flagged loudly in a debug build"
        );
        let chunks: Vec<usize> = (0..store.len())
            .map(|i| store.frame(i).unwrap().0.chunk)
            .collect();
        assert!(
            chunks.windows(2).all(|w| w[0] < w[1]),
            "chunk order must stay monotonic even when finish() is misused, got {chunks:?}"
        );
    }
}
