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
    /// Span reserved for the always-kept final frame.
    last_span: usize,
    last_index: Option<FrameIndex>,
    seen: usize,
}

impl FrameStore {
    pub fn new(n_frames: usize, n_calibrants: usize, budget_bytes: usize) -> Self {
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

    /// Copy this chunk's heap contents into the slab. Returns whether the frame
    /// was retained on the stride. The final frame is handled by `finish`.
    pub fn record(
        &mut self,
        chunk: usize,
        range: Range<usize>,
        points: impl Iterator<Item = CalibrantPoint>,
    ) -> bool {
        self.seen += 1;
        let on_stride =
            chunk.is_multiple_of(self.keep_every) && self.index.len() < self.index.capacity();
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
    pub fn finish(&mut self, last_chunk: usize) {
        if self.index.last().map(|f| f.chunk) == Some(last_chunk) {
            return;
        }
        if let Some(idx) = self.last_index.take() {
            self.index.push(idx);
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
        assert_eq!(store.dropped(), 12 - store.len());
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
        let cap = store.slab_capacity();
        for chunk in 0..12 {
            store.record(chunk, chunk..chunk + 1, (0..4).map(pt));
        }
        assert_eq!(store.slab_capacity(), cap, "recording must not reallocate");
    }
}
