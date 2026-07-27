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

/// Default slab budget for a live Phase 1 run, and the value
/// `CALIB_DASH_FRAME_BUDGET_MB` overrides. A run has as many frames as it has
/// scoring chunks, so this is the budget that decides `keep_every`: too small
/// and the history scrubber only has a handful of batches to step through.
pub const DEFAULT_RUN_BUDGET_BYTES: usize = 64 * 1024 * 1024;

/// Slab budget for replaying a single saved batch (`calib_dash <file>`).
/// Deliberately far smaller than [`DEFAULT_RUN_BUDGET_BYTES`]: there is
/// exactly one frame to hold, so nothing here ever strides and the budget only
/// has to clear one frame's worth of points.
pub const REPLAY_BUDGET_BYTES: usize = 1 << 20;

/// One heap entry, flattened. `speclib_index` is carried because churn diffing
/// needs a stable identity for a calibrant: RT coordinates are not unique and
/// cannot distinguish "same peptide, re-scored" from "different peptide, same
/// RT".
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibrantPoint {
    pub library_rt: f64,
    pub observed_rt: f64,
    pub speclib_index: usize,
}

struct FrameIndex {
    chunk: usize,
    offset: usize,
    len: usize,
}

/// What the slab kept, and at what cost. A named struct because all three are
/// `usize` and a tuple would let a caller swap two of them silently.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameSummary {
    /// Frames currently in the index.
    pub retained: usize,
    pub keep_every: usize,
    pub dropped: usize,
}

pub struct FrameStore {
    slab: Vec<CalibrantPoint>,
    index: Vec<FrameIndex>,
    n_calibrants: usize,
    keep_every: usize,
    /// Number of on-stride spans the slab has room for. This, not
    /// `index.capacity()`, is the bound `record` gates on, since
    /// `Vec::capacity()` may over-allocate past what the slab has room for.
    retained: usize,
    /// Span reserved for the always-kept final frame.
    last_span: usize,
    last_index: Option<FrameIndex>,
    seen: usize,
}

impl FrameStore {
    pub fn new(n_frames: usize, n_calibrants: usize, budget_bytes: usize) -> Self {
        let n_calibrants = n_calibrants.max(1);
        let frame_bytes = n_calibrants * std::mem::size_of::<CalibrantPoint>();
        let max_frames = (budget_bytes / frame_bytes).max(1);
        let keep_every = n_frames.max(1).div_ceil(max_frames).max(1);
        let retained = n_frames.max(1).div_ceil(keep_every);
        // +1 span for the always-kept last frame, which need not land on the
        // stride.
        let slab = vec![
            CalibrantPoint {
                library_rt: 0.0,
                observed_rt: 0.0,
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

    pub fn summary(&self) -> FrameSummary {
        FrameSummary {
            retained: self.index.len(),
            keep_every: self.keep_every,
            dropped: self.seen.saturating_sub(self.index.len()),
        }
    }

    /// Copy this chunk's heap contents into the slab. The final frame is
    /// handled by `finish`.
    pub fn record(&mut self, chunk: usize, points: impl Iterator<Item = CalibrantPoint>) {
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
        let idx = FrameIndex { chunk, offset, len };
        if on_stride {
            self.index.push(idx);
        } else {
            self.last_index = Some(idx);
        }
    }

    /// Promote the reserved span into the index if the final chunk was not
    /// already retained on the stride. A reserved entry older than the index
    /// tail is stale and dropped: the index must stay in chunk order for the
    /// history scrubber.
    pub fn finish(&mut self) {
        if let Some(idx) = self.last_index.take()
            && self.index.last().is_none_or(|f| f.chunk < idx.chunk)
        {
            self.index.push(idx);
        }
    }

    /// The `i`th retained frame as `(chunk, points)`.
    pub fn frame(&self, i: usize) -> Option<(usize, &[CalibrantPoint])> {
        let idx = self.index.get(i)?;
        Some((idx.chunk, &self.slab[idx.offset..idx.offset + idx.len]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(i: usize) -> CalibrantPoint {
        CalibrantPoint {
            library_rt: i as f64,
            observed_rt: i as f64 * 2.0,
            speclib_index: i,
        }
    }

    /// `pt`, tagged with the chunk that recorded it, so one frame's span cannot
    /// pass for another's.
    fn pt_in(chunk: usize, i: usize) -> CalibrantPoint {
        CalibrantPoint {
            speclib_index: chunk * 100 + i,
            ..pt(i)
        }
    }

    fn budget_for(frames: usize, n_cal: usize) -> usize {
        frames * n_cal * std::mem::size_of::<CalibrantPoint>()
    }

    #[test]
    fn keep_every_strides_to_fit_the_budget() {
        // (n_frames, n_calibrants, budget in whole frames, expected keep_every)
        let cases = [
            (10, 4, 100, 1),
            (12, 4, 3, 4),
            // A budget below one frame still keeps one, striding over all ten.
            (10, 100, 0, 10),
        ];
        for (n_frames, n_cal, budget_frames, expected) in cases {
            let store = FrameStore::new(n_frames, n_cal, budget_for(budget_frames, n_cal));
            assert_eq!(
                store.summary().keep_every,
                expected,
                "{n_frames} frames of {n_cal} calibrants in a {budget_frames}-frame budget"
            );
        }
    }

    /// Retained frames exceed the *stride* budget by exactly one — the reserved
    /// final span, which is what makes the last chunk always replayable however
    /// the stride falls. Points are asserted per frame, and made
    /// chunk-distinguishable to do it: with the same `pt(i)` in every frame the
    /// whole slab looks alike, so a `record` offset off by a frame would hand
    /// the scrubber a different frame's points unnoticed.
    #[test]
    fn the_last_frame_is_always_retained_with_its_own_points() {
        let mut store = FrameStore::new(12, 4, budget_for(3, 4));
        for chunk in 0..12 {
            store.record(chunk, (0..4).map(|i| pt_in(chunk, i)));
        }
        store.finish();

        // Hand-computed for this fixture: chunks 0, 4, 8 land on the stride
        // (keep_every = 4), chunk 11 is promoted as the reserved last frame,
        // and the other 8 of the 12 recorded chunks are dropped.
        let summary = store.summary();
        assert_eq!(summary.retained, 4, "3 strided + 1 reserved last");
        assert_eq!(summary.dropped, 8);

        let chunks: Vec<usize> = (0..summary.retained)
            .map(|i| store.frame(i).expect("retained frames are readable").0)
            .collect();
        assert_eq!(
            chunks,
            vec![0, 4, 8, 11],
            "the final chunk is the last frame"
        );
        for i in 0..summary.retained {
            let (chunk, pts) = store.frame(i).unwrap();
            let expected: Vec<CalibrantPoint> = (0..4).map(|j| pt_in(chunk, j)).collect();
            assert_eq!(
                pts,
                expected.as_slice(),
                "frame {i} must replay chunk {chunk}'s own points"
            );
        }
    }

    /// The recorded length is the *iterator's*, clamped to `n_calibrants`:
    /// early chunks have not filled the heap yet, and an overlong batch must be
    /// truncated rather than overrun the frame's span.
    #[test]
    fn recorded_points_round_trip() {
        // (points offered, expected recorded length)
        let cases = [(4, 4), (2, 2), (9, 4)];
        for (offered, expected) in cases {
            let mut store = FrameStore::new(1, 4, budget_for(1, 4));
            store.record(0, (0..offered).map(pt));
            let (chunk, pts) = store.frame(0).unwrap();
            assert_eq!(chunk, 0);
            let want: Vec<CalibrantPoint> = (0..expected).map(pt).collect();
            assert_eq!(
                pts,
                want.as_slice(),
                "{offered} points offered must record as {expected}"
            );
        }
    }

    /// `n_calibrants == 0` is a real shape (`calib_dash <file>` derives it from
    /// the snapshot's point count, and a snapshot can be empty), and so is a
    /// zero budget. `n_calibrants.max(1)` is what keeps a frame from being
    /// zero-width: every frame holds at least one point, and the frame-size
    /// divisor `new` computes from it stays nonzero.
    #[test]
    fn a_store_asked_for_no_calibrant_slots_still_keeps_one() {
        let mut store = FrameStore::new(4, 0, 0);
        store.record(0, (0..4).map(pt));
        store.finish();
        let (chunk, pts) = store.frame(0).expect("chunk 0 lands on the stride");
        assert_eq!(chunk, 0);
        assert_eq!(
            pts,
            &[pt(0)],
            "a zero-calibrant store still records into one slot per frame"
        );
    }

    /// The dashboard's "allocate once at startup, never during a run"
    /// promise, for the one allocation big enough to matter. Pinned by
    /// pointer identity rather than `Vec::capacity()` — capacity equality
    /// cannot distinguish "reused the same allocation" from "reallocated and
    /// landed on the same capacity anyway", and `frame()`'s own slice is
    /// already a view into the slab, so no test-only accessor is needed to
    /// see it.
    #[test]
    fn the_slab_is_allocated_once() {
        let mut store = FrameStore::new(12, 4, budget_for(3, 4));
        store.record(0, (0..4).map(pt));
        let base = store
            .frame(0)
            .expect("chunk 0 lands on the stride")
            .1
            .as_ptr();
        for chunk in 1..12 {
            store.record(chunk, (0..4).map(pt));
        }
        store.finish();
        assert_eq!(
            store.frame(0).expect("frame 0 is still there").1.as_ptr(),
            base,
            "the slab must not reallocate"
        );
    }

    #[test]
    fn a_stale_reserved_slot_is_dropped_not_appended_out_of_order() {
        // keep_every = 4, retained = 3: chunks 0, 4, 8 land on the stride and
        // fill the index outright, leaving chunk 7 in the reserved slot. It is
        // older than the index tail, so promoting it would break chunk order.
        let mut store = FrameStore::new(12, 4, budget_for(3, 4));
        for chunk in 0..9 {
            store.record(chunk, (0..4).map(pt));
        }
        store.finish();
        let chunks: Vec<usize> = (0..store.summary().retained)
            .map(|i| store.frame(i).unwrap().0)
            .collect();
        assert_eq!(chunks, vec![0, 4, 8], "chunk order must stay monotonic");
    }
}
