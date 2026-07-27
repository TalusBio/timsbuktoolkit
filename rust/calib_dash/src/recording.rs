//! Owned copies of the borrowed `FitEvent`s, filled in place.
//!
//! `weights` are `f32` for storage and display only — every metric is computed
//! from the `f64` values in the live events, never from this downcast copy.
//! At `bins = 100` a recording is ~50 KB; the `Grid` the fit already holds is
//! ~40 bytes per cell, so this is roughly a tenth of what the fit itself spends.

use calibrt::{
    FitEvent,
    FitObserver,
    GridGeom,
    Point,
    RidgeMeasurement,
};

#[derive(Debug, Clone)]
pub struct DpDecision {
    pub i: usize,
    pub library: f64,
    pub observed: f64,
    pub chose: Option<usize>,
    pub acc_weight: f64,
    pub considered: Vec<(usize, f64)>,
}

/// `Clone` so the batch scrubber (`App::set_scrub_recording`, driven by
/// `CalibDash::sync_scrub`) can hand the Fit tab its own owned copy of a
/// replayed frame's recording without fighting `refit_recording`'s
/// reuse-in-place allocation, which stays live for the *next* scrub — an
/// extra allocation, but only once per keypress while a user is actively
/// browsing history, never on the per-batch hot path.
#[derive(Clone)]
pub struct FitRecording {
    geom: GridGeom,
    weights: Vec<f32>,
    suppressed: Vec<bool>,
    /// Grid indices of the assembled path (DP chain plus greedy tails). The
    /// `bins` capacity is a hint, not a bound — weight ties can push the
    /// survivor count past it, as `calibrt`'s survivor debug assert spells
    /// out. Growing past a hinted capacity is a reallocation, not a bug.
    path_indices: Vec<usize>,
    /// The DP-chosen segment within `path_indices`: `path_indices[..dp_range.start]`
    /// and `path_indices[dp_range.end..]` were greedily attached by Pass 2,
    /// not scored by the DP recurrence.
    dp_range: std::ops::Range<usize>,
    curve: Vec<Point>,
    ridge: Vec<RidgeMeasurement>,
    /// One entry per DP node visited. Same `bins`-capacity hint as
    /// `path_indices` above.
    dp: Vec<DpDecision>,
}

impl FitRecording {
    pub fn new(bins: usize) -> Self {
        Self {
            geom: GridGeom {
                bins,
                x_range: (0.0, 1.0),
                y_range: (0.0, 1.0),
            },
            weights: vec![0.0; bins * bins],
            suppressed: vec![false; bins * bins],
            // Capacity hint only — see the field doc comment.
            path_indices: Vec::with_capacity(bins),
            dp_range: 0..0,
            curve: Vec::with_capacity(bins),
            ridge: Vec::with_capacity(bins),
            dp: Vec::with_capacity(bins),
        }
    }

    /// Clears everything except `weights`/`suppressed` — split out so
    /// `FitStarted` can reuse it after (conditionally) resizing those two for
    /// a changed `bins`, without also paying for the fill/clear this helper
    /// skips.
    fn reset_body_only(&mut self) {
        self.path_indices.clear();
        self.dp_range = 0..0;
        self.curve.clear();
        self.ridge.clear();
        self.dp.clear();
    }

    pub(crate) fn geom(&self) -> GridGeom {
        self.geom
    }

    pub(crate) fn path_indices(&self) -> &[usize] {
        &self.path_indices
    }

    /// The DP-chosen slice of `path_indices`. Indices outside it were attached
    /// by Pass 2's greedy extension.
    pub(crate) fn dp_range(&self) -> std::ops::Range<usize> {
        self.dp_range.clone()
    }

    pub(crate) fn curve(&self) -> &[Point] {
        &self.curve
    }

    pub(crate) fn ridge(&self) -> &[RidgeMeasurement] {
        &self.ridge
    }

    pub(crate) fn dp(&self) -> &[DpDecision] {
        &self.dp
    }

    pub(crate) fn weight(&self, row: usize, col: usize) -> f32 {
        self.weights
            .get(row * self.geom.bins + col)
            .copied()
            .unwrap_or(0.0)
    }

    pub(crate) fn is_suppressed(&self, row: usize, col: usize) -> bool {
        self.suppressed
            .get(row * self.geom.bins + col)
            .copied()
            .unwrap_or(false)
    }
}

/// Grid-bin index of `v` within `range`, at `bins` bins: the one place the
/// value-to-bin arithmetic is written, shared by `FitRecording`'s own
/// row/column placement and by `ui.rs`'s overlay placement, so a mark can
/// never land in a different bin than the weight it is marking.
///
/// Callers own the domain: `bins` must be nonzero and `range` a finite,
/// positive-width span. `ui.rs`'s `bin_of` wrapper is what screens for those
/// before calling in.
pub(crate) fn bin_index(v: f64, range: (f64, f64), bins: usize) -> usize {
    let (lo, hi) = range;
    let span = hi - lo;
    (((v - lo) / span * bins as f64) as usize).min(bins - 1)
}

impl FitObserver for FitRecording {
    fn on_event(&mut self, ev: FitEvent<'_>) {
        match ev {
            FitEvent::FitStarted { geom } => {
                // A geometry change means the caller re-fit at a different
                // `bins`; resize rather than silently mis-indexing.
                if geom.bins != self.geom.bins {
                    self.weights = vec![0.0; geom.bins * geom.bins];
                    self.suppressed = vec![false; geom.bins * geom.bins];
                }
                self.geom = geom;
                // `Suppressed` below only ever sets bits for cells that are
                // suppressed *this* fit; it never clears a bit for a cell
                // that used to be suppressed and now survives. Without this,
                // a re-fit at the same `bins` (the common case — no
                // reallocation above) would carry forward stale suppression
                // flags from the previous fit.
                self.suppressed.fill(false);
                self.reset_body_only();
            }
            FitEvent::GridReady { cells } => {
                // Every cell is rewritten unconditionally: `cells` always has
                // length `bins * bins`, so this alone keeps `weights` correct
                // across re-fits without a separate clear.
                for (i, n) in cells.iter().enumerate() {
                    self.weights[i] = n.center.weight as f32;
                }
            }
            FitEvent::Suppressed { cells } => {
                for (slot, n) in self.suppressed.iter_mut().zip(cells) {
                    if n.suppressed {
                        *slot = true;
                    }
                }
            }
            FitEvent::DpNode {
                i,
                node,
                chose,
                acc_weight,
                considered,
            } => {
                self.dp.push(DpDecision {
                    i,
                    library: node.center.library,
                    observed: node.center.observed,
                    chose,
                    acc_weight,
                    considered: considered.to_vec(),
                });
            }
            FitEvent::PathFound { path, dp_range, .. } => {
                debug_assert!(
                    dp_range.end <= path.len(),
                    "dp_range must fall within the assembled path"
                );
                self.dp_range = dp_range;
                // path_indices are grid indices, derived from the point's cell.
                for p in path {
                    let col = bin_index(p.library, self.geom.x_range, self.geom.bins);
                    let row = bin_index(p.observed, self.geom.y_range, self.geom.bins);
                    self.path_indices.push(row * self.geom.bins + col);
                }
            }
            FitEvent::CurveFit { curve } => self.curve.extend_from_slice(curve.points()),
            FitEvent::RidgeMeasured { widths } => self.ridge.extend_from_slice(widths),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use calibrt::{
        CalibrationState,
        LibraryRT,
        ObserveOpts,
        ObservedRTSeconds,
    };

    fn diagonal_state(bins: usize) -> CalibrationState {
        let mut s = CalibrationState::new(bins, (0.0, bins as f64), (0.0, bins as f64), 5).unwrap();
        let pts: Vec<_> = (0..bins)
            .map(|i| {
                let v = i as f64 + 0.5;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0 + i as f64)
            })
            .collect();
        s.update(pts.into_iter()).unwrap();
        s
    }

    #[test]
    fn recording_captures_geometry_weights_and_path() {
        let mut s = diagonal_state(10);
        // A genuinely off-diagonal point (library != observed): every other
        // point in this fixture has library == observed, so a row/col swap in
        // `col_of`/`row_of` would pass unnoticed without this one. Weight 0.5
        // is too light to be a row/column max next to the diagonal's weight-8
        // entry at row 7 and weight-3 entry at col 2, so it is suppressed and
        // does not disturb the path/curve length assertions below — it only
        // needs to show up in the raw weight grid.
        s.update(std::iter::once((
            LibraryRT(2.5),
            ObservedRTSeconds(7.5),
            0.5,
        )))
        .unwrap();
        let mut rec = FitRecording::new(10);
        s.fit_with(&mut rec, ObserveOpts::NONE);

        assert_eq!(rec.geom().bins, 10);
        // The diagonal cell (row i, col i) carries weight 1 + i.
        assert!((rec.weight(3, 3) - 4.0).abs() < 1e-6);
        assert_eq!(rec.weight(0, 5), 0.0, "off-diagonal cells are empty");
        assert!(
            (rec.weight(7, 2) - 0.5).abs() < 1e-6,
            "row = observed's bin, col = library's bin"
        );
        assert_eq!(
            rec.weight(2, 7),
            0.0,
            "a row/col swap would put the extra point here instead"
        );
        assert_eq!(rec.path_indices().len(), 10);
        assert_eq!(rec.curve().len(), 10);
    }

    #[test]
    fn recording_captures_the_suppression_mask() {
        let mut s = diagonal_state(10);
        let mut rec = FitRecording::new(10);
        s.fit_with(&mut rec, ObserveOpts::NONE);
        // Every diagonal cell survives; every empty cell is suppressed.
        assert!(!rec.is_suppressed(4, 4));
        assert!(rec.is_suppressed(0, 5));
    }

    #[test]
    fn dp_decisions_are_recorded_only_when_enabled() {
        let mut s = diagonal_state(10);

        let mut off = FitRecording::new(10);
        s.fit_with(&mut off, ObserveOpts::NONE);
        assert!(off.dp().is_empty());

        let mut on = FitRecording::new(10);
        s.reset();
        let pts: Vec<_> = (0..10)
            .map(|i| {
                let v = i as f64 + 0.5;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0 + i as f64)
            })
            .collect();
        s.update(pts.into_iter()).unwrap();
        s.fit_with(&mut on, ObserveOpts { dp_nodes: true });
        assert_eq!(on.dp().len(), 10);
        assert!(
            on.dp().iter().any(|d| d.chose.is_some()),
            "some node picked a predecessor"
        );
    }

    #[test]
    fn a_stale_suppression_bit_does_not_survive_a_refit_at_unchanged_bins() {
        // Fit 1: A at (row 1, col 1) is dominated by B in the same row (B's
        // weight 9 beats A's weight 2), so A fails the row-max check and is
        // suppressed regardless of its column.
        let mut s = CalibrationState::new(3, (0.0, 3.0), (0.0, 3.0), 2).unwrap();
        s.update(
            [
                (LibraryRT(1.5), ObservedRTSeconds(1.5), 2.0),
                (LibraryRT(2.5), ObservedRTSeconds(1.5), 9.0),
            ]
            .into_iter(),
        )
        .unwrap();
        let mut rec = FitRecording::new(3);
        s.fit_with(&mut rec, ObserveOpts::NONE);
        assert!(
            rec.is_suppressed(1, 1),
            "A is dominated by B in its row on the first fit"
        );

        // Fit 2, same `CalibrationState` (so `bins` is unchanged and
        // `FitRecording` does not reallocate its mask): B is gone and A is
        // now the sole occupant of its row and column, so it survives.
        s.reset();
        s.update(std::iter::once((
            LibraryRT(1.5),
            ObservedRTSeconds(1.5),
            5.0,
        )))
        .unwrap();
        s.fit_with(&mut rec, ObserveOpts::NONE);
        assert!(
            !rec.is_suppressed(1, 1),
            "A survives alone in its row/col on the second fit — a stale bit \
             from the first fit must not linger"
        );
    }

    #[test]
    fn a_failed_suppression_records_an_empty_path() {
        // All weights below 1.0 trip calibrt's `max_in_row` initialization, so
        // everything is suppressed and no path exists. The recording must
        // survive it rather than panic or half-fill.
        let mut s = CalibrationState::new(4, (0.0, 4.0), (0.0, 4.0), 2).unwrap();
        s.update((0..4).map(|i| {
            (
                LibraryRT(i as f64 + 0.5),
                ObservedRTSeconds(i as f64 + 0.5),
                0.1,
            )
        }))
        .unwrap();
        let mut rec = FitRecording::new(4);
        s.fit_with(&mut rec, ObserveOpts::NONE);
        assert!(rec.path_indices().is_empty());
        assert!(rec.curve().is_empty());
    }
}
