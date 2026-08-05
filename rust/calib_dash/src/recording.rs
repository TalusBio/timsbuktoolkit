//! Owned copies of the borrowed `FitEvent`s, filled in place, plus the fit's
//! products copied off the state by [`FitRecording::set_fit`].
//!
//! `weights` are `f32` for storage and display only — every metric is computed
//! from the `f64` values in the live events, never from this downcast copy.

use calibrt::{
    CalibrationCurve,
    CalibrationState,
    FitEvent,
    FitObserver,
    GridGeom,
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

#[derive(Clone)]
pub struct FitRecording {
    geom: GridGeom,
    weights: Vec<f32>,
    suppressed: Vec<bool>,
    /// See [`calibrt::CalibrationState::path_indices`]. The `bins` capacity is a
    /// hint, not a bound — weight ties can push the survivor count past it.
    path_indices: Vec<usize>,
    /// See [`calibrt::CalibrationState::dp_range`].
    dp_range: std::ops::Range<usize>,
    /// The curve itself, so the overlay predicts through it instead of
    /// re-deriving calibrt's interpolation.
    curve: Option<CalibrationCurve>,
    ridge: Vec<RidgeMeasurement>,
    /// One entry per DP node visited. Same `bins`-capacity hint as
    /// `path_indices` above.
    dp: Vec<DpDecision>,
    /// Set when a fit starts, cleared by [`Self::set_fit`]. The fit's products are
    /// not events, so a caller that fits without calling `set_fit` leaves the four
    /// product buffers holding the *previous* fit — a wrong panel, not a blank
    /// one. Every product accessor refuses to read while this is set.
    unrecorded_fit: bool,
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
            curve: None,
            ridge: Vec::with_capacity(bins),
            dp: Vec::with_capacity(bins),
            // A recording that has seen no fit has nothing outstanding: Phase 1
            // draws an empty Fit tab before the first fit runs.
            unrecorded_fit: false,
        }
    }

    pub(crate) fn geom(&self) -> GridGeom {
        self.geom
    }

    /// Panics when a fit ran without a following [`Self::set_fit`]. Guards each
    /// accessor below, so the fault surfaces at the reader that would have drawn
    /// the previous fit's data as this one's.
    #[track_caller]
    fn expect_fit_recorded(&self) {
        assert!(
            !self.unrecorded_fit,
            "FitRecording holds an unrecorded fit: `fit_with` ran without a following \
             `set_fit`, so the path, curve and ridge are the previous fit's"
        );
    }

    #[track_caller]
    pub(crate) fn path_indices(&self) -> &[usize] {
        self.expect_fit_recorded();
        &self.path_indices
    }

    /// See [`calibrt::CalibrationState::dp_range`].
    #[track_caller]
    pub(crate) fn dp_range(&self) -> std::ops::Range<usize> {
        self.expect_fit_recorded();
        self.dp_range.clone()
    }

    #[track_caller]
    pub(crate) fn curve(&self) -> Option<&CalibrationCurve> {
        self.expect_fit_recorded();
        self.curve.as_ref()
    }

    #[track_caller]
    pub(crate) fn ridge(&self) -> &[RidgeMeasurement] {
        self.expect_fit_recorded();
        &self.ridge
    }

    /// Copy the fit's products in: the path, its DP segment, the curve and the
    /// ridge are not events, the state carries them after `fit_with` returns.
    /// All four are replaced, so this is safe to call after any fit.
    pub fn set_fit(&mut self, state: &CalibrationState) {
        self.path_indices.clear();
        self.path_indices.extend_from_slice(state.path_indices());
        self.dp_range = state.dp_range();
        self.curve = state.curve().cloned();
        self.ridge.clear();
        self.ridge.extend_from_slice(state.ridge_widths());
        self.unrecorded_fit = false;
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

impl FitObserver for FitRecording {
    fn on_event(&mut self, ev: FitEvent<'_>) {
        match ev {
            FitEvent::FitStarted { geom, cells } => {
                // A geometry change means the caller re-fit at a different
                // `bins`; resize rather than silently mis-indexing.
                if geom.bins != self.geom.bins {
                    self.weights = vec![0.0; geom.bins * geom.bins];
                    self.suppressed = vec![false; geom.bins * geom.bins];
                }
                self.geom = geom;
                // The products this fit will leave on the state are not events;
                // they arrive through `set_fit` or not at all.
                self.unrecorded_fit = true;
                // `Suppressed` only ever sets bits, never clears them.
                self.suppressed.fill(false);
                self.dp.clear();
                // `cells` is always `bins * bins` long, so every cell is
                // rewritten and `weights` needs no separate clear.
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

    /// Mirrors `calibrt`'s own `diagonal_state` test fixture, parameterized on `bins`.
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
        rec.set_fit(&s);

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
        assert_eq!(rec.curve().unwrap().points().len(), 10);
    }

    /// `set_fit` is manual, so forgetting it has to be loud: the products would
    /// otherwise read as the previous fit's, which draws a plausible wrong panel.
    #[test]
    #[should_panic(expected = "unrecorded fit")]
    fn reading_the_products_of_an_unrecorded_fit_panics() {
        let mut s = diagonal_state(10);
        let mut rec = FitRecording::new(10);
        s.fit_with(&mut rec, ObserveOpts::NONE);
        rec.set_fit(&s);
        assert_eq!(rec.path_indices().len(), 10, "the recorded fit reads fine");

        // A second fit, this time without the `set_fit` that copies its products.
        s.fit_with(&mut rec, ObserveOpts::NONE);
        let _ = rec.path_indices();
    }

    /// Nothing has been fit, so nothing is outstanding — Phase 1 draws an empty
    /// Fit tab before the first fit runs.
    #[test]
    fn a_recording_with_no_fit_reads_as_empty() {
        let rec = FitRecording::new(10);
        assert!(rec.path_indices().is_empty());
        assert!(rec.curve().is_none());
        assert!(rec.ridge().is_empty());
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

    /// `CalibDash` reuses one `FitRecording` for every refit, and a scrubbed
    /// frame or a replayed snapshot can arrive at a different `bins`: the grid
    /// buffers must be resized (an unresized one reads as a blank heatmap, not a
    /// panic) and the fit's products must hold the second fit's entries only.
    #[test]
    fn a_refit_at_a_different_bins_resizes_the_grid_and_shows_only_the_second_fit() {
        // Fit 1, a 3x3 grid: 3 path indices, 3 curve points, 3 measurements.
        let mut small = diagonal_state(3);
        let mut rec = FitRecording::new(3);
        small.fit_with(&mut rec, ObserveOpts::NONE);
        rec.set_fit(&small);
        assert_eq!(
            (
                rec.geom().bins,
                rec.path_indices().len(),
                rec.curve().unwrap().points().len(),
                rec.ridge().len()
            ),
            (3, 3, 3, 3),
            "sanity: the small fit filled every buffer"
        );

        // Fit 2, a 10x10 grid through the same recording.
        let mut big = diagonal_state(10);
        big.fit_with(&mut rec, ObserveOpts::NONE);
        rec.set_fit(&big);

        assert_eq!(rec.geom().bins, 10, "the geometry follows the refit");
        // Row 9 exists only in a resized buffer: at 3x3 the flat index 9*10+9
        // is past the end of a 9-cell grid.
        assert!(
            (rec.weight(9, 9) - 10.0).abs() < 1e-6,
            "the diagonal's last cell carries weight 1 + 9, got {}",
            rec.weight(9, 9)
        );
        assert!(
            rec.is_suppressed(9, 5),
            "an empty cell in the resized mask is suppressed"
        );
        assert!(
            !rec.is_suppressed(9, 9),
            "the diagonal's last cell survives suppression"
        );
        assert_eq!(
            (
                rec.path_indices().len(),
                rec.curve().unwrap().points().len(),
                rec.ridge().len()
            ),
            (10, 10, 10),
            "the refit's own entries only — 13 of each would mean the first \
             fit's leaked through"
        );
    }
}
