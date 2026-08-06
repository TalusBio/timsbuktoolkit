//! An owned copy of everything a finished [`calibrt::CalibrationState`] holds
//! that a panel draws: the grid it fit on, and the fit's products.
//!
//! `weights` are `f32` for storage and display only — every metric is computed
//! from the `f64` values on the state, never from this downcast copy.

use calibrt::{
    CalibrationCurve,
    CalibrationState,
    GridGeom,
    RidgeMeasurement,
};

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
        }
    }

    /// Everything the panels draw, read off a state whose fit has finished. A
    /// state whose fit failed carries an empty path, curve and ridge, and this
    /// copies that emptiness through.
    pub fn from_state(state: &CalibrationState) -> Self {
        let cells = state.grid_cells();
        Self {
            geom: GridGeom {
                bins: state.grid_bins(),
                x_range: state.grid_x_range(),
                y_range: state.grid_y_range(),
            },
            weights: cells.iter().map(|n| n.center.weight as f32).collect(),
            suppressed: cells.iter().map(|n| n.suppressed).collect(),
            path_indices: state.path_indices().to_vec(),
            dp_range: state.dp_range(),
            curve: state.curve().cloned(),
            ridge: state.ridge_widths().to_vec(),
        }
    }

    pub(crate) fn geom(&self) -> GridGeom {
        self.geom
    }

    pub(crate) fn path_indices(&self) -> &[usize] {
        &self.path_indices
    }

    /// See [`calibrt::CalibrationState::dp_range`].
    pub(crate) fn dp_range(&self) -> std::ops::Range<usize> {
        self.dp_range.clone()
    }

    pub(crate) fn curve(&self) -> Option<&CalibrationCurve> {
        self.curve.as_ref()
    }

    pub(crate) fn ridge(&self) -> &[RidgeMeasurement] {
        &self.ridge
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

#[cfg(test)]
mod tests {
    use super::*;
    use calibrt::{
        CalibrationState,
        LibraryRT,
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
        s.fit();
        let rec = FitRecording::from_state(&s);

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
        s.fit();
        assert!(
            FitRecording::from_state(&s).is_suppressed(1, 1),
            "A is dominated by B in its row on the first fit"
        );

        // Fit 2, same `CalibrationState` at the same `bins`: B is gone and A is
        // now the sole occupant of its row and column, so it survives.
        s.reset();
        s.update(std::iter::once((
            LibraryRT(1.5),
            ObservedRTSeconds(1.5),
            5.0,
        )))
        .unwrap();
        s.fit();
        assert!(
            !FitRecording::from_state(&s).is_suppressed(1, 1),
            "A survives alone in its row/col on the second fit — a stale bit \
             from the first fit must not linger"
        );
    }

    /// A scrubbed frame or a replayed snapshot can arrive at a different `bins`
    /// than the last one drawn, so the grid buffers must be sized to the state
    /// they came from (an undersized one reads as a blank heatmap, not a panic)
    /// and the products must hold that state's entries only.
    #[test]
    fn a_refit_at_a_different_bins_sizes_the_grid_and_shows_only_the_second_fit() {
        // Fit 1, a 3x3 grid: 3 path indices, 3 curve points, 3 measurements.
        let mut small = diagonal_state(3);
        small.fit();
        let rec = FitRecording::from_state(&small);
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

        // Fit 2, a 10x10 grid.
        let mut big = diagonal_state(10);
        big.fit();
        let rec = FitRecording::from_state(&big);

        assert_eq!(rec.geom().bins, 10, "the geometry follows the refit");
        // Row 9 exists only in a 10x10 buffer: at 3x3 the flat index 9*10+9
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
