//! Core implementation of the Calib-RT algorithm.
//!
//! Original description: <https://doi.org/10.1093/bioinformatics/btae417>

pub mod grid;
mod pathfinding;
pub mod plotting;
pub use grid::Grid;
use tracing::{
    info,
    warn,
};

/// Minimum denominator for slope calculations to avoid division by zero.
const MIN_SLOPE_DENOMINATOR: f64 = 1e-9;

/// Default width for calibration curve plots.
const CALIBRATION_PLOT_WIDTH: usize = 40;

/// Default height for calibration curve plots.
const CALIBRATION_PLOT_HEIGHT: usize = 20;

/// Custom error types for the Calib-RT library.
#[derive(Debug, Clone)]
pub enum CalibRtError {
    /// Returned when calibration is attempted with no input points.
    NoPoints,
    /// Returned when calibration is attempted with not enough points to interpolate.
    InsufficientPoints,
    /// Returned when the grid is created with a zero-width or zero-height range.
    ZeroRange,
    /// Returned when prediction is attempted for a value outside the calibrated range.
    OutOfBounds(f64),
    /// Returned when the weight of a point is invalid (e.g., Nan, infinite ...).
    UnsupportedWeight(f64),
}

/// Represents a single data point on the library-measured-RT plane.
#[derive(Debug, Clone, Copy, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub weight: f64,
}

/// Represents the final calibration curve.
/// It holds the sorted points from the optimal path and can be used
/// to predict calibrated RTs.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CalibrationCurve {
    points: Vec<Point>,
    slopes: Vec<f64>,
}

impl CalibrationCurve {
    /// Creates a new CalibrationCurve from a slice of points.
    /// Precomputes slopes for faster prediction.
    pub(crate) fn new(mut points: Vec<Point>) -> Result<Self, CalibRtError> {
        if points.is_empty() {
            return Err(CalibRtError::NoPoints);
        }
        if points.len() < 2 {
            return Err(CalibRtError::InsufficientPoints);
        }

        points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());

        let slopes = points
            .windows(2)
            .map(|p| (p[1].y - p[0].y) / (p[1].x - p[0].x).max(MIN_SLOPE_DENOMINATOR))
            .collect();

        Ok(Self { points, slopes })
    }

    /// Read access to the sorted calibration points.
    pub fn points(&self) -> &[Point] {
        &self.points
    }

    pub fn wrmse<'a>(&self, test_points: impl Iterator<Item = &'a Point> + 'a) -> f64 {
        let mut total_error = 0.0;
        let mut weight: f64 = 0.0;

        for p in test_points {
            match self.predict(p.x) {
                Ok(predicted_y) => {
                    let error = predicted_y - p.y;
                    total_error += (error * error) * p.weight;
                    weight += p.weight;
                }
                Err(_) => {
                    // Ignore out-of-bounds points for MSE calculation
                }
            }
        }

        if weight == 0.0 {
            f64::NAN // No valid predictions
        } else {
            (total_error / weight).sqrt()
        }
    }

    /// Predicts a calibrated measured RT (Y) for a given library RT (X).
    /// Returns an error if the value is outside the bounds of the calibration curve.
    pub fn predict(&self, x_val: f64) -> Result<f64, CalibRtError> {
        let first_x = self.points.first().unwrap().x;
        let last_x = self.points.last().unwrap().x;
        if x_val < first_x {
            return Err(CalibRtError::OutOfBounds(self.predict_with_index(x_val, 1)));
        }

        if x_val > last_x {
            return Err(CalibRtError::OutOfBounds(
                self.predict_with_index(x_val, self.slopes.len()),
            ));
        }

        // Find the partition point; first element >= x_val.
        let i = self.points.partition_point(|p| p.x < x_val);
        // Clamp to [1, slopes.len()] — partition_point can return 0 when x_val == first_x
        let i = i.max(1).min(self.slopes.len());
        Ok(self.predict_with_index(x_val, i))
    }

    /// Internal prediction function that performs linear interpolation using a precomputed slope.
    ///
    /// # Arguments
    /// * `x_val` - The x-coordinate to predict y for
    /// * `i` - The partition index from the sorted points array (must satisfy: 1 <= i <= points.len())
    ///
    /// # Panics
    /// Panics if `i == 0` or `i > slopes.len()`, as these violate the interpolation invariants.
    ///
    /// # Implementation Note
    /// Uses slope between points[i-1] and points[i] to interpolate.
    /// When called from `predict()`, `i` is guaranteed valid via partition_point().
    /// When called for out-of-bounds extrapolation, caller must ensure valid index.
    fn predict_with_index(&self, x_val: f64, i: usize) -> f64 {
        assert!(
            i > 0 && i <= self.slopes.len(),
            "Index {} out of valid range [1, {}] for interpolation",
            i,
            self.slopes.len()
        );
        let p1 = self.points[i - 1];
        let slope = self.slopes[i - 1];
        p1.y + (x_val - p1.x) * slope
    }
}

/// Serializable snapshot of calibration data — points + config.
/// Used for save/load. Does not include the fitted curve (reconstructed on load).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CalibrationSnapshot {
    pub points: Vec<[f64; 3]>,  // [x, y, weight]
    pub grid_size: usize,
    pub lookback: usize,
}

/// Reusable calibration state for incremental fitting. Owns all allocations.
pub struct CalibrationState {
    grid: grid::Grid,
    path_indices: Vec<usize>,
    dp_max_weights: Vec<f64>,
    dp_prev_indices: Vec<Option<usize>>,
    curve: Option<CalibrationCurve>,
    stale: bool,
    lookback: usize,
}

impl CalibrationState {
    pub fn new(
        grid_size: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
        lookback: usize,
    ) -> Result<Self, CalibRtError> {
        Ok(Self {
            grid: grid::Grid::new(grid_size, x_range, y_range)?,
            path_indices: Vec::new(),
            dp_max_weights: Vec::new(),
            dp_prev_indices: Vec::new(),
            curve: None,
            stale: false,
            lookback,
        })
    }

    pub fn update(&mut self, points: impl Iterator<Item = (f64, f64, f64)>) {
        for (x, y, w) in points {
            let _ = self.grid.add_point(&Point { x, y, weight: w });
        }
        self.stale = true;
    }

    pub fn fit(&mut self) {
        if self.grid.suppress_nonmax().is_err() {
            self.curve = None;
            self.path_indices.clear();
            self.stale = false;
            return;
        }

        // Collect non-suppressed nodes for pathfinding
        let mut filtered: Vec<grid::Node> = self.grid.grid_cells()
            .iter()
            .filter(|n| !n.suppressed && n.center.weight > 0.0)
            .copied()
            .collect();

        // Pathfinding with reused buffers
        let path_points = pathfinding::find_optimal_path(
            &mut filtered,
            self.lookback,
            &mut self.dp_max_weights,
            &mut self.dp_prev_indices,
        );

        // Store path indices by matching path points back to grid cells
        self.path_indices.clear();
        for pp in &path_points {
            if let Some(idx) = self.grid.grid_cells().iter().position(|n| {
                (n.center.x - pp.x).abs() < 1e-9 && (n.center.y - pp.y).abs() < 1e-9
            }) {
                self.path_indices.push(idx);
            }
        }

        self.curve = CalibrationCurve::new(path_points).ok();
        self.stale = false;
    }

    pub fn reset(&mut self) {
        self.grid.reset();
        self.curve = None;
        self.path_indices.clear();
        self.stale = false;
    }

    pub fn grid_cells(&self) -> &[grid::Node] {
        self.grid.grid_cells()
    }

    pub fn grid_bins(&self) -> usize {
        self.grid.bins
    }

    pub fn grid_x_range(&self) -> (f64, f64) {
        self.grid.x_range
    }

    pub fn grid_y_range(&self) -> (f64, f64) {
        self.grid.y_range
    }

    pub fn path_indices(&self) -> &[usize] {
        &self.path_indices
    }

    pub fn curve(&self) -> Option<&CalibrationCurve> {
        self.curve.as_ref()
    }

    /// Bundle current config into a snapshot (caller provides the points).
    pub fn save_snapshot(&self, points: &[(f64, f64, f64)]) -> CalibrationSnapshot {
        CalibrationSnapshot {
            points: points.iter().map(|&(x, y, w)| [x, y, w]).collect(),
            grid_size: self.grid.bins,
            lookback: self.lookback,
        }
    }

    /// Reconstruct a CalibrationState from a snapshot.
    pub fn from_snapshot(snapshot: &CalibrationSnapshot) -> Result<Self, CalibRtError> {
        if snapshot.points.is_empty() {
            return Err(CalibRtError::NoPoints);
        }
        let x_range = compute_range(snapshot.points.iter().map(|p| p[0]))?;
        let y_range = compute_range(snapshot.points.iter().map(|p| p[1]))?;

        let mut state = Self::new(snapshot.grid_size, x_range, y_range, snapshot.lookback)?;
        state.update(snapshot.points.iter().map(|p| (p[0], p[1], p[2])));
        state.fit();
        Ok(state)
    }

    pub fn is_stale(&self) -> bool {
        self.stale
    }
}

#[cfg(test)]
mod calibration_state_tests {
    use super::*;

    #[test]
    fn test_update_fit_cycle() {
        let mut state = CalibrationState::new(10, (0.0, 100.0), (0.0, 100.0), 30).unwrap();
        let points: Vec<(f64, f64, f64)> = (0..10)
            .map(|i| {
                let v = (i as f64) * 10.0 + 5.0;
                (v, v, 1.0)
            })
            .collect();

        state.update(points.into_iter());
        assert!(state.is_stale());

        state.fit();
        assert!(!state.is_stale());
        assert!(state.curve().is_some());

        let curve = state.curve().unwrap();
        let pred = curve.predict(50.0).unwrap();
        assert!((pred - 50.0).abs() < 5.0, "predicted {} expected ~50.0", pred);
    }

    #[test]
    fn test_reset_clears_state() {
        let mut state = CalibrationState::new(10, (0.0, 100.0), (0.0, 100.0), 30).unwrap();
        let points = vec![(25.0, 25.0, 1.0), (75.0, 75.0, 1.0)];
        state.update(points.into_iter());
        state.fit();
        assert!(state.curve().is_some());

        state.reset();
        assert!(state.curve().is_none());
        assert!(state.path_indices().is_empty());
        assert!(!state.is_stale());
    }

    #[test]
    fn test_refit_after_reset_update() {
        let mut state = CalibrationState::new(10, (0.0, 100.0), (0.0, 100.0), 30).unwrap();

        // First fit: y = x
        let points1: Vec<_> = (0..10).map(|i| ((i as f64) * 10.0 + 5.0, (i as f64) * 10.0 + 5.0, 1.0)).collect();
        state.update(points1.into_iter());
        state.fit();
        let curve1_pred = state.curve().unwrap().predict(50.0).unwrap();

        // Reset and refit: y = 2x
        state.reset();
        let points2: Vec<_> = (0..10).map(|i| ((i as f64) * 10.0 + 5.0, (i as f64) * 20.0 + 5.0, 1.0)).collect();
        state.update(points2.into_iter());
        state.fit();
        let curve2_pred = state.curve().unwrap().predict(50.0).unwrap();

        assert!((curve2_pred - curve1_pred).abs() > 10.0);
    }
}

/// Computes the min and max values from an iterator of f64 values.
///
/// # Returns
/// - `Ok((min, max))` if at least one valid value exists
/// - `Err(CalibRtError::NoPoints)` if no valid values exist
fn compute_range(values: impl Iterator<Item = f64>) -> Result<(f64, f64), CalibRtError> {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut count = 0;

    for val in values {
        if val.is_finite() {
            min = min.min(val);
            max = max.max(val);
            count += 1;
        }
    }

    if count == 0 || !min.is_finite() || !max.is_finite() {
        return Err(CalibRtError::NoPoints);
    }

    Ok((min, max))
}

/// Calibrates retention times using the Calib-RT algorithm with explicit ranges.
///
/// This is the lower-level API that requires you to specify the x and y ranges explicitly.
/// Consider using [`calibrate`] for automatic range detection.
///
/// # Arguments
/// * `points` - A slice of `Point` structs representing the data.
/// * `x_range` - The min and max values for the X dimension.
/// * `y_range` - The min and max values for the Y dimension.
/// * `grid_size` - The size of the grid for initial filtering (e.g., 100).
///
/// # Returns
/// A `Result` containing a `CalibrationCurve` or a `CalibRtError`.
pub fn calibrate_with_ranges(
    points: &[Point],
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid_size: usize,
    lookback: usize,
) -> Result<CalibrationCurve, CalibRtError> {
    // Module 1: Grid data and apply nonmaximal suppression
    let mut grid = Grid::new(grid_size, x_range, y_range)?;

    grid.extend_points(points)?;
    grid.suppress_nonmax()?;
    grid.display_heatmap();

    let mut filtered_nodes: Vec<grid::Node> = grid
        .nodes
        .into_iter()
        .filter(|n| !n.suppressed && n.center.weight > 0.0)
        .collect();

    // Module 2: Find the optimal ascending path
    let mut max_weights = Vec::new();
    let mut prev_indices = Vec::new();
    let optimal_path_points = pathfinding::find_optimal_path(
        &mut filtered_nodes, lookback, &mut max_weights, &mut prev_indices,
    );
    // Module 3: Fit the final points and prepare for extrapolation
    let calcurve = CalibrationCurve::new(optimal_path_points);
    match &calcurve {
        Ok(c) => {
            let wrmse = c.wrmse(points.iter());
            info!("Calibration successful, WRMSE: {}", wrmse);
            plotting::plot_function(
                |x| {
                    c.predict(x).map_err(|e| match e {
                        CalibRtError::OutOfBounds(y) => y,
                        _ => panic!("Unexpected error during plotting"),
                    })
                },
                (x_range.0, x_range.1),
                CALIBRATION_PLOT_WIDTH,
                CALIBRATION_PLOT_HEIGHT,
            );
        }
        Err(e) => {
            warn!("Calibration failed: {:?}", e);
        }
    }

    calcurve
}

/// Calibrates retention times using the Calib-RT algorithm with automatic range detection.
///
/// This is a convenience wrapper that automatically computes the x and y ranges from the input points.
/// If you need explicit control over the ranges, use [`calibrate_with_ranges`] instead.
///
/// # Arguments
/// * `points` - A slice of `Point` structs representing the data.
/// * `grid_size` - The size of the grid for initial filtering (e.g., 100).
///
/// # Returns
/// A `Result` containing a `CalibrationCurve` or a `CalibRtError`.
///
/// # Example
/// ```
/// use calibrt::{Point, calibrate};
///
/// let points = vec![
///     Point { x: 1.0, y: 1.5, weight: 1.0 },
///     Point { x: 2.0, y: 2.5, weight: 1.0 },
///     Point { x: 3.0, y: 3.5, weight: 1.0 },
/// ];
///
/// let curve = calibrate(&points, 100).expect("Calibration failed");
/// ```
pub fn calibrate(points: &[Point], grid_size: usize) -> Result<CalibrationCurve, CalibRtError> {
    if points.is_empty() {
        return Err(CalibRtError::NoPoints);
    }

    let x_range = compute_range(points.iter().map(|p| p.x))?;
    let y_range = compute_range(points.iter().map(|p| p.y))?;

    calibrate_with_ranges(points, x_range, y_range, grid_size, 30)
}
