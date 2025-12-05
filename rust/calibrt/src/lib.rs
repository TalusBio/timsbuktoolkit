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
    fn new(mut points: Vec<Point>) -> Result<Self, CalibRtError> {
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
    let optimal_path_points = pathfinding::find_optimal_path(&mut filtered_nodes);
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

    calibrate_with_ranges(points, x_range, y_range, grid_size)
}
