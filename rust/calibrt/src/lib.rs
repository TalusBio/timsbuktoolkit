//! Core implementation of the Calib-RT algorithm.
//!
//! Original description
//! https://doi.org/10.1093/bioinformatics/btae417

pub mod grid;
mod pathfinding;
pub mod plotting;
pub use grid::Grid;
use tracing::info;

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
            .map(|p| (p[1].y - p[0].y) / (p[1].x - p[0].x).max(1e-9))
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

    /// Internal prediction function that assumes x_val is within bounds.
    fn predict_with_index(&self, x_val: f64, i: usize) -> f64 {
        // `i` must be > 0, because if `i` were 0, it would either be an exact
        // match on the first element (handled above) or x_val would be smaller
        // than the first element (caught by the bounds check).
        // Therefore, `i - 1` is always a valid index here.
        let p1 = self.points[i - 1];
        let slope = self.slopes[i - 1];
        p1.y + (x_val - p1.x) * slope
    }
}

/// Calibrates retention times using the Calib-RT algorithm.
///
/// # Arguments
/// * `points` - An iterator over `Point` structs representing the data.
/// * `x_range` - The min and max values for the X dimension.
/// * `y_range` - The min and max values for the Y dimension.
/// * `grid_size` - The size of the grid for initial filtering (e.g., 100).
///
/// # Returns
/// A `Result` containing a `CalibrationCurve` or a `CalibRtError`.
pub fn calibrate(
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
    match calcurve {
        Ok(ref c) => {
            let wrmse = c.wrmse(points.iter());
            info!("RMSE: {}", wrmse);
            plotting::plot_function(
                |x| {
                    c.predict(x).map_err(|e| match e {
                        CalibRtError::OutOfBounds(y) => y,
                        _ => panic!("Unexpected error during plotting"),
                    })
                },
                (x_range.0, x_range.1),
                40,
                20,
            );
        }
        Err(ref _e) => (),
    };

    calcurve
}
