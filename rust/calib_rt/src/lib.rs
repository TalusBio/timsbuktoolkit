/// Core implementation of the Calib-RT algorithm.
/// Original description
/// https://doi.org/10.1093/bioinformatics/btae417

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
}

/// Represents a single data point on the library-measured-RT plane.
#[derive(Debug, Clone, Copy, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct Point {
    pub x: f64,
    pub y: f64,
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

    /// Predicts a calibrated measured RT (Y) for a given library RT (X).
    /// Returns an error if the value is outside the bounds of the calibration curve.
    pub fn predict(&self, x_val: f64) -> Result<f64, CalibRtError> {
        let first_x = self.points.first().unwrap().x;
        let last_x = self.points.last().unwrap().x;
        if x_val <= first_x {
            return Err(CalibRtError::OutOfBounds(self.predict_with_index(x_val, 1)));
        }

        if x_val >= last_x {
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

/// Represents a node (cell) in the grid.
#[derive(Debug, Clone, Copy, Default)]
pub struct Node {
    center: Point,
    frequency: u32,
    suppressed: bool,
}

/// Represents the gridded data after initial filtering.
pub struct Grid {
    nodes: Vec<Node>,
    x_range: (f64, f64),
    y_range: (f64, f64),
    x_span: f64,
    y_span: f64,
    bins: usize,
}

impl Grid {
    /// Creates a new, empty grid with a fixed geometry.
    /// The center of each node is constant based on the grid resolution.
    pub fn new(
        bins: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
    ) -> Result<Self, CalibRtError> {
        if bins == 0 {
            return Err(CalibRtError::ZeroRange);
        };
        let x_span = x_range.1 - x_range.0;
        let y_span = y_range.1 - y_range.0;

        if x_span <= 0.0 || y_span <= 0.0 {
            return Err(CalibRtError::ZeroRange);
        }

        let mut nodes = Vec::with_capacity(bins * bins);
        for r in 0..bins {
            for c in 0..bins {
                let center_x = x_range.0 + (c as f64 + 0.5) * (x_span / bins as f64);
                let center_y = y_range.0 + (r as f64 + 0.5) * (y_span / bins as f64);
                nodes.push(Node {
                    center: Point {
                        x: center_x,
                        y: center_y,
                    },
                    frequency: 0,
                    suppressed: false,
                });
            }
        }

        Ok(Self {
            nodes,
            x_range,
            y_range,
            x_span,
            y_span,
            bins,
        })
    }

    /// Adds a single point to the grid, incrementing the frequency of the corresponding cell.
    pub fn add_point(&mut self, point: &Point) {
        let Point { x, y } = point;

        let gx = (((x - self.x_range.0) / self.x_span) * self.bins as f64) as usize;
        let gy = (((y - self.y_range.0) / self.y_span) * self.bins as f64) as usize;

        let gx = gx.min(self.bins - 1);
        let gy = gy.min(self.bins - 1);

        let index = gy * self.bins + gx;
        if let Some(node) = self.nodes.get_mut(index) {
            node.frequency += 1;
        }
    }

    /// Applies nonmaximal suppression to the grid nodes.
    pub fn suppress_nonmax(&mut self) {
        // We start with 1s to prevent the max being empty in all ...
        let mut max_in_row = vec![1; self.bins];
        let mut max_in_col = vec![1; self.bins];

        for r in 0..self.bins {
            for c in 0..self.bins {
                let index = r * self.bins + c;
                let freq = self.nodes[index].frequency;
                if freq > max_in_row[r] {
                    max_in_row[r] = freq;
                }
                if freq > max_in_col[c] {
                    max_in_col[c] = freq;
                }
            }
        }

        for r in 0..self.bins {
            for c in 0..self.bins {
                let index = r * self.bins + c;
                let node = &mut self.nodes[index];
                node.suppressed = true;
                if node.frequency == max_in_row[r] && node.frequency == max_in_col[c] {
                    node.suppressed = false;
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------
// Module 2: Optimal Ascending Path Identification
// --------------------------------------------------------------------------------

/// Finds the highest-weight path through the nodes that satisfies the monotonic constraint.
fn find_optimal_path(nodes: &mut [Node]) -> Vec<Point> {
    if nodes.is_empty() {
        return Vec::new();
    }

    // Sort nodes primarily by x, then by y to process them in order for DAG pathfinding.
    nodes.sort_by(|a, b| {
        a.center
            .x
            .partial_cmp(&b.center.x)
            .unwrap()
            .then_with(|| a.center.y.partial_cmp(&b.center.y).unwrap())
    });

    let n = nodes.len();
    let mut max_weights = vec![0.0; n];
    let mut prev_node_indices = vec![None; n];

    for i in 0..n {
        max_weights[i] = nodes[i].frequency as f64; // Path can start at any node

        for j in 0..i {
            // 2.1 & 2.2: Check for monotonic edge and calculate weight
            if nodes[i].center.x > nodes[j].center.x && nodes[i].center.y > nodes[j].center.y {
                let dx = nodes[i].center.x - nodes[j].center.x;
                let dy = nodes[i].center.y - nodes[j].center.y;
                let dist = (dx * dx + dy * dy).sqrt();

                if dist > 1e-6 {
                    // Avoid division by zero
                    let edge_weight =
                        (nodes[i].frequency as f64 * nodes[j].frequency as f64) / dist;
                    let new_weight = max_weights[j] + edge_weight;

                    if new_weight > max_weights[i] {
                        max_weights[i] = new_weight;
                        prev_node_indices[i] = Some(j);
                    }
                }
            }
        }
    }

    // 2.3 Path Finding: Find the path with the maximum weight sum
    let mut max_path_weight = 0.0;
    let mut end_of_path_idx = 0;
    for i in 0..n {
        if max_weights[i] > max_path_weight {
            max_path_weight = max_weights[i];
            end_of_path_idx = i;
        }
    }

    // Reconstruct the path
    let mut path = Vec::new();
    let mut current_idx_opt = Some(end_of_path_idx);
    while let Some(current_idx) = current_idx_opt {
        path.push(nodes[current_idx].center);
        current_idx_opt = prev_node_indices[current_idx];
    }
    path.reverse();

    path
}

// --------------------------------------------------------------------------------
// Main Calib-RT Function
// --------------------------------------------------------------------------------

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
pub fn calibrate<'a>(
    points: impl IntoIterator<Item = &'a Point>,
    x_range: (f64, f64),
    y_range: (f64, f64),
    grid_size: usize,
) -> Result<CalibrationCurve, CalibRtError> {
    // Module 1: Grid data and apply nonmaximal suppression
    let mut grid = Grid::new(grid_size, x_range, y_range)?;

    let mut point_count = 0;
    for point in points {
        grid.add_point(point);
        point_count += 1;
    }

    if point_count == 0 {
        return Err(CalibRtError::NoPoints);
    }

    grid.suppress_nonmax();

    let mut filtered_nodes: Vec<Node> = grid
        .nodes
        .into_iter()
        .filter(|n| !n.suppressed && n.frequency > 0)
        .collect();

    // Module 2: Find the optimal ascending path
    let optimal_path_points = find_optimal_path(&mut filtered_nodes);

    // Module 3: Fit the final points and prepare for extrapolation
    CalibrationCurve::new(optimal_path_points)
}
