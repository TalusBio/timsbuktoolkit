//! Core implementation of the Calib-RT algorithm.
//!
//! Original description: <https://doi.org/10.1093/bioinformatics/btae417>

pub mod grid;
mod pathfinding;
pub mod plotting;
pub mod types;
pub use grid::Grid;
use tracing::{
    info,
    warn,
};
pub use types::{
    LibraryRT,
    ObservedRTSeconds,
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
    pub library: f64,
    pub observed: f64,
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

        points.sort_by(|a, b| a.library.partial_cmp(&b.library).unwrap());

        let slopes = points
            .windows(2)
            .map(|p| {
                (p[1].observed - p[0].observed)
                    / (p[1].library - p[0].library).max(MIN_SLOPE_DENOMINATOR)
            })
            .collect();

        Ok(Self { points, slopes })
    }

    /// Read access to the sorted calibration points.
    pub fn points(&self) -> &[Point] {
        &self.points
    }

    pub fn wrmse(
        &self,
        test_points: impl Iterator<Item = (LibraryRT<f64>, ObservedRTSeconds<f64>, f64)>,
    ) -> f64 {
        let mut total_error = 0.0;
        let mut weight: f64 = 0.0;

        for (lib_rt, obs_rt, w) in test_points {
            match self.predict(lib_rt) {
                Ok(predicted_y) => {
                    let error = predicted_y.0 - obs_rt.0;
                    total_error += (error * error) * w;
                    weight += w;
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
    pub fn predict(&self, lib_rt: LibraryRT<f64>) -> Result<ObservedRTSeconds<f64>, CalibRtError> {
        let x_val = lib_rt.0;
        let first_x = self.points.first().unwrap().library;
        let last_x = self.points.last().unwrap().library;
        if x_val < first_x {
            return Err(CalibRtError::OutOfBounds(self.predict_with_index(x_val, 1)));
        }

        if x_val > last_x {
            return Err(CalibRtError::OutOfBounds(
                self.predict_with_index(x_val, self.slopes.len()),
            ));
        }

        // Find the partition point; first element >= x_val.
        let i = self.points.partition_point(|p| p.library < x_val);
        // Clamp to [1, slopes.len()] — partition_point can return 0 when x_val == first_x
        let i = i.max(1).min(self.slopes.len());
        Ok(ObservedRTSeconds(self.predict_with_index(x_val, i)))
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
        p1.observed + (x_val - p1.library) * slope
    }
}

/// Default `fraction` for [`CalibrationState::measure_ridge_width`]: expand
/// away from the path cell until the weight drops below 10% of it (FW@10%max).
/// Lives here because calibrt owns the measurement — every consumer that wants
/// "the" ridge width must read this rather than spell `0.1` again.
pub const DEFAULT_RIDGE_FRACTION: f64 = 0.1;

/// Measurement of the evidence ridge width at one grid column.
#[derive(Debug, Clone)]
pub struct RidgeMeasurement {
    /// Center library RT position.
    pub library: LibraryRT<f64>,
    /// Half-width of the ridge in y-units (seconds).
    pub half_width: f64,
    /// Weight inside the ridge bounds.
    pub ridge_weight: f64,
    /// Total weight in the full column (all rows at this x).
    pub column_weight: f64,
}

/// Serializable snapshot of calibration data — points + config.
/// Used for save/load. Does not include the fitted curve (reconstructed on load).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CalibrationSnapshot {
    pub points: Vec<[f64; 3]>, // [library, observed, weight]
    pub grid_size: usize,
    pub lookback: usize,
}

/// Grid geometry, emitted once per fit so a consumer can lay out `cells`
/// without inferring it from the node coordinates (which is impossible for an
/// empty or single-occupied grid).
#[derive(Debug, Clone, Copy)]
pub struct GridGeom {
    pub bins: usize,
    pub x_range: (f64, f64),
    pub y_range: (f64, f64),
}

/// One step of the fit, borrowed from the state that produced it.
///
/// `cells` is `bins * bins` ROW-MAJOR: `index = row * bins + col`, where `row`
/// indexes the observed-RT axis and `col` the library-RT axis. This matches
/// `Grid::add_point`'s `gy * bins + gx`. Nothing in the type system enforces
/// it and every consumer depends on it.
pub enum FitEvent<'a> {
    FitStarted {
        geom: GridGeom,
    },
    GridReady {
        cells: &'a [grid::Node],
    },
    Suppressed {
        cells: &'a [grid::Node],
    },
    /// Emitted once per DP node, only when `ObserveOpts::dp_nodes` is set.
    /// `considered` holds every `(predecessor_index, edge_weight)` the node
    /// evaluated, including the ones it rejected.
    DpNode {
        i: usize,
        node: &'a grid::Node,
        chose: Option<usize>,
        acc_weight: f64,
        considered: &'a [(usize, f64)],
    },
    PathFound {
        path: &'a [Point],
        /// The DP-chosen segment within `path`: `path[..dp_range.start]` is a
        /// greedily attached prefix and `path[dp_range.end..]` a greedily
        /// attached suffix (Pass 2's monotonic extension beyond what the DP
        /// itself scored).
        dp_range: std::ops::Range<usize>,
        /// The DP recurrence's objective value at the chosen end node —
        /// covers only `path[dp_range]`, not the greedily-attached prefix or
        /// suffix.
        dp_weight: f64,
    },
    CurveFit {
        curve: &'a CalibrationCurve,
    },
    RidgeMeasured {
        widths: &'a [RidgeMeasurement],
    },
}

pub trait FitObserver {
    fn on_event(&mut self, ev: FitEvent<'_>);
}

/// The no-op observer: `fit_with` is generic, so this monomorphizes away.
impl FitObserver for () {
    fn on_event(&mut self, _: FitEvent<'_>) {}
}

#[derive(Debug, Clone, Copy)]
pub struct ObserveOpts {
    /// Emit `DpNode` from the DP's inner loop. Off by default: it fires once
    /// per node.
    pub dp_nodes: bool,
}

impl ObserveOpts {
    pub const NONE: Self = Self { dp_nodes: false };
}

/// Reusable calibration state for incremental fitting. Owns all allocations.
pub struct CalibrationState {
    grid: grid::Grid,
    path_indices: Vec<usize>,
    /// Survivors of suppression, refilled per fit. Sized at `bins`, which is a
    /// hint and not a bound: `suppress_nonmax` keeps every node tied for a
    /// row/column max, so ties can legitimately exceed `bins`.
    filtered: Vec<grid::Node>,
    scratch: pathfinding::PathfindingScratch,
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
            filtered: Vec::with_capacity(grid_size),
            scratch: pathfinding::PathfindingScratch::default(),
            curve: None,
            stale: false,
            lookback,
        })
    }

    /// Feed points into the grid. Returns an error if any point has
    /// non-finite coordinates or weight (NaN/Inf), indicating a bug upstream.
    pub fn update(
        &mut self,
        points: impl Iterator<Item = (LibraryRT<f64>, ObservedRTSeconds<f64>, f64)>,
    ) -> Result<(), CalibRtError> {
        for (lib_rt, obs_rt, w) in points {
            self.grid.add_point(&Point {
                library: lib_rt.0,
                observed: obs_rt.0,
                weight: w,
            })?;
        }
        self.stale = true;
        Ok(())
    }

    pub fn fit(&mut self) {
        self.fit_with(&mut (), ObserveOpts::NONE)
    }

    pub fn fit_with<O: FitObserver>(&mut self, obs: &mut O, opts: ObserveOpts) {
        obs.on_event(FitEvent::FitStarted {
            geom: GridGeom {
                bins: self.grid.bins,
                x_range: self.grid.x_range,
                y_range: self.grid.y_range,
            },
        });
        obs.on_event(FitEvent::GridReady {
            cells: self.grid.grid_cells(),
        });

        let suppression_failed = self.grid.suppress_nonmax().is_err();
        obs.on_event(FitEvent::Suppressed {
            cells: self.grid.grid_cells(),
        });
        if suppression_failed {
            self.clear_fit();
            return;
        }

        self.filtered.clear();
        self.filtered.extend(
            self.grid
                .grid_cells()
                .iter()
                .filter(|n| !n.suppressed && n.center.weight > 0.0)
                .copied(),
        );

        let mut path_points = Vec::new();
        let (dp_range, dp_weight) = pathfinding::find_optimal_path(
            &mut self.filtered,
            self.lookback,
            &mut self.scratch,
            &mut path_points,
            obs,
            opts,
        );

        self.path_indices.clear();
        for pp in &path_points {
            if let Some(idx) = self.grid.grid_cells().iter().position(|n| {
                (n.center.library - pp.library).abs() < 1e-9
                    && (n.center.observed - pp.observed).abs() < 1e-9
            }) {
                self.path_indices.push(idx);
            }
        }

        obs.on_event(FitEvent::PathFound {
            path: &path_points,
            dp_range,
            dp_weight,
        });

        self.curve = CalibrationCurve::new(path_points).ok();
        self.stale = false;

        if let Some(c) = self.curve.as_ref() {
            obs.on_event(FitEvent::CurveFit { curve: c });
        }
    }

    /// Drop the previous fit's results, keeping the grid and the buffers.
    fn clear_fit(&mut self) {
        self.curve = None;
        self.path_indices.clear();
        self.stale = false;
    }

    pub fn reset(&mut self) {
        self.grid.reset();
        self.clear_fit();
    }

    /// Re-point `self` at a new geometry and clear the previous fit — see
    /// [`grid::Grid::reconfigure`] for what stays allocated.
    pub fn reconfigure(
        &mut self,
        bins: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
    ) -> Result<(), CalibRtError> {
        self.grid.reconfigure(bins, x_range, y_range)?;
        self.clear_fit();
        Ok(())
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

    /// Measure the width of the evidence "mountain" around the fitted path.
    ///
    /// For each grid column that contains a path cell, expands up and down
    /// from the path cell until cell weight drops below `fraction` of the
    /// path cell's weight. Returns `(column_center_x, half_width_y, total_weight)`.
    ///
    /// `fraction`: weight threshold as a fraction of the path cell's weight
    /// (e.g., 0.1 = expand until weight < 10% of path cell).
    /// `total_weight`: sum of all cell weights in the expanded range — heavier
    /// columns should carry more authority in tolerance estimation.
    pub fn measure_ridge_width(&mut self, fraction: f64) -> Vec<RidgeMeasurement> {
        self.measure_ridge_width_with(fraction, &mut ())
    }

    /// As [`Self::measure_ridge_width`], but reports the measurements through
    /// `obs` once they're computed.
    pub fn measure_ridge_width_with<O: FitObserver>(
        &mut self,
        fraction: f64,
        obs: &mut O,
    ) -> Vec<RidgeMeasurement> {
        let bins = self.grid.bins;
        let y_span = self.grid.y_range.1 - self.grid.y_range.0;
        let cell_h = y_span / bins as f64;

        // Sync node weights into buffer A, then blur into buffer B
        self.grid.sync_weights();
        self.grid.blur_weights();

        let mut widths = Vec::new();

        for &path_idx in &self.path_indices {
            let gx = path_idx % bins;
            let gy = path_idx / bins;
            let path_weight = self.grid.blurred_weight(gy, gx);
            if path_weight <= 0.0 {
                continue;
            }

            let threshold = path_weight * fraction;

            // Expand upward (increasing gy) from path cell
            let mut upper_gy = gy;
            let mut total_weight = path_weight;
            for dy in 1..bins {
                let check_gy = gy + dy;
                if check_gy >= bins {
                    break;
                }
                let w = self.grid.blurred_weight(check_gy, gx);
                if w < threshold {
                    break;
                }
                total_weight += w;
                upper_gy = check_gy;
            }

            // Expand downward (decreasing gy) from path cell
            let mut lower_gy = gy;
            for dy in 1..bins {
                if dy > gy {
                    break;
                }
                let check_gy = gy - dy;
                let w = self.grid.blurred_weight(check_gy, gx);
                if w < threshold {
                    break;
                }
                total_weight += w;
                lower_gy = check_gy;
            }

            // Sum all weights in this column for the in-ridge ratio
            let column_weight: f64 = (0..bins).map(|row| self.grid.blurred_weight(row, gx)).sum();

            let half_width = ((upper_gy - lower_gy) as f64 + 1.0) * cell_h * 0.5;

            widths.push(RidgeMeasurement {
                library: LibraryRT(self.grid.grid_cells()[path_idx].center.library),
                half_width,
                ridge_weight: total_weight,
                column_weight,
            });
        }

        obs.on_event(FitEvent::RidgeMeasured { widths: &widths });
        widths
    }

    /// Bundle current config into a snapshot (caller provides the points).
    pub fn save_snapshot(
        &self,
        points: &[(LibraryRT<f64>, ObservedRTSeconds<f64>, f64)],
    ) -> CalibrationSnapshot {
        CalibrationSnapshot {
            points: points
                .iter()
                .map(|&(lib, obs, w)| [lib.0, obs.0, w])
                .collect(),
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
        state.update(
            snapshot
                .points
                .iter()
                .map(|p| (LibraryRT(p[0]), ObservedRTSeconds(p[1]), p[2])),
        )?;
        state.fit();
        Ok(state)
    }

    pub fn is_stale(&self) -> bool {
        self.stale
    }

    #[cfg(test)]
    pub fn filtered_ptr(&self) -> *const grid::Node {
        self.filtered.as_ptr()
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
    let mut state = CalibrationState::new(grid_size, x_range, y_range, lookback)?;
    state.update(points.iter().map(|p| {
        (
            LibraryRT(p.library),
            ObservedRTSeconds(p.observed),
            p.weight,
        )
    }))?;

    state.fit();
    state.grid.display_heatmap();

    // No path at all — including the suppression short-circuit, which never
    // builds one — is `NoPoints`; a path too short to interpolate is not.
    let calcurve = state
        .curve()
        .cloned()
        .ok_or(if state.path_indices().is_empty() {
            CalibRtError::NoPoints
        } else {
            CalibRtError::InsufficientPoints
        });
    match &calcurve {
        Ok(c) => {
            let wrmse = c.wrmse(points.iter().map(|p| {
                (
                    LibraryRT(p.library),
                    ObservedRTSeconds(p.observed),
                    p.weight,
                )
            }));
            info!("Calibration successful, WRMSE: {}", wrmse);
            plotting::plot_function(
                |x| {
                    c.predict(LibraryRT(x))
                        .map(|obs| obs.0)
                        .map_err(|e| match e {
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
///     Point { library: 1.0, observed: 1.5, weight: 1.0 },
///     Point { library: 2.0, observed: 2.5, weight: 1.0 },
///     Point { library: 3.0, observed: 3.5, weight: 1.0 },
/// ];
///
/// let curve = calibrate(&points, 100).expect("Calibration failed");
/// ```
pub fn calibrate(points: &[Point], grid_size: usize) -> Result<CalibrationCurve, CalibRtError> {
    if points.is_empty() {
        return Err(CalibRtError::NoPoints);
    }

    let x_range = compute_range(points.iter().map(|p| p.library))?;
    let y_range = compute_range(points.iter().map(|p| p.observed))?;

    calibrate_with_ranges(points, x_range, y_range, grid_size, 30)
}

#[cfg(test)]
mod observer_tests {
    use super::*;

    /// Records event names plus the payloads the assertions need.
    #[derive(Default)]
    struct Recorder {
        names: Vec<&'static str>,
        geom: Option<GridGeom>,
        dp_edges: Vec<(usize, Option<usize>)>,
        /// `(library, observed)` of each DP node, indexed by `i`. `i` runs
        /// 0..n in order (one `DpNode` event per node), so `push`ing here
        /// keeps `dp_coords[i]` aligned with the node the DP loop saw at `i`.
        dp_coords: Vec<(f64, f64)>,
    }

    impl FitObserver for Recorder {
        fn on_event(&mut self, ev: FitEvent<'_>) {
            match ev {
                FitEvent::FitStarted { geom, .. } => {
                    self.names.push("start");
                    self.geom = Some(geom);
                }
                FitEvent::GridReady { .. } => self.names.push("grid"),
                FitEvent::Suppressed { .. } => self.names.push("suppressed"),
                FitEvent::DpNode { i, node, chose, .. } => {
                    self.names.push("dp");
                    self.dp_edges.push((i, chose));
                    self.dp_coords
                        .push((node.center.library, node.center.observed));
                }
                FitEvent::PathFound { .. } => self.names.push("path"),
                FitEvent::CurveFit { .. } => self.names.push("curve"),
                FitEvent::RidgeMeasured { .. } => self.names.push("ridge"),
            }
        }
    }

    /// A clean diagonal ridge: 10 points on y = x, one per grid column.
    fn diagonal_state() -> CalibrationState {
        let mut s = CalibrationState::new(10, (0.0, 10.0), (0.0, 10.0), 5).unwrap();
        let pts: Vec<_> = (0..10)
            .map(|i| {
                let v = i as f64 + 0.5;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0 + i as f64)
            })
            .collect();
        s.update(pts.into_iter()).unwrap();
        s
    }

    #[test]
    fn events_arrive_in_pipeline_order() {
        let mut s = diagonal_state();
        let mut rec = Recorder::default();
        s.fit_with(&mut rec, ObserveOpts::NONE);
        assert_eq!(
            rec.names,
            vec!["start", "grid", "suppressed", "path", "curve"],
            "no dp events when dp_nodes is off"
        );
        let g = rec.geom.expect("FitStarted must be emitted");
        assert_eq!(g.bins, 10);
        assert_eq!(g.x_range, (0.0, 10.0));
        assert_eq!(g.y_range, (0.0, 10.0));
    }

    #[test]
    fn dp_events_appear_only_when_enabled_and_respect_monotonicity() {
        let mut s = diagonal_state();
        let mut rec = Recorder::default();
        s.fit_with(&mut rec, ObserveOpts { dp_nodes: true });
        assert!(rec.names.contains(&"dp"), "dp events must be emitted");
        // The monotonic constraint the DP is supposed to enforce:
        // a chosen predecessor's library AND observed RT must both be
        // strictly less than the successor's.
        for (i, chose) in &rec.dp_edges {
            if let Some(j) = chose {
                let (lib_i, obs_i) = rec.dp_coords[*i];
                let (lib_j, obs_j) = rec.dp_coords[*j];
                assert!(
                    lib_i > lib_j && obs_i > obs_j,
                    "node {i} at ({lib_i}, {obs_i}) chose predecessor {j} at ({lib_j}, {obs_j}), \
                     which does not strictly precede it on both RT axes"
                );
            }
        }
        assert_eq!(rec.dp_edges.len(), 10, "one event per DP node");
    }

    #[test]
    fn ridge_events_come_from_the_ridge_call_not_the_fit() {
        let mut s = diagonal_state();
        let mut rec = Recorder::default();
        s.fit_with(&mut rec, ObserveOpts::NONE);
        assert!(!rec.names.contains(&"ridge"));
        s.measure_ridge_width_with(0.1, &mut rec);
        assert_eq!(rec.names.last(), Some(&"ridge"));
    }
}

#[cfg(test)]
mod calibration_state_tests {
    use super::*;

    #[test]
    fn test_update_fit_cycle() {
        let mut state = CalibrationState::new(10, (0.0, 100.0), (0.0, 100.0), 30).unwrap();
        let points: Vec<(LibraryRT<f64>, ObservedRTSeconds<f64>, f64)> = (0..10)
            .map(|i| {
                let v = (i as f64) * 10.0 + 5.0;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0)
            })
            .collect();

        state.update(points.into_iter()).unwrap();
        assert!(state.is_stale());

        state.fit();
        assert!(!state.is_stale());
        assert!(state.curve().is_some());

        let curve = state.curve().unwrap();
        let pred = curve.predict(LibraryRT(50.0)).unwrap();
        assert!(
            (pred.0 - 50.0).abs() < 5.0,
            "predicted {} expected ~50.0",
            pred.0
        );
    }

    #[test]
    fn test_reset_clears_state() {
        let mut state = CalibrationState::new(10, (0.0, 100.0), (0.0, 100.0), 30).unwrap();
        let points = vec![
            (LibraryRT(25.0), ObservedRTSeconds(25.0), 1.0),
            (LibraryRT(75.0), ObservedRTSeconds(75.0), 1.0),
        ];
        state.update(points.into_iter()).unwrap();
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
        let points1: Vec<_> = (0..10)
            .map(|i| {
                (
                    LibraryRT((i as f64) * 10.0 + 5.0),
                    ObservedRTSeconds((i as f64) * 10.0 + 5.0),
                    1.0,
                )
            })
            .collect();
        state.update(points1.into_iter()).unwrap();
        state.fit();
        let curve1_pred = state.curve().unwrap().predict(LibraryRT(50.0)).unwrap();

        // Reset and refit: y = 2x
        state.reset();
        let points2: Vec<_> = (0..10)
            .map(|i| {
                (
                    LibraryRT((i as f64) * 10.0 + 5.0),
                    ObservedRTSeconds((i as f64) * 20.0 + 5.0),
                    1.0,
                )
            })
            .collect();
        state.update(points2.into_iter()).unwrap();
        state.fit();
        let curve2_pred = state.curve().unwrap().predict(LibraryRT(50.0)).unwrap();

        assert!((curve2_pred.0 - curve1_pred.0).abs() > 10.0);
    }

    #[test]
    fn repeated_fits_do_not_grow_the_scratch_buffers() {
        let mut s = CalibrationState::new(20, (0.0, 20.0), (0.0, 20.0), 5).unwrap();
        let pts: Vec<_> = (0..20)
            .map(|i| {
                let v = i as f64 + 0.5;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0 + i as f64)
            })
            .collect();
        s.update(pts.iter().copied()).unwrap();
        s.fit();
        let filtered_ptr = s.filtered_ptr();

        for _ in 0..5 {
            s.reset();
            s.update(pts.iter().copied()).unwrap();
            s.fit();
        }

        assert_eq!(
            s.filtered_ptr(),
            filtered_ptr,
            "filtered must be the same allocation, not a same-sized new one"
        );
    }

    #[test]
    fn reconfigure_reuses_state_across_batches_with_changing_ranges() {
        let mut s = CalibrationState::new(10, (0.0, 10.0), (0.0, 10.0), 3).unwrap();
        let pts1: Vec<_> = (0..10)
            .map(|i| {
                let v = i as f64 + 0.5;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0 + i as f64)
            })
            .collect();
        s.update(pts1.iter().copied()).unwrap();
        s.fit();
        let filtered_ptr = s.filtered_ptr();

        // Same bins, a completely different (shifted, wider) range — the case
        // `reconfigure` exists to keep allocation-free.
        s.reconfigure(10, (100.0, 200.0), (100.0, 200.0)).unwrap();
        let pts2: Vec<_> = (0..10)
            .map(|i| {
                let v = 100.0 + i as f64 * 10.0 + 0.5;
                (LibraryRT(v), ObservedRTSeconds(v), 1.0 + i as f64)
            })
            .collect();
        s.update(pts2.iter().copied()).unwrap();
        s.fit();

        assert_eq!(s.grid_x_range(), (100.0, 200.0));
        assert_eq!(s.grid_y_range(), (100.0, 200.0));
        assert!(
            s.curve().unwrap().predict(LibraryRT(150.0)).is_ok(),
            "the new fit must be defined over the new range"
        );
        assert_eq!(
            s.filtered_ptr(),
            filtered_ptr,
            "reconfigure at unchanged bins must not reallocate `filtered`"
        );
    }
}
