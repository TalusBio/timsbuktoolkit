//! Core implementation of the Calib-RT algorithm.
//!
//! Original description: <https://doi.org/10.1093/bioinformatics/btae417>

pub mod grid;
mod pathfinding;
pub mod types;
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

/// The grid weight one calibrant contributes. Every producer weighs calibrants
/// equally: weight decides which nodes survive `suppress_nonmax` and scales every
/// DP edge, so a per-calibrant weight would fit a curve no other consumer computes.
pub const CALIBRANT_WEIGHT: f64 = 1.0;

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

/// Everything a consumer reports about a fit's [`RidgeMeasurement`]s, folded in
/// the one place the arithmetic is written — the dashboard, the search's derived
/// tolerances and the CLI's log line all read this rather than each summing the
/// slice their own way.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RidgeSummary {
    /// Mean half-width with each column weighted by its `ridge_weight`, so a
    /// heavy column carries more authority than a lonely one. NaN when the
    /// total weight is not positive: a weighted mean over no weight is not a
    /// number, and 0.0 would read as "the ridge is infinitely tight".
    pub weighted_half_width: f64,
    pub min_half_width: f64,
    pub max_half_width: f64,
    pub n_columns: usize,
    /// Fraction of the columns' total weight that falls inside the ridge bounds
    /// (0.0–1.0). Higher = better agreement between library and raw file.
    pub in_ridge_ratio: f64,
}

impl RidgeSummary {
    /// `None` for an empty slice: there is no column count, minimum or maximum
    /// to report, and no fold over nothing produces one.
    ///
    /// The mean divides by the total ridge weight itself. A zero total divides
    /// 0.0 by 0.0 (every term carries the same weight) and lands on the NaN the
    /// field documents.
    pub fn of(widths: &[RidgeMeasurement]) -> Option<Self> {
        if widths.is_empty() {
            return None;
        }
        let sum = |f: fn(&RidgeMeasurement) -> f64| widths.iter().map(f).sum::<f64>();
        let ridge_weight = sum(|m| m.ridge_weight);
        let column_weight = sum(|m| m.column_weight);
        let hw = || widths.iter().map(|m| m.half_width);
        Some(Self {
            weighted_half_width: sum(|m| m.half_width * m.ridge_weight) / ridge_weight,
            min_half_width: hw().fold(f64::INFINITY, f64::min),
            max_half_width: hw().fold(f64::NEG_INFINITY, f64::max),
            n_columns: widths.len(),
            in_ridge_ratio: if column_weight > 0.0 {
                ridge_weight / column_weight
            } else {
                0.0
            },
        })
    }
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
/// `Grid::add_point`'s `gy * bins + gx`.
pub enum FitEvent<'a> {
    FitStarted {
        geom: GridGeom,
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
        /// `path`'s cells as row-major grid indices, one per point, from the
        /// grid's own [`GridGeom`] arithmetic — so an overlay never has to
        /// re-derive them and risk landing in a different cell.
        indices: &'a [usize],
        /// The DP-chosen segment within `path`: `path[..dp_range.start]` is a
        /// greedily attached prefix and `path[dp_range.end..]` a greedily
        /// attached suffix (Pass 2's monotonic extension beyond what the DP
        /// itself scored).
        dp_range: std::ops::Range<usize>,
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

/// The no-op observer.
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

/// Reusable calibration state for incremental fitting. Keeps the grid and the
/// pathfinding buffers across fits; the path itself is allocated per fit.
pub struct CalibrationState {
    grid: grid::Grid,
    path_indices: Vec<usize>,
    /// Survivors of suppression, refilled per fit. Sized at `bins`, which is a
    /// hint and not a bound: `suppress_nonmax` keeps every node tied for a
    /// row/column max, so ties can legitimately exceed `bins`.
    filtered: Vec<grid::Node>,
    scratch: pathfinding::PathfindingScratch,
    curve: Option<CalibrationCurve>,
    /// The fit's ridge widths at [`DEFAULT_RIDGE_FRACTION`], measured by
    /// [`Self::fit_with`] so that reading them takes `&self`. A consumer that
    /// wants another fraction calls [`Self::measure_ridge_width`] instead; that
    /// path leaves this alone.
    ridge_widths: Vec<RidgeMeasurement>,
    /// The points fed through [`Self::update`] since the grid was last cleared.
    /// Kept because the grid bins them on the way in: a fitted grid cannot say
    /// what it was fit on, and that is what [`Self::snapshot`] has to persist.
    fit_points: Vec<Point>,
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
            ridge_widths: Vec::new(),
            fit_points: Vec::new(),
            lookback,
        })
    }

    /// A state whose geometry is not known yet. The unit range is a placeholder
    /// that the first [`Self::refit`] replaces with the points' own extents;
    /// fitting before then fits an empty grid.
    pub fn deferred(grid_size: usize, lookback: usize) -> Result<Self, CalibRtError> {
        Self::new(grid_size, (0.0, 1.0), (0.0, 1.0), lookback)
    }

    /// The whole re-fit sequence, in the one place every consumer shares:
    /// derive the grid geometry from the points, `reconfigure` onto it, `update`,
    /// then `fit_with`. Returns the [`GridRanges`] the fit actually ran on.
    ///
    /// The geometry is derived *here*, by [`point_ranges`], rather than passed in:
    /// a caller that supplies the acquisition RT range instead would clamp an
    /// iRT-scaled library — whose RTs fall entirely outside it — into one edge
    /// column. Every calibrant weighs [`CALIBRANT_WEIGHT`].
    ///
    /// On `Err` the previous fit is left alone: a later, larger point set may well
    /// span a usable range.
    pub fn refit<O: FitObserver>(
        &mut self,
        bins: usize,
        points: impl Iterator<Item = (f64, f64)> + Clone,
        obs: &mut O,
        opts: ObserveOpts,
    ) -> Result<GridRanges, CalibRtError> {
        let (x_range, y_range) = point_ranges(points.clone())?;
        self.reconfigure(bins, x_range, y_range)?;
        self.update(
            points.map(|(lib, obs)| (LibraryRT(lib), ObservedRTSeconds(obs), CALIBRANT_WEIGHT)),
        )?;
        self.fit_with(obs, opts);
        Ok((x_range, y_range))
    }

    /// Feed points into the grid. Returns an error if any point has
    /// non-finite coordinates or weight (NaN/Inf), indicating a bug upstream.
    pub fn update(
        &mut self,
        points: impl Iterator<Item = (LibraryRT<f64>, ObservedRTSeconds<f64>, f64)>,
    ) -> Result<(), CalibRtError> {
        for (lib_rt, obs_rt, w) in points {
            let point = Point {
                library: lib_rt.0,
                observed: obs_rt.0,
                weight: w,
            };
            self.grid.add_point(&point)?;
            self.fit_points.push(point);
        }
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

        let (path_points, dp_range) = pathfinding::find_optimal_path(
            &mut self.filtered,
            self.lookback,
            &mut self.scratch,
            obs,
            opts,
        );

        self.path_indices.clear();
        self.path_indices.extend(
            path_points
                .iter()
                .map(|p| self.grid.cell_of(p.library, p.observed)),
        );

        obs.on_event(FitEvent::PathFound {
            path: &path_points,
            indices: &self.path_indices,
            dp_range,
        });

        self.curve = CalibrationCurve::new(path_points).ok();

        if let Some(c) = self.curve.as_ref() {
            obs.on_event(FitEvent::CurveFit { curve: c });
        }

        // Measure the ridge here so consumers read it through `&self`; the
        // per-query RT tolerance needs it on every scoring thread. Deliberately
        // silent: `FitEvent::RidgeMeasured` stays the explicit measurement's
        // event, so a consumer that also measures at its own fraction records
        // one ridge per fit rather than two.
        let widths = self.measure_ridge_width_core(DEFAULT_RIDGE_FRACTION);
        self.ridge_widths = widths;
    }

    /// Drop the previous fit's results, keeping the grid and the buffers.
    fn clear_fit(&mut self) {
        self.curve = None;
        self.path_indices.clear();
        self.ridge_widths.clear();
    }

    pub fn reset(&mut self) {
        self.grid.reset();
        self.fit_points.clear();
        self.clear_fit();
    }

    /// Re-point `self` at a new geometry and clear the previous fit — see
    /// `grid::Grid::reconfigure` for what stays allocated.
    pub fn reconfigure(
        &mut self,
        bins: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
    ) -> Result<(), CalibRtError> {
        self.grid.reconfigure(bins, x_range, y_range)?;
        self.fit_points.clear();
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

    /// The points the grid currently holds, in the order [`Self::update`] took
    /// them.
    pub fn fit_points(&self) -> &[Point] {
        &self.fit_points
    }

    /// Everything needed to rebuild an equal state: the points, and the grid
    /// geometry they were binned under. Refitting a snapshot under its own
    /// `(grid_size, lookback)` reproduces the curve, so this is the state's
    /// persistent form — see [`Self::from_snapshot`].
    pub fn snapshot(&self) -> CalibrationSnapshot {
        CalibrationSnapshot {
            points: self
                .fit_points
                .iter()
                .map(|p| [p.library, p.observed, p.weight])
                .collect(),
            grid_size: self.grid_bins(),
            lookback: self.lookback,
        }
    }

    /// The current fit's ridge widths at [`DEFAULT_RIDGE_FRACTION`], sorted by
    /// library RT. Empty until a fit succeeds, and again once one fails, which is
    /// the "no ridge data" case consumers fall back from.
    pub fn ridge_widths(&self) -> &[RidgeMeasurement] {
        &self.ridge_widths
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
        let widths = self.measure_ridge_width_core(fraction);
        obs.on_event(FitEvent::RidgeMeasured { widths: &widths });
        widths
    }

    /// The measurement itself. Reads the nodes into the blur buffers and walks
    /// out from each path cell, so it is idempotent: the accumulated weights it
    /// derives from are never written, and a later `update` can still add points
    /// to the same grid.
    ///
    /// `&mut self` is for the two derived blur buffers, not for the fit.
    fn measure_ridge_width_core(&mut self, fraction: f64) -> Vec<RidgeMeasurement> {
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

        // Sorted by library RT because every consumer interpolates between
        // adjacent columns, which is only meaningful in x order. The DP path is
        // not guaranteed to arrive that way.
        widths.sort_by(|a, b| {
            a.library
                .0
                .partial_cmp(&b.library.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        widths
    }

    /// Reconstruct a CalibrationState from a snapshot.
    pub fn from_snapshot(snapshot: &CalibrationSnapshot) -> Result<Self, CalibRtError> {
        if snapshot.points.is_empty() {
            return Err(CalibRtError::NoPoints);
        }
        let (x_range, y_range) = point_ranges(snapshot.points.iter().map(|p| (p[0], p[1])))?;
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
}

/// A grid geometry's extents: `(x_range, y_range)`, each a `(min, max)`, in the
/// order [`CalibrationState::new`] and [`CalibrationState::reconfigure`] take.
pub type GridRanges = ((f64, f64), (f64, f64));

/// The [`GridRanges`] a set of `(library, observed)` pairs spans.
///
/// A point with a non-finite coordinate contributes to neither axis (it would
/// poison both bounds); `Err(NoPoints)` when that leaves nothing. A range that
/// comes out empty or inverted on either axis is `Err(ZeroRange)` here rather
/// than later out of `Grid::new`, so a caller that only wants to know whether
/// a grid is configurable never has to build one.
pub fn point_ranges(
    points: impl IntoIterator<Item = (f64, f64)>,
) -> Result<GridRanges, CalibRtError> {
    let mut x = (f64::INFINITY, f64::NEG_INFINITY);
    let mut y = x;
    for (px, py) in points {
        if !px.is_finite() || !py.is_finite() {
            continue;
        }
        x = (x.0.min(px), x.1.max(px));
        y = (y.0.min(py), y.1.max(py));
    }
    if !x.0.is_finite() {
        return Err(CalibRtError::NoPoints);
    }
    if x.0 >= x.1 || y.0 >= y.1 {
        return Err(CalibRtError::ZeroRange);
    }
    Ok((x, y))
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
    let (x_range, y_range) = point_ranges(points.iter().map(|p| (p.library, p.observed)))?;
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
    }

    impl FitObserver for Recorder {
        fn on_event(&mut self, ev: FitEvent<'_>) {
            match ev {
                FitEvent::FitStarted { geom, .. } => {
                    self.names.push("start");
                    self.geom = Some(geom);
                }
                FitEvent::Suppressed { .. } => self.names.push("suppressed"),
                FitEvent::DpNode { i, chose, .. } => {
                    self.names.push("dp");
                    self.dp_edges.push((i, chose));
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

    /// The grid bins its input, so a fitted grid cannot report what it was fit
    /// on. The state keeps the points alongside it, and must clear them exactly
    /// when the grid's accumulation is cleared — otherwise a snapshot persists
    /// points the grid no longer holds.
    #[test]
    fn fit_points_track_the_grids_accumulation() {
        let mut s = diagonal_state();
        assert_eq!(s.fit_points().len(), 10);

        // A second batch accumulates, matching `update`'s effect on the grid.
        s.update(std::iter::once((LibraryRT(2.5), ObservedRTSeconds(2.5), 1.0)))
            .unwrap();
        assert_eq!(s.fit_points().len(), 11);

        // A fit reads the points; it must not consume them.
        s.fit();
        assert!(s.curve().is_some(), "the diagonal must fit");
        assert_eq!(s.fit_points().len(), 11, "fitting is not consuming");

        let snap = s.snapshot();
        assert_eq!(snap.points.len(), 11);
        assert_eq!((snap.grid_size, snap.lookback), (10, 5));
        assert_eq!(
            snap.points[0],
            [0.5, 0.5, 1.0],
            "weights survive into the snapshot"
        );

        // Refitting the snapshot reproduces the curve it was taken from.
        let rebuilt = CalibrationState::from_snapshot(&snap).unwrap();
        let (orig, new) = (s.curve().unwrap(), rebuilt.curve().unwrap());
        assert_eq!(orig.points().len(), new.points().len());
        for x in [1.0, 3.3, 7.7] {
            assert_eq!(
                orig.predict(LibraryRT(x)).unwrap().0,
                new.predict(LibraryRT(x)).unwrap().0,
                "prediction diverged at {x}"
            );
        }

        s.reconfigure(10, (0.0, 20.0), (0.0, 20.0)).unwrap();
        assert!(
            s.fit_points().is_empty(),
            "reconfigure clears the grid, so it must clear the points"
        );

        let mut s = diagonal_state();
        s.reset();
        assert!(s.fit_points().is_empty(), "reset likewise");
    }

    #[test]
    fn events_arrive_in_pipeline_order() {
        let mut s = diagonal_state();
        let mut rec = Recorder::default();
        s.fit_with(&mut rec, ObserveOpts::NONE);
        assert_eq!(
            rec.names,
            vec!["start", "suppressed", "path", "curve"],
            "no dp events when dp_nodes is off, and no ridge event from the fit"
        );
        let g = rec.geom.expect("FitStarted must be emitted");
        assert_eq!(g.bins, 10);
        assert_eq!(g.x_range, (0.0, 10.0));
        assert_eq!(g.y_range, (0.0, 10.0));

        // And the ridge call is where it does come from.
        s.measure_ridge_width_with(0.1, &mut rec);
        assert_eq!(rec.names.last(), Some(&"ridge"));
    }

    #[test]
    fn dp_events_appear_once_per_node_when_enabled() {
        let mut s = diagonal_state();
        let mut rec = Recorder::default();
        s.fit_with(&mut rec, ObserveOpts { dp_nodes: true });
        assert!(rec.names.contains(&"dp"), "dp events must be emitted");
        assert_eq!(rec.dp_edges.len(), 10, "one event per DP node");
    }
}

#[cfg(test)]
mod ridge_summary_tests {
    use super::*;

    /// A column carrying twice its ridge weight in total, so `in_ridge_ratio`
    /// is 0.5 for any mix of these.
    fn m(half_width: f64, ridge_weight: f64) -> RidgeMeasurement {
        RidgeMeasurement {
            library: LibraryRT(1.0),
            half_width,
            ridge_weight,
            column_weight: ridge_weight * 2.0,
        }
    }

    #[test]
    fn the_mean_is_weighted_by_ridge_weight_and_the_spread_is_not() {
        let s = RidgeSummary::of(&[m(10.0, 1.0), m(20.0, 3.0)]).unwrap();
        // (10*1 + 20*3) / 4 = 17.5, not the unweighted 15.0.
        assert!((s.weighted_half_width - 17.5).abs() < 1e-9);
        assert_eq!((s.min_half_width, s.max_half_width), (10.0, 20.0));
        assert_eq!(s.n_columns, 2);
        assert!((s.in_ridge_ratio - 0.5).abs() < 1e-9);
        assert!(RidgeSummary::of(&[]).is_none(), "nothing to count or bound");
    }

    /// A total ridge weight of 0.25 must report the half-width it measured, not
    /// a fraction of it. And a weightless column keeps its count and bounds —
    /// only the mean goes NaN, where a 0.0 would read as a perfectly tight
    /// ridge.
    #[test]
    fn a_total_weight_below_one_does_not_shrink_the_mean() {
        let s = RidgeSummary::of(&[m(30.0, 0.25)]).unwrap();
        assert!(
            (s.weighted_half_width - 30.0).abs() < 1e-9,
            "got {}",
            s.weighted_half_width
        );
        let zero = RidgeSummary::of(&[m(5.0, 0.0)]).unwrap();
        assert!(zero.weighted_half_width.is_nan());
        assert_eq!((zero.n_columns, zero.min_half_width), (1, 5.0));
    }

    /// The two axes are bounded independently, a non-finite point drops out of
    /// both, and the three refusals are distinguishable: nothing left to bound
    /// is `NoPoints`, while a single collapsed axis is `ZeroRange` — which is
    /// what a grid built from these ranges would have said anyway.
    #[test]
    fn point_ranges_bounds_both_axes_and_names_each_refusal() {
        let pts = [(2.0, 30.0), (1.0, 10.0), (f64::NAN, 99.0), (3.0, 20.0)];
        assert_eq!(
            point_ranges(pts).unwrap(),
            ((1.0, 3.0), (10.0, 30.0)),
            "the NaN point must not reach the observed axis either"
        );
        for (pts, want) in [
            (vec![], "NoPoints"),
            (vec![(f64::INFINITY, 1.0)], "NoPoints"),
            (vec![(1.0, 1.0), (1.0, 2.0)], "ZeroRange"),
            (vec![(1.0, 1.0), (2.0, 1.0)], "ZeroRange"),
        ] {
            let got = format!("{:?}", point_ranges(pts).unwrap_err());
            assert_eq!(got, want);
        }
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

        state.fit();
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
    }
}
