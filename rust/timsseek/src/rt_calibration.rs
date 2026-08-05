use crate::ScorerQueriable;
use crate::scoring::pipeline::Scorer;
pub use calibrt::{
    CALIBRANT_WEIGHT,
    CalibRtError,
    CalibrationCurve as RTCalibration,
    CalibrationSnapshot,
    CalibrationState as CalibratedGrid,
    FitEvent,
    FitObserver,
    LibraryRT,
    ObserveOpts,
    ObservedRTSeconds,
    Point,
    RidgeMeasurement,
    RidgeSummary,
};
use serde::{
    Deserialize,
    Serialize,
};
use timsquery::Tolerance;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};

/// The only version [`SavedCalibration`] reads or writes. Named once so the
/// writer and the reader's gate cannot disagree.
pub const CALIBRATION_FORMAT_VERSION: &str = "v3";

/// RT tolerance floor used when nothing configured one. Keeps a window from
/// closing to nothing where the ridge is narrow.
pub const DEFAULT_RT_FLOOR_MINUTES: f32 = 0.5;

/// The RT tolerance a ridge half-width implies, in minutes, floored at
/// `floor_minutes` — the same [`FloorsTriplet::rt_minutes`] that
/// [`DimensionErrors::derive_windows`] applies to the uniform fallback, so a
/// query gets one configured floor whichever path it takes.
///
/// The one place the ridge-to-tolerance rule is written. The search applies it
/// per query at the interpolated half-width; a UI showing a single number
/// applies it to [`RidgeSummary::weighted_half_width`]. Both must agree, or the
/// window a user is shown is not the window the search would open.
pub fn rt_tolerance_from_ridge(half_width_seconds: f64, floor_minutes: f32) -> f32 {
    ((half_width_seconds / 60.0) as f32).max(floor_minutes)
}

/// Per-dimension residual statistics. `1.4826 * mad` is a robust stdev
/// estimator (equals stdev for Gaussian data, resists outliers otherwise).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ErrorStats {
    pub mean: f32,
    pub stdev: f32,
    pub median: f32,
    pub mad: f32,
    pub n: usize,
}

impl ErrorStats {
    pub fn from_slice(xs: &[f32]) -> Self {
        if xs.is_empty() {
            return Self::default();
        }
        let n = xs.len();
        let mean = xs.iter().sum::<f32>() / n as f32;
        let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / n as f32;
        let stdev = var.sqrt();

        // Reuse the sorted buffer for the median, then the MAD computation.
        let mut buf: Vec<f32> = xs.to_vec();
        buf.sort_by(f32::total_cmp);
        let median = buf[n / 2];

        for v in buf.iter_mut() {
            *v = (*v - median).abs();
        }
        buf.sort_by(f32::total_cmp);
        let mad = buf[n / 2];

        Self {
            mean,
            stdev,
            median,
            mad,
            n,
        }
    }
}

/// Residual statistics per measurement dimension.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DimensionErrors {
    pub mz_ppm: ErrorStats,
    pub mobility_pct: ErrorStats,
    pub rt_seconds: ErrorStats,
}

impl DimensionErrors {
    /// The tolerance windows these residuals imply under `params`.
    ///
    /// Lives here rather than in whichever binary happens to measure the
    /// residuals: [`DerivationParams`] describes this derivation and is written
    /// into every saved calibration, so a second implementation of it would put
    /// two different meanings behind one recorded `method`.
    pub fn derive_windows(&self, params: &DerivationParams) -> DerivedWindows {
        let (mz_left, mz_right) =
            mad_symmetric_bounds(&self.mz_ppm, params.sigma.mz, params.floors.mz_ppm);
        DerivedWindows {
            // RT takes the half-width alone: these residuals are what is left
            // after the curve, so they carry no offset for the median to preserve.
            rt_minutes: (params.sigma.rt * 1.4826 * self.rt_seconds.mad / 60.0)
                .max(params.floors.rt_minutes),
            mz_ppm: (mz_left as f64, mz_right as f64),
            mobility_pct: mad_symmetric_bounds(
                &self.mobility_pct,
                params.sigma.mobility,
                params.floors.mobility_pct,
            ),
        }
    }
}

/// The tolerance windows [`DimensionErrors::derive_windows`] produces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DerivedWindows {
    pub rt_minutes: f32,
    pub mz_ppm: (f64, f64),
    pub mobility_pct: (f32, f32),
}

/// `median ± n_sigma * 1.4826 * MAD`, floored — the `"mad_symmetric"` method
/// [`DerivationParams`] names. Robust to tails: equal to `mean ± n_sigma * stdev`
/// for a Gaussian population, and it resists the outlier inflation that would
/// widen a window around a handful of bad matches.
fn mad_symmetric_bounds(stats: &ErrorStats, n_sigma: f32, min_val: f32) -> (f32, f32) {
    if stats.n == 0 {
        return (min_val, min_val);
    }
    let sigma = 1.4826 * stats.mad;
    let left = (-(stats.median - n_sigma * sigma)).max(min_val);
    let right = (stats.median + n_sigma * sigma).max(min_val);
    (left, right)
}

/// How tolerances are derived from `DimensionErrors`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerivationParams {
    /// Bound method, e.g. "mad_symmetric" = `median ± n_sigma * 1.4826 * MAD`.
    pub method: String,
    pub sigma: SigmaTriplet,
    pub floors: FloorsTriplet,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SigmaTriplet {
    pub mz: f32,
    pub mobility: f32,
    pub rt: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FloorsTriplet {
    pub mz_ppm: f32,
    pub mobility_pct: f32,
    pub rt_minutes: f32,
}

impl Default for DerivationParams {
    fn default() -> Self {
        Self {
            method: "mad_symmetric".to_string(),
            sigma: SigmaTriplet {
                mz: 1.5,
                mobility: 3.0,
                rt: 3.0,
            },
            floors: FloorsTriplet {
                mz_ppm: 0.1,
                mobility_pct: 0.1,
                rt_minutes: DEFAULT_RT_FLOOR_MINUTES,
            },
        }
    }
}

/// A fitted RT calibration plus the tolerance windows measured alongside it.
///
/// The grid is the source of truth for the curve, the ridge widths and the points
/// they came from. Taken by `&` on every scoring thread; nothing here refits.
pub struct CalibrationResult {
    state: CalibratedGrid,
    /// Fallback uniform RT tolerance, used where the grid measured no ridge.
    rt_tolerance_minutes: f32,
    mz_tolerance_ppm: (f64, f64),
    mobility_tolerance_pct: (f32, f32),
    errors: DimensionErrors,
    derivation: Option<DerivationParams>,
}

impl CalibrationResult {
    /// Wrap a fitted grid. `Err` when the grid has no curve: a result that
    /// cannot convert an RT is not one Phase 3 can use, and rejecting it here
    /// keeps [`Self::convert_irt`] infallible.
    pub fn new(
        state: CalibratedGrid,
        rt_tolerance_minutes: f32,
        mz_tolerance_ppm: (f64, f64),
        mobility_tolerance_pct: (f32, f32),
    ) -> Result<Self, CalibRtError> {
        if state.curve().is_none() {
            return Err(CalibRtError::NoPoints);
        }
        Ok(Self {
            state,
            rt_tolerance_minutes,
            mz_tolerance_ppm,
            mobility_tolerance_pct,
            errors: DimensionErrors::default(),
            derivation: None,
        })
    }

    /// The grid this was fit on — the source for the curve, the ridge widths and
    /// the calibrant points.
    pub fn state(&self) -> &CalibratedGrid {
        &self.state
    }

    pub fn with_error_stats(mut self, errors: DimensionErrors) -> Self {
        self.errors = errors;
        self
    }

    pub fn with_derivation(mut self, derivation: DerivationParams) -> Self {
        self.derivation = Some(derivation);
        self
    }

    pub fn errors(&self) -> &DimensionErrors {
        &self.errors
    }

    /// Interpolate ridge half-width at a given library RT (seconds).
    /// Returns the half-width in seconds, or None if no ridge data.
    fn ridge_half_width_at(&self, library_rt: LibraryRT<f64>) -> Option<f64> {
        ridge_half_width_interp(self.state.ridge_widths(), library_rt.0)
    }

    /// The fitted curve. Present for the lifetime of `self` — [`Self::new`]
    /// rejects a grid without one.
    fn curve(&self) -> &RTCalibration {
        self.state
            .curve()
            .expect("CalibrationResult::new rejects a grid with no curve")
    }

    /// Convert indexed RT to calibrated absolute RT (seconds).
    pub fn convert_irt(&self, irt: LibraryRT<f32>) -> ObservedRTSeconds<f32> {
        match self.curve().predict(LibraryRT(irt.0 as f64)) {
            Ok(rt) => ObservedRTSeconds(rt.0 as f32),
            Err(CalibRtError::OutOfBounds(rt)) => ObservedRTSeconds(rt as f32),
            Err(_) => ObservedRTSeconds(irt.0),
        }
    }

    /// Derived RT tolerance in minutes.
    pub fn rt_tolerance_minutes(&self) -> f32 {
        self.rt_tolerance_minutes
    }

    /// The configured RT floor. `rt_tolerance_minutes` was already floored at it
    /// by [`DimensionErrors::derive_windows`]; the ridge path has to apply it
    /// itself.
    fn rt_floor_minutes(&self) -> f32 {
        self.derivation
            .as_ref()
            .map_or(DEFAULT_RT_FLOOR_MINUTES, |d| d.floors.rt_minutes)
    }

    /// Get per-query tolerance. Uses position-dependent ridge width when available,
    /// falls back to uniform `rt_tolerance_minutes` otherwise.
    /// `rt` is the library RT in seconds (pre-calibration).
    pub fn get_tolerance(&self, _mz: f64, _mobility: f32, rt: LibraryRT<f32>) -> Tolerance {
        let rt_tol_minutes = match self.ridge_half_width_at(LibraryRT(rt.0 as f64)) {
            Some(half_width) => rt_tolerance_from_ridge(half_width, self.rt_floor_minutes()),
            None => self.rt_tolerance_minutes,
        };

        Tolerance {
            ms: MzTolerance::Ppm(self.mz_tolerance_ppm),
            rt: RtTolerance::Minutes((rt_tol_minutes, rt_tol_minutes)),
            mobility: MobilityTolerance::Pct(self.mobility_tolerance_pct),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }

    pub fn mz_tolerance(&self) -> (f64, f64) {
        self.mz_tolerance_ppm
    }

    pub fn mobility_tolerance(&self) -> (f32, f32) {
        self.mobility_tolerance_pct
    }

    /// Summary of ridge width measurements for reporting.
    pub fn ridge_width_summary(&self) -> Option<RidgeSummary> {
        RidgeSummary::of(self.state.ridge_widths())
    }

    /// Tolerance for the secondary spectral query at a detected apex.
    pub fn get_spectral_tolerance(&self) -> Tolerance {
        Tolerance {
            ms: MzTolerance::Ppm(self.mz_tolerance_ppm),
            rt: RtTolerance::Minutes((0.5 / 60.0, 0.5 / 60.0)), // ~0.5 seconds
            mobility: MobilityTolerance::Pct(self.mobility_tolerance_pct),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }

    /// Tight mobility tolerance for isotope pattern matching.
    pub fn get_isotope_tolerance(&self) -> Tolerance {
        self.get_spectral_tolerance()
            .with_mobility_tolerance(MobilityTolerance::Pct((3.0, 3.0)))
    }

    /// Write the calibration where a viewer or the dashboard can read it.
    ///
    /// The grid's snapshot is the record of the fit — a reader refits it to get
    /// the curve and the ridge widths back. Everything a search measured that the
    /// grid does not know goes in [`ResidualBlock`].
    pub fn save_json(
        &self,
        rt_range_seconds: [f64; 2],
        n_scored: usize,
        path: &std::path::Path,
    ) -> Result<(), String> {
        SavedCalibration::new(
            rt_range_seconds,
            self.state.snapshot(),
            self.rt_tolerance_minutes,
            Some(ResidualBlock {
                errors: self.errors.clone(),
                derivation: self.derivation.clone().unwrap_or_default(),
                mz_ppm: [self.mz_tolerance_ppm.0, self.mz_tolerance_ppm.1],
                mobility_pct: [
                    self.mobility_tolerance_pct.0 as f64,
                    self.mobility_tolerance_pct.1 as f64,
                ],
            }),
            n_scored,
        )
        .write(path)
    }

    /// Fallback when calibration fails: identity RT mapping, secondary tolerance.
    pub fn fallback<I: ScorerQueriable>(pipeline: &Scorer<I>) -> Self {
        let range = pipeline.index.ms1_cycle_mapping().range_milis();
        let start = range.0 as f64 / 1000.0;
        let end = range.1 as f64 / 1000.0;
        let mut state = CalibratedGrid::new(10, (start, end), (start, end), 10)
            .expect("the run's own RT range is a valid grid geometry");
        state
            .update(
                [start, end]
                    .into_iter()
                    .map(|rt| (LibraryRT(rt), ObservedRTSeconds(rt), CALIBRANT_WEIGHT)),
            )
            .expect("the identity endpoints are finite");
        state.fit();

        Self::new(state, 1.0, (10.0, 10.0), (5.0, 5.0))
            .expect("two identity endpoints always fit a two-point curve")
    }
}

/// JSON v3 calibration file format — shared between CLI and viewer.
///
/// `calibration` is the grid's own snapshot and the only record of the fit: the
/// curve and the ridge widths are recomputed by refitting it, so the file cannot
/// carry a curve that disagrees with the points that produced it.
#[derive(Debug, Serialize, Deserialize)]
pub struct SavedCalibration {
    pub version: String,
    pub rt_range_seconds: [f64; 2],
    pub calibration: CalibrationSnapshot,
    /// The uniform RT tolerance. Every writer has one — it is what a query falls
    /// back to where the grid measured no ridge.
    pub rt_tolerance_minutes: f32,
    /// What a search measured beyond the curve. `None` for a writer that measures
    /// no residuals — zeros there would read as "measured, and tight".
    #[serde(default)]
    pub residuals: Option<ResidualBlock>,
    pub n_scored: usize,
}

/// The residual statistics a search measures at its calibrant apexes, and the
/// m/z and mobility windows derived from them. One block because the windows are
/// meaningless without the statistics and the parameters that shaped them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidualBlock {
    pub errors: DimensionErrors,
    pub derivation: DerivationParams,
    pub mz_ppm: [f64; 2],
    pub mobility_pct: [f64; 2],
}

impl SavedCalibration {
    /// Assemble a file record. The version is stamped here, not by callers.
    pub fn new(
        rt_range_seconds: [f64; 2],
        calibration: CalibrationSnapshot,
        rt_tolerance_minutes: f32,
        residuals: Option<ResidualBlock>,
        n_scored: usize,
    ) -> Self {
        Self {
            version: CALIBRATION_FORMAT_VERSION.to_string(),
            rt_range_seconds,
            calibration,
            rt_tolerance_minutes,
            residuals,
            n_scored,
        }
    }

    /// Serialize to `path` in the layout [`Self::read`] expects.
    pub fn write(&self, path: &std::path::Path) -> Result<(), String> {
        let json = serde_json::to_string_pretty(self).map_err(|e| e.to_string())?;
        std::fs::write(path, json).map_err(|e| e.to_string())
    }

    /// Parse a calibration file and check its provenance. The `Option<String>`
    /// is a reason to distrust the file, not an error: a calibration is only
    /// valid for the run it was fit on, and `raw_rt_range` — the RT span of the
    /// run it is about to be used on — is the one cheap way to catch the wrong
    /// file. `None` there means nothing verifies it, which also warns.
    ///
    /// `calib_dash` is the deliberate exception to reading through here: it pulls
    /// `calibration` out with plain `serde_json` rather than depend on this crate.
    pub fn read(
        path: &std::path::Path,
        raw_rt_range: Option<[f64; 2]>,
    ) -> Result<(Self, Option<String>), String> {
        let json = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let saved: Self = serde_json::from_str(&json).map_err(|e| e.to_string())?;
        if saved.version != CALIBRATION_FORMAT_VERSION {
            return Err(format!(
                "Unsupported calibration version: {} (expected {CALIBRATION_FORMAT_VERSION})",
                saved.version
            ));
        }
        let warning = saved.provenance_warning(raw_rt_range);
        Ok((saved, warning))
    }

    /// Number of calibrant points the curve was fit on.
    pub fn n_calibrants(&self) -> usize {
        self.calibration.points.len()
    }

    fn provenance_warning(&self, raw_rt_range: Option<[f64; 2]>) -> Option<String> {
        let Some(raw) = raw_rt_range else {
            return Some(
                "No raw RT range to check the calibration against — nothing verifies it was \
                 fit on this run"
                    .to_string(),
            );
        };
        let overlap_lo = self.rt_range_seconds[0].max(raw[0]);
        let overlap_hi = self.rt_range_seconds[1].min(raw[1]);
        let overlap = (overlap_hi - overlap_lo).max(0.0);
        let span = self.rt_range_seconds[1] - self.rt_range_seconds[0];
        if span <= 0.0 || overlap / span >= 0.5 {
            return None;
        }
        Some(format!(
            "Calibration RT range [{:.1}, {:.1}]s overlaps the raw file's [{:.1}, {:.1}]s by \
             {:.0}% — it may have been fit on a different run",
            self.rt_range_seconds[0],
            self.rt_range_seconds[1],
            raw[0],
            raw[1],
            (overlap / span) * 100.0,
        ))
    }
}

/// Linearly interpolate the ridge half-width at a given library RT (seconds).
/// Endpoints clamp; empty input returns None.
pub fn ridge_half_width_interp(widths: &[RidgeMeasurement], library_rt_s: f64) -> Option<f64> {
    if widths.is_empty() {
        return None;
    }
    if library_rt_s <= widths[0].library.0 {
        return Some(widths[0].half_width);
    }
    if library_rt_s >= widths[widths.len() - 1].library.0 {
        return Some(widths[widths.len() - 1].half_width);
    }
    let pos = widths.partition_point(|m| m.library.0 < library_rt_s);
    if pos == 0 {
        return Some(widths[0].half_width);
    }
    let left = &widths[pos - 1];
    let right = &widths[pos];
    let t = (library_rt_s - left.library.0) / (right.library.0 - left.library.0).max(1e-9);
    Some(left.half_width + t * (right.half_width - left.half_width))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    /// A curve with enough spread that the grid fit is non-degenerate, plus ridge
    /// widths and tolerances distinct from every default, so a field silently
    /// dropped by the round-trip shows up as a mismatch.
    fn sample_calibration() -> CalibrationResult {
        let (grid_size, lookback) = (16, 4);
        let points: Vec<Point> = (0..40)
            .map(|i| {
                let library = i as f64 * 10.0;
                Point {
                    library,
                    observed: 30.0 + library * 1.7,
                    weight: CALIBRANT_WEIGHT,
                }
            })
            .collect();
        let mut state = CalibratedGrid::deferred(grid_size, lookback).unwrap();
        state
            .refit(
                grid_size,
                points.iter().map(|p| (p.library, p.observed)),
                &mut (),
                ObserveOpts::NONE,
            )
            .unwrap();
        assert!(
            !state.ridge_widths().is_empty(),
            "test fixture must produce ridge widths for the round-trip to be meaningful"
        );

        let errors = DimensionErrors {
            mz_ppm: ErrorStats::from_slice(&[1.0, 2.0, 3.0]),
            mobility_pct: ErrorStats::from_slice(&[0.5, 1.5]),
            rt_seconds: ErrorStats::from_slice(&[4.0, 6.0, 8.0]),
        };
        CalibrationResult::new(state, 1.25, (7.5, 8.5), (2.5, 3.5))
            .unwrap()
            .with_error_stats(errors)
            .with_derivation(DerivationParams::default())
    }

    fn save_fixture(calibration: &CalibrationResult) -> TempDir {
        let dir = tempfile::tempdir().unwrap();
        calibration
            .save_json([0.0, 1200.0], 999, &dir.path().join("calibration.json"))
            .unwrap();
        dir
    }

    /// The snapshot is the whole record of the fit, so refitting it must reproduce
    /// the curve and the ridge widths. Walks a reader's real path:
    /// `read` -> `from_snapshot` -> `new`.
    #[test]
    fn a_saved_snapshot_refits_to_the_same_calibration() {
        let original = sample_calibration();
        let dir = save_fixture(&original);
        let path = dir.path().join("calibration.json");

        let (saved, warning) = SavedCalibration::read(&path, Some([0.0, 1200.0])).unwrap();
        assert!(
            warning.is_none(),
            "the fixture's RT range is the one it claims: {warning:?}"
        );

        let residuals = saved
            .residuals
            .expect("a search's calibration carries its residuals");
        assert_eq!(residuals.mz_ppm, [7.5, 8.5]);
        assert_eq!(residuals.mobility_pct, [2.5, 3.5]);
        assert_eq!(
            residuals.errors.rt_seconds.n,
            original.errors().rt_seconds.n
        );

        let refit = CalibratedGrid::from_snapshot(&saved.calibration).unwrap();
        assert_eq!(refit.fit_points(), original.state().fit_points());

        let reloaded = CalibrationResult::new(
            refit,
            saved.rt_tolerance_minutes,
            (residuals.mz_ppm[0], residuals.mz_ppm[1]),
            (
                residuals.mobility_pct[0] as f32,
                residuals.mobility_pct[1] as f32,
            ),
        )
        .unwrap();

        // The ridge was refit, not read back, so matching it is the claim that the
        // persisted points and geometry reproduce the grid they came from.
        let (new_ridge, orig_ridge) = (
            reloaded.state().ridge_widths(),
            original.state().ridge_widths(),
        );
        assert_eq!(new_ridge.len(), orig_ridge.len(), "ridge width count");
        for (a, b) in new_ridge.iter().zip(orig_ridge.iter()) {
            assert_eq!(a.library.0, b.library.0);
            assert_eq!(a.half_width, b.half_width);
            assert_eq!(a.ridge_weight, b.ridge_weight);
        }

        // The tolerance getter folds the curve, the ridge widths and the floors
        // together, so matching it across the library RT range covers all three.
        for i in 0..=40 {
            let rt = LibraryRT(i as f32 * 10.0);
            assert_eq!(
                reloaded.get_tolerance(500.0, 1.0, rt),
                original.get_tolerance(500.0, 1.0, rt),
                "tolerance mismatch at library RT {}",
                rt.0
            );
            assert_eq!(
                reloaded.convert_irt(rt).0,
                original.convert_irt(rt).0,
                "RT conversion mismatch at library RT {}",
                rt.0
            );
        }
    }

    /// Pins the three per-dimension rules `derive_windows` folds together.
    #[test]
    fn derive_windows_applies_each_dimensions_rule() {
        let params = DerivationParams {
            method: "mad_symmetric".to_string(),
            sigma: SigmaTriplet {
                mz: 2.0,
                mobility: 2.0,
                rt: 3.0,
            },
            floors: FloorsTriplet {
                mz_ppm: 0.1,
                mobility_pct: 0.1,
                rt_minutes: 0.25,
            },
        };

        // Symmetric about zero: median 0, MAD 1.
        let symmetric = ErrorStats::from_slice(&[-2.0, -1.0, 0.0, 1.0, 2.0]);
        assert_eq!((symmetric.median, symmetric.mad), (0.0, 1.0));

        // Offset by 10: the median rides along, so the window is not centred on
        // zero — that offset is the systematic error the window has to cover.
        let offset = ErrorStats::from_slice(&[8.0, 9.0, 10.0, 11.0, 12.0]);
        assert_eq!((offset.median, offset.mad), (10.0, 1.0));

        let rt_seconds = ErrorStats::from_slice(&[-60.0, 0.0, 60.0]);
        assert_eq!(rt_seconds.mad, 60.0);

        let w = DimensionErrors {
            mz_ppm: offset,
            mobility_pct: symmetric,
            rt_seconds,
        }
        .derive_windows(&params);

        // Expectations are literals, not the implementation's own expression: at
        // MAD 1 and 2 sigma the half-spread is 2 * 1.4826 = 2.9652.
        //
        // The offset dimension's left bound lands under its floor — its median is
        // further from zero than the spread, so the window does not reach back.
        assert_eq!(w.mz_ppm.0, params.floors.mz_ppm as f64);
        assert!((w.mz_ppm.1 - 12.9652).abs() < 1e-4, "{:?}", w.mz_ppm);
        assert!(
            (w.mobility_pct.0 - 2.9652).abs() < 1e-4,
            "{:?}",
            w.mobility_pct
        );
        assert!(
            (w.mobility_pct.1 - 2.9652).abs() < 1e-4,
            "{:?}",
            w.mobility_pct
        );
        // RT ignores the median: 3 sigma over a 60s MAD, reported in minutes.
        assert!((w.rt_minutes - 4.4478).abs() < 1e-4, "{}", w.rt_minutes);

        // A dimension with nothing measured falls to its floor rather than zero,
        // which would be a window no query could match inside.
        let nothing = DimensionErrors::default().derive_windows(&params);
        let mz_floor = params.floors.mz_ppm as f64;
        assert_eq!(nothing.mz_ppm, (mz_floor, mz_floor));
        assert_eq!(
            nothing.mobility_pct,
            (params.floors.mobility_pct, params.floors.mobility_pct)
        );
        assert_eq!(nothing.rt_minutes, params.floors.rt_minutes);
    }

    /// A file with no `residuals` block, as a tool that measures none writes it.
    /// Spelled out rather than round-tripped through our own writer, so the test
    /// pins the layout a reader must accept.
    fn foreign_file(version: &str) -> String {
        format!(
            r#"{{
              "version": "{version}",
              "rt_range_seconds": [0.0, 1200.0],
              "calibration": {{
                "points": [[0.0, 30.0, 1.0], [100.0, 200.0, 1.0],
                           [200.0, 370.0, 1.0], [300.0, 540.0, 1.0]],
                "grid_size": 16,
                "lookback": 4
              }},
              "rt_tolerance_minutes": 1.25,
              "n_scored": 999
            }}"#
        )
    }

    fn write_foreign(version: &str) -> TempDir {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("calibration.json"), foreign_file(version)).unwrap();
        dir
    }

    #[test]
    fn read_accepts_a_calibration_with_no_residuals() {
        let dir = write_foreign(CALIBRATION_FORMAT_VERSION);
        let (saved, _) =
            SavedCalibration::read(&dir.path().join("calibration.json"), Some([0.0, 1200.0]))
                .unwrap();

        assert!(saved.residuals.is_none());
        assert_eq!(saved.n_calibrants(), 4);
        assert!(CalibratedGrid::from_snapshot(&saved.calibration).is_ok());
    }

    /// A foreign version must be refused rather than reinterpreted: the reader
    /// refits the snapshot, so a layout whose fields mean something else would be
    /// silently taken for this one.
    #[test]
    fn read_refuses_a_foreign_version() {
        let dir = write_foreign("v2");
        let err = SavedCalibration::read(&dir.path().join("calibration.json"), Some([0.0, 1200.0]))
            .unwrap_err();
        assert!(err.contains("v2") && err.contains("v3"), "{err}");
    }
}
