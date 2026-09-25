use crate::ScorerQueriable;
use crate::scoring::pipeline::Scorer;
pub use calibrt::{
    CALIBRANT_WEIGHT,
    CALIBRATION_FORMAT_VERSION,
    CalibRtError,
    CalibrationCurve as RTCalibration,
    CalibrationSnapshot,
    CalibrationState as CalibratedGrid,
    LibraryRT,
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

/// RT tolerance floor used when nothing configured one. Keeps a window from
/// closing to nothing where the ridge is narrow.
const DEFAULT_RT_FLOOR_MINUTES: f32 = 0.5;

/// The RT tolerance a ridge half-width implies, in minutes, floored at
/// `floor_minutes` -- the same [`FloorsTriplet::rt_minutes`] that
/// [`DimensionErrors::derive_windows`] applies to the uniform fallback. The
/// per-query search path and a single-number UI must both come through here, or
/// the window a user is shown is not the window the search would open.
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

/// `median ± n_sigma * 1.4826 * MAD`, floored -- the `"mad_symmetric"` method
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

/// Independently calibrated dimensions, with an optional RT fit.
///
/// The grid is the source of truth for the curve, the ridge widths and the points
/// they came from. Taken by `&` on every scoring thread; nothing here refits.
pub struct CalibrationResult {
    state: CalibratedGrid,
    /// Fallback uniform RT tolerance, used where the grid measured no ridge.
    rt_tolerance_minutes: f32,
    tolerance: Tolerance,
    errors: DimensionErrors,
    derivation: Option<DerivationParams>,
    library_rt_axis: timsquery::RtAxis,
}

impl CalibrationResult {
    /// Wrap a fitted grid. `Err` when the grid has no curve.
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
            tolerance: Tolerance {
                ms: MzTolerance::Ppm(mz_tolerance_ppm),
                rt: RtTolerance::Unrestricted,
                mobility: MobilityTolerance::Pct(mobility_tolerance_pct),
                quad: QuadTolerance::Absolute((0.1, 0.1)),
            },
            errors: DimensionErrors::default(),
            derivation: None,
            library_rt_axis: timsquery::RtAxis::Unspecified,
        })
    }

    /// Record the coordinate domain used when fitting this calibration.
    pub fn with_rt_axis(mut self, axis: &timsquery::RtAxis) -> Self {
        self.library_rt_axis = axis.clone();
        self
    }

    pub fn has_rt_calibration(&self) -> bool {
        self.state.curve().is_some()
    }

    /// Each axis with no measurements retains its configured extraction window.
    pub fn from_measurements(
        mut state: CalibratedGrid,
        configured: &Tolerance,
        errors: DimensionErrors,
        derivation: DerivationParams,
    ) -> Self {
        if state.curve().is_none() {
            state.reset();
        }
        let windows = errors.derive_windows(&derivation);
        let mut tolerance = configured
            .clone()
            .with_rt_tolerance(RtTolerance::Unrestricted);
        if errors.mz_ppm.n > 0 {
            let mut bounds = windows.mz_ppm;
            if let MzTolerance::Ppm(broad) = configured.ms {
                bounds = (bounds.0.min(broad.0), bounds.1.min(broad.1));
            }
            tolerance.ms = MzTolerance::Ppm(bounds);
        }
        if errors.mobility_pct.n > 0 {
            tolerance.mobility = MobilityTolerance::Pct(windows.mobility_pct);
        }
        Self {
            state,
            tolerance,
            rt_tolerance_minutes: windows.rt_minutes,
            errors,
            derivation: Some(derivation),
            library_rt_axis: timsquery::RtAxis::Unspecified,
        }
    }

    /// The grid this was fit on -- the source for the curve, the ridge widths and
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

    /// Convert indexed RT to calibrated absolute RT (seconds).
    pub fn convert_irt(&self, irt: LibraryRT<f32>) -> Result<ObservedRTSeconds<f32>, CalibRtError> {
        match self
            .state
            .curve()
            .ok_or(CalibRtError::NoPoints)?
            .predict(LibraryRT(irt.0 as f64))
        {
            Ok(rt) => Ok(ObservedRTSeconds(rt.0 as f32)),
            // calibrt reports its extrapolated observed value in this variant.
            Err(CalibRtError::OutOfBounds(seconds)) => Ok(ObservedRTSeconds(seconds as f32)),
            Err(err) => Err(err),
        }
    }

    /// Keep calibrated m/z and mobility windows without imposing an RT constraint.
    pub fn unrestricted_rt_tolerance(&self) -> Tolerance {
        self.tolerance.clone()
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
    /// `rt` is a coordinate on the original library axis.
    pub fn get_tolerance(&self, _mz: f64, _mobility: f32, rt: LibraryRT<f32>) -> Tolerance {
        let rt_tol_minutes = match self.ridge_half_width_at(LibraryRT(rt.0 as f64)) {
            Some(half_width) => rt_tolerance_from_ridge(half_width, self.rt_floor_minutes()),
            None => self.rt_tolerance_minutes,
        };

        if !self.has_rt_calibration() {
            return self.unrestricted_rt_tolerance();
        }
        self.tolerance
            .clone()
            .with_rt_tolerance(RtTolerance::Minutes((rt_tol_minutes, rt_tol_minutes)))
    }

    pub fn mz_tolerance(&self) -> &MzTolerance {
        &self.tolerance.ms
    }

    pub fn mobility_tolerance(&self) -> &MobilityTolerance {
        &self.tolerance.mobility
    }

    /// Summary of ridge width measurements for reporting.
    pub fn ridge_width_summary(&self) -> Option<RidgeSummary> {
        RidgeSummary::of(self.state.ridge_widths())
    }

    /// Tolerance for the secondary spectral query at a detected apex.
    pub fn get_spectral_tolerance(&self) -> Tolerance {
        self.tolerance
            .clone()
            .with_rt_tolerance(RtTolerance::Minutes((0.5 / 60.0, 0.5 / 60.0)))
    }

    /// Tight mobility refinement only when mobility calibration succeeded.
    pub fn get_isotope_tolerance(&self) -> Tolerance {
        let tolerance = self.get_spectral_tolerance();
        if self.errors.mobility_pct.n > 0 {
            tolerance.with_mobility_tolerance(MobilityTolerance::Pct((3.0, 3.0)))
        } else {
            tolerance
        }
    }

    /// Write the calibration where a viewer or the dashboard can read it.
    ///
    /// The grid's snapshot is the record of the fit -- a reader refits it to get
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
                tolerance: self.unrestricted_rt_tolerance(),
            }),
            n_scored,
        )
        .with_library_rt_axis(self.library_rt_axis.clone())
        .write(path)
    }

    /// No measured dimensions: preserve configured windows and search all RT.
    pub fn fallback<I: ScorerQueriable>(pipeline: &Scorer<I>) -> Self {
        Self::from_measurements(
            CalibratedGrid::deferred(10, 10).expect("valid default grid geometry"),
            &pipeline.broad_tolerance,
            DimensionErrors::default(),
            DerivationParams::default(),
        )
    }
}

/// The calibration file format, carrying the residual block a search measures.
pub type SavedCalibration = calibrt::SavedCalibration<ResidualBlock>;

/// The residual statistics a search measures at its calibrant apexes, and the
/// m/z and mobility windows derived from them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidualBlock {
    pub errors: DimensionErrors,
    pub derivation: DerivationParams,
    pub tolerance: Tolerance,
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

    #[test]
    fn dimensions_calibrate_independently_and_missing_axes_keep_configuration() {
        let configured = Tolerance {
            ms: MzTolerance::Absolute((0.01, 0.02)),
            mobility: MobilityTolerance::Unrestricted,
            rt: RtTolerance::Minutes((1.0, 2.0)),
            quad: QuadTolerance::Absolute((0.3, 0.4)),
        };
        for (mz, mob) in [(false, false), (true, false), (false, true), (true, true)] {
            let errors = DimensionErrors {
                mz_ppm: ErrorStats::from_slice(if mz { &[2.0, 3.0, 4.0] } else { &[] }),
                mobility_pct: ErrorStats::from_slice(if mob { &[0.2, 0.3, 0.4] } else { &[] }),
                ..Default::default()
            };
            let result = CalibrationResult::from_measurements(
                CalibratedGrid::deferred(10, 3).unwrap(),
                &configured,
                errors,
                DerivationParams::default(),
            );
            assert!(!result.has_rt_calibration());
            assert!(result.convert_irt(LibraryRT(20.0)).is_err());
            let primary = result.get_tolerance(500.0, 1.0, LibraryRT(20.0));
            assert_eq!(primary.rt, RtTolerance::Unrestricted);
            assert_eq!(primary.quad, configured.quad);
            if mz {
                assert!(matches!(primary.ms, MzTolerance::Ppm(_)));
            } else {
                assert_eq!(primary.ms, configured.ms);
            }
            if mob {
                assert!(matches!(primary.mobility, MobilityTolerance::Pct(_)));
            } else {
                assert_eq!(primary.mobility, configured.mobility);
            }
            let secondary = result.get_spectral_tolerance();
            assert_eq!(secondary.ms, primary.ms);
            assert_eq!(secondary.mobility, primary.mobility);
            assert_eq!(secondary.quad, primary.quad);
            if !mob {
                assert_eq!(result.get_isotope_tolerance().mobility, configured.mobility);
            }
            let dir = TempDir::new().unwrap();
            let path = dir.path().join("calibration.json");
            result.save_json([0.0, 600.0], 10, &path).unwrap();
            let (saved, _) = SavedCalibration::read(&path, None).unwrap();
            assert!(saved.calibration.points.is_empty());
            let residuals = saved.residuals.unwrap();
            assert_eq!(residuals.tolerance.ms, primary.ms);
            assert_eq!(residuals.tolerance.mobility, primary.mobility);
            assert_eq!(residuals.errors.mz_ppm.n, if mz { 3 } else { 0 });
        }
    }

    #[test]
    fn measured_mz_window_stays_inside_configured_ppm_window() {
        let configured = Tolerance {
            ms: MzTolerance::Ppm((8.0, 15.0)),
            ..Tolerance::default()
        };
        let result = CalibrationResult::from_measurements(
            CalibratedGrid::deferred(10, 3).unwrap(),
            &configured,
            DimensionErrors {
                mz_ppm: ErrorStats::from_slice(&[-60.0, 0.0, 60.0]),
                ..Default::default()
            },
            DerivationParams::default(),
        );
        assert_eq!(*result.mz_tolerance(), configured.ms);
    }

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
            .refit(grid_size, points.iter().map(|p| (p.library, p.observed)))
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
    fn unavailable_calibration_never_reinterprets_the_input_as_seconds() {
        let mut calibration = sample_calibration();
        calibration.state = CalibratedGrid::deferred(10, 3).unwrap();
        assert!(matches!(
            calibration.convert_irt(LibraryRT(-20.0)),
            Err(CalibRtError::NoPoints)
        ));
        assert!(matches!(
            calibration.convert_irt(LibraryRT(0.0)),
            Err(CalibRtError::NoPoints)
        ));
    }

    #[test]
    fn a_saved_snapshot_refits_to_the_same_calibration() {
        let axis = timsquery::RtAxis::NormalizedIndex {
            scale: Some("test-reference".into()),
        };
        let original = sample_calibration().with_rt_axis(&axis);
        let dir = save_fixture(&original);
        let path = dir.path().join("calibration.json");

        let (saved, warning) = SavedCalibration::read(&path, Some([0.0, 1200.0])).unwrap();
        assert!(
            warning.is_none(),
            "the fixture's RT range is the one it claims: {warning:?}"
        );

        assert_eq!(saved.library_rt_axis, axis);

        let residuals = saved
            .residuals
            .expect("a search's calibration carries its residuals");
        assert_eq!(residuals.tolerance.ms, MzTolerance::Ppm((7.5, 8.5)));
        assert_eq!(
            residuals.tolerance.mobility,
            MobilityTolerance::Pct((2.5, 3.5))
        );
        assert_eq!(
            residuals.errors.rt_seconds.n,
            original.errors().rt_seconds.n
        );

        let refit = CalibratedGrid::from_snapshot(&saved.calibration).unwrap();
        assert_eq!(refit.fit_points(), original.state().fit_points());

        let reloaded =
            CalibrationResult::new(refit, saved.rt_tolerance_minutes, (7.5, 8.5), (2.5, 3.5))
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
                reloaded
                    .convert_irt(rt)
                    .map(|v| v.0)
                    .map_err(|e| format!("{e:?}")),
                original
                    .convert_irt(rt)
                    .map(|v| v.0)
                    .map_err(|e| format!("{e:?}")),
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
        // zero -- that offset is the systematic error the window has to cover.
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
        // The offset dimension's left bound lands under its floor -- its median is
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
              "library_rt_axis": {{"kind": "seconds"}},
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
        assert!(
            err.contains("v2") && err.contains(CALIBRATION_FORMAT_VERSION),
            "{err}"
        );
    }
}
