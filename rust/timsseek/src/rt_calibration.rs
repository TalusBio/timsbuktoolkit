use crate::ScorerQueriable;
use crate::scoring::pipeline::Scorer;
pub use calibrt::{
    CALIBRANT_WEIGHT,
    CalibRtError,
    CalibrationCurve as RTCalibration,
    CalibrationSnapshot,
    CalibrationState as CalibratedGrid,
    DEFAULT_RIDGE_FRACTION,
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

/// Multiplier applied to the ridge half-width to get the query tolerance.
/// 1.0 = use the FW@10%max directly (already generous).
/// Increase for more conservative searches.
const RIDGE_WIDTH_MULTIPLIER: f64 = 1.0;

/// Minimum RT tolerance in minutes (prevents pathologically tight windows).
const MIN_RT_TOLERANCE_MINUTES: f32 = 0.5;

/// The only version [`SavedCalibration`] reads or writes. Named once so the
/// writer and the reader's gate cannot disagree.
pub const CALIBRATION_FORMAT_VERSION: &str = "v3";

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
                rt_minutes: 0.5,
            },
        }
    }
}

/// A fitted RT calibration plus the tolerance windows measured alongside it.
///
/// Owns the [`CalibratedGrid`] rather than copies of its outputs, so the curve,
/// the ridge widths and the points it was fit on cannot drift from each other or
/// go stale. Taken by `&` on every scoring thread; nothing here refits.
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

    pub fn fit_points(&self) -> &[Point] {
        self.state.fit_points()
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

    /// Get per-query tolerance. Uses position-dependent ridge width when available,
    /// falls back to uniform `rt_tolerance_minutes` otherwise.
    /// `rt` is the library RT in seconds (pre-calibration).
    pub fn get_tolerance(&self, _mz: f64, _mobility: f32, rt: LibraryRT<f32>) -> Tolerance {
        let rt_tol_minutes = self
            .ridge_half_width_at(LibraryRT(rt.0 as f64))
            .map(|hw| (hw * RIDGE_WIDTH_MULTIPLIER / 60.0) as f32)
            .unwrap_or(self.rt_tolerance_minutes)
            .max(MIN_RT_TOLERANCE_MINUTES);

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

    /// Save calibration to JSON v2 format.
    ///
    /// The grid's own snapshot carries the points and their geometry, so
    /// [`Self::from_saved`] can refit an equal grid; the residual statistics and
    /// tolerances are read back verbatim.
    pub fn save_json(
        &self,
        rt_range_seconds: [f64; 2],
        n_scored: usize,
        path: &std::path::Path,
    ) -> Result<(), String> {
        let derivation = self.derivation.clone().unwrap_or_default();
        let saved = SavedCalibration {
            version: CALIBRATION_FORMAT_VERSION.to_string(),
            rt_range_seconds,
            calibration: self.state.snapshot(),
            errors: self.errors.clone(),
            derivation,
            tolerances: SavedTolerances {
                rt_minutes: self.rt_tolerance_minutes,
                mz_ppm: [self.mz_tolerance_ppm.0, self.mz_tolerance_ppm.1],
                mobility_pct: [
                    self.mobility_tolerance_pct.0 as f64,
                    self.mobility_tolerance_pct.1 as f64,
                ],
            },
            n_scored,
        };
        let json = serde_json::to_string_pretty(&saved).map_err(|e| e.to_string())?;
        std::fs::write(path, json).map_err(|e| e.to_string())
    }

    /// Rebuild a `CalibrationResult` from a saved file. Refitting the persisted
    /// points under their persisted geometry reproduces the grid, and with it the
    /// curve and the ridge widths; every other field is read back verbatim.
    pub fn from_saved(saved: &SavedCalibration) -> Result<Self, CalibRtError> {
        let state = CalibratedGrid::from_snapshot(&saved.calibration)?;
        Ok(Self::new(
            state,
            saved.tolerances.rt_minutes,
            (saved.tolerances.mz_ppm[0], saved.tolerances.mz_ppm[1]),
            (
                saved.tolerances.mobility_pct[0] as f32,
                saved.tolerances.mobility_pct[1] as f32,
            ),
        )?
        .with_error_stats(saved.errors.clone())
        .with_derivation(saved.derivation.clone()))
    }

    /// Read a saved calibration off disk and rebuild it, logging any
    /// provenance warning. See [`SavedCalibration::read`] and
    /// [`Self::from_saved`].
    pub fn read_json(
        path: &std::path::Path,
        raw_rt_range: Option<[f64; 2]>,
    ) -> Result<Self, String> {
        let (saved, warning) = SavedCalibration::read(path, raw_rt_range)?;
        if let Some(warning) = warning {
            tracing::warn!("{warning}");
        }
        Self::from_saved(&saved).map_err(|e| format!("{e:?}"))
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
    pub errors: DimensionErrors,
    pub derivation: DerivationParams,
    pub tolerances: SavedTolerances,
    pub n_scored: usize,
}

impl SavedCalibration {
    /// Parse a calibration file and check its provenance. The `Option<String>`
    /// is a human-readable reason to distrust the file, not an error: a
    /// calibration is only valid for the run it was fit on, and comparing RT
    /// spans is the one cheap way to catch the wrong file.
    ///
    /// `raw_rt_range` is the RT span of the file the calibration is about to be
    /// used on; `None` means there is nothing to check against, which is itself
    /// worth warning about.
    ///
    /// Every reader goes through here so the version gate and the provenance
    /// check cannot drift between the CLI and the viewer.
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

#[derive(Debug, Serialize, Deserialize)]
pub struct SavedTolerances {
    pub rt_minutes: f32,
    pub mz_ppm: [f64; 2],
    pub mobility_pct: [f64; 2],
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

    /// Save the fixture to a fresh temp dir, returning the dir so it outlives
    /// the path (dropping it removes the tree, including after a panic).
    fn save_fixture(calibration: &CalibrationResult) -> TempDir {
        let dir = tempfile::tempdir().unwrap();
        calibration
            .save_json([0.0, 1200.0], 999, &dir.path().join("calibration.json"))
            .unwrap();
        dir
    }

    #[test]
    fn save_json_then_from_saved_is_lossless() {
        let original = sample_calibration();
        let dir = save_fixture(&original);
        let path = dir.path().join("calibration.json");

        let loaded = CalibrationResult::read_json(&path, Some([0.0, 1200.0])).unwrap();

        assert_eq!(
            loaded.rt_tolerance_minutes(),
            original.rt_tolerance_minutes()
        );
        assert_eq!(loaded.mz_tolerance(), original.mz_tolerance());
        assert_eq!(loaded.mobility_tolerance(), original.mobility_tolerance());
        assert_eq!(loaded.fit_points(), original.fit_points());

        // The ridge is refit, not read back, so matching it is the claim that
        // the persisted points and geometry reproduce the grid they came from.
        let (loaded_ridge, orig_ridge) =
            (loaded.state().ridge_widths(), original.state().ridge_widths());
        assert_eq!(loaded_ridge.len(), orig_ridge.len(), "ridge width count");
        for (loaded_w, orig_w) in loaded_ridge.iter().zip(orig_ridge.iter()) {
            assert_eq!(loaded_w.library.0, orig_w.library.0);
            assert_eq!(loaded_w.half_width, orig_w.half_width);
            assert_eq!(loaded_w.ridge_weight, orig_w.ridge_weight);
        }

        // The tolerance getter folds the curve, the ridge widths and the floors
        // together, so matching it across the library RT range covers all three.
        for i in 0..=40 {
            let rt = LibraryRT(i as f32 * 10.0);
            assert_eq!(
                loaded.get_tolerance(500.0, 1.0, rt),
                original.get_tolerance(500.0, 1.0, rt),
                "tolerance mismatch at library RT {}",
                rt.0
            );
            assert_eq!(
                loaded.convert_irt(rt).0,
                original.convert_irt(rt).0,
                "RT conversion mismatch at library RT {}",
                rt.0
            );
        }
    }

    /// A file from an older writer must be refused, not silently reinterpreted:
    /// v2 stored ridge widths alongside the points, and reading one as v3 would
    /// take the stored widths for refit ones.
    #[test]
    fn read_refuses_a_foreign_version() {
        let dir = save_fixture(&sample_calibration());
        let path = dir.path().join("calibration.json");

        let text = std::fs::read_to_string(&path).unwrap();
        let mut value: serde_json::Value = serde_json::from_str(&text).unwrap();
        value.as_object_mut().unwrap()["version"] = serde_json::json!("v2");
        std::fs::write(&path, serde_json::to_string(&value).unwrap()).unwrap();

        let err = SavedCalibration::read(&path, Some([0.0, 1200.0])).unwrap_err();
        assert!(err.contains("v2") && err.contains("v3"), "{err}");
    }
}
