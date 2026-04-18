use crate::ScorerQueriable;
use crate::scoring::pipeline::Scorer;
pub use calibrt::{
    CalibRtError,
    CalibrationCurve as RTCalibration,
    CalibrationSnapshot,
    CalibrationState as CalibratedGrid,
    LibraryRT,
    ObservedRTSeconds,
    Point,
    RidgeMeasurement,
    calibrate_with_ranges,
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

/// Immutable calibration result. Provides RT conversion and per-query tolerance
/// without mutating the speclib.
pub struct CalibrationResult {
    cal_curve: RTCalibration,
    /// Fallback uniform RT tolerance (used when no ridge data available).
    rt_tolerance_minutes: f32,
    mz_tolerance_ppm: (f64, f64),
    mobility_tolerance_pct: (f32, f32),
    /// Position-dependent ridge widths (sorted by x). Empty = use uniform fallback.
    ridge_widths: Vec<RidgeMeasurement>,
    errors: DimensionErrors,
    derivation: Option<DerivationParams>,
}

impl CalibrationResult {
    pub fn new(
        cal_curve: RTCalibration,
        rt_tolerance_minutes: f32,
        mz_tolerance_ppm: (f64, f64),
        mobility_tolerance_pct: (f32, f32),
    ) -> Self {
        Self {
            cal_curve,
            rt_tolerance_minutes,
            mz_tolerance_ppm,
            mobility_tolerance_pct,
            ridge_widths: Vec::new(),
            errors: DimensionErrors::default(),
            derivation: None,
        }
    }

    pub fn with_ridge_widths(mut self, mut widths: Vec<RidgeMeasurement>) -> Self {
        widths.sort_by(|a, b| {
            a.library
                .0
                .partial_cmp(&b.library.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.ridge_widths = widths;
        self
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
        ridge_half_width_interp(&self.ridge_widths, library_rt.0)
    }

    /// Convert indexed RT to calibrated absolute RT (seconds).
    pub fn convert_irt(&self, irt: LibraryRT<f32>) -> ObservedRTSeconds<f32> {
        match self.cal_curve.predict(LibraryRT(irt.0 as f64)) {
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
    pub fn ridge_width_summary(&self) -> Option<RidgeWidthSummary> {
        if self.ridge_widths.is_empty() {
            return None;
        }
        let total_weight: f64 = self.ridge_widths.iter().map(|m| m.ridge_weight).sum();
        let weighted_avg = self
            .ridge_widths
            .iter()
            .map(|m| m.half_width * m.ridge_weight)
            .sum::<f64>()
            / total_weight.max(1.0);
        let min = self
            .ridge_widths
            .iter()
            .map(|m| m.half_width)
            .fold(f64::MAX, f64::min);
        let max = self
            .ridge_widths
            .iter()
            .map(|m| m.half_width)
            .fold(0.0f64, f64::max);
        let total_column_weight: f64 = self.ridge_widths.iter().map(|m| m.column_weight).sum();
        let in_ridge_ratio = if total_column_weight > 0.0 {
            total_weight / total_column_weight
        } else {
            0.0
        };
        Some(RidgeWidthSummary {
            weighted_avg,
            min,
            max,
            n_columns: self.ridge_widths.len(),
            in_ridge_ratio,
        })
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
    pub fn save_json(
        &self,
        calibrant_points: &[(f64, f64, f64)],
        rt_range_seconds: [f64; 2],
        grid_size: usize,
        lookback: usize,
        n_scored: usize,
        path: &std::path::Path,
    ) -> Result<(), String> {
        let derivation = self.derivation.clone().unwrap_or_default();
        let saved = SavedCalibration {
            version: "v2".to_string(),
            rt_range_seconds,
            calibration: CalibrationSnapshot {
                points: calibrant_points
                    .iter()
                    .map(|&(x, y, w)| [x, y, w])
                    .collect(),
                grid_size,
                lookback,
            },
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
            n_calibrants: calibrant_points.len(),
            n_scored,
        };
        let json = serde_json::to_string_pretty(&saved).map_err(|e| e.to_string())?;
        std::fs::write(path, json).map_err(|e| e.to_string())
    }

    /// Load calibration from JSON v2 format.
    pub fn load_json(
        path: &std::path::Path,
        raw_rt_range: Option<[f64; 2]>,
    ) -> Result<LoadedCalibration, String> {
        let json = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let saved: SavedCalibration = serde_json::from_str(&json).map_err(|e| e.to_string())?;
        if saved.version != "v2" {
            return Err(format!(
                "Unsupported calibration version: {} (expected v2)",
                saved.version
            ));
        }

        let warning = match raw_rt_range {
            Some(raw) => {
                let overlap_lo = saved.rt_range_seconds[0].max(raw[0]);
                let overlap_hi = saved.rt_range_seconds[1].min(raw[1]);
                let overlap = (overlap_hi - overlap_lo).max(0.0);
                let span = saved.rt_range_seconds[1] - saved.rt_range_seconds[0];
                if span > 0.0 && overlap / span < 0.5 {
                    Some(
                        "RT range mismatch — calibration may not be valid for this file"
                            .to_string(),
                    )
                } else {
                    None
                }
            }
            None => Some("No raw file loaded — cannot verify RT range compatibility".to_string()),
        };

        Ok(LoadedCalibration {
            snapshot: saved.calibration,
            tolerances: saved.tolerances,
            errors: saved.errors,
            derivation: saved.derivation,
            n_calibrants: saved.n_calibrants,
            n_scored: saved.n_scored,
            warning,
        })
    }

    /// Fallback when calibration fails: identity RT mapping, secondary tolerance.
    pub fn fallback<I: ScorerQueriable>(pipeline: &Scorer<I>) -> Self {
        let range = pipeline.index.ms1_cycle_mapping().range_milis();
        let start = range.0 as f64 / 1000.0;
        let end = range.1 as f64 / 1000.0;
        let points = vec![
            Point {
                library: start,
                observed: start,
                weight: 1.0,
            },
            Point {
                library: end,
                observed: end,
                weight: 1.0,
            },
        ];
        let cal_curve = calibrate_with_ranges(&points, (start, end), (start, end), 10, 10)
            .expect("Identity calibration should not fail");

        Self {
            cal_curve,
            rt_tolerance_minutes: 1.0,
            mz_tolerance_ppm: (10.0, 10.0),
            mobility_tolerance_pct: (5.0, 5.0),
            ridge_widths: Vec::new(),
            errors: DimensionErrors::default(),
            derivation: None,
        }
    }
}

/// Summary of ridge width measurements for reporting.
pub struct RidgeWidthSummary {
    pub weighted_avg: f64,
    pub min: f64,
    pub max: f64,
    pub n_columns: usize,
    /// Fraction of total column weight that falls inside the ridge bounds (0.0–1.0).
    /// Higher = better fit between spectral library and raw file.
    pub in_ridge_ratio: f64,
}

/// JSON v2 calibration file format — shared between CLI and viewer.
///
/// Schema breaking change from v1: residual statistics (`errors`) and
/// derivation parameters (`derivation`) are now persisted, and `tolerances`
/// is a derived/convenience view rather than the source of truth.
#[derive(Debug, Serialize, Deserialize)]
pub struct SavedCalibration {
    pub version: String,
    pub rt_range_seconds: [f64; 2],
    pub calibration: CalibrationSnapshot,
    pub errors: DimensionErrors,
    pub derivation: DerivationParams,
    pub tolerances: SavedTolerances,
    pub n_calibrants: usize,
    pub n_scored: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SavedTolerances {
    pub rt_minutes: f32,
    pub mz_ppm: [f64; 2],
    pub mobility_pct: [f64; 2],
}

/// Result of loading a calibration JSON file.
pub struct LoadedCalibration {
    pub snapshot: CalibrationSnapshot,
    pub tolerances: SavedTolerances,
    pub errors: DimensionErrors,
    pub derivation: DerivationParams,
    pub n_calibrants: usize,
    pub n_scored: usize,
    pub warning: Option<String>,
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
