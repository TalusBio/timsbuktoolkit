use crate::ScorerQueriable;
use crate::scoring::pipeline::Scorer;
pub use calibrt::{
    CalibRtError,
    CalibrationCurve as RTCalibration,
    Point,
    calibrate_with_ranges,
};
use timsquery::Tolerance;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};

/// Immutable calibration result. Provides RT conversion and per-query tolerance
/// without mutating the speclib.
pub struct CalibrationResult {
    cal_curve: RTCalibration,
    rt_tolerance_minutes: f32,
    mz_tolerance_ppm: (f64, f64),
    mobility_tolerance_pct: (f32, f32),
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
        }
    }

    /// Convert indexed RT to calibrated absolute RT (seconds).
    pub fn convert_irt(&self, irt_seconds: f32) -> f32 {
        match self.cal_curve.predict(irt_seconds as f64) {
            Ok(rt) => rt as f32,
            Err(CalibRtError::OutOfBounds(rt)) => rt as f32,
            Err(_) => irt_seconds,
        }
    }

    /// Get per-query tolerance. Initially uniform; future: position-dependent.
    pub fn get_tolerance(&self, _mz: f64, _mobility: f32, _rt: f32) -> Tolerance {
        Tolerance {
            ms: MzTolerance::Ppm(self.mz_tolerance_ppm),
            rt: RtTolerance::Minutes((self.rt_tolerance_minutes, self.rt_tolerance_minutes)),
            mobility: MobilityTolerance::Pct(self.mobility_tolerance_pct),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
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

    /// Fallback when calibration fails: identity RT mapping, secondary tolerance.
    pub fn fallback<I: ScorerQueriable>(pipeline: &Scorer<I>) -> Self {
        let range = pipeline.index.ms1_cycle_mapping().range_milis();
        let start = range.0 as f64 / 1000.0;
        let end = range.1 as f64 / 1000.0;
        let points = vec![
            Point {
                x: start,
                y: start,
                weight: 1.0,
            },
            Point {
                x: end,
                y: end,
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
        }
    }
}

