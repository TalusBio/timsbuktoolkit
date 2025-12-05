use timsquery::{
    MzMobilityStatsCollector,
    SpectralCollector,
};

use crate::IonAnnot;
use crate::utils::top_n_array::TopNArray;

/// Number of top precursor isotopes to track for error calculation.
/// Typically captures M+0, M+1, M+2 isotopes for calibration purposes.
const PRECURSOR_TOP_N: usize = 3;

/// Number of top fragment ions to track for error calculation.
/// Balances statistical power with outlier resistance.
const FRAGMENT_TOP_N: usize = 7;

/// Container for measured m/z and mobility offsets from top ions.
///
/// Tracks the highest-intensity precursors and fragments to calculate
/// systematic measurement errors for calibration and quality assessment.
#[derive(Debug)]
pub struct MzMobilityOffsets {
    /// Top precursor measurements (highest intensity first).
    pub(crate) ms1: TopNArray<PRECURSOR_TOP_N, ObsIonWithError>,

    /// Top fragment measurements (highest intensity first).
    pub(crate) ms2: TopNArray<FRAGMENT_TOP_N, ObsIonWithError>,

    /// Reference mobility used for error calculation.
    pub(crate) ref_mobility: f64,
}

/// Observed ion measurement with deviations from theoretical values.
///
/// Represents a single ion measurement including its observed intensity and
/// the calculated errors relative to theoretical m/z and mobility predictions.
/// Sorted by intensity to prioritize high-quality measurements.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ObsIonWithError {
    /// Observed ion intensity (weight from aggregator).
    pub(crate) intensity: f64,

    /// Mass-to-charge ratio error in parts per million (ppm).
    /// Positive means observed m/z is higher than theoretical.
    pub(crate) mz_error_ppm: f32,

    /// Ion mobility error as percentage of reference mobility.
    /// Positive means observed mobility is higher than theoretical.
    pub(crate) mobility_error_pct: f32,
}

impl PartialOrd for ObsIonWithError {
    fn partial_cmp(&self, other: &ObsIonWithError) -> Option<std::cmp::Ordering> {
        self.intensity.partial_cmp(&other.intensity)
    }
}

impl Default for ObsIonWithError {
    fn default() -> Self {
        Self {
            intensity: 0.0,
            mz_error_ppm: f32::NAN,
            mobility_error_pct: f32::NAN,
        }
    }
}

impl MzMobilityOffsets {
    pub fn new(
        item: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
        ref_mobility: f64,
    ) -> Self {
        let mut ms1 = TopNArray::new();
        let mut ms2 = TopNArray::new();

        for ((key, ref_mz), val) in item.iter_precursors() {
            if key < 0i8 {
                continue;
            }
            let intensity = val.weight();
            let mz_error_ppm = (val.mean_mz().unwrap_or(f64::NAN) - ref_mz) as f32;
            // Convert to PPM
            let mz_error_ppm = mz_error_ppm / (ref_mz as f32) * 1e6;
            let mobility_error_pct =
                (val.mean_mobility().unwrap_or(f64::NAN) - ref_mobility) as f32;
            // Convert to percentage
            let mobility_error_pct = mobility_error_pct / (ref_mobility as f32) * 1e2;
            ms1.push(ObsIonWithError {
                intensity,
                mz_error_ppm,
                mobility_error_pct,
            });
        }

        for ((_key, ref_mz), val) in item.iter_fragments() {
            let intensity = val.weight();
            let mz_error_ppm = (val.mean_mz().unwrap_or(f64::NAN) - ref_mz) as f32;
            let mz_error_ppm = mz_error_ppm / (*ref_mz as f32) * 1e6;
            let mobility_error_pct =
                (val.mean_mobility().unwrap_or(f64::NAN) - ref_mobility) as f32;
            let mobility_error_pct = mobility_error_pct / (ref_mobility as f32) * 1e2;
            ms2.push(ObsIonWithError {
                intensity,
                mz_error_ppm,
                mobility_error_pct,
            });
        }

        Self {
            ms1,
            ms2,
            ref_mobility,
        }
    }

    pub fn ms1_mz_errors(&self) -> [f32; PRECURSOR_TOP_N] {
        let mut out = [0.0; PRECURSOR_TOP_N];
        let vals = self.ms1.get_values();
        for i in 0..PRECURSOR_TOP_N {
            out[i] = vals[i].mz_error_ppm;
        }
        out
    }

    pub fn ms2_mz_errors(&self) -> [f32; FRAGMENT_TOP_N] {
        let mut out = [0.0; FRAGMENT_TOP_N];
        let vals = self.ms2.get_values();
        for i in 0..FRAGMENT_TOP_N {
            out[i] = vals[i].mz_error_ppm;
        }
        out
    }

    pub fn ms1_mobility_errors(&self) -> [f32; PRECURSOR_TOP_N] {
        let mut out = [0.0; PRECURSOR_TOP_N];
        let vals = self.ms1.get_values();
        for i in 0..PRECURSOR_TOP_N {
            out[i] = vals[i].mobility_error_pct;
        }
        out
    }

    pub fn ms2_mobility_errors(&self) -> [f32; FRAGMENT_TOP_N] {
        let mut out = [0.0; FRAGMENT_TOP_N];
        let vals = self.ms2.get_values();
        for i in 0..FRAGMENT_TOP_N {
            out[i] = vals[i].mobility_error_pct;
        }
        out
    }

    pub fn avg_delta_mobs(&self) -> (MzMobilityStatsCollector, MzMobilityStatsCollector) {
        let mut ms2 = MzMobilityStatsCollector::default();
        let mut ms1 = MzMobilityStatsCollector::default();
        let vals = self.ms2.get_values();
        for v in vals.iter().take(3) {
            if v.mobility_error_pct.is_nan() {
                continue;
            }
            ms2.add(
                v.intensity,
                v.mz_error_ppm as f64,
                v.mobility_error_pct as f64,
            );
        }

        let vals = self.ms1.get_values();
        for v in vals.iter() {
            if v.mobility_error_pct.is_nan() {
                continue;
            }
            ms1.add(
                v.intensity,
                v.mz_error_ppm as f64,
                v.mobility_error_pct as f64,
            );
        }

        (ms1, ms2)
    }
}
