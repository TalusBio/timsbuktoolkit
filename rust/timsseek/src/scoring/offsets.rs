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

/// Number of top fragments contributing to the obs-mobility estimate.
/// Subset of FRAGMENT_TOP_N to bias toward the highest-confidence ions.
const FRAGMENT_OBS_MOB_TOP_N: usize = 3;

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

        // Guard the percent-error divide: a zero/absent reference mobility
        // (no-IM library, or an mzML sentinel) would yield `inf`, which slips
        // past the `is_nan()`-only accumulator guards in `weighted_ms1/ms2`
        // (Phase-2 calibration). NaN is the correct "no information" sentinel.
        let ref_mob_f32 = ref_mobility as f32;
        let ref_mob_valid = ref_mob_f32.is_finite() && ref_mob_f32 != 0.0;

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
            let mobility_error_pct = if ref_mob_valid {
                mobility_error_pct / ref_mob_f32 * 1e2
            } else {
                f32::NAN
            };
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
            let mobility_error_pct = if ref_mob_valid {
                mobility_error_pct / ref_mob_f32 * 1e2
            } else {
                f32::NAN
            };
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

    /// Intensity-weighted mean of precursor mz / mobility errors across the
    /// top-N (label >= 0) channels. Returns None when no channel has signal.
    /// Used by Phase-2 calibration to estimate population-level offsets;
    /// rescoring uses the per-ion arrays directly.
    pub fn weighted_ms1(&self) -> Option<(f32, f32)> {
        let (mut w_mz, mut mz) = (0.0f64, 0.0f64);
        let (mut w_mob, mut mob) = (0.0f64, 0.0f64);
        for v in self.ms1.get_values_sorted() {
            if v.intensity <= 0.0 {
                continue;
            }
            if !v.mz_error_ppm.is_nan() {
                w_mz += v.intensity;
                mz += v.intensity * v.mz_error_ppm as f64;
            }
            if !v.mobility_error_pct.is_nan() {
                w_mob += v.intensity;
                mob += v.intensity * v.mobility_error_pct as f64;
            }
        }
        // m/z and mobility are INDEPENDENT: a run with no searchable mobility
        // (mzML against a no-IM library) has all-NaN mobility errors, but its m/z
        // errors are perfectly good and MUST still calibrate. Gate each on its
        // own weight; return NaN mobility rather than dropping the whole calibrant.
        if w_mz > 0.0 {
            let mz_avg = (mz / w_mz) as f32;
            let mob_avg = if w_mob > 0.0 {
                (mob / w_mob) as f32
            } else {
                f32::NAN
            };
            Some((mz_avg, mob_avg))
        } else {
            None
        }
    }

    pub fn ms1_mz_errors(&self) -> [f32; PRECURSOR_TOP_N] {
        let mut out = [0.0; PRECURSOR_TOP_N];
        let vals = self.ms1.get_values_sorted();
        for i in 0..PRECURSOR_TOP_N {
            out[i] = vals[i].mz_error_ppm;
        }
        out
    }

    pub fn ms2_mz_errors(&self) -> [f32; FRAGMENT_TOP_N] {
        let mut out = [0.0; FRAGMENT_TOP_N];
        let vals = self.ms2.get_values_sorted();
        for i in 0..FRAGMENT_TOP_N {
            out[i] = vals[i].mz_error_ppm;
        }
        out
    }

    pub fn ms1_mobility_errors(&self) -> [f32; PRECURSOR_TOP_N] {
        let mut out = [0.0; PRECURSOR_TOP_N];
        let vals = self.ms1.get_values_sorted();
        for i in 0..PRECURSOR_TOP_N {
            out[i] = vals[i].mobility_error_pct;
        }
        out
    }

    pub fn ms2_mobility_errors(&self) -> [f32; FRAGMENT_TOP_N] {
        let mut out = [0.0; FRAGMENT_TOP_N];
        let vals = self.ms2.get_values_sorted();
        for i in 0..FRAGMENT_TOP_N {
            out[i] = vals[i].mobility_error_pct;
        }
        out
    }

    /// Intensity-weighted absolute mobility deltas (1/k0 units) for MS1 and MS2.
    /// Converts the stored percent error back to absolute via `ref_mobility`.
    /// The collector's "mz" slot carries the ppm error (unused by current
    /// consumers); only `mean_mobility()` is meaningful here.
    pub fn avg_delta_mobs(&self) -> (MzMobilityStatsCollector, MzMobilityStatsCollector) {
        let mut ms1 = MzMobilityStatsCollector::default();
        let mut ms2 = MzMobilityStatsCollector::default();
        let pct_to_abs = self.ref_mobility / 100.0;

        for v in self
            .ms2
            .get_values_sorted()
            .iter()
            .take(FRAGMENT_OBS_MOB_TOP_N)
        {
            if v.mobility_error_pct.is_nan() {
                continue;
            }
            ms2.add(
                v.intensity,
                v.mz_error_ppm as f64,
                v.mobility_error_pct as f64 * pct_to_abs,
            );
        }

        for v in self.ms1.get_values_sorted().iter() {
            if v.mobility_error_pct.is_nan() {
                continue;
            }
            ms1.add(
                v.intensity,
                v.mz_error_ppm as f64,
                v.mobility_error_pct as f64 * pct_to_abs,
            );
        }

        (ms1, ms2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ion(intensity: f64, mz_err: f32, mob_err: f32) -> ObsIonWithError {
        ObsIonWithError {
            intensity,
            mz_error_ppm: mz_err,
            mobility_error_pct: mob_err,
        }
    }

    #[test]
    fn weighted_ms1_calibrates_mz_when_mobility_absent() {
        // mzML vs no-IM library: finite m/z error, NaN mobility error. m/z MUST
        // still calibrate — the calibrant must not be dropped just because
        // mobility is absent (regression guard for the `w_mz && w_mob` gate).
        let mut ms1 = TopNArray::new();
        ms1.push(ion(10.0, 5.0, f32::NAN));
        ms1.push(ion(20.0, 7.0, f32::NAN));
        let offsets = MzMobilityOffsets {
            ms1,
            ms2: TopNArray::new(),
            ref_mobility: 0.0,
        };
        let (mz, mob) = offsets
            .weighted_ms1()
            .expect("m/z must calibrate even when mobility is absent");
        assert!(mz.is_finite() && mz > 0.0, "mz avg should be finite: {mz}");
        assert!(mob.is_nan(), "mobility must stay NaN, not fabricated: {mob}");
    }

    #[test]
    fn weighted_ms1_none_without_mz_signal() {
        let offsets = MzMobilityOffsets {
            ms1: TopNArray::new(),
            ms2: TopNArray::new(),
            ref_mobility: 0.9,
        };
        assert!(offsets.weighted_ms1().is_none());
    }
}
