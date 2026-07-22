//! Cross-field / cross-block ML features. These reference more than one field,
//! so they are NOT `#[feat]` attrs (which are strictly unary) and NOT
//! compile-exhaustive — they are open-ended interactions. Names are codegen
//! from the `NAMES` table (no golden). Emitted after the per-block features,
//! before the sequence block.

use crate::scoring::blocks::{
    FeatSink,
    NameSink,
};
use crate::scoring::results::ScoringFields;

/// Feature names emitted by [`features`], in the same order. Single source for
/// both the value walk and the set-level name walk ([`feature_names`]).
const NAMES: [&str; 9] = [
    "main_over_delta_next",
    "rt_err",
    "ms2_intensity_ratios_max",
    "main_times_delta_next",
    "split_product_x_coverage",
    "ms2_mz_mean_abs_error",
    "ms2_mob_mean_abs_error",
    "ms1_mz_mean_abs_error",
    "ms1_mob_mean_abs_error",
];

/// Mean of the absolute, finite, non-zero errors (NaN if none qualify).
fn mean_abs_error(errs: &[f32]) -> f64 {
    let (sum, n) = errs
        .iter()
        .filter(|e| e.is_finite() && **e != 0.0)
        .fold((0.0f64, 0u32), |(s, n), &e| (s + (e as f64).abs(), n + 1));
    if n > 0 { sum / n as f64 } else { f64::NAN }
}

/// Emit the cross-field feature *values*, in `NAMES` order.
pub fn features(s: &ScoringFields, o: &mut FeatSink) {
    o.push((s.primary.main_score / s.primary.delta_next) as f64);
    o.push((s.rt.obs_rt_seconds - s.rt.calibrated_rt_seconds) as f64);
    o.push({
        let ratios = &s.rel_intensities.ms2_intensity_ratios;
        ratios
            .iter()
            .filter(|r| r.is_finite())
            .fold(f32::NEG_INFINITY, |a, &b| a.max(b)) as f64
    });
    o.push((s.primary.main_score * s.primary.delta_next) as f64);
    o.push((s.split.split_product_score * s.features.fragment_coverage) as f64);
    o.push(mean_abs_error(&s.ion_errors.ms2_mz_errors));
    o.push(mean_abs_error(&s.ion_errors.ms2_mobility_errors));
    o.push(mean_abs_error(&s.ion_errors.ms1_mz_errors));
    o.push(mean_abs_error(&s.ion_errors.ms1_mobility_errors));
}

/// Emit the cross-field feature *names* (set-level).
pub fn feature_names(o: &mut NameSink) {
    for n in NAMES {
        o.push(n);
    }
}
