//! Cross-field / cross-block ML features. These reference more than one field,
//! so they are NOT `#[feat]` attrs (which are strictly unary) and NOT
//! compile-exhaustive — they are open-ended interactions, guarded by the golden
//! feature-name-set test. Emitted after the per-block features, before the
//! sequence block.

use crate::scoring::blocks::FeatSink;
use crate::scoring::results::ScoringFields;

/// Mean of the absolute, finite, non-zero errors (NaN if none qualify).
fn mean_abs_error(errs: &[f32]) -> f64 {
    let (sum, n) = errs
        .iter()
        .filter(|e| e.is_finite() && **e != 0.0)
        .fold((0.0f64, 0u32), |(s, n), &e| (s + (e as f64).abs(), n + 1));
    if n > 0 { sum / n as f64 } else { f64::NAN }
}

pub fn features(s: &ScoringFields, o: &mut FeatSink) {
    o.push(
        "main_over_delta_next",
        (s.primary.main_score / s.primary.delta_next) as f64,
    );
    o.push(
        "rt_err",
        (s.rt.obs_rt_seconds - s.rt.calibrated_rt_seconds) as f64,
    );
    o.push("ms2_intensity_ratios_max", {
        let ratios = &s.rel_intensities.ms2_intensity_ratios;
        ratios
            .iter()
            .filter(|r| r.is_finite())
            .fold(f32::NEG_INFINITY, |a, &b| a.max(b)) as f64
    });
    o.push(
        "main_times_delta_next",
        (s.primary.main_score * s.primary.delta_next) as f64,
    );
    o.push(
        "split_product_x_coverage",
        (s.split.split_product_score * s.features.fragment_coverage) as f64,
    );
    o.push(
        "ms2_mz_mean_abs_error",
        mean_abs_error(&s.ion_errors.ms2_mz_errors),
    );
    o.push(
        "ms2_mob_mean_abs_error",
        mean_abs_error(&s.ion_errors.ms2_mobility_errors),
    );
    o.push(
        "ms1_mz_mean_abs_error",
        mean_abs_error(&s.ion_errors.ms1_mz_errors),
    );
    o.push(
        "ms1_mob_mean_abs_error",
        mean_abs_error(&s.ion_errors.ms1_mobility_errors),
    );
}
