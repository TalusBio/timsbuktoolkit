//! Cross-field / cross-block ML features. These reference more than one field,
//! so they are NOT `#[feat]` attrs (which are strictly unary) — instead a
//! hand-written [`Derived::compute`] reads a view of the finalized
//! [`ScoringFields`] and fills a normal `score_block!` struct. The block emits
//! ML features/names through the same walks as every other block (so value and
//! name can't desync), but is NOT part of `compose_scoring_fields!`, so
//! `columns()` is never called and these stay ML-only (no parquet column).
//! Emitted after the per-block features, before the sequence block.

use timsseek_macros::ScoreBlock;

use crate::scoring::results::ScoringFields;

/// Cross-field interaction features (ML-only). Field names ARE the feature
/// names; order is the emission order.
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct Derived {
    #[feat(log2)]
    pub main_over_delta_next: f64,
    #[feat(raw, abs)]
    pub rt_err: f64,
    #[feat(raw)]
    pub ms2_intensity_ratios_max: f64,
    #[feat(log2)]
    pub main_times_delta_next: f64,
    #[feat(log2)]
    pub split_product_x_coverage: f64,
    #[feat(raw)]
    pub ms2_mz_mean_abs_error: f64,
    #[feat(raw)]
    pub ms2_mob_mean_abs_error: f64,
    #[feat(raw)]
    pub ms1_mz_mean_abs_error: f64,
    #[feat(raw)]
    pub ms1_mob_mean_abs_error: f64,
}

impl Derived {
    /// Derive the cross-field features from a finalized [`ScoringFields`].
    pub fn compute(s: &ScoringFields) -> Self {
        Self {
            main_over_delta_next: (s.primary.main_score / s.primary.delta_next) as f64,
            rt_err: (s.rt.obs_rt_seconds - s.rt.calibrated_rt_seconds) as f64,
            ms2_intensity_ratios_max: s
                .rel_intensities
                .ms2_intensity_ratios
                .iter()
                .filter(|r| r.is_finite())
                .fold(f32::NEG_INFINITY, |a, &b| a.max(b))
                as f64,
            main_times_delta_next: (s.primary.main_score * s.primary.delta_next) as f64,
            split_product_x_coverage: (s.split.split_product_score * s.features.fragment_coverage)
                as f64,
            ms2_mz_mean_abs_error: mean_abs_error(&s.ion_errors.ms2_mz_errors),
            ms2_mob_mean_abs_error: mean_abs_error(&s.ion_errors.ms2_mobility_errors),
            ms1_mz_mean_abs_error: mean_abs_error(&s.ion_errors.ms1_mz_errors),
            ms1_mob_mean_abs_error: mean_abs_error(&s.ion_errors.ms1_mobility_errors),
        }
    }
}

/// Mean of the absolute, finite, non-zero errors (NaN if none qualify).
fn mean_abs_error(errs: &[f32]) -> f64 {
    let (sum, n) = errs
        .iter()
        .filter(|e| e.is_finite() && **e != 0.0)
        .fold((0.0f64, 0u32), |(s, n), &e| (s + (e as f64).abs(), n + 1));
    if n > 0 { sum / n as f64 } else { f64::NAN }
}
