use super::cv::{
    CrossValidatedScorer,
    DataBuffer,
    FeatureLike,
    GBMConfig,
};
use super::{
    LabelledScore,
    TargetDecoy,
};
use rand::prelude::*;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use tracing::debug;

/// Assign q_values in place.
///
/// # Invariants
/// * `scores` must be sorted in descending order (e.g. best PSM is first)
///
/// Implementation derived from the Sage implementation of qval (Thanks Mike) github.com/lazear/sage
fn assign_qval<T: LabelledScore>(scores: &mut [T], key: impl Fn(&T) -> f32) {
    assert!(
        scores.windows(2).all(|w| key(&w[0]) >= key(&w[1])),
        "Expecting scores to be sorted in descending order",
    );

    let mut decoy = 1;
    let mut target = 0;

    for score in scores.iter_mut() {
        match score.get_label() {
            TargetDecoy::Decoy => decoy += 1,
            TargetDecoy::Target => target += 1,
        }
        score.assign_qval((decoy as f32) / (target as f32));
    }

    // Reverse slice, and calculate the cumulative minimum
    let mut q_min = 1.0f32;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.get_qval());
        score.assign_qval(q_min);
    }

    // We do a third pass to ensure that values with the same score have the same q-value
    let mut last_score = f32::NAN;
    let mut last_qval = 1.0f32;
    for score in scores.iter_mut() {
        let current_score = key(score);
        if current_score != last_score {
            last_score = current_score;
            last_qval = score.get_qval();
            continue;
        }
        score.assign_qval(last_qval);
    }
}

pub fn report_qvalues_at_thresholds<T: LabelledScore + std::fmt::Debug>(
    scores: &[T],
    thresholds: &[f32],
) -> Vec<(f32, usize, usize, usize)> {
    let mut out = Vec::new();

    for &thresh in thresholds {
        let n_below_thresh = scores.iter().filter(|s| s.get_qval() <= thresh).count();
        let n_targets = scores
            .iter()
            .filter(|s| s.get_qval() <= thresh && matches!(s.get_label(), TargetDecoy::Target))
            .count();
        let n_decoys = scores
            .iter()
            .filter(|s| s.get_qval() <= thresh && matches!(s.get_label(), TargetDecoy::Decoy))
            .count();
        out.push((thresh, n_below_thresh, n_targets, n_decoys));
    }

    out
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
/// Fixed shuffle seed used by `rescore`. Makes the pre-rescore shuffle
/// (and therefore the fold assignment + downstream target counts)
/// reproducible across runs, eliminating RNG-driven noise in benches.
/// `GBMConfig::default().seed == 0` already makes the boosting itself
/// deterministic; this seals the only remaining entropy source.
const RESCORE_SHUFFLE_SEED: u64 = 42;

pub fn rescore(mut data: Vec<CompetedCandidate>) -> Vec<FinalResult> {
    let config = GBMConfig::default();

    // Canonicalize input order before the seeded shuffle. Upstream
    // stages can emit candidates in an order that drifts with
    // floating-point accumulation quirks (e.g. different peak-bucket
    // layouts produce identical features but different vec orderings).
    // Without a stable sort here, the seeded shuffle sees different
    // inputs across runs -> different fold assignment -> different
    // q-values -> drifting target counts across equivalent configs.
    // (library_id, precursor_charge) is a non-FP composite key that
    // should be unique per candidate after target-decoy competition.
    data.sort_unstable_by_key(|c| (c.scoring.library_id, c.scoring.precursor_charge));

    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(RESCORE_SHUFFLE_SEED);
    data.shuffle(&mut rng);

    let mut scorer = CrossValidatedScorer::<CompetedCandidate>::new_from_shuffled(3, data, config);
    scorer
        .fit(&mut DataBuffer::default(), &mut DataBuffer::default())
        .unwrap();

    let mut scored = scorer.score();
    // Sort by score descending
    #[cfg(feature = "rayon")]
    scored.par_sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    #[cfg(not(feature = "rayon"))]
    scored.sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut scored, |x| CompetedCandidate::get_score(x) as f32);
    debug!("Best:\n{:#?}", scored.first());
    debug!("Worst:\n{:#?}", scored.last());

    scored.into_iter().map(|c| c.into_final()).collect()
}

use crate::scoring::results::{
    CompetedCandidate,
    FinalResult,
};

fn mean_abs_error(errs: &[f32]) -> f64 {
    let (sum, n) = errs
        .iter()
        .filter(|e| e.is_finite() && **e != 0.0)
        .fold((0.0f64, 0u32), |(s, n), &e| (s + (e as f64).abs(), n + 1));
    if n > 0 { sum / n as f64 } else { f64::NAN }
}

// ---------------------------------------------------------------------------
// CompetedCandidate: FeatureLike + LabelledScore
// ---------------------------------------------------------------------------

impl FeatureLike for CompetedCandidate {
    fn as_feature(&self) -> impl IntoIterator<Item = f64> + '_ {
        let s = &self.scoring;

        vec![
            (s.precursor_mz / 5.0).round(),
            s.precursor_charge as f64,
            s.precursor_mobility as f64,
            s.calibrated_rt_seconds.round() as f64,
            s.n_scored_fragments as f64,
            // Combined
            s.main_score as f64,
            (s.main_score / s.delta_next) as f64,
            s.delta_next as f64,
            s.delta_second_next as f64,
            s.obs_rt_seconds as f64,
            s.obs_mobility as f64,
            (s.obs_rt_seconds - s.calibrated_rt_seconds) as f64,
            s.calibrated_sq_delta_rt as f64,
            s.delta_ms1_ms2_mobility as f64,
            s.sq_delta_ms1_ms2_mobility as f64,
            s.rising_cycles as f64,
            s.falling_cycles as f64,
            // MS2
            s.npeaks as f64,
            s.apex_lazyscore as f64,
            (s.ms2_summed_intensity as f64).ln_1p(),
            s.ms2_lazyscore as f64,
            s.ms2_isotope_lazyscore as f64,
            s.ms2_isotope_lazyscore_ratio as f64,
            s.lazyscore_z as f64,
            s.lazyscore_vs_baseline as f64,
            // Split product & apex features
            (s.split_product_score as f64).ln_1p(),
            (s.cosine_au as f64).ln_1p(),
            (s.scribe_au as f64).ln_1p(),
            s.cosine_cg as f64,
            s.scribe_cg as f64,
            s.cosine_weighted_coelution as f64,
            s.cosine_gradient_consistency as f64,
            s.scribe_weighted_coelution as f64,
            s.scribe_gradient_consistency as f64,
            s.peak_shape as f64,
            s.ratio_cv as f64,
            s.centered_apex as f64,
            s.precursor_coelution as f64,
            s.fragment_coverage as f64,
            s.precursor_apex_match as f64,
            s.xic_quality as f64,
            s.fragment_apex_agreement as f64,
            s.isotope_correlation as f64,
            s.gaussian_correlation as f64,
            s.per_frag_gaussian_corr as f64,
            // MS2 per-ion errors
            s.ms2_mz_errors[0] as f64,
            s.ms2_mz_errors[1] as f64,
            s.ms2_mz_errors[2] as f64,
            s.ms2_mz_errors[3] as f64,
            s.ms2_mz_errors[4] as f64,
            s.ms2_mz_errors[5] as f64,
            s.ms2_mz_errors[6] as f64,
            s.ms2_mobility_errors[0] as f64,
            s.ms2_mobility_errors[1] as f64,
            s.ms2_mobility_errors[2] as f64,
            s.ms2_mobility_errors[3] as f64,
            s.ms2_mobility_errors[4] as f64,
            s.ms2_mobility_errors[5] as f64,
            s.ms2_mobility_errors[6] as f64,
            // MS1
            (s.ms1_summed_intensity as f64).ln_1p(),
            // MS1 per-ion errors
            s.ms1_mz_errors[0] as f64,
            s.ms1_mz_errors[1] as f64,
            s.ms1_mz_errors[2] as f64,
            s.ms1_mobility_errors[0] as f64,
            s.ms1_mobility_errors[1] as f64,
            s.ms1_mobility_errors[2] as f64,
            // Relative intensities
            s.ms1_intensity_ratios[0] as f64,
            s.ms1_intensity_ratios[1] as f64,
            s.ms1_intensity_ratios[2] as f64,
            s.ms2_intensity_ratios[0] as f64,
            s.ms2_intensity_ratios[1] as f64,
            s.ms2_intensity_ratios[2] as f64,
            s.ms2_intensity_ratios[3] as f64,
            s.ms2_intensity_ratios[4] as f64,
            s.ms2_intensity_ratios[5] as f64,
            s.ms2_intensity_ratios[6] as f64,
            self.delta_group as f64,
            self.delta_group_ratio as f64,
            s.calibrated_rt_seconds as f64,
            s.calibrated_sq_delta_rt as f64,
            // Derived intensity features
            {
                let ratios = &s.ms2_intensity_ratios;
                ratios
                    .iter()
                    .filter(|r| r.is_finite())
                    .fold(f32::NEG_INFINITY, |a, &b| a.max(b)) as f64
            },
            // Interaction features
            (s.main_score * s.delta_next) as f64,
            (s.split_product_score * s.fragment_coverage) as f64,
            // Summary error features
            mean_abs_error(&s.ms2_mz_errors),
            mean_abs_error(&s.ms2_mobility_errors),
            mean_abs_error(&s.ms1_mz_errors),
            mean_abs_error(&s.ms1_mobility_errors),
        ]
    }

    fn get_y(&self) -> f64 {
        if self.scoring.is_target { 1.0 } else { 0.0 }
    }

    fn assign_score(&mut self, score: f64) {
        self.discriminant_score = score as f32;
    }

    fn get_score(&self) -> f64 {
        self.discriminant_score as f64
    }
}

impl LabelledScore for CompetedCandidate {
    fn get_label(&self) -> TargetDecoy {
        if self.scoring.is_target {
            TargetDecoy::Target
        } else {
            TargetDecoy::Decoy
        }
    }

    fn assign_qval(&mut self, q: f32) {
        self.qvalue = q;
    }

    fn get_qval(&self) -> f32 {
        self.qvalue
    }
}

impl LabelledScore for FinalResult {
    fn get_label(&self) -> TargetDecoy {
        if self.scoring.is_target {
            TargetDecoy::Target
        } else {
            TargetDecoy::Decoy
        }
    }

    fn assign_qval(&mut self, q: f32) {
        self.qvalue = q;
    }

    fn get_qval(&self) -> f32 {
        self.qvalue
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // scores = np.array([10, 10, 9, 8, 7, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1])
    // target = np.array([1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0])
    // qvals = np.array([
    //     1 / 4,
    //     1 / 4,
    //     1 / 4,
    //     1 / 4,
    //     2 / 6,
    //     2 / 6,
    //     2 / 6,
    //     3 / 7,
    //     3 / 7,
    //     4 / 7,
    //     5 / 8,
    //     5 / 8,
    //     1,
    //     1,
    //     1,
    //     1,
    // ])

    #[test]
    fn test_assign_qval() {
        let scores = vec![10, 10, 9, 8, 7, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1];
        let target = vec![1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0];
        let qvals = vec![
            1. / 4.,
            1. / 4.,
            1. / 4.,
            1. / 4.,
            2. / 6.,
            2. / 6.,
            2. / 6.,
            3. / 7.,
            3. / 7.,
            4. / 8.,   // 4. / 7.,
            4. / 8.,   // 5. / 8.,
            4. / 8.,   // 5. / 8.,
            6.0 / 8.0, // 1.,
            6.0 / 8.0, // 1.,
            6.0 / 8.0, // 1.,
            6.0 / 8.0, // 1.,
        ];
        struct TestScore {
            score: f64,
            label: TargetDecoy,
            qval: f32,
        }

        impl LabelledScore for TestScore {
            fn get_label(&self) -> TargetDecoy {
                self.label
            }

            fn assign_qval(&mut self, q: f32) {
                self.qval = q
            }

            fn get_qval(&self) -> f32 {
                self.qval
            }
        }

        let mut data = scores
            .iter()
            .zip(target.iter())
            .map(|(&s, &l)| TestScore {
                score: s as f64,
                label: if l == 1 {
                    TargetDecoy::Target
                } else {
                    TargetDecoy::Decoy
                },
                qval: 0.0,
            })
            .collect::<Vec<_>>();

        assign_qval(&mut data, |x| x.score as f32);

        for i in 0..qvals.len() {
            let real = qvals[i];
            let model = data[i].get_qval();

            assert!(
                (real - model).abs() < 1e-6,
                "At index {}: expected {}, got {}; REAL: {:?} MODEL {:?}",
                i,
                real,
                model,
                qvals,
                data.iter().map(|x| x.get_qval()).collect::<Vec<_>>(),
            );
        }
    }
}
