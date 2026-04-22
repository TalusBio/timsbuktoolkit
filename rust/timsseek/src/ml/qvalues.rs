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

/// Canonical 20-dim AA-count feature names, aligned with
/// `ParsedSequence::aa_counts` output order (see `CANONICAL_AA_LETTERS`):
/// A C D E F G H I K L M N P Q R S T V W Y.
pub const AA_COUNT_NAMES: [&str; 20] = [
    "aa_count_A",
    "aa_count_C",
    "aa_count_D",
    "aa_count_E",
    "aa_count_F",
    "aa_count_G",
    "aa_count_H",
    "aa_count_I",
    "aa_count_K",
    "aa_count_L",
    "aa_count_M",
    "aa_count_N",
    "aa_count_P",
    "aa_count_Q",
    "aa_count_R",
    "aa_count_S",
    "aa_count_T",
    "aa_count_V",
    "aa_count_W",
    "aa_count_Y",
];

impl CompetedCandidate {
    /// Single source of truth for feature values and names.
    ///
    /// Base block: 88 dims always present.
    /// Sequence block: 22 dims appended when `peptide.parsed.is_some()`.
    /// The gate is speclib-wide: either all candidates in a run have Some,
    /// or all have None — so vector length is stable within a single fit.
    fn named_features(&self) -> Vec<(f64, &'static str)> {
        let s = &self.scoring;
        let mut v: Vec<(f64, &'static str)> = vec![
            // Identity / precursor
            ((s.precursor_mz / 5.0).round(), "precursor_mz_round5"),
            (s.precursor_charge as f64, "precursor_charge"),
            (s.precursor_mobility as f64, "precursor_mobility"),
            (
                s.calibrated_rt_seconds.round() as f64,
                "calibrated_rt_seconds_round",
            ),
            (s.n_scored_fragments as f64, "n_scored_fragments"),
            // Combined
            (s.main_score as f64, "main_score"),
            ((s.main_score / s.delta_next) as f64, "main_over_delta_next"),
            (s.delta_next as f64, "delta_next"),
            (s.delta_second_next as f64, "delta_second_next"),
            (s.obs_rt_seconds as f64, "obs_rt_seconds"),
            (s.obs_mobility as f64, "obs_mobility"),
            (
                (s.obs_rt_seconds - s.calibrated_rt_seconds) as f64,
                "rt_err",
            ),
            (s.calibrated_sq_delta_rt as f64, "calibrated_sq_delta_rt"),
            (s.delta_ms1_ms2_mobility as f64, "delta_ms1_ms2_mobility"),
            (
                s.sq_delta_ms1_ms2_mobility as f64,
                "sq_delta_ms1_ms2_mobility",
            ),
            (s.rising_cycles as f64, "rising_cycles"),
            (s.falling_cycles as f64, "falling_cycles"),
            // MS2
            (s.npeaks as f64, "npeaks"),
            (s.apex_lazyscore as f64, "apex_lazyscore"),
            (
                (s.ms2_summed_intensity as f64).ln_1p(),
                "ms2_summed_intensity_ln1p",
            ),
            (s.ms2_lazyscore as f64, "ms2_lazyscore"),
            (s.ms2_isotope_lazyscore as f64, "ms2_isotope_lazyscore"),
            (
                s.ms2_isotope_lazyscore_ratio as f64,
                "ms2_isotope_lazyscore_ratio",
            ),
            (s.lazyscore_z as f64, "lazyscore_z"),
            (s.lazyscore_vs_baseline as f64, "lazyscore_vs_baseline"),
            // Split product & apex features
            (
                (s.split_product_score as f64).ln_1p(),
                "split_product_score_ln1p",
            ),
            ((s.cosine_au as f64).ln_1p(), "cosine_au_ln1p"),
            ((s.scribe_au as f64).ln_1p(), "scribe_au_ln1p"),
            (s.cosine_cg as f64, "cosine_cg"),
            (s.scribe_cg as f64, "scribe_cg"),
            (
                s.cosine_weighted_coelution as f64,
                "cosine_weighted_coelution",
            ),
            (
                s.cosine_gradient_consistency as f64,
                "cosine_gradient_consistency",
            ),
            (
                s.scribe_weighted_coelution as f64,
                "scribe_weighted_coelution",
            ),
            (
                s.scribe_gradient_consistency as f64,
                "scribe_gradient_consistency",
            ),
            (s.peak_shape as f64, "peak_shape"),
            (s.ratio_cv as f64, "ratio_cv"),
            (s.centered_apex as f64, "centered_apex"),
            (s.precursor_coelution as f64, "precursor_coelution"),
            (s.fragment_coverage as f64, "fragment_coverage"),
            (s.precursor_apex_match as f64, "precursor_apex_match"),
            (s.xic_quality as f64, "xic_quality"),
            (s.fragment_apex_agreement as f64, "fragment_apex_agreement"),
            (s.isotope_correlation as f64, "isotope_correlation"),
            (s.gaussian_correlation as f64, "gaussian_correlation"),
            (s.per_frag_gaussian_corr as f64, "per_frag_gaussian_corr"),
            // MS2 per-ion errors (7 mz + 7 mobility)
            (s.ms2_mz_errors[0] as f64, "ms2_mz_err_0"),
            (s.ms2_mz_errors[1] as f64, "ms2_mz_err_1"),
            (s.ms2_mz_errors[2] as f64, "ms2_mz_err_2"),
            (s.ms2_mz_errors[3] as f64, "ms2_mz_err_3"),
            (s.ms2_mz_errors[4] as f64, "ms2_mz_err_4"),
            (s.ms2_mz_errors[5] as f64, "ms2_mz_err_5"),
            (s.ms2_mz_errors[6] as f64, "ms2_mz_err_6"),
            (s.ms2_mobility_errors[0] as f64, "ms2_mob_err_0"),
            (s.ms2_mobility_errors[1] as f64, "ms2_mob_err_1"),
            (s.ms2_mobility_errors[2] as f64, "ms2_mob_err_2"),
            (s.ms2_mobility_errors[3] as f64, "ms2_mob_err_3"),
            (s.ms2_mobility_errors[4] as f64, "ms2_mob_err_4"),
            (s.ms2_mobility_errors[5] as f64, "ms2_mob_err_5"),
            (s.ms2_mobility_errors[6] as f64, "ms2_mob_err_6"),
            // MS1
            (
                (s.ms1_summed_intensity as f64).ln_1p(),
                "ms1_summed_intensity_ln1p",
            ),
            // MS1 per-ion errors (3 mz + 3 mobility)
            (s.ms1_mz_errors[0] as f64, "ms1_mz_err_0"),
            (s.ms1_mz_errors[1] as f64, "ms1_mz_err_1"),
            (s.ms1_mz_errors[2] as f64, "ms1_mz_err_2"),
            (s.ms1_mobility_errors[0] as f64, "ms1_mob_err_0"),
            (s.ms1_mobility_errors[1] as f64, "ms1_mob_err_1"),
            (s.ms1_mobility_errors[2] as f64, "ms1_mob_err_2"),
            // Relative intensities
            (s.ms1_intensity_ratios[0] as f64, "ms1_intensity_ratio_0"),
            (s.ms1_intensity_ratios[1] as f64, "ms1_intensity_ratio_1"),
            (s.ms1_intensity_ratios[2] as f64, "ms1_intensity_ratio_2"),
            (s.ms2_intensity_ratios[0] as f64, "ms2_intensity_ratio_0"),
            (s.ms2_intensity_ratios[1] as f64, "ms2_intensity_ratio_1"),
            (s.ms2_intensity_ratios[2] as f64, "ms2_intensity_ratio_2"),
            (s.ms2_intensity_ratios[3] as f64, "ms2_intensity_ratio_3"),
            (s.ms2_intensity_ratios[4] as f64, "ms2_intensity_ratio_4"),
            (s.ms2_intensity_ratios[5] as f64, "ms2_intensity_ratio_5"),
            (s.ms2_intensity_ratios[6] as f64, "ms2_intensity_ratio_6"),
            (self.delta_group as f64, "delta_group"),
            (self.delta_group_ratio as f64, "delta_group_ratio"),
            (s.calibrated_rt_seconds as f64, "calibrated_rt_seconds"),
            (
                s.calibrated_sq_delta_rt as f64,
                "calibrated_sq_delta_rt_dup",
            ),
            // Derived intensity features
            (
                {
                    let ratios = &s.ms2_intensity_ratios;
                    ratios
                        .iter()
                        .filter(|r| r.is_finite())
                        .fold(f32::NEG_INFINITY, |a, &b| a.max(b)) as f64
                },
                "ms2_intensity_ratios_max",
            ),
            // Interaction features
            (
                (s.main_score * s.delta_next) as f64,
                "main_times_delta_next",
            ),
            (
                (s.split_product_score * s.fragment_coverage) as f64,
                "split_product_x_coverage",
            ),
            // Summary error features
            (mean_abs_error(&s.ms2_mz_errors), "ms2_mz_mean_abs_error"),
            (
                mean_abs_error(&s.ms2_mobility_errors),
                "ms2_mob_mean_abs_error",
            ),
            (mean_abs_error(&s.ms1_mz_errors), "ms1_mz_mean_abs_error"),
            (
                mean_abs_error(&s.ms1_mobility_errors),
                "ms1_mob_mean_abs_error",
            ),
        ];

        // Sequence-derived block. Gated — all-or-none per speclib load.
        if let Some(counts) = s.peptide.aa_counts() {
            let length = s.peptide.length().unwrap() as f64;
            let n_mods = s.peptide.n_mods().unwrap() as f64;
            v.push((length, "peptide_length"));
            for (i, c) in counts.iter().enumerate() {
                v.push((*c, AA_COUNT_NAMES[i]));
            }
            v.push((n_mods, "peptide_n_mods"));
        }
        v
    }

    pub fn feature_names(&self) -> Vec<&'static str> {
        self.named_features().into_iter().map(|(_, n)| n).collect()
    }
}

impl FeatureLike for CompetedCandidate {
    fn as_feature(&self) -> impl IntoIterator<Item = f64> + '_ {
        self.named_features()
            .into_iter()
            .map(|(v, _)| v)
            .collect::<Vec<_>>()
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

#[cfg(test)]
mod feature_tests {
    use super::*;
    use crate::models::DecoyMarking;
    use crate::models::sequence::{
        AminoAcid,
        ParsedSequence,
        Peptide,
    };
    use crate::scoring::results::{
        CompetedCandidate,
        ScoringFields,
    };
    use smallvec::smallvec;
    use std::sync::Arc;

    fn base_scoring_fields(peptide: Peptide) -> ScoringFields {
        ScoringFields {
            peptide,
            library_id: 0,
            decoy_group_id: 0,
            precursor_mz: 500.0,
            precursor_charge: 2,
            precursor_mobility: 0.9,
            is_target: true,
            library_rt: 60.0,
            calibrated_rt_seconds: 3600.0,
            obs_rt_seconds: 3601.0,
            calibrated_sq_delta_rt: 1.0,
            obs_mobility: 0.91,
            delta_ms1_ms2_mobility: 0.01,
            sq_delta_ms1_ms2_mobility: 0.0001,
            main_score: 10.0,
            delta_next: 2.0,
            delta_second_next: 1.0,
            apex_lazyscore: 5.0,
            ms2_lazyscore: 4.0,
            ms2_isotope_lazyscore: 3.0,
            ms2_isotope_lazyscore_ratio: 0.5,
            lazyscore_z: 2.0,
            lazyscore_vs_baseline: 1.5,
            split_product_score: 0.8,
            cosine_au: 0.7,
            scribe_au: 0.6,
            cosine_cg: 0.5,
            scribe_cg: 0.4,
            cosine_weighted_coelution: 0.9,
            cosine_gradient_consistency: 0.85,
            scribe_weighted_coelution: 0.88,
            scribe_gradient_consistency: 0.82,
            peak_shape: 0.95,
            ratio_cv: 0.1,
            centered_apex: 0.5,
            precursor_coelution: 0.9,
            fragment_coverage: 0.8,
            precursor_apex_match: 0.7,
            xic_quality: 0.75,
            fragment_apex_agreement: 0.85,
            isotope_correlation: 0.9,
            gaussian_correlation: 0.88,
            per_frag_gaussian_corr: 0.87,
            rising_cycles: 3,
            falling_cycles: 2,
            npeaks: 5,
            n_scored_fragments: 6,
            ms2_summed_intensity: 1000.0,
            ms1_summed_intensity: 500.0,
            ms2_mz_errors: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            ms2_mobility_errors: [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07],
            ms1_mz_errors: [0.1, 0.2, 0.3],
            ms1_mobility_errors: [0.01, 0.02, 0.03],
            ms2_intensity_ratios: [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
            ms1_intensity_ratios: [0.9, 0.8, 0.7],
        }
    }

    fn sample_competed_candidate_parsed() -> CompetedCandidate {
        let parsed = ParsedSequence {
            // PEPTIDEK — 8 residues
            residues: smallvec![
                AminoAcid::from_ascii(b'P'),
                AminoAcid::from_ascii(b'E'),
                AminoAcid::from_ascii(b'P'),
                AminoAcid::from_ascii(b'T'),
                AminoAcid::from_ascii(b'I'),
                AminoAcid::from_ascii(b'D'),
                AminoAcid::from_ascii(b'E'),
                AminoAcid::from_ascii(b'K'),
            ],
            mods: smallvec![],
        };
        let peptide = Peptide {
            raw: Arc::from("PEPTIDEK"),
            parsed: Some(parsed),
            decoy: DecoyMarking::Target,
            decoy_group: 0,
        };
        CompetedCandidate {
            scoring: base_scoring_fields(peptide),
            delta_group: 1.0,
            delta_group_ratio: 0.5,
            discriminant_score: 0.0,
            qvalue: 1.0,
        }
    }

    fn sample_competed_candidate_unparsed() -> CompetedCandidate {
        let peptide = Peptide {
            raw: Arc::from("PEPTIDEK"),
            parsed: None,
            decoy: DecoyMarking::Target,
            decoy_group: 0,
        };
        CompetedCandidate {
            scoring: base_scoring_fields(peptide),
            delta_group: 1.0,
            delta_group_ratio: 0.5,
            discriminant_score: 0.0,
            qvalue: 1.0,
        }
    }

    #[test]
    fn features_and_names_same_length() {
        let cand = sample_competed_candidate_parsed();
        let feats: Vec<f64> = cand.as_feature().into_iter().collect();
        let names = cand.feature_names();
        assert_eq!(feats.len(), names.len());
    }

    #[test]
    fn sequence_block_present_when_gate_on() {
        let cand = sample_competed_candidate_parsed();
        let names = cand.feature_names();
        assert!(names.contains(&"peptide_length"));
        assert!(names.contains(&"aa_count_A"));
        assert!(names.contains(&"aa_count_Y"));
        assert!(names.contains(&"peptide_n_mods"));
        assert_eq!(names[names.len() - 22], "peptide_length");
        assert_eq!(names[names.len() - 1], "peptide_n_mods");
    }

    #[test]
    fn sequence_block_absent_when_gate_off() {
        let cand = sample_competed_candidate_unparsed();
        let names = cand.feature_names();
        assert!(!names.contains(&"peptide_length"));
        assert!(!names.contains(&"peptide_n_mods"));
    }

    #[test]
    fn gate_delta_is_22_dims() {
        let on = sample_competed_candidate_parsed().feature_names().len();
        let off = sample_competed_candidate_unparsed().feature_names().len();
        assert_eq!(on - off, 22);
    }

    #[test]
    fn base_feature_count_locked() {
        let off = sample_competed_candidate_unparsed().feature_names().len();
        assert_eq!(
            off, 87,
            "base (gate-off) feature count is locked at 87; update this test if the set intentionally changes"
        );
        let on = sample_competed_candidate_parsed().feature_names().len();
        assert_eq!(on, 109, "gate-on total is 87 base + 22 sequence = 109");
    }
}
