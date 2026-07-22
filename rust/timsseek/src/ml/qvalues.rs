use super::cv::{
    CrossValidatedScorer,
    DataBuffer,
    FeatureLike,
    GBMConfig,
    RescoreFeatureStats,
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

pub fn rescore(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
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
    // NOTE: `library_id` is now the POSITIONAL target index (lazy arena) /
    // materialized eg id — a target and its ±decoy variants share it, but
    // target-decoy competition leaves exactly ONE candidate per
    // (target, charge), so the key stays unique among survivors here.
    data.sort_unstable_by_key(|c| {
        (
            c.scoring.identity.library_id,
            c.scoring.identity.precursor_charge,
        )
    });

    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(RESCORE_SHUFFLE_SEED);
    data.shuffle(&mut rng);

    let mut scorer = CrossValidatedScorer::<CompetedCandidate>::new_from_shuffled(3, data, config);
    scorer
        .fit(&mut DataBuffer::default(), &mut DataBuffer::default())
        .unwrap();

    let names: Vec<&'static str> = scorer
        .data()
        .first()
        .map(|c| c.feature_names())
        .unwrap_or_default();
    let stats = scorer.feature_stats(&names);

    let mut scored = scorer.score();
    // Sort by score descending
    #[cfg(feature = "rayon")]
    scored.par_sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    #[cfg(not(feature = "rayon"))]
    scored.sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut scored, |x| CompetedCandidate::get_score(x) as f32);
    debug!("Best:\n{:#?}", scored.first());
    debug!("Worst:\n{:#?}", scored.last());

    (scored.into_iter().map(|c| c.into_final()).collect(), stats)
}

use crate::scoring::blocks::{
    FeatSink,
    ScoreBlock,
    derived,
    sequence_counts,
};
use crate::scoring::results::{
    CompetedCandidate,
    FinalResult,
};

// ---------------------------------------------------------------------------
// CompetedCandidate: FeatureLike + LabelledScore
// ---------------------------------------------------------------------------

impl CompetedCandidate {
    /// Single source of truth for feature values and names.
    ///
    /// Assembled by walking each block's `features()` (codegen-derived), then
    /// the post-model meta block, then the cross-field derived features, then
    /// the conditional sequence block LAST.
    ///
    /// Base set: 86 dims always present. Sequence block: 22 dims appended when
    /// `peptide.aa_counts()` is `Some`. The gate is speclib-wide, so the vector
    /// length is stable within a single fit.
    fn named_features(&self) -> Vec<(f64, &'static str)> {
        // Single assembly point for the ML feature vector. The full set is
        // spread across four sources: per-block projections
        // (`ScoringFields::push_features`), post-model meta (`result_meta`),
        // cross-field interactions (`scoring::blocks::derived`), and the
        // conditional sequence block (`scoring::blocks::sequence_counts`). The
        // `feature_name_set_matches_golden` test is what pins the resulting set.
        let mut sink = FeatSink::new();
        self.scoring.push_features(&mut sink);
        self.result_meta().features(&mut sink);
        derived::features(&self.scoring, &mut sink);
        sequence_counts::features(&self.scoring.identity.peptide, &mut sink);
        sink.into_named()
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
        if self.scoring.identity.is_target {
            1.0
        } else {
            0.0
        }
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
        if self.scoring.identity.is_target {
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
        if self.scoring.identity.is_target {
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
        ScoringFields::sample(peptide)
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

    /// The 86 base feature names (gate off). This is the GBM contract: the
    /// *set* must stay stable across refactors, so a rename that preserves the
    /// count (which `base_feature_count_locked` alone would miss) is caught here.
    const GOLDEN_BASE_FEATURES: &[&str] = &[
        "precursor_mz_round5",
        "precursor_charge",
        "precursor_mobility",
        "calibrated_rt_seconds_round",
        "calibrated_rt_seconds",
        "obs_rt_seconds",
        "calibrated_sq_delta_rt",
        "obs_mobility",
        "delta_ms1_ms2_mobility",
        "sq_delta_ms1_ms2_mobility",
        "main_score",
        "delta_next",
        "delta_second_next",
        "split_product_score_ln1p",
        "cosine_au_ln1p",
        "scribe_au_ln1p",
        "cosine_cg",
        "scribe_cg",
        "cosine_weighted_coelution",
        "cosine_gradient_consistency",
        "scribe_weighted_coelution",
        "scribe_gradient_consistency",
        "peak_shape",
        "ratio_cv",
        "centered_apex",
        "precursor_coelution",
        "fragment_coverage",
        "precursor_apex_match",
        "xic_quality",
        "fragment_apex_agreement",
        "isotope_correlation",
        "gaussian_correlation",
        "per_frag_gaussian_corr",
        "apex_lazyscore",
        "lazyscore_z",
        "lazyscore_vs_baseline",
        "ms2_lazyscore",
        "ms2_isotope_lazyscore",
        "ms2_isotope_lazyscore_ratio",
        "rising_cycles",
        "falling_cycles",
        "npeaks",
        "n_scored_fragments",
        "ms2_summed_intensity_ln1p",
        "ms1_summed_intensity_ln1p",
        "ms2_mz_error_0",
        "ms2_mz_error_1",
        "ms2_mz_error_2",
        "ms2_mz_error_3",
        "ms2_mz_error_4",
        "ms2_mz_error_5",
        "ms2_mz_error_6",
        "ms2_mobility_error_0",
        "ms2_mobility_error_1",
        "ms2_mobility_error_2",
        "ms2_mobility_error_3",
        "ms2_mobility_error_4",
        "ms2_mobility_error_5",
        "ms2_mobility_error_6",
        "ms1_mz_error_0",
        "ms1_mz_error_1",
        "ms1_mz_error_2",
        "ms1_mobility_error_0",
        "ms1_mobility_error_1",
        "ms1_mobility_error_2",
        "ms1_intensity_ratio_0",
        "ms1_intensity_ratio_1",
        "ms1_intensity_ratio_2",
        "ms2_intensity_ratio_0",
        "ms2_intensity_ratio_1",
        "ms2_intensity_ratio_2",
        "ms2_intensity_ratio_3",
        "ms2_intensity_ratio_4",
        "ms2_intensity_ratio_5",
        "ms2_intensity_ratio_6",
        "delta_group",
        "delta_group_ratio",
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

    /// The 22 sequence feature names appended when the gate is on.
    const GOLDEN_SEQUENCE_FEATURES: &[&str] = &[
        "peptide_length",
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
        "peptide_n_mods",
    ];

    #[test]
    fn feature_name_set_matches_golden() {
        use std::collections::BTreeSet;
        let off: BTreeSet<&str> = sample_competed_candidate_unparsed()
            .feature_names()
            .into_iter()
            .collect();
        let golden_base: BTreeSet<&str> = GOLDEN_BASE_FEATURES.iter().copied().collect();
        assert_eq!(
            off, golden_base,
            "base feature name-set drifted; update GOLDEN_BASE_FEATURES only if the change is intended (and bump the workspace version + retrain GBM)"
        );

        let on: BTreeSet<&str> = sample_competed_candidate_parsed()
            .feature_names()
            .into_iter()
            .collect();
        let golden_on: BTreeSet<&str> = GOLDEN_BASE_FEATURES
            .iter()
            .chain(GOLDEN_SEQUENCE_FEATURES.iter())
            .copied()
            .collect();
        assert_eq!(
            on, golden_on,
            "gate-on feature name-set drifted; see GOLDEN_SEQUENCE_FEATURES"
        );
    }

    #[test]
    fn no_duplicate_feature_names() {
        use std::collections::BTreeSet;
        let names = sample_competed_candidate_parsed().feature_names();
        let uniq: BTreeSet<&str> = names.iter().copied().collect();
        assert_eq!(names.len(), uniq.len(), "duplicate feature name emitted");
    }

    #[test]
    fn features_and_names_same_length() {
        let cand = sample_competed_candidate_parsed();
        let feats: Vec<f64> = cand.as_feature().into_iter().collect();
        let names = cand.feature_names();
        assert_eq!(feats.len(), names.len());
    }

    #[test]
    fn neutralize_mobility_nans_every_mobility_feature() {
        // For a run with no searchable mobility axis (mzML/FAIMS), all
        // mobility-derived GBM features must become NaN (forust missing), so
        // they cannot bias the score with sentinel-derived constants.
        let mut cand = sample_competed_candidate_parsed();
        cand.scoring.neutralize_mobility();

        let feats = cand.named_features();
        let mob: Vec<_> = feats.iter().filter(|(_, n)| n.contains("mob")).collect();
        assert_eq!(
            mob.len(),
            16,
            "mobility feature count changed: {:?}",
            mob.iter().map(|(_, n)| *n).collect::<Vec<_>>()
        );
        for (v, n) in &mob {
            assert!(v.is_nan(), "mobility feature {n} should be NaN, got {v}");
        }
        // Non-mobility features are untouched (at least one stays finite).
        assert!(
            feats
                .iter()
                .any(|(v, n)| !n.contains("mob") && v.is_finite()),
            "non-mobility features must remain finite"
        );
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
            off, 86,
            "base (gate-off) feature count is locked at 86; update this test if the set intentionally changes"
        );
        let on = sample_competed_candidate_parsed().feature_names().len();
        assert_eq!(on, 108, "gate-on total is 86 base + 22 sequence = 108");
    }
}
