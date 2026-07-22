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
use crate::scoring::blocks::derived::Derived;
use crate::scoring::blocks::result_meta::ResultMeta;
use crate::scoring::blocks::{
    FeatSink,
    NameSink,
    ScoreBlock,
    sequence_counts,
};
use crate::scoring::results::{
    CompetedCandidate,
    FinalResult,
    ScoringFields,
};
use rand::prelude::*;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use std::sync::Arc;
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

/// Fixed shuffle seed used by `rescore`. Makes the pre-rescore shuffle
/// (and therefore the fold assignment + downstream target counts)
/// reproducible across runs, eliminating RNG-driven noise in benches.
/// `GBMConfig::default().seed == 0` already makes the boosting itself
/// deterministic; this seals the only remaining entropy source.
const RESCORE_SHUFFLE_SEED: u64 = 42;

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore(
    mut data: Vec<CompetedCandidate>,
    feature_names: &[Arc<str>],
) -> (Vec<FinalResult>, RescoreFeatureStats) {
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

    let stats = scorer.feature_stats(feature_names);

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

// ---------------------------------------------------------------------------
// CompetedCandidate: FeatureLike + LabelledScore
// ---------------------------------------------------------------------------

/// The set-level ML feature-name list, built ONCE per fit. Names are a property
/// of the feature *set*, not of any record — the value walk (`feature_values`)
/// and this name walk share the same ordered sources, so they align by
/// construction (macro blocks) or by adjacent hand-written pairs.
///
/// The four sources, in order: per-block names
/// ([`ScoringFields::push_feature_names`]), post-model meta ([`ResultMeta`]),
/// cross-field interactions ([`Derived`]), then the conditional sequence block
/// ([`sequence_counts`]) LAST. `gate_on` = the run has parsed sequences
/// (speclib-wide), appending the 22 sequence names.
pub fn feature_name_set(gate_on: bool) -> Vec<Arc<str>> {
    let mut n = NameSink::new();
    ScoringFields::push_feature_names(&mut n);
    <ResultMeta as ScoreBlock>::feature_names(&mut n);
    <Derived as ScoreBlock>::feature_names(&mut n);
    if gate_on {
        sequence_counts::feature_names(&mut n);
    }
    n.into_names()
}

/// The set-level feature names for a batch of candidates. The sequence gate is
/// speclib-wide, so it is read once from the first candidate.
pub fn feature_name_set_for(candidates: &[CompetedCandidate]) -> Vec<Arc<str>> {
    let gate_on = candidates
        .first()
        .map(|c| c.scoring.identity.peptide.aa_counts().is_some())
        .unwrap_or(false);
    feature_name_set(gate_on)
}

impl CompetedCandidate {
    /// This record's ML feature *values*, in the set-level name order (see
    /// [`feature_name_set`]). Names are not carried per record.
    ///
    /// Assembled from the same four ordered sources as the name walk: per-block
    /// values, post-model meta, cross-field derived, then the conditional
    /// sequence block LAST.
    fn feature_values(&self) -> Vec<f64> {
        let mut sink = FeatSink::new();
        self.scoring.push_features(&mut sink);
        self.result_meta().features(&mut sink);
        Derived::compute(&self.scoring).features(&mut sink);
        sequence_counts::features(&self.scoring.identity.peptide, &mut sink);
        sink.into_values()
    }
}

impl FeatureLike for CompetedCandidate {
    fn as_feature(&self) -> impl IntoIterator<Item = f64> + '_ {
        self.feature_values()
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

    #[test]
    fn no_duplicate_feature_names() {
        use std::collections::BTreeSet;
        let names = feature_name_set(true);
        let uniq: BTreeSet<&str> = names.iter().map(|n| n.as_ref()).collect();
        assert_eq!(names.len(), uniq.len(), "duplicate feature name emitted");
    }

    #[test]
    fn features_and_names_same_length() {
        // The load-bearing invariant: a record's value vector aligns 1:1 with
        // the set-level name list built independently. Both gates checked.
        let on: Vec<f64> = sample_competed_candidate_parsed()
            .as_feature()
            .into_iter()
            .collect();
        assert_eq!(on.len(), feature_name_set(true).len());
        let off: Vec<f64> = sample_competed_candidate_unparsed()
            .as_feature()
            .into_iter()
            .collect();
        assert_eq!(off.len(), feature_name_set(false).len());
    }

    #[test]
    fn neutralize_mobility_nans_every_mobility_feature() {
        // For a run with no searchable mobility axis (mzML/FAIMS), all
        // mobility-derived GBM features must become NaN (forust missing), so
        // they cannot bias the score with sentinel-derived constants.
        let mut cand = sample_competed_candidate_parsed();
        cand.scoring.neutralize_mobility();

        let names = feature_name_set(true);
        let vals: Vec<f64> = cand.as_feature().into_iter().collect();
        assert_eq!(names.len(), vals.len());
        let mob: Vec<(&Arc<str>, &f64)> = names
            .iter()
            .zip(vals.iter())
            .filter(|(n, _)| n.contains("mob"))
            .collect();
        assert_eq!(
            mob.len(),
            16,
            "mobility feature count changed: {:?}",
            mob.iter().map(|(n, _)| n.as_ref()).collect::<Vec<_>>()
        );
        for (n, v) in &mob {
            assert!(v.is_nan(), "mobility feature {n} should be NaN, got {v}");
        }
        // Non-mobility features are untouched (at least one stays finite).
        assert!(
            names
                .iter()
                .zip(vals.iter())
                .any(|(n, v)| !n.contains("mob") && v.is_finite()),
            "non-mobility features must remain finite"
        );
    }

    #[test]
    fn sequence_block_present_when_gate_on() {
        let names = feature_name_set(true);
        let has = |t: &str| names.iter().any(|n| n.as_ref() == t);
        assert!(has("peptide_length"));
        assert!(has("aa_count_A"));
        assert!(has("aa_count_Y"));
        assert!(has("peptide_n_mods"));
        assert_eq!(&*names[names.len() - 22], "peptide_length");
        assert_eq!(&*names[names.len() - 1], "peptide_n_mods");
    }

    #[test]
    fn sequence_block_absent_when_gate_off() {
        let names = feature_name_set(false);
        assert!(!names.iter().any(|n| n.as_ref() == "peptide_length"));
        assert!(!names.iter().any(|n| n.as_ref() == "peptide_n_mods"));
    }

    #[test]
    fn gate_delta_is_22_dims() {
        let on = sample_competed_candidate_parsed()
            .as_feature()
            .into_iter()
            .count();
        let off = sample_competed_candidate_unparsed()
            .as_feature()
            .into_iter()
            .count();
        assert_eq!(on - off, 22);
    }
}
