use super::cv::{
    CrossValidatedScorer,
    DataBuffer,
    FeatureLike,
    FeatureStat,
    FoldStats,
    GBMConfig,
    RescoreFeatureStats,
};
use super::lda::{
    DEFAULT_LDA_SHRINKAGE,
    LdaModel,
    inverse_normal_transform_columns,
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
    let config = GBMConfig::from_env();

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

/// Sage-style shrinkage-LDA rescorer. Single closed-form linear fit over all
/// candidates (no CV, no boosting): ~100x cheaper than the GBM path. The FDR
/// machinery (`assign_qval`, target-decoy competition) is untouched — only the
/// discriminant score source changes.
///
/// Selected at runtime via `TIMSSEEK_RESCORE_MODEL=lda`; the GBM `rescore`
/// remains the default. See `ml::lda` for the fit details.
pub fn rescore_lda(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    use std::time::Instant;

    // Same canonical ordering as the GBM path for run-to-run determinism.
    data.sort_unstable_by_key(|c| {
        (
            c.scoring.identity.library_id,
            c.scoring.identity.precursor_charge,
        )
    });

    let base_names: Vec<Arc<str>> = feature_name_set_for(&data);
    let base_ncols = base_names.len();
    let nrows = data.len();

    // Materialize the base feature matrix once as a single flat row-major buffer
    // (`feat[i*ncols + j]`) — a `Vec<Vec<f64>>` here means one heap allocation
    // per candidate, which dominates at tens of millions of rows. Track finite
    // counts per column so we can encode informative missingness below.
    let t = Instant::now();
    let mut base: Vec<f64> = Vec::with_capacity(nrows * base_ncols);
    for c in data.iter() {
        base.extend(c.as_feature());
    }
    debug_assert_eq!(base.len(), nrows * base_ncols);
    let is_decoy: Vec<bool> = data.iter().map(|c| c.get_y() < 0.5).collect();

    // Optional raw-matrix dump for offline feature-engineering iteration.
    // `TIMSSEEK_LDA_DUMP=/prefix` writes `<prefix>.f64` (row-major matrix),
    // `<prefix>.labels` (u8, 1=target), `<prefix>.names.txt`. Pre-transform.
    if let Ok(prefix) = std::env::var("TIMSSEEK_LDA_DUMP") {
        dump_feature_matrix(&prefix, &base, nrows, &base_names, &is_decoy);
    }

    // Center-is-better -> monotone-is-better. Signed mass/mobility/RT errors are
    // best near zero and bad in *either* tail. A linear discriminant assigns one
    // sign per feature, so it cannot express "close to zero is good" and drives
    // these ~22 high-value features to ~0 weight (GBM splits on |err| instead).
    // Fold each to its magnitude so small error becomes a monotone extreme.
    // Default on; `TIMSSEEK_LDA_ABS=none` disables.
    let use_abs = std::env::var("TIMSSEEK_LDA_ABS")
        .map(|v| !v.eq_ignore_ascii_case("none"))
        .unwrap_or(true);
    if use_abs {
        absolutize_center_better(&mut base, base_ncols, &base_names);
    }

    // Missingness indicators. GBM handles NaN natively (a missing branch), and
    // absence of higher-index fragment ions correlates with target/decoy. LDA
    // imputes NaN -> mean, discarding that signal, so we hand it back as binary
    // `<feat>_isna` columns for every feature whose missing rate is informative
    // (neither ~always-present nor ~always-absent). Default on;
    // `TIMSSEEK_LDA_INDICATORS=none` disables.
    let use_indicators = std::env::var("TIMSSEEK_LDA_INDICATORS")
        .map(|v| !v.eq_ignore_ascii_case("none"))
        .unwrap_or(true);
    let (feat, names, ncols) = if use_indicators {
        let (feat, names) = augment_with_missingness(&base, nrows, &base_names);
        let ncols = names.len();
        (feat, names, ncols)
    } else {
        (base, base_names, base_ncols)
    };
    eprintln!(
        "  LDA: extracted {nrows} x {ncols} feature matrix ({} missingness cols) in {:.2?}",
        ncols - base_ncols,
        t.elapsed()
    );
    let mut feat = feat;

    // Gaussianize each feature toward the LDA normality assumption. Default on;
    // `TIMSSEEK_LDA_TRANSFORM=none` disables for A/B comparison.
    let use_int = std::env::var("TIMSSEEK_LDA_TRANSFORM")
        .map(|v| !v.eq_ignore_ascii_case("none"))
        .unwrap_or(true);
    if use_int {
        let t = Instant::now();
        inverse_normal_transform_columns(&mut feat, nrows, ncols);
        eprintln!(
            "  LDA: inverse-normal transform ({ncols} cols) in {:.2?}",
            t.elapsed()
        );
    }

    let stats = match LdaModel::fit(&feat, nrows, ncols, &is_decoy, DEFAULT_LDA_SHRINKAGE) {
        Some(model) => {
            let t = Instant::now();
            let mut scores = vec![0.0f64; nrows];
            model.score_all(&feat, nrows, &mut scores);
            for (cand, &s) in data.iter_mut().zip(scores.iter()) {
                cand.assign_score(s);
            }
            eprintln!("  LDA: scored {nrows} candidates in {:.2?}", t.elapsed());
            lda_feature_stats(&names, &feat, nrows, model.coef())
        }
        None => {
            tracing::error!("LDA fit failed (singular or empty class); scores left at zero");
            vec![FoldStats {
                fold: 0,
                feature_stats: Vec::new(),
                feature_importance: Vec::new(),
            }]
        }
    };

    let mut scored = data;
    #[cfg(feature = "rayon")]
    scored.par_sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    #[cfg(not(feature = "rayon"))]
    scored.sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut scored, |x| CompetedCandidate::get_score(x) as f32);

    (scored.into_iter().map(|c| c.into_final()).collect(), stats)
}

/// Dump the raw feature matrix + labels for offline analysis. Best-effort:
/// logs and returns on any I/O error rather than aborting the run.
fn dump_feature_matrix(
    prefix: &str,
    base: &[f64],
    nrows: usize,
    names: &[Arc<str>],
    is_decoy: &[bool],
) {
    use std::io::Write;
    let write = || -> std::io::Result<()> {
        let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{prefix}.f64"))?);
        // Header: nrows, ncols as u64 little-endian, then row-major f64.
        f.write_all(&(nrows as u64).to_le_bytes())?;
        f.write_all(&(names.len() as u64).to_le_bytes())?;
        // SAFETY: transmuting &[f64] to &[u8] for a bulk write.
        let bytes =
            unsafe { std::slice::from_raw_parts(base.as_ptr() as *const u8, base.len() * 8) };
        f.write_all(bytes)?;
        f.flush()?;

        let labels: Vec<u8> = is_decoy.iter().map(|&d| u8::from(!d)).collect();
        std::fs::write(format!("{prefix}.labels"), &labels)?;
        let names_txt = names.iter().map(|s| s.as_ref()).collect::<Vec<_>>().join("\n");
        std::fs::write(format!("{prefix}.names.txt"), names_txt)?;
        Ok(())
    };
    match write() {
        Ok(()) => eprintln!("  LDA: dumped feature matrix to {prefix}.{{f64,labels,names.txt}}"),
        Err(e) => tracing::error!("feature dump failed: {e}"),
    }
}

/// Predicate: does this feature's discriminative structure put the *best*
/// value at zero, with both tails bad? Such signed features are useless to a
/// single-sign linear discriminant until folded to their magnitude.
fn is_center_better(name: &str) -> bool {
    // Signed per-ion / precursor mass + mobility errors, RT residual, and the
    // MS1-vs-MS2 mobility offset. Excludes already-magnitude features
    // (`*_sq_*`, `*_mean_abs_error`, `calibrated_sq_delta_rt`).
    name == "rt_err"
        || name == "delta_ms1_ms2_mobility"
        || name.starts_with("ms2_mz_err_")
        || name.starts_with("ms2_mob_err_")
        || name.starts_with("ms1_mz_err_")
        || name.starts_with("ms1_mob_err_")
}

/// Replace each center-is-better column with its magnitude (`|x|`, NaN
/// preserved) so a linear model can weight it monotonically.
fn absolutize_center_better(base: &mut [f64], base_ncols: usize, base_names: &[Arc<str>]) {
    let abs_cols: Vec<usize> = (0..base_ncols)
        .filter(|&j| is_center_better(&base_names[j]))
        .collect();
    if abs_cols.is_empty() {
        return;
    }
    for row in base.chunks_exact_mut(base_ncols) {
        for &j in &abs_cols {
            row[j] = row[j].abs();
        }
    }
}

/// Fraction of rows that must be missing (and present) for a feature's
/// missingness to earn its own indicator column. Below this it carries no
/// signal; above `1 - this` it is ~always absent (equally useless).
const MISSINGNESS_MIN_RATE: f64 = 0.01;

/// Append a binary `<name>_isna` column (1.0 = the source value was non-finite)
/// for every base feature whose missing rate lies in
/// `[MISSINGNESS_MIN_RATE, 1 - MISSINGNESS_MIN_RATE]`. Returns the augmented
/// row-major matrix and the extended name list.
fn augment_with_missingness(
    base: &[f64],
    nrows: usize,
    base_names: &[Arc<str>],
) -> (Vec<f64>, Vec<Arc<str>>) {
    let base_ncols = base_names.len();
    // Per-column NaN counts.
    let mut nan = vec![0u64; base_ncols];
    for row in base.chunks_exact(base_ncols) {
        for (j, &v) in row.iter().enumerate() {
            if !v.is_finite() {
                nan[j] += 1;
            }
        }
    }
    let n = nrows.max(1) as f64;
    let sel: Vec<usize> = (0..base_ncols)
        .filter(|&j| {
            let r = nan[j] as f64 / n;
            r >= MISSINGNESS_MIN_RATE && r <= 1.0 - MISSINGNESS_MIN_RATE
        })
        .collect();
    if sel.is_empty() {
        return (base.to_vec(), base_names.to_vec());
    }

    let ncols = base_ncols + sel.len();
    let mut names: Vec<Arc<str>> = base_names.to_vec();
    for &j in &sel {
        names.push(Arc::from(format!("{}_isna", base_names[j])));
    }

    let mut out = vec![0.0f64; nrows * ncols];
    for (i, row) in base.chunks_exact(base_ncols).enumerate() {
        let obase = i * ncols;
        out[obase..obase + base_ncols].copy_from_slice(row);
        for (t, &j) in sel.iter().enumerate() {
            out[obase + base_ncols + t] = if row[j].is_finite() { 0.0 } else { 1.0 };
        }
    }
    (out, names)
}

/// Single-"fold" feature stats for the LDA path: per-feature finite means +
/// NaN ratios, plus `|coef|` as the importance ranking (LDA weights in
/// standardized space are directly interpretable as importance).
fn lda_feature_stats(
    names: &[Arc<str>],
    feat: &[f64],
    nrows: usize,
    coef: &[f64],
) -> RescoreFeatureStats {
    let ncols = names.len();
    let mut sums = vec![0.0f64; ncols];
    let mut finite = vec![0u32; ncols];
    let mut nan = vec![0u32; ncols];
    for row in feat.chunks_exact(ncols) {
        for j in 0..ncols {
            let v = row[j];
            if v.is_finite() {
                sums[j] += v;
                finite[j] += 1;
            } else {
                nan[j] += 1;
            }
        }
    }
    let n = nrows.max(1) as f32;
    let feature_stats: Vec<FeatureStat> = (0..ncols)
        .map(|j| FeatureStat {
            name: names[j].clone(),
            mean: if finite[j] > 0 {
                (sums[j] / finite[j] as f64) as f32
            } else {
                f32::NAN
            },
            nan_ratio: nan[j] as f32 / n,
        })
        .collect();

    let mut feature_importance: Vec<(Arc<str>, f32)> = names
        .iter()
        .zip(coef.iter())
        .map(|(nm, c)| (nm.clone(), c.abs() as f32))
        .collect();
    feature_importance.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    vec![FoldStats {
        fold: 0,
        feature_stats,
        feature_importance,
    }]
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
