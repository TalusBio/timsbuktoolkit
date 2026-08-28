use super::cv::{
    CrossValidatedScorer,
    FeatureLike,
    FoldDataset,
    FoldModel,
    GBMConfig,
    GbmFoldModel,
    PrecomputedFeatures,
    RescoreFeatureStats,
    StreamingDataset,
    fold_feature_stats,
};
use super::lda::{
    LdaConfig,
    LdaModel,
};
use super::mlp::MlpConfig;
use super::mlp_fold::MlpFoldModel;
use super::{
    LabelledScore,
    N_RESCORE_FOLDS,
    TargetDecoy,
};
use crate::scoring::blocks::derived::Derived;
use crate::scoring::blocks::result_meta::ResultMeta;
use crate::scoring::blocks::{
    NameSink,
    ScoreBlock,
    sequence_counts,
};
use crate::scoring::results::{
    CompetedCandidate,
    FinalResult,
    ScoringFields,
};
use crate::utils::maybe_par;
use rand::prelude::*;
use std::sync::Arc;

#[cfg(test)]
use super::cv::RowMajorDataset;
use tracing::debug;

/// Failure to produce a valid discriminant score for every rescore fold.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RescoreError {
    /// A cross-fit fold failed while fitting or predicting.
    CrossFit {
        model: &'static str,
        fold: usize,
        folds: usize,
        stage: &'static str,
        reason: String,
    },
    /// A model failed outside the dedicated cross-fit driver.
    Model { model: &'static str, reason: String },
    /// One or more folds completed without producing a usable model.
    UntrainedFolds { model: &'static str, folds: Vec<u8> },
}

impl std::fmt::Display for RescoreError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CrossFit {
                model,
                fold,
                folds,
                stage,
                reason,
            } => write!(
                f,
                "cross-fit {model} {stage} failed on fold {}/{}: {reason}",
                fold + 1,
                folds
            ),
            Self::Model { model, reason } => write!(f, "{model} rescore failed: {reason}"),
            Self::UntrainedFolds { model, folds } => {
                write!(f, "{model} produced no trained model for folds {folds:?}")
            }
        }
    }
}

impl std::error::Error for RescoreError {}

/// Scored results and their per-fold feature statistics.
pub type RescoreResult = Result<(Vec<FinalResult>, RescoreFeatureStats), RescoreError>;

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

/// How many targets and decoys pass at one q-value cutoff.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThresholdCounts {
    pub q: f32,
    pub targets: usize,
    pub decoys: usize,
}

/// Count the targets and decoys at or below each q-value threshold, one
/// [`ThresholdCounts`] per threshold in the order given.
pub fn report_qvalues_at_thresholds<T: LabelledScore + std::fmt::Debug>(
    scores: &[T],
    thresholds: &[f32],
) -> Vec<ThresholdCounts> {
    let mut out = Vec::new();

    for &q in thresholds {
        let targets = scores
            .iter()
            .filter(|s| s.get_qval() <= q && matches!(s.get_label(), TargetDecoy::Target))
            .count();
        let decoys = scores
            .iter()
            .filter(|s| s.get_qval() <= q && matches!(s.get_label(), TargetDecoy::Decoy))
            .count();
        out.push(ThresholdCounts { q, targets, decoys });
    }

    out
}

/// Fixed shuffle seed used by `rescore`. Makes the pre-rescore shuffle
/// (and therefore the fold assignment + downstream target counts)
/// reproducible across runs, eliminating RNG-driven noise in benches.
/// `GBMConfig::default().seed == 0` already makes the boosting itself
/// deterministic; this seals the only remaining entropy source.
const RESCORE_SHUFFLE_SEED: u64 = 42;

/// Canonicalize candidate order, then apply the fixed seeded shuffle.
///
/// Shared by all four rescorers so the row order — and therefore the positional
/// fold partition `fold(i) = i % N_RESCORE_FOLDS` derived from it — has exactly
/// one definition.
///
/// Canonicalizing first matters because upstream stages can emit candidates in
/// an order that drifts with floating-point accumulation quirks (e.g. different
/// peak-bucket layouts produce identical features but different vec orderings).
/// Without a stable sort here, the seeded shuffle sees different inputs across
/// runs -> different fold assignment -> different q-values -> drifting target
/// counts across equivalent configs. `(library_id, precursor_charge)` is a
/// non-FP composite key that should be unique per candidate after target-decoy
/// competition.
///
/// NOTE: `library_id` is the POSITIONAL target index (lazy arena) / materialized
/// eg id — a target and its ±decoy variants share it, but target-decoy
/// competition leaves exactly ONE candidate per (target, charge), so the key
/// stays unique among survivors here.
///
/// The shuffle is NOT cosmetic: every fold partition below is positional, so an
/// unshuffled (library-ordered) input would make folds correlate with library
/// layout — e.g. protein-grouped or target/decoy-blocked speclibs would give
/// systematically unbalanced folds.
///
/// Feature rows must be projected after this call because fold assignment and
/// labels are positional. Building a frame first would misalign its rows after
/// the shuffle without changing its shape.
fn canonicalize_and_shuffle(data: &mut [CompetedCandidate]) {
    data.sort_unstable_by_key(|c| (c.scoring.identity.row, c.scoring.identity.precursor_charge));

    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(RESCORE_SHUFFLE_SEED);
    data.shuffle(&mut rng);
}

/// The shared tail of all four rescorers: rank by discriminant score, derive
/// q-values from that ranking, and convert to [`FinalResult`].
///
/// Rank direction is load-bearing and must stay DESCENDING: `assign_qval`
/// asserts it, and every q-value is a running target/decoy ratio over exactly
/// this walk, so the key here and the key handed to `assign_qval` have to be
/// the same function of a candidate.
///
/// Tied scores change the q-values they produce, so the sort is stable: ties
/// keep the seeded-shuffle input order and the result is reproducible on any
/// thread count. Reproducible, not correct — breaking ties on a real key is a
/// separate decision, not made here.
fn finalize(
    mut scored: Vec<CompetedCandidate>,
    stats: RescoreFeatureStats,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    maybe_par::sort_by(&mut scored, |a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut scored, |x| CompetedCandidate::get_score(x) as f32);
    debug!("Best:\n{:#?}", scored.first());
    debug!("Worst:\n{:#?}", scored.last());

    (scored.into_iter().map(|c| c.into_final()).collect(), stats)
}

/// A cross-fit (leak-free) model over the canonical rescore fold partition.
///
/// See [`crossfit`] for the partition contract.
struct CrossFit<M: FoldModel> {
    /// Held-out score per row, row-aligned with the input matrix.
    scores: Vec<f64>,
    /// The per-fold fitted models, indexed by fold. Always
    /// `N_RESCORE_FOLDS` long — a `CrossFit` only exists if every fold fit.
    models: Vec<M>,
    /// `fold_rows[f]` = the rows `models[f]` SCORED (ascending), i.e. exactly
    /// fold `f`. Kept so the sidecar stats summarize exactly those rows without
    /// re-deriving the partition a second time.
    fold_rows: Vec<Vec<usize>>,
}

impl<M: FoldModel> CrossFit<M> {
    /// Per-fold sidecar stats, on the GBM path's shared implementation
    /// ([`fold_feature_stats`]): for fold `f` the means / NaN ratios cover the
    /// rows `models[f]` scored, and the importance is that model's own
    /// per-column measurement.
    ///
    /// `data` must be the dataset this was cross-fit over.
    fn feature_stats<D: FoldDataset>(&self, data: &D) -> RescoreFeatureStats {
        let models: Vec<Option<&M>> = self.models.iter().map(Some).collect();
        fold_feature_stats(data, &self.fold_rows, &models)
    }
}

/// Cross-fit a [`FoldModel`] over the canonical rescore fold partition and
/// return each row's HELD-OUT score.
///
/// The shared statement of leak-freedom for [`crossfit_lda`], `rescore_lda`,
/// and `rescore_hybrid`.
///
/// Generic over the model because the partition is independent of the fitted
/// model. `what` names the model in failure logs.
///
/// # Why a row may never be scored by a model that saw it
/// These models are label-aware: they fit on the target/decoy labels. One
/// in-sample fit over all rows, scoring those same rows, lets every row's
/// discriminant peek at its own label; the target/decoy separation then looks
/// better than it is, and `assign_qval` derives the q-values from exactly that
/// separation. The result is an FDR that is wrong in the flattering direction
/// with nothing downstream to catch it.
///
/// The hybrid needs this even more sharply: there the held-out score is not the
/// final score but one `lda_score` column fed to a
/// cross-validated GBM. An in-sample column would smuggle a row's own label
/// into a feature the GBM reads while that row is held out of its GBM fold — so
/// the GBM's own CV cannot notice, and the leak surfaces only as an optimistic
/// FDR.
///
/// # Partition
/// Fold membership comes from [`FoldDataset::get_fold`] — for the
/// production datasets, that is `i % N_RESCORE_FOLDS`. For
/// fold `f` the model is fitted on every row with `get_fold(i) != f` and then
/// scores only the rows with `get_fold(i) == f`, so no row contributes to the
/// model scoring it.
///
/// This is NOT `CrossValidatedScorer`'s partition and must not be unified with
/// it: that one fits on fold `f` alone, early-stops on `f + 1`, and scores the
/// remaining `n_folds - 2` folds, so a row is scored by several models that
/// each saw a fraction of the data. Here a row is scored exactly once, by a
/// model that saw everything else. Both satisfy leak-freedom, which is the only
/// property either has to satisfy; how many rows a model sees and how many
/// models score a row are free to differ, and do.
///
/// What the two DO have to agree on is the fold ASSIGNMENT `get_fold`, or a
/// hybrid row's cross-fit column can come from a model trained on rows the GBM
/// is holding out — leak restored, silently. Both sides read that assignment
/// from [`FoldDataset::get_fold`], so it has one definition; the
/// `crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on` test pins the
/// two partitions against each other rather than against a re-typed modulo.
///
/// Returns an error if any fold cannot fit or score all of its held-out rows.
/// Callers therefore cannot accidentally consume a partially filled score
/// column.
fn crossfit<D: FoldDataset, M: FoldModel>(
    data: &D,
    cfg: &M::Config,
    what: &'static str,
) -> Result<CrossFit<M>, RescoreError>
where
    M::Error: std::fmt::Display,
{
    let n_folds = N_RESCORE_FOLDS as usize;
    let nrows = data.nrows();
    let mut scores = vec![0.0f64; nrows];
    let mut models = Vec::with_capacity(n_folds);
    let mut fold_rows: Vec<Vec<usize>> = Vec::with_capacity(n_folds);
    debug_assert!(
        (0..nrows).all(|i| data.get_fold(i) < n_folds),
        "dataset partitions into more folds than the cross-fit walks"
    );

    for f in 0..n_folds {
        // TRAIN rows = every row NOT in fold f. HELD = exactly fold f.
        // Ascending in both cases, so the row order the fit reduces over is the
        // dataset's own order.
        let train: Vec<usize> = (0..nrows).filter(|&i| data.get_fold(i) != f).collect();
        let held: Vec<usize> = (0..nrows).filter(|&i| data.get_fold(i) == f).collect();

        // `val` is EMPTY and the partition is why: this walk trains on
        // all-but-fold-`f` and scores fold `f`, so every row is already spoken
        // for and there is no third slice to hand over. The LDA is closed-form
        // and ignores it. `MlpFoldModel` DOES early-stop, and handles the empty
        // slice by carving a deterministic inner validation set out of `train`
        // — see its `fit` for the rule; nothing about it reaches fold `f`.
        let model = M::fit(cfg, data, f, &train, &[]).map_err(|e| RescoreError::CrossFit {
            model: what,
            fold: f,
            folds: n_folds,
            stage: "fit",
            reason: e.to_string(),
        })?;
        // Score ONLY the held-out fold, with a model that never saw it.
        let preds = model
            .predict(data, &held)
            .map_err(|e| RescoreError::CrossFit {
                model: what,
                fold: f,
                folds: n_folds,
                stage: "predict",
                reason: e.to_string(),
            })?;
        if preds.len() != held.len() {
            return Err(RescoreError::CrossFit {
                model: what,
                fold: f,
                folds: n_folds,
                stage: "predict",
                reason: format!(
                    "returned {} scores for {} held-out rows",
                    preds.len(),
                    held.len()
                ),
            });
        }
        for (&i, s) in held.iter().zip(preds) {
            scores[i] = s;
        }
        fold_rows.push(held);
        models.push(model);
    }

    Ok(CrossFit {
        scores,
        models,
        fold_rows,
    })
}

/// Cross-fit an [`LdaModel`] on [`LdaConfig::default`] — see [`crossfit`] for
/// the partition and the failure policy.
///
/// `N_RESCORE_FOLDS` fits instead of 1 is affordable precisely because the LDA
/// is closed-form and cheap.
fn crossfit_lda<D: FoldDataset>(data: &D) -> Result<CrossFit<LdaModel>, RescoreError> {
    crossfit::<D, LdaModel>(data, &LdaConfig::default(), "LDA")
}

/// Reject tree folds that never produced a split.
fn ensure_tree_splits(
    model: &'static str,
    stats: &RescoreFeatureStats,
) -> Result<(), RescoreError> {
    let folds: Vec<u8> = stats
        .iter()
        .filter(|fs| fs.feature_importance.is_empty())
        .map(|fs| fs.fold)
        .collect();
    if folds.is_empty() {
        Ok(())
    } else {
        Err(RescoreError::UntrainedFolds { model, folds })
    }
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore(mut data: Vec<CompetedCandidate>) -> RescoreResult {
    let config = GBMConfig::default();

    canonicalize_and_shuffle(&mut data);

    // The ALL-lane matrix (linear ++ nonlinear), built over `data` in its
    // post-shuffle order — see `canonicalize_and_shuffle` for why that order is
    // load-bearing.
    let names = all_feature_name_set();
    debug_assert_eq!(names.len(), ALL_NCOLS);
    let feat = build_all_matrix(competed_rows(&data));
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let precomputed = PrecomputedFeatures::from_row_major(feat, ALL_NCOLS, responses);

    let mut scorer =
        CrossValidatedScorer::<CompetedCandidate, GbmFoldModel>::new_from_shuffled_with_precomputed(
            N_RESCORE_FOLDS,
            data,
            config,
            precomputed,
            names,
        );
    scorer.fit().map_err(|e| RescoreError::Model {
        model: "GBM",
        reason: e.to_string(),
    })?;

    let stats = scorer.feature_stats();
    ensure_tree_splits("GBM", &stats)?;

    Ok(finalize(scorer.score(), stats))
}

/// Sage-style shrinkage-LDA rescorer: closed-form linear fits, no boosting. It
/// is generally much cheaper and less sensitive than GBM, with neither gap
/// constant across candidate counts.
///
/// The FDR machinery (`assign_qval`, target-decoy competition) is untouched —
/// only the discriminant score source changes.
///
/// CROSS-FIT, not a single in-sample fit: every row's score comes from an LDA
/// fitted without that row, via [`crossfit_lda`] — see [`crossfit`] for the
/// partition and why it is mandatory.
///
/// Returns PER-FOLD [`FoldStats`] (one per fold, like the GBM path): feature
/// means/NaN ratios over each fold's held-out rows, `|coef|` importance from
/// that fold's model.
///
/// Selected via the `rescore_model` config field / `--rescore-model` CLI flag
/// ([`crate::ml::RescoreModel::Lda`]).
/// See `ml::lda` for the fit details.
pub fn rescore_lda(mut data: Vec<CompetedCandidate>) -> RescoreResult {
    // Canonical sort + seeded shuffle — the same helper, key and seed as every
    // other rescorer.
    canonicalize_and_shuffle(&mut data);

    // LDA trains on the LINEAR lane only: fields that are approx-Gaussian after
    // their declared per-row transform (raw/log2/ln1p). Skew-taming is done at
    // emit time by the grammar, so there is no data-dependent normalization step
    // here — the only remaining data-dependent op is LDA's own standardization.
    // Built after the shuffle, per `canonicalize_and_shuffle`.
    let names: Vec<Arc<str>> = linear_feature_name_set();
    let nrows = data.len();
    let ncols = LINEAR_NCOLS;
    debug_assert_eq!(names.len(), ncols);

    let dataset = StreamingDataset::new(
        &data,
        names.clone(),
        N_RESCORE_FOLDS as usize,
        write_competed_linear_row,
    );

    let cf = crossfit_lda(&dataset)?;
    debug!("LDA cross-fit scored {nrows} candidates across {N_RESCORE_FOLDS} folds");
    let stats = cf.feature_stats(&dataset);
    let scores = cf.scores;
    drop(dataset);
    for (cand, s) in data.iter_mut().zip(scores) {
        cand.assign_score(s);
    }

    Ok(finalize(data, stats))
}

/// The LINEAR-lane dataset a hybrid cross-fits its extra column over.
///
/// `data` must already be through [`canonicalize_and_shuffle`].
///
/// The fold count is `N_RESCORE_FOLDS`, i.e. the SAME
/// `get_fold` the GBM's [`CrossValidatedScorer`] derives its partition from
/// below. That shared definition is what makes "the same fold assignment on both
/// sides" structural rather than a comment; see [`crossfit`].
fn hybrid_linear_dataset(data: &[CompetedCandidate]) -> StreamingDataset<'_, CompetedCandidate> {
    StreamingDataset::new(
        data,
        linear_feature_name_set(),
        N_RESCORE_FOLDS as usize,
        write_competed_linear_row,
    )
}

/// The frame the hybrid hands the GBM: the NONLINEAR lane matrix with one
/// cross-fit score appended as the TRAILING column, plus the matching names.
///
/// Matrix values and names are built together so the appended column cannot
/// drift from its name or position.
///
/// `data` must be post-[`canonicalize_and_shuffle`], and `score` / `responses`
/// must be in that same row order.
fn hybrid_frame(
    data: &[CompetedCandidate],
    responses: Vec<f64>,
    score_name: &str,
    score: Vec<f64>,
) -> (PrecomputedFeatures, Vec<Arc<str>>) {
    let nrows = data.len();
    // Hard assert, not `debug_assert`: the zip below would silently TRUNCATE the
    // matrix to the shorter of the two and every width check downstream would
    // still pass on the truncated frame.
    assert_eq!(
        score.len(),
        nrows,
        "the appended {score_name} column must carry one value per row"
    );
    assert_eq!(responses.len(), nrows, "responses must be one per row");

    let ncols = NONLINEAR_NCOLS + 1;
    let nl = build_nonlinear_matrix(data);
    let mut feat = Vec::with_capacity(nrows * ncols);
    for (row, s) in nl.chunks_exact(NONLINEAR_NCOLS).zip(score) {
        feat.extend_from_slice(row);
        feat.push(s);
    }
    let mut names = nonlinear_feature_name_set();
    names.push(Arc::from(score_name));
    debug_assert_eq!(names.len(), ncols);
    debug_assert_eq!(feat.len(), nrows * ncols);

    (
        PrecomputedFeatures::from_row_major(feat, ncols, responses),
        names,
    )
}

/// Hybrid rescorer: cross-fit an LDA on the LINEAR lane, push its (leak-free)
/// `lda_score` as one extra column into the NONLINEAR lane, then train the GBM
/// CV on `nonlinear + lda_score` instead of the full feature frame.
///
/// LEAK-FREEDOM: `lda_score` is cross-fit via [`crossfit`] — see there for the
/// partition, why a label-aware feature fed to a CV'd GBM in particular must be
/// leak-free, and why the fold ASSIGNMENT has to match the one
/// `CrossValidatedScorer` derives its own partition from. ASSIGNMENT, not
/// partition: the two train/score splits differ deliberately and both are
/// leak-free, so "unifying" them is the refactor to not make.
///
/// Selected via the `rescore_model` config field / `--rescore-model` CLI flag
/// ([`crate::ml::RescoreModel::Hybrid`]).
#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore_hybrid(mut data: Vec<CompetedCandidate>) -> RescoreResult {
    let config = GBMConfig::default();

    // Canonical sort + seeded shuffle — IDENTICAL to `rescore` (same helper,
    // same key, same seed) so fold assignment and downstream q-values are
    // reproducible.
    canonicalize_and_shuffle(&mut data);

    // Both lane projections read `data` in its post-shuffle order (see
    // `canonicalize_and_shuffle`), so row `i` is the same candidate in the
    // streamed linear source, nonlinear frame, lda_score, responses, and moved data.
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let lin_dataset = hybrid_linear_dataset(&data);

    let lda_score = crossfit::<_, LdaModel>(&lin_dataset, &LdaConfig::default(), "LDA")?.scores;

    let (precomputed, names) = hybrid_frame(&data, responses, "lda_score", lda_score);

    let mut scorer =
        CrossValidatedScorer::<CompetedCandidate, GbmFoldModel>::new_from_shuffled_with_precomputed(
            N_RESCORE_FOLDS,
            data,
            config,
            precomputed,
            names,
        );
    scorer.fit().map_err(|e| RescoreError::Model {
        model: "hybrid GBM",
        reason: e.to_string(),
    })?;

    let stats = scorer.feature_stats();
    ensure_tree_splits("hybrid GBM", &stats)?;

    Ok(finalize(scorer.score(), stats))
}

/// `config` is a parameter so tests can train a smaller net than
/// [`MlpConfig::default`] without tuning the production default for test runtime.
/// Fold assignment and per-fold RNGs are deterministic after
/// [`canonicalize_and_shuffle`]. A failed fold returns [`RescoreError`].
fn rescore_mlp_with(mut data: Vec<CompetedCandidate>, config: MlpConfig) -> RescoreResult {
    canonicalize_and_shuffle(&mut data);

    // Fold assignment remains positional after this shuffle. Feature values are
    // streamed one row at a time into the MLP's two reusable batch buffers;
    // there is deliberately no raw or transformed fold matrix on this path.
    let names = all_feature_name_set();
    let ncols = ALL_NCOLS;
    debug_assert_eq!(names.len(), ncols);

    let mut scorer =
        CrossValidatedScorer::<CompetedCandidate, MlpFoldModel>::new_from_shuffled_streaming(
            N_RESCORE_FOLDS,
            data,
            config,
            names,
            write_competed_all_row,
        );
    scorer.fit_parallel().map_err(|e| RescoreError::Model {
        model: "MLP",
        reason: e.to_string(),
    })?;

    let stats = scorer.feature_stats();

    Ok(finalize(scorer.score(), stats))
}

/// MLP rescorer over the ALL lane (linear ++ nonlinear) — the same feature set
/// [`rescore`] gives the GBM, so the two are directly comparable. The one
/// [`crate::ml::RescoreModel::Mlp`] selects.
///
/// See [`rescore_mlp_with`] for the cross-fit and determinism contracts.
///
/// Runtime and sensitivity comparisons are not constant across candidate counts.
pub fn rescore_mlp(data: Vec<CompetedCandidate>) -> RescoreResult {
    rescore_mlp_with(data, MlpConfig::default())
}

// ---------------------------------------------------------------------------
// Lane feature projections — shared ordering for streaming and matrix consumers
// ---------------------------------------------------------------------------
//
// GBM consumers retain a flat row-major `Vec<f64>`; LDA and MLP project the same
// rows into one-row scratch and retain no raw frame. Every contributing block's
// width is an inherent const, so the lane widths below are compile-time constants.
//
// The per-row walk order is fixed — scoring blocks (composition order) ->
// `ResultMeta` -> `Derived` -> (nonlinear only) `sequence_counts` — and the
// name-set walks further down repeat it exactly, which is what keeps column
// `j` and name `j` the same feature.

/// LINEAR-lane width (LDA): fields approx-Gaussian after their declared
/// per-row transform.
const LINEAR_NCOLS: usize =
    ScoringFields::LINEAR_LEN + ResultMeta::LINEAR_LEN + Derived::LINEAR_LEN;

/// NONLINEAR-lane width: the rest (context, counts, sequence).
const NONLINEAR_NCOLS: usize = ScoringFields::NONLINEAR_LEN
    + ResultMeta::NONLINEAR_LEN
    + Derived::NONLINEAR_LEN
    + sequence_counts::LEN;

/// ALL-lane width (GBM) = linear ++ nonlinear.
const ALL_NCOLS: usize = LINEAR_NCOLS + NONLINEAR_NCOLS;

// Neither lane may collapse to nothing. Now that the widths are consts this is
// a build failure rather than a test failure — a lane that lost every feature
// would otherwise train a model on a zero-column matrix.
const _: () = assert!(LINEAR_NCOLS > 0, "linear lane collapsed");
const _: () = assert!(NONLINEAR_NCOLS > 0, "nonlinear lane collapsed");

trait ValueSink {
    fn push(&mut self, values: &[f64]);
}

impl ValueSink for Vec<f64> {
    fn push(&mut self, values: &[f64]) {
        self.extend_from_slice(values);
    }
}

struct SliceSink<'a> {
    out: &'a mut [f64],
    written: usize,
}

impl<'a> SliceSink<'a> {
    fn new(out: &'a mut [f64]) -> Self {
        Self { out, written: 0 }
    }

    fn finish(self) {
        assert_eq!(self.written, self.out.len());
    }
}

impl ValueSink for SliceSink<'_> {
    fn push(&mut self, values: &[f64]) {
        let end = self
            .written
            .checked_add(values.len())
            .expect("feature row width overflow");
        assert!(end <= self.out.len(), "feature row exceeds its lane width");
        self.out[self.written..end].copy_from_slice(values);
        self.written = end;
    }
}

/// Project one row's linear values into a streaming or retained sink.
fn project_linear_row(
    scoring: &ScoringFields,
    meta: &ResultMeta,
    derived: &Derived,
    out: &mut impl ValueSink,
) {
    out.push(&scoring.linear_feature_array());
    out.push(&meta.linear_feature_array());
    out.push(&derived.linear_feature_array());
}

/// Project one row's nonlinear values into a streaming or retained sink.
fn project_nonlinear_row(
    scoring: &ScoringFields,
    meta: &ResultMeta,
    derived: &Derived,
    out: &mut impl ValueSink,
) {
    out.push(&scoring.nonlinear_feature_array());
    out.push(&meta.nonlinear_feature_array());
    out.push(&derived.nonlinear_feature_array());
    out.push(&sequence_counts::nonlinear_feature_array(
        &scoring.identity.peptide,
    ));
}

/// Write one competed candidate's all-lane row into caller-owned scratch.
/// This is the MLP boundary: raw `f64` values live for one row only and are
/// transformed into a reusable `f32` batch buffer by the consumer.
fn write_competed_all_row(candidate: &CompetedCandidate, out: &mut [f64]) {
    assert_eq!(out.len(), ALL_NCOLS);
    let scoring = &candidate.scoring;
    let meta = candidate.result_meta();
    let derived = Derived::compute(scoring);
    let mut sink = SliceSink::new(out);
    project_linear_row(scoring, &meta, &derived, &mut sink);
    project_nonlinear_row(scoring, &meta, &derived, &mut sink);
    sink.finish();
}

/// Matrix-free LINEAR-lane projection for the standalone and hybrid LDA paths.
fn write_competed_linear_row(candidate: &CompetedCandidate, out: &mut [f64]) {
    assert_eq!(out.len(), LINEAR_NCOLS);
    let scoring = &candidate.scoring;
    let meta = candidate.result_meta();
    let derived = Derived::compute(scoring);
    let mut sink = SliceSink::new(out);
    project_linear_row(scoring, &meta, &derived, &mut sink);
    sink.finish();
}

/// The LINEAR-lane matrix for `data` in its CURRENT order (call AFTER any
/// shuffle, so row `i` aligns with `data[i]`). `LINEAR_NCOLS` wide.
#[cfg(test)]
fn build_linear_matrix(data: &[CompetedCandidate]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * LINEAR_NCOLS);
    for c in data {
        let meta = c.result_meta();
        project_linear_row(&c.scoring, &meta, &Derived::compute(&c.scoring), &mut out);
    }
    out
}

/// The NONLINEAR-lane matrix for `data`, `NONLINEAR_NCOLS` wide.
fn build_nonlinear_matrix(data: &[CompetedCandidate]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * NONLINEAR_NCOLS);
    for c in data {
        let meta = c.result_meta();
        project_nonlinear_row(&c.scoring, &meta, &Derived::compute(&c.scoring), &mut out);
    }
    out
}

/// The ALL-lane matrix (linear then nonlinear, per row) — the GBM feature set,
/// `ALL_NCOLS` wide, matching [`all_feature_name_set`]'s order.
///
/// ONE pass over `rows` with ONE `Derived::compute` per row: the two lanes are
/// adjacent within a row, so there is nothing to gain from walking twice.
///
/// Takes `(scoring, meta)` pairs rather than a row type, because both sides of
/// rescoring feed it and they agree on nothing else. See [`competed_rows`] for
/// the pre-rescore side and [`feature_frame`] for the post-rescore one.
fn build_all_matrix<'a>(
    rows: impl ExactSizeIterator<Item = (&'a ScoringFields, ResultMeta)>,
) -> Vec<f64> {
    let mut out = Vec::with_capacity(rows.len() * ALL_NCOLS);
    for (s, meta) in rows {
        let derived = Derived::compute(s);
        project_linear_row(s, &meta, &derived, &mut out);
        project_nonlinear_row(s, &meta, &derived, &mut out);
    }
    out
}

/// Competed candidates in the shape [`build_all_matrix`] consumes.
fn competed_rows(
    data: &[CompetedCandidate],
) -> impl ExactSizeIterator<Item = (&ScoringFields, ResultMeta)> {
    data.iter().map(|c| (&c.scoring, c.result_meta()))
}

/// The ALL-lane feature names + matrix for post-rescore rows: the dashboard's
/// entry point.
///
/// Row-major: value `j` of row `i` is at `matrix[i * names.len() + j]`.
pub fn feature_frame(data: &[FinalResult]) -> (Vec<Arc<str>>, Vec<f64>) {
    let rows = data.iter().map(|r| (&r.scoring, r.result_meta()));
    (all_feature_name_set(), build_all_matrix(rows))
}

/// LINEAR-lane feature names (LDA), in [`project_linear_row`]'s order.
pub fn linear_feature_name_set() -> Vec<Arc<str>> {
    let mut n = NameSink::new();
    <ScoringFields as ScoreBlock>::linear_feature_names(&mut n);
    <ResultMeta as ScoreBlock>::linear_feature_names(&mut n);
    <Derived as ScoreBlock>::linear_feature_names(&mut n);
    n.into_names()
}

/// NONLINEAR-lane feature names, in [`project_nonlinear_row`]'s order. The
/// `sequence_counts` names are unconditional — a peptide with no parsed
/// sequence contributes NaN values under them, not a shorter row.
pub fn nonlinear_feature_name_set() -> Vec<Arc<str>> {
    let mut n = NameSink::new();
    <ScoringFields as ScoreBlock>::nonlinear_feature_names(&mut n);
    <ResultMeta as ScoreBlock>::nonlinear_feature_names(&mut n);
    <Derived as ScoreBlock>::nonlinear_feature_names(&mut n);
    sequence_counts::nonlinear_feature_names(&mut n);
    n.into_names()
}

/// The ALL-lane feature names (GBM) = linear ++ nonlinear, matching
/// [`build_all_matrix`]'s column order.
pub fn all_feature_name_set() -> Vec<Arc<str>> {
    let mut v = linear_feature_name_set();
    v.extend(nonlinear_feature_name_set());
    v
}

// ---------------------------------------------------------------------------
// Label / score plumbing. `FeatureLike` is the label+score half only — the
// GBM/LDA feature values reach the scorer as a prebuilt lane frame
// (`PrecomputedFeatures::from_row_major`); the MLP supplies the same projection
// through `write_competed_all_row` one row at a time. Neither route needs
// `FeatureVector` (the self-walking supertrait used only by
// `CrossValidatedScorer::new_from_shuffled`).
// ---------------------------------------------------------------------------

impl FeatureLike for CompetedCandidate {
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
            4. / 8.,
            4. / 8.,
            4. / 8.,
            6. / 8.,
            6. / 8.,
            6. / 8.,
            6. / 8.,
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
    use crate::ml::lda::DEFAULT_SHRINKAGE;
    use crate::models::DecoyMarking;
    use crate::models::sequence::Peptide;
    use crate::scoring::results::{
        CompetedCandidate,
        ScoringFields,
    };
    use std::sync::Arc;

    fn base_scoring_fields(peptide: Peptide) -> ScoringFields {
        ScoringFields::sample(peptide)
    }

    /// A target candidate over `PEPTIDEK` — 8 residues, no mods. With
    /// `sequence_features` set the sequence lanes carry those counts; without it
    /// they stay NaN, which is the only difference between the two cases.
    fn sample_competed_candidate(sequence_features: bool) -> CompetedCandidate {
        let peptide = Peptide {
            raw: Arc::from("PEPTIDEK"),
            decoy: DecoyMarking::Target,
            decoy_group: 0,
            sequence_features,
        };
        CompetedCandidate {
            scoring: base_scoring_fields(peptide),
            delta_group_ln1p_diff: 1.0,
            delta_group_ln1p_ratio: 0.5,
            discriminant_score: 0.0,
            qvalue: 1.0,
        }
    }

    // --- Lane walks (the live ML path) ---

    /// LANE WIDTH PARITY: a lane matrix's row width MUST equal that lane's
    /// name-set length. The two are built by separate walks (one positional
    /// over `[f64; N]` arrays, one over `NameSink`), so this is the assertion
    /// that keeps column `j` and name `j` the same feature; a drift silently
    /// misattributes feature importances and stats.
    ///
    /// Checked with an UNPARSED peptide as well as a parsed one: the sequence
    /// block is unconditional now, so a row must be exactly as wide either way
    /// — the width being a compile-time const is what makes that structural.
    #[test]
    fn lane_matrix_widths_match_name_sets() {
        for data in [
            vec![sample_competed_candidate(true)],
            vec![sample_competed_candidate(false)],
        ] {
            assert_eq!(linear_feature_name_set().len(), LINEAR_NCOLS);
            assert_eq!(nonlinear_feature_name_set().len(), NONLINEAR_NCOLS);
            assert_eq!(all_feature_name_set().len(), ALL_NCOLS);

            assert_eq!(build_linear_matrix(&data).len(), LINEAR_NCOLS);
            assert_eq!(build_nonlinear_matrix(&data).len(), NONLINEAR_NCOLS);
            assert_eq!(build_all_matrix(competed_rows(&data)).len(), ALL_NCOLS);
        }
    }

    /// The all-lane matrix is exactly the linear row followed by the nonlinear
    /// row, per row — the one-pass build must not reorder or drop anything
    /// relative to the two single-lane builds.
    #[test]
    fn all_matrix_is_linear_then_nonlinear_per_row() {
        let data = vec![
            sample_competed_candidate(true),
            sample_competed_candidate(false),
        ];
        let lin = build_linear_matrix(&data);
        let nl = build_nonlinear_matrix(&data);
        let all = build_all_matrix(competed_rows(&data));

        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        for i in 0..data.len() {
            let row = &all[i * ALL_NCOLS..(i + 1) * ALL_NCOLS];
            let mut streamed_linear = vec![0.0; LINEAR_NCOLS];
            write_competed_linear_row(&data[i], &mut streamed_linear);
            assert_eq!(
                bits(&lin[i * LINEAR_NCOLS..(i + 1) * LINEAR_NCOLS]),
                bits(&streamed_linear),
                "the matrix-free LDA row writer must preserve the canonical linear projection"
            );
            let mut streamed = vec![0.0; ALL_NCOLS];
            write_competed_all_row(&data[i], &mut streamed);
            assert_eq!(
                bits(row),
                bits(&streamed),
                "the matrix-free MLP row writer must preserve the canonical lane projection"
            );
            assert_eq!(
                bits(&row[..LINEAR_NCOLS]),
                bits(&lin[i * LINEAR_NCOLS..(i + 1) * LINEAR_NCOLS])
            );
            assert_eq!(
                bits(&row[LINEAR_NCOLS..]),
                bits(&nl[i * NONLINEAR_NCOLS..(i + 1) * NONLINEAR_NCOLS])
            );
        }
    }

    /// An unparsed peptide preserves the feature width by emitting NaNs.
    #[test]
    fn unparsed_sequence_emits_nan_not_a_narrower_row() {
        let names = nonlinear_feature_name_set();
        let seq_start = NONLINEAR_NCOLS - sequence_counts::LEN;
        assert_eq!(&*names[seq_start], "peptide_length");
        assert_eq!(&*names[NONLINEAR_NCOLS - 1], "peptide_n_mods");

        let unparsed = build_nonlinear_matrix(&[sample_competed_candidate(false)]);
        assert!(
            unparsed[seq_start..].iter().all(|v| v.is_nan()),
            "unparsed sequence features must all be NaN: {:?}",
            &unparsed[seq_start..]
        );

        let parsed = build_nonlinear_matrix(&[sample_competed_candidate(true)]);
        // PEPTIDEK: 8 residues, no mods.
        assert_eq!(parsed[seq_start], 8.0);
        assert_eq!(parsed[NONLINEAR_NCOLS - 1], 0.0);
    }

    #[test]
    fn feature_names_are_unique() {
        // Exact counts may change, but duplicate names make reports ambiguous.
        let mut seen = std::collections::HashSet::new();
        for n in all_feature_name_set() {
            assert!(seen.insert(n.clone()), "dup feature name: {n}");
        }
    }

    /// Build a non-degenerate synthetic candidate set: `n` rows, alternating
    /// target/decoy, distinct `library_id`, with the LINEAR-lane count fields
    /// varied by label + row so the cross-fit LDA has real within-class scatter
    /// and a class-mean gap (i.e. it actually fits, exercising the score path).
    ///
    /// Every field varied here lands in the linear lane, including the two
    /// log-space group-delta features. The nonlinear lane is therefore constant.
    ///
    /// The nonlinear lane is intentionally constant; use
    /// [`synthetic_competed_nonlinear_signal`] when a nonlinear split is needed.
    fn synthetic_competed(n: u32) -> Vec<CompetedCandidate> {
        (0..n)
            .map(|i| {
                let mut c = sample_competed_candidate(true);
                c.scoring.identity.library_id = i as u64;
                let is_target = i % 2 == 0;
                c.scoring.identity.is_target = is_target;
                let base: u8 = if is_target { 20 } else { 8 };
                let jitter = (i % 5) as u8;
                c.scoring.counts.rising_cycles = base + jitter;
                c.scoring.counts.falling_cycles = base.saturating_sub(jitter);
                c.scoring.counts.npeaks = base + (i % 3) as u8;
                c.scoring.finalize_counts.n_scored_fragments = base + (i % 4) as u8;
                c.delta_group_ln1p_diff = if is_target { 2.0 } else { 0.5 } + (i % 7) as f32 * 0.1;
                c.delta_group_ln1p_ratio = if is_target { 0.8 } else { 0.3 };
                c
            })
            .collect()
    }

    fn synthetic_competed_linear_only(n: u32) -> Vec<CompetedCandidate> {
        synthetic_competed(n)
            .into_iter()
            .map(|mut c| {
                let sample = sample_competed_candidate(true);
                c.delta_group_ln1p_diff = sample.delta_group_ln1p_diff;
                c.delta_group_ln1p_ratio = sample.delta_group_ln1p_ratio;
                c
            })
            .collect()
    }

    #[test]
    fn crossfit_lda_scores_are_held_out() {
        const N: usize = 60;
        let n_folds = N_RESCORE_FOLDS as usize;
        let is_target = |i: usize| i.is_multiple_of(2);

        let mut feat = Vec::with_capacity(N * 3);
        let mut is_decoy = Vec::with_capacity(N);
        for i in 0..N {
            let label = if is_target(i) { 1.0 } else { 0.0 };
            let in_fold0 = i % n_folds == 0;
            feat.push(if in_fold0 { label } else { 0.0 }); // col0
            feat.push(if in_fold0 { 0.0 } else { label }); // col1
            feat.push(((i / (2 * n_folds)) % 5) as f64); // col2
            is_decoy.push(!is_target(i));
        }

        let responses: Vec<f64> = is_decoy
            .iter()
            .map(|&d| if d { 0.0 } else { 1.0 })
            .collect();
        let names: Vec<Arc<str>> = ["col0", "col1", "col2"].map(Arc::from).to_vec();
        let dataset = RowMajorDataset::new(
            PrecomputedFeatures::from_row_major(feat.clone(), 3, responses),
            names,
            n_folds,
        );

        let cf = crossfit_lda(&dataset).expect("cross-fit must succeed");
        assert_eq!(
            cf.models.len(),
            n_folds,
            "cross-fit must produce one model per fold"
        );

        let fold0: Vec<usize> = (0..N).step_by(n_folds).collect();
        assert!(fold0.iter().any(|&i| is_target(i)) && fold0.iter().any(|&i| !is_target(i)));

        let held: Vec<f64> = fold0.iter().map(|&i| cf.scores[i]).collect();
        let spread = held.iter().cloned().fold(f64::MIN, f64::max)
            - held.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            spread < 1e-6,
            "held-out fold-0 scores must not separate; spread={spread} scores={held:?}"
        );

        let insample = LdaModel::fit_matrix(&feat, N, 3, &is_decoy, DEFAULT_SHRINKAGE).unwrap();
        let (mut min_t, mut max_d) = (f64::MAX, f64::MIN);
        for &i in &fold0 {
            let s = insample.score(&feat[i * 3..(i + 1) * 3]);
            if is_target(i) {
                min_t = min_t.min(s);
            } else {
                max_d = max_d.max(s);
            }
        }
        let insample_gap = min_t - max_d;
        assert!(
            insample_gap > 1e-3,
            "in-sample control must separate fold-0 rows; gap={insample_gap}"
        );
        assert!(
            insample_gap > 1e3 * spread.max(f64::MIN_POSITIVE),
            "hold-out separation ({spread}) must be negligible next to in-sample ({insample_gap})"
        );
    }

    #[test]
    fn rescore_lda_is_cross_fit_and_deterministic() {
        let n = 90;
        let (out_a, stats_a) = rescore_lda(synthetic_competed(n)).unwrap();
        let (out_b, _stats_b) = rescore_lda(synthetic_competed(n)).unwrap();

        assert_eq!(out_a.len(), n as usize);

        assert_eq!(stats_a.len(), N_RESCORE_FOLDS as usize);
        for (f, fs) in stats_a.iter().enumerate() {
            assert_eq!(fs.fold, f as u8);
            assert!(
                !fs.feature_stats.is_empty(),
                "fold {f} has no feature stats"
            );
            assert!(
                !fs.feature_importance.is_empty(),
                "fold {f} has no importance"
            );
        }

        for r in &out_a {
            assert!(r.discriminant_score.is_finite());
            assert!((0.0..=1.0).contains(&r.qvalue));
        }

        let key = |out: &[FinalResult]| -> Vec<(u64, u32)> {
            let mut v: Vec<(u64, u32)> = out
                .iter()
                .map(|r| {
                    (
                        r.scoring.identity.library_id,
                        r.discriminant_score.to_bits(),
                    )
                })
                .collect();
            v.sort_unstable();
            v
        };
        assert_eq!(key(&out_a), key(&out_b), "lda rescore not deterministic");
    }

    #[test]
    fn lda_sidecar_reports_every_linear_feature_including_zero_weights() {
        let (_out, stats) = rescore_lda(synthetic_competed(90)).unwrap();
        assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);

        for fs in &stats {
            assert_eq!(
                fs.feature_importance.len(),
                LINEAR_NCOLS,
                "fold {}: LDA must report one importance per linear-lane feature",
                fs.fold
            );
            assert!(
                fs.feature_importance.iter().all(|(_, g)| g.is_finite()),
                "fold {}: LDA importance is never the 'unreported' NAN sentinel",
                fs.fold
            );
            assert!(
                fs.feature_importance.iter().any(|(_, g)| *g == 0.0),
                "fold {}: this fixture has constant linear features, whose \
                 zero weights must survive to the sidecar",
                fs.fold
            );
        }

        let lane: std::collections::HashSet<Arc<str>> =
            linear_feature_name_set().into_iter().collect();
        let reported: std::collections::HashSet<Arc<str>> = stats[0]
            .feature_importance
            .iter()
            .map(|(n, _)| n.clone())
            .collect();
        assert_eq!(reported, lane);
    }

    fn mlp_test_cfg(seed: u64) -> MlpConfig {
        MlpConfig {
            hidden: vec![16, 8],
            lr: 1e-3,
            epochs: 150,
            batch_size: 16,
            seed,
            ..MlpConfig::default()
        }
    }

    fn pair_auc(scores: &[f64], is_target: &[bool]) -> f64 {
        let pick = |want: bool| -> Vec<f64> {
            scores
                .iter()
                .zip(is_target)
                .filter(|&(_, &b)| b == want)
                .map(|(&s, _)| s)
                .collect()
        };
        let (t, d) = (pick(true), pick(false));
        assert!(!t.is_empty() && !d.is_empty(), "AUC needs both classes");
        let mut acc = 0.0;
        for a in &t {
            for b in &d {
                acc += match a.partial_cmp(b) {
                    Some(std::cmp::Ordering::Greater) => 1.0,
                    Some(std::cmp::Ordering::Equal) => 0.5,
                    _ => 0.0,
                };
            }
        }
        acc / (t.len() * d.len()) as f64
    }

    struct LabelRow {
        y: f64,
        score: f64,
    }

    impl FeatureLike for LabelRow {
        fn get_y(&self) -> f64 {
            self.y
        }

        fn assign_score(&mut self, score: f64) {
            self.score = score;
        }

        fn get_score(&self) -> f64 {
            self.score
        }
    }

    #[test]
    fn mlp_crossfit_scores_are_held_out() {
        const N: usize = 90;
        let n_folds = N_RESCORE_FOLDS as usize;
        let is_target = |i: usize| i.is_multiple_of(2);

        let mut feat = Vec::with_capacity(N * 3);
        let mut responses = Vec::with_capacity(N);
        for i in 0..N {
            let label = if is_target(i) { 1.0 } else { 0.0 };
            let in_fold0 = i % n_folds == 0;
            feat.push(if in_fold0 { label } else { 0.0 }); // col0
            feat.push(if in_fold0 { 0.0 } else { label }); // col1
            feat.push(((i / (2 * n_folds)) % 5) as f64); // col2
            responses.push(label);
        }
        let names: Vec<Arc<str>> = ["col0", "col1", "col2"].map(Arc::from).to_vec();

        let fold0: Vec<usize> = (0..N).step_by(n_folds).collect();
        let fold0_is_target: Vec<bool> = fold0.iter().map(|&i| is_target(i)).collect();
        assert!(fold0_is_target.iter().any(|&b| b) && fold0_is_target.iter().any(|&b| !b));

        for seed in [7u64, 13, 42, 1234] {
            let rows: Vec<LabelRow> = responses
                .iter()
                .map(|&y| LabelRow { y, score: 0.0 })
                .collect();
            let mut scorer =
                CrossValidatedScorer::<LabelRow, MlpFoldModel>::new_from_shuffled_with_precomputed(
                    N_RESCORE_FOLDS,
                    rows,
                    mlp_test_cfg(seed),
                    PrecomputedFeatures::from_row_major(feat.clone(), 3, responses.clone()),
                    names.clone(),
                );
            scorer.fit().expect("cross-fit must succeed");
            let scores = scorer.get_scores().unwrap();

            let held: Vec<f64> = fold0.iter().map(|&i| scores[i]).collect();
            let held_auc = pair_auc(&held, &fold0_is_target);
            assert!(
                (held_auc - 0.5).abs() < 1e-12,
                "seed {seed}: held-out fold-0 rows must not separate — the only \
                 column that could separate them was culled from the model that \
                 scores them. AUC={held_auc} scores={held:?}"
            );

            let dataset = RowMajorDataset::new(
                PrecomputedFeatures::from_row_major(feat.clone(), 3, responses.clone()),
                names.clone(),
                n_folds,
            );
            let all: Vec<usize> = (0..N).collect();
            let insample = MlpFoldModel::fit(&mlp_test_cfg(seed), &dataset, 0, &all, &[]).unwrap();
            assert!(
                insample.transform().culled().is_empty(),
                "seed {seed}: col0 must survive the cull in sample, or the \
                 control cannot separate and proves nothing"
            );
            let insample_auc = pair_auc(
                &insample.predict(&dataset, &fold0).unwrap(),
                &fold0_is_target,
            );
            assert!(
                insample_auc > 0.99,
                "seed {seed}: in-sample control must separate fold-0 rows \
                 (AUC={insample_auc}); without it the hold-out assertion is vacuous"
            );

            let leaky = MlpFoldModel::fit(&mlp_test_cfg(seed), &dataset, 0, &fold0, &[]).unwrap();
            let leaky_auc = pair_auc(&leaky.predict(&dataset, &fold0).unwrap(), &fold0_is_target);
            assert!(
                leaky_auc > 0.99,
                "seed {seed}: a model trained ON fold 0 must separate fold 0 \
                 (AUC={leaky_auc}), else (a) would hold even for a leaking scorer"
            );
        }
    }

    #[test]
    fn rescore_mlp_is_deterministic() {
        let n = 90;
        let key = |out: &[FinalResult]| -> Vec<(u64, u32)> {
            let mut v: Vec<(u64, u32)> = out
                .iter()
                .map(|r| {
                    (
                        r.scoring.identity.library_id,
                        r.discriminant_score.to_bits(),
                    )
                })
                .collect();
            v.sort_unstable();
            v
        };

        for seed in [7u64, 13, 42, 1234] {
            let run = || rescore_mlp_with(synthetic_competed(n), mlp_test_cfg(seed));
            let (out_a, stats_a) = run().unwrap();
            let (out_b, _) = run().unwrap();

            assert_eq!(out_a.len(), n as usize);
            for r in &out_a {
                assert!(
                    r.discriminant_score.is_finite(),
                    "seed {seed}: non-finite score {}",
                    r.discriminant_score
                );
                assert!((0.0..=1.0).contains(&r.qvalue));
            }
            assert_eq!(
                key(&out_a),
                key(&out_b),
                "seed {seed}: mlp rescore not deterministic"
            );

            assert_eq!(stats_a.len(), N_RESCORE_FOLDS as usize);
            for (f, fs) in stats_a.iter().enumerate() {
                assert_eq!(fs.fold, f as u8);
                assert_eq!(
                    fs.feature_stats.len(),
                    ALL_NCOLS,
                    "seed {seed}: fold {f} stats must be all-feature width"
                );
                assert_eq!(
                    fs.feature_importance.len(),
                    ALL_NCOLS,
                    "seed {seed}: fold {f} — the MLP measures every lane \
                         column (0.0 for culled ones), so none may be dropped"
                );
                assert!(
                    fs.feature_importance.iter().all(|(_, g)| g.is_finite()),
                    "seed {seed}: fold {f} must never emit the NAN sentinel"
                );
                assert!(
                    fs.feature_importance.iter().any(|(_, g)| *g > 0.0),
                    "seed {seed}: fold {f} reported no weight at all"
                );
            }
        }
    }

    #[test]
    fn mlp_streamed_rows_are_aligned_with_the_shuffled_data() {
        for seed in [7u64, 13, 42, 1234] {
            let (out, _) = rescore_mlp_with(synthetic_competed(120), mlp_test_cfg(seed)).unwrap();
            let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
            let is_target: Vec<bool> = out.iter().map(|r| r.scoring.identity.is_target).collect();
            let auc = pair_auc(&scores, &is_target);
            assert!(
                auc > 0.9,
                "seed {seed}: AUC {auc} is near chance; streamed rows may not align with \
                 shuffled candidates"
            );
        }
    }

    #[test]
    fn rescore_hybrid_smoke_and_determinism() {
        let n = 360;
        let (out_a, stats_a) = rescore_hybrid(synthetic_competed(n)).unwrap();
        let (out_b, _stats_b) = rescore_hybrid(synthetic_competed(n)).unwrap();

        assert_eq!(stats_a.len(), N_RESCORE_FOLDS as usize);
        for fs in &stats_a {
            assert!(
                !fs.feature_importance.is_empty(),
                "fold {}: the GBM reported no importance at all, so it built no trees and \
                 the score is a per-fold constant — every other assertion in this test \
                 then passes without measuring anything (this is what n = 90 did)",
                fs.fold
            );
        }

        assert_eq!(out_a.len(), n as usize);
        assert_eq!(out_b.len(), n as usize);

        for r in &out_a {
            assert!(
                r.discriminant_score.is_finite(),
                "non-finite discriminant score: {}",
                r.discriminant_score
            );
            assert!(
                (0.0..=1.0).contains(&r.qvalue),
                "qvalue out of [0,1]: {}",
                r.qvalue
            );
        }

        let key = |out: &[FinalResult]| -> Vec<(u64, u32)> {
            let mut v: Vec<(u64, u32)> = out
                .iter()
                .map(|r| {
                    (
                        r.scoring.identity.library_id,
                        r.discriminant_score.to_bits(),
                    )
                })
                .collect();
            v.sort_unstable();
            v
        };
        assert_eq!(
            key(&out_a),
            key(&out_b),
            "hybrid rescore not deterministic across runs"
        );
    }

    struct FoldSpy {
        train: Vec<usize>,
        ncols: usize,
    }

    impl FoldModel for FoldSpy {
        type Config = ();
        type Error = std::convert::Infallible;

        fn fit<D: FoldDataset>(
            _cfg: &(),
            data: &D,
            _fold: usize,
            train: &[usize],
            _val: &[usize],
        ) -> Result<Self, std::convert::Infallible> {
            Ok(FoldSpy {
                train: train.to_vec(),
                ncols: data.column_names().len(),
            })
        }

        fn predict<D: FoldDataset>(
            &self,
            _data: &D,
            rows: &[usize],
        ) -> Result<Vec<f64>, std::convert::Infallible> {
            Ok(vec![0.0; rows.len()])
        }

        fn importance(&self) -> Vec<f32> {
            vec![0.0; self.ncols]
        }
    }

    struct ShortPredict;

    impl FoldModel for ShortPredict {
        type Config = ();
        type Error = std::convert::Infallible;

        fn fit<D: FoldDataset>(
            _cfg: &(),
            _data: &D,
            _fold: usize,
            _train: &[usize],
            _val: &[usize],
        ) -> Result<Self, Self::Error> {
            Ok(Self)
        }

        fn predict<D: FoldDataset>(
            &self,
            _data: &D,
            rows: &[usize],
        ) -> Result<Vec<f64>, Self::Error> {
            Ok(vec![0.0; rows.len().saturating_sub(1)])
        }

        fn importance(&self) -> Vec<f32> {
            Vec::new()
        }
    }

    #[test]
    fn crossfit_rejects_incomplete_prediction_vectors() {
        let mut data = synthetic_competed(12);
        canonicalize_and_shuffle(&mut data);
        let dataset = hybrid_linear_dataset(&data);

        let error = match crossfit::<_, ShortPredict>(&dataset, &(), "short predictor") {
            Err(error) => error,
            Ok(_) => panic!("short predictions must fail the cross-fit"),
        };
        assert!(matches!(
            error,
            RescoreError::CrossFit {
                fold: 0,
                stage: "predict",
                ..
            }
        ));
    }

    #[test]
    fn crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on() {
        const N: usize = 47;

        let mut data = synthetic_competed(N as u32);
        canonicalize_and_shuffle(&mut data);
        let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();

        let lin_dataset = hybrid_linear_dataset(&data);
        let cf = crossfit::<_, FoldSpy>(&lin_dataset, &(), "spy").expect("the spy cannot fail");

        let (precomputed, names) = hybrid_frame(&data, responses, "spy_score", cf.scores.clone());
        let scorer =
            CrossValidatedScorer::<CompetedCandidate, GbmFoldModel>::new_from_shuffled_with_precomputed(
                N_RESCORE_FOLDS,
                data,
                GBMConfig::default(),
                precomputed,
                names,
            );

        assert_eq!(
            cf.fold_rows,
            scorer.fold_rows(),
            "the cross-fit and the GBM scorer partition the rows differently, so a \
             hybrid row's cross-fit column can come from a model trained on rows \
             the GBM is holding out"
        );

        for (f, model) in cf.models.iter().enumerate() {
            let scored = &cf.fold_rows[f];
            let expected: Vec<usize> = (0..N).filter(|i| !scored.contains(i)).collect();
            assert_eq!(model.train, expected, "fold {f} trained on the wrong rows");
            assert!(
                !scored.is_empty() && !expected.is_empty(),
                "fold {f} is empty, so the equality above is vacuous"
            );
        }

        let mut seen: Vec<usize> = cf.fold_rows.iter().flatten().copied().collect();
        seen.sort_unstable();
        assert_eq!(seen, (0..N).collect::<Vec<_>>());
    }

    fn indistinguishable_competed(n: u32) -> Vec<CompetedCandidate> {
        (0..n)
            .map(|i| {
                let mut c = sample_competed_candidate(true);
                c.scoring.identity.library_id = i as u64;
                c.scoring.identity.is_target = i % 2 == 0;
                c
            })
            .collect()
    }

    #[test]
    fn standalone_mlp_rejects_an_untrainable_frame() {
        let error = rescore_mlp_with(indistinguishable_competed(24), mlp_test_cfg(7)).unwrap_err();
        assert!(matches!(error, RescoreError::Model { model: "MLP", .. }));
    }

    #[test]
    fn lda_and_hybrid_reject_an_untrainable_linear_frame() {
        for error in [
            rescore_lda(indistinguishable_competed(24)).unwrap_err(),
            rescore_hybrid(indistinguishable_competed(24)).unwrap_err(),
        ] {
            assert!(matches!(
                error,
                RescoreError::CrossFit {
                    model: "LDA",
                    stage: "fit",
                    ..
                }
            ));
        }
    }

    #[test]
    fn untrained_gbm_folds_are_errors() {
        let stats = vec![crate::ml::cv::FoldStats {
            fold: 2,
            feature_stats: Vec::new(),
            feature_importance: Vec::new(),
        }];

        assert_eq!(
            ensure_tree_splits("GBM", &stats),
            Err(RescoreError::UntrainedFolds {
                model: "GBM",
                folds: vec![2],
            })
        );
    }

    fn assert_nonlinear_lane_is_flat(fixture: &[CompetedCandidate], column: &str) {
        let nl = build_nonlinear_matrix(fixture);
        for (j, name) in nonlinear_feature_name_set().iter().enumerate() {
            let first = nl[j].to_bits();
            assert!(
                (0..fixture.len()).all(|i| nl[i * NONLINEAR_NCOLS + j].to_bits() == first),
                "fixture premise broken: nonlinear feature {name} varies across rows, so \
                 the GBM could separate without {column} and the assertions that follow \
                 would prove nothing"
            );
        }
    }

    #[test]
    fn hybrid_lda_score_carries_the_linear_lane_into_the_gbm() {
        assert_nonlinear_lane_is_flat(&synthetic_competed_linear_only(360), "lda_score");

        let (out, stats) = rescore_hybrid(synthetic_competed_linear_only(360)).unwrap();

        let split_on_it = stats
            .iter()
            .filter(|fs| {
                fs.feature_importance
                    .iter()
                    .any(|(n, _)| &**n == "lda_score")
            })
            .count();
        assert!(
            split_on_it > 0,
            "no fold split on lda_score even though it is the only non-constant \
             column, so the appended column never reached the GBM or is constant \
             (a silently failed cross-fit)"
        );

        let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
        let is_target: Vec<bool> = out.iter().map(|r| r.scoring.identity.is_target).collect();
        let auc = pair_auc(&scores, &is_target);
        assert!(
            auc > 0.8,
            "AUC {auc} — the GBM sees a constant nonlinear lane here, so anything above \
             chance has to come through lda_score, and it did not"
        );
    }

    #[test]
    fn neutralize_mobility_nans_every_mobility_feature() {
        // For a run with no searchable mobility axis (mzML/FAIMS), every
        // mobility-derived GBM feature must become NaN (`forust` missing), so it
        // cannot bias the score with sentinel-derived constants. Walked over the
        // LIVE lane path (`build_all_matrix` == linear ++ nonlinear), not a
        // hand-listed field set, so a new mobility field that forgets to
        // neutralize fails here.
        //
        // Two things are pinned, and both matter:
        //   (a) every mobility feature is NaN (or, for the `_isna` companions,
        //       flipped to 1.0 == "missing" — an isna column is BY CONSTRUCTION
        //       never NaN, so demanding NaN there would be wrong);
        //   (b) every non-mobility feature is bit-for-bit unchanged. Without (b)
        //       an impl that NaN'd the whole record would pass (a).
        let before = build_all_matrix(competed_rows(&[sample_competed_candidate(true)]));

        let mut cand = sample_competed_candidate(true);
        cand.scoring.neutralize_mobility();
        let after = build_all_matrix(competed_rows(&[cand]));

        let names = all_feature_name_set();
        assert_eq!(names.len(), ALL_NCOLS);
        assert_eq!(before.len(), ALL_NCOLS);
        assert_eq!(after.len(), ALL_NCOLS);

        let is_mobility = |n: &str| n.contains("mob");
        let mut nan_mob: Vec<&str> = Vec::new();
        let mut isna_mob: Vec<&str> = Vec::new();
        let mut finite_non_mob = 0usize;

        for (j, name) in names.iter().enumerate() {
            let (b, a) = (before[j], after[j]);
            if !is_mobility(name) {
                assert!(
                    a == b || (a.is_nan() && b.is_nan()),
                    "non-mobility feature {name} changed: {b} -> {a}"
                );
                if a.is_finite() {
                    finite_non_mob += 1;
                }
            } else if name.ends_with("_isna") {
                assert_eq!(a, 1.0, "mobility isna companion {name} should be 1.0");
                isna_mob.push(name);
            } else {
                assert!(a.is_nan(), "mobility feature {name} should be NaN, got {a}");
                nan_mob.push(name);
            }
        }

        // No exact counts: adding a mobility-derived feature must not break
        // this test, since the assertions above are already per-feature. The
        // only thing worth pinning is that the `mob` name filter still matches
        // something — otherwise the loop above would pass vacuously.
        assert!(
            !nan_mob.is_empty(),
            "no mobility features matched the name filter"
        );
        assert!(
            !isna_mob.is_empty(),
            "no mobility _isna companions matched the name filter"
        );
        assert!(
            finite_non_mob > 20,
            "non-mobility features must stay finite, got {finite_non_mob}"
        );
    }

    /// Row-major layout: value `j` of row `i` lives at `matrix[i * nf + j]`.
    ///
    /// Every consumer indexes on that contract — `rescore_dash` sweeps the
    /// matrix a row at a time with every column's accumulator live, so an
    /// interleaved or transposed write would silently mix features together
    /// rather than fail. Rows carry distinct `delta_group_ln1p_diff` values so the
    /// assertion can tell them apart; a length check alone passes even when
    /// the layout is wrong.
    #[test]
    fn feature_frame_rows_are_contiguous() {
        let rows: Vec<_> = [1.0f32, 2.0, 3.0]
            .into_iter()
            .map(|delta_group_ln1p_diff| {
                let mut c = sample_competed_candidate(true);
                c.delta_group_ln1p_diff = delta_group_ln1p_diff;
                c.into_final()
            })
            .collect();

        let (names, got) = feature_frame(&rows);
        let nf = names.len();
        assert_eq!(got.len(), rows.len() * nf, "one value per name per row");

        let j = names
            .iter()
            .position(|n| &**n == "delta_group_ln1p_diff")
            .expect("delta_group_ln1p_diff is an ALL-lane feature");
        for (i, expected) in [1.0f64, 2.0, 3.0].into_iter().enumerate() {
            assert_eq!(
                got[i * nf + j],
                expected,
                "row {i}'s delta_group_ln1p_diff is not at matrix[{i} * {nf} + {j}]"
            );
        }
    }
}
