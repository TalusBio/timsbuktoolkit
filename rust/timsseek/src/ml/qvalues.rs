use super::cv::{
    CrossValidatedScorer,
    FeatureLike,
    FoldDataset,
    FoldModel,
    FoldStats,
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
use super::mlp_fold::{
    MlpFoldError,
    MlpFoldModel,
};
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
/// # Feature rows are projected AFTER this call — THE statement of why
/// Fold assignment ([`FoldDataset::get_fold`]) and labels are both POSITIONAL,
/// so projected row `i` has to come from `data[i]`. A frame built BEFORE the shuffle
/// compiles, comes out the right width, and runs to completion — it just pairs
/// every feature row with a different candidate's label and fold, which no type
/// or width check can notice and which surfaces only as a model scoring at
/// chance. `names` comes from the same lane walk as the values, so columns and
/// names align by construction either way, and the width assertions at the call
/// sites cannot see this at all.
///
/// Every rescorer below projects its lane rows after calling this and points
/// here instead of restating the argument; the rationale lived in five copies
/// and copies drift.
/// `mlp_streamed_rows_are_aligned_with_the_shuffled_data` is the test that
/// would catch an inversion.
fn canonicalize_and_shuffle(data: &mut [CompetedCandidate]) {
    data.sort_unstable_by_key(|c| {
        (
            c.scoring.identity.library_id,
            c.scoring.identity.precursor_charge,
        )
    });

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
/// THE canonical statement of leak-freedom — [`crossfit_lda`],
/// [`crossfit_score_column`] and therefore `rescore_lda` / `rescore_hybrid`
/// route through here and their docs just point back at this one.
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
/// # THE partition
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
/// # Failure policy — all-or-nothing
/// Returns `None` if ANY fold's fit or scoring fails (singular scatter / empty
/// class for the LDA, a fully culled train feature set for the MLP). Per-fold
/// fallback was rejected on purpose: leaving one fold's rows at 0.0 while the
/// others carry real discriminant values makes the score distribution
/// FOLD-DEPENDENT, i.e. a function of a row's position in the shuffle rather
/// than of its evidence. That silently corrupts the q-value ranking (a whole
/// third of the rows is pushed into an arbitrary slab of the sort) in a way that
/// no downstream check would catch. A uniform failure is louder and cannot
/// reintroduce fold identity as a splittable feature, so callers degrade the
/// WHOLE run instead. Callers log the failure; this function logs which fold
/// failed. "Uniform" does NOT mean the ranking survives — see
/// [`abort_standalone_mlp`] for what each caller does about that.
///
/// `scores` is local to this call and never leaves it through the `None` arm, so
/// no caller can obtain a partially filled column by ignoring the failure. The
/// `Some` arm has one other route to one, and it is not closed by scoping: a
/// [`FoldModel::predict`] that returned fewer scores than it was handed rows
/// would leave the tail of `held` at `0.0`, because the `zip` below stops at the
/// shorter side. All three impls return `rows.len()`, so there is no live bug;
/// the `debug_assert_eq!` is there so a fourth cannot quietly not. The same
/// hazard is a HARD assert one layer up in [`hybrid_frame`], where the shorter
/// side would truncate the matrix instead.
fn crossfit<D: FoldDataset, M: FoldModel>(
    data: &D,
    cfg: &M::Config,
    what: &str,
) -> Option<CrossFit<M>>
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
        let model = match M::fit(cfg, data, f, &train, &[]) {
            Ok(m) => m,
            Err(e) => {
                tracing::error!("cross-fit {what}: fold {f}/{n_folds} fit failed ({e})");
                return None;
            }
        };
        // Score ONLY the held-out fold, with a model that never saw it.
        let preds = match model.predict(data, &held) {
            Ok(p) => p,
            Err(e) => {
                tracing::error!("cross-fit {what}: fold {f}/{n_folds} scoring failed ({e})");
                return None;
            }
        };
        // The zip truncates to the shorter side, which for a short `preds` means
        // a PARTIALLY FILLED column — the exact shape the failure policy above
        // exists to prevent, and one no width check downstream would notice.
        debug_assert_eq!(
            preds.len(),
            held.len(),
            "cross-fit {what}: fold {f}'s model returned {} scores for {} held-out rows",
            preds.len(),
            held.len()
        );
        for (&i, s) in held.iter().zip(preds) {
            scores[i] = s;
        }
        fold_rows.push(held);
        models.push(model);
    }

    Some(CrossFit {
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
fn crossfit_lda<D: FoldDataset>(data: &D) -> Option<CrossFit<LdaModel>> {
    crossfit::<D, LdaModel>(data, &LdaConfig::default(), "LDA")
}

/// The extra COLUMN the hybrid appends to the nonlinear lane: each row's HELD-OUT
/// `M` score over the canonical rescore partition (see [`crossfit`]), or a
/// UNIFORMLY ZERO column if the cross-fit failed — plus which of the two it is.
///
/// `model` names the model in [`crossfit`]'s per-fold failure logs; `column`
/// names the frame column in this function's own.
///
/// Always `data.nrows()` long, and never partially filled — the fallback is a
/// fresh vector rather than whatever [`crossfit`] had accumulated before the
/// failing fold, so a run where folds 0 and 1 fitted and fold 2 did not still
/// yields all zeros. A partially filled column is the failure mode worth naming:
/// it makes a row's score depend on its position in the shuffle, gives the GBM a
/// feature it can split on to recover fold identity, and corrupts the q-value
/// ranking with nothing downstream to catch it. All-zero is constant -> zero
/// split gain -> the hybrid degrades to "GBM on the nonlinear lane".
///
/// Whether THAT degradation is benign is argued in exactly one place —
/// [`abort_standalone_mlp`] — and the short version is "only while the nonlinear
/// lane is itself trainable". [`report_untrained_folds`] is the runtime check the
/// callers run afterwards, and the flag returned here is what lets it tell "the
/// column degraded AND nothing trained" from "nothing trained at all".
///
/// The models are discarded: a hybrid's sidecar reports the GBM's importance over
/// `nonlinear + <column>`, and this model's own per-column weights belong to a
/// different feature set (the LINEAR lane) that the sidecar has no column for.
/// [`rescore_lda`] is the entry point that reports those.
fn crossfit_score_column<D: FoldDataset, M: FoldModel>(
    data: &D,
    cfg: &M::Config,
    model: &str,
    column: &str,
) -> (Vec<f64>, bool)
where
    M::Error: std::fmt::Display,
{
    match crossfit::<D, M>(data, cfg, model) {
        Some(cf) => (cf.scores, false),
        None => {
            tracing::error!(
                "hybrid: cross-fit {model} failed; {column} is uniformly 0, so the GBM now trains \
                 on the nonlinear lane ALONE. That is a weaker model, not a broken one — \
                 PROVIDED the nonlinear lane is trainable at this row count. If it is not, the \
                 GBM builds no trees and the q-values are meaningless; the check after the fit \
                 reports that case explicitly."
            );
            (vec![0.0f64; data.nrows()], true)
        }
    }
}

/// How many folds produced NO MODEL AT ALL: empty importance means no tree was
/// ever built, so every row that fold scored carries the same value.
///
/// A COUNT, and the callers act on `> 0` rather than on `== stats.len()`, because
/// the partial case is the worse one. One untrained fold out of three leaves a
/// third of the rows on real discriminant values and two thirds on a per-fold
/// constant: the score distribution then depends on a row's POSITION IN THE
/// SHUFFLE, which is exactly what [`crossfit`]'s own failure policy calls out as
/// silently corrupting the q-value ranking with nothing downstream to catch it.
/// A uniform failure is at least uniform.
///
/// No false positives on the healthy side: `GbmFoldModel::importance` reports
/// `NaN` for the columns forust never split on (dropped at the sidecar boundary)
/// and a finite value for the ones it did, `0.0` included, so a fold that built
/// even one split has at least one entry.
///
/// It does NOT see a fold that trained on noise: trees built off uninformative
/// columns report importance like any other. The claim is strictly "some fold
/// produced no model".
fn untrained_folds(stats: &RescoreFeatureStats) -> usize {
    stats
        .iter()
        .filter(|fs| fs.feature_importance.is_empty())
        .count()
}

/// Log the untrained-frame state, if that is what happened. `frame` names what
/// the GBM was handed; `degraded_column` is `Some(name)` only when a hybrid's
/// appended column is the all-zero fallback.
///
/// Called UNCONDITIONALLY by all three GBM-bearing rescorers, including the
/// default [`rescore`]. The check used to be gated behind the degraded flag, so
/// the path most users run never ran it — and a small real search (below
/// `GBMConfig::default`'s `min_leaf_weight` floor of roughly 34 rows per fold)
/// therefore emitted per-fold constants and q-values that rank by shuffle order
/// with nothing said about it. That is pre-existing behavior for [`rescore`], not
/// a regression, but the state is exactly the one this function calls "NOT usable
/// for FDR" and there is now machinery to detect it.
///
/// The two states stay DISTINGUISHABLE in the message, because they are different
/// diagnoses. "The fallback column degraded AND nothing trained" points at two
/// failures whose combination is what makes the run useless (see
/// [`crossfit_score_column`] and [`abort_standalone_mlp`] for why either alone is
/// survivable). "Nothing trained at all" points at the frame or the row count,
/// with no cross-fit involved — the only actionable thing there is more rows, a
/// lane with live variance, or a different `rescore_model`.
///
/// This is a detector rather than an abort because by this point the run has a
/// full result set, and the same condition is benign on a real search (27 live
/// nonlinear features, rows in the hundreds of thousands) while pathological on a
/// small or degenerate input. The caller's job is to make it VISIBLE.
fn report_untrained_folds(frame: &str, degraded_column: Option<&str>, stats: &RescoreFeatureStats) {
    let untrained = untrained_folds(stats);
    if untrained == 0 {
        return;
    }
    // Identical in both arms, and the sentence the operator has to act on.
    let unusable = "Every row those folds scored carries the SAME value, so its q-value ranks by \
                    the internal shuffle order rather than by evidence and is NOT usable for FDR.";
    let partial = if untrained < stats.len() {
        "The remaining folds DID train, so the score distribution is fold-dependent — worse than \
         uniformly degenerate, not better. "
    } else {
        ""
    };
    match degraded_column {
        Some(column) => tracing::error!(
            "rescore: {column} was zeroed by a failed cross-fit AND {untrained} of {} GBM folds \
             then built no trees on {frame}. {unusable} {partial}Treat this run as failed.",
            stats.len(),
        ),
        None => tracing::error!(
            "rescore: {untrained} of {} GBM folds built no trees on {frame}, with no failed \
             cross-fit involved — the frame itself is untrainable here, typically too few rows \
             per fold for GBMConfig's min_leaf_weight or a lane with no live variance. \
             {unusable} {partial}Treat this run as failed.",
            stats.len(),
        ),
    }
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
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
    scorer.fit().unwrap();

    let stats = scorer.feature_stats();
    // No cross-fit column here, so there is nothing to have degraded — but a
    // frame that trained nothing is still unusable, and this is the path most
    // users run. See `report_untrained_folds`.
    report_untrained_folds("the ALL lane", None, &stats);

    finalize(scorer.score(), stats)
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
pub fn rescore_lda(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    use std::time::Instant;

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

    let t = Instant::now();
    let (scores, stats) = match crossfit_lda(&dataset) {
        Some(cf) => {
            eprintln!(
                "  LDA: cross-fit ({N_RESCORE_FOLDS} folds) + scored {nrows} candidates in {:.2?}",
                t.elapsed()
            );
            let stats = cf.feature_stats(&dataset);
            (cf.scores, stats)
        }
        None => {
            // All-or-nothing (see `crossfit`). Here that means every score
            // stays 0.0, so the q-values collapse to a single uninformative
            // value — uniformly useless rather than silently mis-ranked.
            tracing::error!(
                "cross-fit LDA failed; ALL scores left at zero (uniform, so the ranking is \
                 degenerate rather than fold-dependent)"
            );
            (
                vec![0.0; nrows],
                vec![FoldStats {
                    fold: 0,
                    feature_stats: Vec::new(),
                    feature_importance: Vec::new(),
                }],
            )
        }
    };
    drop(dataset);
    for (cand, s) in data.iter_mut().zip(scores) {
        cand.assign_score(s);
    }

    finalize(data, stats)
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
/// Row-major, so "extra trailing column" is literally one more push at the end
/// of each row. Matrix and names are built by the same call precisely so the
/// appended value and the appended name cannot get out of step — the two
/// Keeping the value and name together prevents their width and order from
/// drifting apart.
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
/// CV on `nonlinear + lda_score`. The GBM re-sees `NONLINEAR_NCOLS + 1` features
/// (28 today) instead of the full ALL lane (128) — the compression play.
///
/// The compression buys real time and costs real sensitivity; neither trade-off
/// is constant across candidate counts.
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
pub fn rescore_hybrid(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
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

    // --- CROSS-FIT lda_score (leak-free), via the shared partition ---
    // Uniformly zero if the cross-fit failed, never partially filled, and the
    // See `crossfit_score_column` for the fallback's shape and
    // `abort_standalone_mlp` for when degrading is benign.
    // `LdaConfig::default()` is the same default `crossfit_lda` pins.
    let (lda_score, degraded) = crossfit_score_column::<_, LdaModel>(
        &lin_dataset,
        &LdaConfig::default(),
        "LDA",
        "lda_score",
    );

    let (precomputed, names) = hybrid_frame(&data, responses, "lda_score", lda_score);

    let mut scorer =
        CrossValidatedScorer::<CompetedCandidate, GbmFoldModel>::new_from_shuffled_with_precomputed(
            N_RESCORE_FOLDS,
            data,
            config,
            precomputed,
            names,
        );
    scorer.fit().unwrap();

    let stats = scorer.feature_stats();
    report_untrained_folds(
        "nonlinear + lda_score",
        degraded.then_some("lda_score"),
        &stats,
    );

    finalize(scorer.score(), stats)
}

/// FAILURE POLICY for the standalone MLP rescorer: abort the run, loudly and
/// with the fold-level diagnostic, rather than degrade.
///
/// Diverging (`-> !`) rather than returning a sentinel, because there is no
/// value a standalone rescorer could return here that means anything. This
/// differs from the hybrid LDA fallback:
///
/// * In [`rescore_hybrid`] a failed cross-fit zeroes ONE COLUMN of a frame that
///   still carries `NONLINEAR_NCOLS` real features. A constant column has zero
///   split gain, so the GBM ignores it and the run degrades to "GBM on the
///   nonlinear lane" — a weaker model, but one whose q-values still rank rows
///   by their evidence.
///
///   THIS HOLDS ONLY WHILE THE NONLINEAR LANE IS ITSELF TRAINABLE, and that
///   condition is load-bearing, not a caveat. If the GBM finds no admissible
///   split on the remaining columns — too few rows per fold against
///   `GBMConfig::default`'s `min_leaf_weight`, or a lane with no live variance —
///   it builds no trees, emits one constant score per fold, and the q-values are
///   exactly as meaningless as the standalone degradation this function rejects.
///   Both sides are measured and pinned by
///   `degraded_hybrid_ranks_only_when_the_rest_of_the_frame_is_trainable`, which
///   also records that `synthetic_competed` — the suite's main fixture — falls on
///   the WRONG side of it at every row count, because its nonlinear lane is
///   entirely constant.
///
///   What keeps the asymmetry defensible is therefore not that the hybrid
///   degradation is always benign — it is not — but that it is benign in the
///   regime a real search runs in (27 live nonlinear features, rows in the
///   hundreds of thousands), and that the way it fails LOUDEST is the way it is
///   most likely to fail: every GBM-bearing rescorer runs
///   [`report_untrained_folds`] after fitting and, if any fold produced no model,
///   logs that the run is not usable for FDR. An operator can see that in the log;
///   the standalone degradation would have been invisible.
///
///   The detector's reach is narrower than "we would notice a bad degraded run".
///   It sees a fold that built NO trees. It cannot see a fold that built trees on
///   columns carrying no real signal — that run is degraded, uninformative, and
///   silent. What bounds the damage there is not this check but the fact that
///   such a lane would have to be uninformative on a real search's worth of rows,
///   and that the sidecar reports the per-fold importance an operator can read.
/// * Here the failed fit IS the discriminant. Degrading to a uniform score makes
///   every row's score identical, so `assign_qval` walks a ranking that is
///   nothing but the seeded shuffle order and every q-value it emits is a
///   fiction. Nothing downstream can tell that apart from a real result: the
///   output file has the same columns, the same row count, and a q-value
///   distribution that looks plausible. That is strictly worse than not
///   finishing, because a search whose FDR is meaningless but unmarked gets used.
///
/// [`rescore_lda`] is the awkward case for this rule: its scores are also the
/// discriminant, and it degrades. It is deliberately left alone — it is a
/// shipping path, its failure condition is a different one (a singular scatter
/// matrix, not a culled lane), and re-deciding its policy is not what this task
/// was scoped to. The asymmetry to defend is not "MLP vs LDA" but "the score
/// itself vs one column of many", and on that axis this path and `rescore_lda`
/// belong on the same side. If `rescore_lda`'s degradation is ever revisited,
/// this is the argument to revisit it with.
///
/// A `Result` on the public rescorers was the other candidate and was rejected:
/// all four share one signature so [`crate::ml::rescore_with`] can dispatch on a
/// config value, the only caller is the CLI's search pipeline, and the only
/// decision that caller could make is "abort the run" — so the `Result` would
/// buy an unwrap at the call site and a signature split across the six.
///
/// # Three variants, three audiences
/// `scorer.fit()` surfaces ALL of [`MlpFoldError`]'s variants — the fit's
/// `NoUsableColumns`, and via `assign_scores` -> `get_scores` -> `predict` both
/// `WidthMismatch` and `NonFiniteScore`. They need different advice, so this
/// matches rather than printing one message for all three.
///
/// * `NoUsableColumns` is a statement about the INPUT (this lane is dead at this
///   row count), which a different `rescore_model` can work around.
/// * `WidthMismatch` is an internal invariant violation — the scorer predicts
///   against the very dataset it was constructed from, so the widths cannot
///   disagree unless the lane wiring is broken — and sending that operator to
///   their config would waste their time on a code defect.
/// * `NonFiniteScore` is a NUMERIC failure of the model itself, and the reason it
///   is an error at all is that its natural downstream abort names the wrong
///   thing: `assign_qval`'s only assertion is about sort order, so a `NaN`
///   discriminant used to abort the run blaming the SORT. See the variant's doc.
fn abort_standalone_mlp(nrows: usize, e: MlpFoldError) -> ! {
    let where_ = format!("all-feature lane ({ALL_NCOLS} columns, {nrows} candidates)");
    match e {
        MlpFoldError::NoUsableColumns { .. } => panic!(
            "MLP rescore aborted: {where_}: {e}\n\
             A standalone MLP rescorer has no usable fallback — degrading would score every \
             row identically, which makes every q-value meaningless while still looking like a \
             finished search. Re-run with a different `rescore_model` (e.g. `gbm`)."
        ),
        MlpFoldError::WidthMismatch { .. } => panic!(
            "MLP rescore aborted: {where_}: {e}\n\
             This is an INTERNAL INVARIANT VIOLATION, not a problem with the input or the \
             config: `rescore_mlp_with` hands the feature matrix to the scorer, and the scorer \
             predicts against that same dataset, so a width disagreement means the lane \
             wiring in `ml::qvalues` / `ml::mlp_fold` is broken. Changing `rescore_model` \
             will not help — please report this with the two widths above."
        ),
        MlpFoldError::NonFiniteScore { .. } => panic!(
            "MLP rescore aborted: {where_}: {e}\n\
             THE MLP'S RAW LOGIT IS THE DISCRIMINANT SCORE on this path, so a non-finite one \
             cannot be ranked: it would reach `assign_qval`, whose only assertion is that \
             scores are sorted in descending order, and abort there blaming the SORT for a \
             model failure. This is a numeric divergence in the fit, not a config error and \
             not a dead lane — lower `MlpConfig::lr`, or re-run with `gbm` / `lda`, whose \
             discriminants are finite by construction."
        ),
    }
}

/// The MLP rescorer body: the GBM's pipeline
/// ([`canonicalize_and_shuffle`] -> lane matrix -> [`CrossValidatedScorer`] ->
/// [`finalize`]) with [`MlpFoldModel`] swapped in for `GbmFoldModel`.
///
/// `config` is a parameter so tests can train a smaller net than
/// [`MlpConfig::default`] without tuning the production default for test runtime.
///
/// # Cross-fit semantics
/// UNCHANGED from the GBM path, because it is literally the same scorer — see
/// [`CrossValidatedScorer`] for the partition. Two MLP-specific notes on top of
/// it: the early-stopping fold (`f + 1`) IS used, by [`MlpFoldModel`]'s
/// patience rule exactly as the GBM uses it through forust's
/// `early_stopping_rounds` (unless the fold is under that fit's minimum
/// held-out size, which no production fold is — see `MIN_INNER_VAL_ROWS`), and
/// every fitted statistic (cull set,
/// standardization moments, imputation means, weights) is train-fold-only
/// inside [`MlpFoldModel::fit`]. The early-stopping fold reaches the loss
/// measurement and nothing else, and the scorer never asks a model to score
/// either the fold it trained on or the fold it stopped on.
///
/// # Determinism
/// [`canonicalize_and_shuffle`] pins the row order and therefore the positional
/// fold partition; each fold's RNG is `MlpConfig::rng_for_fold(fold)`, a pure
/// function of the config seed and the fold index. Nothing else is stochastic,
/// so two runs of the same build on the same input produce bit-identical
/// scores.
///
/// # Failure policy
/// A fold that cannot be fitted ABORTS the run — see [`abort_standalone_mlp`]
/// for why this path does not degrade the way the hybrid does.
fn rescore_mlp_with(
    mut data: Vec<CompetedCandidate>,
    config: MlpConfig,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    canonicalize_and_shuffle(&mut data);

    // Fold assignment remains positional after this shuffle. Feature values are
    // streamed one row at a time into the MLP's two reusable batch buffers;
    // there is deliberately no raw or transformed fold matrix on this path.
    let names = all_feature_name_set();
    let ncols = ALL_NCOLS;
    debug_assert_eq!(names.len(), ncols);

    let nrows = data.len();
    let mut scorer =
        CrossValidatedScorer::<CompetedCandidate, MlpFoldModel>::new_from_shuffled_streaming(
            N_RESCORE_FOLDS,
            data,
            config,
            names,
            write_competed_all_row,
        );
    if let Err(e) = scorer.fit_parallel() {
        abort_standalone_mlp(nrows, e);
    }

    let stats = scorer.feature_stats();

    finalize(scorer.score(), stats)
}

/// MLP rescorer over the ALL lane (linear ++ nonlinear) — the same feature set
/// [`rescore`] gives the GBM, so the two are directly comparable. The one
/// [`crate::ml::RescoreModel::Mlp`] selects.
///
/// See [`rescore_mlp_with`] for the cross-fit and determinism contracts.
///
/// Runtime and sensitivity comparisons are not constant across candidate counts.
pub fn rescore_mlp(data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
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

/// Append one row's LINEAR-lane values. `derived` is passed in rather than
/// recomputed so the all-lane walk pays for `Derived::compute` once per row.
/// Takes the blocks rather than a candidate so that both `CompetedCandidate`
/// (pre-rescore) and `FinalResult` (post-rescore) feed the SAME builder —
/// there is exactly one definition of a feature row.
fn push_linear_row(
    scoring: &ScoringFields,
    meta: &ResultMeta,
    derived: &Derived,
    out: &mut Vec<f64>,
) {
    out.extend_from_slice(&scoring.linear_feature_array());
    out.extend_from_slice(&meta.linear_feature_array());
    out.extend_from_slice(&derived.linear_feature_array());
    // sequence_counts has NO linear lane (context features).
}

/// Append one row's NONLINEAR-lane values (see [`push_linear_row`]).
fn push_nonlinear_row(
    scoring: &ScoringFields,
    meta: &ResultMeta,
    derived: &Derived,
    out: &mut Vec<f64>,
) {
    out.extend_from_slice(&scoring.nonlinear_feature_array());
    out.extend_from_slice(&meta.nonlinear_feature_array());
    out.extend_from_slice(&derived.nonlinear_feature_array());
    out.extend_from_slice(&sequence_counts::nonlinear_feature_array(
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
    let mut at = 0usize;
    macro_rules! write {
        ($values:expr) => {{
            let values = $values;
            let end = at + values.len();
            out[at..end].copy_from_slice(&values);
            at = end;
        }};
    }
    write!(scoring.linear_feature_array());
    write!(meta.linear_feature_array());
    write!(derived.linear_feature_array());
    write!(scoring.nonlinear_feature_array());
    write!(meta.nonlinear_feature_array());
    write!(derived.nonlinear_feature_array());
    write!(sequence_counts::nonlinear_feature_array(
        &scoring.identity.peptide
    ));
    debug_assert_eq!(at, ALL_NCOLS);
}

/// Matrix-free LINEAR-lane projection for the standalone and hybrid LDA paths.
fn write_competed_linear_row(candidate: &CompetedCandidate, out: &mut [f64]) {
    assert_eq!(out.len(), LINEAR_NCOLS);
    let scoring = &candidate.scoring;
    let meta = candidate.result_meta();
    let derived = Derived::compute(scoring);
    let scoring_values = scoring.linear_feature_array();
    let meta_values = meta.linear_feature_array();
    let derived_values = derived.linear_feature_array();
    let mut at = 0usize;
    for values in [
        scoring_values.as_slice(),
        meta_values.as_slice(),
        derived_values.as_slice(),
    ] {
        let end = at + values.len();
        out[at..end].copy_from_slice(values);
        at = end;
    }
    debug_assert_eq!(at, LINEAR_NCOLS);
}

/// The LINEAR-lane matrix for `data` in its CURRENT order (call AFTER any
/// shuffle, so row `i` aligns with `data[i]`). `LINEAR_NCOLS` wide.
#[cfg(test)]
fn build_linear_matrix(data: &[CompetedCandidate]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * LINEAR_NCOLS);
    for c in data {
        let meta = c.result_meta();
        push_linear_row(&c.scoring, &meta, &Derived::compute(&c.scoring), &mut out);
    }
    out
}

/// The NONLINEAR-lane matrix for `data`, `NONLINEAR_NCOLS` wide.
fn build_nonlinear_matrix(data: &[CompetedCandidate]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * NONLINEAR_NCOLS);
    for c in data {
        let meta = c.result_meta();
        push_nonlinear_row(&c.scoring, &meta, &Derived::compute(&c.scoring), &mut out);
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
        push_linear_row(s, &meta, &derived, &mut out);
        push_nonlinear_row(s, &meta, &derived, &mut out);
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

/// LINEAR-lane feature names (LDA), in [`push_linear_row`]'s order.
pub fn linear_feature_name_set() -> Vec<Arc<str>> {
    let mut n = NameSink::new();
    <ScoringFields as ScoreBlock>::linear_feature_names(&mut n);
    <ResultMeta as ScoreBlock>::linear_feature_names(&mut n);
    <Derived as ScoreBlock>::linear_feature_names(&mut n);
    n.into_names()
}

/// NONLINEAR-lane feature names, in [`push_nonlinear_row`]'s order. The
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
            vec![sample_competed_candidate_parsed()],
            vec![sample_competed_candidate_unparsed()],
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
            sample_competed_candidate_parsed(),
            sample_competed_candidate_unparsed(),
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

    /// The sequence block used to be gated off entirely for an unparsed
    /// peptide, changing the feature-set width. Now it is always present and
    /// contributes NaN — which forust reads as missing — so the width is
    /// label-independent and the "missing" signal is explicit.
    #[test]
    fn unparsed_sequence_emits_nan_not_a_narrower_row() {
        let names = nonlinear_feature_name_set();
        let seq_start = NONLINEAR_NCOLS - sequence_counts::LEN;
        assert_eq!(&*names[seq_start], "peptide_length");
        assert_eq!(&*names[NONLINEAR_NCOLS - 1], "peptide_n_mods");

        let unparsed = build_nonlinear_matrix(&[sample_competed_candidate_unparsed()]);
        assert!(
            unparsed[seq_start..].iter().all(|v| v.is_nan()),
            "unparsed sequence features must all be NaN: {:?}",
            &unparsed[seq_start..]
        );

        let parsed = build_nonlinear_matrix(&[sample_competed_candidate_parsed()]);
        // PEPTIDEK: 8 residues, no mods.
        assert_eq!(parsed[seq_start], 8.0);
        assert_eq!(parsed[NONLINEAR_NCOLS - 1], 0.0);
    }

    #[test]
    fn feature_names_are_unique() {
        // No exact counts here on purpose: adding or removing a score must not
        // break this test. Two features sharing a name would make the
        // importance/stat reports ambiguous. (That neither lane collapsed is
        // now a `const _: () = assert!(..)` next to the width consts.)
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
    /// EVERY field this varies lands in the LINEAR lane — `delta_group` and
    /// `delta_group_ratio` included. They are `ResultMeta`'s two LINEAR features
    /// and `ResultMeta::NONLINEAR_LEN == 0`, so this fixture's NONLINEAR lane is
    /// bit-constant at EVERY row count. The property is stated by MEASUREMENT,
    /// not by this comment: [`assert_nonlinear_lane_is_flat`] checks it and
    /// [`assert_nonlinear_lane_varies`] is its mirror, so a future feature that
    /// changes the answer fails a test rather than rotting a doc.
    ///
    /// The consequence is load-bearing and has been misread repeatedly: a GBM
    /// handed the nonlinear lane of this fixture has NOTHING to split on, so
    /// nothing here exercises "the model keeps working off the nonlinear lane".
    /// [`synthetic_competed_nonlinear_signal`] is the fixture that does.
    fn synthetic_competed(n: u32) -> Vec<CompetedCandidate> {
        (0..n)
            .map(|i| {
                let mut c = sample_competed_candidate_parsed();
                c.scoring.identity.library_id = i;
                let is_target = i % 2 == 0;
                c.scoring.identity.is_target = is_target;
                let base: u8 = if is_target { 20 } else { 8 };
                let jitter = (i % 5) as u8;
                c.scoring.counts.rising_cycles = base + jitter;
                c.scoring.counts.falling_cycles = base.saturating_sub(jitter);
                c.scoring.counts.npeaks = base + (i % 3) as u8;
                c.scoring.finalize_counts.n_scored_fragments = base + (i % 4) as u8;
                c.delta_group = if is_target { 2.0 } else { 0.5 } + (i % 7) as f32 * 0.1;
                c.delta_group_ratio = if is_target { 0.8 } else { 0.3 };
                c
            })
            .collect()
    }

    fn synthetic_competed_linear_only(n: u32) -> Vec<CompetedCandidate> {
        synthetic_competed(n)
            .into_iter()
            .map(|mut c| {
                let sample = sample_competed_candidate_parsed();
                c.delta_group = sample.delta_group;
                c.delta_group_ratio = sample.delta_group_ratio;
                c
            })
            .collect()
    }

    fn peptide_of_len(len: usize) -> Peptide {
        const ALPHABET: &[u8] = b"PEPTIDEKAVLGSR";
        let raw: String = (0..len)
            .map(|i| ALPHABET[i % ALPHABET.len()] as char)
            .collect();
        let residues = raw.bytes().map(AminoAcid::from_ascii).collect();
        Peptide {
            raw: Arc::from(raw.as_str()),
            parsed: Some(ParsedSequence {
                residues,
                mods: smallvec![],
            }),
            decoy: DecoyMarking::Target,
            decoy_group: 0,
        }
    }

    fn synthetic_competed_nonlinear_signal(n: u32) -> Vec<CompetedCandidate> {
        (0..n)
            .map(|i| {
                let is_target = i % 2 == 0;
                let len = if is_target { 16 } else { 7 } + (i % 3) as usize;
                let mut c = sample_competed_candidate_parsed();
                c.scoring.identity.peptide = peptide_of_len(len);
                c.scoring.identity.library_id = i;
                c.scoring.identity.is_target = is_target;
                c
            })
            .collect()
    }

    fn degraded_hybrid(
        mut data: Vec<CompetedCandidate>,
    ) -> (Vec<FinalResult>, RescoreFeatureStats) {
        canonicalize_and_shuffle(&mut data);
        let n = data.len();
        let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
        let (precomputed, names) = hybrid_frame(&data, responses, "lda_score", vec![0.0f64; n]);
        let mut scorer =
            CrossValidatedScorer::<CompetedCandidate, GbmFoldModel>::new_from_shuffled_with_precomputed(
                N_RESCORE_FOLDS,
                data,
                GBMConfig::default(),
                precomputed,
                names,
            );
        scorer.fit().unwrap();
        let stats = scorer.feature_stats();
        finalize(scorer.score(), stats)
    }

    fn assert_nonlinear_lane_varies(fixture: &[CompetedCandidate]) {
        let nl = build_nonlinear_matrix(fixture);
        let varying = nonlinear_feature_name_set()
            .iter()
            .enumerate()
            .filter(|(j, _)| {
                (0..fixture.len()).any(|i| {
                    let v = nl[i * NONLINEAR_NCOLS + j];
                    v.to_bits() != nl[*j].to_bits() && !(v.is_nan() && nl[*j].is_nan())
                })
            })
            .count();
        assert!(
            varying > 0,
            "fixture premise broken: EVERY nonlinear-lane column is constant, so the GBM \
             has nothing to split on and 'the GBM keeps working off the rest of the frame' \
             cannot be observed here"
        );
    }

    #[test]
    fn degraded_hybrid_ranks_only_when_the_rest_of_the_frame_is_trainable() {
        let trainable = synthetic_competed_nonlinear_signal(360);
        assert_nonlinear_lane_varies(&trainable);
        let (out, stats) = degraded_hybrid(trainable);
        assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);
        for fs in &stats {
            assert!(
                !fs.feature_importance.is_empty(),
                "fold {}: the degraded frame built no trees even though the nonlinear lane \
                 is trainable, so the hybrid's whole justification for degrading is false",
                fs.fold
            );
        }
        let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
        let is_target: Vec<bool> = out.iter().map(|r| r.scoring.identity.is_target).collect();
        let auc = pair_auc(&scores, &is_target);
        assert!(
            auc > 0.8,
            "AUC {auc}: with the appended column zeroed the GBM must still rank on the \
             nonlinear lane. If this is at chance, degrading is NOT a weaker model — it is \
             a meaningless one, and the hybrid should abort like the standalone path does"
        );
        assert_eq!(
            untrained_folds(&stats),
            0,
            "the runtime detector must not fire when every fold of the degraded frame trained"
        );

        let flat = synthetic_competed(360);
        assert_nonlinear_lane_is_flat(&flat, "lda_score");
        let (out, stats) = degraded_hybrid(flat);
        assert!(
            stats.iter().all(|fs| fs.feature_importance.is_empty()),
            "measured: on this fixture the degraded frame builds no trees at all. If that \
             has changed, the boundary this test documents has moved and the doc on \
             `abort_standalone_mlp` needs re-checking"
        );
        assert_eq!(
            untrained_folds(&stats),
            N_RESCORE_FOLDS as usize,
            "the runtime detector must fire on exactly this state, or an operator gets \
             meaningless q-values with nothing in the log"
        );
        let distinct: std::collections::HashSet<u32> =
            out.iter().map(|r| r.discriminant_score.to_bits()).collect();
        assert!(
            distinct.len() <= N_RESCORE_FOLDS as usize,
            "no trees means one constant score per fold; got {} distinct scores",
            distinct.len()
        );
        let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
        let is_target: Vec<bool> = out.iter().map(|r| r.scoring.identity.is_target).collect();
        let auc = pair_auc(&scores, &is_target);
        assert!(
            (auc - 0.5).abs() < 0.05,
            "a per-fold constant cannot rank better than chance; AUC {auc}"
        );
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
        let (out_a, stats_a) = rescore_lda(synthetic_competed(n));
        let (out_b, _stats_b) = rescore_lda(synthetic_competed(n));

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

        let key = |out: &[FinalResult]| -> Vec<(u32, u32)> {
            let mut v: Vec<(u32, u32)> = out
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
        let (_out, stats) = rescore_lda(synthetic_competed(90));
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
        let key = |out: &[FinalResult]| -> Vec<(u32, u32)> {
            let mut v: Vec<(u32, u32)> = out
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
            let (out_a, stats_a) = run();
            let (out_b, _) = run();

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
            let (out, _) = rescore_mlp_with(synthetic_competed(120), mlp_test_cfg(seed));
            let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
            let is_target: Vec<bool> = out.iter().map(|r| r.scoring.identity.is_target).collect();
            let auc = pair_auc(&scores, &is_target);
            assert!(
                auc > 0.9,
                "seed {seed}: AUC {auc} is near chance. The cause this test was \
                     written for is a feature frame built BEFORE the shuffle rather than \
                     after it, which pairs every feature row with another candidate's \
                     label and fold — but any MLP training regression lands here too. \
                     Check `mlp_fold`'s `fitted_model_separates_a_separable_toy_on_held_\
                     out_rows` first: if that still passes, the fit is fine and the \
                     ordering is the problem"
            );
        }
    }

    #[test]
    fn rescore_with_dispatches_every_variant() {
        use crate::ml::{
            RescoreModel,
            rescore_with,
        };

        const N: u32 = 360;
        type Rescorer = fn(Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats);
        let key = |out: &[FinalResult]| -> Vec<(u32, u32)> {
            let mut v: Vec<(u32, u32)> = out
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
            RescoreModel::default(),
            RescoreModel::Mlp,
            "the library default must match the shipped configuration"
        );

        let mut seen: Vec<(RescoreModel, Vec<(u32, u32)>)> = Vec::new();
        for (variant, direct, frame_width) in [
            (RescoreModel::Gbm, rescore as Rescorer, ALL_NCOLS),
            (RescoreModel::Lda, rescore_lda as Rescorer, LINEAR_NCOLS),
            (
                RescoreModel::Hybrid,
                rescore_hybrid as Rescorer,
                NONLINEAR_NCOLS + 1,
            ),
            (RescoreModel::Mlp, rescore_mlp as Rescorer, ALL_NCOLS),
        ] {
            let (dispatched, stats) = rescore_with(variant, synthetic_competed(N));
            assert_eq!(dispatched.len(), N as usize);
            for r in &dispatched {
                assert!(
                    r.discriminant_score.is_finite(),
                    "{variant:?}: non-finite discriminant score {}",
                    r.discriminant_score
                );
                assert!(
                    (0.0..=1.0).contains(&r.qvalue),
                    "{variant:?}: qvalue out of [0,1]: {}",
                    r.qvalue
                );
            }

            assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);
            for fs in &stats {
                assert_eq!(
                    fs.feature_stats.len(),
                    frame_width,
                    "rescore_with({variant:?}) fold {}: trained on a {}-column frame, not the \
                     {frame_width} this variant's lane has",
                    fs.fold,
                    fs.feature_stats.len(),
                );
            }

            let (direct_out, _) = direct(synthetic_competed(N));
            assert_eq!(
                key(&dispatched),
                key(&direct_out),
                "rescore_with({variant:?}) does not reproduce its own entry point, so the \
                 arm dispatches to a different model"
            );

            seen.push((variant, key(&dispatched)));
        }

        for i in 0..seen.len() {
            for j in (i + 1)..seen.len() {
                assert_ne!(
                    seen[i].1, seen[j].1,
                    "{:?} and {:?} produced bit-identical scores, so (b) would hold just as \
                     well for a permutation of the arms",
                    seen[i].0, seen[j].0
                );
            }
        }
    }

    #[test]
    fn rescore_hybrid_smoke_and_determinism() {
        let n = 360;
        let (out_a, stats_a) = rescore_hybrid(synthetic_competed(n));
        let (out_b, _stats_b) = rescore_hybrid(synthetic_competed(n));

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

        let key = |out: &[FinalResult]| -> Vec<(u32, u32)> {
            let mut v: Vec<(u32, u32)> = out
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
                let mut c = sample_competed_candidate_parsed();
                c.scoring.identity.library_id = i;
                c.scoring.identity.is_target = i % 2 == 0;
                c
            })
            .collect()
    }

    #[test]
    #[should_panic(expected = "MLP rescore aborted: all-feature lane")]
    fn standalone_mlp_aborts_rather_than_scoring_every_row_the_same() {
        rescore_mlp_with(indistinguishable_competed(24), mlp_test_cfg(7));
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

        let (out, stats) = rescore_hybrid(synthetic_competed_linear_only(360));

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
        // mobility-derived GBM feature must become NaN (forust missing), so it
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
        let before = build_all_matrix(competed_rows(&[sample_competed_candidate_parsed()]));

        let mut cand = sample_competed_candidate_parsed();
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
    /// rather than fail. Rows carry distinct `delta_group` values so the
    /// assertion can tell them apart; a length check alone passes even when
    /// the layout is wrong.
    #[test]
    fn feature_frame_rows_are_contiguous() {
        let rows: Vec<_> = [1.0f32, 2.0, 3.0]
            .into_iter()
            .map(|delta_group| {
                let mut c = sample_competed_candidate_parsed();
                c.delta_group = delta_group;
                c.into_final()
            })
            .collect();

        let (names, got) = feature_frame(&rows);
        let nf = names.len();
        assert_eq!(got.len(), rows.len() * nf, "one value per name per row");

        let j = names
            .iter()
            .position(|n| &**n == "delta_group")
            .expect("delta_group is an ALL-lane feature");
        for (i, expected) in [1.0f64, 2.0, 3.0].into_iter().enumerate() {
            assert_eq!(
                got[i * nf + j],
                expected,
                "row {i}'s delta_group is not at matrix[{i} * {nf} + {j}]"
            );
        }
    }
}
