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
    RowMajorDataset,
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
/// Shared by ALL SIX rescorers so the row order — and therefore the positional
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
/// # Lane matrices are built AFTER this call — THE statement of why
/// Fold assignment ([`FoldDataset::get_fold`]) and labels are both POSITIONAL,
/// so matrix row `i` has to be `data[i]`. A matrix built BEFORE the shuffle
/// compiles, comes out the right width, and runs to completion — it just pairs
/// every feature row with a different candidate's label and fold, which no type
/// or width check can notice and which surfaces only as a model scoring at
/// chance. `names` comes from the same lane walk as the values, so columns and
/// names align by construction either way, and the width assertions at the call
/// sites cannot see this at all.
///
/// Every rescorer below builds its lane matrix after calling this and points
/// here instead of restating the argument; the rationale lived in five copies
/// and copies drift.
/// `mlp_lane_matrix_is_row_aligned_with_the_shuffled_data` is the test that
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

/// The shared tail of ALL SIX rescorers: rank by discriminant score, derive
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
/// [`crossfit_score_column`] and therefore `rescore_lda` / `rescore_hybrid` /
/// `rescore_hybrid_mlp` all route through here and their docs just point back
/// at this one.
///
/// Generic over the model because the partition is the whole content of this
/// function and it does not vary with what is being fitted: the LDA path and the
/// MLP hybrid's `mlp_score` column want the SAME walk, and writing it twice
/// would be two chances to get the fold arithmetic subtly different in a way
/// only an optimistic FDR would reveal. `what` names the model in the failure
/// logs, which is the only thing the two calls differ in.
///
/// # Why a row may never be scored by a model that saw it
/// These models are label-aware: they fit on the target/decoy labels. One
/// in-sample fit over all rows, scoring those same rows, lets every row's
/// discriminant peek at its own label; the target/decoy separation then looks
/// better than it is, and `assign_qval` derives the q-values from exactly that
/// separation. The result is an FDR that is wrong in the flattering direction
/// with nothing downstream to catch it.
///
/// The hybrids need this even more sharply: there the held-out score is not the
/// final score but ONE COLUMN (`lda_score` / `mlp_score`) fed to a
/// cross-validated GBM. An in-sample column would smuggle a row's own label
/// into a feature the GBM reads while that row is held out of its GBM fold — so
/// the GBM's own CV cannot notice, and the leak surfaces only as an optimistic
/// FDR.
///
/// # THE partition
/// Fold membership comes from [`FoldDataset::get_fold`] — for the
/// [`RowMajorDataset`] every caller builds, that is `i % N_RESCORE_FOLDS`. For
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
/// from [`RowMajorDataset::get_fold`], so it has one definition; the
/// `crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on` test pins the
/// two partitions against each other rather than against a re-typed modulo.
///
/// # Failure policy — all-or-nothing
/// Returns `None` if ANY fold's fit or scoring fails (singular scatter / empty
/// class for the LDA, a fully culled train matrix for the MLP). Per-fold
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

/// The extra COLUMN a hybrid appends to the nonlinear lane: each row's HELD-OUT
/// `M` score over the canonical rescore partition (see [`crossfit`]), or a
/// UNIFORMLY ZERO column if the cross-fit failed — plus which of the two it is.
///
/// ONE function for BOTH hybrids' degradation arms, so `lda_score` and
/// `mlp_score` cannot drift apart. They used to differ in nothing but the model
/// name and the column name, with a byte-identical five-line prose block in each,
/// and only one of the two was even factored out of its caller.
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
/// [`rescore_lda`] and [`rescore_mlp_linear`] are the entry points that report
/// those.
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
/// A uniform failure is at least uniform. An earlier version of this used `all`
/// and was therefore blind to the state it most needed to report — measured on
/// `degraded_hybrid(synthetic_competed_nonlinear_signal(n))`: `n = 198` gives
/// `[0, 0, 0]`, `n = 204` gives `[1, 0, 0]`, `n = 208` gives `[1, 0, 1]`.
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

/// Sage-style shrinkage-LDA rescorer: closed-form linear fits, no boosting — so
/// ~20x cheaper than the GBM path, measured as Phase 5 wall time on a
/// 114138-candidate mzML run (0.48s against the GBM's 9.59s). One dataset, one
/// machine, and that run was FAIMS so a quarter of the linear lane was culled as
/// dead: the order of magnitude is the claim, not the digits. The figure used to
/// read "~100x" with nothing said about what it was measured against, which is
/// exactly why nobody noticed it going stale.
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
/// ([`crate::ml::RescoreModel::Lda`]); the GBM `rescore` remains the default.
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

    // Materialize the linear-lane matrix once as a single flat row-major buffer
    // (`feat[i*ncols + j]`) — the layout LDA wants.
    let t = Instant::now();
    let feat = build_linear_matrix(&data);
    debug_assert_eq!(feat.len(), nrows * ncols);
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();

    // Optional raw-matrix dump for offline feature-engineering iteration.
    // `TIMSSEEK_LDA_DUMP=/prefix` writes `<prefix>.f64` (row-major matrix),
    // `<prefix>.labels` (u8, 1=target), `<prefix>.names.txt`. This is the
    // LINEAR-lane matrix + linear names (offline python reads this).
    if let Ok(prefix) = std::env::var("TIMSSEEK_LDA_DUMP") {
        let is_decoy: Vec<bool> = responses.iter().map(|&y| y < 0.5).collect();
        dump_feature_matrix(&prefix, &feat, nrows, &names, &is_decoy);
    }
    eprintln!(
        "  LDA: extracted {nrows} x {ncols} linear-lane matrix in {:.2?}",
        t.elapsed()
    );

    // The matrix, its names and the fold partition become one dataset from here
    // on; the cross-fit and the sidecar stats both read rows through it, so
    // neither can disagree with the other about a column or a fold.
    let dataset = RowMajorDataset::new(
        PrecomputedFeatures::from_row_major(feat, ncols, responses),
        names,
        N_RESCORE_FOLDS as usize,
    );

    let t = Instant::now();
    let stats = match crossfit_lda(&dataset) {
        Some(cf) => {
            for (cand, &s) in data.iter_mut().zip(cf.scores.iter()) {
                cand.assign_score(s);
            }
            eprintln!(
                "  LDA: cross-fit ({N_RESCORE_FOLDS} folds) + scored {nrows} candidates in {:.2?}",
                t.elapsed()
            );
            cf.feature_stats(&dataset)
        }
        None => {
            // All-or-nothing (see `crossfit`). Here that means every score
            // stays 0.0, so the q-values collapse to a single uninformative
            // value — uniformly useless rather than silently mis-ranked.
            tracing::error!(
                "cross-fit LDA failed; ALL scores left at zero (uniform, so the ranking is \
                 degenerate rather than fold-dependent)"
            );
            vec![FoldStats {
                fold: 0,
                feature_stats: Vec::new(),
                feature_importance: Vec::new(),
            }]
        }
    };

    finalize(data, stats)
}

/// The LINEAR-lane dataset a hybrid cross-fits its extra column over.
///
/// `data` must already be through [`canonicalize_and_shuffle`] and `responses`
/// must come from the same order.
///
/// The fold count is `N_RESCORE_FOLDS` on a [`RowMajorDataset`], i.e. the SAME
/// `get_fold` the GBM's [`CrossValidatedScorer`] derives its partition from
/// below. That shared definition is what makes "the same fold assignment on both
/// sides" structural rather than a comment; see [`crossfit`].
fn hybrid_linear_dataset(data: &[CompetedCandidate], responses: Vec<f64>) -> RowMajorDataset {
    let lane = Lane::Linear;
    let feat = lane.matrix(data);
    debug_assert_eq!(feat.len(), data.len() * lane.ncols());
    RowMajorDataset::new(
        PrecomputedFeatures::from_row_major(feat, lane.ncols(), responses),
        lane.names(),
        N_RESCORE_FOLDS as usize,
    )
}

/// The frame BOTH hybrids hand the GBM: the NONLINEAR lane matrix with one
/// cross-fit score appended as the TRAILING column, plus the matching names.
///
/// Row-major, so "extra trailing column" is literally one more push at the end
/// of each row. Matrix and names are built by the same call precisely so the
/// appended value and the appended name cannot get out of step — the two
/// hybrids differ in nothing here but the string, and a copy per hybrid would be
/// a second place for the width, the walk and the name to drift apart.
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
/// (28 today) instead of the full ALL lane (128) — the compression play — at
/// ~parity.
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

    // BOTH lane matrices are built from `data` in its post-shuffle order (see
    // `canonicalize_and_shuffle`), so row `i` is the same candidate in lin, nl,
    // lda_score, responses, and the moved `data`.
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let lin_dataset = hybrid_linear_dataset(&data, responses.clone());

    // --- CROSS-FIT lda_score (leak-free), via the shared partition ---
    // Uniformly zero if the cross-fit failed, never partially filled, and the
    // same helper the MLP hybrid uses — see `crossfit_score_column` for the
    // fallback's shape and `abort_standalone_mlp` for when degrading is benign.
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

/// WHICH lane feature set a rescorer trains on — the one axis the MLP path is
/// parameterized over.
///
/// Each variant bundles the three things that have to agree for a lane to be
/// coherent: the width const, the name set, and the matrix builder. They are
/// selected together here rather than at each call site precisely because
/// picking two of the three from one lane and the third from another is the
/// silent-misattribution bug `lane_matrix_widths_match_name_sets` exists to
/// catch, and the enum removes the opportunity.
///
/// Deliberately NOT public: the lane is an internal knob of the MLP rescorers
/// (see [`rescore_mlp_linear`] / [`rescore_mlp_all`], each of which pins one
/// variant), not part of the model selection surface.
#[derive(Debug, Clone, Copy)]
enum Lane {
    /// Fields approx-Gaussian after their declared per-row transform — what the
    /// LDA trains on.
    Linear,
    /// Linear ++ nonlinear, the full feature set the GBM trains on.
    All,
}

impl Lane {
    /// Row width. A compile-time const per lane (see [`LINEAR_NCOLS`]).
    fn ncols(self) -> usize {
        match self {
            Lane::Linear => LINEAR_NCOLS,
            Lane::All => ALL_NCOLS,
        }
    }

    /// Column names, in the matrix's column order.
    fn names(self) -> Vec<Arc<str>> {
        match self {
            Lane::Linear => linear_feature_name_set(),
            Lane::All => all_feature_name_set(),
        }
    }

    /// The row-major lane matrix for `data` IN ITS CURRENT ORDER — so call it
    /// AFTER [`canonicalize_and_shuffle`], never before.
    fn matrix(self, data: &[CompetedCandidate]) -> Vec<f64> {
        match self {
            Lane::Linear => build_linear_matrix(data),
            Lane::All => build_all_matrix(competed_rows(data)),
        }
    }
}

/// FAILURE POLICY for the STANDALONE MLP rescorers: abort the run, loudly and
/// with the fold-level diagnostic, rather than degrade.
///
/// Diverging (`-> !`) rather than returning a sentinel, because there is no
/// value a standalone rescorer could return here that means anything. This
/// deliberately DIFFERS from what the same failure does in the hybrids, and the
/// asymmetry is the whole point:
///
/// * In [`rescore_hybrid_mlp_with`] a failed cross-fit zeroes ONE COLUMN of a
///   frame that still carries `NONLINEAR_NCOLS` real features. A constant column
///   has zero split gain, so the GBM ignores it and the run degrades to "GBM on
///   the nonlinear lane" — a weaker model, but one whose q-values still rank rows
///   by their evidence. Same for `lda_score` in [`rescore_hybrid`].
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
/// all six share one signature so [`crate::ml::rescore_with`] can dispatch on a
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
fn abort_standalone_mlp(lane: Lane, ncols: usize, nrows: usize, e: MlpFoldError) -> ! {
    let where_ = format!("{lane:?} lane ({ncols} columns, {nrows} candidates)");
    match e {
        MlpFoldError::NoUsableColumns { .. } => panic!(
            "MLP rescore aborted: {where_}: {e}\n\
             A standalone MLP rescorer has no usable fallback — degrading would score every \
             row identically, which makes every q-value meaningless while still looking like a \
             finished search. Re-run with a different `rescore_model` (e.g. `gbm`, or the \
             `hybrid_mlp` variant, which tolerates this failure by dropping one column)."
        ),
        MlpFoldError::WidthMismatch { .. } => panic!(
            "MLP rescore aborted: {where_}: {e}\n\
             This is an INTERNAL INVARIANT VIOLATION, not a problem with the input or the \
             config: `rescore_mlp_lane` hands the lane matrix to the scorer, and the scorer \
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

/// The body of BOTH MLP rescorers: the GBM's pipeline
/// ([`canonicalize_and_shuffle`] -> lane matrix -> [`CrossValidatedScorer`] ->
/// [`finalize`]) with [`MlpFoldModel`] swapped in for `GbmFoldModel`.
///
/// `lane` and `config` are parameters rather than two near-identical copies of
/// this function: the lane changes only which of the three matched constants
/// [`Lane`] hands back, and `config` is what lets tests train a smaller net than
/// [`MlpConfig::default`] without the default itself being tuned for test
/// runtime.
///
/// # Cross-fit semantics
/// UNCHANGED from the GBM path, because it is literally the same scorer — see
/// [`CrossValidatedScorer`] for the partition. Two MLP-specific notes on top of
/// it: the early-stopping fold (`f + 1`) IS used, by [`MlpFoldModel`]'s
/// patience rule exactly as the GBM uses it through forust's
/// `early_stopping_rounds`, and every fitted statistic (cull set,
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
/// for why this path does not degrade the way the hybrids do.
fn rescore_mlp_lane(
    mut data: Vec<CompetedCandidate>,
    lane: Lane,
    config: MlpConfig,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    canonicalize_and_shuffle(&mut data);

    // Built here, AFTER the shuffle, over `data` in its current order — see
    // `canonicalize_and_shuffle` for why that ordering is load-bearing. The
    // three matched constants come from `lane` together, so a lane matrix can
    // never be paired with another lane's width or names.
    let names = lane.names();
    let ncols = lane.ncols();
    debug_assert_eq!(names.len(), ncols);
    let feat = lane.matrix(&data);
    debug_assert_eq!(feat.len(), data.len() * ncols);

    let nrows = data.len();
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let precomputed = PrecomputedFeatures::from_row_major(feat, ncols, responses);

    let mut scorer =
        CrossValidatedScorer::<CompetedCandidate, MlpFoldModel>::new_from_shuffled_with_precomputed(
            N_RESCORE_FOLDS,
            data,
            config,
            precomputed,
            names,
        );
    if let Err(e) = scorer.fit() {
        abort_standalone_mlp(lane, ncols, nrows, e);
    }

    let stats = scorer.feature_stats();

    finalize(scorer.score(), stats)
}

// BOTH entry points NAME THEIR LANE, and there is deliberately no bare
// `rescore_mlp` for either to be mistaken for.
//
// The lane is not a detail an entry point may leave implicit: it selects the
// feature set the model trains on, nothing type-checks which one a caller
// meant, and a wrong pick still compiles, still runs, and still produces
// plausible q-values. The concrete hazard is the model-selection arm a later
// task writes — `RescoreModel::Mlp => rescore_mlp(data)` reads as correct
// against ANY mapping, so a bare name would make the inversion invisible at
// exactly the site that decides it. Spelling the lane in the name is what
// forces that arm to state which one it wants.
//
// The realized mapping, pinned by
// `rescore_with_dispatches_each_mlp_variant_to_its_lane`:
//   `RescoreModel::Mlp`    -> `rescore_mlp_linear` (the default MLP shape)
//   `RescoreModel::MlpAll` -> `rescore_mlp_all`

/// MLP rescorer over the LINEAR lane only — the same feature set
/// [`rescore_lda`] uses, so it isolates "nonlinear model" from "more features".
/// THE DEFAULT MLP SHAPE, i.e. the one [`crate::ml::RescoreModel::Mlp`] selects.
///
/// See [`rescore_mlp_lane`] for the cross-fit and determinism contracts.
pub fn rescore_mlp_linear(data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    rescore_mlp_lane(data, Lane::Linear, MlpConfig::default())
}

/// MLP rescorer over the ALL lane (linear ++ nonlinear) — the same feature set
/// [`rescore`] gives the GBM, so the two are directly comparable. The one
/// [`crate::ml::RescoreModel::MlpAll`] selects.
///
/// See [`rescore_mlp_lane`] for the cross-fit and determinism contracts.
pub fn rescore_mlp_all(data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    rescore_mlp_lane(data, Lane::All, MlpConfig::default())
}

/// MLP counterpart of [`rescore_hybrid`]: cross-fit the MLP on the LINEAR lane,
/// push its (leak-free) `mlp_score` as one extra column into the NONLINEAR lane,
/// then train the GBM CV on `nonlinear + mlp_score`. Same shape as the LDA
/// hybrid, with a nonlinear compressor in place of the linear one, so the two are
/// directly comparable.
///
/// `mlp_config` is a parameter for the same reason [`rescore_mlp_lane`]'s is:
/// tests train a smaller net than [`MlpConfig::default`] without the default
/// being tuned for test runtime. The GBM half uses [`GBMConfig::default`], as
/// [`rescore_hybrid`] does.
///
/// # Leak-freedom
/// `mlp_score` is NOT the output — it is a FEATURE the GBM reads, which makes it
/// the sharpest leak surface here. If a row's `mlp_score` came from an MLP that
/// trained on that row, the row's own label rides into a column the GBM sees
/// while the GBM is holding that row out of its own fold: the GBM's
/// cross-validation structurally cannot detect it, and it surfaces only as an
/// optimistically wrong FDR. [`crossfit`] is therefore the only producer of this
/// column, and the fold ASSIGNMENT it walks is the same
/// [`RowMajorDataset::get_fold`] the [`CrossValidatedScorer`] below derives its
/// own partition from.
///
/// Deliberately NOT public: the config parameter is a test knob, and
/// [`rescore_hybrid_mlp`] is the entry point.
fn rescore_hybrid_mlp_with(
    mut data: Vec<CompetedCandidate>,
    mlp_config: MlpConfig,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    let config = GBMConfig::default();

    // Canonical sort + seeded shuffle — IDENTICAL to every other rescorer (same
    // helper, same key, same seed).
    canonicalize_and_shuffle(&mut data);

    // BOTH lane matrices come from `data` in its post-shuffle order (see
    // `canonicalize_and_shuffle`), so row `i` is the same candidate in the linear
    // dataset, the nonlinear matrix, `mlp_score`, `responses`, and the moved
    // `data`.
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let lin_dataset = hybrid_linear_dataset(&data, responses.clone());

    // Uniformly zero if the cross-fit failed, never partially filled — see
    // `crossfit_score_column`, which is also the LDA hybrid's arm.
    let (mlp_score, degraded) =
        crossfit_score_column::<_, MlpFoldModel>(&lin_dataset, &mlp_config, "MLP", "mlp_score");

    let (precomputed, names) = hybrid_frame(&data, responses, "mlp_score", mlp_score);

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
        "nonlinear + mlp_score",
        degraded.then_some("mlp_score"),
        &stats,
    );

    finalize(scorer.score(), stats)
}

/// The MLP hybrid on [`MlpConfig::default`] — see [`rescore_hybrid_mlp_with`].
///
/// The one [`crate::ml::RescoreModel::HybridMlp`] selects.
#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore_hybrid_mlp(data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    rescore_hybrid_mlp_with(data, MlpConfig::default())
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
        // Header: nrows, ncols as u64 little-endian, then row-major f64, also
        // little-endian. The reader is offline python, which assumes LE.
        f.write_all(&(nrows as u64).to_le_bytes())?;
        f.write_all(&(names.len() as u64).to_le_bytes())?;
        for v in base {
            f.write_all(&v.to_le_bytes())?;
        }
        f.flush()?;

        let labels: Vec<u8> = is_decoy.iter().map(|&d| u8::from(!d)).collect();
        std::fs::write(format!("{prefix}.labels"), &labels)?;
        let names_txt = names
            .iter()
            .map(|s| s.as_ref())
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(format!("{prefix}.names.txt"), names_txt)?;
        Ok(())
    };
    match write() {
        Ok(()) => eprintln!("  LDA: dumped feature matrix to {prefix}.{{f64,labels,names.txt}}"),
        Err(e) => tracing::error!("feature dump failed: {e}"),
    }
}

// ---------------------------------------------------------------------------
// Lane feature matrices — the ML consumer's live path
// ---------------------------------------------------------------------------
//
// A lane matrix is a flat row-major `Vec<f64>` (`feat[i * ncols + j]`), the
// layout `LdaModel::fit_matrix` and `PrecomputedFeatures::from_row_major` both index
// directly. `ncols` is not carried alongside it: every contributing block's
// width is an inherent const, so the lane widths below are compile-time
// constants.
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

/// The LINEAR-lane matrix for `data` in its CURRENT order (call AFTER any
/// shuffle, so row `i` aligns with `data[i]`). `LINEAR_NCOLS` wide.
fn build_linear_matrix(data: &[CompetedCandidate]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * LINEAR_NCOLS);
    for c in data {
        let meta = c.result_meta();
        push_linear_row(&c.scoring, &meta, &Derived::compute(&c.scoring), &mut out);
    }
    out
}

/// The NONLINEAR-lane matrix for `data`, `NONLINEAR_NCOLS` wide (see
/// [`build_linear_matrix`] for the ordering contract).
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
// feature values reach the scorer as a prebuilt lane frame
// (`PrecomputedFeatures::from_row_major`), so neither type here implements
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
                // Give the linear lane (counts) label-correlated variance.
                let base: u8 = if is_target { 20 } else { 8 };
                let jitter = (i % 5) as u8;
                c.scoring.counts.rising_cycles = base + jitter;
                c.scoring.counts.falling_cycles = base.saturating_sub(jitter);
                c.scoring.counts.npeaks = base + (i % 3) as u8;
                c.scoring.finalize_counts.n_scored_fragments = base + (i % 4) as u8;
                // The two `ResultMeta` fields. Both are LINEAR-lane features
                // (all-lane columns 90 and 91, inside `LINEAR_NCOLS`) — this is
                // MORE linear signal, not nonlinear signal. See the doc above.
                c.delta_group = if is_target { 2.0 } else { 0.5 } + (i % 7) as f32 * 0.1;
                c.delta_group_ratio = if is_target { 0.8 } else { 0.3 };
                c
            })
            .collect()
    }

    /// [`synthetic_competed`] with its two `ResultMeta` fields reset to the
    /// sample defaults, which removes LINEAR-lane signal — NOT nonlinear signal.
    ///
    /// `delta_group` / `delta_group_ratio` are all-lane columns 90 and 91, i.e.
    /// inside `LINEAR_NCOLS`, so this fixture's nonlinear lane is bit-constant
    /// for exactly the reason `synthetic_competed`'s already is: nothing either
    /// one varies reaches that lane. What resetting them buys is a NARROWER
    /// linear lane — only the four count fields vary — so a hybrid's cross-fit
    /// score column is the sole channel by which any signal at all can reach the
    /// GBM.
    ///
    /// Built for the hybrid tests, where "the GBM cannot separate these rows
    /// without the appended column" has to be a fact about the fixture rather
    /// than a hope. That the nonlinear lane is flat is stated by MEASUREMENT at
    /// the point of use ([`assert_nonlinear_lane_is_flat`]), not by this comment,
    /// so a new nonlinear feature that varies with something set here fails
    /// loudly instead of quietly weakening the test.
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

    /// A peptide of exactly `len` residues, so a fixture can put label signal in
    /// the `sequence_counts` block — the tail of the NONLINEAR lane.
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

    /// The fixture [`synthetic_competed`] is NOT: label signal in the NONLINEAR
    /// lane, carried by peptide length (and therefore by the residue counts) so
    /// that a GBM handed the nonlinear lane alone has something to split on.
    ///
    /// Needed because `synthetic_competed`'s nonlinear lane is entirely CONSTANT
    /// — measured, not assumed: `delta_group` and `delta_group_ratio`, the fields
    /// it varies for "GBM signal", are all-lane columns 90 and 91, i.e. inside
    /// `LINEAR_NCOLS`. So no fixture in this file could exercise "the GBM keeps
    /// working off the rest of the frame" before this one, at any row count.
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

    /// EXACTLY what a hybrid's failure arm hands the GBM: the nonlinear lane plus
    /// an all-zero score column, through the real [`hybrid_frame`] and the real
    /// [`CrossValidatedScorer`]. The only thing not reproduced is the failing
    /// cross-fit that would have produced those zeros.
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

    /// Not one nonlinear-lane column is constant across `fixture`'s rows — the
    /// mirror of [`assert_nonlinear_lane_is_flat`], and the premise of the
    /// benign-degradation half below.
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

    /// THE PREMISE of the hybrids' failure policy: when a cross-fit fails, the
    /// appended column goes uniformly zero and the GBM carries on off the rest of
    /// the frame — so the run degrades to a weaker model rather than to a
    /// meaningless one.
    ///
    /// That claim was shipped in `abort_standalone_mlp`'s doc as the whole reason
    /// the hybrids may degrade where the standalone rescorers abort, and it was
    /// UNTESTED. This test pins both halves of it, because the claim is
    /// CONDITIONAL and the condition is not decorative:
    ///
    /// * (a) with a trainable nonlinear lane, the degraded frame builds trees and
    ///   ranks targets above decoys — the premise holds;
    /// * (b) with `synthetic_competed`, whose nonlinear lane is entirely constant,
    ///   the SAME degraded frame reports no importance on any fold and emits one
    ///   score per fold. The q-values are then pure shuffle order — precisely the
    ///   outcome the standalone abort exists to prevent.
    ///
    /// (b) is not a bug being documented as a feature: it is the boundary of the
    /// premise, and `rescore_hybrid` / `rescore_hybrid_mlp_with` now detect that
    /// exact state at runtime (see [`untrained_folds`]) and say so
    /// in the log instead of returning it quietly. Do not delete (b) — without it
    /// the next reader will believe the degradation is unconditionally benign,
    /// which is the error this test was written to correct.
    #[test]
    fn degraded_hybrid_ranks_only_when_the_rest_of_the_frame_is_trainable() {
        // (a) THE PREMISE, on a fixture whose nonlinear lane can be trained.
        let trainable = synthetic_competed_nonlinear_signal(360);
        assert_nonlinear_lane_varies(&trainable);
        let (out, stats) = degraded_hybrid(trainable);
        assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);
        for fs in &stats {
            assert!(
                !fs.feature_importance.is_empty(),
                "fold {}: the degraded frame built no trees even though the nonlinear lane \
                 is trainable, so the hybrids' whole justification for degrading is false",
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
             a meaningless one, and the hybrids should abort like the standalone path does"
        );
        assert_eq!(
            untrained_folds(&stats),
            0,
            "the runtime detector must not fire when every fold of the degraded frame trained"
        );

        // (b) THE CONDITION, on the suite's main fixture: nonlinear lane constant,
        // so the same degradation yields one constant score per fold.
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

    /// THE WORST STATE, and the one a `all`-based detector was blind to: SOME
    /// folds of a degraded frame build trees and some build none.
    ///
    /// Two thirds of the rows then carry a per-fold constant while one third
    /// carries real discriminant values, so a row's score distribution is a
    /// function of its POSITION IN THE SEEDED SHUFFLE. [`crossfit`]'s failure
    /// policy names that exact shape as the one that "silently corrupts the
    /// q-value ranking … in a way that no downstream check would catch", and it is
    /// worse than the uniform case precisely because it is not uniform: uniform
    /// scores are visibly degenerate, a fold-dependent mixture looks like a
    /// working result.
    ///
    /// The row count is chosen from a measured plateau, not guessed:
    /// `synthetic_competed_nonlinear_signal` through the degraded frame gives
    /// `[0, 0, 0]` at 192-198, `[1, 0, 0]` at 200-206, `[1, 0, 1]` at 208 and
    /// `[1, 2, 1]` at 210. 204 sits in the middle of the partial plateau. The
    /// partial-ness is ASSERTED rather than assumed, so a shift in that boundary
    /// fails here instead of quietly turning this into a second copy of the
    /// uniform test.
    #[test]
    fn detector_fires_when_only_some_folds_of_a_degraded_frame_trained() {
        let (_out, stats) = degraded_hybrid(synthetic_competed_nonlinear_signal(204));
        assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);

        let untrained = untrained_folds(&stats);
        assert!(
            untrained > 0 && untrained < stats.len(),
            "fixture premise broken: this must be the PARTIAL state (some folds trained, some \
             did not), or it does not distinguish `any` from `all`. importance lengths: {:?}",
            stats
                .iter()
                .map(|fs| fs.feature_importance.len())
                .collect::<Vec<_>>()
        );

        // THE ASSERTION, on the REPORTER and not just on the count: a threshold of
        // `== stats.len()` would keep the worst outcome on the degraded path off
        // the operator's log entirely.
        let capture = |degraded: Option<&str>| -> String {
            let log = ErrorLog::default();
            tracing::subscriber::with_default(log.clone(), || {
                report_untrained_folds("nonlinear + lda_score", degraded, &stats)
            });
            let msgs = log.messages();
            msgs.iter()
                .find(|m| m.contains("built no trees"))
                .cloned()
                .unwrap_or_else(|| {
                    panic!(
                        "a PARTIALLY untrained frame was not reported at all — the worst state \
                         on this path reached an operator in silence. Captured: {msgs:#?}"
                    )
                })
        };

        // Both diagnoses must reach the log, and both must carry the count and
        // the partial-ness — those are what separate the fold-dependent case
        // from the uniform one.
        let with_degradation = capture(Some("lda_score"));
        let untrained_only = capture(None);
        for hit in [&with_degradation, &untrained_only] {
            assert!(
                hit.contains(&format!("{untrained} of {} GBM folds", stats.len())),
                "the message must say HOW MANY folds failed, or the operator cannot tell the \
                 fold-dependent case from the uniform one: {hit}"
            );
            assert!(
                hit.contains("fold-dependent"),
                "a partial failure must be named as such — it is worse than the uniform case, \
                 not better: {hit}"
            );
        }

        // …and they must stay DISTINGUISHABLE: "the fallback column degraded AND
        // nothing trained" is a different diagnosis from "nothing trained at
        // all", and an operator acts on them differently.
        assert!(
            with_degradation.contains("lda_score")
                && with_degradation.contains("zeroed by a failed cross-fit"),
            "the degraded arm must name the zeroed column and say it was zeroed: \
             {with_degradation}"
        );
        assert!(
            !untrained_only.contains("zeroed by a failed cross-fit"),
            "with no cross-fit failure the message must not invent one — that sends the reader \
             looking at the wrong half of the pipeline: {untrained_only}"
        );
        assert!(
            untrained_only.contains("no failed cross-fit"),
            "the untrained-only arm has to say the cross-fit was fine, or it reads as the \
             degraded one with a word missing: {untrained_only}"
        );
    }

    /// A `tracing::Subscriber` that records ERROR-level messages, so a test can
    /// assert an operator would actually SEE something rather than that a
    /// predicate returned true.
    ///
    /// Hand-rolled: `tracing-subscriber` is not a dependency of this crate and
    /// this work may not add one.
    #[derive(Clone, Default)]
    struct ErrorLog(Arc<std::sync::Mutex<Vec<String>>>);

    impl ErrorLog {
        fn messages(&self) -> Vec<String> {
            self.0.lock().unwrap().clone()
        }
    }

    impl tracing::Subscriber for ErrorLog {
        fn enabled(&self, md: &tracing::Metadata<'_>) -> bool {
            *md.level() == tracing::Level::ERROR
        }

        fn new_span(&self, _: &tracing::span::Attributes<'_>) -> tracing::span::Id {
            tracing::span::Id::from_u64(1)
        }

        fn record(&self, _: &tracing::span::Id, _: &tracing::span::Record<'_>) {}

        fn record_follows_from(&self, _: &tracing::span::Id, _: &tracing::span::Id) {}

        fn event(&self, event: &tracing::Event<'_>) {
            if *event.metadata().level() != tracing::Level::ERROR {
                return;
            }
            struct Msg(String);
            impl tracing::field::Visit for Msg {
                fn record_debug(
                    &mut self,
                    field: &tracing::field::Field,
                    value: &dyn std::fmt::Debug,
                ) {
                    if field.name() == "message" {
                        self.0 = format!("{value:?}");
                    }
                }
            }
            let mut msg = Msg(String::new());
            event.record(&mut msg);
            self.0.lock().unwrap().push(msg.0);
        }

        fn enter(&self, _: &tracing::span::Id) {}

        fn exit(&self, _: &tracing::span::Id) {}
    }

    /// THE WIRING, for BOTH hybrids: a degraded frame that trained nothing must
    /// reach the LOG, not just a predicate.
    ///
    /// [`untrained_folds`] being right is worth nothing if a rescorer never calls
    /// it, and the two arms compute their `degraded` flag independently — a
    /// hardcoded `false` at either call site would otherwise be caught by nothing.
    /// So this drives the real public entry points and asserts on captured ERROR
    /// output, once per arm, keyed on the COLUMN NAME so one arm's message cannot
    /// stand in for the other's.
    ///
    /// `indistinguishable_competed(24)` degrades BOTH cross-fits: the LDA's
    /// scatter matrix is singular even after shrinkage when every column is
    /// constant, and the MLP's cull leaves nothing (see the fixture's doc). The GBM
    /// then trains on a nonlinear lane that is constant too, so no fold builds a
    /// tree and the detector must fire in both runs.
    #[test]
    fn both_hybrids_log_a_degraded_frame_that_trained_nothing() {
        for (column, run) in [
            (
                "lda_score",
                (|| {
                    rescore_hybrid(indistinguishable_competed(24));
                }) as fn(),
            ),
            ("mlp_score", || {
                rescore_hybrid_mlp_with(indistinguishable_competed(24), mlp_test_cfg(7));
            }),
        ] {
            let log = ErrorLog::default();
            tracing::subscriber::with_default(log.clone(), run);

            let msgs = log.messages();
            let hit = msgs.iter().find(|m| {
                m.contains(column) && m.contains("built no trees") && m.contains("NOT usable")
            });
            assert!(
                hit.is_some(),
                "{column}: the degraded-and-untrained state never reached the log, so the \
                 hybrid either did not flag its own degradation or never ran the check. \
                 Captured ERROR messages: {msgs:#?}"
            );
            // The count is in the message, so the operator learns whether it was
            // the uniform case or the (worse) fold-dependent one.
            let hit = hit.unwrap();
            assert!(
                hit.contains(&format!(
                    "{} of {} GBM folds",
                    N_RESCORE_FOLDS, N_RESCORE_FOLDS
                )),
                "{column}: every fold was untrained here, so the message must say so: {hit}"
            );
        }
    }

    /// THE WIRING for the DEFAULT rescorer, which is the path most users run and
    /// the one the detector was originally not pointed at.
    ///
    /// [`rescore`] has no cross-fit column, so under the old `degraded`-gated
    /// check it could never report anything: a search below
    /// `GBMConfig::default`'s `min_leaf_weight` floor emitted per-fold constants
    /// and q-values that rank by shuffle order, silently. That is pre-existing
    /// behavior rather than a regression, which is exactly why it needed a test —
    /// nothing else on this branch would have noticed the omission.
    ///
    /// `indistinguishable_competed(24)` makes every ALL-lane column bit-constant,
    /// so no fold builds a tree (see the fixture's doc). The message must ALSO not
    /// blame a cross-fit failure, because there was no cross-fit: the two
    /// diagnoses are different and an operator acts on them differently.
    #[test]
    fn plain_rescore_reports_a_frame_that_trained_nothing() {
        let log = ErrorLog::default();
        tracing::subscriber::with_default(log.clone(), || {
            rescore(indistinguishable_competed(24));
        });

        let msgs = log.messages();
        let hit = msgs
            .iter()
            .find(|m| m.contains("built no trees") && m.contains("NOT usable"))
            .unwrap_or_else(|| {
                panic!(
                    "the default rescorer trained nothing and said nothing — the detector is \
                     not wired into the path most users run. Captured ERROR messages: {msgs:#?}"
                )
            });
        assert!(
            hit.contains(&format!(
                "{} of {} GBM folds",
                N_RESCORE_FOLDS, N_RESCORE_FOLDS
            )),
            "every fold was untrained here, so the message must say so: {hit}"
        );
        assert!(
            !hit.contains("zeroed by a failed cross-fit"),
            "`rescore` has no cross-fit column to zero, so blaming one sends the reader to a \
             stage this path does not have: {hit}"
        );
    }

    /// The cross-fit must actually HOLD OUT: a row's score has to come from a
    /// model that never saw that row.
    ///
    /// Construction (60 rows, alternating target/decoy, 3 columns) — every
    /// column is a label copy or label-balanced noise, arranged so the only
    /// thing that could separate FOLD 0 is a column that is constant on every
    /// other fold:
    ///
    /// * `col0` — a copy of the label on fold-0 rows (`i % N_RESCORE_FOLDS == 0`)
    ///   and exactly `0.0` everywhere else.
    /// * `col1` — the mirror: a copy of the label on non-fold-0 rows, `0.0` on
    ///   fold-0 rows. Keeps the training folds separable so the fit is not
    ///   degenerate.
    /// * `col2` — varies per row but is label-balanced inside EVERY fold (each
    ///   value appears once as a target and once as a decoy per fold), so it
    ///   supplies within-class scatter — needed for a non-singular `Sw` — with
    ///   no discriminative signal in any subset of folds.
    ///
    /// The model for fold 0 is fitted on folds 1..N only, where `col0` is
    /// constant -> zero weight; and fold-0 rows all have `col1 == 0`. So every
    /// held-out fold-0 score must collapse to the SAME value: zero separation.
    /// An in-sample fit over all rows does see `col0` track the label and
    /// separates those very rows cleanly. Both halves are asserted — the
    /// in-sample control is what stops the first assertion from passing
    /// vacuously (e.g. for a helper that returned all zeros).
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
            // col2: constant within each block of 6 consecutive rows. A block of
            // 6 holds exactly one target and one decoy of each fold, so this
            // column's class means are equal in every fold and in every union of
            // folds.
            feat.push(((i / (2 * n_folds)) % 5) as f64); // col2
            is_decoy.push(!is_target(i));
        }

        // Same matrix, handed over as the dataset the cross-fit now reads
        // through. `RowMajorDataset`'s `get_fold` IS `i % N_RESCORE_FOLDS`, so
        // the fold-0 construction above still describes the partition exactly.
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

        // (a) HELD OUT: fold-0 scores carry no information -> they are all equal.
        let held: Vec<f64> = fold0.iter().map(|&i| cf.scores[i]).collect();
        let spread = held.iter().cloned().fold(f64::MIN, f64::max)
            - held.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            spread < 1e-6,
            "held-out fold-0 scores must not separate; spread={spread} scores={held:?}"
        );

        // (b) CONTROL: the same rows, scored by an in-sample fit over ALL rows,
        // separate perfectly — so (a) is a property of the hold-out, not of the
        // data being unlearnable.
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

    /// 90 rows, like `rescore_mlp_is_deterministic_on_both_lanes` and for the same
    /// reason: `min_leaf_weight`'s no-trees trap is a property of `GBMConfig` and
    /// this path trains no GBM. An LDA reports a finite `|coef|` for every linear
    /// -lane column at any row count, so the non-empty assertions below cannot go
    /// vacuous the way a GBM's sparse importance can. The row count is free here,
    /// so it is small. (It sat at 360 for a while with no comment saying why, which
    /// read as if the GBM trap applied.)
    #[test]
    fn rescore_lda_is_cross_fit_and_deterministic() {
        let n = 90;
        let (out_a, stats_a) = rescore_lda(synthetic_competed(n));
        let (out_b, _stats_b) = rescore_lda(synthetic_competed(n));

        assert_eq!(out_a.len(), n as usize);

        // Per-fold stats, one FoldStats per fold (the GBM path's shape), which
        // also proves every fold's LDA fit succeeded — the single-fit path used
        // to emit exactly one `fold: 0` entry.
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

    /// The LDA sidecar is FULL WIDTH: one importance row per linear-lane
    /// feature, per fold, including the features whose `|coef|` is exactly 0.0.
    ///
    /// An LDA looks at every column, so it has a measurement for every column —
    /// a zero weight means "this feature is dead or constant in this fold",
    /// which is a finding, not a gap. The sidecar used to lose those rows to a
    /// `!= 0.0` filter meant for the tree model's "never split on this"
    /// sentinel; under `synthetic_competed` most linear fields are constant, so
    /// that filter dropped the large majority of the lane here.
    ///
    /// Both directions are pinned: the width, and that the zeros are really
    /// present (a full-width vector of all-nonzero values would satisfy the
    /// length check alone).
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

        // Names are the linear lane exactly — no duplicates, no strays.
        let lane: std::collections::HashSet<Arc<str>> =
            linear_feature_name_set().into_iter().collect();
        let reported: std::collections::HashSet<Arc<str>> = stats[0]
            .feature_importance
            .iter()
            .map(|(n, _)| n.clone())
            .collect();
        assert_eq!(reported, lane);
    }

    // -----------------------------------------------------------------
    // MLP rescorer
    // -----------------------------------------------------------------

    /// Test config for the MLP rescorers: same architecture family as
    /// [`MlpConfig::default`], shrunk so the seed sweeps below stay fast, and
    /// with the batch size cut so the fixtures here (90-360 rows) actually get
    /// several minibatch steps per epoch instead of one full-batch step. The
    /// default is left alone on purpose — tuning a production default for test
    /// runtime is how a default stops meaning anything.
    ///
    /// # Why `lr` is `1e-3` and not the default's `3e-4`
    /// [`MlpConfig::default`] records that `1e-2` and `1e-3` both "walked a
    /// 2-input ReLU net into a dead-unit plateau in testing", so raising it here
    /// looks like contradicting a warning written 400 lines away. It is not: that
    /// warning is about the DEFAULT's step budget (30 epochs at batch 256, i.e.
    /// roughly one full-batch step per epoch), where a large step has nothing to
    /// recover from a dead unit with. This config runs 150 epochs at batch 16 —
    /// one to two orders of magnitude more steps on the same fixtures — and the
    /// plateau does not reproduce across the four-seed sweeps below.
    ///
    /// `3e-4` was tried and MEASURABLY UNDER-TRAINS at this budget:
    /// `mlp_crossfit_scores_are_held_out`'s leak control has to reach AUC > 0.99
    /// on a 30-row fold and only reaches 0.92, which would make the hold-out
    /// assertion it guards vacuous. So `1e-3` is a measured requirement of these
    /// fixtures, not an oversight — do not "align it with the default" without
    /// re-running that test.
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

    /// Fraction of (target, decoy) pairs ranked correctly, ties counted as half
    /// a win. Threshold-free and invariant to the logit's offset and scale,
    /// which is what makes it comparable between a held-out fit and an
    /// in-sample one.
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

    /// A bare [`FeatureLike`] row: a label and a score slot, which is all
    /// [`CrossValidatedScorer`] reads from `T` once the matrix is supplied
    /// externally. Lets the leak test below drive the REAL scorer over a
    /// hand-built matrix instead of re-deriving its partition.
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

    /// The MLP cross-fit must actually HOLD OUT: a row's score has to come from
    /// a model that never saw that row.
    ///
    /// This re-points the `crossfit_lda_scores_are_held_out` construction at
    /// [`CrossValidatedScorer`] + [`MlpFoldModel`], i.e. the exact pair
    /// [`rescore_mlp_lane`] uses. The re-pointing is not cosmetic: the two
    /// partitions differ, and the fixture has to be built for the one under
    /// test. Under `CrossValidatedScorer` with `N_RESCORE_FOLDS == 3`, model `f`
    /// trains on fold `f`, early-stops on `f + 1` and scores fold `f + 2`, so
    /// FOLD 0 IS SCORED BY MODEL 1 ALONE — trained on fold-1 rows only. (Under
    /// `crossfit_lda` it would instead be scored by a model trained on folds 1
    /// AND 2.)
    ///
    /// Construction (90 rows, alternating target/decoy, 3 columns):
    ///
    /// * `col0` — a copy of the label on fold-0 rows, exactly `0.0` elsewhere.
    ///   Constant on fold-1 rows, so `ColumnTransform` CULLS it from model 1
    ///   outright: the leaking column is not even an input to the model that
    ///   scores fold 0.
    /// * `col1` — the mirror: the label on non-fold-0 rows, `0.0` on fold-0
    ///   rows. Keeps model 1 learnable (so the fit is not degenerate) while
    ///   being a constant input on every row it scores.
    /// * `col2` — constant within each block of `2 * n_folds` consecutive rows.
    ///   A block holds exactly one target and one decoy of each fold, so within
    ///   fold 0 the targets and the decoys carry the IDENTICAL multiset of
    ///   values. It supplies genuine variance without discriminative signal.
    ///
    /// So the only input that varies across the fold-0 rows model 1 scores is
    /// `col2`, whose values pair up exactly between the classes: the held-out
    /// AUC must be exactly 0.5.
    ///
    /// NON-VACUITY is the in-sample control, and it is doing real work here:
    /// fitted over ALL rows, `col0` is no longer constant, survives the cull,
    /// and tracks the label on precisely the rows being scored — so the same
    /// model class on the same data separates them perfectly. Without it, a
    /// `predict` that returned a constant, a scorer that never fitted, or a
    /// fixture whose data was simply unlearnable would all sail through the
    /// first assertion.
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

        // Several seeds: the MLP's initialization is the one thing that varies
        // run to run, and a single seed reports an init-dependent training trap
        // as a clean pass or a leak purely on luck.
        for seed in [7u64, 13, 42, 1234] {
            // (a) HELD OUT, through the real scorer and the real partition.
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

            // (b) CONTROL: same model class, same data, fitted IN SAMPLE over
            // every row, scoring those same fold-0 rows.
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

            // (c) THE LEAK, made concrete: this is the model a partition bug
            // would hand fold 0 — trained on fold-0 rows, where `col0` is the
            // label. It separates them perfectly. So (a) is measuring the
            // partition, not a fixture that happens to be unlearnable.
            let leaky = MlpFoldModel::fit(&mlp_test_cfg(seed), &dataset, 0, &fold0, &[]).unwrap();
            let leaky_auc = pair_auc(&leaky.predict(&dataset, &fold0).unwrap(), &fold0_is_target);
            assert!(
                leaky_auc > 0.99,
                "seed {seed}: a model trained ON fold 0 must separate fold 0 \
                 (AUC={leaky_auc}), else (a) would hold even for a leaking scorer"
            );
        }
    }

    /// Determinism + the per-fold sidecar shape for the MLP rescorers, over
    /// BOTH lanes and several seeds. Bit-identical is the assertion, as in
    /// `rescore_lda_is_cross_fit_and_deterministic`: the whole q-value pipeline
    /// is a ranking of these scores, so "close" is not a property worth having.
    ///
    /// 90 rows, not the 360 the GBM-bearing tests need. `min_leaf_weight`'s
    /// no-trees trap is a property of `GBMConfig`, and this path trains no GBM:
    /// `MlpFoldModel` reports a finite importance for every lane column at any
    /// row count, so the full-width assertions below cannot go vacuous the way a
    /// GBM's sparse importance can. The row count is free here, so it is small.
    #[test]
    fn rescore_mlp_is_deterministic_on_both_lanes() {
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

        for lane in [Lane::Linear, Lane::All] {
            for seed in [7u64, 13, 42, 1234] {
                let run = || rescore_mlp_lane(synthetic_competed(n), lane, mlp_test_cfg(seed));
                let (out_a, stats_a) = run();
                let (out_b, _) = run();

                assert_eq!(out_a.len(), n as usize);
                for r in &out_a {
                    assert!(
                        r.discriminant_score.is_finite(),
                        "{lane:?}/{seed}: non-finite score {}",
                        r.discriminant_score
                    );
                    assert!((0.0..=1.0).contains(&r.qvalue));
                }
                assert_eq!(
                    key(&out_a),
                    key(&out_b),
                    "{lane:?}/{seed}: mlp rescore not deterministic"
                );

                // One FoldStats per fold, FULL WIDTH on both halves. Stronger
                // than the LDA path's non-empty checks, and it can be: the MLP
                // measures every lane column (0.0 for the culled ones), so
                // anything short of the lane width is a dropped row.
                assert_eq!(stats_a.len(), N_RESCORE_FOLDS as usize);
                for (f, fs) in stats_a.iter().enumerate() {
                    assert_eq!(fs.fold, f as u8);
                    assert_eq!(
                        fs.feature_stats.len(),
                        lane.ncols(),
                        "{lane:?}/{seed}: fold {f} stats must be lane width"
                    );
                    assert_eq!(
                        fs.feature_importance.len(),
                        lane.ncols(),
                        "{lane:?}/{seed}: fold {f} — the MLP measures every lane \
                         column (0.0 for culled ones), so none may be dropped"
                    );
                    assert!(
                        fs.feature_importance.iter().all(|(_, g)| g.is_finite()),
                        "{lane:?}/{seed}: fold {f} must never emit the NAN sentinel"
                    );
                    assert!(
                        fs.feature_importance.iter().any(|(_, g)| *g > 0.0),
                        "{lane:?}/{seed}: fold {f} reported no weight at all"
                    );
                }
            }
        }
    }

    /// ROW ALIGNMENT: the lane matrix must be built AFTER
    /// [`canonicalize_and_shuffle`], so matrix row `i` is `data[i]`.
    ///
    /// Nothing in the type system says so — both orders compile, both produce a
    /// correctly-shaped matrix, and both run to completion. What a pre-shuffle
    /// build changes is which LABEL and which FOLD each feature row is paired
    /// with: the shuffle is a permutation, so every row would carry another
    /// candidate's label. The model then has nothing to learn and the scores
    /// land at chance.
    ///
    /// `synthetic_competed` gives targets and decoys clearly different counts,
    /// so a correctly-aligned fit ranks them apart; the assertion is that it
    /// does, well above the 0.5 a misaligned build is pinned to. The upper
    /// bound is not asserted — how WELL it separates is the model's business,
    /// not this test's.
    #[test]
    fn mlp_lane_matrix_is_row_aligned_with_the_shuffled_data() {
        for seed in [7u64, 13, 42, 1234] {
            for lane in [Lane::Linear, Lane::All] {
                let (out, _) = rescore_mlp_lane(synthetic_competed(120), lane, mlp_test_cfg(seed));
                let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
                let is_target: Vec<bool> =
                    out.iter().map(|r| r.scoring.identity.is_target).collect();
                let auc = pair_auc(&scores, &is_target);
                assert!(
                    auc > 0.9,
                    "{lane:?}/{seed}: AUC {auc} is near chance. The cause this test was \
                     written for is a lane matrix built BEFORE the shuffle rather than \
                     after it, which pairs every feature row with another candidate's \
                     label and fold — but any MLP training regression lands here too. \
                     Check `mlp_fold`'s `fitted_model_separates_a_separable_toy_on_held_\
                     out_rows` first: if that still passes, the fit is fine and the \
                     ordering is the problem"
                );
            }
        }
    }

    /// Each public entry point runs on the default [`MlpConfig`] (not the shrunk
    /// test one), stays deterministic, and — the load-bearing part — TRAINS ON
    /// THE LANE ITS NAME CLAIMS.
    ///
    /// The lane assertion is what makes the naming a fact rather than a
    /// convention. Nothing type-checks it: both entry points have the same
    /// signature, both run to completion on either lane, and both produce
    /// plausible q-values either way, so an inverted wiring is invisible
    /// everywhere except in which features the model actually saw. That is
    /// recoverable from the output: the sidecar's importance rows are the lane's
    /// name set exactly (`MlpFoldModel` measures every column), so comparing
    /// them against `Lane::names()` pins the lane from the outside.
    ///
    /// A length check alone would be weaker but would still work here, since
    /// `LINEAR_NCOLS != ALL_NCOLS`; the name-set comparison is used instead
    /// because it does not depend on that inequality holding as features are
    /// added and removed.
    ///
    /// This test is also the one that would have caught the inversion this
    /// function's names were shipped with (`rescore_mlp` = ALL,
    /// `rescore_mlp_linear` = LINEAR, against a spec that maps
    /// `RescoreModel::Mlp` to the LINEAR lane): the old version iterated both
    /// entry points but asserted nothing about which lane either one used.
    #[test]
    fn public_mlp_entry_points_train_on_the_lane_they_name() {
        type Rescorer = fn(Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats);

        for (name, rescorer, lane) in [
            (
                "rescore_mlp_linear",
                rescore_mlp_linear as Rescorer,
                Lane::Linear,
            ),
            ("rescore_mlp_all", rescore_mlp_all as Rescorer, Lane::All),
        ] {
            let (out_a, stats) = rescorer(synthetic_competed(90));
            let (out_b, _) = rescorer(synthetic_competed(90));
            assert_eq!(out_a.len(), 90);
            assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);

            // THE LANE, read back out of the sidecar.
            let expected: std::collections::HashSet<Arc<str>> = lane.names().into_iter().collect();
            for fs in &stats {
                let reported: std::collections::HashSet<Arc<str>> = fs
                    .feature_importance
                    .iter()
                    .map(|(n, _)| n.clone())
                    .collect();
                assert_eq!(
                    reported,
                    expected,
                    "{name} fold {}: trained on the wrong lane — it reported {} features, \
                     but {lane:?} has {}. The entry point's name and the lane it passes to \
                     `rescore_mlp_lane` have been inverted.",
                    fs.fold,
                    reported.len(),
                    expected.len(),
                );
            }

            let bits = |out: &[FinalResult]| -> Vec<(u32, u32)> {
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
            assert_eq!(bits(&out_a), bits(&out_b), "{name} is not deterministic");
            assert!(out_a.iter().all(|r| r.discriminant_score.is_finite()));
        }
    }

    /// The enum -> lane mapping, pinned through the dispatcher users actually
    /// reach (`rescore_model` config field / `--rescore-model`).
    ///
    /// `public_mlp_entry_points_train_on_the_lane_they_name` pins each entry
    /// point against the lane in its own name; it says nothing about which entry
    /// point `rescore_with` hands a variant to. That arm is a separate,
    /// untype-checked mapping: `RescoreModel::Mlp => rescore_mlp_all(data)`
    /// compiles, runs, and yields plausible q-values off the wrong feature set.
    /// So the two halves are pinned separately, and this is the half that pins
    /// the dispatch.
    ///
    /// Technique for the two pure-MLP variants is the sibling test's: the
    /// sidecar's importance rows ARE the trained lane's name set, because
    /// `MlpFoldModel` reports a finite importance for every lane column
    /// (culled ones included). `HybridMlp` trains a GBM, whose importance is
    /// sparse at this fixture size (no tree ever splits), so it is pinned on
    /// `feature_stats` instead — those names come from the dataset's columns,
    /// not from the model, and therefore identify the frame the GBM was handed
    /// (`nonlinear ++ mlp_score`) regardless of what it learned.
    ///
    /// A name set pins the LANE, not the MODEL: `RescoreModel::Mlp =>
    /// rescore_lda(data)` would satisfy the checks above, since the LDA trains on
    /// the linear lane too and also reports dense full-width importance. So each
    /// arm is ALSO compared BIT-FOR-BIT against the entry point it names. Every
    /// path here is deterministic (pinned by
    /// `public_mlp_entry_points_train_on_the_lane_they_name` and
    /// `rescore_hybrid_mlp_smoke_and_determinism`), so bit-equality is a fair
    /// test, and it is what makes the mapping a statement about the model rather
    /// than about the feature set.
    #[test]
    fn rescore_with_dispatches_each_mlp_variant_to_its_lane() {
        use crate::ml::{
            RescoreModel,
            rescore_with,
        };

        let expected_set = |names: Vec<Arc<str>>| -> std::collections::HashSet<Arc<str>> {
            names.into_iter().collect()
        };
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

        // --- The two pure-MLP variants, read off the MLP's own importance ---
        for (variant, lane) in [
            (RescoreModel::Mlp, Lane::Linear),
            (RescoreModel::MlpAll, Lane::All),
        ] {
            let (out, stats) = rescore_with(variant, synthetic_competed(90));
            assert_eq!(out.len(), 90);
            assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);

            let expected = expected_set(lane.names());
            for fs in &stats {
                let reported: std::collections::HashSet<Arc<str>> = fs
                    .feature_importance
                    .iter()
                    .map(|(n, _)| n.clone())
                    .collect();
                assert_eq!(
                    reported,
                    expected,
                    "rescore_with({variant:?}) fold {}: dispatched to the wrong lane — it \
                     reported {} features, but {lane:?} has {}. The `rescore_with` arm for \
                     {variant:?} calls the other lane's entry point.",
                    fs.fold,
                    reported.len(),
                    expected.len(),
                );
            }

            // ...and it is the MLP entry point, not some other model on the same
            // lane. See the doc.
            let direct = match variant {
                RescoreModel::Mlp => rescore_mlp_linear(synthetic_competed(90)).0,
                _ => rescore_mlp_all(synthetic_competed(90)).0,
            };
            assert_eq!(
                key(&out),
                key(&direct),
                "rescore_with({variant:?}) does not produce what its own entry point does, \
                 so the arm calls a DIFFERENT MODEL that happens to train on {lane:?}"
            );
        }

        // --- HybridMlp: the GBM's frame is `nonlinear ++ mlp_score` ---
        let (out, stats) = rescore_with(RescoreModel::HybridMlp, synthetic_competed(90));
        assert_eq!(out.len(), 90);
        assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);

        let mut hybrid_names = nonlinear_feature_name_set();
        hybrid_names.push(Arc::from("mlp_score"));
        let expected = expected_set(hybrid_names);
        for fs in &stats {
            let reported: std::collections::HashSet<Arc<str>> =
                fs.feature_stats.iter().map(|f| f.name.clone()).collect();
            // The LDA hybrid's frame is the SAME WIDTH (nonlinear + one score
            // column), so a count comparison would not separate the two — the
            // discriminating column is the score's NAME.
            assert!(
                reported.contains("mlp_score"),
                "rescore_with(HybridMlp) fold {}: the GBM's frame has no `mlp_score` column. \
                 The arm does not call `rescore_hybrid_mlp` (extra columns: {:?}).",
                fs.fold,
                reported.difference(&expected).collect::<Vec<_>>(),
            );
            assert_eq!(
                reported, expected,
                "rescore_with(HybridMlp) fold {}: the GBM was handed something other than the \
                 nonlinear lane plus `mlp_score`.",
                fs.fold,
            );
        }
    }

    /// THE DEFAULT PATH, plus the `Gbm` / `Lda` / `Hybrid` arms of
    /// [`rescore_with`].
    ///
    /// [`rescore`] is what runs when nobody sets `rescore_model`, and nothing
    /// called it before this test. Its three arms had nothing pinning their
    /// mapping either: the same untype-checked dispatch that
    /// `rescore_with_dispatches_each_mlp_variant_to_its_lane` covers for the MLP
    /// variants, and it matters more here, because `Gbm` is the one that runs by
    /// default.
    ///
    /// Pinned on BIT-IDENTICAL scores against the entry point each arm is
    /// supposed to call. That is stronger than a sidecar-name comparison, which
    /// pins only the FEATURE SET: `Lda => rescore_mlp_linear` reports the same
    /// dense linear-lane name set the LDA does and would satisfy one. Every path
    /// here is deterministic run to run — pinned separately by
    /// `rescore_lda_is_cross_fit_and_deterministic` and
    /// `rescore_hybrid_smoke_and_determinism` — so bit-equality is a fair test.
    ///
    /// The frame WIDTH is checked too, off `feature_stats` (whose names come from
    /// the dataset's columns rather than from the model, so they hold with or
    /// without trees). The three widths are distinct — `ALL_NCOLS`,
    /// `LINEAR_NCOLS`, `NONLINEAR_NCOLS + 1` — so this alone separates the three
    /// arms' lanes.
    ///
    /// NON-VACUITY is the last block: the three arms must produce three DIFFERENT
    /// score vectors, or "each arm calls its own entry point" would also hold for
    /// a permutation of them. `n = 360` is what makes that true.
    /// `GBMConfig::default`'s `min_leaf_weight` of 5.0 exceeds a 30-row fold's
    /// total logloss hessian, so at 90 rows both GBM-bearing arms emit one
    /// constant per fold and would be far harder to tell apart. Do not shrink it.
    #[test]
    fn rescore_with_dispatches_the_default_and_the_non_mlp_variants() {
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

        // The variant the config file advertises as the default really is the
        // one `#[default]` selects — the flag being absent has to land on `Gbm`.
        assert_eq!(
            RescoreModel::default(),
            RescoreModel::Gbm,
            "the default rescorer moved; `default_config.toml` and this test's whole \
             premise both say it is the GBM"
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

            // (a) THE FRAME the model was handed.
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

            // (b) THE MODEL: bit-identical to the entry point this arm names.
            let (direct_out, _) = direct(synthetic_competed(N));
            assert_eq!(
                key(&dispatched),
                key(&direct_out),
                "rescore_with({variant:?}) does not reproduce its own entry point, so the \
                 arm dispatches to a different model"
            );

            seen.push((variant, key(&dispatched)));
        }

        // (c) NON-VACUITY: three distinct models, three distinct score vectors.
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

    /// Smoke + determinism for the LDA hybrid.
    ///
    /// # Why the importance assertion is here
    /// Every OTHER assertion in this test holds in the regime where the GBM
    /// builds no trees at all: a per-fold constant score is finite, its q-values
    /// are in `[0, 1]`, and it is bit-identical run to run. `GBMConfig::default`'s
    /// `min_leaf_weight` of 5.0 exceeds a 30-row fold's total logloss hessian, so
    /// at `n = 90` this test would pass while measuring nothing but the
    /// reproducibility of a constant — measured: `n = 90` gives per-fold
    /// importance lengths `[0, 0, 0]`, `n = 360` gives `[1, 1, 1]`.
    ///
    /// `n = 360` is therefore load-bearing, and the non-empty importance check is
    /// what makes it so: an empty importance means no tree was ever built, i.e.
    /// the output does not depend on the features, the appended `lda_score`
    /// column, or anything else this test names. Do not drop it, and do not
    /// shrink `n` without watching it fail.
    ///
    /// It does NOT claim the output depends on `lda_score` specifically — that is
    /// `hybrid_lda_score_carries_the_linear_lane_into_the_gbm`'s job, on a fixture
    /// built so that nothing else could carry the signal.
    #[test]
    fn rescore_hybrid_smoke_and_determinism() {
        let n = 360;
        let (out_a, stats_a) = rescore_hybrid(synthetic_competed(n));
        let (out_b, _stats_b) = rescore_hybrid(synthetic_competed(n));

        // NON-VACUITY, first: if this fails, every assertion below is about a
        // per-fold constant rather than about a trained model.
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

        // (a) output length preserved
        assert_eq!(out_a.len(), n as usize);
        assert_eq!(out_b.len(), n as usize);

        // (b) all discriminant scores finite, (c) qvalues in [0, 1]
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

        // (d) determinism: same seed + sort key -> identical scores per candidate.
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

    // -----------------------------------------------------------------
    // MLP hybrid
    // -----------------------------------------------------------------

    /// A [`FoldModel`] that fits nothing and records the row slice it was handed.
    /// Lets a test read a partition's TRAIN set back out of the driver that built
    /// it, without a real fit's cost or its failure modes.
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

    /// THE assertion the hybrids' leak-freedom rests on: [`crossfit`] and
    /// [`CrossValidatedScorer`] must agree on the fold ASSIGNMENT.
    ///
    /// They deliberately DISAGREE on the train/score split — `crossfit` fits on
    /// all-but-fold-`f` and scores fold `f`, the scorer fits on fold `f` alone
    /// and scores the folds that are neither it nor its early-stopping
    /// neighbour. Both are leak-free on their own. What is not safe is the two
    /// slicing the rows into DIFFERENT folds: then a row's `mlp_score` can come
    /// from an MLP trained on rows the GBM is holding out, and because the
    /// offending value sits in a FEATURE rather than in the score, the GBM's own
    /// cross-validation cannot see it — the only symptom is an optimistically
    /// wrong FDR.
    ///
    /// So the check is a set equality between the two sides, not each side
    /// against a re-typed `i % N_RESCORE_FOLDS`: a repeated modulo is exactly
    /// the arrangement that let them drift before, and this test would pass
    /// against a drifted pair if it re-derived the answer itself.
    ///
    /// `N` is deliberately NOT a multiple of `N_RESCORE_FOLDS`, so a fold-size
    /// off-by-one has somewhere to show up.
    ///
    /// Both sides go through THE HYBRIDS' OWN HELPERS rather than a hand-built
    /// `RowMajorDataset`: [`hybrid_linear_dataset`] for the cross-fit's dataset
    /// and [`hybrid_frame`] for the scorer's. Each of those spells
    /// `N_RESCORE_FOLDS as usize` in its own body, and an earlier version of this
    /// test pinned the two DRIVERS while leaving those two literals uncompared —
    /// so a hybrid could have handed one side a differently-folded dataset and
    /// this test would still have passed.
    #[test]
    fn crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on() {
        const N: usize = 47;

        // The real pre-rescore sequence, so the row order under test is the one a
        // hybrid actually partitions.
        let mut data = synthetic_competed(N as u32);
        canonicalize_and_shuffle(&mut data);
        let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();

        // Side A: the cross-fit that produces `mlp_score` / `lda_score`, over the
        // dataset `rescore_hybrid` / `rescore_hybrid_mlp_with` build for it.
        let lin_dataset = hybrid_linear_dataset(&data, responses.clone());
        let cf = crossfit::<_, FoldSpy>(&lin_dataset, &(), "spy").expect("the spy cannot fail");

        // Side B: the GBM scorer the hybrid hands that column to, over the frame
        // `hybrid_frame` builds. No fit needed — the partition is built in the
        // constructor.
        let (precomputed, names) = hybrid_frame(&data, responses, "spy_score", cf.scores.clone());
        let scorer =
            CrossValidatedScorer::<CompetedCandidate, GbmFoldModel>::new_from_shuffled_with_precomputed(
                N_RESCORE_FOLDS,
                data,
                GBMConfig::default(),
                precomputed,
                names,
            );

        // (a) SAME ASSIGNMENT. `cf.fold_rows[f]` is the rows model `f` SCORED,
        // which under `crossfit` is exactly fold `f` — the same set the scorer
        // TRAINS its model `f` on.
        assert_eq!(
            cf.fold_rows,
            scorer.fold_rows(),
            "the cross-fit and the GBM scorer partition the rows differently, so a \
             hybrid row's cross-fit column can come from a model trained on rows \
             the GBM is holding out"
        );

        // (b) and the split really is the complement: model `f` trained on every
        // row it did NOT score, so no row is scored by a model that saw it.
        for (f, model) in cf.models.iter().enumerate() {
            let scored = &cf.fold_rows[f];
            let expected: Vec<usize> = (0..N).filter(|i| !scored.contains(i)).collect();
            assert_eq!(model.train, expected, "fold {f} trained on the wrong rows");
            assert!(
                !scored.is_empty() && !expected.is_empty(),
                "fold {f} is empty, so the equality above is vacuous"
            );
        }

        // (c) non-vacuity for (a): the partition is a real partition (every row
        // in exactly one fold), so (a) is not comparing two empty shapes.
        let mut seen: Vec<usize> = cf.fold_rows.iter().flatten().copied().collect();
        seen.sort_unstable();
        assert_eq!(seen, (0..N).collect::<Vec<_>>());
    }

    /// FAILURE POLICY: if ANY fold's MLP fit fails, `mlp_score` is uniformly
    /// zero — not the partial column the earlier folds already produced.
    ///
    /// The fixture makes the LAST fold the failing one, which is the only
    /// arrangement that can tell the two apart: both columns are the label-free
    /// value `0.0` on every fold-0 and fold-1 row and vary only on fold-2 rows,
    /// so
    ///
    /// * fold 0's model (trains on folds 1+2) sees varying columns — fits;
    /// * fold 1's model (trains on folds 0+2) likewise — fits;
    /// * fold 2's model (trains on folds 0+1) sees two constant columns, both
    ///   culled, and returns [`MlpFoldError::NoUsableColumns`].
    ///
    /// So a driver that returned its accumulator on failure would hand back real
    /// scores on two thirds of the rows and zeros on the rest. That column is
    /// worse than no column: it is a function of a row's position in the seeded
    /// shuffle, the GBM can split on it to recover fold identity, and the
    /// resulting q-value ranking is corrupt in a way nothing downstream checks.
    /// A uniform zero is constant, so it has zero split gain and the hybrid
    /// degrades to "GBM on the nonlinear lane".
    ///
    /// Asserted on RAW BITS, so a `-0.0` or a denormal cannot pass as zero.
    #[test]
    fn crossfit_score_column_is_uniformly_zero_when_any_fold_fails() {
        const N: usize = 60;
        const NCOLS: usize = 2;
        let n_folds = N_RESCORE_FOLDS as usize;
        let last = n_folds - 1;

        let mut feat = Vec::with_capacity(N * NCOLS);
        let mut responses = Vec::with_capacity(N);
        for i in 0..N {
            let live = i % n_folds == last;
            feat.push(if live { i as f64 } else { 0.0 });
            feat.push(if live { -(i as f64) } else { 0.0 });
            responses.push((i % 2) as f64);
        }
        let names: Vec<Arc<str>> = ["col0", "col1"].map(Arc::from).to_vec();
        let dataset = RowMajorDataset::new(
            PrecomputedFeatures::from_row_major(feat, NCOLS, responses),
            names,
            n_folds,
        );
        let complement = |f: usize| -> Vec<usize> { (0..N).filter(|i| i % n_folds != f).collect() };

        for seed in [7u64, 13, 42, 1234] {
            let cfg = mlp_test_cfg(seed);

            // The fixture does what its doc says: the last fold's model is the
            // one that cannot fit, and the earlier ones can.
            for f in 0..last {
                assert!(
                    MlpFoldModel::fit(&cfg, &dataset, f, &complement(f), &[]).is_ok(),
                    "seed {seed}: fold {f} must FIT, or the partial column this test \
                     rules out was never produced in the first place"
                );
            }
            assert_eq!(
                MlpFoldModel::fit(&cfg, &dataset, last, &complement(last), &[]).err(),
                Some(MlpFoldError::NoUsableColumns {
                    fold: last,
                    ncols: NCOLS,
                    train_rows: complement(last).len(),
                }),
                "seed {seed}: fold {last} must fail for this test to exercise the \
                 failure path at all"
            );

            // Non-vacuity: fold 0's model, which the driver DOES fit before
            // hitting the failure, produces scores that are not already zero —
            // so an implementation that leaked its accumulator would be visibly
            // different from the assertion below.
            let fold0: Vec<usize> = (0..N).step_by(n_folds).collect();
            let partial = MlpFoldModel::fit(&cfg, &dataset, 0, &complement(0), &[])
                .unwrap()
                .predict(&dataset, &fold0)
                .unwrap();
            assert!(
                partial.iter().any(|s| s.to_bits() != 0.0f64.to_bits()),
                "seed {seed}: the surviving folds score everything at exactly zero, \
                 so 'uniformly zero' and 'partially filled' are indistinguishable here"
            );

            // THE ASSERTION.
            let (col, degraded) =
                crossfit_score_column::<_, MlpFoldModel>(&dataset, &cfg, "MLP", "mlp_score");
            assert_eq!(
                col.len(),
                N,
                "seed {seed}: the column must still be full length"
            );
            assert!(
                col.iter().all(|s| s.to_bits() == 0.0f64.to_bits()),
                "seed {seed}: mlp_score must be uniformly zero after a failed fold, got {col:?}"
            );
            // The caller has to be TOLD it degraded, or it cannot run the
            // "…and then the GBM learned nothing" check that makes this state
            // visible in the log.
            assert!(
                degraded,
                "seed {seed}: a zeroed column must be reported as degraded to the caller"
            );
        }
    }

    /// `n` candidates that are IDENTICAL in every feature — only `library_id`
    /// (a positional index, not a feature) and the label differ.
    ///
    /// Every lane column is therefore bit-constant, and at `n = 24` — the only
    /// size any caller passes, in this module and above this definition as well as
    /// below it — the variance `ColumnTransform` computes as `sumsq/n - mean^2`
    /// comes out EXACTLY zero, so it culls all 101 of them and
    /// `MlpFoldModel::fit` has nothing left to train on. A GBM handed the same
    /// fixture likewise builds no trees, which is what the detector tests use it
    /// for.
    ///
    /// THE ROW COUNT IS LOAD-BEARING and cannot be scaled freely: that expression
    /// keeps a hair of floating-point residue on a constant column of realistic
    /// magnitude once the row count gets large enough, and the residue clears
    /// `MIN_STD` (1e-12). Measured on this fixture: 8 and 16 train rows cull all
    /// 101 columns, 30 and 32 leave 10 of them alive and the fit SUCCEEDS. So
    /// `n = 24` is chosen to put both paths under the cull —
    /// `CrossValidatedScorer` trains fold `f` on fold `f` alone (8 rows), while
    /// [`crossfit`] trains on all-but-fold-`f` (16 rows) — and raising it to 48
    /// silently stops exercising the hybrid's failure path.
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

    /// FAILURE POLICY, standalone half: a fold that cannot be fitted ABORTS the
    /// run. It does not fall back to a uniform score.
    ///
    /// This is the decision the panic encodes, so it is pinned as a test rather
    /// than left to a comment: the degradation the hybrids use
    /// (`crossfit_score_column_is_uniformly_zero_when_any_fold_fails`) would, on
    /// this path, give every row the same discriminant and therefore q-values
    /// that are a function of the seeded shuffle order alone — meaningless, and
    /// indistinguishable from a real result in the output file. A future reader
    /// "unifying" the two policies has to delete this test to do it.
    ///
    /// The `expected` substring covers the policy (aborted), the entry point and
    /// the LANE. The diagnostic tail is `MlpFoldError`'s own `Display`, pinned in
    /// `mlp_fold`'s `fit_rejects_a_fully_culled_matrix` and
    /// `an_empty_train_slice_reports_zero_rows_rather_than_blaming_the_cull`; the
    /// abort only interpolates it.
    ///
    /// The row count is not free — see [`indistinguishable_competed`].
    #[test]
    #[should_panic(expected = "MLP rescore aborted: Linear lane")]
    fn standalone_mlp_aborts_rather_than_scoring_every_row_the_same() {
        rescore_mlp_lane(
            indistinguishable_competed(24),
            Lane::Linear,
            mlp_test_cfg(7),
        );
    }

    /// The abort's three arms must say DIFFERENT things, and the difference is the
    /// whole point of matching on the variant.
    ///
    /// `scorer.fit()` surfaces all three: `NoUsableColumns` from the fit, and
    /// `WidthMismatch` / `NonFiniteScore` through `assign_scores` -> `get_scores`
    /// -> `predict`. Only the first is a statement about the INPUT that a different
    /// `rescore_model` could work around. A single message for all of them sent an
    /// operator to their config over an internal invariant violation.
    ///
    /// Driven straight at [`abort_standalone_mlp`] because neither predict-side arm
    /// is reachable from a fixture: [`rescore_mlp_lane`] hands the scorer the very
    /// dataset the scorer then predicts against (so a width mismatch would require
    /// the bug the message exists for), and no input can make the net diverge on
    /// demand — `mlp_fold`'s `predict_rejects_a_non_finite_logit_instead_of_
    /// returning_it` covers that half by poisoning a fitted net. The
    /// `NoUsableColumns` arm is covered end to end by
    /// `standalone_mlp_aborts_rather_than_scoring_every_row_the_same`.
    ///
    /// Every panic reaches stderr while this runs — the default hook is left in
    /// place on purpose rather than swapped out, since a `set_hook` here would
    /// silence any concurrently-running test's panic message too.
    #[test]
    fn the_abort_blames_the_config_only_for_the_input_failure() {
        let message = |e: MlpFoldError| -> String {
            let payload = std::panic::catch_unwind(move || {
                abort_standalone_mlp(Lane::Linear, LINEAR_NCOLS, 24, e)
            })
            .expect_err("abort_standalone_mlp diverges");
            payload
                .downcast_ref::<String>()
                .cloned()
                .expect("panic! with a format string carries a String payload")
        };

        let dead_lane = message(MlpFoldError::NoUsableColumns {
            fold: 2,
            ncols: LINEAR_NCOLS,
            train_rows: 16,
        });
        let bad_width = message(MlpFoldError::WidthMismatch {
            fitted: LINEAR_NCOLS,
            got: ALL_NCOLS,
        });
        let diverged = message(MlpFoldError::NonFiniteScore {
            row: 17,
            n_bad: 8,
            scored: 24,
        });

        // Shared: all name the abort, the lane, and the error's own diagnostic.
        for (what, msg) in [
            ("dead lane", &dead_lane),
            ("width", &bad_width),
            ("diverged", &diverged),
        ] {
            assert!(
                msg.contains("MLP rescore aborted: Linear lane")
                    && msg.contains(&format!("{LINEAR_NCOLS} columns"))
                    && msg.contains("24 candidates"),
                "{what}: the abort must place itself (lane, width, row count): {msg}"
            );
        }

        // The INPUT failure points at the config, and names the fold — that comes
        // from `MlpFoldError`'s own Display, which is why the variant carries it.
        assert!(
            dead_lane.contains("rescore_model") && dead_lane.contains("fold 2"),
            "a dead lane is an input problem: say which fold, and offer the way out: {dead_lane}"
        );

        // The INVARIANT violation must NOT, and must say so in as many words.
        assert!(
            bad_width.contains("INTERNAL INVARIANT VIOLATION")
                && bad_width.contains("will not help"),
            "a width mismatch is a code defect, not a config problem, and the message has to \
             say so rather than leaving the reader to infer it: {bad_width}"
        );
        assert!(
            !bad_width.contains("no usable fallback"),
            "the width arm must not reuse the dead-lane advice, which is what sent an \
             operator to their config over a code defect: {bad_width}"
        );
        assert!(
            bad_width.contains(&format!("{ALL_NCOLS} lane columns"))
                && bad_width.contains(&format!("fitted on {LINEAR_NCOLS}")),
            "the two widths are the only actionable content here: {bad_width}"
        );

        // The NUMERIC failure must name the discriminant — the whole reason this
        // is an error rather than a value is that `assign_qval` would have blamed
        // the sort instead.
        assert!(
            diverged.contains("DISCRIMINANT SCORE") && diverged.contains("sorted"),
            "a non-finite logit has to say it IS the score and that the downstream abort \
             would have blamed the sort, or the next reader re-diagnoses it from scratch: \
             {diverged}"
        );
        assert!(
            diverged.contains("dataset row 17") && diverged.contains("8 of 24"),
            "the row and the count come from `MlpFoldError`'s Display, which is why the \
             variant carries them: {diverged}"
        );
        assert!(
            !diverged.contains("INTERNAL INVARIANT VIOLATION")
                && !diverged.contains("no usable fallback"),
            "a diverged fit is neither a code defect nor a dead lane, and reusing either \
             arm's advice would misdirect: {diverged}"
        );
    }

    /// The other half of that decision: the SAME fixture that aborts the
    /// standalone rescorer runs to completion through the hybrid, because there a
    /// failed MLP cross-fit costs one column out of `NONLINEAR_NCOLS + 1` and the
    /// GBM keeps the rest.
    ///
    /// The two paths are deliberately asymmetric and this is the pair of tests
    /// that says so; `abort_standalone_mlp`'s docs carry the reasoning.
    ///
    /// The cross-fit really does fail on this fixture at this size (16 train rows
    /// per fold, all columns culled — see [`indistinguishable_competed`]), so the
    /// `mlp_score` column reaching the GBM below is the all-zero fallback.
    ///
    /// NOTE for the `min_leaf_weight` audit: this runs a GBM on 8 rows per fold,
    /// so it builds NO TREES. That is fine, and deliberate — nothing here claims
    /// the GBM learned anything. Every assertion is about the run COMPLETING and
    /// about the frame it was handed (`feature_stats` names come from the dataset,
    /// not from the model), all of which hold in the no-trees regime. Do not add
    /// an assertion about scores or importance to this test; put it at 360 rows in
    /// `hybrid_mlp_score_carries_the_linear_lane_into_the_gbm` instead.
    #[test]
    fn hybrid_mlp_degrades_where_the_standalone_mlp_aborts() {
        const N: u32 = 24;
        let data = indistinguishable_competed(N);

        // The premise: the MLP cross-fit this hybrid depends on cannot fit here.
        // Without it, the test would pass on a fixture where nothing degraded.
        let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
        let lin = hybrid_linear_dataset(&data, responses);
        assert!(
            crossfit::<_, MlpFoldModel>(&lin, &mlp_test_cfg(7), "MLP").is_none(),
            "fixture premise broken: the MLP cross-fit SUCCEEDS here, so this test is not \
             exercising the degraded path at all"
        );

        let (out, stats) = rescore_hybrid_mlp_with(data, mlp_test_cfg(7));
        assert_eq!(
            out.len(),
            N as usize,
            "the hybrid must still return every row"
        );
        assert!(
            out.iter().all(|r| r.discriminant_score.is_finite()),
            "a failed cross-fit must not put non-finite values in the score column"
        );
        assert!(
            out.iter().all(|r| (0.0..=1.0).contains(&r.qvalue)),
            "q-values must stay in range on the degraded path"
        );
        assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);
        // The frame the GBM was handed still carries the appended column — the
        // degradation zeroes its VALUES, it does not drop the column.
        assert!(
            stats
                .iter()
                .all(|fs| fs.feature_stats.iter().any(|f| &*f.name == "mlp_score")),
            "the degraded hybrid must still hand the GBM an mlp_score column"
        );
    }

    /// Smoke + determinism for the MLP hybrid, mirroring
    /// `rescore_hybrid_smoke_and_determinism` and swept over seeds because a
    /// model is trained: the MLP's initialization is the one thing that varies
    /// between configurations, and a single seed reports an init-dependent
    /// training trap as a clean pass purely on luck.
    ///
    /// Bit-identical is the assertion, not approximately-equal: the whole
    /// q-value pipeline is a ranking of these scores.
    ///
    /// # Why 360 rows, and why the seed-sensitivity check is mandatory
    /// This test was written at 90 and was BLIND. `GBMConfig::default`'s
    /// `min_leaf_weight` of 5.0 is most of a 30-row fold's total logloss
    /// hessian, so the GBM finds no admissible split, builds no trees, and emits
    /// one constant score per fold — a value that does not depend on
    /// `mlp_score`, or on the MLP, at all. Two runs at different MLP SEEDS
    /// produced identical score vectors, which means the bit-identity assertion
    /// below would have passed just as happily if `crossfit_score_column` returned
    /// `rand::random()` per row. That is the one global constraint this test
    /// exists to hold, so being unable to fail was the whole problem.
    ///
    /// 360 rows puts 120 in each training fold, trees get built, and the score
    /// becomes a function of the appended column. The `assert_ne!` at the end is
    /// what keeps that true: it fails the moment the output stops depending on
    /// the MLP, i.e. the moment every assertion above goes vacuous again. Do not
    /// drop it, and do not shrink `n` without checking it still fails.
    #[test]
    fn rescore_hybrid_mlp_smoke_and_determinism() {
        let n = 360;
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

        // Kept from the sweep so the non-vacuity check below can reuse this run
        // instead of training a thirteenth net for a result the loop already
        // computed bit-identically.
        let mut first: Option<(u64, Vec<(u32, u32)>)> = None;
        for seed in [7u64, 13, 42, 1234] {
            let run = || rescore_hybrid_mlp_with(synthetic_competed(n), mlp_test_cfg(seed));
            let (out_a, stats_a) = run();
            let (out_b, _) = run();

            assert_eq!(out_a.len(), n as usize);
            for r in &out_a {
                assert!(
                    r.discriminant_score.is_finite(),
                    "seed {seed}: non-finite discriminant score {}",
                    r.discriminant_score
                );
                assert!(
                    (0.0..=1.0).contains(&r.qvalue),
                    "seed {seed}: qvalue out of [0,1]: {}",
                    r.qvalue
                );
            }
            assert_eq!(
                key(&out_a),
                key(&out_b),
                "seed {seed}: mlp hybrid not deterministic across runs"
            );
            assert_eq!(stats_a.len(), N_RESCORE_FOLDS as usize);
            first.get_or_insert_with(|| (seed, key(&out_a)));
        }

        // NON-VACUITY, and the load-bearing assertion of this test: the output
        // must actually DEPEND on the MLP. Two different MLP seeds train
        // different nets, so they produce a different `mlp_score` column and
        // must produce different final scores. If this passes, the bit-identity
        // assertions above are measuring reproducibility of the MLP; if it
        // fails, they are measuring reproducibility of a constant.
        //
        // Any two of the sweep's seeds would do; this pairs the sweep's first
        // against a seed the sweep did not use, so a hypothetical collision
        // between two swept seeds cannot mask the check.
        let (swept_seed, swept) = first.expect("the sweep above ran at least once");
        let novel = key(&rescore_hybrid_mlp_with(synthetic_competed(n), mlp_test_cfg(999)).0);
        assert_ne!(
            swept, novel,
            "seed {swept_seed} and seed 999 gave bit-identical hybrid scores, so the output \
             does not depend on mlp_score and every determinism assertion in this test is \
             vacuous (at 90 rows this was exactly the case: no trees, one constant \
             score per fold)"
        );

        // The public entry point, on the DEFAULT `MlpConfig` rather than the
        // shrunk test one — the config the production path actually uses.
        let (out_a, _) = rescore_hybrid_mlp(synthetic_competed(n));
        let (out_b, _) = rescore_hybrid_mlp(synthetic_competed(n));
        assert_eq!(out_a.len(), n as usize);
        assert!(out_a.iter().all(|r| r.discriminant_score.is_finite()));
        assert_eq!(
            key(&out_a),
            key(&out_b),
            "rescore_hybrid_mlp is not deterministic on the default config"
        );
    }

    /// The GBM in the MLP hybrid trains on `nonlinear + mlp_score` — NOT on the
    /// all lane, and not on the nonlinear lane alone.
    ///
    /// Nothing type-checks the feature set: every lane produces a
    /// correctly-shaped matrix and plausible q-values, so a hybrid that quietly
    /// dropped the appended column, appended it in the wrong position, or handed
    /// the GBM the full lane set would still run. What is recoverable from the
    /// outside is the sidecar, and three things are pinned through it:
    ///
    /// * the per-fold stats are FULL WIDTH at `NONLINEAR_NCOLS + 1` and the
    ///   trailing column is named `mlp_score`, i.e. the column is appended, not
    ///   inserted;
    /// * no fold reports an IMPORTANCE name outside `nonlinear ++ mlp_score`.
    ///
    /// THE REAL GUARD IS THE WIDTH CHECK, and the distinction matters. The stats
    /// names come from the dataset's columns, so they describe the frame the GBM
    /// was HANDED whatever it learned; the importance names come from the model,
    /// and on this fixture the model reports exactly one column. Measured at 360
    /// rows: every fold's importance is `[("mlp_score", 1.0)]` and nothing else.
    /// So the importance loop never inspects a lane name at all — it is a cheap
    /// tripwire for a future fixture, not the assertion carrying this test.
    ///
    /// That is a consequence of `synthetic_competed`'s NONLINEAR LANE BEING
    /// BIT-CONSTANT (see its doc): `mlp_score` is the only column here the GBM
    /// can split on, so it is the only one it reports. Whether that column is
    /// load-bearing rather than merely present is still
    /// `hybrid_mlp_score_carries_the_linear_lane_into_the_gbm`'s job — that test
    /// asserts on the resulting AUC, which this one does not look at.
    #[test]
    fn hybrid_mlp_gbm_trains_on_the_nonlinear_lane_plus_mlp_score() {
        let nonlinear: std::collections::HashSet<Arc<str>> =
            nonlinear_feature_name_set().into_iter().collect();
        let linear: std::collections::HashSet<Arc<str>> =
            linear_feature_name_set().into_iter().collect();

        // 360 rows, not 90: `GBMConfig::default`'s `min_leaf_weight` of 5.0 is
        // most of a small fold's total logloss hessian, so a GBM trained on 30
        // rows finds no admissible split, builds no trees, and reports an EMPTY
        // importance — under which the importance-name assertions below would all
        // pass vacuously. `at_least_one_reported` pins that. (It does not make
        // them non-vacuous in the stronger sense: the only name they ever see on
        // this fixture is `mlp_score`. See the doc above.)
        let mut at_least_one_reported = false;
        for seed in [7u64, 13, 42, 1234] {
            let (_out, stats) =
                rescore_hybrid_mlp_with(synthetic_competed(360), mlp_test_cfg(seed));
            assert_eq!(stats.len(), N_RESCORE_FOLDS as usize);

            for fs in &stats {
                at_least_one_reported |= !fs.feature_importance.is_empty();
                assert_eq!(
                    fs.feature_stats.len(),
                    NONLINEAR_NCOLS + 1,
                    "seed {seed}: fold {} must be the nonlinear lane plus ONE column",
                    fs.fold
                );
                assert_eq!(
                    &*fs.feature_stats[NONLINEAR_NCOLS].name, "mlp_score",
                    "seed {seed}: fold {} — mlp_score must be the TRAILING column",
                    fs.fold
                );

                for (name, gain) in &fs.feature_importance {
                    assert!(
                        gain.is_finite(),
                        "seed {seed}: fold {} emitted the NAN sentinel",
                        fs.fold
                    );
                    assert!(
                        &**name == "mlp_score" || nonlinear.contains(name),
                        "seed {seed}: fold {} reported {name}, which is neither a \
                         nonlinear-lane feature nor mlp_score — the GBM was handed the \
                         wrong feature set{}",
                        fs.fold,
                        if linear.contains(name) {
                            " (it is a LINEAR-lane feature, i.e. the all lane leaked in)"
                        } else {
                            ""
                        }
                    );
                }
            }
        }
        assert!(
            at_least_one_reported,
            "no fold reported any importance at all, so every name assertion above \
             passed vacuously"
        );
    }

    /// The premise BOTH `*_carries_the_linear_lane_into_the_gbm` tests rest on:
    /// not one nonlinear-lane column varies across `fixture`'s rows, so a GBM
    /// handed the nonlinear lane has literally nothing to split on and cannot
    /// rank better than chance. Anything above chance downstream therefore had to
    /// arrive through the appended `column` and nowhere else.
    ///
    /// Checked rather than assumed, and checked per-feature rather than by a
    /// count, so a future nonlinear feature that happens to vary with something
    /// `synthetic_competed_linear_only` sets fails loudly here instead of quietly
    /// giving the GBM a second signal and turning both tests into "the GBM can
    /// separate something".
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

    /// `lda_score` is not decoration: it is the channel through which the LINEAR
    /// lane reaches a GBM that is never shown it.
    ///
    /// The `rescore_hybrid` counterpart of
    /// [`hybrid_mlp_score_carries_the_linear_lane_into_the_gbm`], and the only
    /// test that observes whether that (shipping) path's `lda_score` column is
    /// CORRECT rather than merely present.
    ///
    /// `rescore_hybrid_smoke_and_determinism` gets part of the way and cannot
    /// finish the job. It runs at 360 rows on `synthetic_competed`, whose
    /// nonlinear lane is bit-constant, so `lda_score` is the GBM's only
    /// splittable column and its non-vacuity check does catch a fully DEAD one:
    /// measured, a zeroed column takes that fixture's per-fold importance from
    /// `[1, 1, 1]` to `[0, 0, 0]` and the assertion fires. What it cannot catch is
    /// a column that is alive but WRONG — row-permuted, truncated, or filled with
    /// noise — because its only other assertions are that two runs of the same
    /// function agree with each other. This fixture makes `lda_score` the sole
    /// non-constant column and then asserts on the resulting AUC, so it is a
    /// statement about the column's VALUES.
    ///
    /// No seed sweep: the LDA is a closed-form solve with no initialization to
    /// vary, so there is no seed for the "sweep seeds where a model is trained"
    /// rule to sweep — a rerun is bit-identical by construction, which
    /// `rescore_lda_is_cross_fit_and_deterministic` already pins.
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

    /// `mlp_score` is not decoration: it is the channel through which the LINEAR
    /// lane reaches a GBM that is never shown it.
    ///
    /// `synthetic_competed_linear_only` puts ALL the label signal in the linear
    /// lane and leaves the nonlinear lane bit-for-bit constant across rows (see
    /// [`assert_nonlinear_lane_is_flat`], which checks it). A GBM handed the
    /// nonlinear lane alone therefore has nothing to split on and cannot rank
    /// better than chance. So on this fixture:
    ///
    /// * `mlp_score` MUST appear in some fold's importance — forust omits the
    ///   columns it never split on, and a constant column (all-zero from a
    ///   silently failed cross-fit, say) has zero split gain and would be absent
    ///   from every fold;
    /// * and the hybrid's output must actually SEPARATE, which is only reachable
    ///   through that one column.
    ///
    /// The second assertion is the one that cannot be satisfied by accident: it
    /// fails for a dropped column, a constant column, a column filled from the
    /// wrong row order, and a cross-fit whose scores are noise.
    #[test]
    fn hybrid_mlp_score_carries_the_linear_lane_into_the_gbm() {
        assert_nonlinear_lane_is_flat(&synthetic_competed_linear_only(360), "mlp_score");

        for seed in [7u64, 13, 42, 1234] {
            let (out, stats) =
                rescore_hybrid_mlp_with(synthetic_competed_linear_only(360), mlp_test_cfg(seed));

            let split_on_it = stats
                .iter()
                .filter(|fs| {
                    fs.feature_importance
                        .iter()
                        .any(|(n, _)| &**n == "mlp_score")
                })
                .count();
            assert!(
                split_on_it > 0,
                "seed {seed}: no fold split on mlp_score even though it is the only \
                 non-constant column, so the appended column never reached the GBM \
                 or is constant (a silently failed cross-fit)"
            );

            let scores: Vec<f64> = out.iter().map(|r| r.discriminant_score as f64).collect();
            let is_target: Vec<bool> = out.iter().map(|r| r.scoring.identity.is_target).collect();
            let auc = pair_auc(&scores, &is_target);
            assert!(
                auc > 0.8,
                "seed {seed}: AUC {auc} — the GBM sees a constant nonlinear lane here, so \
                 anything above chance has to come through mlp_score, and it did not"
            );
        }
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
