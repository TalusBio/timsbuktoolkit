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
/// Shared by ALL three rescorers so the row order — and therefore the
/// positional fold partition `fold(i) = i % N_RESCORE_FOLDS` derived from it —
/// has exactly one definition.
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

/// The shared tail of ALL three rescorers: rank by discriminant score, derive
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
/// [`crossfit_mlp_scores`] and therefore `rescore_lda` / `rescore_hybrid` /
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
/// is holding out — leak restored, silently. That agreement used to be the same
/// modulo written out in two modules with nothing enforcing the match; both
/// sides now read it from `RowMajorDataset::get_fold`, so there is one
/// definition rather than two that happen to coincide. The
/// `crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on` test pins the
/// two partitions against each other rather than against a repeated modulo.
///
/// # Failure policy — all-or-nothing
/// Returns `None` if ANY fold's fit or scoring fails (singular scatter / empty
/// class for the LDA, a fully culled train matrix for the MLP). Per-fold
/// fallback was rejected on purpose: leaving one fold's rows at 0.0 while the
/// others carry real discriminant values makes the score distribution
/// FOLD-DEPENDENT, i.e. a function of a row's position in the shuffle rather
/// than of its evidence. That silently corrupts the q-value ranking (a whole
/// third of the rows is pushed into an arbitrary slab of the sort) in a way that
/// no downstream check would catch. A uniform failure is both louder and
/// harmless to the ranking, so callers degrade the WHOLE run instead. Callers
/// log the failure; this function logs which fold failed.
///
/// `scores` is local to this call and only reaches the caller through the
/// `Some` arm, so a partially filled column is not something a caller can
/// obtain even by ignoring the `None`.
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

        // `val` is empty: neither model here early-stops (the LDA is
        // closed-form, the MLP runs a fixed epoch count), and `FoldModel`
        // documents that such a model ignores the slice.
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

/// The `mlp_score` column [`rescore_hybrid_mlp`] appends to the nonlinear lane:
/// each row's HELD-OUT MLP score over the canonical rescore partition (see
/// [`crossfit`]), or a UNIFORMLY ZERO column if the cross-fit failed.
///
/// Always `data.nrows()` long, and never partially filled — the fallback is a
/// fresh vector rather than whatever [`crossfit`] had accumulated before the
/// failing fold, so a run where folds 0 and 1 fitted and fold 2 did not still
/// yields all zeros. A partially filled column is the failure mode worth naming:
/// it makes a row's `mlp_score` depend on its position in the shuffle, gives the
/// GBM a feature it can split on to recover fold identity, and corrupts the
/// q-value ranking with nothing downstream to catch it. All-zero is constant ->
/// zero split gain -> the hybrid degrades to "GBM on the nonlinear lane", which
/// is sane and loud.
///
/// The models are discarded: the hybrid's sidecar reports the GBM's importance
/// over `nonlinear + mlp_score`, and the MLP's own per-column weights belong to
/// a different feature set (the LINEAR lane) that the sidecar has no column for.
/// [`rescore_mlp_linear`] is the entry point that reports those.
fn crossfit_mlp_scores<D: FoldDataset>(data: &D, cfg: &MlpConfig) -> Vec<f64> {
    match crossfit::<D, MlpFoldModel>(data, cfg, "MLP") {
        Some(cf) => cf.scores,
        None => {
            tracing::error!(
                "hybrid: cross-fit MLP failed; mlp_score is uniformly 0 (GBM falls back to the \
                 nonlinear lane alone)"
            );
            vec![0.0f64; data.nrows()]
        }
    }
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    let config = GBMConfig::default();

    canonicalize_and_shuffle(&mut data);

    // Build the ALL-lane matrix (linear ++ nonlinear) over the shuffled `data`,
    // so row `i` aligns with `data[i]` — fold assignment + labels are
    // positional, so the matrix MUST be built AFTER the shuffle. The GBM sees
    // the full lane feature set; feature names come from the same walk order,
    // so columns and `names` align by construction (width asserted below + the
    // lane-parity test).
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

    finalize(scorer.score(), stats)
}

/// Sage-style shrinkage-LDA rescorer: closed-form linear fits, no boosting —
/// ~100x cheaper than the GBM path. The FDR machinery (`assign_qval`,
/// target-decoy competition) is untouched — only the discriminant score source
/// changes.
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

    // Canonical sort + seeded shuffle — IDENTICAL to `rescore` / `rescore_hybrid`.
    // The shuffle is load-bearing here now that the fit is cross-fit: the fold
    // partition below is positional.
    canonicalize_and_shuffle(&mut data);

    // LDA trains on the LINEAR lane only: fields that are approx-Gaussian after
    // their declared per-row transform (raw/log2/ln1p). Skew-taming is done at
    // emit time by the grammar, so there is no data-dependent normalization step
    // here — the only remaining data-dependent op is LDA's own standardization.
    // The matrix is built AFTER the shuffle over `data` in its current order, so
    // row `i` aligns with `data[i]` — fold assignment and labels are positional;
    // `names` comes from the same walk order, so columns and names align by
    // construction (width asserted below + the parity test).
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
/// must come from the same order — labels and folds are positional.
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
/// CV on `nonlinear + lda_score`. The GBM re-sees ~30 features instead of the
/// full ~131 (the compression play) at ~parity.
///
/// LEAK-FREEDOM: `lda_score` is cross-fit via [`crossfit`] — see there for the
/// partition, why a label-aware feature fed to a CV'd GBM in particular must be
/// leak-free, and why that partition has to match the one
/// `CrossValidatedScorer` builds internally.
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

    let nrows = data.len();

    // Build BOTH lane matrices AFTER the shuffle over the SAME `data`, so row
    // `i` is the same candidate in lin, nl, lda_score, responses, and the moved
    // `data`.
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let lin_dataset = hybrid_linear_dataset(&data, responses.clone());

    // --- CROSS-FIT lda_score (leak-free), via the shared partition ---
    // On failure the column is left uniformly zero rather than partially
    // filled: a fold-dependent `lda_score` is a feature the GBM can split on to
    // recover fold identity, which is strictly worse than no feature at all.
    // An all-zero column is constant -> zero split gain -> hybrid degrades to
    // "GBM on the nonlinear lane", which is a sane, loud degradation.
    let lda_score = match crossfit_lda(&lin_dataset) {
        Some(cf) => cf.scores,
        None => {
            tracing::error!(
                "hybrid: cross-fit LDA failed; lda_score is uniformly 0 (GBM falls back to the \
                 nonlinear lane alone)"
            );
            vec![0.0f64; nrows]
        }
    };

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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
/// UNCHANGED from the GBM path, because it is literally the same scorer: fold
/// `f`'s model trains on fold `f`, early-stops on `f + 1` (which
/// [`MlpFoldModel`] ignores — it has no early stopping), and scores the
/// remaining `n_folds - 2` folds; each row's score is the mean over the models
/// that saw neither its fold nor the fold before it. No row is ever scored by a
/// model that trained on it, and every fitted statistic (cull set,
/// standardization moments, imputation means, weights) is train-fold-only
/// inside [`MlpFoldModel::fit`].
///
/// # Determinism
/// [`canonicalize_and_shuffle`] pins the row order and therefore the positional
/// fold partition; each fold's RNG is `MlpConfig::rng_for_fold(fold)`, a pure
/// function of the config seed and the fold index. Nothing else is stochastic,
/// so two runs of the same build on the same input produce bit-identical
/// scores.
fn rescore_mlp_lane(
    mut data: Vec<CompetedCandidate>,
    lane: Lane,
    config: MlpConfig,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    canonicalize_and_shuffle(&mut data);

    // ORDER IS LOAD-BEARING: the matrix is built here, AFTER the shuffle, over
    // `data` in its current order, so row `i` is `data[i]`. Fold assignment
    // (`i % N_RESCORE_FOLDS`) and labels are both positional, so building it
    // before the shuffle would pair every row's features with another
    // candidate's label and fold — a silent correctness bug with no type or
    // width check that could notice. `names` comes from the same lane, so
    // columns and names align by construction (width asserted below + the lane
    // parity test).
    let names = lane.names();
    let ncols = lane.ncols();
    debug_assert_eq!(names.len(), ncols);
    let feat = lane.matrix(&data);
    debug_assert_eq!(feat.len(), data.len() * ncols);

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
    // Same failure policy as the GBM path: a rescorer that cannot fit has no
    // meaningful score to fall back to, and a partially-scored run would be
    // fold-dependent (see `crossfit`'s all-or-nothing note).
    if let Err(e) = scorer.fit() {
        panic!("MLP rescore ({lane:?} lane) failed: {e}");
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
// Intended mapping (the enum variants themselves are NOT this task's scope):
//   `RescoreModel::Mlp`    -> `rescore_mlp_linear` (the default MLP shape)
//   `RescoreModel::MlpAll` -> `rescore_mlp_all`

/// MLP rescorer over the LINEAR lane only — the same feature set
/// [`rescore_lda`] uses, so it isolates "nonlinear model" from "more features".
/// THE DEFAULT MLP SHAPE, i.e. the one `RescoreModel::Mlp` is meant to select.
///
/// See [`rescore_mlp_lane`] for the cross-fit and determinism contracts.
pub fn rescore_mlp_linear(data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    rescore_mlp_lane(data, Lane::Linear, MlpConfig::default())
}

/// MLP rescorer over the ALL lane (linear ++ nonlinear) — the same feature set
/// [`rescore`] gives the GBM, so the two are directly comparable. The one
/// `RescoreModel::MlpAll` is meant to select.
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

    // BOTH lane matrices are built AFTER the shuffle, over the SAME `data`, so
    // row `i` is the same candidate in the linear dataset, the nonlinear matrix,
    // `mlp_score`, `responses`, and the moved `data`. Fold assignment and labels
    // are positional, so building either one earlier would pair every feature
    // row with another candidate's label and fold.
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let lin_dataset = hybrid_linear_dataset(&data, responses.clone());

    // Uniformly zero if the cross-fit failed, never partially filled — see
    // `crossfit_mlp_scores`.
    let mlp_score = crossfit_mlp_scores(&lin_dataset, &mlp_config);

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

    finalize(scorer.score(), stats)
}

/// The MLP hybrid on [`MlpConfig::default`] — see [`rescore_hybrid_mlp_with`].
///
/// The one `RescoreModel::HybridMlp` is meant to select (the variant itself is
/// not this function's scope).
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
// layout `LdaModel::fit` and `PrecomputedFeatures::from_row_major` both index
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
                // Also vary the meta delta-group (nonlinear lane) for GBM signal.
                c.delta_group = if is_target { 2.0 } else { 0.5 } + (i % 7) as f32 * 0.1;
                c.delta_group_ratio = if is_target { 0.8 } else { 0.3 };
                c
            })
            .collect()
    }

    /// [`synthetic_competed`] with the NONLINEAR lane flattened: the counts keep
    /// their label-correlated variance, but `delta_group` / `delta_group_ratio`
    /// are left at the sample defaults, so every nonlinear-lane column is
    /// constant across rows.
    ///
    /// Built for the hybrid tests, where "the GBM cannot separate these rows
    /// without the appended column" has to be a fact about the fixture rather
    /// than a hope. The nonlinear lane being FLAT is asserted at the point of use,
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
        let insample = LdaModel::fit(&feat, N, 3, &is_decoy, DEFAULT_SHRINKAGE).unwrap();
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
        let n = 360;
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
    /// with the batch size cut so the fixtures here (90-120 rows) actually get
    /// several minibatch steps per epoch instead of one full-batch step. The
    /// default is left alone on purpose — tuning a production default for test
    /// runtime is how a default stops meaning anything.
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
    #[test]
    fn rescore_mlp_is_deterministic_on_both_lanes() {
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

                // One FoldStats per fold, full width, non-empty on both halves
                // — the LDA path's assertions.
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
                    "{lane:?}/{seed}: AUC {auc} is near chance — the feature rows are \
                     paired with the wrong candidates' labels, i.e. the lane matrix was \
                     built before the shuffle rather than after it"
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
    #[test]
    fn rescore_with_dispatches_each_mlp_variant_to_its_lane() {
        use crate::ml::{
            RescoreModel,
            rescore_with,
        };

        let expected_set = |names: Vec<Arc<str>>| -> std::collections::HashSet<Arc<str>> {
            names.into_iter().collect()
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
    #[test]
    fn crossfit_holds_out_exactly_the_rows_the_gbm_scorer_trains_on() {
        const N: usize = 47;
        const NCOLS: usize = 2;
        let n_folds = N_RESCORE_FOLDS as usize;

        let feat: Vec<f64> = (0..N).flat_map(|i| [i as f64, (i % 2) as f64]).collect();
        let responses: Vec<f64> = (0..N).map(|i| (i % 2) as f64).collect();
        let names: Vec<Arc<str>> = ["col0", "col1"].map(Arc::from).to_vec();

        // Side A: the cross-fit that produces `mlp_score` / `lda_score`.
        let dataset = RowMajorDataset::new(
            PrecomputedFeatures::from_row_major(feat.clone(), NCOLS, responses.clone()),
            names.clone(),
            n_folds,
        );
        let cf = crossfit::<_, FoldSpy>(&dataset, &(), "spy").expect("the spy cannot fail");

        // Side B: the GBM scorer the hybrid hands that column to. No fit needed
        // — the partition is built in the constructor.
        let rows: Vec<LabelRow> = responses
            .iter()
            .map(|&y| LabelRow { y, score: 0.0 })
            .collect();
        let scorer =
            CrossValidatedScorer::<LabelRow, GbmFoldModel>::new_from_shuffled_with_precomputed(
                N_RESCORE_FOLDS,
                rows,
                GBMConfig::default(),
                PrecomputedFeatures::from_row_major(feat, NCOLS, responses),
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
    fn crossfit_mlp_scores_is_uniformly_zero_when_any_fold_fails() {
        use crate::ml::mlp_fold::MlpFoldError;

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
                Some(MlpFoldError::NoUsableColumns),
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
            let col = crossfit_mlp_scores(&dataset, &cfg);
            assert_eq!(
                col.len(),
                N,
                "seed {seed}: the column must still be full length"
            );
            assert!(
                col.iter().all(|s| s.to_bits() == 0.0f64.to_bits()),
                "seed {seed}: mlp_score must be uniformly zero after a failed fold, got {col:?}"
            );
        }
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
    /// below would have passed just as happily if `crossfit_mlp_scores` returned
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
        }

        // NON-VACUITY, and the load-bearing assertion of this test: the output
        // must actually DEPEND on the MLP. Two different MLP seeds train
        // different nets, so they produce a different `mlp_score` column and
        // must produce different final scores. If this passes, the bit-identity
        // assertions above are measuring reproducibility of the MLP; if it
        // fails, they are measuring reproducibility of a constant.
        let by_seed =
            |seed: u64| key(&rescore_hybrid_mlp_with(synthetic_competed(n), mlp_test_cfg(seed)).0);
        assert_ne!(
            by_seed(7),
            by_seed(999),
            "two different MLP seeds gave bit-identical hybrid scores, so the output \
             does not depend on mlp_score and every determinism assertion in this \
             test is vacuous (at 90 rows this was exactly the case: no trees, one \
             constant score per fold)"
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
    /// * no LINEAR-lane name is reported, i.e. the GBM did not get the all lane
    ///   (the MLP consumed the linear lane and handed over one column);
    /// * no fold reports a name outside `nonlinear ++ mlp_score` at all.
    ///
    /// That the appended column is not merely PRESENT but load-bearing is
    /// `hybrid_mlp_score_carries_the_linear_lane_into_the_gbm`'s job — this
    /// fixture's nonlinear lane separates targets from decoys on its own, so the
    /// GBM here has no need to split on `mlp_score` and its absence from the
    /// importance would prove nothing.
    #[test]
    fn hybrid_mlp_gbm_trains_on_the_nonlinear_lane_plus_mlp_score() {
        let nonlinear: std::collections::HashSet<Arc<str>> =
            nonlinear_feature_name_set().into_iter().collect();
        let linear: std::collections::HashSet<Arc<str>> =
            linear_feature_name_set().into_iter().collect();

        // 360 rows, not 90: `GBMConfig::default`'s `min_leaf_weight` of 5.0 is
        // most of a small fold's total logloss hessian, so a GBM trained on 30
        // rows finds no admissible split, builds no trees, and reports an EMPTY
        // importance — under which the name assertions below would all pass
        // vacuously. `at_least_one_reported` pins that.
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
    /// `rescore_hybrid_smoke_and_determinism` cannot: its fixture separates on the
    /// nonlinear lane too, so it notices a dead `lda_score` only as a thinner
    /// importance report, and a row-permuted or truncated column not at all. This
    /// fixture leaves `lda_score` as the only non-constant column, so the AUC
    /// assertion below is a statement about the column's VALUES.
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
