use super::cv::{
    CrossValidatedScorer,
    DataBuffer,
    FeatureLike,
    FeatureStat,
    FoldStats,
    GBMConfig,
    PrecomputedFeatures,
    RescoreFeatureStats,
};
use super::lda::{
    DEFAULT_SHRINKAGE,
    LdaModel,
};
use super::{
    LabelledScore,
    N_RESCORE_FOLDS,
    TargetDecoy,
};
use crate::scoring::blocks::derived::Derived;
use crate::scoring::blocks::result_meta::ResultMeta;
use crate::scoring::blocks::{
    FeatFrame,
    FrameSink,
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

/// A cross-fit (leak-free) LDA over the canonical rescore fold partition.
///
/// See `crossfit_lda` for the partition contract.
struct CrossFitLda {
    /// Held-out score per row, row-aligned with the input matrix.
    scores: Vec<f64>,
    /// The per-fold fitted models, indexed by fold. Always
    /// `N_RESCORE_FOLDS` long — a `CrossFitLda` only exists if every fold fit.
    models: Vec<LdaModel>,
}

/// Cross-fit an [`LdaModel`] over the canonical rescore fold partition and
/// return each row's HELD-OUT score.
///
/// THE partition: `fold(i) = i % N_RESCORE_FOLDS`. This is the single
/// definition used by `rescore_lda` and `rescore_hybrid`, and it MUST equal the
/// partition `CrossValidatedScorer` builds internally
/// (`assigned_fold[i] = i % n_folds`, see
/// `ml::cv::CrossValidatedScorer::new_from_shuffled`). For fold `f` the model is
/// fitted on every row `i` with `i % N_RESCORE_FOLDS != f` and then scores only
/// the rows with `i % N_RESCORE_FOLDS == f`, so no row ever contributes to the
/// model that scores it.
///
/// # Failure policy — all-or-nothing
/// Returns `None` if ANY fold's fit fails (singular scatter / empty class).
/// Per-fold fallback was rejected on purpose: leaving one fold's rows at 0.0
/// while the others carry real discriminant values makes the score
/// distribution FOLD-DEPENDENT, i.e. a function of a row's position in the
/// shuffle rather than of its evidence. That silently corrupts the q-value
/// ranking (a whole third of the rows is pushed into an arbitrary slab of the
/// sort) in a way that no downstream check would catch. A uniform failure is
/// both louder and harmless to the ranking, so callers degrade the WHOLE run
/// instead. Callers log the failure; this function logs which fold failed.
fn crossfit_lda(
    feat: &[f64],
    nrows: usize,
    ncols: usize,
    is_decoy: &[bool],
) -> Option<CrossFitLda> {
    let n_folds = N_RESCORE_FOLDS as usize;
    let mut scores = vec![0.0f64; nrows];
    let mut models = Vec::with_capacity(n_folds);

    // Reused across folds: the gathered training matrix is the only sizeable
    // allocation here and its capacity is fold-invariant.
    let mut train: Vec<f64> = Vec::with_capacity(feat.len());
    let mut train_is_decoy: Vec<bool> = Vec::with_capacity(nrows);

    for f in 0..n_folds {
        // TRAIN rows = every row NOT in fold f, gathered contiguously.
        train.clear();
        train_is_decoy.clear();
        for i in 0..nrows {
            if i % n_folds != f {
                train.extend_from_slice(&feat[i * ncols..(i + 1) * ncols]);
                train_is_decoy.push(is_decoy[i]);
            }
        }
        let train_nrows = train_is_decoy.len();
        let model = match LdaModel::fit(
            &train,
            train_nrows,
            ncols,
            &train_is_decoy,
            DEFAULT_SHRINKAGE,
        ) {
            Some(m) => m,
            None => {
                tracing::error!(
                    "cross-fit LDA: fold {f}/{n_folds} fit failed (singular or empty class)"
                );
                return None;
            }
        };
        // Score ONLY the held-out fold, with a model that never saw it.
        for i in (f..nrows).step_by(n_folds) {
            scores[i] = model.score(&feat[i * ncols..(i + 1) * ncols]);
        }
        models.push(model);
    }

    Some(CrossFitLda { scores, models })
}

#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub fn rescore(mut data: Vec<CompetedCandidate>) -> (Vec<FinalResult>, RescoreFeatureStats) {
    let config = GBMConfig::default();

    canonicalize_and_shuffle(&mut data);

    // Build the ALL-lane frame (linear ++ nonlinear) over the shuffled `data`,
    // so frame row `i` aligns with `data[i]` — fold assignment + labels are
    // positional, so the frame MUST be built AFTER the shuffle. The GBM sees
    // the full lane feature set; feature names come from the same walk, so the
    // frame columns and `names` align by construction (asserted below + the
    // lane-parity test).
    let names = all_feature_name_set_for(&data);
    let frame = build_all_frame(&data);
    debug_assert_eq!(
        frame.names(),
        names.as_slice(),
        "ALL-lane frame columns must match all_feature_name_set order"
    );
    let ncols = frame.ncols();
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let precomputed = PrecomputedFeatures::from_row_major(frame.row_major(), ncols, responses);

    let mut scorer = CrossValidatedScorer::<CompetedCandidate>::new_from_shuffled_with_precomputed(
        N_RESCORE_FOLDS,
        data,
        config,
        precomputed,
    );
    scorer
        .fit(&mut DataBuffer::default(), &mut DataBuffer::default())
        .unwrap();

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

/// Sage-style shrinkage-LDA rescorer: closed-form linear fits, no boosting —
/// ~100x cheaper than the GBM path. The FDR machinery (`assign_qval`,
/// target-decoy competition) is untouched — only the discriminant score source
/// changes.
///
/// CROSS-FIT, not a single in-sample fit. Every row's score comes from an LDA
/// fitted WITHOUT that row, via the shared `crossfit_lda` partition
/// (`fold(i) = i % N_RESCORE_FOLDS`). An LDA is label-aware, so a single fit
/// over all rows followed by scoring those same rows would let each row's
/// discriminant peek at its own label — exactly the leak `rescore_hybrid`
/// avoids — and the resulting target/decoy separation (hence the q-values that
/// `assign_qval` derives from it) would be optimistically biased. `N_RESCORE_FOLDS`
/// fits instead of 1 is affordable precisely because the LDA is cheap.
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
    // The frame is built AFTER the shuffle over `data` in its current order, so
    // row `i` aligns with `data[i]` — fold assignment and labels are positional;
    // `names` comes from the same walk, so frame columns and names align by
    // construction (asserted below + the parity test).
    let names: Vec<Arc<str>> = linear_feature_name_set_for(&data);
    let nrows = data.len();

    // Materialize the linear-lane matrix once as a single flat row-major buffer
    // (`feat[i*ncols + j]`) — the layout LDA wants; the column store transposes
    // in one pass.
    let t = Instant::now();
    let frame = build_lane_frame(&data, Lane::Linear);
    debug_assert_eq!(
        frame.names(),
        names.as_slice(),
        "linear-lane frame columns must match linear_feature_name_set order"
    );
    let ncols = frame.ncols();
    let feat = frame.row_major();
    debug_assert_eq!(feat.len(), nrows * ncols);
    let is_decoy: Vec<bool> = data.iter().map(|c| c.get_y() < 0.5).collect();

    // Optional raw-matrix dump for offline feature-engineering iteration.
    // `TIMSSEEK_LDA_DUMP=/prefix` writes `<prefix>.f64` (row-major matrix),
    // `<prefix>.labels` (u8, 1=target), `<prefix>.names.txt`. This is the
    // LINEAR-lane matrix + linear names (offline python reads this).
    if let Ok(prefix) = std::env::var("TIMSSEEK_LDA_DUMP") {
        dump_feature_matrix(&prefix, &feat, nrows, &names, &is_decoy);
    }
    eprintln!(
        "  LDA: extracted {nrows} x {ncols} linear-lane matrix in {:.2?}",
        t.elapsed()
    );

    let t = Instant::now();
    let stats = match crossfit_lda(&feat, nrows, ncols, &is_decoy) {
        Some(cf) => {
            for (cand, &s) in data.iter_mut().zip(cf.scores.iter()) {
                cand.assign_score(s);
            }
            eprintln!(
                "  LDA: cross-fit ({N_RESCORE_FOLDS} folds) + scored {nrows} candidates in {:.2?}",
                t.elapsed()
            );
            lda_feature_stats(&names, &feat, nrows, ncols, &cf.models)
        }
        None => {
            // All-or-nothing (see `crossfit_lda`): a partial cross-fit would
            // leave one fold's rows at 0.0 and the rest at real discriminant
            // values, making the score distribution depend on shuffle position.
            // Zeroing every row keeps the failure uniform — q-values collapse to
            // a single uninformative value instead of silently mis-ranking a
            // third of the data.
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

    let mut scored = data;
    #[cfg(feature = "rayon")]
    scored.par_sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    #[cfg(not(feature = "rayon"))]
    scored.sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut scored, |x| CompetedCandidate::get_score(x) as f32);

    (scored.into_iter().map(|c| c.into_final()).collect(), stats)
}

/// Hybrid rescorer: cross-fit an LDA on the LINEAR lane, push its (leak-free)
/// `lda_score` as one extra column into the NONLINEAR lane, then train the GBM
/// CV on `nonlinear + lda_score`. The GBM re-sees ~30 features instead of the
/// full ~131 (the compression play) at ~parity.
///
/// LEAK-FREEDOM: `lda_score` is label-aware (LDA fits on target/decoy labels).
/// A single LDA fit on all rows fed to the CV'd GBM would let each held-out
/// fold's `lda_score` peek at its own labels -> optimistic FDR. Instead the
/// column comes from the shared `crossfit_lda` helper, which owns the ONE
/// definition of the partition (`fold(i) = i % N_RESCORE_FOLDS`) — the same
/// helper `rescore_lda` uses. That partition MUST match `CrossValidatedScorer`'s
/// internal `assigned_fold[i] = i % n_folds`
/// (`ml::cv::CrossValidatedScorer::new_from_shuffled_with_precomputed`) so that
/// a row's `lda_score` never saw its own label in either the LDA fit or the
/// GBM fold it is held out of; keeping the partition in one function is what
/// stops the two sides from drifting.
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

    // Build BOTH lane frames AFTER the shuffle over the SAME `data`, so row `i`
    // is the same candidate in lin, nl, lda_score, responses, and the moved
    // `data`.
    let lin = build_lane_frame(&data, Lane::Linear);
    let lin_ncols = lin.ncols();
    let lin_rm = lin.row_major();
    let is_decoy: Vec<bool> = data.iter().map(|c| c.get_y() < 0.5).collect();

    // --- CROSS-FIT lda_score (leak-free), via the shared partition ---
    // On failure the column is left uniformly zero rather than partially
    // filled: a fold-dependent `lda_score` is a feature the GBM can split on to
    // recover fold identity, which is strictly worse than no feature at all.
    // An all-zero column is constant -> zero split gain -> hybrid degrades to
    // "GBM on the nonlinear lane", which is a sane, loud degradation.
    let lda_score = match crossfit_lda(&lin_rm, nrows, lin_ncols, &is_decoy) {
        Some(cf) => cf.scores,
        None => {
            tracing::error!(
                "hybrid: cross-fit LDA failed; lda_score is uniformly 0 (GBM falls back to the \
                 nonlinear lane alone)"
            );
            vec![0.0f64; nrows]
        }
    };

    // --- Build the NONLINEAR frame + push lda_score as the LAST column ---
    // `nl.push_column` appends lda_score after the nonlinear columns, so `names`
    // must be the nonlinear names THEN "lda_score" to match `row_major` order.
    let mut nl = build_lane_frame(&data, Lane::Nonlinear);
    let mut names = nonlinear_feature_name_set_for(&data);
    nl.push_column("lda_score", lda_score);
    names.push(Arc::from("lda_score"));
    debug_assert_eq!(
        nl.names(),
        names.as_slice(),
        "hybrid frame columns must match nonlinear names ++ lda_score"
    );

    let ncols = nl.ncols();
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let precomputed = PrecomputedFeatures::from_row_major(nl.row_major(), ncols, responses);

    let mut scorer = CrossValidatedScorer::<CompetedCandidate>::new_from_shuffled_with_precomputed(
        N_RESCORE_FOLDS,
        data,
        config,
        precomputed,
    );
    scorer
        .fit(&mut DataBuffer::default(), &mut DataBuffer::default())
        .unwrap();

    let stats = scorer.feature_stats(&names);

    let mut scored = scorer.score();
    // Sort by score descending — identical tail to `rescore`.
    #[cfg(feature = "rayon")]
    scored.par_sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    #[cfg(not(feature = "rayon"))]
    scored.sort_unstable_by(|a, b| b.get_score().total_cmp(&a.get_score()));
    assign_qval(&mut scored, |x| CompetedCandidate::get_score(x) as f32);
    debug!("Best:\n{:#?}", scored.first());
    debug!("Worst:\n{:#?}", scored.last());

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

/// Per-fold feature stats for the cross-fit LDA path, mirroring the GBM's
/// [`CrossValidatedScorer::feature_stats`] semantics exactly: for fold `f` the
/// means / NaN ratios are computed over that fold's HELD-OUT rows
/// (`i % N_RESCORE_FOLDS == f`, the rows `models[f]` actually scored), and the
/// importance ranking is `|coef|` of that fold's model (LDA weights in
/// standardized space are directly interpretable as importance).
///
/// `models` is indexed by fold and must be `N_RESCORE_FOLDS` long — the
/// invariant `crossfit_lda` upholds.
fn lda_feature_stats(
    names: &[Arc<str>],
    feat: &[f64],
    nrows: usize,
    ncols: usize,
    models: &[LdaModel],
) -> RescoreFeatureStats {
    debug_assert_eq!(names.len(), ncols);
    let n_folds = N_RESCORE_FOLDS as usize;
    let mut out: RescoreFeatureStats = Vec::with_capacity(n_folds);

    for (f, model) in models.iter().enumerate() {
        let mut sums = vec![0.0f64; ncols];
        let mut finite = vec![0u32; ncols];
        let mut nan = vec![0u32; ncols];
        let mut n = 0usize;
        for i in (f..nrows).step_by(n_folds) {
            let row = &feat[i * ncols..(i + 1) * ncols];
            for j in 0..ncols {
                let v = row[j];
                if v.is_finite() {
                    sums[j] += v;
                    finite[j] += 1;
                } else {
                    nan[j] += 1;
                }
            }
            n += 1;
        }
        let feature_stats: Vec<FeatureStat> = if n > 0 {
            (0..ncols)
                .map(|j| FeatureStat {
                    name: names[j].clone(),
                    mean: if finite[j] > 0 {
                        (sums[j] / finite[j] as f64) as f32
                    } else {
                        f32::NAN
                    },
                    nan_ratio: nan[j] as f32 / n as f32,
                })
                .collect()
        } else {
            Vec::new()
        };

        let mut feature_importance: Vec<(Arc<str>, f32)> = names
            .iter()
            .zip(model.coef().iter())
            .map(|(nm, c)| (nm.clone(), c.abs() as f32))
            .collect();
        feature_importance
            .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        out.push(FoldStats {
            fold: f as u8,
            feature_stats,
            feature_importance,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// Lane feature walks (FeatFrame) — the ML consumer's live path
// ---------------------------------------------------------------------------

/// Which disjoint feature lane a walk targets. `Linear` = LDA (fields
/// approx-Gaussian after their declared transform); `Nonlinear` = the rest
/// (context, counts, sequence). `all` = linear ++ nonlinear (GBM).
#[derive(Clone, Copy)]
enum Lane {
    Linear,
    Nonlinear,
}

/// Build a batch [`FeatFrame`] for `data` in its CURRENT order (call AFTER any
/// shuffle so frame row `i` aligns with `data[i]`). The walk order per lane is
/// fixed — scoring blocks (composition order) -> `ResultMeta` -> `Derived` ->
/// (nonlinear only) `sequence_counts` — identical to the name-set walk below,
/// so frame columns and names cannot desync.
fn build_lane_frame(data: &[CompetedCandidate], lane: Lane) -> FeatFrame {
    let mut frame = FeatFrame::with_capacity(0, data.len());
    {
        let mut s = FrameSink::new(&mut frame, data.len());
        for c in data {
            s.begin_row();
            match lane {
                Lane::Linear => {
                    c.scoring.push_linear_features(&mut s);
                    c.result_meta().linear_features(&mut s);
                    Derived::compute(&c.scoring).linear_features(&mut s);
                    // sequence_counts has NO linear lane (context features).
                }
                Lane::Nonlinear => {
                    c.scoring.push_nonlinear_features(&mut s);
                    c.result_meta().nonlinear_features(&mut s);
                    Derived::compute(&c.scoring).nonlinear_features(&mut s);
                    sequence_counts::nonlinear_features(&c.scoring.identity.peptide, &mut s);
                }
            }
        }
    }
    frame
}

/// The ALL-lane frame (linear then nonlinear, per row) — the GBM feature set.
/// Column order = linear columns ++ nonlinear columns, matching
/// [`all_feature_name_set`].
fn build_all_frame(data: &[CompetedCandidate]) -> FeatFrame {
    let mut frame = build_lane_frame(data, Lane::Linear);
    frame.extend(build_lane_frame(data, Lane::Nonlinear));
    frame
}

/// LINEAR-lane feature names (LDA), in the `build_lane_frame` walk order.
///
/// Takes NO gate, unlike the other two lanes: `sequence_counts` is the only
/// gated block and it has no linear lane, so this set is gate-invariant. If a
/// gated LINEAR feature is ever added, this signature has to grow the parameter
/// back — the absence of it is the reminder.
pub fn linear_feature_name_set() -> Vec<Arc<str>> {
    let mut n = NameSink::new();
    ScoringFields::push_linear_feature_names(&mut n);
    <ResultMeta as ScoreBlock>::linear_feature_names(&mut n);
    <Derived as ScoreBlock>::linear_feature_names(&mut n);
    n.into_names()
}

/// NONLINEAR-lane feature names, in the `build_lane_frame` walk order. The
/// gate appends the 22 trailing `sequence_counts` names when the run has parsed
/// sequences (speclib-wide).
pub fn nonlinear_feature_name_set(gate_on: bool) -> Vec<Arc<str>> {
    let mut n = NameSink::new();
    ScoringFields::push_nonlinear_feature_names(&mut n);
    <ResultMeta as ScoreBlock>::nonlinear_feature_names(&mut n);
    <Derived as ScoreBlock>::nonlinear_feature_names(&mut n);
    if gate_on {
        sequence_counts::nonlinear_feature_names(&mut n);
    }
    n.into_names()
}

/// The ALL-lane feature names (GBM) = linear ++ nonlinear, matching
/// `build_all_frame`'s column order.
pub fn all_feature_name_set(gate_on: bool) -> Vec<Arc<str>> {
    // Only the nonlinear half is gate-sensitive.
    let mut v = linear_feature_name_set();
    v.extend(nonlinear_feature_name_set(gate_on));
    v
}

/// Read the speclib-wide sequence gate off the first candidate.
fn gate_for(candidates: &[CompetedCandidate]) -> bool {
    candidates
        .first()
        .map(|c| c.scoring.identity.peptide.aa_counts().is_some())
        .unwrap_or(false)
}

/// [`linear_feature_name_set`], kept as a `*_for` for call-site symmetry with
/// the gated lanes (the linear lane itself has nothing to gate on).
pub fn linear_feature_name_set_for(_candidates: &[CompetedCandidate]) -> Vec<Arc<str>> {
    linear_feature_name_set()
}

/// [`nonlinear_feature_name_set`] with the gate read from `candidates`.
pub fn nonlinear_feature_name_set_for(candidates: &[CompetedCandidate]) -> Vec<Arc<str>> {
    nonlinear_feature_name_set(gate_for(candidates))
}

/// [`all_feature_name_set`] with the gate read from `candidates`.
pub fn all_feature_name_set_for(candidates: &[CompetedCandidate]) -> Vec<Arc<str>> {
    all_feature_name_set(gate_for(candidates))
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

    #[test]
    fn lane_frame_columns_match_name_sets() {
        // LANE ORDER PARITY: the FeatFrame column order MUST equal the name-set
        // order for each lane. A desync here silently misattributes feature
        // importances / stats and (for the value/label join) corrupts nothing
        // by itself, but a column/name drift is exactly the class of bug this
        // asserts against. Gate on (parsed sequence).
        let data = vec![sample_competed_candidate_parsed()];

        let lin = build_lane_frame(&data, Lane::Linear);
        assert_eq!(lin.names(), linear_feature_name_set().as_slice());

        let nonlin = build_lane_frame(&data, Lane::Nonlinear);
        assert_eq!(nonlin.names(), nonlinear_feature_name_set(true).as_slice());

        let all = build_all_frame(&data);
        assert_eq!(all.names(), all_feature_name_set(true).as_slice());
        assert_eq!(all.ncols(), lin.ncols() + nonlin.ncols());
    }

    #[test]
    fn lane_frame_columns_match_name_sets_gate_off() {
        // With no parsed sequence the nonlinear lane drops the 22 sequence names.
        let data = vec![sample_competed_candidate_unparsed()];
        let nonlin = build_lane_frame(&data, Lane::Nonlinear);
        assert_eq!(nonlin.names(), nonlinear_feature_name_set(false).as_slice());
        let all = build_all_frame(&data);
        assert_eq!(all.names(), all_feature_name_set(false).as_slice());
    }

    #[test]
    fn lane_feature_sets_are_nontrivial() {
        // No exact counts here on purpose: adding or removing a score must not
        // break this test. What it does pin is that no lane collapsed to a
        // no-op, and that the full set has no duplicate names.
        let lin = linear_feature_name_set().len();
        let nonlin = nonlinear_feature_name_set(true).len();
        let all = all_feature_name_set(true).len();
        assert!(lin > 0, "linear lane collapsed");
        assert!(nonlin > 0, "nonlinear lane collapsed");
        assert!(all > 0, "feature set collapsed");
        // No duplicate names across the full set.
        let mut seen = std::collections::HashSet::new();
        for n in all_feature_name_set(true) {
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

        let cf = crossfit_lda(&feat, N, 3, &is_decoy).expect("cross-fit must succeed");
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

    #[test]
    fn rescore_hybrid_smoke_and_determinism() {
        let n = 90;
        let (out_a, _stats_a) = rescore_hybrid(synthetic_competed(n));
        let (out_b, _stats_b) = rescore_hybrid(synthetic_competed(n));

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

    #[test]
    fn gate_delta_is_the_sequence_count_block() {
        // Structural, not an arbitrary total: the sequence-count gate adds one
        // column per canonical amino acid plus `peptide_length` and
        // `peptide_n_mods`. Derived from the source of truth so adding an
        // amino-acid column updates both sides at once.
        let expected = crate::models::AA_COUNT_NAMES.len() + 2;
        let on = nonlinear_feature_name_set(true).len();
        let off = nonlinear_feature_name_set(false).len();
        assert_eq!(on - off, expected);
    }

    #[test]
    fn neutralize_mobility_nans_every_mobility_feature() {
        // For a run with no searchable mobility axis (mzML/FAIMS), every
        // mobility-derived GBM feature must become NaN (forust missing), so it
        // cannot bias the score with sentinel-derived constants. Walked over the
        // LIVE lane path (`build_all_frame` == linear ++ nonlinear), not a
        // hand-listed field set, so a new mobility field that forgets to
        // neutralize fails here.
        //
        // Two things are pinned, and both matter:
        //   (a) every mobility feature is NaN (or, for the `_isna` companions,
        //       flipped to 1.0 == "missing" — an isna column is BY CONSTRUCTION
        //       never NaN, so demanding NaN there would be wrong);
        //   (b) every non-mobility feature is bit-for-bit unchanged. Without (b)
        //       an impl that NaN'd the whole record would pass (a).
        let before = build_all_frame(&[sample_competed_candidate_parsed()]);

        let mut cand = sample_competed_candidate_parsed();
        cand.scoring.neutralize_mobility();
        let after = build_all_frame(&[cand]);

        let names = before.names().to_vec();
        assert_eq!(names.as_slice(), after.names());
        assert_eq!(names, all_feature_name_set(true));

        let is_mobility = |n: &str| n.contains("mob");
        let mut nan_mob: Vec<&str> = Vec::new();
        let mut isna_mob: Vec<&str> = Vec::new();
        let mut finite_non_mob = 0usize;

        for (j, name) in names.iter().enumerate() {
            let (b, a) = (before.column(j)[0], after.column(j)[0]);
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
}
