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
    NameSink,
    ScoreBlock,
    sequence_counts,
};
use crate::scoring::results::{
    CompetedCandidate,
    FeatureRow,
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

pub fn report_qvalues_at_thresholds<T: LabelledScore + std::fmt::Debug>(
    scores: &[T],
    thresholds: &[f32],
) -> Vec<(f32, usize, usize, usize)> {
    // One pass for all thresholds. Three filter-counts per threshold walked
    // the whole slice 3N times for N cutoffs and re-derived `n_below` from a
    // third scan instead of the two class counts.
    let mut counts = vec![(0usize, 0usize); thresholds.len()];
    for s in scores {
        let q = s.get_qval();
        let is_target = matches!(s.get_label(), TargetDecoy::Target);
        for (&thresh, c) in thresholds.iter().zip(counts.iter_mut()) {
            if q <= thresh {
                if is_target {
                    c.0 += 1;
                } else {
                    c.1 += 1;
                }
            }
        }
    }
    thresholds
        .iter()
        .zip(counts)
        .map(|(&thresh, (n_targets, n_decoys))| (thresh, n_targets + n_decoys, n_targets, n_decoys))
        .collect()
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
/// THE canonical statement of leak-freedom — `rescore_lda` and `rescore_hybrid`
/// both route through here and their docs just point back at this one.
///
/// # Why a row may never be scored by a model that saw it
/// An LDA is label-aware: it fits on the target/decoy labels. One in-sample fit
/// over all rows, scoring those same rows, lets every row's discriminant peek at
/// its own label; the target/decoy separation then looks better than it is, and
/// `assign_qval` derives the q-values from exactly that separation. The result
/// is an FDR that is wrong in the flattering direction with nothing downstream
/// to catch it. `N_RESCORE_FOLDS` fits instead of 1 is affordable precisely
/// because the LDA is closed-form and cheap.
///
/// `rescore_hybrid` needs this even more sharply: there the held-out score is
/// not the final score but ONE COLUMN (`lda_score`) fed to a cross-validated
/// GBM. An in-sample `lda_score` would smuggle a row's own label into a feature
/// the GBM reads while that row is held out of its GBM fold — so the GBM's own
/// CV cannot notice, and the leak surfaces only as an optimistic FDR.
///
/// # THE partition
/// `fold(i) = i % N_RESCORE_FOLDS`, the single definition used by both callers.
/// For fold `f` the model is fitted on every row `i` with
/// `i % N_RESCORE_FOLDS != f` and then scores only the rows with
/// `i % N_RESCORE_FOLDS == f`, so no row contributes to the model scoring it.
///
/// It MUST equal the partition `CrossValidatedScorer` builds internally
/// (`assigned_fold[i] = i % n_folds`, see
/// `ml::cv::CrossValidatedScorer::new_from_shuffled`), or a hybrid row's
/// `lda_score` can come from an LDA trained on rows the GBM is holding out —
/// leak restored, silently. NOTHING ENFORCES THAT MATCH: it is the same modulo
/// written in two modules, kept together only by each side having exactly one
/// definition. Changing either means changing both.
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

    // Build the ALL-lane matrix (linear ++ nonlinear) over the shuffled `data`,
    // so row `i` aligns with `data[i]` — fold assignment + labels are
    // positional, so the matrix MUST be built AFTER the shuffle. The GBM sees
    // the full lane feature set; feature names come from the same walk order,
    // so columns and `names` align by construction (width asserted below + the
    // lane-parity test).
    let names = all_feature_name_set();
    debug_assert_eq!(names.len(), ALL_NCOLS);
    let feat = build_all_matrix(&data);
    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let precomputed = PrecomputedFeatures::from_row_major(feat, ALL_NCOLS, responses);

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

    finalize(scorer.score(), stats)
}

/// Sage-style shrinkage-LDA rescorer: closed-form linear fits, no boosting —
/// ~100x cheaper than the GBM path. The FDR machinery (`assign_qval`,
/// target-decoy competition) is untouched — only the discriminant score source
/// changes.
///
/// CROSS-FIT, not a single in-sample fit: every row's score comes from an LDA
/// fitted without that row, via [`crossfit_lda`] — see there for the partition
/// and why it is mandatory.
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
            // All-or-nothing (see `crossfit_lda`). Here that means every score
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

/// Hybrid rescorer: cross-fit an LDA on the LINEAR lane, push its (leak-free)
/// `lda_score` as one extra column into the NONLINEAR lane, then train the GBM
/// CV on `nonlinear + lda_score`. The GBM re-sees ~30 features instead of the
/// full ~131 (the compression play) at ~parity.
///
/// LEAK-FREEDOM: `lda_score` is cross-fit via [`crossfit_lda`] — see there for
/// the partition, why a label-aware feature fed to a CV'd GBM in particular
/// must be leak-free, and why that partition has to match the one
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
    let lin = build_linear_matrix(&data);
    let is_decoy: Vec<bool> = data.iter().map(|c| c.get_y() < 0.5).collect();

    // --- CROSS-FIT lda_score (leak-free), via the shared partition ---
    // On failure the column is left uniformly zero rather than partially
    // filled: a fold-dependent `lda_score` is a feature the GBM can split on to
    // recover fold identity, which is strictly worse than no feature at all.
    // An all-zero column is constant -> zero split gain -> hybrid degrades to
    // "GBM on the nonlinear lane", which is a sane, loud degradation.
    let lda_score = match crossfit_lda(&lin, nrows, LINEAR_NCOLS, &is_decoy) {
        Some(cf) => cf.scores,
        None => {
            tracing::error!(
                "hybrid: cross-fit LDA failed; lda_score is uniformly 0 (GBM falls back to the \
                 nonlinear lane alone)"
            );
            vec![0.0f64; nrows]
        }
    };

    // --- NONLINEAR matrix with lda_score appended as the LAST column ---
    // Row-major, so "extra trailing column" is literally "one more push at the
    // end of each row"; `names` is the nonlinear names THEN "lda_score" to
    // match.
    let ncols = NONLINEAR_NCOLS + 1;
    let nl = build_nonlinear_matrix(&data);
    let mut feat = Vec::with_capacity(nrows * ncols);
    for (row, s) in nl.chunks_exact(NONLINEAR_NCOLS).zip(lda_score) {
        feat.extend_from_slice(row);
        feat.push(s);
    }
    let mut names = nonlinear_feature_name_set();
    names.push(Arc::from("lda_score"));
    debug_assert_eq!(names.len(), ncols);
    debug_assert_eq!(feat.len(), nrows * ncols);

    let responses: Vec<f64> = data.iter().map(|c| c.get_y()).collect();
    let precomputed = PrecomputedFeatures::from_row_major(feat, ncols, responses);

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

    finalize(scorer.score(), stats)
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
fn build_linear_matrix<R: FeatureRow>(data: &[R]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * LINEAR_NCOLS);
    for r in data {
        let (s, meta) = (r.scoring(), r.result_meta());
        push_linear_row(s, &meta, &Derived::compute(s), &mut out);
    }
    out
}

/// The NONLINEAR-lane matrix for `data`, `NONLINEAR_NCOLS` wide (see
/// [`build_linear_matrix`] for the ordering contract).
fn build_nonlinear_matrix<R: FeatureRow>(data: &[R]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * NONLINEAR_NCOLS);
    for r in data {
        let (s, meta) = (r.scoring(), r.result_meta());
        push_nonlinear_row(s, &meta, &Derived::compute(s), &mut out);
    }
    out
}

/// The ALL-lane matrix (linear then nonlinear, per row) — the GBM feature set,
/// `ALL_NCOLS` wide, matching [`all_feature_name_set`]'s order.
///
/// ONE pass over `data` with ONE `Derived::compute` per row: the two lanes are
/// adjacent within a row, so there is nothing to gain from walking twice.
///
/// Generic over [`FeatureRow`] so the matrix trained on `CompetedCandidate`s
/// and the one the dashboard displays for `FinalResult`s are the same code.
fn build_all_matrix<R: FeatureRow>(data: &[R]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len() * ALL_NCOLS);
    for r in data {
        let (s, meta) = (r.scoring(), r.result_meta());
        let derived = Derived::compute(s);
        push_linear_row(s, &meta, &derived, &mut out);
        push_nonlinear_row(s, &meta, &derived, &mut out);
    }
    out
}

/// The ALL-lane feature names + matrix for post-rescore rows: the dashboard's
/// entry point.
///
/// Row-major: value `j` of row `i` is at `matrix[i * names.len() + j]`.
pub fn feature_frame(data: &[FinalResult]) -> (Vec<Arc<str>>, Vec<f64>) {
    (all_feature_name_set(), build_all_matrix(data))
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
            assert_eq!(build_all_matrix(&data).len(), ALL_NCOLS);
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
        let all = build_all_matrix(&data);

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
        let before = build_all_matrix(&[sample_competed_candidate_parsed()]);

        let mut cand = sample_competed_candidate_parsed();
        cand.scoring.neutralize_mobility();
        let after = build_all_matrix(&[cand]);

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

    /// `build_all_matrix` over `CompetedCandidate` and over `FinalResult` must
    /// produce byte-identical rows. They share the row builders, so what this
    /// guards is `into_final`: if that conversion ever drops, rounds or
    /// reorders a value the feature row reads, the dashboard would display a
    /// matrix the model never trained on, silently.
    ///
    /// The candidate is the PARSED one on purpose. `sample_default()` has no
    /// sequence, so all 22 `sequence_counts` columns come out NaN and the
    /// NaN-equals-NaN escape below passes them without comparing anything.
    #[test]
    fn into_final_preserves_every_feature_value() {
        let competed = sample_competed_candidate_parsed();
        let expected = build_all_matrix(std::slice::from_ref(&competed));

        let final_result = competed.clone().into_final();
        let (names, got) = feature_frame(std::slice::from_ref(&final_result));

        assert_eq!(names, all_feature_name_set());
        assert_eq!(got.len(), names.len(), "one row, one value per name");
        let compared = expected
            .iter()
            .zip(&got)
            .filter(|(a, b)| !a.is_nan() && !b.is_nan())
            .count();
        assert!(
            compared > 22,
            "only {compared} of {} columns were non-NaN; the fixture is not \
             exercising the sequence_counts block",
            names.len()
        );
        for (j, (a, b)) in expected.iter().zip(got.iter()).enumerate() {
            assert!(
                a == b || (a.is_nan() && b.is_nan()),
                "column {j} ({}) differs: competed={a} final={b}",
                names[j]
            );
        }
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
