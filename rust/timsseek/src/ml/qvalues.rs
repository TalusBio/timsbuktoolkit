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
    data.sort_unstable_by_key(|c| (c.scoring.library_id, c.scoring.precursor_charge));

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
    data.sort_unstable_by_key(|c| (c.scoring.library_id, c.scoring.precursor_charge));

    let base_names: Vec<&'static str> = data.first().map(|c| c.feature_names()).unwrap_or_default();
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
    names: &[&'static str],
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
        std::fs::write(format!("{prefix}.names.txt"), names.join("\n"))?;
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
fn absolutize_center_better(base: &mut [f64], base_ncols: usize, base_names: &[&'static str]) {
    let abs_cols: Vec<usize> = (0..base_ncols)
        .filter(|&j| is_center_better(base_names[j]))
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
///
/// Indicator names are interned via `Box::leak` (a handful of small strings,
/// once per run) to satisfy the `&'static str` name contract shared with the
/// GBM feature-stats path.
fn augment_with_missingness(
    base: &[f64],
    nrows: usize,
    base_names: &[&'static str],
) -> (Vec<f64>, Vec<&'static str>) {
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
    let mut names: Vec<&'static str> = base_names.to_vec();
    for &j in &sel {
        let nm: &'static str = Box::leak(format!("{}_isna", base_names[j]).into_boxed_str());
        names.push(nm);
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
    names: &[&'static str],
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
            name: names[j],
            mean: if finite[j] > 0 {
                (sums[j] / finite[j] as f64) as f32
            } else {
                f32::NAN
            },
            nan_ratio: nan[j] as f32 / n,
        })
        .collect();

    let mut feature_importance: Vec<(&'static str, f32)> = names
        .iter()
        .zip(coef.iter())
        .map(|(nm, c)| (*nm, c.abs() as f32))
        .collect();
    feature_importance.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    vec![FoldStats {
        fold: 0,
        feature_stats,
        feature_importance,
    }]
}

use crate::models::AA_COUNT_NAMES;
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
            off, 86,
            "base (gate-off) feature count is locked at 86; update this test if the set intentionally changes"
        );
        let on = sample_competed_candidate_parsed().feature_names().len();
        assert_eq!(on, 108, "gate-on total is 86 base + 22 sequence = 108");
    }
}
