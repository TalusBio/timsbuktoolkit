//! Adapter from a finished rescoring run to the `rescore_dash` TUI.
//!
//! Nothing here is pipeline logic; it exists so `processing` carries two
//! `#[cfg(feature = "dashboard")]` call sites and no feature plumbing.

use std::collections::HashMap;
use std::sync::Arc;
use timsseek::ml::RescoreFeatureStats;
use timsseek::scoring::results::FinalResult;
use timsseek::scoring::timings::TimedStep;
use tracing::info;

/// Dev-only, gated on an environment variable rather than a CLI flag: a
/// debugging tool, not part of the documented interface.
///
/// Any value except an explicit `0`/`false`. Matching only `"1"` made
/// `TIMSSEEK_RESCORE_DASHBOARD=true` a silent no-op.
fn enabled() -> bool {
    std::env::var("TIMSSEEK_RESCORE_DASHBOARD")
        .is_ok_and(|v| !matches!(v.trim().to_ascii_lowercase().as_str(), "" | "0" | "false"))
}

/// Fold-averaged GBM gain, aligned to `feature_names`.
///
/// A feature no fold reported gets `0.0`; the LDA path reports coefficients
/// for the linear lane only, so that is a normal outcome, not a missing
/// value. Averaged here so the dashboard indexes gain by column instead of
/// searching the per-fold lists by name.
fn mean_gain_per_feature(stats: &RescoreFeatureStats, feature_names: &[Arc<str>]) -> Vec<f32> {
    let mut sum = vec![0.0f32; feature_names.len()];
    let mut n = vec![0u32; feature_names.len()];
    let index: HashMap<&str, usize> = feature_names
        .iter()
        .enumerate()
        .map(|(j, name)| (&**name, j))
        .collect();
    for fold in stats {
        for (name, gain) in &fold.feature_importance {
            if let Some(&j) = index.get(&**name) {
                sum[j] += gain;
                n[j] += 1;
            }
        }
    }
    sum.iter()
        .zip(&n)
        .map(|(&s, &k)| if k == 0 { 0.0 } else { s / k as f32 })
        .collect()
}

/// Sweep `data` into a `Dashboard`, or `None` if the dashboard is off or the
/// view is rejected.
///
/// The ~1 GB feature matrix is built and dropped inside this function, so
/// nothing matrix-sized outlives it. `enabled()` is checked first, so a run
/// without the env var pays nothing; whether stdout is a terminal is
/// `rescore_dash::run`'s call, so that it is the one place that warns.
pub fn build(
    data: &[FinalResult],
    feature_stats: &RescoreFeatureStats,
    qval_report: &[(f32, usize, usize, usize)],
) -> Option<rescore_dash::Dashboard> {
    if !enabled() {
        return None;
    }
    let step = TimedStep::begin("Rescore dashboard precompute");
    let (feature_names, matrix) = timsseek::ml::qvalues::feature_frame(data);
    let is_target: Vec<bool> = data.iter().map(|r| r.scoring.identity.is_target).collect();
    let score: Vec<f32> = data.iter().map(|r| r.discriminant_score).collect();
    let qvalue: Vec<f32> = data.iter().map(|r| r.qvalue).collect();
    let gain = mean_gain_per_feature(feature_stats, &feature_names);
    // The same counts already in the run log, so the panel and the log cannot
    // disagree and the dashboard does not re-scan for them.
    let thresholds: Vec<rescore_dash::ThresholdRow> = qval_report
        .iter()
        .map(|&(q, _total, targets, decoys)| rescore_dash::ThresholdRow { q, targets, decoys })
        .collect();
    let matrix_mb = (matrix.len() * size_of::<f64>()) as f64 / (1024.0 * 1024.0);
    info!(
        "rescore dashboard: {} rows x {} features ({matrix_mb:.1} MB matrix, freed after precompute)",
        is_target.len(),
        feature_names.len()
    );
    let built = rescore_dash::Dashboard::build(
        &rescore_dash::RescoreView {
            feature_names: &feature_names,
            features: &matrix,
            is_target: &is_target,
            score: &score,
            qvalue: &qvalue,
            thresholds: &thresholds,
            gain: &gain,
        },
        rescore_dash::DEFAULT_SAMPLE,
    );
    step.finish();
    built
        .inspect_err(|e| tracing::warn!("rescore dashboard input rejected: {e}"))
        .ok()
}
