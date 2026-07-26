pub mod cv;
pub mod lda;
pub mod qvalues;
pub use cv::RescoreFeatureStats;
pub use qvalues::{
    rescore,
    rescore_hybrid,
    rescore_lda,
};

use crate::scoring::results::{
    CompetedCandidate,
    FinalResult,
};
use serde::{
    Deserialize,
    Serialize,
};

/// Fold count for every rescorer (`rescore`, `rescore_lda`, `rescore_hybrid`).
///
/// Shared rather than repeated per call site because the LDA cross-fit
/// (`qvalues::crossfit_lda`, used by both the LDA and hybrid paths) requires its
/// partition (`i % N_RESCORE_FOLDS`) to be the SAME partition
/// `CrossValidatedScorer` builds internally — if the two drift, a row's
/// `lda_score` can peek at its own label.
///
/// Must be >= 3: `CrossValidatedScorer::get_scores` averages over the
/// `n_folds - 2` folds that are neither the training nor the early-stopping
/// fold, which degenerates (division by zero) at 2.
pub const N_RESCORE_FOLDS: u8 = 3;

/// Which rescorer produces the final discriminant score.
///
/// Lives here, next to the three rescorers it selects between, so that
/// integration tests, benches and the pyo3 bindings can pick a model without
/// going through the CLI. The CLI keeps a `clap::ValueEnum` mirror of this enum
/// (`timsseek_cli::cli::CliRescoreModel`) purely because `ValueEnum` is a
/// foreign trait it cannot implement for a lib-owned type; that mirror converts
/// into this via an exhaustive `From`, so adding a variant here is a compile
/// error there rather than a silent drift.
///
/// Selected via the `rescore_model` config field / the `--rescore-model` CLI
/// flag (CLI wins). Dispatch through [`rescore_with`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum RescoreModel {
    /// Gradient-boosted trees over the ALL lane (linear ++ nonlinear, ~131
    /// features), cross-validated over `N_RESCORE_FOLDS`. The default.
    #[default]
    Gbm,
    /// Sage-style shrinkage LDA on the LINEAR lane only (see [`lda`]). ~100x
    /// cheaper than the GBM; cross-fit over the same fold partition.
    Lda,
    /// Cross-fit LDA (linear lane) -> `lda_score` column -> GBM CV on
    /// `nonlinear + lda_score`: the same fold partition on both sides, so a
    /// row's `lda_score` never saw its own label.
    Hybrid,
}

/// Dispatch to the rescorer named by `model`.
///
/// All three rescorers share this signature (candidates in, scored+q-valued
/// results plus per-fold feature stats out), so this is the single place model
/// selection is resolved.
pub fn rescore_with(
    model: RescoreModel,
    data: Vec<CompetedCandidate>,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    match model {
        RescoreModel::Gbm => rescore(data),
        RescoreModel::Lda => rescore_lda(data),
        RescoreModel::Hybrid => rescore_hybrid(data),
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TargetDecoy {
    Target,
    Decoy,
}

pub trait LabelledScore {
    fn get_label(&self) -> TargetDecoy;
    fn assign_qval(&mut self, q: f32);
    fn get_qval(&self) -> f32;
}

impl LabelledScore for (f64, TargetDecoy, f32) {
    fn get_label(&self) -> TargetDecoy {
        self.1
    }

    fn assign_qval(&mut self, q: f32) {
        self.2 = q
    }

    fn get_qval(&self) -> f32 {
        self.2
    }
}
