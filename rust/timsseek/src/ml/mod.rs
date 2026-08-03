pub mod cv;
pub mod lda;
mod mlp;
mod mlp_fold;
pub mod qvalues;
pub use cv::RescoreFeatureStats;
pub use qvalues::{
    RescoreError,
    RescoreResult,
    rescore,
    rescore_hybrid,
    rescore_lda,
    rescore_mlp,
};

use crate::scoring::results::CompetedCandidate;
use serde::{
    Deserialize,
    Serialize,
};

/// Fold count shared by every rescorer.
///
/// Must be >= 3: `CrossValidatedScorer::get_scores` averages over the
/// `n_folds - 2` folds that are neither the training nor the early-stopping
/// fold, which degenerates (division by zero) at 2.
pub const N_RESCORE_FOLDS: u8 = 3;

/// Which rescorer produces the final discriminant score.
///
/// Selected via the `rescore_model` config field / the `--rescore-model` CLI
/// flag (CLI wins). Dispatch through [`rescore_with`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum RescoreModel {
    /// Gradient-boosted trees over all linear and nonlinear features.
    Gbm,
    /// Sage-style shrinkage LDA on the linear features only.
    Lda,
    /// Cross-fit LDA followed by GBM on nonlinear features plus the LDA score.
    Hybrid,
    /// MLP over all linear and nonlinear features.
    #[default]
    Mlp,
}

/// Dispatch to the rescorer named by `model`.
///
/// Tests pin every arm against its direct entry point.
pub fn rescore_with(model: RescoreModel, data: Vec<CompetedCandidate>) -> RescoreResult {
    match model {
        RescoreModel::Gbm => rescore(data),
        RescoreModel::Lda => rescore_lda(data),
        RescoreModel::Hybrid => rescore_hybrid(data),
        RescoreModel::Mlp => rescore_mlp(data),
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
