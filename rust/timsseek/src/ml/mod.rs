pub mod cv;
pub mod lda;
pub mod mlp;
pub mod mlp_fold;
pub mod qvalues;
pub use cv::RescoreFeatureStats;
pub use qvalues::{
    rescore,
    rescore_hybrid,
    rescore_hybrid_mlp,
    rescore_lda,
    rescore_mlp_all,
    rescore_mlp_linear,
};

use crate::scoring::results::{
    CompetedCandidate,
    FinalResult,
};
use serde::{
    Deserialize,
    Serialize,
};

/// Fold count for every rescorer (`rescore`, `rescore_lda`, `rescore_hybrid`,
/// `rescore_mlp_linear`, `rescore_mlp_all`, `rescore_hybrid_mlp`).
///
/// Shared rather than repeated per call site because the cross-fit
/// (`qvalues::crossfit`, used by the LDA path and by both hybrids' extra column)
/// requires the same fold ASSIGNMENT that `CrossValidatedScorer` uses, or a
/// hybrid row's `lda_score` / `mlp_score` can peek at its own label. Both sides
/// read that assignment from `RowMajorDataset::get_fold`, so it has one
/// definition; the fold COUNT living here is what keeps it fed the same argument.
///
/// The two do NOT use the same train/score split on top of that assignment, and
/// must not be made to: `CrossValidatedScorer` fits on one fold and scores the
/// rest, `crossfit` fits on all but one fold and scores that one. Both are
/// leak-free; see `qvalues::crossfit` for why the difference is deliberate.
///
/// Must be >= 3: `CrossValidatedScorer::get_scores` averages over the
/// `n_folds - 2` folds that are neither the training nor the early-stopping
/// fold, which degenerates (division by zero) at 2.
pub const N_RESCORE_FOLDS: u8 = 3;

/// Which rescorer produces the final discriminant score.
///
/// Lives here, next to the six rescorers it selects between, so that
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
#[serde(rename_all = "snake_case")]
pub enum RescoreModel {
    /// Gradient-boosted trees over the ALL lane (linear ++ nonlinear, 128
    /// features), cross-validated over `N_RESCORE_FOLDS`. The default.
    #[default]
    Gbm,
    /// Sage-style shrinkage LDA on the LINEAR lane only (101 features, see
    /// [`lda`]). ~20x cheaper than the GBM in Phase 5 wall time on a
    /// 114138-candidate mzML run — see [`rescore_lda`] for what that measurement
    /// does and does not cover. Cross-fit over the same fold ASSIGNMENT.
    Lda,
    /// Cross-fit LDA (linear lane) -> `lda_score` column -> GBM CV on
    /// `nonlinear + lda_score`: the same fold ASSIGNMENT on both sides, so a
    /// row's `lda_score` never saw its own label.
    Hybrid,
    /// MLP on the LINEAR lane only (101 features) — the same feature set
    /// [`Lda`](RescoreModel::Lda) trains on, so it isolates "nonlinear model"
    /// from "more features". THE DEFAULT MLP SHAPE, hence the unqualified name.
    Mlp,
    /// MLP on the ALL lane (linear ++ nonlinear, 128 features) — the same
    /// feature set [`Gbm`](RescoreModel::Gbm) trains on, so the two are
    /// directly comparable.
    MlpAll,
    /// MLP counterpart of [`Hybrid`](RescoreModel::Hybrid): cross-fit MLP
    /// (linear lane) -> `mlp_score` column -> GBM CV on
    /// `nonlinear + mlp_score`, same fold ASSIGNMENT on both sides.
    HybridMlp,
}

/// Dispatch to the rescorer named by `model`.
///
/// All six rescorers share this signature (candidates in, scored+q-valued
/// results plus per-fold feature stats out), so this is the single place model
/// selection is resolved.
///
/// NOTHING TYPE-CHECKS ANY OF THESE ARMS: every rescorer has this signature, so
/// an inverted arm compiles, runs, and produces plausible q-values off the wrong
/// model or the wrong feature set. The MLP entry points spell their lane in their
/// name for exactly that reason (see `qvalues`'s note above
/// `rescore_mlp_linear`), and both halves of the mapping are pinned as tests:
/// `rescore_with_dispatches_each_mlp_variant_to_its_lane` for the three MLP
/// arms, `rescore_with_dispatches_the_default_and_the_non_mlp_variants` for
/// `Gbm` / `Lda` / `Hybrid` — the latter matters most, since `Gbm` is what runs
/// when nobody sets the flag.
pub fn rescore_with(
    model: RescoreModel,
    data: Vec<CompetedCandidate>,
) -> (Vec<FinalResult>, RescoreFeatureStats) {
    match model {
        RescoreModel::Gbm => rescore(data),
        RescoreModel::Lda => rescore_lda(data),
        RescoreModel::Hybrid => rescore_hybrid(data),
        RescoreModel::Mlp => rescore_mlp_linear(data),
        RescoreModel::MlpAll => rescore_mlp_all(data),
        RescoreModel::HybridMlp => rescore_hybrid_mlp(data),
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
