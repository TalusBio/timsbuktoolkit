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
///
/// # Cost — THE CANONICAL TABLE, AND WHY IT IS A TABLE
/// **No variant here has a single cost multiplier.** Every ratio below moves
/// with the candidate count, in both directions and by more than a factor of
/// two, so a lone "~Nx" is wrong at every scale but one. Successive versions of
/// this documentation each picked a number and each rotted; the table exists so
/// the next reader interpolates between anchors instead.
///
/// Phase 5 wall time / targets at 1% FDR, two runs of one Astral mzML:
///
/// | model | 114138 candidates | 2174837 candidates |
/// |---|---|---|
/// | `lda` | 0.476s / 24141 | 9.667s / 102843 |
/// | `hybrid` | 3.455s / 24407 | 148.062s / 110425 |
/// | `mlp` | 6.141s / 25077 | 323.012s / 125242 — OLD DEFAULTS |
/// | `mlp_all` | 7.342s / 25451 | 281.953s / 134655 |
/// | `gbm` | 9.588s / 25240 | 574.315s / 134710 |
/// | `hybrid_mlp` | 11.258s / 25454 | 773.560s / 133977 — OLD DEFAULTS |
///
/// EVERY 114k CELL IS AT THE CURRENT [`mlp::MlpConfig`] DEFAULTS. The two cells
/// marked OLD DEFAULTS are not: `mlp` and `hybrid_mlp` have not been re-timed at
/// 2.17M since the retune, and the retune cut both of their 114k times by more
/// than half (`mlp` 14.143s -> 6.141s, `hybrid_mlp` 37.831s -> 11.258s at that
/// scale), so read those two mid-size numbers as loose upper bounds, not as
/// current cost. `lda`, `hybrid` and `gbm` never construct an `MlpConfig`, so the
/// retune cannot have moved them at either scale.
///
/// ## The speed ratio GROWS with candidate count
/// `lda` against `gbm`: ~20x cheaper at 114k, ~59x cheaper at 2.17M, and ~100x
/// on a ~28M-entry library measured earlier in a different regime. Three anchors
/// on one curve, not three attempts at one constant.
///
/// The MLP variants can cross the GBM entirely, and `mlp` is the sharpest case
/// on record: at 114k candidates it took 14.143s against the GBM's 9.588s at the
/// pre-retune defaults (1.5x SLOWER) and takes 6.141s at the current ones (1.56x
/// FASTER). Same variant, same input, same scale — the comparison changed sign
/// from hyperparameter tuning alone. `mlp_all` shows the scale axis of the same
/// effect: 1.73x the GBM's Phase 5 time at 114k (slower) and 0.85x at 2.17M
/// pre-retune, 0.77x and 0.49x retuned. So "the MLP is slower than the GBM" is a
/// question about your library size AND which defaults you are on, never a
/// property of the model.
///
/// ## The LDA's ACCURACY cost grows with scale too, and it is the bigger trap
/// `lda` gives up 4.4% of the GBM's targets at 114k candidates but **23.7%** at
/// 2.17M; `hybrid` -3.3% then -18.0%. Someone who picks `lda` for its speed
/// after a small-run comparison is giving up several times what that comparison
/// showed them.
///
/// ## Caveats that survive every number above
/// * Both runs are the same FAIMS acquisition, i.e. NO ion mobility: 25 of the
///   101 linear and 27 of the 128 all-lane columns were culled as dead. A TIMS
///   run exercises more features and may time differently.
/// * Two datasets, one machine. Orders of magnitude and directions of movement,
///   not a benchmark suite.
/// * Phase 5 is a fraction of a whole search, so 2x here is not 2x end to end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum RescoreModel {
    /// Gradient-boosted trees over the ALL lane (linear ++ nonlinear, 128
    /// features), cross-validated over `N_RESCORE_FOLDS`. The default, and the
    /// most sensitive option measured at both scales (25240 and 134710 targets at
    /// 1% FDR) — though [`MlpAll`](RescoreModel::MlpAll) now lands within 0.04% of
    /// it for half the Phase 5 time on the larger of the two. See the cost table
    /// above.
    #[default]
    Gbm,
    /// Sage-style shrinkage LDA on the LINEAR lane only (101 features, see
    /// [`lda`]). Cross-fit over the same fold ASSIGNMENT.
    ///
    /// THE CHEAPEST AND THE LEAST SENSITIVE, and both gaps widen with the
    /// candidate count: 20x cheaper for -4.4% targets at 114k candidates, 59x
    /// cheaper for **-23.7%** at 2.17M. See the type-level cost table, and
    /// [`rescore_lda`] for what the measurement does and does not cover.
    Lda,
    /// Cross-fit LDA (linear lane) -> `lda_score` column -> GBM CV on
    /// `nonlinear + lda_score`: the same fold ASSIGNMENT on both sides, so a
    /// row's `lda_score` never saw its own label.
    ///
    /// Between `lda` and `gbm` on both axes, and NOT at parity with the GBM the
    /// way this was once documented: -3.3% targets at 114k candidates, where it
    /// took 0.36x the GBM's Phase 5 wall time (3.455s against 9.588s), and -18.0%
    /// at 2.17M for 0.26x. Unaffected by the [`mlp::MlpConfig`] retune — the
    /// compressor here is the LDA, so both figures are current.
    Hybrid,
    /// MLP on the LINEAR lane only (101 features) — the same feature set
    /// [`Lda`](RescoreModel::Lda) trains on, so it isolates "nonlinear model"
    /// from "more features". THE DEFAULT MLP SHAPE, hence the unqualified name.
    ///
    /// THE MOST DRAMATIC BENEFICIARY OF THE RETUNE: at 114k candidates it went
    /// from 14.143s to 6.141s of Phase 5 time, i.e. from 1.5x SLOWER than
    /// [`Gbm`](RescoreModel::Gbm) to 1.56x faster at that scale, on tuning alone.
    /// Its 2.17M figure (323.012s / 125242) is still an OLD-DEFAULTS number and
    /// has not been re-measured. Sensitivity is the trade: -0.6% of the GBM's
    /// targets at 114k, and -7.0% at 2.17M before the retune. See the cost table
    /// above.
    Mlp,
    /// MLP on the ALL lane (linear ++ nonlinear, 128 features) — the same
    /// feature set [`Gbm`](RescoreModel::Gbm) trains on, so the two are
    /// directly comparable.
    ///
    /// THE ARM THE [`mlp::MlpConfig`] DEFAULTS WERE TUNED ON, and the only MLP
    /// variant timed at those defaults at BOTH scales — on those two runs the best
    /// speed/sensitivity point measured: 1.31x faster than `gbm` for +211 targets
    /// at 114k candidates, 2.04x faster for -55 (-0.04%) at 2.17M.
    MlpAll,
    /// MLP counterpart of [`Hybrid`](RescoreModel::Hybrid): cross-fit MLP
    /// (linear lane) -> `mlp_score` column -> GBM CV on
    /// `nonlinear + mlp_score`, same fold ASSIGNMENT on both sides.
    ///
    /// THE MOST EXPENSIVE VARIANT AT BOTH SCALES MEASURED, structurally: it pays
    /// for two models, so the GBM's Phase 5 time is a floor under it however the
    /// MLP half is tuned. At 114k candidates it costs 11.258s for 25454 targets
    /// (+214 over `gbm`, the most of any variant, for 1.17x its time) — retuning
    /// took it from 37.831s to that, so the "structurally hopeless" verdict it
    /// once carried was really a verdict on the old hyperparameters. Its 2.17M
    /// figure (773.560s / 133977) is an OLD-DEFAULTS number and will improve; do
    /// not judge it on that row.
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
