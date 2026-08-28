//! Result-meta family -- post-model, output-only fields (filled in place after
//! rescoring, not computed via Inputs). `columns` emits all four; only the two
//! delta-group fields carry `#[feat(raw)]`, so they are the only ones in a
//! feature lane -- the LINEAR one, since `linear = true` is the default
//! (`discriminant_score`/`qvalue` are Parquet-only).

use timsseek_macros::ScoreBlock;

/// Stage: post-model (output-only).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct ResultMeta {
    /// `ln_1p(best_score) - ln_1p(runner_up_score)` within the competition group.
    #[feat(raw)]
    pub delta_group_ln1p_diff: f32,
    /// `ln_1p(runner_up_score) / ln_1p(best_score)` within the competition group.
    #[feat(raw)]
    pub delta_group_ln1p_ratio: f32,
    pub discriminant_score: f32,
    pub qvalue: f32,
}
