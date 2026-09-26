//! Result-meta family -- post-model fields written to Parquet.
//! Exploration E1 retains the group margins in output but omits all four fields
//! from rescorer feature lanes.

use timsseek_macros::ScoreBlock;

/// Stage: post-model (output-only).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct ResultMeta {
    /// `ln_1p(best_score) - ln_1p(runner_up_score)` within the competition group.
    pub delta_group_ln1p_diff: f32,
    /// `ln_1p(runner_up_score) / ln_1p(best_score)` within the competition group.
    pub delta_group_ln1p_ratio: f32,
    pub discriminant_score: f32,
    pub qvalue: f32,
}
