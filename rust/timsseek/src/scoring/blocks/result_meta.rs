//! Result metadata written to Parquet after rescoring.

use timsseek_macros::ScoreBlock;

/// Stage: post-model (output-only).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct ResultMeta {
    pub discriminant_score: f32,
    pub qvalue: f32,
}
