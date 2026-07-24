//! Result-meta family — post-model, output-only fields (filled in place after
//! the GBM, not computed via Inputs). `columns` emits all four; `features`
//! emits only the two delta-group fields (`discriminant_score`/`qvalue` are
//! Parquet-only).

use timsseek_macros::ScoreBlock;

/// Stage: post-model (output-only).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct ResultMeta {
    #[feat(raw)]
    pub delta_group: f32,
    #[feat(raw)]
    pub delta_group_ratio: f32,
    pub discriminant_score: f32,
    pub qvalue: f32,
}
