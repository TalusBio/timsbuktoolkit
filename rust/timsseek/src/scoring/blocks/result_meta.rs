//! Result-meta family — post-model, output-only fields (filled in place after
//! the GBM, not computed via Inputs). `columns` emits all four; `features`
//! emits only the two delta-group fields (`discriminant_score`/`qvalue` are
//! Parquet-only).

use crate::score_block;

score_block! {
    /// Stage: post-model (output-only).
    pub struct ResultMeta {
        #[raw] pub delta_group: f32,
        #[raw] pub delta_group_ratio: f32,
        #[col_only] pub discriminant_score: f32,
        #[col_only] pub qvalue: f32,
    }
}
