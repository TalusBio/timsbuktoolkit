//! Summed-intensity family (MS1/MS2 apex sums). Both are ML-projected as
//! `ln1p` only (the raw sums are Parquet columns).

use timsseek_macros::ScoreBlock;

/// Stage: apex (intensity sums at the joint apex).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct Intensities {
    #[feat(ln1p)]
    pub ms2_summed_intensity: f32,
    #[feat(ln1p)]
    pub ms1_summed_intensity: f32,
}
