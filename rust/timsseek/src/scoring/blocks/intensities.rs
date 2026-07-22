//! Summed-intensity family (MS1/MS2 apex sums). Both are ML-projected as
//! `ln1p` only (the raw sums are Parquet columns).

use crate::score_block;

score_block! {
    /// Stage: apex (intensity sums at the joint apex).
    pub struct Intensities {
        #[feat(ln1p)] pub ms2_summed_intensity: f32,
        #[feat(ln1p)] pub ms1_summed_intensity: f32,
    }
}
