//! Split-product family: cosine/scribe apex agreement scores.
//!
//! The struct is macro-generated; `compute` lives in
//! [`crate::scoring::scores::apex_features::compute_split_product`] (needs the
//! chromatogram buffers, so it runs at the apex stage). This block is the
//! typed projection of that computation, built via [`From`].

use crate::score_block;
use crate::scoring::scores::apex_features::SplitProductScore;

score_block! {
    /// Stage: apex (built from `SplitProductScore`, computed while chromatogram
    /// buffers are live).
    pub struct SplitProduct {
        #[feat(ln1p)] pub split_product_score: f32,
        #[feat(ln1p)] pub cosine_au: f32,
        #[feat(ln1p)] pub scribe_au: f32,
        #[raw] pub cosine_cg: f32,
        #[raw] pub scribe_cg: f32,
        #[raw] pub cosine_weighted_coelution: f32,
        #[raw] pub cosine_gradient_consistency: f32,
        #[raw] pub scribe_weighted_coelution: f32,
        #[raw] pub scribe_gradient_consistency: f32,
    }
}

impl From<&SplitProductScore> for SplitProduct {
    fn from(s: &SplitProductScore) -> Self {
        Self {
            split_product_score: s.base_score,
            cosine_au: s.cosine_au,
            scribe_au: s.scribe_au,
            cosine_cg: s.cosine_cg,
            scribe_cg: s.scribe_cg,
            cosine_weighted_coelution: s.cosine_weighted_coelution,
            cosine_gradient_consistency: s.cosine_gradient_consistency,
            scribe_weighted_coelution: s.scribe_weighted_coelution,
            scribe_gradient_consistency: s.scribe_gradient_consistency,
        }
    }
}
