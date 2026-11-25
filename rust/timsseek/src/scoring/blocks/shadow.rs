//! Single-cycle fragment evidence after subtracting a co-isolated candidate.

use timsseek_macros::ScoreBlock;

#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct ShadowScores {
    /// Cosine of predicted and residual fragment intensities at the apex.
    #[feat(raw)]
    pub shadow_cosine: f32,
    /// Sum of residual matched intensity divided by the original matched sum.
    #[feat(raw)]
    pub shadow_retained: f32,
    /// Fraction of observed subject channels whose residual intensity is zero.
    #[feat(raw)]
    pub shadow_lost: f32,
}

impl Default for ShadowScores {
    fn default() -> Self {
        Self {
            shadow_cosine: f32::NAN,
            shadow_retained: f32::NAN,
            shadow_lost: f32::NAN,
        }
    }
}
