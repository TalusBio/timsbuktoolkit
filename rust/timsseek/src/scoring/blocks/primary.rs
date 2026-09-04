//! Primary-scores family: main ranking score and peak-separation deltas.

use timsseek_macros::ScoreBlock;

/// Stage: apex (main_score is the weighted product at the joint apex;
/// deltas compare against the 2nd/3rd global peaks).
#[derive(Debug, Clone, Copy, ::serde::Serialize, ScoreBlock)]
pub struct PrimaryScores {
    #[feat(log2)]
    pub main_score: f32,
    #[feat(log2)]
    pub delta_next: f32,
    #[feat(log2)]
    pub delta_second_next: f32,
}
