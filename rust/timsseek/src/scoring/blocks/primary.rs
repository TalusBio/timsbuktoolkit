//! Primary-scores family: main ranking score and peak-separation deltas.

use crate::score_block;

score_block! {
    /// Stage: apex (main_score is the weighted product at the joint apex;
    /// deltas compare against the 2nd/3rd global peaks).
    pub struct PrimaryScores {
        #[raw] pub main_score: f32,
        #[raw] pub delta_next: f32,
        #[raw] pub delta_second_next: f32,
    }
}
