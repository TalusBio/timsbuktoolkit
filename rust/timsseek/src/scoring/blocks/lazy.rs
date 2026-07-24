//! Lazyscore families. Split by data-death: `ApexLazyScores` is computed from
//! the chromatogram traces (apex stage); `SecondaryLazyScores` is computed from
//! the secondary spectral-query collectors (finalize stage).
//!
//! Canary block for the `#[derive(ScoreBlock)]` migration (see
//! `timsseek_macros`): both structs route their fields through the new
//! `column_schema`/`linear_features`/`linear_feature_names` lane methods
//! instead of the legacy `score_block!`-generated `features`/`feature_names`
//! walk, which the derive leaves defaulted (no-op) for these two blocks.

use serde::Serialize;
use timsseek_macros::ScoreBlock;

use crate::scoring::pipeline::SecondaryLazyScoresRaw;

/// Stage: apex (traces are dead by finalize).
#[derive(Debug, Clone, Copy, Serialize, ScoreBlock)]
pub struct ApexLazyScores {
    #[feat(raw)]
    pub apex_lazyscore: f32,
    #[feat(raw)]
    pub lazyscore_z: f32,
    #[feat(raw)]
    pub lazyscore_vs_baseline: f32,
}

/// Stage: finalize (from the secondary-query inner/isotope collectors).
#[derive(Debug, Clone, Copy, Serialize, ScoreBlock)]
pub struct SecondaryLazyScores {
    #[feat(raw)]
    pub ms2_lazyscore: f32,
    #[feat(raw)]
    pub ms2_isotope_lazyscore: f32,
    #[feat(ln1p)]
    pub ms2_isotope_lazyscore_ratio: f32,
}

impl From<SecondaryLazyScoresRaw> for SecondaryLazyScores {
    fn from(s: SecondaryLazyScoresRaw) -> Self {
        Self {
            ms2_lazyscore: s.lazyscore,
            ms2_isotope_lazyscore: s.iso_lazyscore,
            ms2_isotope_lazyscore_ratio: s.ratio,
        }
    }
}
