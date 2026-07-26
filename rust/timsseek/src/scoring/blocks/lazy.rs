//! Lazyscore families. Split by data-death: `ApexLazyScores` is computed from
//! the chromatogram traces (apex stage); `SecondaryLazyScores` is computed from
//! the secondary spectral-query collectors (finalize stage).
//!
//! Both structs are ordinary `#[derive(ScoreBlock)]` blocks (see
//! `timsseek_macros`): the derive generates all six projection methods
//! (`column_schema`/`columns` plus both feature lanes' values and names) from
//! the field list. Every field here is `#[feat(..)]`-annotated without
//! `linear = false`, so they all land in the LINEAR lane — the live path the
//! ML consumer reads via `FeatFrame`.

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
