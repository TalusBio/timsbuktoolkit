//! Lazyscore families. Split by data-death: `ApexLazyScores` is computed from
//! the chromatogram traces (apex stage); `SecondaryLazyScores` is computed from
//! the secondary spectral-query collectors (finalize stage).

use crate::score_block;
use crate::scoring::pipeline::SecondaryLazyScoresRaw;

score_block! {
    /// Stage: apex (traces are dead by finalize).
    pub struct ApexLazyScores {
        #[raw] pub apex_lazyscore: f32,
        #[raw] pub lazyscore_z: f32,
        #[raw] pub lazyscore_vs_baseline: f32,
    }
}

score_block! {
    /// Stage: finalize (from the secondary-query inner/isotope collectors).
    pub struct SecondaryLazyScores {
        #[raw] pub ms2_lazyscore: f32,
        #[raw] pub ms2_isotope_lazyscore: f32,
        #[raw] pub ms2_isotope_lazyscore_ratio: f32,
    }
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
