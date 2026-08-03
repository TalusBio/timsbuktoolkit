//! Lazyscore families. Split by data-death: `ApexLazyScores` is computed from
//! the chromatogram traces (apex stage); `SecondaryLazyScores` is computed from
//! the secondary spectral-query collectors (finalize stage).
//!
//! Both structs are ordinary `#[derive(ScoreBlock)]` blocks (see
//! `timsseek_macros`): the derive generates every projection
//! (`column_schema`/`columns` plus both feature lanes' value arrays and names)
//! from the field list. Every field here is `#[feat(..)]`-annotated without
//! `linear = false`, so they all land in the LINEAR lane — the live path the
//! ML consumer reads.

use serde::Serialize;
use timsseek_macros::ScoreBlock;

use crate::scoring::pipeline::SecondaryLazyScoresRaw;

/// Stage: apex (traces are dead by finalize).
#[derive(Debug, Clone, Copy, Serialize, ScoreBlock)]
pub struct ApexLazyScores {
    #[feat(ln1p)]
    pub apex_lazyscore: f32,
    #[feat(raw)]
    pub lazyscore_z: f32,
    #[feat(ln1p)]
    pub lazyscore_vs_baseline: f32,
}

/// Stage: finalize (from the secondary-query inner/isotope collectors).
#[derive(Debug, Clone, Copy, Serialize, ScoreBlock)]
pub struct SecondaryLazyScores {
    #[feat(ln1p)]
    pub ms2_lazyscore: f32,
    #[feat(ln1p)]
    pub ms2_isotope_lazyscore: f32,
    #[feat(raw)]
    pub ms2_isotope_lazyscore_log_diff: f32,
}

impl From<SecondaryLazyScoresRaw> for SecondaryLazyScores {
    fn from(s: SecondaryLazyScoresRaw) -> Self {
        Self {
            ms2_lazyscore: s.lazyscore,
            ms2_isotope_lazyscore: s.iso_lazyscore,
            ms2_isotope_lazyscore_log_diff: s.isotope_lazyscore_log_diff,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::blocks::ScoreBlock;

    #[test]
    fn isotope_log_difference_is_not_transformed_twice() {
        let raw = SecondaryLazyScoresRaw {
            lazyscore: 2.0,
            iso_lazyscore: 100.0,
            isotope_lazyscore_log_diff: 2.0f32.ln_1p() - 100.0f32.ln_1p(),
        };
        let expected = raw.isotope_lazyscore_log_diff as f64;
        assert!(expected < -1.0, "fixture must expose ln1p's invalid domain");

        let values = SecondaryLazyScores::from(raw).linear_feature_array();
        assert_eq!(values[2], expected);
        assert!(values[2].is_finite());

        let mut names = crate::scoring::blocks::NameSink::new();
        SecondaryLazyScores::linear_feature_names(&mut names);
        assert_eq!(&*names.into_names()[2], "ms2_isotope_lazyscore_log_diff");
    }
}
