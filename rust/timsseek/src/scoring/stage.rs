//! What each pipeline stage ranks its candidates by.
//!
//! The apex finder says where an apex is ([`super::apex_finding::ApexLocation`]);
//! a stage says how candidates compare across precursors. The two are kept
//! apart in the types: only [`super::apex_finding::TraceScorer::rank_at`] makes
//! a [`RankScore`], so a stage cannot be fed a value that was never meant to
//! be compared across extractions.

use std::cmp::Ordering;

use crate::errors::DataProcessingError;

/// A pipeline stage that ranks candidates by one scalar.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    /// Phase 1: calibrant selection over the broad, full-run extraction.
    Calibrant,
    /// Phase 3: final scoring over the narrow, calibrated extraction.
    Final,
}

impl Stage {
    pub const ALL: [Stage; 2] = [Stage::Calibrant, Stage::Final];

    /// The only label a report prints for a stage.
    pub fn name(self) -> &'static str {
        match self {
            Self::Calibrant => "calibrant",
            Self::Final => "final",
        }
    }
}

/// The scalar a [`Stage`] ranks by. Finite, and ordered across the candidates
/// of one stage in one run. Built only by the scorer.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct RankScore(f32);

impl RankScore {
    pub(super) fn new(value: f32) -> Result<Self, DataProcessingError> {
        if value.is_finite() {
            Ok(Self(value))
        } else {
            Err(DataProcessingError::ExpectedFiniteNonNanData {
                context: format!("rank score is not finite: {value}"),
            })
        }
    }

    pub fn get(self) -> f32 {
        self.0
    }

    /// Total order; equal to `partial_cmp` because the value is finite.
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }

    /// A fixed score for tests that need a candidate but not a scorer.
    ///
    /// # Panics
    ///
    /// Panics if `value` is not finite.
    #[cfg(any(test, feature = "test-support"))]
    pub fn for_tests(value: f32) -> Self {
        Self::new(value).expect("test rank score is finite")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_score_rejects_non_finite_values() {
        assert!(RankScore::new(f32::NAN).is_err());
        assert!(RankScore::new(f32::INFINITY).is_err());
        assert!(RankScore::new(f32::NEG_INFINITY).is_err());
        assert_eq!(RankScore::new(-2.5).unwrap().get(), -2.5);
    }

    #[test]
    fn stage_names_are_distinct() {
        assert_eq!(Stage::ALL.map(Stage::name), ["calibrant", "final"]);
    }
}
