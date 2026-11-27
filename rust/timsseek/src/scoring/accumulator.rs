//! Accumulator for parallel scoring results.
//!
//! Aggregates successful `ScoredCandidate`s, per-reason skip counts, and
//! phase timings across Rayon workers using fold-reduce.

use super::results::ScoredCandidate;
use super::skip::{
    SkipCounts,
    SkipReason,
};
use super::timings::ScoreTimings;
pub(super) type ScoreItem = (Result<ScoredCandidate, SkipReason>, ScoreTimings);

#[derive(Default)]
pub(super) struct IonSearchAccumulator {
    pub(super) res: Vec<ScoredCandidate>,
    pub(super) skips: SkipCounts,
    pub(super) timings: ScoreTimings,
}

impl IonSearchAccumulator {
    pub(super) fn merge(&mut self, other: Self) {
        self.res.extend(other.res);
        self.skips += other.skips;
        self.timings += other.timings;
    }

    pub(super) fn fold(&mut self, item: ScoreItem) {
        match item.0 {
            Ok(elem) => self.res.push(elem),
            Err(reason) => self.skips.bump(reason),
        }
        self.timings += item.1;
    }
}
