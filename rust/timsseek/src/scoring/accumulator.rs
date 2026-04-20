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
#[cfg(feature = "rayon")]
use rayon::iter::{
    FromParallelIterator,
    IntoParallelIterator,
    ParallelIterator,
};

pub(super) type ScoreItem = (Result<ScoredCandidate, SkipReason>, ScoreTimings);

#[derive(Default)]
pub(super) struct IonSearchAccumulator {
    pub(super) res: Vec<ScoredCandidate>,
    pub(super) skips: SkipCounts,
    pub(super) timings: ScoreTimings,
}

impl IonSearchAccumulator {
    pub(super) fn reduce(mut self, other: Self) -> Self {
        self.res.extend(other.res);
        self.skips += other.skips;
        self.timings += other.timings;
        self
    }

    pub(super) fn fold(mut self, item: ScoreItem) -> Self {
        match item.0 {
            Ok(elem) => self.res.push(elem),
            Err(reason) => self.skips.bump(reason),
        }
        self.timings += item.1;
        self
    }
}

impl FromIterator<ScoreItem> for IonSearchAccumulator {
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = ScoreItem>,
    {
        iter.into_iter()
            .fold(IonSearchAccumulator::default(), IonSearchAccumulator::fold)
    }
}

#[cfg(feature = "rayon")]
impl FromParallelIterator<ScoreItem> for IonSearchAccumulator {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = ScoreItem>,
    {
        par_iter
            .into_par_iter()
            .fold(IonSearchAccumulator::default, IonSearchAccumulator::fold)
            .reduce(IonSearchAccumulator::default, IonSearchAccumulator::reduce)
    }
}
