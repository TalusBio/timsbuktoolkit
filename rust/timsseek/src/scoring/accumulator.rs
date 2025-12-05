//! Accumulator for parallel scoring results.
//!
//! This module provides efficient collection of scoring results from parallel iterators.
//! It aggregates both successful search results and timing measurements across threads.

use super::search_results::IonSearchResults;
use super::timings::ScoreTimings;
use rayon::iter::{
    FromParallelIterator,
    IntoParallelIterator,
    ParallelIterator,
};

/// Accumulator for collecting scoring results and timings from parallel operations.
///
/// This struct is used internally by `score_iter()` to efficiently collect results
/// from Rayon's parallel iterator. It implements both sequential and parallel collection
/// via `FromIterator` and `FromParallelIterator`.
///
/// # Design
///
/// The accumulator uses the fold-reduce pattern for parallel aggregation:
/// 1. **Fold**: Each thread accumulates results into a local accumulator
/// 2. **Reduce**: Local accumulators are merged pairwise to produce final result
///
/// This minimizes contention and allows efficient parallel collection of results.
#[derive(Default)]
pub(super) struct IonSearchAccumulator {
    pub(super) res: Vec<IonSearchResults>,
    pub(super) timings: ScoreTimings,
}

impl IonSearchAccumulator {
    /// Merges two accumulators by extending results and summing timings.
    ///
    /// Used in the reduce phase of parallel collection.
    pub(super) fn reduce(mut self, other: Self) -> Self {
        self.res.extend(other.res);
        self.timings += other.timings;
        self
    }

    /// Adds a single scoring result (if present) to the accumulator.
    ///
    /// Used in the fold phase of parallel collection. Successful results are collected,
    /// while `None` results (failed scoring) are discarded.
    pub(super) fn fold(mut self, item: (Option<IonSearchResults>, ScoreTimings)) -> Self {
        if let Some(elem) = item.0 {
            self.res.push(elem);
        }
        self.timings += item.1;
        self
    }
}

impl FromIterator<(Option<IonSearchResults>, ScoreTimings)> for IonSearchAccumulator {
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = (Option<IonSearchResults>, ScoreTimings)>,
    {
        iter.into_iter()
            .fold(IonSearchAccumulator::default(), IonSearchAccumulator::fold)
    }
}

impl FromParallelIterator<(Option<IonSearchResults>, ScoreTimings)> for IonSearchAccumulator {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = (Option<IonSearchResults>, ScoreTimings)>,
    {
        par_iter
            .into_par_iter()
            .fold(IonSearchAccumulator::default, IonSearchAccumulator::fold)
            .reduce(IonSearchAccumulator::default, IonSearchAccumulator::reduce)
    }
}
