use crate::data_sources::speclib::ExpectedIntensities;
use crate::errors::{
    DataProcessingError,
    TimsSeekError,
};
use crate::fragment_mass::IonAnnot;
use crate::models::DigestSlice;
use rayon::prelude::*;
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};
use timsquery::{
    EGCAggregator,
    ElutionGroup,
    GenerallyQueriable,
    QuadSplittedTransposedIndex,
    Tolerance,
};

use super::calculate_scores::{
    IntensityArrays,
    LocalizedPreScore,
    PreScore,
};
use super::search_results::{
    IonSearchResults,
    SearchResultBuilder,
};

#[derive(Debug, Clone)]
pub struct QueryItemToScore {
    pub digest: DigestSlice,
    pub charge: u8,
    pub querie: Arc<ElutionGroup<IonAnnot>>,
    pub expected_intensity: ExpectedIntensities,
}

// Metrics specific to the scoring process
#[derive(Debug, Clone, Default, Copy)]
pub struct ScoringMetrics {
    pub num_processed: usize,
    pub num_skipped_localization: usize,
    pub num_skipped_scoring: usize,
    pub time_localizing: Duration,
    pub time_scoring: Duration,
}

pub struct Scorer<I: GenerallyQueriable<IonAnnot>> {
    index_cycle_rt_ms: Arc<[u32]>,
    index: I,
    tolerance: Tolerance,
}

// Private helper methods for individual scoring steps
impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    #[inline] // Potentially inline small helper methods
    fn _build_prescore(&self, item: QueryItemToScore) -> PreScore {
        let mut agg =
            EGCAggregator::new(item.querie.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        self.index.add_query(&mut agg, &self.tolerance);

        let prescore = PreScore {
            charge: item.charge,
            digest: item.digest,
            expected_intensities: item.expected_intensity,
            query_values: agg,
        };
        prescore
    }

    // This internal version takes the buffer for optimal reuse in parallel contexts
    #[inline]
    fn _localize_step(
        &self,
        prescore: PreScore,
        buffer: &mut IntensityArrays,
    ) -> Result<LocalizedPreScore, TimsSeekError> {
        prescore.localize_with_buffer(buffer)
        // Error handling can be done here or propagated
    }

    #[inline]
    fn _finalize_step(
        &self,
        localized_score: &LocalizedPreScore,
    ) -> Result<IonSearchResults, DataProcessingError> {
        let builder = SearchResultBuilder::default();
        builder.with_localized_pre_score(localized_score).finalize()
    }
}

impl Scorer<QuadSplittedTransposedIndex> {
    pub fn from_quad_splitted(index: QuadSplittedTransposedIndex, tolerance: Tolerance) -> Self {
        let index_cycle_rt_ms = index.cycle_rt_ms.clone();
        Self {
            tolerance,
            index,
            index_cycle_rt_ms,
        }
    }
}

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    /// Scores a single query item by orchestrating the internal steps.
    /// Useful for testing or single-item processing scenarios.
    pub fn score(&self, item: QueryItemToScore) -> Option<IonSearchResults> {
        // Allocate buffer internally for this single use case
        let mut buffer = match IntensityArrays::new_empty(5, 10, self.index_cycle_rt_ms.clone()) {
            Ok(b) => b,
            Err(_) => return None, // Or handle error appropriately
        };

        let prescore = self._build_prescore(item);
        match self._localize_step(prescore, &mut buffer) {
            Ok(localized) => self._finalize_step(&localized).ok(),
            Err(_e) => {
                // Log error if needed
                None
            }
        }
    }

    /// Scores a collection of query items efficiently in parallel.
    /// Handles parallelization and buffer management internally.
    pub fn score_iter(
        &self,
        items_to_score: Vec<QueryItemToScore>, // Takes Vec for easy parallel iteration
    ) -> (Vec<IonSearchResults>, ScoringMetrics) {
        // Returns results and metrics

        let mut metrics = ScoringMetrics::default();
        let num_input_items = items_to_score.len();
        let loc_score_start = Instant::now(); // Combined timing for iter

        let init_fn = || IntensityArrays::new_empty(5, 10, self.index_cycle_rt_ms.clone()).unwrap();

        let results: Vec<IonSearchResults> = items_to_score
            .into_par_iter()
            .map_init(
                // Use map_init for efficient buffer reuse
                init_fn,
                |buffer, item| {
                    // Closure takes buffer and item
                    let prescore = self._build_prescore(item);
                    match self._localize_step(prescore, buffer) {
                        // Reuse buffer
                        Ok(localized) => {
                            match self._finalize_step(&localized) {
                                Ok(final_result) => Some(final_result),
                                Err(_e) => None, // Filter out finalization errors
                            }
                        }
                        Err(_e) => None, // Filter out localization errors
                    }
                },
            )
            .flatten() // Remove the None values (errors)
            .collect();

        let elapsed = loc_score_start.elapsed();
        metrics.num_processed = results.len();
        // Note: Calculating skipped counts accurately might need more atomic counters
        // within the closure if you need distinct localization vs finalization skips.
        // This gives the total number skipped for any reason.
        let num_skipped = num_input_items.saturating_sub(results.len());

        // Assign time - might need refinement if separate timings are critical
        metrics.time_localizing = elapsed / 2; // Crude split example
        metrics.time_scoring = elapsed / 2; // Needs better measurement if important

        (results, metrics)
    }

    // Optional: score_full_iter(...) similar to score_iter but for full results
}

// // --- Calling Code ---
// // Remains simple!
//
// // ... setup ...
// let scorer = Scorer::new(index.cycle_rt_ms.clone());
// // ... query index ...
// let items_to_score: Vec<QueryItemToScore> = ... ; // Prepare items
//
// // Single call to process the whole chunk efficiently
// let (ion_search_results, scoring_metrics) = scorer.score_iter(items_to_score);
//
// // Combine metrics
// let chunk_metrics = RuntimeMetrics {
//     num_skipped: num_skipped_in_query + (scoring_metrics.num_processed - num_input_items), // Example calculation
//     time_localizing: scoring_metrics.time_localizing,
//     time_scoring: scoring_metrics.time_scoring,
//     // ... other fields ...
// };
//
