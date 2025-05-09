use crate::errors::DataProcessingError;
use crate::{
    IonAnnot,
    QueryItemToScore,
};
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;
use timsquery::{
    ChromatogramCollector,
    GenerallyQueriable,
    Tolerance,
};

use super::calculate_scores::{
    IntensityArrays,
    LocalizedPreScore,
    LongitudinalMainScoreElements,
    PreScore,
};
use super::full_results::FullQueryResult;
use super::search_results::{
    IonSearchResults,
    SearchResultBuilder,
};
use tracing::info;

pub struct Scorer<I: GenerallyQueriable<IonAnnot>> {
    pub index_cycle_rt_ms: Arc<[u32]>,
    pub index: I,
    pub tolerance: Tolerance,
}

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    // does inlining do anything here?
    #[inline]
    fn _build_prescore(&self, item: &QueryItemToScore) -> PreScore {
        let mut agg =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        self.index.add_query(&mut agg, &self.tolerance);

        PreScore {
            charge: item.charge,
            digest: item.digest.clone(),
            expected_intensities: item.expected_intensity.clone(),
            query_values: agg,
        }
    }

    // This internal version takes the buffer for optimal reuse in parallel contexts
    #[inline]
    fn _localize_step(
        &self,
        prescore: PreScore,
        buffer: &mut IntensityArrays,
    ) -> Result<LocalizedPreScore, DataProcessingError> {
        prescore.localize_with_buffer(buffer)
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

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    /// Scores a single query item by orchestrating the internal steps.
    /// Useful for testing or single-item processing scenarios.
    pub fn score(&self, item: QueryItemToScore) -> Option<IonSearchResults> {
        let mut buffer = IntensityArrays::new_empty(5, 10, self.index_cycle_rt_ms.clone()).ok()?;
        let prescore = self._build_prescore(&item);
        let localized = self._localize_step(prescore, &mut buffer).ok()?;
        self._finalize_step(&localized).ok()
    }

    pub fn score_full(
        &self,
        queries: QueryItemToScore,
    ) -> Result<FullQueryResult, DataProcessingError> {
        let mut res =
            ChromatogramCollector::new(queries.query.clone(), self.index_cycle_rt_ms.clone())
                .unwrap();
        self.index.add_query(&mut res, &self.tolerance);
        let builder = SearchResultBuilder::default();
        let int_arrs = IntensityArrays::new(&res, &queries.expected_intensity)?;
        let prescore = PreScore {
            charge: queries.charge,
            digest: queries.digest,
            expected_intensities: queries.expected_intensity,
            query_values: res.clone(),
        };

        let longitudinal_main_score_elements = LongitudinalMainScoreElements::try_new(&int_arrs)?;

        let res2 = builder
            .with_localized_pre_score(&prescore.localize()?)
            .finalize()?;
        let longitudinal_main_score = longitudinal_main_score_elements.main_score_iter().collect();

        Ok(FullQueryResult {
            main_score_elements: longitudinal_main_score_elements,
            longitudinal_main_score,
            extractions: res.clone(),
            search_results: res2,
        })
    }

    /// Scores a collection of query items efficiently in parallel.
    /// Handles parallelization and buffer management internally.
    pub fn score_iter(
        &self,
        items_to_score: &[QueryItemToScore], // Takes Vec for easy parallel iteration
    ) -> Vec<IonSearchResults> {
        let num_input_items = items_to_score.len();
        let loc_score_start = Instant::now(); // Combined timing for iter

        let init_fn = || IntensityArrays::new_empty(5, 10, self.index_cycle_rt_ms.clone()).unwrap();

        let results: Vec<IonSearchResults> = items_to_score
            .into_par_iter()
            .map_init(init_fn, |buffer, item| {
                let prescore = self._build_prescore(item);
                match self._localize_step(prescore, buffer) {
                    // Reuse buffer
                    Ok(localized) => self._finalize_step(&localized).ok(),
                    Err(_e) => None, // TODO: LOG
                }
            })
            .flatten() // Remove the None values (errors)
            .collect();

        let elapsed = loc_score_start.elapsed();
        info!("Scoring {} items took: {:?}", num_input_items, elapsed);

        results
    }
}
