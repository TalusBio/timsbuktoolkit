use crate::errors::DataProcessingError;
use crate::{
    IonAnnot,
    QueryItemToScore,
};
use parquet::format::NanoSeconds;
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;
use timsquery::{
    ChromatogramCollector,
    GenerallyQueriable,
    IncludedRange,
    MzMobilityStatsCollector,
    SpectralCollector,
    Tolerance,
};

use super::calculate_scores::{
    IntensityArrays,
    LongitudinalMainScoreElements,
    MainScore,
    PreScore,
    RelativeIntensities,
};
use super::full_results::FullQueryResult;
use super::offsets::MzMobilityOffsets;
use super::search_results::{
    IonSearchResults,
    SearchResultBuilder,
};
use tracing::info;

pub struct Scorer<I: GenerallyQueriable<IonAnnot>> {
    pub index_cycle_rt_ms: Arc<[u32]>,
    pub index: I,
    pub tolerance: Tolerance,
    // The secondsty tolerance is used for ...
    // the secondary query and is meant to be
    // essentually the same but with a narrower retention time
    // range. (instead of the full range)
    pub secondary_tolerance: Tolerance,
    pub fragmented_range: IncludedRange<f64>,
}

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    // does inlining do anything here?
    #[inline]
    fn _build_prescore(&self, item: &QueryItemToScore) -> PreScore {
        let mut agg =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        self.index.add_query(&mut agg, &self.tolerance);

        // TODO: implement an 'is_empty' or equivalent to check if anything was added
        // and if not do a quick exit.

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
        prescore: &PreScore,
        buffer: &mut IntensityArrays,
    ) -> Result<MainScore, DataProcessingError> {
        prescore.localize_with_buffer(buffer)
    }

    #[inline]
    fn _secondary_query(
        &self,
        item: &QueryItemToScore,
        prescore: &PreScore,
        main_score: &MainScore,
    ) -> SpectralCollector<IonAnnot, MzMobilityStatsCollector> {
        // TODO: add the change in the tolerance target rt
        let new_rt = main_score.retention_time_ms;
        let new_query = Arc::new(item.query.with_rt_seconds(new_rt as f32 / 1000.0));
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);
        // TODO use a specific peak width for the secondary query
        self.index.add_query(&mut agg, &self.secondary_tolerance);
        agg
    }

    #[inline]
    fn _finalize_step(
        &self,
        main_score: &MainScore,
        pre_score: &PreScore,
        offsets: &MzMobilityOffsets,
        rel_inten: RelativeIntensities,
    ) -> Result<IonSearchResults, DataProcessingError> {
        let builder = SearchResultBuilder::default();
        builder
            .with_pre_score(pre_score)
            .with_sorted_offsets(offsets)
            .with_relative_intensities(rel_inten)
            .with_main_score(*main_score)
            .finalize()
    }

    #[inline]
    pub fn _build_buffer(
        &self,
        item: &QueryItemToScore,
    ) -> Result<IntensityArrays, DataProcessingError> {
        IntensityArrays::new_empty(
            item.query.precursors.len(),
            item.query.fragments.len(),
            self.index_cycle_rt_ms.clone(),
        )
    }
}

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    /// Scores a single query item by orchestrating the internal steps.
    /// Useful for testing or single-item processing scenarios.
    pub fn buffered_score(
        &self,
        item: QueryItemToScore,
        buffer: &mut IntensityArrays,
    ) -> Option<IonSearchResults> {
        let prescore = self._build_prescore(&item);
        let main_score = match self._localize_step(&prescore, buffer) {
            Ok(score) => score,
            Err(e) => {
                // info!("Error in localization: {:?}", e);
                // info!("Query id: {:#?}", item.query.id);
                return None;
            }
        };
        let secondary_query = self._secondary_query(&item, &prescore, &main_score);
        let offsets = MzMobilityOffsets::new(&secondary_query, item.query.mobility as f64);
        let rel_inten = RelativeIntensities::new(&secondary_query);
        match self._finalize_step(&main_score, &prescore, &offsets, rel_inten) {
            Ok(res) => Some(res),
            Err(e) => {
                info!("Error in scoring: {:?}", e);
                None
            }
        }
    }

    pub fn score(&self, item: QueryItemToScore) -> Option<IonSearchResults> {
        let mut buffer = self._build_buffer(&item).ok()?;
        self.buffered_score(item, &mut buffer)
    }

    pub fn score_full(
        &self,
        item: QueryItemToScore,
    ) -> Result<FullQueryResult, DataProcessingError> {
        let mut res =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        self.index.add_query(&mut res, &self.tolerance);
        let builder = SearchResultBuilder::default();
        let int_arrs = IntensityArrays::new(&res, &item.expected_intensity)?;
        let pre_score = PreScore {
            charge: item.charge,
            digest: item.digest.clone(),
            expected_intensities: item.expected_intensity.clone(),
            query_values: res.clone(),
        };

        let longitudinal_main_score_elements = LongitudinalMainScoreElements::try_new(&int_arrs)?;
        let mut buffer = self._build_buffer(&item)?;
        let main_score = self._localize_step(&pre_score, &mut buffer)?;
        let secondary_query = self._secondary_query(&item, &pre_score, &main_score);
        let offsets = MzMobilityOffsets::new(&secondary_query, item.query.mobility as f64);
        let rel_inten = RelativeIntensities::new(&secondary_query);

        let res2 = builder
            .with_pre_score(&pre_score)
            .with_sorted_offsets(&offsets)
            .with_relative_intensities(rel_inten)
            .with_main_score(main_score)
            .finalize()?;
        let longitudinal_main_score = longitudinal_main_score_elements.main_score_iter().collect();

        Ok(FullQueryResult {
            main_score_elements: longitudinal_main_score_elements,
            longitudinal_main_score,
            extractions: res.clone(),
            search_results: res2,
        })
    }

    pub fn score_iter(&self, items_to_score: &[QueryItemToScore]) -> Vec<IonSearchResults> {
        let num_input_items = items_to_score.len();
        let loc_score_start = Instant::now();

        // There are still A LOT of allocations that can be saved with this pattern
        // but this is a good start.
        let init_fn = || IntensityArrays::new_empty(5, 10, self.index_cycle_rt_ms.clone()).unwrap();

        let results: Vec<IonSearchResults> = items_to_score
            .into_par_iter()
            .filter(|x| {
                let lims = IncludedRange::from(x.query.get_precursor_mz_limits());
                self.fragmented_range.intersects(lims)
            })
            .map_init(init_fn, |buffer, item| {
                // The pre-score is essentually just the collected intensities.
                // TODO: Figure out how to remove this clone ...
                self.buffered_score(item.clone(), buffer)
            })
            .flatten() // Remove the None values (errors)
            .collect();

        let elapsed = loc_score_start.elapsed();
        info!("Scoring {} items took: {:?}", num_input_items, elapsed);

        results
    }
}
