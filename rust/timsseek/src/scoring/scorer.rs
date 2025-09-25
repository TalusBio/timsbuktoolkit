use crate::errors::DataProcessingError;
use crate::utils::elution_group_ops::isotope_offset_fragments;
use crate::{
    IonAnnot,
    QueryItemToScore,
};
use rayon::prelude::*;
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    GenerallyQueriable,
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
use super::hyperscore::single_lazyscore;
use super::offsets::MzMobilityOffsets;
use super::search_results::{
    IonSearchResults,
    SearchResultBuilder,
};
use tracing::{
    debug,
    info,
    warn,
};

#[derive(Debug, Clone)]
pub struct SecondaryQuery {
    inner: SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
    isotope: SpectralCollector<IonAnnot, f32>,
}

pub struct SecondaryLazyScores {
    pub lazyscore: f32,
    pub iso_lazyscore: f32,
    pub ratio: f32,
}

impl SecondaryQuery {
    fn lazyscores(&self) -> SecondaryLazyScores {
        let lazyscore = single_lazyscore(
            self.inner
                .iter_fragments()
                .map(|((_k, _mz), v)| v.weight() as f32),
        );
        let iso_lazyscore =
            single_lazyscore(self.isotope.iter_fragments().map(|((_k, _mz), v)| *v));
        let ratio = iso_lazyscore / lazyscore.max(1.0);
        SecondaryLazyScores {
            lazyscore,
            iso_lazyscore,
            ratio,
        }
    }
}

pub struct Scorer<I: GenerallyQueriable<IonAnnot>> {
    pub index_cycle_rt_ms: Arc<[u32]>,
    pub index: I,
    pub tolerance: Tolerance,
    // The secondsty tolerance is used for ...
    // the secondary query and is meant to be
    // essentually the same but with a narrower retention time
    // range. (instead of the full range)
    pub secondary_tolerance: Tolerance,
    pub fragmented_range: TupleRange<f64>,
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
        inten_buffer: &mut IntensityArrays,
        score_buffer: &mut LongitudinalMainScoreElements,
    ) -> Result<MainScore, DataProcessingError> {
        prescore.localize_with_buffer(inten_buffer, score_buffer)
    }

    #[inline]
    fn _secondary_query(
        &self,
        item: &QueryItemToScore,
        _prescore: &PreScore,
        main_score: &MainScore,
    ) -> SecondaryQuery {
        // TODO: add the change in the tolerance target rt
        // proably usig the information in the prescore
        let new_rt = main_score.retention_time_ms;
        let new_query = Arc::new(item.query.as_ref().with_rt_seconds(new_rt as f32 / 1000.0));
        let mut isotope_agg: SpectralCollector<_, f32> =
            SpectralCollector::new(isotope_offset_fragments(&new_query, 1i8).into());
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);

        // TODO use a specific peak width for the secondary query
        self.index.add_query(&mut agg, &self.secondary_tolerance);
        self.index
            .add_query(&mut isotope_agg, &self.secondary_tolerance);
        SecondaryQuery {
            inner: agg,
            isotope: isotope_agg,
        }
    }

    #[inline]
    fn _finalize_step(
        &self,
        item: &QueryItemToScore,
        main_score: &MainScore,
        pre_score: &PreScore,
        secondary_query: &SecondaryQuery,
    ) -> Result<IonSearchResults, DataProcessingError> {
        let offsets = MzMobilityOffsets::new(&secondary_query.inner, item.query.mobility as f64);
        let rel_inten = RelativeIntensities::new(&secondary_query.inner);
        let builder = SearchResultBuilder::default();
        builder
            .with_pre_score(pre_score)
            .with_sorted_offsets(&offsets)
            .with_relative_intensities(rel_inten)
            .with_secondary_lazyscores(secondary_query.lazyscores())
            .with_main_score(*main_score)
            .finalize()
    }

    #[inline]
    pub fn _build_buffers(
        &self,
        item: &QueryItemToScore,
    ) -> Result<(IntensityArrays, LongitudinalMainScoreElements), DataProcessingError> {
        let inten_buffer = IntensityArrays::new_empty(
            item.query.precursors.len(),
            item.query.fragments.len(),
            self.index_cycle_rt_ms.clone(),
        )?;
        let score_buffer =
            LongitudinalMainScoreElements::new_with_capacity(self.index_cycle_rt_ms.len());
        Ok((inten_buffer, score_buffer))
    }
}

#[derive(Debug, Default)]
pub struct ScoreTimings {
    pub prescore: Duration,
    pub localize: Duration,
    pub secondary_query: Duration,
    pub finalization: Duration,
}

impl std::ops::AddAssign for ScoreTimings {
    fn add_assign(&mut self, rhs: Self) {
        self.prescore += rhs.prescore;
        self.localize += rhs.localize;
        self.secondary_query += rhs.secondary_query;
        self.finalization += rhs.finalization;
    }
}

#[derive(Default)]
struct IonSearchAccumulator {
    res: Vec<IonSearchResults>,
    timings: ScoreTimings,
}

impl IonSearchAccumulator {
    fn reduce(mut self, other: Self) -> Self {
        self.res.extend(other.res);
        self.timings += other.timings;
        self
    }

    fn fold(mut self, item: (Option<IonSearchResults>, ScoreTimings)) -> Self {
        if let Some(elem) = item.0 {
            self.res.push(elem);
        }
        self.timings += item.1;
        self
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

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    /// Scores a single query item by orchestrating the internal steps.
    /// Useful for testing or single-item processing scenarios.
    pub fn buffered_score(
        &self,
        item: QueryItemToScore,
        int_buffer: &mut IntensityArrays,
        score_buffer: &mut LongitudinalMainScoreElements,
        timings: &mut ScoreTimings,
    ) -> Option<IonSearchResults> {
        let st = Instant::now();
        let prescore = self._build_prescore(&item);
        timings.prescore += st.elapsed();

        let st = Instant::now();
        let main_score = match self._localize_step(&prescore, int_buffer, score_buffer) {
            Ok(score) => score,
            Err(_e) => {
                // TODO: trim what errors I should log and which are
                // accepted here
                // info!("Error in localization: {:?}", e);
                // info!("Query id: {:#?}", item.query.id);
                return None;
            }
        };
        timings.localize += st.elapsed();

        let st = Instant::now();
        let secondary_query = self._secondary_query(&item, &prescore, &main_score);
        timings.secondary_query += st.elapsed();

        let st = Instant::now();
        let out = self._finalize_step(&item, &main_score, &prescore, &secondary_query);
        timings.finalization += st.elapsed();
        match out {
            Ok(res) => Some(res),
            Err(e) => {
                // I think I can panic here ... since the only
                // possible errors from finalizing are code errors, not recoverable
                // or data ones ...
                // But since this happens in a non-main thread ...
                // not sure if it really matters ...
                warn!("Error in scoring: {:?}", e);
                None
            }
        }
    }

    pub fn score(&self, item: QueryItemToScore) -> Option<IonSearchResults> {
        let mut buffers = self._build_buffers(&item).ok()?;
        let mut timings = ScoreTimings::default();
        let maybe_score = self.buffered_score(item, &mut buffers.0, &mut buffers.1, &mut timings);
        debug!("{:?}", timings);
        maybe_score
    }

    pub fn score_full(
        &self,
        item: QueryItemToScore,
    ) -> Result<FullQueryResult, DataProcessingError> {
        let mut res =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        self.index.add_query(&mut res, &self.tolerance);
        let pre_score = PreScore {
            charge: item.charge,
            digest: item.digest.clone(),
            expected_intensities: item.expected_intensity.clone(),
            query_values: res.clone(),
        };

        // This is a pretty big allocation ..
        let (mut inten_buffer, mut score_buffer) = self._build_buffers(&item)?;
        let main_score = self._localize_step(&pre_score, &mut inten_buffer, &mut score_buffer)?;
        let secondary_query = self._secondary_query(&item, &pre_score, &main_score);
        let res2 = self._finalize_step(&item, &main_score, &pre_score, &secondary_query)?;
        let longitudinal_main_score = score_buffer.main_score_iter().collect();

        Ok(FullQueryResult {
            main_score_elements: score_buffer,
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
        let init_fn = || {
            let int_buffer =
                IntensityArrays::new_empty(5, 10, self.index_cycle_rt_ms.clone()).unwrap();
            let score_buffer =
                LongitudinalMainScoreElements::new_with_capacity(self.index_cycle_rt_ms.len());
            (int_buffer, score_buffer)
        };

        let results: IonSearchAccumulator = items_to_score
            .into_par_iter()
            .with_min_len(512)
            .filter(|x| {
                let tmp = x.query.get_precursor_mz_limits();
                let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should alredy be ordered");
                self.fragmented_range.intersects(lims)
            })
            .map_init(init_fn, |(int_buffer, score_buff), item| {
                // The pre-score is essentually just the collected intensities.
                // TODO: Figure out how to remove this clone ...
                let mut timings = ScoreTimings::default();
                let maybe_score =
                    self.buffered_score(item.clone(), int_buffer, score_buff, &mut timings);
                (maybe_score, timings)
            })
            .collect();

        let elapsed = loc_score_start.elapsed();
        let avg_speed =
            std::time::Duration::from_nanos(elapsed.as_nanos() as u64 / num_input_items as u64);
        let throughput = num_input_items as f64 / elapsed.as_secs_f64();
        info!(
            "Scoring {} items took: {:?} throughput: {}/s, avg: {:?}",
            num_input_items, elapsed, throughput, avg_speed
        );

        info!("{:?}", results.timings);

        results.res
    }
}
