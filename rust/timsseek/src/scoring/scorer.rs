use crate::errors::DataProcessingError;
use crate::utils::elution_group_ops::isotope_offset_fragments;
use crate::{
    IonAnnot,
    QueryItemToScore,
};
use rayon::prelude::*;
use serde::Serialize;
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};
use timsquery::models::tolerance::MobilityTolerance;
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    GenerallyQueriable,
    MzMobilityStatsCollector,
    SpectralCollector,
    Tolerance,
};

use super::calculate_scores::{
    MainScore,
    PeptideScorer,
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
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, item), level = "trace")
    )]
    fn _build_prescore(&self, item: &QueryItemToScore) -> PreScore {
        let span = tracing::span!(tracing::Level::TRACE, "build_prescore::new_collector").entered();
        // TODO: pass the collector as a buffer!!!
        // Note: This is surprisingly cheap to do compared to adding the query.
        let mut agg =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        span.exit();

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

    fn get_mobility(item: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>) -> f64 {
        // Get the weighted mean mobility of the fragments
        let (sum_weights, sum_mobilities) = item
            .iter_fragments()
            .filter_map(|((_k, _mz), v)| {
                v.mean_mobility()
                    .ok()
                    .map(|mobility| (v.weight().log2(), mobility * v.weight().log2()))
            })
            .fold((0.0, 0.0), |(other_w, other_m), (this_w, this_m)| {
                (this_w + other_w, this_m + other_m)
            });

        // Same for precursors
        let (p_sum_weights, p_sum_mobilities) = item
            .iter_precursors()
            .filter_map(|((_k, _mz), v)| {
                v.mean_mobility()
                    .ok()
                    .map(|mobility| (v.weight().log2(), mobility * v.weight().log2()))
            })
            .fold((0.0, 0.0), |(other_w, other_m), (this_w, this_m)| {
                (this_w + other_w, this_m + other_m)
            });

        let total_weight = sum_weights + p_sum_weights;
        let out = if total_weight == 0.0 {
            // No mobility information, return 0
            0.0
        } else {
            (sum_mobilities + p_sum_mobilities) / total_weight
        };
        assert!(out >= 0.0);
        // Mobility should be in that ballpark
        assert!(out <= 2.0);
        out
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, item, main_score), level = "trace")
    )]
    fn _secondary_query(&self, item: &QueryItemToScore, main_score: &MainScore) -> SecondaryQuery {
        // TODO: add the change in the tolerance target rt
        // proably usig the information in the prescore
        let new_rt_seconds = main_score.retention_time_ms as f32 / 1000.0;

        // Make a search here and then use that as a filter for mobility.
        let new_query = item.query.clone().with_rt_seconds(new_rt_seconds);
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);
        self.index.add_query(&mut agg, &self.secondary_tolerance);

        // Get the mobility from the main query
        // and use that to make an isotope queries
        let mobility = Self::get_mobility(&agg);
        let new_query = item
            .query
            .clone()
            .with_rt_seconds(new_rt_seconds)
            .with_mobility(mobility as f32);

        let mut isotope_agg: SpectralCollector<_, f32> =
            SpectralCollector::new(isotope_offset_fragments(&new_query, 1i8));
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);
        // I can probably cache this ...
        // But I dont want to store a "tertiary" tolerance in the scorer
        let tol_use = self
            .secondary_tolerance
            .clone()
            .with_mobility_tolerance(MobilityTolerance::Pct((5.0, 5.0)));

        self.index.add_query(&mut agg, &tol_use);
        self.index.add_query(&mut isotope_agg, &tol_use);
        SecondaryQuery {
            inner: agg,
            isotope: isotope_agg,
        }
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
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
}

#[derive(Debug, Default)]
pub struct ScoreTimings {
    pub prescore: Duration,
    pub localize: Duration,
    pub secondary_query: Duration,
    pub finalization: Duration,
}

impl Serialize for ScoreTimings {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ScoreTimings", 4)?;
        state.serialize_field("prescore_ms", &self.prescore.as_millis())?;
        state.serialize_field("localize_ms", &self.localize.as_millis())?;
        state.serialize_field("secondary_query_ms", &self.secondary_query.as_millis())?;
        state.serialize_field("finalization_ms", &self.finalization.as_millis())?;
        state.end()
    }
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

impl<I: GenerallyQueriable<IonAnnot>> Scorer<I> {
    pub fn buffered_score(
        &self,
        item: QueryItemToScore,
        scorer: &mut PeptideScorer,
        timings: &mut ScoreTimings,
    ) -> Option<IonSearchResults> {
        let st = Instant::now();
        let prescore = self._build_prescore(&item);
        timings.prescore += st.elapsed();

        let st = Instant::now();
        let main_score = match scorer.score(&prescore) {
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
        let secondary_query = self._secondary_query(&item, &main_score);
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
        let mut scorer = PeptideScorer::new(
            item.query.precursors.len(),
            item.query.fragments.len(),
            self.index_cycle_rt_ms.clone(),
        )
        .ok()?;
        let mut timings = ScoreTimings::default();
        let maybe_score = self.buffered_score(item, &mut scorer, &mut timings);
        debug!("{:?}", timings);
        maybe_score
    }

    pub fn score_full(
        &self,
        item: QueryItemToScore,
    ) -> Result<FullQueryResult, DataProcessingError> {
        let mut chromatogram_collector =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        self.index
            .add_query(&mut chromatogram_collector, &self.tolerance);
        let pre_score = PreScore {
            charge: item.charge,
            digest: item.digest.clone(),
            expected_intensities: item.expected_intensity.clone(),
            query_values: chromatogram_collector.clone(),
        };

        let mut scorer = PeptideScorer::new(
            item.query.precursors.len(),
            item.query.fragments.len(),
            self.index_cycle_rt_ms.clone(),
        )?;
        let main_score = scorer.score(&pre_score)?;

        let secondary_query = self._secondary_query(&item, &main_score);
        let search_results =
            self._finalize_step(&item, &main_score, &pre_score, &secondary_query)?;

        let time_resolved_scores = scorer.last_time_resolved_scores();
        let longitudinal_main_score = time_resolved_scores.main_score_iter().collect();

        Ok(FullQueryResult {
            main_score_elements: time_resolved_scores.clone(),
            longitudinal_main_score,
            extractions: chromatogram_collector,
            search_results,
        })
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, items_to_score), level = "trace")
    )]
    pub fn score_iter(
        &self,
        items_to_score: &[QueryItemToScore],
    ) -> (Vec<IonSearchResults>, ScoreTimings) {
        let num_input_items = items_to_score.len();
        let loc_score_start = Instant::now();

        let init_fn = || {
            PeptideScorer::new(
                // We are assuming that 5 precursors and 10 fragments are sufficient
                5,
                10,
                self.index_cycle_rt_ms.clone(),
            )
            .unwrap()
        };

        let filter_fn = |x: &&QueryItemToScore| {
            let tmp = x.query.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };

        #[cfg(not(feature = "serial_scoring"))]
        let results: IonSearchAccumulator = {
            items_to_score
                .into_par_iter()
                .with_min_len(512)
                .filter(filter_fn)
                .map_init(init_fn, |scorer, item| {
                    let mut timings = ScoreTimings::default();
                    let maybe_score = self.buffered_score(item.clone(), scorer, &mut timings);
                    (maybe_score, timings)
                })
                .collect()
        };

        #[cfg(feature = "serial_scoring")]
        let results: IonSearchAccumulator = {
            let mut scorer = init_fn();
            items_to_score
                .into_iter()
                .filter(filter_fn)
                .map(|item| {
                    let mut timings = ScoreTimings::default();
                    let maybe_score = self.buffered_score(item.clone(), &mut scorer, &mut timings);
                    (maybe_score, timings)
                })
                .collect()
        };
        // let results: IonSearchAccumulator = items_to_score
        //     .into_par_iter()
        //     .with_min_len(512)
        //     .filter(|x| {
        //         let tmp = x.query.get_precursor_mz_limits();
        //         let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should alredy be ordered");
        //         self.fragmented_range.intersects(lims)
        //     })
        //     .map_init(init_fn, |scorer, item| {
        //         // The pre-score is essentually just the collected intensities.
        //         // TODO: Figure out how to remove this clone ...
        //         let mut timings = ScoreTimings::default();
        //         let maybe_score = self.buffered_score(item.clone(), scorer, &mut timings);
        //         (maybe_score, timings)
        //     })
        //     .collect();

        let elapsed = loc_score_start.elapsed();
        let avg_speed =
            std::time::Duration::from_nanos(elapsed.as_nanos() as u64 / num_input_items as u64);
        let throughput = num_input_items as f64 / elapsed.as_secs_f64();
        info!(
            "Scoring {} items took: {:?} throughput: {}/s, avg: {:?}",
            num_input_items, elapsed, throughput, avg_speed
        );

        info!("{:?}", results.timings);

        (results.res, results.timings)
    }
}
