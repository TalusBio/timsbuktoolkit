//! Peptide scoring pipeline.
//!
//! # Performance-Critical Design: Buffer Reuse
//!
//! The scoring pipeline processes thousands of peptide queries per second. To achieve this
//! throughput, it's critical to minimize allocations in the hot path.
//!
//! ## Why `process_query()` Exists
//!
//! Each scoring operation requires several buffers:
//! - `ApexFinder` holds time-series feature buffers (size varies with query)
//! - Chromatogram collectors (size varies with query complexity)
//!
//! Creating a new `ApexFinder` for every query incurs allocation overhead. The solution
//! is **buffer reuse**:
//!
//! 1. `process_query()` accepts a mutable `&mut ApexFinder` reference
//! 2. `process_batch()` uses Rayon's `map_init()` to create one buffer per thread
//! 3. Each thread reuses its buffer across thousands of queries
//!
//! ## Scoring Pipeline
//!
//! The scoring process has four stages with separated metadata and scoring data:
//!
//! 1. **Context Build** (27% of time): Extract chromatograms, separate metadata from scoring data
//! 2. **Apex Finding** (62% of time): Find peak apex using time-series features (bottleneck)
//! 3. **Refinement** (11% of time): Two-pass query at detected apex with narrow tolerances
//! 4. **Finalization** (<1% of time): Assemble final results with calculated offsets and metadata

use crate::errors::DataProcessingError;
use crate::utils::elution_group_ops::isotope_offset_fragments;
use crate::{
    IonAnnot,
    QueryItemToScore,
    ScorerQueriable,
};
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;
use timsquery::models::tolerance::MobilityTolerance;
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    MzMobilityStatsCollector,
    OptionallyRestricted,
    SpectralCollector,
    Tolerance,
};

use super::accumulator::IonSearchAccumulator;
use super::apex_finding::{
    ApexFinder,
    ApexScore,
    RelativeIntensities,
};
// use super::calculate_scores::RelativeIntensities; // Removed
use super::full_results::FullQueryResult;
use super::hyperscore::single_lazyscore;
use super::offsets::MzMobilityOffsets;
use super::search_results::{
    IonSearchResults,
    SearchResultBuilder,
};
use super::timings::ScoreTimings;
use tracing::{
    info,
    warn,
};

/// Hierarchical tolerance configuration for the scoring pipeline.
///
/// Progressive refinement: prescore (broad RT) -> secondary (narrow RT) -> tertiary (tight mobility).
#[derive(Debug, Clone)]
pub struct ToleranceHierarchy {
    /// Broad search for prescore phase (full RT range, e.g., +/- 5 min).
    pub prescore: Tolerance,

    /// Refined search at detected apex (narrow RT window, e.g., +/- 30-60s).
    pub secondary: Tolerance,
}

impl ToleranceHierarchy {
    /// Creates tertiary tolerance with 3% mobility constraints for isotope matching.
    pub fn tertiary_tolerance(&self) -> Tolerance {
        self.secondary
            .clone()
            .with_mobility_tolerance(MobilityTolerance::Pct((3.0, 3.0)))
    }
}

/// Filter out zero-intensity ions and update expected intensities in one pass.
///
/// This maintains index alignment by removing ions from the chromatogram collector
/// and expected intensities simultaneously in a single loop.
#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn filter_zero_intensity_ions(
    agg: &mut ChromatogramCollector<IonAnnot, f32>,
    expected: &mut crate::ExpectedIntensities,
) {
    // Early-exit predicate: stop at first non-zero value (much faster than summing)
    let predicate = |chrom: &[f32]| chrom.iter().any(|&x| x > 0.0);

    // Filter precursors
    for (k, _mz) in agg.drain_nonmatching_precursors(predicate) {
        expected.precursor_intensities.remove(&k);
    }

    // Filter fragments
    for (k, _mz) in agg.drain_nonmatching_fragments(predicate) {
        expected.fragment_intensities.remove(&k);
    }

    // Assert all lengths match
    assert_eq!(
        agg.precursors.num_ions(),
        agg.eg.precursor_count(),
        "Precursor count mismatch after filtering"
    );
    assert_eq!(
        agg.precursors.num_ions(),
        expected.precursor_intensities.len(),
        "Precursor expected intensities count mismatch after filtering"
    );
    assert_eq!(
        agg.fragments.num_ions(),
        agg.eg.fragment_count(),
        "Fragment count mismatch after filtering"
    );
}

#[derive(Debug, Clone, Copy)]
pub struct SecondaryLazyScores {
    pub lazyscore: f32,
    pub iso_lazyscore: f32,
    pub ratio: f32,
}

/// Compute lazyscores from inner and isotope collectors.
fn compute_secondary_lazyscores(
    inner: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
    isotope: &SpectralCollector<IonAnnot, f32>,
) -> SecondaryLazyScores {
    let lazyscore = single_lazyscore(
        inner
            .iter_fragments()
            .map(|((_k, _mz), v)| v.weight() as f32),
    );
    let iso_lazyscore =
        single_lazyscore(isotope.iter_fragments().map(|((_k, _mz), v)| *v));
    let ratio = iso_lazyscore / lazyscore.max(1.0);
    SecondaryLazyScores {
        lazyscore,
        iso_lazyscore,
        ratio,
    }
}

/// Multi-stage pipeline for scoring peptide queries against indexed MS data.
///
/// Pipeline stages: build context → find apex → refine → finalize.
/// Uses progressive tolerance refinement, metadata separation, and buffer reuse for high throughput.
pub struct ScoringPipeline<I: ScorerQueriable> {
    /// Retention time in milliseconds for each index cycle.
    pub index_cycle_rt_ms: Arc<[u32]>,

    /// Indexed peak data that implements the required query aggregators.
    pub index: I,

    /// Hierarchical tolerance configuration for progressive refinement.
    pub tolerances: ToleranceHierarchy,

    /// m/z range where peptides were fragmented.
    /// Queries with precursors outside this range are filtered out.
    pub fragmented_range: TupleRange<f64>,
}

enum SkippingReason {
    // TODO: Implement more options and a counter ...
    RetentionTimeOutOfBounds,
}

impl<I: ScorerQueriable> ScoringPipeline<I> {
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn build_candidate_context(
        &self,
        item: &QueryItemToScore,
    ) -> Result<(super::apex_finding::PeptideMetadata, super::apex_finding::ScoringContext), SkippingReason> {
        let max_range = TupleRange::try_new(
            *self.index_cycle_rt_ms.first().unwrap(),
            *self.index_cycle_rt_ms.last().unwrap(),
        )
        .expect("Reference RTs should be sorted and valid");
        let rt_range = match self
            .tolerances
            .prescore
            .rt_range_as_milis(item.query.rt_seconds())
        {
            OptionallyRestricted::Unrestricted => max_range,
            OptionallyRestricted::Restricted(r) => r,
        };

        if !max_range.intersects(rt_range) {
            return Err(SkippingReason::RetentionTimeOutOfBounds);
        }
        let mut agg = tracing::span!(
            tracing::Level::TRACE,
            "build_candidate_context::new_collector"
        ).in_scope(|| {
            match ChromatogramCollector::new(
                item.query.clone(),
                rt_range,
                &self.index_cycle_rt_ms,
            ) {
                Ok(collector) => collector,
                Err(e) => {
                    let tol_range = self.tolerances.prescore.rt_range_as_milis(item.query.rt_seconds());
                    panic!(
                        "Failed to create ChromatogramCollector for query id {:#?}: {:?} with RT tolerance {:#?}",
                        item.query, e, tol_range,
                    )
                }
            }
        });

        tracing::span!(tracing::Level::TRACE, "build_candidate_context::add_query").in_scope(|| {
            self.index.add_query(&mut agg, &self.tolerances.prescore);
        });

        // Filter out zero-intensity ions and update expected intensities in one pass
        let mut expected_intensities = item.expected_intensity.clone();
        filter_zero_intensity_ions(&mut agg, &mut expected_intensities);

        let metadata = super::apex_finding::PeptideMetadata {
            digest: item.digest.clone(),
            charge: item.query.precursor_charge(),
            library_id: agg.eg.id() as u32,
            ref_rt_seconds: item.query.rt_seconds(),
            ref_mobility_ook0: item.query.mobility_ook0(),
            ref_precursor_mz: item.query.mono_precursor_mz(),
        };

        let scoring_ctx = super::apex_finding::ScoringContext {
            expected_intensities,
            query_values: agg,
        };

        Ok((metadata, scoring_ctx))
    }

    /// Calculates the weighted mean ion mobility across fragments and precursors.
    fn get_mobility(item: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>) -> f64 {
        // Calculate weighted mean mobility from fragments
        let (sum_weights, sum_mobilities) = item
            .iter_fragments()
            .filter_map(|((_k, _mz), v)| {
                v.mean_mobility().ok().map(|mobility| {
                    let weight = v.weight().log2();
                    (weight, mobility * weight)
                })
            })
            .fold((0.0, 0.0), |(other_w, other_m), (this_w, this_m)| {
                (this_w + other_w, this_m + other_m)
            });

        // Same calculation for precursors
        let (p_sum_weights, p_sum_mobilities) = item
            .iter_precursors()
            .filter_map(|((_k, _mz), v)| {
                v.mean_mobility().ok().map(|mobility| {
                    let weight = v.weight().log2();
                    (weight, mobility * weight)
                })
            })
            .fold((0.0, 0.0), |(other_w, other_m), (this_w, this_m)| {
                (this_w + other_w, this_m + other_m)
            });

        let total_weight = sum_weights + p_sum_weights;
        let out = if total_weight == 0.0 {
            // No mobility information available
            0.0
        } else {
            (sum_mobilities + p_sum_mobilities) / total_weight
        };

        assert!(out >= 0.0);
        assert!(out <= 2.0, "Mobility {} outside expected range [0, 2]", out);
        out
    }

    /// Performs refined secondary query at detected apex with two-pass strategy.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, item, main_score), level = "trace")
    )]
    fn execute_secondary_query(
        &self,
        item: &QueryItemToScore,
        main_score: &ApexScore,
    ) -> (
        SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
        SpectralCollector<IonAnnot, f32>,
    ) {
        let new_rt_seconds = main_score.retention_time_ms as f32 / 1000.0;

        // **Pass 1**: Query at apex RT to determine observed mobility
        let new_query = item.query.clone().with_rt_seconds(new_rt_seconds);
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);
        self.index.add_query(&mut agg, &self.tolerances.secondary);

        // Calculate weighted mean mobility from observed data
        let mobility = Self::get_mobility(&agg);

        // **Pass 2**: Query at apex RT + observed mobility with tight tolerance
        let new_query = item
            .query
            .clone()
            .with_rt_seconds(new_rt_seconds)
            .with_mobility(mobility as f32);

        // Query isotope pattern (+1 neutron offset) for ratio scoring
        let mut isotope_agg: SpectralCollector<_, f32> =
            SpectralCollector::new(isotope_offset_fragments(&new_query, 1i8));

        // Query main pattern with tight mobility tolerance
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);

        let tol_use = self.tolerances.tertiary_tolerance();

        self.index.add_query(&mut agg, &tol_use);
        self.index.add_query(&mut isotope_agg, &tol_use);

        (agg, isotope_agg)
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn finalize_results(
        &self,
        metadata: &super::apex_finding::PeptideMetadata,
        nqueries: u8,
        main_score: &ApexScore,
        inner_collector: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
        isotope_collector: &SpectralCollector<IonAnnot, f32>,
    ) -> Result<IonSearchResults, DataProcessingError> {
        let offsets = MzMobilityOffsets::new(inner_collector, metadata.ref_mobility_ook0 as f64);
        let rel_inten = RelativeIntensities::new(inner_collector);
        let lazyscores = compute_secondary_lazyscores(inner_collector, isotope_collector);

        let builder = SearchResultBuilder::default();

        builder
            .with_metadata(metadata)
            .with_nqueries(nqueries)
            .with_sorted_offsets(&offsets)
            .with_relative_intensities(rel_inten)
            .with_secondary_lazyscores(lazyscores)
            .with_apex_score(main_score)
            .finalize()
    }
}

impl<I: ScorerQueriable> ScoringPipeline<I> {
    pub fn process_query(
        &self,
        item: QueryItemToScore,
        buffer: &mut ApexFinder,
        timings: &mut ScoreTimings,
    ) -> Option<IonSearchResults> {
        let st = Instant::now();
        let (metadata, scoring_ctx) = match self.build_candidate_context(&item) {
            Ok(result) => result,
            Err(SkippingReason::RetentionTimeOutOfBounds) => {
                return None;
            }
        };
        timings.prescore += st.elapsed();

        if scoring_ctx.expected_intensities.fragment_intensities.is_empty() {
            return None;
        }

        let st = Instant::now();
        let apex_score = match buffer.find_apex(&scoring_ctx) {
            Ok(score) => score,
            Err(_e) => {
                return None;
            }
        };
        timings.localize += st.elapsed();

        let st = Instant::now();
        let (inner_collector, isotope_collector) =
            self.execute_secondary_query(&item, &apex_score);
        timings.secondary_query += st.elapsed();

        let nqueries = scoring_ctx.query_values.fragments.num_ions() as u8;

        let st = Instant::now();
        let out = self.finalize_results(
            &metadata,
            nqueries,
            &apex_score,
            &inner_collector,
            &isotope_collector,
        );
        timings.finalization += st.elapsed();

        match out {
            Ok(res) => Some(res),
            Err(e) => {
                warn!("Error in scoring: {:?}", e);
                None
            }
        }
    }

    pub fn process_query_full(
        &self,
        item: QueryItemToScore,
    ) -> Result<FullQueryResult, DataProcessingError> {
        let mut buffer = ApexFinder::new(self.index_cycle_rt_ms.clone());

        // Re-implementing logic here because process_query consumes `item` and returns `Option`.
        // We want intermediate results for `FullQueryResult`.

        let (metadata, scoring_ctx) = self
            .build_candidate_context(&item)
            .map_err(|_| DataProcessingError::ExpectedNonEmptyData {
                context: Some("RT out of bounds".into()),
            })?;

        let apex_score = buffer.find_apex(&scoring_ctx)?;
        let (inner_collector, isotope_collector) = self.execute_secondary_query(&item, &apex_score);

        let nqueries = scoring_ctx.query_values.fragments.num_ions() as u8;
        let search_results = self.finalize_results(
            &metadata,
            nqueries,
            &apex_score,
            &inner_collector,
            &isotope_collector,
        )?;

        // Extract query_values before it's consumed
        let extractions = scoring_ctx.query_values;

        Ok(FullQueryResult {
            main_score_elements: buffer.traces.clone(),
            longitudinal_main_score: buffer.traces.main_score.clone(),
            extractions,
            search_results,
        })
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, items_to_score), level = "trace")
    )]
    pub fn process_batch(
        &self,
        items_to_score: &[QueryItemToScore],
    ) -> (Vec<IonSearchResults>, ScoreTimings) {
        let num_input_items = items_to_score.len();
        let loc_score_start = Instant::now();

        let init_fn = || {
            ApexFinder::new(self.index_cycle_rt_ms.clone())
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
                    let maybe_score = self.process_query(item.clone(), scorer, &mut timings);
                    (maybe_score, timings)
                })
                .collect()
        };

        #[cfg(feature = "serial_scoring")]
        let results: IonSearchAccumulator = {
            let mut scorer = init_fn();
            items_to_score
                .iter()
                .filter(filter_fn)
                .map(|item| {
                    let mut timings = ScoreTimings::default();
                    let maybe_score = self.process_query(item.clone(), &mut scorer, &mut timings);
                    (maybe_score, timings)
                })
                .collect()
        };

        let elapsed = loc_score_start.elapsed();
        let avg_speed =
            std::time::Duration::from_nanos(elapsed.as_nanos() as u64 / num_input_items as u64);
        let throughput = num_input_items as f64 / elapsed.as_secs_f64();
        let million_per_min = 1e-6 * throughput * 60.0;
        info!(
            "Scoring {} items took: {:?} throughput: {:#.1}/s, million_per_min: {:#.1}, avg: {:?}",
            num_input_items, elapsed, throughput, million_per_min, avg_speed
        );

        info!("{:?}", results.timings);

        (results.res, results.timings)
    }
}
