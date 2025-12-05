//! Peptide scoring engine with buffer reuse optimization.
//!
//! # Performance-Critical Design: Buffer Reuse
//!
//! The scoring pipeline processes thousands of peptide queries per second. To achieve this
//! throughput, it's critical to minimize allocations in the hot path.
//!
//! ## Why `buffered_score()` Exists
//!
//! Each scoring operation requires several buffers:
//! - `PeptideScorer` holds time-series feature buffers (size varies with query)
//! - Chromatogram collectors (size varies with query complexity)
//!
//! Creating a new `PeptideScorer` for every query incurs allocation overhead. The solution
//! is **buffer reuse**:
//!
//! 1. `buffered_score()` accepts a mutable `&mut PeptideScorer` reference
//! 2. `score_iter()` uses Rayon's `map_init()` to create one scorer per thread
//! 3. Each thread reuses its scorer across thousands of queries
//!
//! ```rust,,ignore
//! use timsseek::{Scorer, QueryItemToScore};
//! let scorer: Scorer<_> = unimplemented!();
//! let items: Vec<QueryItemToScore> = vec![];
//! // Each thread gets its own reusable PeptideScorer buffer
//! let (results, timings) = scorer.score_iter(&items);
//! ```
//!
//! ## Scoring Pipeline
//!
//! The scoring process has four stages (percentages from benchmark on 225k queries):
//!
//! 1. **Prescore** (27% of time): Extract chromatograms and compute initial intensity scores
//! 2. **Localization** (62% of time): Find peak apex using time-series features (bottleneck)
//! 3. **Secondary Query** (11% of time): Refine search at detected apex with narrow tolerances
//! 4. **Finalization** (<1% of time): Assemble final results with calculated offsets

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
    SpectralCollector,
    Tolerance,
};

use super::accumulator::IonSearchAccumulator;
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
use super::timings::ScoreTimings;
use tracing::{
    debug,
    info,
    warn,
};

/// Default buffer size for precursor ions in PeptideScorer.
/// Sufficient for most peptides with isotopic envelope (M+0, M+1, M+2, M+3, M+4).
const DEFAULT_PRECURSOR_BUFFER_SIZE: usize = 5;

/// Default buffer size for fragment ions in PeptideScorer.
/// Sufficient for typical y/b ion series in most peptides.
const DEFAULT_FRAGMENT_BUFFER_SIZE: usize = 10;

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

/// Peptide scoring engine that searches indexed MS data for spectral library matches.
///
/// The scorer implements a multi-stage pipeline to identify and score peptide candidates:
/// 1. Prescore: Collect chromatographic data across full RT range
/// 2. Localize: Find peak apex using time-series features
/// 3. Secondary Query: Refine search at detected apex with narrow tolerances
/// 4. Finalize: Calculate offsets and assemble results
///
/// # Tolerance Hierarchy
///
/// The scorer uses three tolerance levels for progressive refinement:
///
/// - **`tolerance`**: Initial broad search (full RT range, e.g., ±5 minutes)
///   - Used in prescore phase to collect chromatographic profiles
///   - Captures the entire elution profile including peak tails
///
/// - **`secondary_tolerance`**: Refined search (narrow RT window, e.g., ±30-60 seconds)
///   - Used after localization identifies the peak apex
///   - Focuses on the apex region to reduce noise and improve S/N
///   - Typically has same m/z tolerance but much narrower RT window
///
/// - **Tertiary tolerance** (created inline in `_secondary_query`): Isotope-specific
///   - Derived from `secondary_tolerance` with tighter mobility constraints (5%)
///   - Used for isotope pattern matching where mobility must be very precise
///
/// # Performance Considerations
///
/// This struct is designed for batch processing via `score_iter()`. Creating a new
/// `Scorer` per query is acceptable, but creating new `PeptideScorer` buffers per
/// query adds allocation overhead (exact cost not recently profiled). Use `score_iter()`
/// for production workloads to benefit from buffer reuse.
pub struct Scorer<I: ScorerQueriable> {
    /// Retention time (in milliseconds) for each index cycle.
    /// Shared across all queries to reduce memory overhead.
    pub index_cycle_rt_ms: Arc<[u32]>,

    /// Indexed peak data that implements the required query aggregators.
    pub index: I,

    /// Initial broad search tolerance (full RT range).
    /// Used in prescore phase to collect chromatographic data.
    pub tolerance: Tolerance,

    /// Refined search tolerance (narrow RT window around detected apex).
    /// Used in secondary query phase after localization. Typically has same
    /// m/z tolerance as `tolerance` but narrower RT window (e.g., ±30-60s vs ±5min).
    pub secondary_tolerance: Tolerance,

    /// m/z range where peptides were fragmented.
    /// Used to filter out queries with precursors outside the acquisition window.
    pub fragmented_range: TupleRange<f64>,
}

impl<I: ScorerQueriable> Scorer<I> {
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn _build_prescore(&self, item: &QueryItemToScore) -> PreScore {
        let span = tracing::span!(tracing::Level::TRACE, "build_prescore::new_collector").entered();
        // TODO: pass the collector as a buffer!!!
        // Note: This is surprisingly cheap to do compared to adding the query.
        let mut agg =
            ChromatogramCollector::new(item.query.clone(), self.index_cycle_rt_ms.clone()).unwrap();
        span.exit();

        self.index.add_query(&mut agg, &self.tolerance);

        // Filter out zero-intensity ions and update expected intensities in one pass
        let mut expected_intensities = item.expected_intensity.clone();
        filter_zero_intensity_ions(&mut agg, &mut expected_intensities);

        PreScore {
            charge: item.charge,
            digest: item.digest.clone(),
            expected_intensities,
            query_values: agg,
        }
    }

    /// Calculates the weighted mean ion mobility across fragments and precursors.
    ///
    /// # Algorithm
    ///
    /// Uses log2 weighting to reduce bias from high-abundance ions:
    /// ```text
    /// mobility = Σ(intensity_log2 × mobility) / Σ(intensity_log2)
    /// ```
    ///
    /// # Why Log2 Weighting?
    ///
    /// Linear weighting by intensity can be dominated by one or two very bright ions,
    /// which may have measurement artifacts or be outliers. Log2 weighting:
    /// - Reduces influence of extremely bright ions
    /// - Preserves relative differences in low-abundance signals
    /// - Improves robustness when ion intensity spans orders of magnitude
    ///
    /// For example, with intensities [1000, 100, 10]:
    /// - Linear weights: [1000, 100, 10] → dominated by first ion
    /// - Log2 weights: [9.97, 6.64, 3.32] → more balanced contribution
    ///
    /// # Returns
    ///
    /// Weighted mean mobility in inverse reduced mobility units (1/K0).
    /// Returns 0.0 if no mobility information is available.
    ///
    /// # Panics
    ///
    /// Asserts that mobility is in valid range [0.0, 2.0]. This should always hold
    /// for real MS data; assertion failures indicate data corruption or indexing bugs.
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
    ///
    /// # Algorithm
    ///
    /// After localization identifies the peak apex, this method refines the search
    /// using progressively tighter tolerances:
    ///
    /// **Pass 1: Determine Observed Mobility**
    /// 1. Query at apex RT using `secondary_tolerance` (narrow RT, relaxed mobility)
    /// 2. Calculate weighted mean mobility from observed fragments/precursors
    ///
    /// **Pass 2: Final Refined Query**
    /// 1. Query at apex RT + observed mobility using tight tolerance (3% mobility)
    /// 2. Also query isotope pattern (+1 neutron offset) for ratio calculation
    ///
    /// # Why Two Passes?
    ///
    /// Mobility predictions from spectral libraries can have systematic errors. By first
    /// querying with relaxed mobility to observe the actual signal, then refining with
    /// tight mobility constraints around the observed value, we:
    /// - Reduce noise from off-target ions
    /// - Improve isotope pattern matching
    /// - Correct for mobility calibration drift
    ///
    /// # Returns
    ///
    /// A `SecondaryQuery` containing:
    /// - Refined spectral data with mobility statistics
    /// - Isotope pattern (+1 neutron) for ratio scoring
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, item, main_score), level = "trace")
    )]
    fn _secondary_query(&self, item: &QueryItemToScore, main_score: &MainScore) -> SecondaryQuery {
        let new_rt_seconds = main_score.retention_time_ms as f32 / 1000.0;

        // **Pass 1**: Query at apex RT to determine observed mobility
        let new_query = item.query.clone().with_rt_seconds(new_rt_seconds);
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);
        self.index.add_query(&mut agg, &self.secondary_tolerance);

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

        // TODO: Consider caching this tolerance to avoid repeated allocation
        let tol_use = self
            .secondary_tolerance
            .clone()
            .with_mobility_tolerance(MobilityTolerance::Pct((3.0, 3.0)));

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
        let offsets =
            MzMobilityOffsets::new(&secondary_query.inner, item.query.mobility_ook0() as f64);
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

impl<I: ScorerQueriable> Scorer<I> {
    pub fn buffered_score(
        &self,
        item: QueryItemToScore,
        scorer: &mut PeptideScorer,
        timings: &mut ScoreTimings,
    ) -> Option<IonSearchResults> {
        let st = Instant::now();
        let prescore = self._build_prescore(&item);
        timings.prescore += st.elapsed();

        if prescore
            .expected_intensities
            .fragment_intensities
            .is_empty()
        {
            // No fragments with intensity - cannot score
            return None;
        }

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
            item.query.precursor_count(),
            item.query.fragment_count(),
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
            item.query.precursor_count(),
            item.query.fragment_count(),
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
                DEFAULT_PRECURSOR_BUFFER_SIZE,
                DEFAULT_FRAGMENT_BUFFER_SIZE,
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
        let million_per_min = 1e-6 * throughput * 60.0;
        info!(
            "Scoring {} items took: {:?} throughput: {:#.1}/s, million_per_min: {:#.1}, avg: {:?}",
            num_input_items, elapsed, throughput, million_per_min, avg_speed
        );

        info!("{:?}", results.timings);

        (results.res, results.timings)
    }
}
