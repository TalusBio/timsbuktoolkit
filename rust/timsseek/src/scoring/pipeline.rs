//! Peptide scoring pipeline.
//!
//! # Performance-Critical Design: Buffer Reuse
//!
//! The scoring pipeline processes thousands of peptide queries per second. To achieve this
//! throughput, it's critical to minimize allocations in the hot path.
//!
//! Each scoring operation requires several buffers:
//! - `ApexFinder` holds time-series feature buffers (size varies with query)
//! - Chromatogram collectors (size varies with query complexity)
//!
//! `prescore_batch` and `score_calibrated_batch` use Rayon's `map_init()` / `fold` to create
//! one `ApexFinder` buffer per thread, which is reused across thousands of queries.
//!
//! ## Scoring Pipeline
//!
//! The scoring process has two phases:
//!
//! 1. **Prescore** (Phase 1): Broad extraction + `find_apex_location` — yields calibrant candidates.
//! 2. **Calibrated scoring** (Phase 3): Narrow calibrated extraction + `find_apex` + secondary query.

use crate::errors::DataProcessingError;
use crate::utils::elution_group_ops::isotope_offset_fragments;
use crate::{
    IonAnnot,
    QueryItemToScore,
    ScorerQueriable,
};
use rayon::prelude::*;
use std::time::Instant;
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    KeyLike,
    MzMobilityStatsCollector,
    OptionallyRestricted,
    SpectralCollector,
    Tolerance,
};

use super::accumulator::IonSearchAccumulator;
use super::apex_finding::{
    ApexFinder,
    ApexLocation,
    ApexScore,
    PeptideMetadata,
    RelativeIntensities,
};
use super::full_results::FullQueryResult;
use super::hyperscore::single_lazyscore;
use super::offsets::MzMobilityOffsets;
use super::results::{
    ScoredCandidate,
    ScoredCandidateBuilder,
};
use super::timings::ScoreTimings;
use crate::rt_calibration::CalibrationResult;
use tracing::warn;

/// Lightweight calibrant candidate — just enough to re-query in Phase 2.
/// Implements Ord by score (ascending) for use in BinaryHeap<Reverse<_>>.
#[derive(Debug, Clone)]
pub struct CalibrantCandidate {
    pub score: f32,
    pub apex_rt_seconds: f32,
    pub speclib_index: usize,
}

impl PartialEq for CalibrantCandidate {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
    }
}

impl Eq for CalibrantCandidate {}

impl PartialOrd for CalibrantCandidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.score.partial_cmp(&other.score)
    }
}

impl Ord for CalibrantCandidate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Bounded min-heap: keeps only the top-N candidates by score.
/// Uses Reverse so the *smallest* score is at the top (ejected first).
pub struct CalibrantHeap {
    heap: std::collections::BinaryHeap<std::cmp::Reverse<CalibrantCandidate>>,
    capacity: usize,
}

impl CalibrantHeap {
    pub fn new(capacity: usize) -> Self {
        Self {
            heap: std::collections::BinaryHeap::with_capacity(capacity + 1),
            capacity,
        }
    }

    pub fn push(&mut self, candidate: CalibrantCandidate) {
        if !candidate.score.is_finite() || candidate.score <= 0.0 {
            return;
        }
        if self.heap.len() < self.capacity {
            self.heap.push(std::cmp::Reverse(candidate));
        } else if let Some(std::cmp::Reverse(min)) = self.heap.peek() {
            if candidate.score > min.score {
                self.heap.pop();
                self.heap.push(std::cmp::Reverse(candidate));
            }
        }
    }

    pub fn merge(mut self, other: Self) -> Self {
        for item in other.heap {
            self.push(item.0);
        }
        self
    }

    pub fn into_vec(self) -> Vec<CalibrantCandidate> {
        self.heap.into_iter().map(|std::cmp::Reverse(c)| c).collect()
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }
}

/// Calibration configuration — all tunable parameters with defaults.
#[derive(Debug, Clone)]
pub struct CalibrationConfig {
    pub n_calibrants: usize,
    pub grid_size: usize,
    pub mz_sigma: f32,
    pub mobility_sigma: f32,
    pub rt_sigma_factor: f32,
    pub min_rt_tolerance_minutes: f32,
    pub lowess_frac: f32,
    pub calibration_query_rt_window_minutes: f32,
    pub dp_lookback: usize,
}

impl Default for CalibrationConfig {
    fn default() -> Self {
        Self {
            n_calibrants: 2000,
            grid_size: 100,
            mz_sigma: 2.0,
            mobility_sigma: 3.0,
            rt_sigma_factor: 3.0,
            min_rt_tolerance_minutes: 0.5,
            lowess_frac: 0.5,
            calibration_query_rt_window_minutes: 0.5,
            dp_lookback: 30,
        }
    }
}

/// Number of top fragments to retain for scoring (by predicted intensity).
const TOP_N_FRAGMENTS: usize = 8;

/// Retain only the top `n` fragments by predicted intensity.
///
/// Removes lower-ranked fragments from the chromatogram collector (fragments array + eg)
/// and from expected intensities, maintaining the invariant that all three agree on count.
fn select_top_n_fragments<T: KeyLike>(
    agg: &mut ChromatogramCollector<T, f32>,
    expected: &mut crate::ExpectedIntensities<T>,
    n: usize,
) {
    let n_frags = agg.fragments.num_ions();
    if n_frags <= n {
        return;
    }

    // Build (positional_index, key, expected_intensity) for all fragments
    let mut indexed: Vec<(usize, T, f32)> = agg
        .fragments
        .iter_mzs()
        .enumerate()
        .map(|(idx, ((key, _mz), _chrom))| {
            let intensity = expected
                .fragment_intensities
                .get(key)
                .copied()
                .unwrap_or(0.0);
            (idx, key.clone(), intensity)
        })
        .collect();

    // Sort descending by expected intensity
    indexed.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    // Collect indices and keys to drop (everything beyond top N)
    let mut to_drop: Vec<(usize, T)> = indexed
        .into_iter()
        .skip(n)
        .map(|(idx, key, _)| (idx, key))
        .collect();

    // Sort drop-indices descending so highest index is removed first (avoids shift)
    to_drop.sort_by(|a, b| b.0.cmp(&a.0));

    for (idx, key) in to_drop {
        agg.fragments
            .drop_row_idx(idx)
            .expect("index should be in bounds");
        agg.eg
            .try_drop_fragment(&key)
            .expect("key should exist in eg");
        expected.fragment_intensities.remove(&key);
    }

    debug_assert_eq!(agg.fragments.num_ions(), agg.eg.fragment_count());
    debug_assert_eq!(
        agg.fragments.num_ions(),
        expected.fragment_intensities.len()
    );
}

/// Filter out zero-intensity ions and update expected intensities in one pass.
///
/// This maintains index alignment by removing ions from the chromatogram collector
/// and expected intensities simultaneously in a single loop.
#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
fn filter_zero_intensity_ions<T: KeyLike>(
    agg: &mut ChromatogramCollector<T, f32>,
    expected: &mut crate::ExpectedIntensities<T>,
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
    let iso_lazyscore = single_lazyscore(isotope.iter_fragments().map(|((_k, _mz), v)| *v));
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
pub struct Scorer<I: ScorerQueriable> {
    /// Indexed peak data that implements the required query aggregators.
    pub index: I,

    /// Broad tolerance used during the prescore phase.
    pub broad_tolerance: Tolerance,

    /// m/z range where peptides were fragmented.
    /// Queries with precursors outside this range are filtered out.
    pub fragmented_range: TupleRange<f64>,
}

enum SkippingReason {
    // TODO: Implement more options and a counter ...
    RetentionTimeOutOfBounds,
}

impl<I: ScorerQueriable> Scorer<I> {
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn build_broad_extraction(
        &self,
        item: &QueryItemToScore,
    ) -> Result<
        (
            super::apex_finding::PeptideMetadata,
            super::apex_finding::Extraction<IonAnnot>,
        ),
        SkippingReason,
    > {
        let max_range = self.index.ms1_cycle_mapping().range_milis();
        let max_range = TupleRange::try_new(max_range.0, max_range.1)
            .expect("Reference RTs should be sorted and valid");
        let rt_range = match self
            .broad_tolerance
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
            "build_broad_extraction::new_collector"
        ).in_scope(|| {
            match ChromatogramCollector::new(
                item.query.clone(),
                rt_range,
                self.index.ms1_cycle_mapping(),
            ) {
                Ok(collector) => collector,
                Err(e) => {
                    let tol_range = self.broad_tolerance.rt_range_as_milis(item.query.rt_seconds());
                    panic!(
                        "Failed to create ChromatogramCollector for query id {:#?}: {:?} with RT tolerance {:#?}",
                        item.query, e, tol_range,
                    )
                }
            }
        });

        tracing::span!(tracing::Level::TRACE, "build_broad_extraction::add_query").in_scope(
            || {
                self.index.add_query(&mut agg, &self.broad_tolerance);
            },
        );

        // Filter out zero-intensity ions and update expected intensities in one pass
        let mut expected_intensities = item.expected_intensity.clone();
        filter_zero_intensity_ions(&mut agg, &mut expected_intensities);

        // Retain only top-N fragments by predicted intensity for scoring
        select_top_n_fragments(&mut agg, &mut expected_intensities, TOP_N_FRAGMENTS);

        let metadata = super::apex_finding::PeptideMetadata {
            digest: item.digest.clone(),
            charge: item.query.precursor_charge(),
            library_id: agg.eg.id() as u32,
            query_rt_seconds: item.query.rt_seconds(),
            ref_mobility_ook0: item.query.mobility_ook0(),
            ref_precursor_mz: item.query.mono_precursor_mz(),
        };

        let scoring_ctx = super::apex_finding::Extraction {
            expected_intensities,
            chromatograms: agg,
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
        tracing::instrument(skip(self, item, main_score, spectral_tol, isotope_tol), level = "trace")
    )]
    fn execute_secondary_query(
        &self,
        item: &QueryItemToScore,
        main_score: &ApexScore,
        spectral_tol: &Tolerance,
        isotope_tol: &Tolerance,
    ) -> (
        SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
        SpectralCollector<IonAnnot, f32>,
    ) {
        let new_rt_seconds = main_score.retention_time_ms as f32 / 1000.0;

        // **Pass 1**: Query at apex RT to determine observed mobility
        let new_query = item.query.clone().with_rt_seconds(new_rt_seconds);
        let mut agg: SpectralCollector<_, MzMobilityStatsCollector> =
            SpectralCollector::new(new_query);
        self.index.add_query(&mut agg, spectral_tol);

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

        self.index.add_query(&mut agg, isotope_tol);
        self.index.add_query(&mut isotope_agg, isotope_tol);

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
    ) -> Result<ScoredCandidate, DataProcessingError> {
        let offsets = MzMobilityOffsets::new(inner_collector, metadata.ref_mobility_ook0 as f64);
        let rel_inten = RelativeIntensities::new(inner_collector);
        let lazyscores = compute_secondary_lazyscores(inner_collector, isotope_collector);

        ScoredCandidateBuilder::default()
            .with_metadata(metadata)
            .with_nqueries(nqueries)
            .with_sorted_offsets(&offsets)
            .with_relative_intensities(rel_inten)
            .with_secondary_lazyscores(lazyscores)
            .with_apex_score(main_score)
            .finalize()
    }
}

impl<I: ScorerQueriable> Scorer<I> {
    pub fn process_query_full(
        &self,
        item: QueryItemToScore,
        calibration: &CalibrationResult,
    ) -> Result<FullQueryResult, DataProcessingError> {
        let mut buffer = ApexFinder::new(self.num_cycles());

        // Re-implementing logic here because process_query consumes `item` and returns `Option`.
        // We want intermediate results for `FullQueryResult`.

        let (metadata, scoring_ctx) = self.build_broad_extraction(&item).map_err(|_| {
            DataProcessingError::ExpectedNonEmptyData {
                context: Some("RT out of bounds".into()),
            }
        })?;

        let apex_score = buffer.find_apex(&scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))?;
        let spectral_tol = calibration.get_spectral_tolerance();
        let isotope_tol = calibration.get_isotope_tolerance();
        let (inner_collector, isotope_collector) =
            self.execute_secondary_query(&item, &apex_score, &spectral_tol, &isotope_tol);

        let nqueries = scoring_ctx.chromatograms.fragments.num_ions() as u8;
        let search_results = self.finalize_results(
            &metadata,
            nqueries,
            &apex_score,
            &inner_collector,
            &isotope_collector,
        )?;

        // Extract chromatograms before it's consumed
        let extractions = scoring_ctx.chromatograms;

        Ok(FullQueryResult {
            main_score_elements: buffer.traces.clone(),
            longitudinal_main_score: buffer.traces.apex_profile.clone(),
            extractions,
            search_results,
        })
    }

    /// Build a chromatogram extraction using calibrated RT and per-query tolerance.
    /// The speclib is NOT mutated — CalibrationResult provides the RT conversion.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn build_calibrated_extraction(
        &self,
        item: &QueryItemToScore,
        calibration: &CalibrationResult,
    ) -> Result<
        (
            super::apex_finding::PeptideMetadata,
            super::apex_finding::Extraction<IonAnnot>,
        ),
        SkippingReason,
    > {
        let original_irt = item.query.rt_seconds();
        let calibrated_rt = calibration.convert_irt(original_irt);
        let tolerance = calibration.get_tolerance(
            item.query.mono_precursor_mz(),
            item.query.mobility_ook0(),
            calibrated_rt,
        );

        let calibrated_query = item.query.clone().with_rt_seconds(calibrated_rt);

        let max_range = self.index.ms1_cycle_mapping().range_milis();
        let max_range = TupleRange::try_new(max_range.0, max_range.1)
            .expect("Reference RTs should be sorted and valid");
        let rt_range = match tolerance.rt_range_as_milis(calibrated_rt) {
            OptionallyRestricted::Unrestricted => max_range,
            OptionallyRestricted::Restricted(r) => r,
        };

        if !max_range.intersects(rt_range) {
            return Err(SkippingReason::RetentionTimeOutOfBounds);
        }

        let mut agg = ChromatogramCollector::new(
            calibrated_query,
            rt_range,
            self.index.ms1_cycle_mapping(),
        )
        .map_err(|_| SkippingReason::RetentionTimeOutOfBounds)?;

        self.index.add_query(&mut agg, &tolerance);

        let mut expected_intensities = item.expected_intensity.clone();
        filter_zero_intensity_ions(&mut agg, &mut expected_intensities);
        select_top_n_fragments(&mut agg, &mut expected_intensities, TOP_N_FRAGMENTS);

        let metadata = super::apex_finding::PeptideMetadata {
            digest: item.digest.clone(),
            charge: item.query.precursor_charge(),
            library_id: agg.eg.id() as u32,
            query_rt_seconds: calibrated_rt,
            ref_mobility_ook0: item.query.mobility_ook0(),
            ref_precursor_mz: item.query.mono_precursor_mz(),
        };

        let scoring_ctx = super::apex_finding::Extraction {
            expected_intensities,
            chromatograms: agg,
        };

        Ok((metadata, scoring_ctx))
    }

    /// Phase 3: Score a peptide using calibrated extraction window.
    /// Expects narrow calibrated extraction (from CalibrationResult).
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn score_calibrated_extraction(
        &self,
        item: &QueryItemToScore,
        calibration: &CalibrationResult,
        buffer: &mut ApexFinder,
        timings: &mut ScoreTimings,
    ) -> Option<ScoredCandidate> {
        let st = Instant::now();
        let (metadata, scoring_ctx) =
            tracing::span!(tracing::Level::TRACE, "score_calibrated::extraction").in_scope(
                || match self.build_calibrated_extraction(item, calibration) {
                    Ok(result) => Some(result),
                    Err(_) => None,
                },
            )?;
        timings.prescore += st.elapsed();

        if scoring_ctx
            .expected_intensities
            .fragment_intensities
            .is_empty()
        {
            return None;
        }

        let st = Instant::now();
        let apex_score =
            tracing::span!(tracing::Level::TRACE, "score_calibrated::apex_scoring").in_scope(
                || {
                    buffer
                        .find_apex(&scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                        .ok()
                },
            )?;
        timings.localize += st.elapsed();

        let st = Instant::now();
        let spectral_tol = calibration.get_spectral_tolerance();
        let isotope_tol = calibration.get_isotope_tolerance();
        let (inner_collector, isotope_collector) =
            tracing::span!(tracing::Level::TRACE, "score_calibrated::secondary_query")
                .in_scope(|| self.execute_secondary_query(item, &apex_score, &spectral_tol, &isotope_tol));
        timings.secondary_query += st.elapsed();

        let nqueries = scoring_ctx.chromatograms.fragments.num_ions() as u8;
        let st = Instant::now();
        let out = tracing::span!(tracing::Level::TRACE, "score_calibrated::finalize").in_scope(
            || {
                self.finalize_results(
                    &metadata,
                    nqueries,
                    &apex_score,
                    &inner_collector,
                    &isotope_collector,
                )
            },
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

    /// Phase 3 batch: Score all peptides with calibrated tolerances.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn score_calibrated_batch(
        &self,
        items_to_score: &[QueryItemToScore],
        calibration: &CalibrationResult,
    ) -> (Vec<ScoredCandidate>, ScoreTimings) {
        let init_fn = || ApexFinder::new(self.num_cycles());
        let filter_fn = |x: &&QueryItemToScore| {
            let tmp = x.query.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };

        #[cfg(not(feature = "serial_scoring"))]
        let results: IonSearchAccumulator = {
            items_to_score
                .par_iter()
                .filter(filter_fn)
                .map_init(init_fn, |scorer, item| {
                    let mut t = ScoreTimings::default();
                    let result =
                        self.score_calibrated_extraction(item, calibration, scorer, &mut t);
                    (result, t)
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
                    let mut t = ScoreTimings::default();
                    let result =
                        self.score_calibrated_extraction(item, calibration, &mut scorer, &mut t);
                    (result, t)
                })
                .collect()
        };

        (results.res, results.timings)
    }

    /// Phase 1: Lightweight prescore — broad extraction + find_apex_location only.
    /// Returns the apex location (with split product score) and metadata.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn prescore(
        &self,
        item: &QueryItemToScore,
        buffer: &mut ApexFinder,
    ) -> Option<(ApexLocation, PeptideMetadata)> {
        let (metadata, scoring_ctx) = tracing::span!(tracing::Level::TRACE, "prescore::extraction")
            .in_scope(|| match self.build_broad_extraction(item) {
                Ok(result) => Some(result),
                Err(SkippingReason::RetentionTimeOutOfBounds) => None,
            })?;

        if scoring_ctx
            .expected_intensities
            .fragment_intensities
            .is_empty()
        {
            return None;
        }

        let apex_location =
            tracing::span!(tracing::Level::TRACE, "prescore::scoring").in_scope(|| {
                buffer
                    .find_apex_location(&scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                    .ok()
            })?;

        Some((apex_location, metadata))
    }

    /// Phase 1 batch: Prescore all peptides, collecting top-N calibrant candidates via bounded heaps.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn prescore_batch(
        &self,
        items_to_score: &[QueryItemToScore],
        speclib_offset: usize,
        config: &CalibrationConfig,
    ) -> CalibrantHeap {
        let filter_fn = |x: &&QueryItemToScore| {
            let tmp = x.query.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };

        #[cfg(not(feature = "serial_scoring"))]
        let heap: CalibrantHeap = {
            let init_fn =
                || (ApexFinder::new(self.num_cycles()), CalibrantHeap::new(config.n_calibrants));

            items_to_score
                .par_iter()
                .enumerate()
                .filter(|(_, x)| filter_fn(x))
                .fold(init_fn, |(mut scorer, mut heap), (chunk_idx, item)| {
                    if let Some((loc, _meta)) = self.prescore(item, &mut scorer) {
                        heap.push(CalibrantCandidate {
                            score: loc.score,
                            apex_rt_seconds: loc.retention_time_ms as f32 / 1000.0,
                            speclib_index: speclib_offset + chunk_idx,
                        });
                    }
                    (scorer, heap)
                })
                .map(|(_, heap)| heap)
                .reduce(
                    || CalibrantHeap::new(config.n_calibrants),
                    |a, b| a.merge(b),
                )
        };

        #[cfg(feature = "serial_scoring")]
        let heap: CalibrantHeap = {
            let mut scorer = ApexFinder::new(self.num_cycles());
            let mut heap = CalibrantHeap::new(config.n_calibrants);
            for (chunk_idx, item) in items_to_score.iter().enumerate().filter(|(_, x)| filter_fn(x))
            {
                if let Some((loc, _meta)) = self.prescore(item, &mut scorer) {
                    heap.push(CalibrantCandidate {
                        score: loc.score,
                        apex_rt_seconds: loc.retention_time_ms as f32 / 1000.0,
                        speclib_index: speclib_offset + chunk_idx,
                    });
                }
            }
            heap
        };

        heap
    }

    fn map_rt_index_to_milis(&self, rt_index: usize) -> u32 {
        self.index
            .ms1_cycle_mapping()
            .rt_milis_for_index(&MS1CycleIndex::new(rt_index as u32))
            .unwrap_or(0)
    }

    fn num_cycles(&self) -> usize {
        self.index.ms1_cycle_mapping().len()
    }
}
