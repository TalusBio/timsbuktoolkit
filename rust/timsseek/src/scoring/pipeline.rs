//! Peptide scoring pipeline.
//!
//! # Performance-Critical Design: Buffer Reuse
//!
//! The scoring pipeline processes thousands of peptide queries per second. To achieve this
//! throughput, it's critical to minimize allocations in the hot path.
//!
//! Each scoring operation requires several buffers:
//! - `TraceScorer` holds time-series feature buffers (size varies with query)
//! - Chromatogram collectors (size varies with query complexity)
//!
//! `prescore_batch` and `score_calibrated_batch` use Rayon's `map_init()` / `fold` to create
//! one `TraceScorer` buffer per thread, which is reused across thousands of queries.
//!
//! ## Scoring Pipeline
//!
//! The scoring process has two phases:
//!
//! 1. **Prescore** (Phase 1): Broad extraction + `find_apex_location` — yields calibrant candidates.
//! 2. **Calibrated scoring** (Phase 3): Narrow calibrated extraction + `find_apex` + secondary query.

use crate::errors::DataProcessingError;
use crate::{
    IonAnnot,
    QueryItemToScore,
    ScorerQueriable,
    timed,
};
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    KeyLike,
    MzMobilityStatsCollector,
    SpectralCollector,
    Tolerance,
};

use super::accumulator::IonSearchAccumulator;
use super::apex_finding::{
    ApexLocation,
    ApexScore,
    Extraction,
    RelativeIntensities,
    TraceScorer,
};
use super::full_results::ViewerResult;
use super::hyperscore::single_lazyscore;
use super::offsets::MzMobilityOffsets;
use super::results::{
    ScoredCandidate,
    ScoredCandidateBuilder,
};
use super::timings::ScoreTimings;
use crate::rt_calibration::{
    CalibrationResult,
    LibraryRT,
    ObservedRTSeconds,
};
use tracing::warn;

/// Per-rayon-worker scoring state. Holds a `TraceScorer` plus reusable scratch:
/// - `Extraction` slot (ChromatogramCollector reset-and-reused, not reallocated)
/// - `inner_collector` / `isotope_collector` reused across Phase 3 secondary queries
/// - `isotope_scratch_eg` holds the neutron-offset-applied eg; `Option<>` because
///   `TimsElutionGroup` has no `Default` (bon builder with required fields) —
///   init lazily on first peptide.
pub struct ScoringWorker {
    pub scorer: TraceScorer,
    pub extraction: Option<Extraction<IonAnnot>>,
    pub inner_collector: Option<SpectralCollector<IonAnnot, MzMobilityStatsCollector>>,
    pub isotope_collector: Option<SpectralCollector<IonAnnot, f32>>,
    pub isotope_scratch_eg: Option<timsquery::TimsElutionGroup<IonAnnot>>,
}

impl ScoringWorker {
    pub fn new(num_cycles: usize, max_frags: usize) -> Self {
        Self {
            scorer: TraceScorer::new(num_cycles, max_frags),
            extraction: None,
            inner_collector: None,
            isotope_collector: None,
            isotope_scratch_eg: None,
        }
    }
}

/// Lightweight calibrant candidate — just enough to re-query in Phase 2.
/// Implements Ord by score (ascending) for use in BinaryHeap<Reverse<_>>.
#[derive(Debug, Clone)]
pub struct CalibrantCandidate {
    pub score: f32,
    pub apex_rt: ObservedRTSeconds<f32>,
    pub speclib_index: usize,
    pub library_rt: LibraryRT<f32>,
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
        self.score.total_cmp(&other.score)
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
        self.heap
            .into_iter()
            .map(|std::cmp::Reverse(c)| c)
            .collect()
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }

    /// Iterate over heap contents. Order is arbitrary (not sorted by score).
    pub fn iter(&self) -> impl Iterator<Item = &CalibrantCandidate> {
        self.heap.iter().map(|r| &r.0)
    }
}

/// Calibration configuration — all tunable parameters with defaults.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CalibrationConfig {
    pub n_calibrants: usize,
    pub grid_size: usize,
    pub mz_sigma: f32,
    pub mobility_sigma: f32,
    pub rt_sigma_factor: f32,
    pub min_rt_tolerance_minutes: f32,
    pub calibration_query_rt_window_minutes: f32,
    pub dp_lookback: usize,
}

/// Phase-2 per-calibrant RT aggregation window. Narrow enough to isolate the
/// apex scan's neighborhood (roughly a handful of cycles) while still averaging
/// across peak width; wider admits neighbor-RT contamination into the stats.
const DEFAULT_CALIBRATION_RT_WINDOW_SECONDS: f32 = 5.0;

impl Default for CalibrationConfig {
    fn default() -> Self {
        Self {
            n_calibrants: 2000,
            grid_size: 100,
            mz_sigma: 1.5,
            mobility_sigma: 3.0,
            rt_sigma_factor: 3.0,
            min_rt_tolerance_minutes: 0.5,
            calibration_query_rt_window_minutes: DEFAULT_CALIBRATION_RT_WINDOW_SECONDS / 60.0,
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
pub(crate) fn select_top_n_fragments<T: KeyLike + Default>(
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
            let intensity = expected.get_fragment(key).unwrap_or(0.0);
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
        expected.remove_fragment(&key);
    }

    debug_assert_eq!(agg.fragments.num_ions(), expected.fragment_len());
}

/// Filter out zero-intensity ions and update expected intensities in one pass.
///
/// This maintains index alignment by removing ions from the chromatogram collector
/// and expected intensities simultaneously in a single loop.
#[cfg_attr(
    feature = "instrumentation",
    tracing::instrument(skip_all, level = "trace")
)]
pub(crate) fn filter_zero_intensity_ions<T: KeyLike + Default>(
    agg: &mut ChromatogramCollector<T, f32>,
    expected: &mut crate::ExpectedIntensities<T>,
) {
    // Early-exit predicate: stop at first non-zero value (much faster than summing)
    let predicate = |chrom: &[f32]| chrom.iter().any(|&x| x > 0.0);

    // Filter precursors
    for (k, _mz) in agg.drain_nonmatching_precursors(predicate) {
        expected.remove_precursor(k);
    }

    // Filter fragments
    for (k, _mz) in agg.drain_nonmatching_fragments(predicate) {
        expected.remove_fragment(&k);
    }

    // Assert arrays and expected masks agree (eg is no longer mutated here).
    assert_eq!(
        agg.precursors.num_ions(),
        expected.precursor_len(),
        "Precursor expected intensities count mismatch after filtering"
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

pub enum SkippingReason {
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
        let extraction = super::extraction::build_extraction(
            &item.query,
            item.expected_intensity.clone(),
            &self.index,
            &self.broad_tolerance,
            Some(TOP_N_FRAGMENTS),
        )?;

        let library_rt = item.query.rt_seconds();
        let metadata = super::apex_finding::PeptideMetadata {
            digest: item.digest.clone(),
            charge: item.query.precursor_charge(),
            library_id: extraction.chromatograms.id as u32,
            library_rt,
            calibrated_rt_seconds: library_rt, // no calibration in broad path
            ref_mobility_ook0: item.query.mobility_ook0(),
            ref_precursor_mz: item.query.mono_precursor_mz(),
        };

        Ok((metadata, extraction))
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
        tracing::instrument(
            skip(self, item, main_score, spectral_tol, isotope_tol),
            level = "trace"
        )
    )]
    /// Performs refined secondary query at detected apex with two-pass strategy.
    /// Results populate `worker.inner_collector` and `worker.isotope_collector` in place.
    fn execute_secondary_query(
        &self,
        item: &QueryItemToScore,
        main_score: &ApexScore,
        spectral_tol: &Tolerance,
        isotope_tol: &Tolerance,
        worker: &mut ScoringWorker,
    ) {
        let new_rt_seconds = main_score.retention_time_ms as f32 / 1000.0;

        // **Pass 1**: query at apex RT to determine observed mobility.
        let inner = worker
            .inner_collector
            .get_or_insert_with(|| SpectralCollector::new(&item.query));
        inner.reset_with_overrides(&item.query, Some(new_rt_seconds), None);
        self.index.add_query(inner, spectral_tol);

        let mobility = Self::get_mobility(inner);

        // **Pass 2**: same collector, now with mobility override.
        let inner = worker.inner_collector.as_mut().expect("init above");
        inner.reset_with_overrides(&item.query, Some(new_rt_seconds), Some(mobility as f32));

        // Isotope scratch eg holds item.query with +1 neutron offset applied
        // (buffer-override — reuses Vec capacity after warm-up).
        let scratch_eg = worker
            .isotope_scratch_eg
            .get_or_insert_with(|| item.query.clone());
        crate::utils::elution_group_ops::apply_isotope_offset_fragments_into(
            scratch_eg,
            &item.query,
            1i8,
        );

        let isotope = worker
            .isotope_collector
            .get_or_insert_with(|| SpectralCollector::new(scratch_eg));
        isotope.reset_with_overrides(scratch_eg, Some(new_rt_seconds), Some(mobility as f32));

        // Both queries share the same isotope_tol per existing logic.
        let inner = worker.inner_collector.as_mut().expect("init above");
        self.index.add_query(inner, isotope_tol);
        let isotope = worker.isotope_collector.as_mut().expect("init above");
        self.index.add_query(isotope, isotope_tol);
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
    pub fn score_for_viewer(
        &self,
        item: QueryItemToScore,
        calibration: &CalibrationResult,
    ) -> Result<ViewerResult, DataProcessingError> {
        // One-shot worker for the viewer path (not hot).
        let mut worker =
            ScoringWorker::new(self.num_cycles(), item.expected_intensity.fragment_len());

        let (metadata, scoring_ctx) = self.build_broad_extraction(&item).map_err(|_| {
            DataProcessingError::ExpectedNonEmptyData {
                context: Some("RT out of bounds".into()),
            }
        })?;

        let apex_score = worker
            .scorer
            .find_apex(&scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))?;
        let spectral_tol = calibration.get_spectral_tolerance();
        let isotope_tol = calibration.get_isotope_tolerance();
        self.execute_secondary_query(&item, &apex_score, &spectral_tol, &isotope_tol, &mut worker);
        let inner_collector = worker.inner_collector.as_ref().expect("set by secondary");
        let isotope_collector = worker.isotope_collector.as_ref().expect("set by secondary");

        let nqueries = scoring_ctx.chromatograms.fragments.num_ions() as u8;
        let scored = self.finalize_results(
            &metadata,
            nqueries,
            &apex_score,
            inner_collector,
            isotope_collector,
        )?;

        // Extract chromatograms before it's consumed
        let chromatograms = scoring_ctx.chromatograms;

        Ok(ViewerResult {
            traces: worker.scorer.traces.clone(),
            longitudinal_apex_profile: worker.scorer.traces.apex_profile.clone(),
            chromatograms,
            scored,
        })
    }

    /// Build a chromatogram extraction using calibrated RT and per-query tolerance.
    /// The speclib is NOT mutated — CalibrationResult provides the RT conversion.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    /// Populate `worker.extraction` with a calibrated-tolerance extraction.
    /// Reuses the worker's backing `ChromatogramCollector` storage.
    fn build_calibrated_extraction_into(
        &self,
        item: &QueryItemToScore,
        calibration: &CalibrationResult,
        worker: &mut ScoringWorker,
    ) -> Result<super::apex_finding::PeptideMetadata, SkippingReason> {
        let original_irt = LibraryRT(item.query.rt_seconds());
        let calibrated_rt = calibration.convert_irt(original_irt);
        let tolerance = calibration.get_tolerance(
            item.query.mono_precursor_mz(),
            item.query.mobility_ook0(),
            original_irt, // library RT — ridge widths are indexed by library RT
        );

        super::extraction::build_extraction_into(
            &mut worker.extraction,
            &item.query,
            Some(calibrated_rt.0),
            item.expected_intensity.clone(),
            &self.index,
            &tolerance,
            Some(TOP_N_FRAGMENTS),
        )?;

        let extr = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        Ok(super::apex_finding::PeptideMetadata {
            digest: item.digest.clone(),
            charge: item.query.precursor_charge(),
            library_id: extr.chromatograms.id as u32,
            library_rt: original_irt.0 as f32,
            calibrated_rt_seconds: calibrated_rt.0,
            ref_mobility_ook0: item.query.mobility_ook0(),
            ref_precursor_mz: item.query.mono_precursor_mz(),
        })
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
        worker: &mut ScoringWorker,
        timings: &mut ScoreTimings,
    ) -> Option<ScoredCandidate> {
        let metadata = timed!(
            timings.extraction,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::extraction").in_scope(|| {
                self.build_calibrated_extraction_into(item, calibration, worker)
                    .ok()
            })
        )?;

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        if scoring_ctx
            .expected_intensities
            .fragment_intensities
            .is_empty()
        {
            return None;
        }

        let apex_score = timed!(
            timings.scoring,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::apex_scoring").in_scope(
                || {
                    worker
                        .scorer
                        .find_apex(scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                        .ok()
                },
            )
        )?;

        timed!(timings.spectral_query, {
            let spectral_tol = calibration.get_spectral_tolerance();
            let isotope_tol = calibration.get_isotope_tolerance();
            tracing::span!(tracing::Level::TRACE, "score_calibrated::secondary_query").in_scope(
                || {
                    self.execute_secondary_query(
                        item,
                        &apex_score,
                        &spectral_tol,
                        &isotope_tol,
                        worker,
                    )
                },
            )
        });
        let inner_collector = worker.inner_collector.as_ref().expect("set by secondary");
        let isotope_collector = worker.isotope_collector.as_ref().expect("set by secondary");

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        let nqueries = scoring_ctx.chromatograms.fragments.num_ions() as u8;
        let out = timed!(
            timings.assembly,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::finalize").in_scope(|| {
                self.finalize_results(
                    &metadata,
                    nqueries,
                    &apex_score,
                    &inner_collector,
                    &isotope_collector,
                )
            },)
        );

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
        let num_cycles = self.num_cycles();
        // Pre-scan so scratch capacity holds every peptide — no realloc in hot path.
        let max_frags = items_to_score
            .iter()
            .map(|i| i.expected_intensity.fragment_len())
            .max()
            .unwrap_or(0);
        let init_fn = || {
            tracing::debug!(
                target: "alloc_track",
                phase = "phase3_score",
                num_cycles,
                max_frags,
                "rayon worker init: ScoringWorker"
            );
            ScoringWorker::new(num_cycles, max_frags)
        };
        let filter_fn = |x: &&QueryItemToScore| {
            let tmp = x.query.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };

        #[cfg(feature = "rayon")]
        let results: IonSearchAccumulator = {
            items_to_score
                .par_iter()
                .filter(filter_fn)
                .map_init(init_fn, |worker, item| {
                    let mut t = ScoreTimings::default();
                    let result =
                        self.score_calibrated_extraction(item, calibration, worker, &mut t);
                    (result, t)
                })
                .collect()
        };

        #[cfg(not(feature = "rayon"))]
        let results: IonSearchAccumulator = {
            let mut scorer = init_fn();
            items_to_score
                .iter()
                .filter(filter_fn)
                .map(|item| {
                    let _span =
                        tracing::span!(tracing::Level::TRACE, "score_calibrated_item").entered();
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
        worker: &mut ScoringWorker,
        timings: &mut super::timings::PrescoreTimings,
    ) -> Option<ApexLocation> {
        let extraction_ok = timed!(
            timings.extraction,
            tracing::span!(tracing::Level::TRACE, "prescore::extraction").in_scope(|| {
                super::extraction::build_extraction_into(
                    &mut worker.extraction,
                    &item.query,
                    None,
                    item.expected_intensity.clone(),
                    &self.index,
                    &self.broad_tolerance,
                    Some(TOP_N_FRAGMENTS),
                )
                .ok()
            })
        );
        extraction_ok?;

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        if scoring_ctx
            .expected_intensities
            .fragment_intensities
            .is_empty()
        {
            return None;
        }

        let result = timed!(
            timings.scoring,
            tracing::span!(tracing::Level::TRACE, "prescore::scoring").in_scope(|| {
                worker
                    .scorer
                    .find_apex_location(scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                    .ok()
            })
        );

        if result.is_some() {
            timings.n_scored += 1;
        }
        result
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
        timings: &mut super::timings::PrescoreTimings,
    ) -> CalibrantHeap {
        let filter_fn = |x: &&QueryItemToScore| {
            let tmp = x.query.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };

        #[cfg(feature = "rayon")]
        let (heap, par_timings): (CalibrantHeap, super::timings::PrescoreTimings) = {
            use super::timings::PrescoreTimings;
            let num_cycles = self.num_cycles();
            let n_calibrants = config.n_calibrants;
            let max_frags = items_to_score
                .iter()
                .map(|i| i.expected_intensity.fragment_len())
                .max()
                .unwrap_or(0);
            let init_fn = || {
                tracing::debug!(
                    target: "alloc_track",
                    phase = "phase1_prescore",
                    num_cycles,
                    n_calibrants,
                    max_frags,
                    "rayon worker init: ScoringWorker + CalibrantHeap"
                );
                (
                    ScoringWorker::new(num_cycles, max_frags),
                    CalibrantHeap::new(n_calibrants),
                    PrescoreTimings::default(),
                )
            };

            items_to_score
                .par_iter()
                .enumerate()
                .filter(|(_, x)| filter_fn(x))
                .fold(
                    init_fn,
                    |(mut worker, mut heap, mut t), (chunk_idx, item)| {
                        t.n_passed_filter += 1;
                        if let Some(loc) = self.prescore(item, &mut worker, &mut t) {
                            heap.push(CalibrantCandidate {
                                score: loc.score,
                                apex_rt: ObservedRTSeconds(loc.retention_time_ms as f32 / 1000.0),
                                speclib_index: speclib_offset + chunk_idx,
                                library_rt: LibraryRT(item.query.rt_seconds()),
                            });
                        }
                        (worker, heap, t)
                    },
                )
                .map(|(_, heap, t)| (heap, t))
                .reduce(
                    || {
                        (
                            CalibrantHeap::new(config.n_calibrants),
                            PrescoreTimings::default(),
                        )
                    },
                    |(a_heap, mut a_t), (b_heap, b_t)| {
                        a_t += b_t;
                        (a_heap.merge(b_heap), a_t)
                    },
                )
        };
        #[cfg(feature = "rayon")]
        {
            *timings += par_timings;
        }

        #[cfg(not(feature = "rayon"))]
        let heap: CalibrantHeap = {
            let max_frags = items_to_score
                .iter()
                .map(|i| i.expected_intensity.fragment_len())
                .max()
                .unwrap_or(0);
            let mut scorer = TraceScorer::new(self.num_cycles(), max_frags);
            let mut heap = CalibrantHeap::new(config.n_calibrants);
            for (chunk_idx, item) in items_to_score.iter().enumerate() {
                if !filter_fn(&item) {
                    continue;
                }
                timings.n_passed_filter += 1;
                if let Some(loc) = self.prescore(item, &mut scorer, timings) {
                    heap.push(CalibrantCandidate {
                        score: loc.score,
                        apex_rt: ObservedRTSeconds(loc.retention_time_ms as f32 / 1000.0),
                        speclib_index: speclib_offset + chunk_idx,
                        library_rt: LibraryRT(item.query.rt_seconds()),
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
