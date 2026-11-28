//! Peptide scoring pipeline.
//!
//! # Performance-Critical Design: Buffer Reuse
//!
//! The scoring pipeline processes thousands of peptide queries per second. To achieve this
//! throughput, it's critical to minimize allocations in the hot path.
//!
//! Each scoring operation requires several buffers:
//! - `ScoringWorker` holds time-series feature buffers (size varies with query)
//! - Chromatogram collectors (size varies with query complexity)
//!
//! The batch scorers reuse one `ScoringWorker` per Rayon thread across
//! blocks and batches. Without Rayon, one worker scores the entire batch.
//!
//! ## Scoring Pipeline
//!
//! The scoring process has two phases:
//!
//! 1. **Prescore** (Phase 1): Broad extraction + `find_apex_location` -- yields calibrant candidates.
//! 2. **Calibrated scoring** (Phase 3): Narrow calibrated extraction + `find_apex` + secondary query.

use crate::data_sources::reference_library::{
    ExpectedIntensity,
    RefQuery,
    ReferenceLibrary,
    RowHandles,
    ScoredIdentity,
};
use crate::errors::DataProcessingError;
use crate::rt_calibration::{
    CalibrationResult,
    LibraryRT,
    ObservedRTSeconds,
};
use crate::utils::maybe_par::BlockWorkers;
use crate::{
    IonAnnot,
    ScorerQueriable,
    timed,
};
use std::collections::HashMap;
use std::hash::{
    Hash,
    Hasher,
};
use std::sync::Mutex;
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timsquery::models::FlatIdx;
use timsquery::traits::Target;
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    ExtractionQuery,
    KeyLike,
    MzMobilityStatsCollector,
    PeakTolerance,
    ResolvedRt,
    RtSelection,
    SpectralCollector,
    Tolerance,
};
use tracing::warn;

use super::accumulator::IonSearchAccumulator;
use super::apex_finding::{
    ApexBlocks,
    ApexLocation,
    Extraction,
    RelativeIntensityCollector,
    TraceScorer,
};
use super::blocks::poisson::{
    PoissonEvidence,
    compute_poisson_evidence,
};
use super::hyperscore::single_lazyscore;
use super::offsets::MzMobilityOffsets;
use super::results::{
    FinalizeInputs,
    ScoredCandidate,
    ScoringFields,
};
pub use super::shadow::ShadowRun;
use super::skip::{
    SkipCounts,
    SkipReason,
    apex_error_to_skip,
    finalize_error_to_skip,
};
use super::timings::{
    PrescoreTimings,
    ScoreTimings,
};

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Phase3Stage {
    #[default]
    PrecursorRange,
    ExpectedFragments,
    Extraction,
    Apex,
    Finalize,
    Scored,
}

/// One selected scoring slot's Phase 3 outcome. The caller preallocates slots
/// so scoring records only requested candidates, including prefilter rejects.
#[derive(Debug, Default)]
pub struct Phase3Trace {
    pub stage: Phase3Stage,
    pub skip: Option<SkipReason>,
    pub rt_window: Option<ResolvedRt>,
    pub main_score: Option<f32>,
    pub apex_rt_seconds: Option<f32>,
    pub matched_fragments: Option<u8>,
}

struct CandidateIdentity {
    flat: FlatIdx,
    handles: RowHandles,
    source_id: timsquery::models::OwnedSourceId,
    is_target: bool,
}

#[derive(Clone, Copy)]
struct BatchScoring<'a, 'index> {
    calibration: Option<&'a CalibrationResult>,
    shadow: Option<&'a ShadowRun<'index>>,
}

/// Scratch and run inputs shared by successive Phase 1 batches.
pub struct PrescoreBatch<'a, I: ScorerQueriable> {
    scorer: &'a Scorer<I>,
    library: &'a ReferenceLibrary,
    capacity_per_band: usize,
    library_rt_range: Option<(f32, f32)>,
    workers: BlockWorkers<ScoringWorker>,
}

type TraceSlots = HashMap<FlatIdx, Mutex<Phase3Trace>>;

/// Scratch and run inputs shared by successive Phase 3 batches.
pub struct BatchScorer<'a, 'index, I: ScorerQueriable> {
    scorer: &'a Scorer<I>,
    library: &'a ReferenceLibrary,
    batch: BatchScoring<'a, 'index>,
    traces: Option<&'a TraceSlots>,
    workers: BlockWorkers<ScoringWorker>,
}

impl<I: ScorerQueriable> BatchScorer<'_, '_, I> {
    pub fn score(&mut self, flats: &[FlatIdx]) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        if self.traces.is_some() {
            self.scorer.score_batch::<true>(
                self.library,
                flats,
                self.batch,
                self.traces,
                &mut self.workers,
            )
        } else {
            self.scorer.score_batch::<false>(
                self.library,
                flats,
                self.batch,
                None,
                &mut self.workers,
            )
        }
    }
}

/// Per-rayon-worker scoring state. Holds a `TraceScorer` plus reusable scratch:
/// - `Extraction` slot (ChromatogramCollector reset-and-reused, not reallocated)
/// - `inner_collector` / `isotope_collector` reused across Phase 3 secondary queries
pub struct ScoringWorker {
    pub scorer: TraceScorer,
    pub extraction: Option<Extraction<IonAnnot>>,
    pub inner_collector: Option<SpectralCollector<IonAnnot, MzMobilityStatsCollector>>,
    pub isotope_collector: Option<SpectralCollector<IonAnnot, f32>>,
}

impl ScoringWorker {
    pub fn new(num_cycles: usize, max_frags: usize) -> Self {
        Self::with_denoising(num_cycles, max_frags, false)
    }

    pub(super) fn with_denoising(
        num_cycles: usize,
        max_frags: usize,
        denoise_traces: bool,
    ) -> Self {
        Self {
            scorer: TraceScorer::with_denoising(num_cycles, max_frags, denoise_traces),
            extraction: None,
            inner_collector: None,
            isotope_collector: None,
        }
    }
}

/// Lightweight calibrant candidate -- just enough to re-query in Phase 2.
/// Implements Ord by score (ascending) for use in BinaryHeap<Reverse<_>>.
#[derive(Debug, Clone)]
pub struct CalibrantCandidate {
    pub score: f32,
    pub apex_rt: ObservedRTSeconds<f32>,
    pub speclib_index: FlatIdx,
    pub library_rt: LibraryRT<f32>,
}

impl PartialEq for CalibrantCandidate {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for CalibrantCandidate {}

impl PartialOrd for CalibrantCandidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CalibrantCandidate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score
            .total_cmp(&other.score)
            .then_with(|| self.speclib_index.cmp(&other.speclib_index))
    }
}

pub const CALIBRANT_HEAP_GROUPS: usize = 2;
pub const CALIBRANT_RT_BANDS: usize = 3;

/// Fixed, implementation-independent hash for assigning opaque library handles
/// to calibration folds. `FlatIdx` remains responsible for exposing only its
/// `Hash` behavior; this code neither constructs nor decodes one.
struct CalibrantSplitHasher(u64);

impl Default for CalibrantSplitHasher {
    fn default() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }
}

impl Hasher for CalibrantSplitHasher {
    fn write(&mut self, bytes: &[u8]) {
        // FNV-1a followed by an avalanche in `finish`. The split needs stable,
        // well-mixed bits, not adversarial collision resistance.
        let mut hash = self.0;
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        self.0 = hash;
    }

    fn write_u32(&mut self, value: u32) {
        self.write(&value.to_le_bytes());
    }

    fn finish(&self) -> u64 {
        let mut value = self.0;
        value ^= value >> 33;
        value = value.wrapping_mul(0xff51_afd7_ed55_8ccd);
        value ^= value >> 33;
        value = value.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        value ^ (value >> 33)
    }
}

fn calibrant_heap_group(index: FlatIdx) -> usize {
    let mut hasher = CalibrantSplitHasher::default();
    index.hash(&mut hasher);
    (hasher.finish() as usize) % CALIBRANT_HEAP_GROUPS
}

/// Two independent calibration folds, each split into equal-width library-RT
/// bands (observed RT when library RT is absent). This is a bounded, mergeable
/// accumulator for the parallel prescore.
pub struct StratifiedCalibrantHeaps {
    heaps: [[CalibrantHeap; CALIBRANT_RT_BANDS]; CALIBRANT_HEAP_GROUPS],
    rt_min: f32,
    band_scale: f32,
    observed_rt_bands: bool,
}

impl StratifiedCalibrantHeaps {
    pub fn new(capacity_per_band: usize, rt_range: (f32, f32)) -> Self {
        let (rt_min, rt_max) = rt_range;
        let band_scale = if rt_max > rt_min {
            CALIBRANT_RT_BANDS as f32 / (rt_max - rt_min)
        } else {
            0.0
        };
        Self {
            heaps: std::array::from_fn(|_| {
                std::array::from_fn(|_| CalibrantHeap::new(capacity_per_band))
            }),
            rt_min,
            band_scale,
            observed_rt_bands: false,
        }
    }

    pub fn for_library(
        capacity: usize,
        library_range: Option<(f32, f32)>,
        observed_range: (f64, f64),
    ) -> Self {
        let mut heaps = Self::new(
            capacity,
            library_range.unwrap_or((observed_range.0 as f32, observed_range.1 as f32)),
        );
        heaps.observed_rt_bands = library_range.is_none();
        heaps
    }

    fn band_of(&self, library_rt: f32) -> usize {
        if self.band_scale == 0.0 {
            return 0;
        }
        (((library_rt - self.rt_min) * self.band_scale) as usize).min(CALIBRANT_RT_BANDS - 1)
    }

    pub fn push(&mut self, candidate: CalibrantCandidate) -> Result<(), SkipReason> {
        let group = calibrant_heap_group(candidate.speclib_index);
        let band = self.band_of(if self.observed_rt_bands {
            candidate.apex_rt.0
        } else {
            candidate.library_rt.0
        });
        self.heaps[group][band].push(candidate)
    }

    pub fn merge(mut self, other: Self) -> Self {
        for (left_group, right_group) in self.heaps.iter_mut().zip(other.heaps) {
            for (left, right) in left_group.iter_mut().zip(right_group) {
                left.merge_in_place(right);
            }
        }
        self
    }

    pub fn merge_in_place(&mut self, other: Self) {
        for (left_group, right_group) in self.heaps.iter_mut().zip(other.heaps) {
            for (left, right) in left_group.iter_mut().zip(right_group) {
                left.merge_in_place(right);
            }
        }
    }

    pub fn group_iter(&self, group: usize) -> impl Iterator<Item = &CalibrantCandidate> + Clone {
        self.heaps[group].iter().flat_map(CalibrantHeap::iter)
    }

    pub fn iter(&self) -> impl Iterator<Item = &CalibrantCandidate> + Clone {
        self.heaps
            .iter()
            .flat_map(|group| group.iter().flat_map(CalibrantHeap::iter))
    }

    pub fn into_vec(self) -> Vec<CalibrantCandidate> {
        self.heaps
            .into_iter()
            .flatten()
            .flat_map(CalibrantHeap::into_vec)
            .collect()
    }
}

/// Bounded top-N set: keeps only the best `capacity` candidates. The
/// candidate order is total (score, then library index), so the set is always
/// sorted and its contents and order do not depend on insertion or merge order.
pub struct CalibrantHeap {
    best: std::collections::BTreeSet<CalibrantCandidate>,
    capacity: usize,
}

impl CalibrantHeap {
    pub fn new(capacity: usize) -> Self {
        Self {
            best: std::collections::BTreeSet::new(),
            capacity,
        }
    }

    /// Returns `Err(SkipReason::CalibrantNanScore)` when `candidate.score` is
    /// NaN so the caller can surface the count in reports. Non-positive or
    /// +Inf scores are silently discarded -- they are either sentinel
    /// (`score <= 0.0` means "no usable signal") or already dominated by any
    /// finite positive candidate in the heap.
    pub fn push(&mut self, candidate: CalibrantCandidate) -> Result<(), SkipReason> {
        if candidate.score.is_nan() {
            return Err(SkipReason::CalibrantNanScore);
        }
        if !candidate.score.is_finite() || candidate.score <= 0.0 {
            return Ok(());
        }
        if self.best.len() < self.capacity {
            self.best.insert(candidate);
        } else if let Some(min) = self.best.first()
            && candidate > *min
        {
            self.best.pop_first();
            self.best.insert(candidate);
        }
        Ok(())
    }

    pub fn merge(mut self, other: Self) -> Self {
        self.merge_in_place(other);
        self
    }

    pub fn merge_in_place(&mut self, other: Self) {
        for item in other.best {
            // Ignore: other heap already filtered on its own push.
            let _ = self.push(item);
        }
    }

    /// Best first.
    pub fn into_vec(self) -> Vec<CalibrantCandidate> {
        self.best.into_iter().rev().collect()
    }

    pub fn len(&self) -> usize {
        self.best.len()
    }

    /// Best first, like [`Self::into_vec`].
    pub fn iter(&self) -> impl Iterator<Item = &CalibrantCandidate> + Clone {
        self.best.iter().rev()
    }
}

/// Calibration configuration -- all tunable parameters with defaults.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CalibrationConfig {
    pub n_calibrants: usize,
    pub grid_size: usize,
    /// Accepted for older configs. The m/z mixture uses a fixed 3 sigma.
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
            // Keep in sync with `timsseek_cli/assets/default_config.toml`
            // ([calibration]) -- that file carries the rationale for these
            // sigmas, and `config::tests::default_template_parses`
            // fails if the two drift apart.
            mz_sigma: 3.0,
            mobility_sigma: 3.0,
            rt_sigma_factor: 3.0,
            min_rt_tolerance_minutes: 0.5,
            calibration_query_rt_window_minutes: DEFAULT_CALIBRATION_RT_WINDOW_SECONDS / 60.0,
            dp_lookback: 30,
        }
    }
}

/// Number of top fragments to retain for scoring (by predicted intensity).
pub const TOP_N_FRAGMENTS: usize = 8;

/// Library-only pre-gate: reject peptides whose library entry carries no
/// predicted fragments before doing any extraction work. Shared by the
/// broad (prescore) and calibrated (score_calibrated) paths so the two
/// don't drift.
#[inline]
fn gate_expected_fragments(query: &RefQuery<'_>) -> Result<(), SkipReason> {
    if query.fragment_count() == 0 {
        return Err(SkipReason::NoExpectedFragments);
    }
    Ok(())
}

/// Retain only the top `n` fragments by predicted intensity.
///
/// Removes lower-ranked fragments from the collector and expected intensities,
/// keeping their fragment counts aligned.
pub fn select_top_n_fragments<T: KeyLike + Default>(
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
    to_drop.sort_by_key(|a| std::cmp::Reverse(a.0));

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
pub fn filter_zero_intensity_ions<T: KeyLike + Default>(
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

    // Assert collector rows and expected intensities agree.
    assert_eq!(
        agg.precursors.num_ions(),
        expected.precursor_len(),
        "Precursor expected intensities count mismatch after filtering"
    );
}

/// Raw secondary lazyscores produced during scoring; projected into the
/// [`crate::scoring::blocks::lazy::SecondaryLazyScores`] block via `From`.
#[derive(Debug, Clone, Copy)]
pub struct SecondaryLazyScoresRaw {
    pub lazyscore: f32,
    pub iso_lazyscore: f32,
    pub isotope_lazyscore_log_diff: f32,
}

/// Compute lazyscores from inner and isotope collectors.
fn compute_secondary_lazyscores(
    inner: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
    isotope: Option<&SpectralCollector<IonAnnot, f32>>,
) -> SecondaryLazyScoresRaw {
    let lazyscore = single_lazyscore(
        inner
            .iter_fragments()
            .map(|((_k, _mz), v)| v.weight() as f32),
    );
    let iso_lazyscore = isotope.map_or(f32::NAN, |isotope| {
        single_lazyscore(isotope.iter_fragments().map(|((_k, _mz), v)| *v))
    });
    let isotope_lazyscore_log_diff = lazyscore.ln_1p() - iso_lazyscore.ln_1p();
    SecondaryLazyScoresRaw {
        lazyscore,
        iso_lazyscore,
        isotope_lazyscore_log_diff,
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

    /// Smooth chromatograms for Phase 3 shape and coelution features.
    pub denoise_traces: bool,
}

impl<I: ScorerQueriable> Scorer<I> {
    pub fn acquisition_rt_range_seconds(&self) -> (f64, f64) {
        let (lo, hi) = self.index.ms1_cycle_mapping().range_milis();
        (lo as f64 / 1000.0, hi as f64 / 1000.0)
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
    /// Results populate `worker.inner_collector` and `worker.isotope_collector` in place.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn execute_secondary_query(
        &self,
        query: &RefQuery<'_>,
        apex: &ApexBlocks,
        spectral_tol: &Tolerance,
        isotope_tol: &Tolerance,
        worker: &mut ScoringWorker,
    ) {
        let apex_rt = RtSelection::Centered(calibrt::ObservedRTSeconds(
            apex.retention_time_ms as f32 / 1000.0,
        ));
        let resolve = |tolerance: &Tolerance| {
            ResolvedRt::from_mapping(apex_rt, &tolerance.rt, self.index.ms1_cycle_mapping())
                .expect("detected apex lies inside acquisition")
        };
        let spectral_query = ExtractionQuery::new(query, resolve(spectral_tol));
        let inner = worker
            .inner_collector
            .get_or_insert_with(|| SpectralCollector::new(&spectral_query));
        inner.reset_with(&spectral_query);
        self.index.add_query(inner, &spectral_tol.peak_tolerance());
        let mobility = Self::get_mobility(inner);
        let isotope_query =
            ExtractionQuery::new(query, resolve(isotope_tol)).with_mobility(mobility as f32);
        inner.reset_with(&isotope_query);
        self.index.add_query(inner, &isotope_tol.peak_tolerance());
        if !query.scoring_plan().fragment_isotopes().enabled {
            worker.isotope_collector = None;
            return;
        }
        let isotope = worker
            .isotope_collector
            .get_or_insert_with(|| SpectralCollector::new(&isotope_query));
        isotope.reset_with(&isotope_query);
        crate::utils::elution_group_ops::shift_fragment_isotopes(isotope, 1);
        self.index.add_query(isotope, &isotope_tol.peak_tolerance());
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn finalize_results(
        &self,
        metadata: &super::apex_finding::CandidateMetadata,
        nqueries: u8,
        apex: ApexBlocks,
        poisson: PoissonEvidence,
        inner_collector: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>,
        isotope_collector: Option<&SpectralCollector<IonAnnot, f32>>,
    ) -> Result<ScoredCandidate, DataProcessingError> {
        let offsets = MzMobilityOffsets::new(inner_collector, metadata.ref_mobility_ook0 as f64);
        let rel_inten = RelativeIntensityCollector::new(inner_collector);
        let secondary_lazy = compute_secondary_lazyscores(inner_collector, isotope_collector);

        let mut scoring = ScoringFields::compute(FinalizeInputs {
            metadata,
            offsets: &offsets,
            rel_inten,
            secondary_lazy,
            nqueries,
            apex,
            poisson,
        });

        // No searchable mobility axis (mzML/FAIMS) → the observed mobility is a
        // sentinel, so drop every mobility feature to NaN (forust-missing).
        if !self.index.mobility_kind().is_scoreable() {
            scoring.neutralize_mobility();
        }
        Ok(ScoredCandidate::new(scoring))
    }
}

impl<I: ScorerQueriable> Scorer<I> {
    /// Build a chromatogram extraction using calibrated RT and per-query tolerance.
    /// The speclib is NOT mutated -- CalibrationResult provides the RT conversion.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    /// Populate `worker.extraction` with a calibrated-tolerance extraction.
    /// Reuses the worker's backing `ChromatogramCollector` storage.
    fn build_calibrated_extraction_into(
        &self,
        query: &RefQuery<'_>,
        identity: CandidateIdentity,
        calibration: &CalibrationResult,
        worker: &mut ScoringWorker,
    ) -> Result<super::apex_finding::CandidateMetadata, SkipReason> {
        let original_rt = query.library_rt();
        let calibrated_rt = if !calibration.has_rt_calibration() {
            None
        } else {
            original_rt
                .map(|rt| {
                    calibration
                        .convert_irt(rt)
                        .map_err(|_| SkipReason::RetentionTimeOutOfBounds)
                })
                .transpose()?
        };
        if let (Some(original_rt), Some(observed_rt)) = (original_rt, calibrated_rt) {
            let tolerance = calibration.get_tolerance(
                query.mono_precursor_mz(),
                query.mobility_ook0(),
                original_rt,
            );
            let observed_query = ExtractionQuery::new(
                query,
                ResolvedRt::from_mapping(
                    RtSelection::Centered(observed_rt),
                    &tolerance.rt,
                    self.index.ms1_cycle_mapping(),
                )
                .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?,
            );
            super::extraction::build_extraction_into(
                &mut worker.extraction,
                &observed_query,
                &self.index,
                &tolerance.peak_tolerance(),
                Some(TOP_N_FRAGMENTS),
            )?;
        } else {
            let tolerance = calibration.unrestricted_rt_tolerance();
            super::extraction::build_extraction_into(
                &mut worker.extraction,
                &ExtractionQuery::new(
                    query,
                    ResolvedRt::from_mapping(
                        RtSelection::FullRun,
                        &self.broad_tolerance.rt,
                        self.index.ms1_cycle_mapping(),
                    )
                    .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?,
                ),
                &self.index,
                &tolerance.peak_tolerance(),
                Some(TOP_N_FRAGMENTS),
            )?;
        }

        Ok(super::apex_finding::CandidateMetadata {
            is_target: identity.is_target,
            charge: query.precursor_charge(),
            handles: identity.handles,
            source_id: identity.source_id,
            library_rt: original_rt.map_or(f32::NAN, |r| r.0),
            calibrated_rt_seconds: calibrated_rt.map_or(f32::NAN, |r| r.0),
            ref_mobility_ook0: query.mobility_ook0(),
            ref_precursor_mz: query.mono_precursor_mz(),
        })
    }

    /// Score with calibrated extraction, or full-RT extraction when calibration
    /// is absent. Both modes share apex scoring and the enabled secondary queries.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn score_extraction<const TRACE: bool>(
        &self,
        query: &RefQuery<'_>,
        identity: CandidateIdentity,
        batch: BatchScoring<'_, '_>,
        worker: &mut ScoringWorker,
        timings: &mut ScoreTimings,
        mut trace: Option<&mut Phase3Trace>,
    ) -> Result<ScoredCandidate, SkipReason> {
        let calibration = batch.calibration;
        let flat = identity.flat;
        if TRACE && let Some(trace) = trace.as_deref_mut() {
            trace.stage = Phase3Stage::ExpectedFragments;
        }
        gate_expected_fragments(query)?;

        if TRACE && let Some(trace) = trace.as_deref_mut() {
            trace.stage = Phase3Stage::Extraction;
        }

        let metadata = timed!(
            timings.extraction,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::extraction").in_scope(|| {
                match calibration {
                    Some(calibration) => {
                        self.build_calibrated_extraction_into(query, identity, calibration, worker)
                    }
                    None => {
                        let tolerance = self.broad_tolerance.clone().with_rt_tolerance(
                            timsquery::models::tolerance::RtTolerance::Unrestricted,
                        );
                        super::extraction::build_extraction_into(
                            &mut worker.extraction,
                            &ExtractionQuery::new(
                                query,
                                ResolvedRt::from_mapping(
                                    RtSelection::FullRun,
                                    &self.broad_tolerance.rt,
                                    self.index.ms1_cycle_mapping(),
                                )
                                .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?,
                            ),
                            &self.index,
                            &tolerance.peak_tolerance(),
                            Some(TOP_N_FRAGMENTS),
                        )?;
                        Ok(super::apex_finding::CandidateMetadata {
                            is_target: identity.is_target,
                            charge: query.precursor_charge(),
                            handles: identity.handles,
                            source_id: identity.source_id,
                            library_rt: query.library_rt().map_or(f32::NAN, |r| r.0),
                            calibrated_rt_seconds: f32::NAN,
                            ref_mobility_ook0: query.mobility_ook0(),
                            ref_precursor_mz: query.mono_precursor_mz(),
                        })
                    }
                }
            })
        );
        if TRACE
            && let Some(trace) = trace.as_deref_mut()
            && matches!(
                metadata,
                Ok(_) | Err(SkipReason::FragmentsOutsideScanRange | SkipReason::NoObservedSignal)
            )
        {
            // These outcomes all reached a freshly reset collector.
            trace.rt_window = Some(
                worker
                    .extraction
                    .as_ref()
                    .expect("extraction set by query")
                    .chromatograms
                    .rt,
            );
        }
        let metadata = metadata?;

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");

        if TRACE && let Some(trace) = trace.as_deref_mut() {
            trace.stage = Phase3Stage::Apex;
        }

        let (apex_score, poisson) = timed!(
            timings.scoring,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::apex_scoring").in_scope(
                || {
                    let apex = worker
                        .scorer
                        .find_apex(scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                        .map_err(|e| apex_error_to_skip(&e))?;
                    let poisson = calibration
                        .and_then(CalibrationResult::poisson_gain)
                        .map_or(PoissonEvidence::MISSING, |gain| {
                            compute_poisson_evidence(
                                &worker.scorer.traces.apex_profile,
                                apex.joint_apex_cycle,
                                (apex.counts.rising_cycles, apex.counts.falling_cycles),
                                &scoring_ctx.chromatograms.fragments,
                                query.iter_expected_fragments(),
                                gain,
                            )
                        });
                    Ok::<_, SkipReason>((apex, poisson))
                },
            )
        )?;

        timed!(timings.spectral_query, {
            let spectral_tol = calibration.map_or_else(
                || {
                    self.broad_tolerance.clone().with_rt_tolerance(
                        timsquery::models::tolerance::RtTolerance::Minutes((
                            0.5 / 60.0,
                            0.5 / 60.0,
                        )),
                    )
                },
                CalibrationResult::get_spectral_tolerance,
            );
            let isotope_tol = calibration.map_or_else(
                || spectral_tol.clone(),
                CalibrationResult::get_isotope_tolerance,
            );
            tracing::span!(tracing::Level::TRACE, "score_calibrated::secondary_query").in_scope(
                || {
                    self.execute_secondary_query(
                        query,
                        &apex_score,
                        &spectral_tol,
                        &isotope_tol,
                        worker,
                    )
                },
            )
        });
        let inner_collector = worker.inner_collector.as_ref().expect("set by secondary");
        let isotope_collector = worker.isotope_collector.as_ref();

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        let nqueries = scoring_ctx.chromatograms.fragments.num_ions() as u8;
        let shadow_context = batch.shadow.and_then(|shadow| {
            shadow.context(
                query,
                flat,
                scoring_ctx.chromatograms.cycle_offset(),
                apex_score.joint_apex_cycle,
                apex_score.peak_width_cycles,
            )
        });
        if TRACE && let Some(trace) = trace {
            trace.stage = Phase3Stage::Finalize;
            trace.main_score = Some(apex_score.primary.main_score);
            trace.apex_rt_seconds = Some(apex_score.retention_time_ms as f32 / 1000.0);
            trace.matched_fragments = Some(nqueries);
        }
        timed!(
            timings.assembly,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::finalize").in_scope(|| {
                self.finalize_results(
                    &metadata,
                    nqueries,
                    apex_score,
                    poisson,
                    inner_collector,
                    isotope_collector,
                )
                .map(|mut candidate| {
                    candidate.shadow_context = shadow_context;
                    candidate
                })
                .map_err(|e| {
                    warn!("Error in scoring: {:?}", e);
                    finalize_error_to_skip(&e)
                })
            })
        )
    }

    /// Phase 3 batch: Score all peptides with calibrated tolerances.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn score_calibrated_batch(
        &self,
        lib: &ReferenceLibrary,
        flats: &[FlatIdx],
        calibration: &CalibrationResult,
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        self.score_shadow_batch(lib, flats, Some(calibration), None, None)
    }

    /// Score spectra over the acquisition's full RT range, without calibration
    /// or assumptions about the library's RT units. No FDR is assigned here.
    pub fn score_raw_batch(
        &self,
        lib: &ReferenceLibrary,
        flats: &[FlatIdx],
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        self.score_shadow_batch(lib, flats, None, None, None)
    }

    pub fn score_traced_batch(
        &self,
        lib: &ReferenceLibrary,
        flats: &[FlatIdx],
        calibration: Option<&CalibrationResult>,
        traces: &HashMap<FlatIdx, Mutex<Phase3Trace>>,
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        self.score_shadow_batch(lib, flats, calibration, Some(traces), None)
    }

    /// Score a Phase 3 batch and capture the run-local shadow context.
    pub fn score_shadow_batch(
        &self,
        lib: &ReferenceLibrary,
        flats: &[FlatIdx],
        calibration: Option<&CalibrationResult>,
        traces: Option<&HashMap<FlatIdx, Mutex<Phase3Trace>>>,
        shadow: Option<&ShadowRun<'_>>,
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        self.batch_scorer(lib, calibration, traces, shadow)
            .score(flats)
    }

    /// Create a scorer that keeps its scratch across Phase 3 batches.
    pub fn batch_scorer<'a, 'index>(
        &'a self,
        library: &'a ReferenceLibrary,
        calibration: Option<&'a CalibrationResult>,
        traces: Option<&'a TraceSlots>,
        shadow: Option<&'a ShadowRun<'index>>,
    ) -> BatchScorer<'a, 'index, I> {
        BatchScorer {
            scorer: self,
            library,
            batch: BatchScoring {
                calibration,
                shadow,
            },
            traces,
            workers: BlockWorkers::new(),
        }
    }

    fn score_batch<const TRACE: bool>(
        &self,
        lib: &ReferenceLibrary,
        flats: &[FlatIdx],
        batch: BatchScoring<'_, '_>,
        traces: Option<&TraceSlots>,
        workers: &mut BlockWorkers<ScoringWorker>,
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        let get_item = |f| lib.item_at(f);
        let num_cycles = self.num_cycles();
        // `fold_reduce` scores contiguous blocks, driving the flyweight by
        // index; the flyweight itself is never stored.
        // Precursor-range gate over the flyweight geometry.
        let filter_fn = |q: &RefQuery<'_>| {
            let tmp = q.precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };
        // All Phase 3 extractions retain at most this many fragment rows.
        let max_frags = TOP_N_FRAGMENTS;

        // One-shot count of pre-filter rejects. filter_fn is cheap (range check);
        // a second pass keeps the rayon hot loop free of synchronised counters.
        let precursor_oofr_count =
            flats.iter().filter(|&&f| !filter_fn(&get_item(f))).count() as u32;

        let mut results = workers.fold_reduce(
            flats,
            || {
                tracing::debug!(
                    target: "alloc_track",
                    phase = "phase3_score",
                    num_cycles,
                    max_frags,
                    "worker init: ScoringWorker"
                );
                ScoringWorker::with_denoising(num_cycles, max_frags, self.denoise_traces)
            },
            |worker, chunk| {
                let mut acc = IonSearchAccumulator::default();
                for &flat in chunk {
                    let q = get_item(flat);
                    let mut trace = if TRACE {
                        traces
                            .expect("traced batch has trace slots")
                            .get(&flat)
                            .map(|entry| entry.lock().expect("trace slot poisoned"))
                    } else {
                        None
                    };
                    if !filter_fn(&q) {
                        if TRACE && let Some(trace) = trace.as_deref_mut() {
                            trace.skip = Some(SkipReason::PrecursorOutOfFragmentedRange);
                        }
                        continue;
                    }
                    let _span =
                        tracing::span!(tracing::Level::TRACE, "score_calibrated_item").entered();
                    let mut t = ScoreTimings::default();
                    let identity = CandidateIdentity {
                        flat,
                        handles: q.handles(),
                        source_id: q.output_id().to_owned_id(),
                        is_target: q.is_target(),
                    };
                    let result = self.score_extraction::<TRACE>(
                        &q,
                        identity,
                        batch,
                        worker,
                        &mut t,
                        trace.as_deref_mut(),
                    );
                    if TRACE && let Some(mut trace) = trace {
                        match &result {
                            Ok(_) => trace.stage = Phase3Stage::Scored,
                            Err(reason) => trace.skip = Some(*reason),
                        }
                    }
                    acc.fold((result, t));
                }
                acc
            },
            IonSearchAccumulator::merge,
        );

        results.skips.precursor_out_of_fragmented_range = precursor_oofr_count;
        (results.res, results.timings, results.skips)
    }

    /// Phase 1: Lightweight prescore -- broad extraction + find_apex_location only.
    /// Returns the [`ApexLocation`] the apex finder picked.
    pub fn prescore(
        &self,
        query: &RefQuery<'_>,
        worker: &mut ScoringWorker,
        timings: &mut PrescoreTimings,
    ) -> Result<ApexLocation, SkipReason> {
        self.prescore_with_windows(
            query,
            worker,
            timings,
            self.prescore_rt(),
            &self.broad_tolerance.peak_tolerance(),
        )
    }

    fn prescore_rt(&self) -> Option<ResolvedRt> {
        ResolvedRt::from_mapping(
            RtSelection::FullRun,
            &self.broad_tolerance.rt,
            self.index.ms1_cycle_mapping(),
        )
        .ok()
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(name = "prescore", skip_all, level = "trace")
    )]
    fn prescore_with_windows(
        &self,
        query: &RefQuery<'_>,
        worker: &mut ScoringWorker,
        timings: &mut PrescoreTimings,
        rt: Option<ResolvedRt>,
        peak_tolerance: &PeakTolerance,
    ) -> Result<ApexLocation, SkipReason> {
        gate_expected_fragments(query)?;

        timed!(
            timings.extraction,
            tracing::span!(tracing::Level::TRACE, "prescore::extraction").in_scope(|| {
                super::extraction::build_extraction_into(
                    &mut worker.extraction,
                    &ExtractionQuery::new(query, rt.ok_or(SkipReason::RetentionTimeOutOfBounds)?),
                    &self.index,
                    peak_tolerance,
                    Some(TOP_N_FRAGMENTS),
                )
            })
        )?;

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");

        let result = timed!(
            timings.scoring,
            tracing::span!(tracing::Level::TRACE, "prescore::scoring").in_scope(|| {
                worker
                    .scorer
                    .find_apex_location(scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                    .map_err(|e| apex_error_to_skip(&e))
            })
        );

        if result.is_ok() {
            timings.n_scored += 1;
        }
        result
    }

    /// Phase 1 batch: prescore all peptides into bounded A/B, library-RT-stratified heaps.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn prescore_partitioned_batch(
        &self,
        lib: &ReferenceLibrary,
        flats: &[FlatIdx],
        capacity_per_band: usize,
        library_rt_range: Option<(f32, f32)>,
        timings: &mut PrescoreTimings,
    ) -> StratifiedCalibrantHeaps {
        self.prescore_batches(lib, capacity_per_band, library_rt_range)
            .score(flats, timings)
    }

    /// Create a scorer that keeps its scratch across Phase 1 batches.
    pub fn prescore_batches<'a>(
        &'a self,
        library: &'a ReferenceLibrary,
        capacity_per_band: usize,
        library_rt_range: Option<(f32, f32)>,
    ) -> PrescoreBatch<'a, I> {
        PrescoreBatch {
            scorer: self,
            library,
            capacity_per_band,
            library_rt_range,
            workers: BlockWorkers::new(),
        }
    }

    fn map_rt_index_to_milis(&self, rt_index: usize) -> u32 {
        self.index
            .ms1_cycle_mapping()
            .rt_milis_for_index(&MS1CycleIndex::new(rt_index as u32))
            .expect("extraction cycles are clamped to the run's cycle mapping")
    }

    fn num_cycles(&self) -> usize {
        self.index.ms1_cycle_mapping().len()
    }
}

impl<I: ScorerQueriable> PrescoreBatch<'_, I> {
    /// Score a batch into bounded A/B, library-RT-stratified heaps.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn score(
        &mut self,
        flats: &[FlatIdx],
        timings: &mut PrescoreTimings,
    ) -> StratifiedCalibrantHeaps {
        let scorer = self.scorer;
        let library = self.library;
        let get_item = |f| library.item_at(f);
        // The flat index IS the global speclib index (see `Speclib::item_at`),
        // so it doubles as `CalibrantCandidate::speclib_index` -- no separate
        // chunk offset needed.
        let filter_fn = |q: &RefQuery<'_>| {
            let tmp = q.precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            scorer.fragmented_range.intersects(lims)
        };

        let num_cycles = scorer.num_cycles();
        let rt = scorer.prescore_rt();
        let peak_tolerance = scorer.broad_tolerance.peak_tolerance();
        let capacity_per_band = self.capacity_per_band;
        let library_rt_range = self.library_rt_range;
        let (heap, par_timings) = self.workers.fold_reduce(
            flats,
            || {
                tracing::debug!(
                    target: "alloc_track",
                    phase = "phase1_prescore",
                    num_cycles,
                    "worker init: ScoringWorker"
                );
                // Phase 1 extractions keep at most `TOP_N_FRAGMENTS` rows, so
                // one worker fits every batch it is reused for.
                ScoringWorker::new(num_cycles, TOP_N_FRAGMENTS)
            },
            |worker, chunk| {
                let mut heap = StratifiedCalibrantHeaps::for_library(
                    capacity_per_band,
                    library_rt_range,
                    scorer.acquisition_rt_range_seconds(),
                );
                let mut t = PrescoreTimings::default();
                for &flat in chunk {
                    let q = get_item(flat);
                    if !filter_fn(&q) {
                        t.skips.bump(SkipReason::PrecursorOutOfFragmentedRange);
                        continue;
                    }
                    t.n_passed_filter += 1;
                    match scorer.prescore_with_windows(&q, worker, &mut t, rt, &peak_tolerance) {
                        Ok(loc) => {
                            let cand = CalibrantCandidate {
                                score: loc.profile_value,
                                apex_rt: ObservedRTSeconds(loc.retention_time_ms as f32 / 1000.0),
                                speclib_index: flat,
                                library_rt: q.library_rt().unwrap_or(LibraryRT(f32::NAN)),
                            };
                            if let Err(reason) = heap.push(cand) {
                                t.skips.bump(reason);
                            }
                        }
                        Err(reason) => t.skips.bump(reason),
                    }
                }
                (heap, t)
            },
            |(heap, timings), (partial_heap, partial_timings)| {
                heap.merge_in_place(partial_heap);
                *timings += partial_timings;
            },
        );
        *timings += par_timings;
        heap
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_sources::reference_library::ReferenceLibrary;
    use timsquery::models::capabilities::TargetCapabilities;
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };

    fn tiny_lazy_lib() -> ReferenceLibrary {
        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        geom.push_row(Row {
            precursor_mz: 900.4,
            charge: 2,
            rt: Some(timsquery::models::RtCoordinate::seconds(1.0)),
            mobility: 1.0,
            frags: &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y8").unwrap(), 800.0),
            ],
            analyte: timsquery::chemistry::analyte::Analyte::from_sequence("PEPTIDEK").as_input(),
            ..Default::default()
        });
        // `Force` so the arena derives the decoy variants these tests score.
        let geom = geom
            .seal(crate::models::DecoyPolicy::Force)
            .expect("fixture ids are usable");
        ReferenceLibrary::from_sealed_arena(timsquery::serde::TargetTable::Mzpaf {
            geom,
            frag_intens: Some(vec![1.0, 0.5]),
        })
        .expect("fixture is a valid reference library")
    }

    #[test]
    fn calibrant_heap_ties_are_reduction_order_independent() {
        let lib = tiny_lazy_lib();
        let mut flats = lib.geometry().flats();
        let first = flats.next().unwrap();
        let second = flats.next().unwrap();
        let candidate = |speclib_index| CalibrantCandidate {
            score: 10.0,
            apex_rt: ObservedRTSeconds(20.0),
            speclib_index,
            library_rt: LibraryRT(10.0),
        };

        let mut forward = CalibrantHeap::new(1);
        forward.push(candidate(first)).unwrap();
        forward.push(candidate(second)).unwrap();
        let mut reverse = CalibrantHeap::new(1);
        reverse.push(candidate(second)).unwrap();
        reverse.push(candidate(first)).unwrap();

        assert_eq!(
            forward.iter().next().unwrap().speclib_index,
            reverse.iter().next().unwrap().speclib_index
        );
    }

    #[test]
    fn calibrant_heap_order_is_independent_of_insertion_order() {
        let lib = tiny_lazy_lib();
        let flats: Vec<_> = lib.geometry().flats().collect();
        let candidate = |rank: usize| CalibrantCandidate {
            score: 1.0 + rank as f32,
            apex_rt: ObservedRTSeconds(20.0),
            speclib_index: flats[rank % flats.len()],
            library_rt: LibraryRT(10.0),
        };
        let filled = |order: &[usize]| {
            let mut heap = CalibrantHeap::new(4);
            for &rank in order {
                heap.push(candidate(rank)).unwrap();
            }
            heap
        };
        let forward = filled(&[0, 1, 2, 3, 4, 5]);
        let shuffled = filled(&[3, 5, 0, 4, 1, 2]);
        let expected: Vec<_> = [5, 4, 3, 2].map(candidate).to_vec();

        assert!(forward.iter().eq(&expected));
        assert!(shuffled.iter().eq(&expected));
        assert_eq!(forward.into_vec(), expected);
        assert_eq!(shuffled.into_vec(), expected);
    }

    #[test]
    fn rt_free_candidates_are_bounded_by_observed_rt_bands() {
        let lib = tiny_lazy_lib();
        let flat = lib.geometry().flats().next().unwrap();
        let make = || StratifiedCalibrantHeaps::for_library(1, None, (0.0, 300.0));
        let mut direct = make();
        let mut left = make();
        let mut right = make();
        for i in 0..30 {
            let candidate = CalibrantCandidate {
                score: i as f32 + 1.0,
                apex_rt: ObservedRTSeconds(i as f32 * 10.0),
                library_rt: LibraryRT(f32::NAN),
                speclib_index: flat,
            };
            direct.push(candidate.clone()).unwrap();
            if i % 2 == 0 {
                left.push(candidate).unwrap();
            } else {
                right.push(candidate).unwrap();
            }
        }
        left.merge_in_place(right);
        let scores = |heap: StratifiedCalibrantHeaps| {
            let mut scores: Vec<_> = heap.into_vec().into_iter().map(|c| c.score).collect();
            scores.sort_by(f32::total_cmp);
            scores
        };
        let retained = scores(direct);
        assert_eq!(retained.len(), CALIBRANT_RT_BANDS);
        assert_eq!(retained, scores(left));
    }

    #[test]
    fn stratified_calibrant_merge_matches_direct_fold() {
        let lib = tiny_lazy_lib();
        let flats: Vec<_> = lib.geometry().flats().collect();
        let candidates: Vec<_> = (0..12)
            .map(|index| CalibrantCandidate {
                score: index as f32 + 1.0,
                apex_rt: ObservedRTSeconds(index as f32),
                speclib_index: flats[index % flats.len()],
                library_rt: LibraryRT((index % CALIBRANT_RT_BANDS) as f32 * 10.0 + 1.0),
            })
            .collect();
        let mut direct = StratifiedCalibrantHeaps::new(2, (0.0, 30.0));
        for candidate in candidates.iter().cloned() {
            direct.push(candidate).unwrap();
        }
        let mut left = StratifiedCalibrantHeaps::new(2, (0.0, 30.0));
        let mut right = StratifiedCalibrantHeaps::new(2, (0.0, 30.0));
        for candidate in candidates[..6].iter().cloned() {
            left.push(candidate).unwrap();
        }
        for candidate in candidates[6..].iter().cloned() {
            right.push(candidate).unwrap();
        }
        left.merge_in_place(right);

        let scores = |heaps: &StratifiedCalibrantHeaps| {
            let mut scores: Vec<_> = heaps.iter().map(|candidate| candidate.score).collect();
            scores.sort_unstable_by(f32::total_cmp);
            scores
        };
        assert_eq!(scores(&direct), scores(&left));
    }
}
