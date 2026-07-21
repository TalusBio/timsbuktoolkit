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
//! `prescore_batch` and `score_calibrated_batch` dispatch through
//! `super::maybe_par::fold_reduce`, which is compile-time gated on the
//! `rayon` feature. With rayon, each worker thread gets its own
//! `ScoringWorker` via the `init` closure and reuses it across thousands
//! of queries; without rayon, a single `ScoringWorker` is built and reused
//! for the whole batch. Either way the same closure bodies run.
//!
//! ## Scoring Pipeline
//!
//! The scoring process has two phases:
//!
//! 1. **Prescore** (Phase 1): Broad extraction + `find_apex_location` — yields calibrant candidates.
//! 2. **Calibrated scoring** (Phase 3): Narrow calibrated extraction + `find_apex` + secondary query.

use crate::data_sources::reference_library::{
    ExpectedIntensity,
    ScoredIdentity,
};
use crate::data_sources::speclib::Speclib;
use crate::errors::DataProcessingError;
use crate::models::sequence::Peptide;
use crate::{
    ExpectedIntensities,
    IonAnnot,
    ScorerQueriable,
    timed,
};
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timsquery::traits::QueryGeom;
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    KeyLike,
    MzMobilityStatsCollector,
    SpectralCollector,
    TimsElutionGroup,
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
use super::hyperscore::single_lazyscore;
use super::offsets::MzMobilityOffsets;
use super::results::{
    ScoredCandidate,
    ScoredCandidateBuilder,
};
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

    /// Returns `Err(SkipReason::CalibrantNanScore)` when `candidate.score` is
    /// NaN so the caller can surface the count in reports. Non-positive or
    /// +Inf scores are silently discarded — they are either sentinel
    /// (`score <= 0.0` means "no usable signal") or already dominated by any
    /// finite positive candidate in the heap.
    pub fn push(&mut self, candidate: CalibrantCandidate) -> Result<(), SkipReason> {
        if candidate.score.is_nan() {
            return Err(SkipReason::CalibrantNanScore);
        }
        if !candidate.score.is_finite() || candidate.score <= 0.0 {
            return Ok(());
        }
        if self.heap.len() < self.capacity {
            self.heap.push(std::cmp::Reverse(candidate));
        } else if let Some(std::cmp::Reverse(min)) = self.heap.peek()
            && candidate.score > min.score
        {
            self.heap.pop();
            self.heap.push(std::cmp::Reverse(candidate));
        }
        Ok(())
    }

    pub fn merge(mut self, other: Self) -> Self {
        for item in other.heap {
            // Ignore: other heap already filtered on its own push.
            let _ = self.push(item.0);
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
pub const TOP_N_FRAGMENTS: usize = 8;

/// Speclib-only pre-gate: reject peptides whose library entry carries no
/// predicted fragments before doing any extraction work. Shared by the
/// broad (prescore) and calibrated (score_calibrated) paths so the two
/// don't drift.
#[inline]
fn gate_expected_fragments(expected: &ExpectedIntensities<IonAnnot>) -> Result<(), SkipReason> {
    if expected.fragment_len() == 0 {
        return Err(SkipReason::NoExpectedFragments);
    }
    Ok(())
}

/// Fill the per-worker scratch elution group in place from a `RefQuery`
/// flyweight (Task 9). `reset_from` copies the per-variant geometry — for a
/// decoy the fragment m/z values are ALREADY shifted by value, so no extra
/// work is needed. It also sets the precursor labels to the isotope-envelope
/// indices via the flyweight's `iter_precursors` (`0..n_isotopes`), which match
/// `expected_precursor_envelope`'s indices, so no separate label pass is needed.
pub fn fill_scratch_from<Q: QueryGeom<Label = IonAnnot>>(
    dst: &mut TimsElutionGroup<IonAnnot>,
    q: &Q,
) {
    dst.reset_from(q);
}

/// Per-worker scratch buffers for the lazy scoring path: one reusable elution
/// group + one reusable expected-intensities set, both refilled from the
/// flyweight per item so the hot loop stays allocation-free after warm-up.
/// Kept OUTSIDE `ScoringWorker` so the filled scratch can be borrowed
/// immutably while the worker is borrowed mutably by the scoring calls.
pub struct ScratchBufs {
    pub eg: TimsElutionGroup<IonAnnot>,
    pub expected: ExpectedIntensities<IonAnnot>,
}

impl ScratchBufs {
    fn new() -> Self {
        Self {
            eg: TimsElutionGroup::empty_like(),
            expected: ExpectedIntensities::default(),
        }
    }

    /// Refill both buffers from the flyweight: geometry into `eg`, and the
    /// expected fragment/precursor intensities into `expected` (deduped by
    /// key via `try_from_pairs` — library keys are unique by construction).
    fn fill_from<Q: QueryGeom<Label = IonAnnot> + ExpectedIntensity>(&mut self, q: &Q) {
        fill_scratch_from(&mut self.eg, q);
        self.expected = ExpectedIntensities::try_from_pairs(
            q.iter_expected_fragments(),
            q.expected_precursor_envelope().into_iter(),
        )
        .expect("library flyweight yields unique fragment/precursor keys");
    }
}

/// Retain only the top `n` fragments by predicted intensity.
///
/// Removes lower-ranked fragments from the chromatogram collector (fragments array + eg)
/// and from expected intensities, maintaining the invariant that all three agree on count.
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

impl<I: ScorerQueriable> Scorer<I> {
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
        tracing::instrument(skip_all, level = "trace")
    )]
    /// Performs refined secondary query at detected apex with two-pass strategy.
    /// Results populate `worker.inner_collector` and `worker.isotope_collector` in place.
    fn execute_secondary_query(
        &self,
        query: &TimsElutionGroup<IonAnnot>,
        main_score: &ApexScore,
        spectral_tol: &Tolerance,
        isotope_tol: &Tolerance,
        worker: &mut ScoringWorker,
    ) {
        let new_rt_seconds = main_score.retention_time_ms as f32 / 1000.0;

        // **Pass 1**: query at apex RT to determine observed mobility.
        let inner = worker
            .inner_collector
            .get_or_insert_with(|| SpectralCollector::new(query));
        inner.reset_with_overrides(query, Some(new_rt_seconds), None);
        self.index.add_query(inner, spectral_tol);

        let mobility = Self::get_mobility(inner);

        // **Pass 2**: same collector, now with mobility override.
        let inner = worker.inner_collector.as_mut().expect("init above");
        inner.reset_with_overrides(query, Some(new_rt_seconds), Some(mobility as f32));

        // Isotope scratch eg holds `query` with +1 neutron offset applied
        // (buffer-override — reuses Vec capacity after warm-up).
        let scratch_eg = worker
            .isotope_scratch_eg
            .get_or_insert_with(|| query.clone());
        crate::utils::elution_group_ops::apply_isotope_offset_fragments_into(
            scratch_eg, query, 1i8,
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

        let mut candidate = ScoredCandidateBuilder::default()
            .with_metadata(metadata)
            .with_nqueries(nqueries)
            .with_sorted_offsets(&offsets)
            .with_relative_intensities(rel_inten)
            .with_secondary_lazyscores(lazyscores)
            .with_apex_score(main_score)
            .finalize()?;

        // No searchable mobility axis (mzML/FAIMS) → the observed mobility is a
        // sentinel, so drop every mobility feature to NaN (forust-missing).
        if !self.index.mobility_kind().is_scoreable() {
            candidate.scoring.neutralize_mobility();
        }
        Ok(candidate)
    }
}

impl<I: ScorerQueriable> Scorer<I> {
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
        query: &TimsElutionGroup<IonAnnot>,
        expected: &ExpectedIntensities<IonAnnot>,
        digest: Peptide,
        calibration: &CalibrationResult,
        worker: &mut ScoringWorker,
    ) -> Result<super::apex_finding::PeptideMetadata, SkipReason> {
        let original_irt = LibraryRT(query.rt_seconds());
        let calibrated_rt = calibration.convert_irt(original_irt);
        let tolerance = calibration.get_tolerance(
            query.mono_precursor_mz(),
            query.mobility_ook0(),
            original_irt, // library RT — ridge widths are indexed by library RT
        );

        super::extraction::build_extraction_into(
            &mut worker.extraction,
            query,
            Some(calibrated_rt.0),
            expected,
            &self.index,
            &tolerance,
            Some(TOP_N_FRAGMENTS),
        )?;

        let extr = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        Ok(super::apex_finding::PeptideMetadata {
            digest,
            charge: query.precursor_charge(),
            library_id: extr.chromatograms.id as u32,
            library_rt: original_irt.0,
            calibrated_rt_seconds: calibrated_rt.0,
            ref_mobility_ook0: query.mobility_ook0(),
            ref_precursor_mz: query.mono_precursor_mz(),
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
        query: &TimsElutionGroup<IonAnnot>,
        expected: &ExpectedIntensities<IonAnnot>,
        digest: Peptide,
        calibration: &CalibrationResult,
        worker: &mut ScoringWorker,
        timings: &mut ScoreTimings,
    ) -> Result<ScoredCandidate, SkipReason> {
        gate_expected_fragments(expected)?;

        let metadata = timed!(
            timings.extraction,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::extraction").in_scope(|| self
                .build_calibrated_extraction_into(query, expected, digest, calibration, worker))
        )?;

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");

        let apex_score = timed!(
            timings.scoring,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::apex_scoring").in_scope(
                || {
                    worker
                        .scorer
                        .find_apex(scoring_ctx, &|idx| self.map_rt_index_to_milis(idx))
                        .map_err(|e| apex_error_to_skip(&e))
                },
            )
        )?;

        timed!(timings.spectral_query, {
            let spectral_tol = calibration.get_spectral_tolerance();
            let isotope_tol = calibration.get_isotope_tolerance();
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
        let isotope_collector = worker.isotope_collector.as_ref().expect("set by secondary");

        let scoring_ctx = worker
            .extraction
            .as_ref()
            .expect("extraction set by build_extraction_into");
        let nqueries = scoring_ctx.chromatograms.fragments.num_ions() as u8;
        timed!(
            timings.assembly,
            tracing::span!(tracing::Level::TRACE, "score_calibrated::finalize").in_scope(|| {
                self.finalize_results(
                    &metadata,
                    nqueries,
                    &apex_score,
                    inner_collector,
                    isotope_collector,
                )
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
        lib: &Speclib,
        flat_range: std::ops::Range<usize>,
        calibration: &CalibrationResult,
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts) {
        // Single columnar store (Task 9 deleted the materialized arm): the
        // flyweight is always a `RefQuery` from the arena, so the loop is
        // monomorphized over one concrete type — statically dispatched, no
        // per-item heap allocation on the scoring hot path.
        self.score_calibrated_batch_impl(|f| lib.item_at(f), flat_range, calibration)
    }

    fn score_calibrated_batch_impl<Q>(
        &self,
        get_item: impl Fn(usize) -> Q + Sync,
        flat_range: std::ops::Range<usize>,
        calibration: &CalibrationResult,
    ) -> (Vec<ScoredCandidate>, ScoreTimings, SkipCounts)
    where
        Q: QueryGeom<Label = IonAnnot> + ExpectedIntensity + ScoredIdentity,
    {
        let num_cycles = self.num_cycles();
        // Materialize the flat index list so `fold_reduce` (which parallelizes
        // over a slice) can drive the flyweight by index. Each item is a
        // `usize`; the flyweight itself is never stored.
        let flats: Vec<usize> = flat_range.collect();
        // Precursor-range gate over the flyweight geometry.
        let filter_fn = |q: &Q| {
            let tmp = q.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };
        // Pre-scan so scratch capacity holds every peptide — no realloc in hot path.
        let max_frags = flats
            .iter()
            .map(|&f| get_item(f).fragment_count())
            .max()
            .unwrap_or(0);

        // One-shot count of pre-filter rejects. filter_fn is cheap (range check);
        // a second pass keeps the rayon hot loop free of synchronised counters.
        let precursor_oofr_count =
            flats.iter().filter(|&&f| !filter_fn(&get_item(f))).count() as u32;

        let (_worker, _scratch, mut results): (ScoringWorker, ScratchBufs, IonSearchAccumulator) =
            super::maybe_par::fold_reduce(
                &flats,
                || {
                    tracing::debug!(
                        target: "alloc_track",
                        phase = "phase3_score",
                        num_cycles,
                        max_frags,
                        "worker init: ScoringWorker + ScratchBufs + IonSearchAccumulator"
                    );
                    (
                        ScoringWorker::new(num_cycles, max_frags),
                        ScratchBufs::new(),
                        IonSearchAccumulator::default(),
                    )
                },
                |(mut worker, mut scratch, acc), (_idx, &flat)| {
                    let q = get_item(flat);
                    if !filter_fn(&q) {
                        return (worker, scratch, acc);
                    }
                    let _span =
                        tracing::span!(tracing::Level::TRACE, "score_calibrated_item").entered();
                    let mut t = ScoreTimings::default();
                    scratch.fill_from(&q);
                    let digest = q.materialize_peptide();
                    let result = self.score_calibrated_extraction(
                        &scratch.eg,
                        &scratch.expected,
                        digest,
                        calibration,
                        &mut worker,
                        &mut t,
                    );
                    (worker, scratch, acc.fold((result, t)))
                },
                |(wa, sa, a), (_wb, _sb, b)| (wa, sa, a.reduce(b)),
            );

        results.skips.precursor_out_of_fragmented_range = precursor_oofr_count;
        (results.res, results.timings, results.skips)
    }

    /// Phase 1: Lightweight prescore — broad extraction + find_apex_location only.
    /// Returns the apex location (with split product score) and metadata.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn prescore(
        &self,
        query: &TimsElutionGroup<IonAnnot>,
        expected: &ExpectedIntensities<IonAnnot>,
        worker: &mut ScoringWorker,
        timings: &mut PrescoreTimings,
    ) -> Result<ApexLocation, SkipReason> {
        gate_expected_fragments(expected)?;

        timed!(
            timings.extraction,
            tracing::span!(tracing::Level::TRACE, "prescore::extraction").in_scope(|| {
                super::extraction::build_extraction_into(
                    &mut worker.extraction,
                    query,
                    None,
                    expected,
                    &self.index,
                    &self.broad_tolerance,
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

    /// Phase 1 batch: Prescore all peptides, collecting top-N calibrant candidates via bounded heaps.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    pub fn prescore_batch(
        &self,
        lib: &Speclib,
        flat_range: std::ops::Range<usize>,
        config: &CalibrationConfig,
        timings: &mut PrescoreTimings,
    ) -> CalibrantHeap {
        // Single columnar store (Task 9): iterate `RefQuery` flyweights from
        // the arena directly — monomorphized, no per-item heap alloc on the
        // prescore hot path (see `score_calibrated_batch`).
        self.prescore_batch_impl(|f| lib.item_at(f), flat_range, config, timings)
    }

    fn prescore_batch_impl<Q>(
        &self,
        get_item: impl Fn(usize) -> Q + Sync,
        flat_range: std::ops::Range<usize>,
        config: &CalibrationConfig,
        timings: &mut PrescoreTimings,
    ) -> CalibrantHeap
    where
        Q: QueryGeom<Label = IonAnnot> + ExpectedIntensity + ScoredIdentity,
    {
        // The flat index IS the global speclib index (see `Speclib::item_at`),
        // so it doubles as `CalibrantCandidate::speclib_index` — no separate
        // chunk offset needed.
        let flats: Vec<usize> = flat_range.collect();
        let filter_fn = |q: &Q| {
            let tmp = q.get_precursor_mz_limits();
            let lims = TupleRange::try_new(tmp.0, tmp.1).expect("Should already be ordered");
            self.fragmented_range.intersects(lims)
        };

        // Count pre-filter rejects once; keeps the rayon loop counter-free.
        let precursor_oofr_count =
            flats.iter().filter(|&&f| !filter_fn(&get_item(f))).count() as u32;
        timings.skips.precursor_out_of_fragmented_range += precursor_oofr_count;

        let num_cycles = self.num_cycles();
        let n_calibrants = config.n_calibrants;
        let max_frags = flats
            .iter()
            .map(|&f| get_item(f).fragment_count())
            .max()
            .unwrap_or(0);

        let (_worker, _scratch, heap, par_timings) = super::maybe_par::fold_reduce(
            &flats,
            || {
                tracing::debug!(
                    target: "alloc_track",
                    phase = "phase1_prescore",
                    num_cycles,
                    n_calibrants,
                    max_frags,
                    "worker init: ScoringWorker + ScratchBufs + CalibrantHeap"
                );
                (
                    ScoringWorker::new(num_cycles, max_frags),
                    ScratchBufs::new(),
                    CalibrantHeap::new(n_calibrants),
                    PrescoreTimings::default(),
                )
            },
            |(mut worker, mut scratch, mut heap, mut t), (_idx, &flat)| {
                let q = get_item(flat);
                if !filter_fn(&q) {
                    return (worker, scratch, heap, t);
                }
                t.n_passed_filter += 1;
                scratch.fill_from(&q);
                match self.prescore(&scratch.eg, &scratch.expected, &mut worker, &mut t) {
                    Ok(loc) => {
                        let cand = CalibrantCandidate {
                            score: loc.score,
                            apex_rt: ObservedRTSeconds(loc.retention_time_ms as f32 / 1000.0),
                            speclib_index: flat,
                            library_rt: LibraryRT(q.rt_seconds()),
                        };
                        if let Err(reason) = heap.push(cand) {
                            t.skips.bump(reason);
                        }
                    }
                    Err(reason) => t.skips.bump(reason),
                }
                (worker, scratch, heap, t)
            },
            |(wa, sa, ha, mut ta), (_wb, _sb, hb, tb)| {
                ta += tb;
                (wa, sa, ha.merge(hb), ta)
            },
        );
        *timings += par_timings;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_sources::reference_library::{
        ExpectedIntensity,
        ReferenceLibrary,
    };
    use timsquery::TimsElutionGroup;
    use timsquery::models::QueryCollection;
    use timsquery::models::capabilities::LibCapabilities;

    fn tiny_lazy_lib() -> ReferenceLibrary {
        let mut geom = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        geom.push_target(
            900.4,
            2,
            1.0,
            1.0,
            &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y8").unwrap(), 800.0),
            ],
            "PEPTIDEK",
            "PEPTIDEK",
            &[],
        );
        geom.seal();
        ReferenceLibrary {
            geom,
            frag_intens: vec![1.0, 0.5],
        }
    }

    #[test]
    fn scratch_eg_filled_from_flyweight_matches_geometry() {
        let lib = tiny_lazy_lib();
        // Variant 1 (+decoy): geometry is mass-shifted, so this exercises the
        // by-value shifted-fragment path through reset_from.
        let q = lib.item_at(1);
        let mut scratch = TimsElutionGroup::<IonAnnot>::empty_like();
        fill_scratch_from(&mut scratch, &q);

        assert!((scratch.mono_precursor_mz() - q.mono_precursor_mz()).abs() < 1e-9);

        let a: Vec<(IonAnnot, f64)> = scratch
            .iter_fragments_refs()
            .map(|(l, m)| (*l, *m))
            .collect();
        let b: Vec<(IonAnnot, f64)> = q.iter_fragments_refs().map(|(l, m)| (*l, m)).collect();
        assert_eq!(a, b);

        // Precursor labels come from the isotope envelope indices.
        let env = q.expected_precursor_envelope();
        let labels: Vec<i8> = scratch.iter_precursors().map(|(iso, _)| iso).collect();
        let expected_labels: Vec<i8> = env.iter().map(|(i, _)| *i).collect();
        assert_eq!(labels, expected_labels);
    }
}
