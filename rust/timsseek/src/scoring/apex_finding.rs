//! Logic for finding the apex of a peptide elution profile.
//!
//! This module replaces the old `LocalizationBuffer` and `calculate_scores.rs`.
//! It implements an efficient, accumulator-based scoring engine that avoids
//! intermediate data transposition.
//!
//! # Usage
//!
//! ```ignore
//! use timsseek::scoring::apex_finding::{ApexFinder, CandidateContext};
//!
//! // 1. Create a reusable finder (one per thread)
//! let mut finder = ApexFinder::new(chromatogram_collector.num_cycles());
//!
//! // 2. Create the context for a specific query
//! let context = CandidateContext {
//!     digest: digest_slice,
//!     charge: 2,
//!     expected_intensities: expected,
//!     query_values: chromatogram_collector,
//! };
//!
//! // 3. Score (reusing the finder's internal buffers)
//! let score = finder.find_apex(&context, rt_mapping_fn).unwrap();
//! println!("Found apex at RT: {} ms with score {}", score.retention_time_ms, score.score);
//! ```

use std::fmt::Display;

use super::{
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};
use crate::IonAnnot;
use crate::errors::DataProcessingError;
use crate::models::{
    DigestSlice,
    ExpectedIntensities,
};
use crate::scoring::scores::apex_features::{
    ApexFeatures,
    SplitProductScore,
    compute_apex_features,
    compute_split_product,
    compute_weighted_score,
    find_joint_apex,
};
use crate::scoring::scores::scribe::SCRIBE_FLOOR;
use crate::utils::top_n_array::TopNArray;
use serde::Serialize;
use timsquery::models::aggregators::ChromatogramCollector;
use timsquery::traits::KeyLike;
use timsquery::{
    MzMobilityStatsCollector,
    SpectralCollector,
};

/// Represents a peptide candidate context required for scoring.
///
/// This bundles the theoretical peptide info (digest, expected intensities)
/// with the raw extracted data (`ChromatogramCollector`).
///
/// The label is in essence anything that identifies the peptide sequence and modifications,
/// so it can be a string, or a more complex struct. (or a simpler one like a smart pointer)
#[derive(Debug)]
pub struct CandidateContext<T: KeyLike, L: Display> {
    /// The peptide sequence and modification information.
    pub label: L,
    /// The precursor charge state.
    pub charge: u8,
    /// The expected theoretical intensities of precursor and fragment ions.
    pub expected_intensities: ExpectedIntensities<T>,
    /// The observed chromatogram data collected from the instrument.
    pub query_values: ChromatogramCollector<T, f32>,
}

/// Immutable peptide metadata that never changes during scoring.
///
/// This separates identity/reference information from the scoring data,
/// making the scoring pipeline cleaner and more efficient.
#[derive(Debug, Clone)]
pub struct PeptideMetadata {
    /// The peptide sequence and modification information.
    pub digest: DigestSlice,

    /// The precursor charge state.
    pub charge: u8,

    /// Library identifier for this peptide.
    pub library_id: u32,

    /// Reference retention time (seconds) from library.
    pub ref_rt_seconds: f32,

    /// Reference ion mobility (ook0) from library.
    pub ref_mobility_ook0: f32,

    /// Reference precursor m/z from library.
    pub ref_precursor_mz: f64,
}

/// Scoring-relevant data that flows through the apex finding pipeline.
///
/// This contains only the data needed for scoring calculations,
/// separated from metadata for clarity and efficiency.
#[derive(Debug)]
pub struct ScoringContext<T: KeyLike> {
    /// The expected theoretical intensities of precursor and fragment ions.
    pub expected_intensities: ExpectedIntensities<T>,

    /// The observed chromatogram data collected from the instrument.
    pub query_values: ChromatogramCollector<T, f32>,
}

/// Lightweight result from Phase 1 apex finding.
/// Contains only the apex location and a basic score for calibrant ranking.
#[derive(Debug, Clone, Copy)]
pub struct ApexLocation {
    /// Basic score (apex profile peak value) for ranking.
    pub score: f32,
    /// Retention time at the apex (ms).
    pub retention_time_ms: u32,
    /// Local cycle index of the apex within the extraction.
    pub apex_cycle: usize,
    /// Peak shape metrics for baseline computation.
    pub raising_cycles: u8,
    pub falling_cycles: u8,
}

/// The result of the full scoring process (Phase 3).
#[derive(Debug, Clone, Copy)]
pub struct ApexScore {
    /// The final weighted product score (higher is better).
    pub score: f32,
    /// Retention time at the apex (ms).
    pub retention_time_ms: u32,
    /// Local cycle index of the joint apex.
    pub joint_apex_cycle: usize,

    // --- Split product components ---
    pub split_product: SplitProductScore,

    // --- 11 features at joint apex ---
    pub features: ApexFeatures,

    // --- Peak discrimination ---
    pub delta_next: f32,
    pub delta_second_next: f32,

    // --- Retained from current (used downstream) ---
    pub lazyscore: f32,
    pub lazyscore_vs_baseline: f32,
    pub lazyscore_z: f32,
    pub npeaks: u8,
    pub ms2_summed_intensity: f32,
    pub ms1_summed_intensity: f32,
    pub raising_cycles: u8,
    pub falling_cycles: u8,
}

/// Stores time-resolved scores for every cycle in the chromatogram.
#[derive(Debug, Clone, Serialize)]
pub struct ScoreTraces {
    /// Per-cycle cosine similarity (sqrt-transformed expected).
    pub ms2_cosine_ref_sim: Vec<f32>,
    /// Per-cycle lazyscore (kept for baseline lambda computation).
    pub ms2_lazyscore: Vec<f32>,
    /// Per-cycle Scribe score.
    pub ms2_scribe: Vec<f32>,
    /// Per-cycle log1p(sum(fragment_intensities)).
    pub ms2_log_intensity: Vec<f32>,
    /// Per-cycle summed precursor intensity (keys >= 0 only).
    pub ms1_precursor_trace: Vec<f32>,
    /// Composite apex profile for peak picking.
    pub main_score: Vec<f32>,
}

impl ScoreTraces {
    pub fn new_with_capacity(capacity: usize) -> Self {
        Self {
            ms2_cosine_ref_sim: Vec::with_capacity(capacity),
            ms2_lazyscore: Vec::with_capacity(capacity),
            ms2_scribe: Vec::with_capacity(capacity),
            ms2_log_intensity: Vec::with_capacity(capacity),
            ms1_precursor_trace: Vec::with_capacity(capacity),
            main_score: Vec::with_capacity(capacity),
        }
    }

    pub fn clear(&mut self) {
        self.ms2_cosine_ref_sim.clear();
        self.ms2_lazyscore.clear();
        self.ms2_scribe.clear();
        self.ms2_log_intensity.clear();
        self.ms1_precursor_trace.clear();
        self.main_score.clear();
    }

    /// Resize all buffers to the specified length (filling with 0.0).
    pub fn resize(&mut self, len: usize) {
        self.ms2_cosine_ref_sim.resize(len, 0.0);
        self.ms2_lazyscore.resize(len, 0.0);
        self.ms2_scribe.resize(len, 0.0);
        self.ms2_log_intensity.resize(len, 0.0);
        self.ms1_precursor_trace.resize(len, 0.0);
        self.main_score.resize(len, 0.0);
    }

    pub fn iter_scores(&self) -> impl Iterator<Item = (&'static str, &[f32])> + '_ {
        vec![
            ("ms2_cosine_ref_sim", &self.ms2_cosine_ref_sim[..]),
            ("ms2_lazyscore", &self.ms2_lazyscore[..]),
            ("ms2_scribe", &self.ms2_scribe[..]),
            ("ms2_log_intensity", &self.ms2_log_intensity[..]),
            ("ms1_precursor_trace", &self.ms1_precursor_trace[..]),
            ("main_score", &self.main_score[..]),
        ]
        .into_iter()
    }
}

/// The core engine for finding peptide apexes.
#[derive(Debug)]
pub struct ApexFinder {
    pub traces: ScoreTraces,
    buffers: ApexFinderBuffers,
}

#[derive(Debug)]
struct ApexFinderBuffers {
    /// Cosine numerator: sum(obs * sqrt(exp)) per cycle.
    temp_ms2_dot_prod: Vec<f32>,
    /// Cosine denominator: sum(obs^2) per cycle.
    temp_ms2_norm_sq_obs: Vec<f32>,
    /// Scribe: sum(sqrt(obs)) per cycle.
    temp_sqrt_sum: Vec<f32>,
    /// Log-intensity: sum(obs) per cycle (finalized as log1p).
    temp_raw_intensity_sum: Vec<f32>,
}

impl ApexFinderBuffers {
    fn new(size: usize) -> Self {
        Self {
            temp_ms2_dot_prod: vec![0.0f32; size],
            temp_ms2_norm_sq_obs: vec![0.0f32; size],
            temp_sqrt_sum: vec![0.0f32; size],
            temp_raw_intensity_sum: vec![0.0f32; size],
        }
    }

    fn clear(&mut self) {
        self.temp_ms2_dot_prod.fill(0.0);
        self.temp_ms2_norm_sq_obs.fill(0.0);
        self.temp_sqrt_sum.fill(0.0);
        self.temp_raw_intensity_sum.fill(0.0);
    }

    fn resize(&mut self, len: usize) {
        self.temp_ms2_dot_prod.resize(len, 0.0);
        self.temp_ms2_norm_sq_obs.resize(len, 0.0);
        self.temp_sqrt_sum.resize(len, 0.0);
        self.temp_raw_intensity_sum.resize(len, 0.0);
    }
}

impl ApexFinder {
    pub fn new(capacity: usize) -> Self {
        Self {
            traces: ScoreTraces::new_with_capacity(capacity),
            buffers: ApexFinderBuffers::new(capacity),
        }
    }

    /// Build cosine and scribe profiles from traces.
    /// cosine_profile[i] = cos^3 * intensity, scribe_profile[i] = scribe * intensity.
    fn build_profiles(&self) -> (Vec<f32>, Vec<f32>) {
        let n = self.traces.ms2_cosine_ref_sim.len();
        let mut cosine_profile = Vec::with_capacity(n);
        let mut scribe_profile = Vec::with_capacity(n);
        for i in 0..n {
            let cos = self.traces.ms2_cosine_ref_sim[i];
            let intensity = self.traces.ms2_log_intensity[i];
            cosine_profile.push(cos * cos * cos * intensity);
            scribe_profile.push(self.traces.ms2_scribe[i] * intensity);
        }
        (cosine_profile, scribe_profile)
    }

    /// Phase 1: Find apex location using broad extraction.
    ///
    /// Returns a lightweight `ApexLocation` with just the peak location and
    /// a basic score (apex profile value). Sufficient for calibrant ranking.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, scoring_ctx, rt_mapper), level = "trace")
    )]
    pub fn find_apex_location<T: KeyLike>(
        &mut self,
        scoring_ctx: &ScoringContext<T>,
        rt_mapper: &dyn Fn(usize) -> u32,
    ) -> Result<ApexLocation, DataProcessingError> {
        let collector = &scoring_ctx.query_values;
        let n_cycles = collector.num_cycles();

        self.traces.clear();
        self.traces.resize(n_cycles);
        self.buffers.clear();
        self.buffers.resize(n_cycles);

        self.compute_pass_1(scoring_ctx)?;
        self.compute_main_score_trace();

        // Peak-pick on apex profile
        let peak_picker = PeakPicker::new(&self.traces.main_score);
        let (max_val, max_loc) = match peak_picker.next_peak() {
            Some(p) => p,
            None => {
                return Err(DataProcessingError::ExpectedNonEmptyData {
                    context: Some("No main score found".into()),
                });
            }
        };
        if max_val == 0.0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("No non-0 main score".into()),
            });
        }

        let (raising_cycles, falling_cycles) = self.calculate_rise_and_fall_cycles(max_loc);
        let cycle_offset = scoring_ctx.query_values.cycle_offset();
        let retention_time_ms = rt_mapper(max_loc + cycle_offset);

        // Compute split product score for calibrant ranking
        let (cosine_profile, scribe_profile) = self.build_profiles();

        let split_product = compute_split_product(
            &cosine_profile,
            &scribe_profile,
            &scoring_ctx.query_values.fragments,
            &scoring_ctx.expected_intensities.fragment_intensities,
        );

        Ok(ApexLocation {
            score: split_product.base_score,
            retention_time_ms,
            apex_cycle: max_loc,
            raising_cycles,
            falling_cycles,
        })
    }

    /// Phase 3: Full scoring on a (narrow) extraction.
    ///
    /// Computes traces, apex profile, split product, 11 features, and weighted score.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip(self, scoring_ctx, rt_mapper), level = "trace")
    )]
    pub fn find_apex<T: KeyLike>(
        &mut self,
        scoring_ctx: &ScoringContext<T>,
        rt_mapper: &dyn Fn(usize) -> u32,
    ) -> Result<ApexScore, DataProcessingError> {
        let collector = &scoring_ctx.query_values;
        let n_cycles = collector.num_cycles();

        // 1. Reset buffers
        self.traces.clear();
        self.traces.resize(n_cycles);
        self.buffers.clear();
        self.buffers.resize(n_cycles);

        // 2. Compute per-cycle scores (single pass)
        self.compute_pass_1(scoring_ctx)?;

        // 3. Compute apex profile (cos^3 * I combined with scribe * I)
        self.compute_main_score_trace();

        // 4. Find apex and extract features
        self.extract_apex_score(scoring_ctx, &rt_mapper)
    }

    /// Single-pass scoring: cosine (sqrt-transformed), scribe, lazyscore,
    /// log-intensity, and precursor trace.
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn compute_pass_1<T: KeyLike>(
        &mut self,
        scoring_ctx: &ScoringContext<T>,
    ) -> Result<(), DataProcessingError> {
        let collector = &scoring_ctx.query_values;

        // --- MS2 (Fragments) ---
        let ms2_dot_prod = &mut self.buffers.temp_ms2_dot_prod;
        let ms2_norm_sq_obs = &mut self.buffers.temp_ms2_norm_sq_obs;
        let sqrt_sum = &mut self.buffers.temp_sqrt_sum;
        let raw_sum = &mut self.buffers.temp_raw_intensity_sum;
        // Sum of sqrt(expected) for cosine norm: ||sqrt(exp)||^2 = sum(exp)
        let mut ms2_sum_exp = 0.0f32;

        // Pre-compute pred_norm for scribe
        let mut pred_norms: Vec<(usize, f32)> = Vec::new();
        let mut pred_sqrt_sum = 0.0f32;

        for (row_idx, ((key, _mz), chrom)) in collector.fragments.iter_mzs().enumerate() {
            let expected = scoring_ctx
                .expected_intensities
                .fragment_intensities
                .get(key)
                .copied()
                .unwrap_or(0.0);

            if expected <= 0.0 {
                continue;
            }

            let sqrt_exp = expected.sqrt();
            ms2_sum_exp += expected; // sqrt_exp * sqrt_exp = expected
            pred_norms.push((row_idx, expected.sqrt()));
            pred_sqrt_sum += sqrt_exp;

            for (i, &intensity) in chrom.iter().enumerate() {
                if intensity > 0.0 {
                    // Lazyscore accumulation
                    self.traces.ms2_lazyscore[i] += intensity.max(1.0).ln();
                    // Cosine: dot(obs, sqrt(exp))
                    ms2_dot_prod[i] += intensity * sqrt_exp;
                    // Cosine: obs norm
                    ms2_norm_sq_obs[i] += intensity * intensity;
                    // Scribe: sqrt(obs) sum
                    sqrt_sum[i] += intensity.sqrt();
                }
                // Raw intensity sum (for log-intensity, includes zeros)
                raw_sum[i] += intensity.max(0.0);
            }
        }

        // Finalize cosine, lazyscore, log-intensity
        let norm_sqrt_exp = ms2_sum_exp.sqrt(); // ||sqrt(exp)|| = sqrt(sum(exp))
        let n = self.traces.ms2_cosine_ref_sim.len();
        for i in 0..n {
            // Lazyscore
            self.traces.ms2_lazyscore[i] =
                crate::utils::math::lnfact_f32(self.traces.ms2_lazyscore[i]);

            // Cosine (sqrt-transformed expected)
            let obs_norm = ms2_norm_sq_obs[i].sqrt();
            if obs_norm > 0.0 && norm_sqrt_exp > 0.0 {
                self.traces.ms2_cosine_ref_sim[i] =
                    (ms2_dot_prod[i] / (obs_norm * norm_sqrt_exp)).clamp(1e-3, 1.0);
            } else {
                self.traces.ms2_cosine_ref_sim[i] = 1e-3;
            }

            // Log-intensity
            self.traces.ms2_log_intensity[i] = raw_sum[i].ln_1p();
        }

        // Finalize scribe (inline, reusing sqrt_sum buffer from accumulation above)
        if pred_sqrt_sum > 0.0 && !pred_norms.is_empty() {
            // Normalize pred_norms in-place
            for entry in pred_norms.iter_mut() {
                entry.1 /= pred_sqrt_sum;
            }

            // Pass B: accumulate SSE
            for &(row_idx, pred_norm_i) in &pred_norms {
                let row = collector
                    .fragments
                    .get_row_idx(row_idx)
                    .expect("row_idx from enumeration must be valid");
                for (t, &intensity) in row.iter().enumerate() {
                    if sqrt_sum[t] == 0.0 {
                        continue;
                    }
                    let obs_norm_i = if intensity > 0.0 {
                        intensity.sqrt() / sqrt_sum[t]
                    } else {
                        0.0
                    };
                    let diff = obs_norm_i - pred_norm_i;
                    self.traces.ms2_scribe[t] += diff * diff;
                }
            }

            // Finalize scribe: -log(sse)
            for t in 0..n {
                if sqrt_sum[t] == 0.0 {
                    self.traces.ms2_scribe[t] = SCRIBE_FLOOR;
                } else {
                    let sse = self.traces.ms2_scribe[t].max(f32::EPSILON);
                    self.traces.ms2_scribe[t] = -sse.ln();
                }
            }
        } else {
            self.traces.ms2_scribe.fill(SCRIBE_FLOOR);
        }

        // --- MS1 Precursor trace ---
        for ((key, _mz), chrom) in collector.precursors.iter_mzs() {
            if *key < 0 {
                continue; // Skip decoy isotope keys
            }
            for (i, &intensity) in chrom.iter().enumerate() {
                if intensity > 0.0 {
                    self.traces.ms1_precursor_trace[i] += intensity;
                }
            }
        }

        Ok(())
    }

    /// Compute the apex profile from cosine and scribe traces.
    ///
    /// apex_profile(t) = C(t) * (0.5 + S_norm(t))
    /// where C(t) = cosine(t)^3 * I(t)
    ///       S(t) = scribe(t) * I(t)
    ///       S_norm = (S - min(S)) / (max(S) - min(S))
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn compute_main_score_trace(&mut self) {
        let len = self.traces.ms2_cosine_ref_sim.len();
        self.traces.main_score.clear();
        self.traces.main_score.reserve(len);

        // Compute S(t) = scribe(t) * I(t), find min/max for normalization
        let mut s_min = f32::INFINITY;
        let mut s_max = f32::NEG_INFINITY;

        // Temp: compute S values inline
        for i in 0..len {
            let s = self.traces.ms2_scribe[i] * self.traces.ms2_log_intensity[i];
            s_min = s_min.min(s);
            s_max = s_max.max(s);
        }

        let s_range = s_max - s_min;

        for i in 0..len {
            let cos = self.traces.ms2_cosine_ref_sim[i];
            let intensity = self.traces.ms2_log_intensity[i];
            let c = cos * cos * cos * intensity; // cos^3 * I

            let s = self.traces.ms2_scribe[i] * intensity;
            let s_norm = if s_range > 0.0 {
                (s - s_min) / s_range
            } else {
                0.5 // Degrade to cosine-only when scribe is constant
            };

            self.traces.main_score.push(c * (0.5 + s_norm));
        }
    }

    fn extract_apex_score<T: KeyLike>(
        &self,
        scoring_ctx: &ScoringContext<T>,
        rt_mapper: &dyn Fn(usize) -> u32,
    ) -> Result<ApexScore, DataProcessingError> {
        let mut peak_picker = PeakPicker::new(&self.traces.main_score);

        // Find best peak
        let (max_val, max_loc) = match peak_picker.next_peak() {
            Some(p) => p,
            None => {
                return Err(DataProcessingError::ExpectedNonEmptyData {
                    context: Some("No main score found".into()),
                });
            }
        };

        if max_val == 0.0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("No non-0 main score".into()),
            });
        }

        // Peak shape (rise/fall) for delta computation
        let (raising_cycles, falling_cycles) = self.calculate_rise_and_fall_cycles(max_loc);

        // Mask and find next peaks for delta scores
        peak_picker.mask_peak(max_loc, raising_cycles as usize, falling_cycles as usize, 2);
        let (next_val, next_loc) = peak_picker.next_peak().unwrap_or((0.0, max_loc));
        let (next_raise, next_fall) = self.calculate_rise_and_fall_cycles(next_loc);
        peak_picker.mask_peak(next_loc, next_raise as usize, next_fall as usize, 1);
        let (second_next_val, _) = peak_picker.next_peak().unwrap_or((0.0, max_loc));

        let delta_next = max_val - next_val;
        let delta_second_next = max_val - second_next_val;

        // Build intermediate profiles for split product and features
        let (cosine_profile, scribe_profile) = self.build_profiles();

        // Split product score (independent cosine/scribe apexes)
        let split_product = compute_split_product(
            &cosine_profile,
            &scribe_profile,
            &scoring_ctx.query_values.fragments,
            &scoring_ctx.expected_intensities.fragment_intensities,
        );

        // Joint precursor-fragment apex
        let joint_apex = find_joint_apex(&cosine_profile, &self.traces.ms1_precursor_trace);

        // 11 features at joint apex
        let n_cycles = cosine_profile.len();
        let features = compute_apex_features(
            &scoring_ctx.query_values.fragments,
            &scoring_ctx.query_values.precursors,
            &scoring_ctx.expected_intensities,
            &cosine_profile,
            &self.traces.ms1_precursor_trace,
            joint_apex,
            n_cycles,
        );

        // Weighted final score
        let score = compute_weighted_score(split_product.base_score, &features);

        // RT at joint apex
        let cycle_offset = scoring_ctx.query_values.cycle_offset();
        let global_loc = joint_apex + cycle_offset;
        let retention_time_ms = rt_mapper(global_loc);

        // Intensity counts at joint apex
        let (ms1_summed_intensity, _) =
            self.sum_intensities_at(&scoring_ctx.query_values.precursors, joint_apex);
        let (ms2_summed_intensity, ms2_npeaks) =
            self.sum_intensities_at(&scoring_ctx.query_values.fragments, joint_apex);

        // Lazyscore baseline stats
        let lambda = self.calculate_baseline_lambda(max_loc, raising_cycles, falling_cycles);
        let k = self.traces.ms2_lazyscore[joint_apex] as f64;
        let norm_lazy_std = lambda.sqrt().max(1.0) as f32;
        let lazyscore_z = self.traces.ms2_lazyscore[joint_apex] / norm_lazy_std;

        if lazyscore_z.is_nan() {
            return Err(DataProcessingError::ExpectedFiniteNonNanData {
                context: format!("Lazy score is NaN {} and {}", k, norm_lazy_std),
            });
        }

        Ok(ApexScore {
            score,
            retention_time_ms,
            joint_apex_cycle: joint_apex,
            split_product,
            features,
            delta_next,
            delta_second_next,
            lazyscore: self.traces.ms2_lazyscore[joint_apex],
            lazyscore_vs_baseline: (k / lambda) as f32,
            lazyscore_z,
            npeaks: ms2_npeaks as u8,
            ms2_summed_intensity,
            ms1_summed_intensity,
            raising_cycles,
            falling_cycles,
        })
    }

    fn calculate_rise_and_fall_cycles(&self, max_loc: usize) -> (u8, u8) {
        let score_slice = &self.traces.ms2_lazyscore;
        let raising = count_falling_steps(max_loc, -1, score_slice);
        let falling = count_falling_steps(max_loc, 1, score_slice);
        (raising, falling)
    }

    fn calculate_baseline_lambda(&self, apex: usize, raise: u8, fall: u8) -> f64 {
        let max_width = self.traces.ms2_lazyscore.len();

        let peak_start = apex.saturating_sub((raise as usize) * 2);
        let peak_end = (apex + (fall as usize) * 2 + 1).min(max_width);

        let base_start = apex.saturating_sub((raise as usize) * 4);
        let base_end = (apex + (fall as usize) * 4 + 1).min(max_width);

        // Sum values in baseline window BUT NOT in peak window
        let mut sum = 0.0;
        let mut count = 0;

        for i in base_start..base_end {
            if i >= peak_start && i < peak_end {
                continue;
            }
            let val = self.traces.ms2_lazyscore[i];
            if !val.is_nan() {
                sum += val as f64;
            }
            count += 1;
        }

        if count == 0 {
            1.0
        } else {
            (sum / count as f64).max(1.0)
        }
    }

    // Helper to sum intensity at a specific time point across all ions in a collector
    fn sum_intensities_at<K: KeyLike>(
        &self,
        mz_array: &timsquery::models::MzMajorIntensityArray<K, f32>,
        idx: usize,
    ) -> (f32, usize) {
        let mut sum = 0.0;
        let mut count = 0;
        for (_, chrom) in mz_array.iter_mzs() {
            // Access the value at `idx`.
            // Assumes Chromatogram supports indexing or we treat it as slice.
            if idx < chrom.len() {
                let val = chrom[idx];
                if val > 0.0 {
                    sum += val;
                    count += 1;
                }
            }
        }
        (sum, count)
    }
}

/// Represents the relative intensities of the MS1 and MS2
/// ions, relative to total intensity of the respective MS level.
#[derive(Debug, Clone, Copy)]
pub struct RelativeIntensities {
    pub ms1: TopNArray<NUM_MS1_IONS, f32>,
    pub ms2: TopNArray<NUM_MS2_IONS, f32>,
}

impl RelativeIntensities {
    pub fn new(agg: &SpectralCollector<IonAnnot, MzMobilityStatsCollector>) -> Self {
        let mut ms1: TopNArray<NUM_MS1_IONS, f32> = TopNArray::new();
        let mut ms2: TopNArray<NUM_MS2_IONS, f32> = TopNArray::new();

        let tot_l1p_ms1: f64 = agg
            .iter_precursors()
            .map(|(_k, v)| v.weight())
            .sum::<f64>()
            .ln_1p();
        let tot_l1p_ms2: f64 = agg
            .iter_fragments()
            .map(|(_k, v)| v.weight())
            .sum::<f64>()
            .ln_1p();

        agg.iter_precursors().for_each(|((k, _mz), v)| {
            // Note isotope keys < 0 mean 'decoy' isotopes
            let weight = v.weight();
            if k >= 0i8 && weight > 0.0 {
                ms1.push((weight.ln_1p() - tot_l1p_ms1) as f32);
            }
        });

        agg.iter_fragments().for_each(|(_k, v)| {
            let weight = v.weight();
            if weight > 0.0 {
                ms2.push((weight.ln_1p() - tot_l1p_ms2) as f32);
            }
        });
        Self { ms1, ms2 }
    }
}

// --- Helpers ---

struct PeakPicker<'a> {
    scores: &'a [f32],
    exclusions: Vec<std::ops::Range<usize>>,
}

impl<'a> PeakPicker<'a> {
    fn new(scores: &'a [f32]) -> Self {
        Self {
            scores,
            exclusions: Vec::new(),
        }
    }

    fn next_peak(&self) -> Option<(f32, usize)> {
        let mut best_val = -f32::INFINITY;
        let mut best_idx = None;

        for (i, &val) in self.scores.iter().enumerate() {
            if val.is_nan() {
                continue;
            }

            // Check exclusions
            let is_excluded = self.exclusions.iter().any(|r| r.contains(&i));
            if is_excluded {
                continue;
            }

            if val > best_val {
                best_val = val;
                best_idx = Some(i);
            }
        }

        best_idx.map(|i| (best_val, i))
    }

    fn mask_peak(&mut self, apex: usize, raise: usize, fall: usize, factor: usize) {
        let start = apex.saturating_sub(raise * factor);
        let end = apex + (fall * factor) + 1;
        // Clamp to len is handled by checks, but good to be safe
        let end = end.min(self.scores.len());
        self.exclusions.push(start..end);
    }
}

fn count_falling_steps(start: usize, step: i32, slc: &[f32]) -> u8 {
    const MAX_WIDTH: u8 = 10;
    let mut count = 0;
    let edge = slc.len() as i32;
    let mut last = slc[start];
    while count < MAX_WIDTH {
        let next = start as i32 + (step * count as i32);
        if next < 0 || next >= edge {
            break;
        }
        let next = next as usize;
        let score = slc[next];
        if score > last {
            break;
        }
        last = score;
        count += 1;
    }
    count
}

/// Applies a 1D Gaussian blur in-place.
pub fn gaussblur_in_place(x: &mut [f32]) {
    let len = x.len();
    if len < 3 {
        return;
    }
    const W_SIDE: f32 = 0.5;
    const W_CENTER: f32 = 1.0;
    const NORM: f32 = 2.0;

    let mut prev_val = x[0];
    x[0] = (x[0] * 1.5 + x[1] * 0.5) / NORM;

    for i in 1..len - 1 {
        let current_val = x[i];
        x[i] = (prev_val * W_SIDE + current_val * W_CENTER + x[i + 1] * W_SIDE) / NORM;
        prev_val = current_val;
    }
    x[len - 1] = (x[len - 1] * 1.5 + prev_val * 0.5) / NORM;
}
