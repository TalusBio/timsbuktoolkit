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
    COELUTION_WINDOW_WIDTH,
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};
use crate::IonAnnot;
use crate::errors::DataProcessingError;
use crate::models::{
    DigestSlice,
    ExpectedIntensities,
};
use crate::scoring::scores::coelution::coelution_score::coelution_vref_score_filter_into;
use crate::scoring::scores::corr_v_ref;
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

/// The result of the apex finding process.
#[derive(Debug, Clone, Copy)]
pub struct ApexScore {
    /// The main composite score (higher is better).
    pub score: f32,
    /// Difference to the next best peak.
    pub delta_next: f32,
    /// Difference to the third best peak.
    pub delta_second_next: f32,
    /// Retention time at the apex (ms).
    pub retention_time_ms: u32,

    // --- MS2 Features ---
    pub ms2_cosine_ref_sim: f32,
    pub ms2_coelution_score: f32,
    pub ms2_summed_intensity: f32,
    pub npeaks: u8,
    pub lazyscore: f32,
    pub lazyscore_vs_baseline: f32,
    pub lazyscore_z: f32,
    pub ms2_corr_v_gauss: f32,

    // --- MS1 Features ---
    pub ms1_corr_v_gauss: f32,
    pub ms1_cosine_ref_sim: f32,
    pub ms1_coelution_score: f32,
    pub ms1_summed_intensity: f32,

    // --- Shape Features ---
    pub raising_cycles: u8,
    pub falling_cycles: u8,
}

/// Stores time-resolved scores for every cycle in the chromatogram.
#[derive(Debug, Clone, Serialize)]
pub struct ScoreTraces {
    pub ms1_cosine_ref_sim: Vec<f32>,
    pub ms1_coelution_score: Vec<f32>,
    pub ms1_corr_v_gauss: Vec<f32>,
    pub ms2_cosine_ref_sim: Vec<f32>,
    pub ms2_coelution_score: Vec<f32>,
    pub ms2_lazyscore: Vec<f32>,
    pub ms2_corr_v_gauss: Vec<f32>,
    pub main_score: Vec<f32>,
}

impl ScoreTraces {
    pub fn new_with_capacity(capacity: usize) -> Self {
        Self {
            ms1_cosine_ref_sim: Vec::with_capacity(capacity),
            ms1_coelution_score: Vec::with_capacity(capacity),
            ms1_corr_v_gauss: Vec::with_capacity(capacity),
            ms2_cosine_ref_sim: Vec::with_capacity(capacity),
            ms2_coelution_score: Vec::with_capacity(capacity),
            ms2_lazyscore: Vec::with_capacity(capacity),
            ms2_corr_v_gauss: Vec::with_capacity(capacity),
            main_score: Vec::with_capacity(capacity),
        }
    }

    pub fn clear(&mut self) {
        self.ms1_cosine_ref_sim.clear();
        self.ms1_coelution_score.clear();
        self.ms1_corr_v_gauss.clear();
        self.ms2_cosine_ref_sim.clear();
        self.ms2_coelution_score.clear();
        self.ms2_lazyscore.clear();
        self.ms2_corr_v_gauss.clear();
        self.main_score.clear();
    }

    /// Resize all buffers to the specified length (filling with 0.0).
    pub fn resize(&mut self, len: usize) {
        self.ms1_cosine_ref_sim.resize(len, 0.0);
        self.ms1_coelution_score.resize(len, 0.0);
        self.ms1_corr_v_gauss.resize(len, 0.0);
        self.ms2_cosine_ref_sim.resize(len, 0.0);
        self.ms2_coelution_score.resize(len, 0.0);
        self.ms2_lazyscore.resize(len, 0.0);
        self.ms2_corr_v_gauss.resize(len, 0.0);
        // main_score is computed later, so we just reserve/clear usually,
        // but for safety we can resize it too.
        self.main_score.resize(len, 0.0);
    }

    pub fn iter_scores(&self) -> impl Iterator<Item = (&'static str, &[f32])> + '_ {
        vec![
            ("ms1_cosine_ref_sim", &self.ms1_cosine_ref_sim[..]),
            ("ms1_coelution_score", &self.ms1_coelution_score[..]),
            ("ms1_corr_v_gauss", &self.ms1_corr_v_gauss[..]),
            ("ms2_cosine_ref_sim", &self.ms2_cosine_ref_sim[..]),
            ("ms2_coelution_score", &self.ms2_coelution_score[..]),
            ("ms2_lazyscore", &self.ms2_lazyscore[..]),
            ("ms2_corr_v_gauss", &self.ms2_corr_v_gauss[..]),
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
    temp_ms2_dot_prod: Vec<f32>,
    temp_ms2_norm_sq_obs: Vec<f32>,
    temp_ms1_dot_prod: Vec<f32>,
    temp_ms1_norm_sq_obs: Vec<f32>,
}

impl ApexFinderBuffers {
    fn new(size: usize) -> Self {
        Self {
            temp_ms2_dot_prod: vec![0.0f32; size],
            temp_ms2_norm_sq_obs: vec![0.0f32; size],
            temp_ms1_dot_prod: vec![0.0f32; size],
            temp_ms1_norm_sq_obs: vec![0.0f32; size],
        }
    }

    fn clear(&mut self) {
        // I can maybe cut some corners by not zeroing the whole vec,
        // but just resizing later.
        // Since every value is over-written anyway.
        // ... This feels safer though.
        self.temp_ms2_dot_prod.fill(0.0);
        self.temp_ms2_norm_sq_obs.fill(0.0);
        self.temp_ms1_dot_prod.fill(0.0);
        self.temp_ms1_norm_sq_obs.fill(0.0);
    }

    fn resize(&mut self, len: usize) {
        self.temp_ms2_dot_prod.resize(len, 0.0);
        self.temp_ms2_norm_sq_obs.resize(len, 0.0);
        self.temp_ms1_dot_prod.resize(len, 0.0);
        self.temp_ms1_norm_sq_obs.resize(len, 0.0);
    }
}

impl ApexFinder {
    pub fn new(capacity: usize) -> Self {
        Self {
            traces: ScoreTraces::new_with_capacity(capacity),
            buffers: ApexFinderBuffers::new(capacity),
        }
    }

    /// Find the peptide apex within the provided scoring context.
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

        // 2. Compute scores (Two-Pass approach)
        self.compute_pass_1(scoring_ctx)?;
        self.compute_pass_2(scoring_ctx)?;

        // 3. Smooth scores
        self.smooth_scores();

        // 4. Compute Main Score (Composite)
        self.compute_main_score_trace();

        // 5. Find Apex and Extract Features
        self.extract_apex_score(scoring_ctx, &rt_mapper)
    }

    /// Pass 1: Scores that depend only on individual ion traces.
    /// - Lazyscore (Hyperscore approximation)
    /// - Cosine Similarity vs Expected Intensities
    /// - Gaussian Shape Correlation
    ///
    /// Can be an error if no valid ions are found for scoring.
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
        // We use accumulators to compute cosine similarity and lazyscore in one go.
        // Lazyscore ~ Sum(ln(1 + intensity))
        // Cosine ~ DotProduct(obs, exp) / (Norm(obs) * Norm(exp))

        let ms2_dot_prod = &mut self.buffers.temp_ms2_dot_prod;
        let ms2_norm_sq_obs = &mut self.buffers.temp_ms2_norm_sq_obs;
        let mut ms2_norm_sq_exp = 0.0f32; // Scalar, since expected is a single vector

        for ((key, _mz), chrom) in collector.fragments.iter_mzs() {
            let expected = scoring_ctx
                .expected_intensities
                .fragment_intensities
                .get(key)
                .copied()
                .unwrap_or(0.0);

            if expected <= 0.0 {
                continue;
            }
            ms2_norm_sq_exp += expected * expected;

            for (i, &intensity) in chrom.iter().enumerate() {
                if intensity > 0.0 {
                    // Lazyscore: lnfact of sum of logs... simplified here to sum of logs
                    // Actual implementation in `hyperscore.rs` uses `lnfact_f32(sum(ln(x)))`.
                    // We'll accumulate the log sums here.
                    let ln_val = intensity.max(1.0).ln();
                    self.traces.ms2_lazyscore[i] += ln_val;

                    // Cosine parts
                    ms2_dot_prod[i] += intensity * expected;
                    ms2_norm_sq_obs[i] += intensity * intensity;
                }
            }
        }

        // Finalize MS2 Cosine & Lazyscore
        let ms2_norm_exp = ms2_norm_sq_exp.sqrt();
        for i in 0..self.traces.ms2_cosine_ref_sim.len() {
            // Finalize Lazyscore
            self.traces.ms2_lazyscore[i] =
                crate::utils::math::lnfact_f32(self.traces.ms2_lazyscore[i]);

            // Finalize Cosine
            let obs_norm = ms2_norm_sq_obs[i].sqrt();
            if obs_norm > 0.0 && ms2_norm_exp > 0.0 {
                self.traces.ms2_cosine_ref_sim[i] = ms2_dot_prod[i] / (obs_norm * ms2_norm_exp);
                // Clip to valid range and min value
                self.traces.ms2_cosine_ref_sim[i] =
                    self.traces.ms2_cosine_ref_sim[i].max(1e-3).min(1.0);
            } else {
                self.traces.ms2_cosine_ref_sim[i] = 1e-3;
            }
        }

        // --- MS2 Gaussian Correlation ---
        // We can reuse the existing function from `corr_v_ref` as it takes `MzMajorIntensityArray`.
        // It iterates ions internally.
        corr_v_ref::calculate_cosine_with_ref_gaussian_into(
            &collector.fragments,
            |_| true, // Filter: take all
            &mut self.traces.ms2_corr_v_gauss,
        )?;

        // --- MS1 (Precursors) ---
        // Similar logic for MS1
        let ms1_dot_prod = &mut self.buffers.temp_ms1_dot_prod;
        let ms1_norm_sq_obs = &mut self.buffers.temp_ms1_norm_sq_obs;
        let mut ms1_norm_sq_exp = 0.0f32;

        for ((key, _mz), chrom) in collector.precursors.iter_mzs() {
            let expected = scoring_ctx
                .expected_intensities
                .precursor_intensities
                .get(key)
                .copied()
                .unwrap_or(0.0);
            if expected <= 0.0 {
                continue;
            }
            ms1_norm_sq_exp += expected * expected;

            for (i, &intensity) in chrom.iter().enumerate() {
                if intensity > 0.0 {
                    ms1_dot_prod[i] += intensity * expected;
                    ms1_norm_sq_obs[i] += intensity * intensity;
                }
            }
        }

        let ms1_norm_exp = ms1_norm_sq_exp.sqrt();
        for i in 0..self.traces.ms1_cosine_ref_sim.len() {
            let obs_norm = ms1_norm_sq_obs[i].sqrt();
            if obs_norm > 0.0 && ms1_norm_exp > 0.0 {
                self.traces.ms1_cosine_ref_sim[i] = ms1_dot_prod[i] / (obs_norm * ms1_norm_exp);
                self.traces.ms1_cosine_ref_sim[i] =
                    self.traces.ms1_cosine_ref_sim[i].max(1e-3).min(1.0);
            } else {
                self.traces.ms1_cosine_ref_sim[i] = 1e-3;
            }
        }

        corr_v_ref::calculate_cosine_with_ref_gaussian_into(
            &collector.precursors,
            |&k| k >= 0, // Filter: Ignore negative keys (decoys/invalid)
            &mut self.traces.ms1_corr_v_gauss,
        )?;

        Ok(())
    }

    /// Pass 2: Scores that depend on the results of Pass 1.
    /// - Coelution Score (compares individual ion traces vs the aggregated Lazyscore trace)
    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn compute_pass_2<T: KeyLike>(
        &mut self,
        scoring_ctx: &ScoringContext<T>,
    ) -> Result<(), DataProcessingError> {
        let collector = &scoring_ctx.query_values;

        // Apply smoothing to Lazyscore BEFORE using it as a reference for Coelution
        // Note: The original code smoothed EVERYTHING at the end.
        // But Coelution uses `ms2_lazyscore` as the reference shape.
        // Does it use the smoothed or raw?
        // Checking `calculate_scores.rs`: `smooth_scores()` happens AFTER `calculate_coelution_scores`.
        // So it uses the RAW lazyscore. Okay, we proceed with raw.

        // MS2 Coelution
        coelution_vref_score_filter_into(
            &collector.fragments,
            &self.traces.ms2_lazyscore,
            COELUTION_WINDOW_WIDTH,
            &|_| true,
            &mut self.traces.ms2_coelution_score,
        )?;

        // MS1 Coelution (vs MS2 Lazyscore reference)
        coelution_vref_score_filter_into(
            &collector.precursors,
            &self.traces.ms2_lazyscore,
            COELUTION_WINDOW_WIDTH,
            &|x: &i8| *x >= 0i8,
            &mut self.traces.ms1_coelution_score,
        )?;

        Ok(())
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn smooth_scores(&mut self) {
        gaussblur_in_place(&mut self.traces.ms2_lazyscore);
        gaussblur_in_place(&mut self.traces.ms1_coelution_score);
        gaussblur_in_place(&mut self.traces.ms2_coelution_score);
        gaussblur_in_place(&mut self.traces.ms2_cosine_ref_sim);
        gaussblur_in_place(&mut self.traces.ms1_cosine_ref_sim);
        gaussblur_in_place(&mut self.traces.ms2_corr_v_gauss);
        gaussblur_in_place(&mut self.traces.ms1_corr_v_gauss);
    }

    #[cfg_attr(
        feature = "instrumentation",
        tracing::instrument(skip_all, level = "trace")
    )]
    fn compute_main_score_trace(&mut self) {
        let len = self.traces.ms1_corr_v_gauss.len();
        self.traces.main_score.clear();
        self.traces.main_score.reserve(len);

        const MS1_SCALING: f32 = 0.75;
        const MS1_OFFSET: f32 = 0.25;

        for i in 0..len {
            let ms1_cos_score =
                MS1_SCALING + (MS1_OFFSET * self.traces.ms1_coelution_score[i].max(1e-3).powi(2));
            let ms1_gauss_score =
                MS1_SCALING + (MS1_OFFSET * self.traces.ms1_corr_v_gauss[i].max(1e-3).powi(2));

            let mut loc_score = 1.0;
            loc_score *= ms1_cos_score;
            loc_score *= ms1_gauss_score;
            loc_score *= self.traces.ms2_cosine_ref_sim[i].max(1e-3).powi(2);
            loc_score *= self.traces.ms2_coelution_score[i].max(1e-3).powi(2);
            loc_score *= self.traces.ms2_corr_v_gauss[i].max(1e-3).powi(2);

            self.traces.main_score.push(loc_score);
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

        // Calculate Peak Shape (Raising/Falling) to determine width
        let (raising_cycles, falling_cycles) = self.calculate_rise_and_fall_cycles(max_loc);

        // Mask the current peak
        peak_picker.mask_peak(max_loc, raising_cycles as usize, falling_cycles as usize, 2);

        // Find next peaks
        let (next_val, next_loc) = peak_picker.next_peak().unwrap_or((0.0, max_loc));

        // For second next, we need to mask the 'next' peak properly.
        let (next_raise, next_fall) = self.calculate_rise_and_fall_cycles(next_loc);
        peak_picker.mask_peak(next_loc, next_raise as usize, next_fall as usize, 1);

        let (second_next_val, _) = peak_picker.next_peak().unwrap_or((0.0, max_loc));

        let delta_next = max_val - next_val;
        let delta_second_next = max_val - second_next_val;

        // Extract features at max_loc
        let cycle_offset = scoring_ctx.query_values.cycle_offset();
        let global_loc = max_loc + cycle_offset;
        let retention_time_ms = rt_mapper(global_loc);

        let (ms1_summed_intensity, _ms1_npeaks) =
            self.sum_intensities_at(&scoring_ctx.query_values.precursors, max_loc);
        let (ms2_summed_intensity, ms2_npeaks) =
            self.sum_intensities_at(&scoring_ctx.query_values.fragments, max_loc);

        // Calculate Lambda (Baseline noise level)
        let lambda = self.calculate_baseline_lambda(max_loc, raising_cycles, falling_cycles);
        let k = self.traces.ms2_lazyscore[max_loc] as f64;
        let norm_lazy_std = lambda.sqrt().max(1.0) as f32;
        let lazyscore_z = self.traces.ms2_lazyscore[max_loc] / norm_lazy_std;

        if lazyscore_z.is_nan() {
            return Err(DataProcessingError::ExpectedFiniteNonNanData {
                context: format!("Lazy score is NaN {} and {}", k, norm_lazy_std),
            });
        }

        Ok(ApexScore {
            score: max_val,
            delta_next,
            delta_second_next,
            retention_time_ms,
            ms2_cosine_ref_sim: self.traces.ms2_cosine_ref_sim[max_loc],
            ms2_coelution_score: self.traces.ms2_coelution_score[max_loc],
            ms2_summed_intensity,
            npeaks: ms2_npeaks as u8,
            lazyscore: self.traces.ms2_lazyscore[max_loc],
            lazyscore_vs_baseline: (k / lambda) as f32,
            lazyscore_z,
            ms2_corr_v_gauss: self.traces.ms2_corr_v_gauss[max_loc],
            ms1_corr_v_gauss: self.traces.ms1_corr_v_gauss[max_loc],
            ms1_cosine_ref_sim: self.traces.ms1_cosine_ref_sim[max_loc],
            ms1_coelution_score: self.traces.ms1_coelution_score[max_loc],
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
