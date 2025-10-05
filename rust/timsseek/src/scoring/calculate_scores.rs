//! # Peptide Scoring Module
//!
//! This module is responsible for scoring peptide candidates against observed mass spectrometry data.
//! It takes theoretical peptide information and raw chromatographic data as input, and it produces
//! a `MainScore` that quantifies the quality of the peptide match.
//!
//! ## Data Flow
//!
//! The scoring process can be broken down into the following stages:
//!
//! 1.  **Input (`PreScore`)**: The process starts with a `PreScore` object, which bundles a peptide
//!     candidate (`DigestSlice`) with its expected theoretical isotope and fragment intensities
//!     (`ExpectedIntensities`), and the raw experimental data collected from the instrument
//!     (`ChromatogramCollector`).
//!
//! 2.  **Data Preparation (`IntensityArrays`)**: The raw data from the `ChromatogramCollector` is
//!     transformed into `IntensityArrays`. These arrays are matrix-like structures that organize
//!     the intensity data by retention time, making it efficient to perform calculations across
//!     the time series.
//!
//! 3.  **Time-Series Feature Calculation (`LongitudinalMainScoreElements`)**: This is the most
//!     complex step. A large struct, `LongitudinalMainScoreElements`, is used as a buffer to
//!     calculate a wide variety of scores for *each point* in the chromatogram's time series.
//!     These "longitudinal" scores include:
//!     - Cosine similarity between observed and theoretical spectra.
//!     - Co-elution scores, measuring how well different fragment ions appear and disappear together.
//!     - Correlation with a Gaussian peak shape.
//!     - "Lazyscore," a measure of local signal-to-noise.
//!
//!     After calculation, these time-series scores are smoothed using a Gaussian blur to reduce noise.
//!
//! 4.  **Apex Detection**: A composite `main_score` is calculated for each time point by multiplying
//!     the various longitudinal scores together. The algorithm then identifies the time points with the
//!     highest `main_score` values, which are considered potential "apexes" of the eluting peptide.
//!
//! 5.  **Final Score Generation (`MainScore`)**: At the best apex found, the system extracts a final
//!     set of features. This includes the value of the main score at the apex, the difference in score
//!     to the next-best apex (`delta_next`), the width of the peak, and various other sub-scores.
//!     These are all collected into the final `MainScore` struct, which is the ultimate output of the
//!     scoring process.
//!
//! ## Planned Refactoring
//!
//! The current implementation has several design issues, including high coupling between structs and
//! an unclear API. A refactoring is planned to introduce a `PeptideScorer` struct that will
//! encapsulate the scoring logic and manage internal buffers, providing a cleaner and more ergonomic
//! interface.

use super::scores::coelution::coelution_score;
use super::scores::corr_v_ref::calculate_cosine_with_ref_gaussian_into;
use super::scores::{
    corr_v_ref,
    hyperscore,
};
use super::{
    COELUTION_WINDOW_WIDTH,
    NUM_MS1_IONS,
    NUM_MS2_IONS,
};
use crate::errors::DataProcessingError;
use crate::models::DigestSlice;
use crate::utils::rolling_calculators::{
    calculate_centered_std,
    calculate_value_vs_baseline_into,
};
use crate::utils::top_n_array::TopNArray;
use crate::{
    ExpectedIntensities,
    IonAnnot,
};
use core::f32;
use serde::Serialize;
use std::sync::Arc;
use timsquery::models::aggregators::ChromatogramCollector;
use timsquery::models::{
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
use timsquery::{
    MzMobilityStatsCollector,
    SpectralCollector,
};
use tracing::warn;

/// Represents a peptide candidate to be scored.
///
/// This struct acts as a data container, bundling the theoretical information about a peptide
/// with the corresponding raw experimental data necessary for scoring. It is the primary input
/// to the [`PeptideScorer::score`] method.
#[derive(Debug)]
pub struct PreScore {
    /// The peptide sequence and modification information.
    pub digest: DigestSlice,
    /// The precursor charge state.
    pub charge: u8,
    /// The expected theoretical intensities of precursor and fragment ions.
    pub expected_intensities: ExpectedIntensities,
    /// The observed chromatogram data collected from the instrument.
    pub query_values: ChromatogramCollector<IonAnnot, f32>,
}

#[derive(Debug, Serialize, Clone)]
pub struct TimeResolvedScores {
    pub ms1_cosine_ref_sim: Vec<f32>,
    pub ms1_coelution_score: Vec<f32>,
    pub ms1_corr_v_gauss: Vec<f32>,
    pub ms2_cosine_ref_sim: Vec<f32>,
    pub ms2_coelution_score: Vec<f32>,
    pub ms2_lazyscore: Vec<f32>,
    pub ms2_lazyscore_vs_baseline: Vec<f32>,
    pub ms2_corr_v_gauss: Vec<f32>,
    pub ms2_lazyscore_vs_baseline_std: f32,
}

#[derive(Debug)]
pub struct IntensityArrays {
    // TODO: Reimplement to use directly the chromatogram collector
    pub ms1_rtmajor: RTMajorIntensityArray<i8, f32>,
    pub ms1_mzmajor: MzMajorIntensityArray<i8, f32>,
    pub ms2_rtmajor: RTMajorIntensityArray<IonAnnot, f32>,
    pub ms2_mzmajor: MzMajorIntensityArray<IonAnnot, f32>,
    pub ms1_expected_intensities: Vec<f32>,
    pub ms2_expected_intensities: Vec<f32>,
}

impl IntensityArrays {
    pub fn new_empty(
        num_ms1: usize,
        num_ms2: usize,
        ref_time_ms: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        // This is mainly used to pre-allocated the memory that will be used later by several.
        // runs.
        let ms1_order: Arc<[(i8, f64)]> = (0..=num_ms1).map(|_o| (-1i8, 867.8309)).collect();
        let ms2_order: Arc<[(IonAnnot, f64)]> = (0..=num_ms2)
            .map(|_o| (IonAnnot::new('p', None, 1, 0).unwrap(), 867.8309))
            .collect();
        let ms2_ref_vec: Vec<f32> = (0..=num_ms2).map(|_| 0.0).collect();
        Ok(Self {
            ms1_rtmajor: RTMajorIntensityArray::try_new_empty(
                ms1_order.clone(),
                ref_time_ms.clone(),
            )?,
            ms1_mzmajor: MzMajorIntensityArray::try_new_empty(
                ms1_order.clone(),
                ref_time_ms.clone(),
            )?,
            ms2_rtmajor: RTMajorIntensityArray::try_new_empty(
                ms2_order.clone(),
                ref_time_ms.clone(),
            )?,
            ms2_mzmajor: MzMajorIntensityArray::try_new_empty(
                ms2_order.clone(),
                ref_time_ms.clone(),
            )?,
            ms1_expected_intensities: vec![0.5; num_ms1],
            ms2_expected_intensities: ms2_ref_vec,
        })
    }

    // Right now I KNOW this is an icredibly dumb operation
    // and I should move this upstream BUT I want to finish the other
    // refactor changes first ... JSPP Apr-18-2025
    pub fn reset_with(
        &mut self,
        intensity_arrays: &ChromatogramCollector<IonAnnot, f32>,
        expected_intensities: &ExpectedIntensities,
    ) -> Result<(), DataProcessingError> {
        self.ms1_rtmajor
            .try_reset_with(&intensity_arrays.precursors)?;
        self.ms1_mzmajor
            .try_reset_with(&intensity_arrays.precursors)?;
        self.ms2_rtmajor
            .try_reset_with(&intensity_arrays.fragments)?;
        self.ms2_mzmajor
            .try_reset_with(&intensity_arrays.fragments)?;
        self.ms1_expected_intensities = expected_intensities.precursor_intensities.clone();
        let ms2_ref_vec: Vec<_> = self
            .ms2_mzmajor
            .mz_order
            .iter()
            .map(|&(k, _)| {
                *expected_intensities
                    .fragment_intensities
                    .get(&k)
                    .expect("Failed to find expected intensity for fragment")
            })
            .collect();
        self.ms2_expected_intensities = ms2_ref_vec;
        Ok(())
    }
}

/// Applies a 1D Gaussian blur in-place with a fixed kernel [0.5, 1.0, 0.5].
///
/// This implementation is memory-efficient, using only a single temporary
/// variable instead of a full-sized buffer.
pub fn gaussblur_in_place(x: &mut [f32]) {
    let len = x.len();
    if len < 3 {
        return;
    }

    // Kernel weights and normalization factor
    const W_SIDE: f32 = 0.5;
    const W_CENTER: f32 = 1.0;
    const NORM: f32 = W_SIDE + W_CENTER + W_SIDE; // This is 2.0

    // Store the original value of the first element before it's modified.
    // This will act as our "previous" value for the first iteration of the loop.
    let mut prev_val = x[0];

    // Handle the first element (index 0) using your custom edge logic.
    // Note: This calculation uses the original x[0] and x[1].
    x[0] = (x[0] * 1.5 + x[1] * 0.5) / NORM;

    // Main convolution loop for the elements in the middle.
    // We iterate from the second element up to the second-to-last.
    for i in 1..len - 1 {
        // `current_val` holds the original value of x[i] before we overwrite it.
        let current_val = x[i];

        // Calculate the new value for x[i].
        // - `prev_val` is the original `x[i-1]`
        // - `current_val` is the original `x[i]`
        // - `x[i+1]` is the original `x[i+1]` (it hasn't been touched yet)
        x[i] = (prev_val * W_SIDE + current_val * W_CENTER + x[i + 1] * W_SIDE) / NORM;

        // For the *next* iteration, the current value will be the previous one.
        prev_val = current_val;
    }

    // Handle the last element (index len - 1).
    // At this point, `prev_val` holds the original value of `x[len-2]`.
    // The value `x[len-1]` is still in its original state.
    x[len - 1] = (x[len - 1] * 1.5 + prev_val * 0.5) / NORM;
}

impl TimeResolvedScores {
    pub fn try_reset_with(
        &mut self,
        intensity_arrays: &IntensityArrays,
    ) -> Result<(), DataProcessingError> {
        self.clear();
        self.calculate_lazyscores(intensity_arrays)?;
        self.calculate_similarity_scores(intensity_arrays);
        self.calculate_coelution_scores(intensity_arrays)?;
        self.calculate_gaussian_correlation_scores(intensity_arrays)?;
        self.smooth_scores();
        self.calculate_baseline_scores(intensity_arrays);
        self.check_lengths()?;
        Ok(())
    }

    fn calculate_lazyscores(
        &mut self,
        intensity_arrays: &IntensityArrays,
    ) -> Result<(), DataProcessingError> {
        self.ms2_lazyscore
            .extend(hyperscore::lazyscore(&intensity_arrays.ms2_rtmajor));

        let max_lzs = self
            .ms2_lazyscore
            .iter()
            .max_by(|x, y| {
                if x.is_nan() {
                    std::cmp::Ordering::Less
                } else if y.is_nan() {
                    std::cmp::Ordering::Greater
                } else {
                    x.partial_cmp(y).unwrap()
                }
            })
            .unwrap_or(&0.0);

        if *max_lzs <= 1e-3 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("No non-0 lazyscore".into()),
            });
        }
        Ok(())
    }

    fn calculate_similarity_scores(&mut self, intensity_arrays: &IntensityArrays) {
        self.ms1_cosine_ref_sim = corr_v_ref::calculate_cosine_with_ref(
            &intensity_arrays.ms1_rtmajor,
            &intensity_arrays.ms1_expected_intensities,
        )
        .unwrap();
        self.ms2_cosine_ref_sim = corr_v_ref::calculate_cosine_with_ref(
            &intensity_arrays.ms2_rtmajor,
            &intensity_arrays.ms2_expected_intensities,
        )
        .unwrap();

        self.ms1_cosine_ref_sim
            .iter_mut()
            .for_each(|x| *x = x.max(1e-3));
        self.ms2_cosine_ref_sim
            .iter_mut()
            .for_each(|x| *x = x.max(1e-3));
    }

    fn calculate_coelution_scores(
        &mut self,
        intensity_arrays: &IntensityArrays,
    ) -> Result<(), DataProcessingError> {
        let filter: Option<fn(&IonAnnot) -> bool> = None;
        self.ms2_coelution_score
            .extend(
                coelution_score::coelution_vref_score_filter::<10, IonAnnot>(
                    &intensity_arrays.ms2_mzmajor,
                    self.ms2_lazyscore.as_slice(),
                    COELUTION_WINDOW_WIDTH,
                    &filter,
                )?,
            );
        self.ms1_coelution_score
            .extend(coelution_score::coelution_vref_score_filter::<6, i8>(
                &intensity_arrays.ms1_mzmajor,
                // Note we DO use the ms2 lazyscore as reference here
                self.ms2_lazyscore.as_slice(),
                7,
                &Some(|x: &i8| *x >= 0i8),
            )?);
        Ok(())
    }

    fn calculate_gaussian_correlation_scores(
        &mut self,
        intensity_arrays: &IntensityArrays,
    ) -> Result<(), DataProcessingError> {
        calculate_cosine_with_ref_gaussian_into(
            &intensity_arrays.ms2_mzmajor,
            |_| true,
            &mut self.ms2_corr_v_gauss,
        )?;
        calculate_cosine_with_ref_gaussian_into(
            &intensity_arrays.ms1_mzmajor,
            |&k| k >= 0,
            &mut self.ms1_corr_v_gauss,
        )?;
        Ok(())
    }

    fn smooth_scores(&mut self) {
        gaussblur_in_place(&mut self.ms2_lazyscore);
        gaussblur_in_place(&mut self.ms1_coelution_score);
        gaussblur_in_place(&mut self.ms2_coelution_score);
        gaussblur_in_place(&mut self.ms2_cosine_ref_sim);
        gaussblur_in_place(&mut self.ms1_cosine_ref_sim);
        gaussblur_in_place(&mut self.ms2_corr_v_gauss);
        gaussblur_in_place(&mut self.ms1_corr_v_gauss);
    }

    fn calculate_baseline_scores(&mut self, intensity_arrays: &IntensityArrays) {
        let rt_len = intensity_arrays.ms1_rtmajor.rts_ms.len();
        let five_pct_index = rt_len * 5 / 100;
        let half_five_pct_index = five_pct_index / 2;

        if five_pct_index > 100 {
            warn!(
                "High five_pct_index: {} (from rt_len: {})",
                five_pct_index, rt_len
            );
        }

        calculate_value_vs_baseline_into(
            &self.ms2_lazyscore,
            five_pct_index,
            &mut self.ms2_lazyscore_vs_baseline,
        );

        let baseline_slice = &self.ms2_lazyscore_vs_baseline
            [half_five_pct_index..(self.ms2_lazyscore_vs_baseline.len() - half_five_pct_index)];
        let lzb_std = calculate_centered_std(baseline_slice);

        self.ms2_lazyscore_vs_baseline_std = lzb_std.max(1.0);
    }

    fn check_lengths(&self) -> Result<(), DataProcessingError> {
        let len = self.ms1_cosine_ref_sim.len();
        macro_rules! check_len {
            ($x:expr) => {
                if $x.len() != len {
                    return Err(DataProcessingError::ExpectedSlicesSameLength {
                        expected: len,
                        other: $x.len(),
                        context: format!("Slices not same length -> {:?}", stringify!($x)).into(),
                    });
                }
            };
        }

        // This is reallty not necessary but exhaustive pattern
        // matching helps me not forget any of the fields.
        let Self {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore,
            ms2_lazyscore_vs_baseline,
            // split_lazyscore,
            ms2_corr_v_gauss,
            ms1_corr_v_gauss,
            ms2_lazyscore_vs_baseline_std: _,
        } = self;

        check_len!(ms1_cosine_ref_sim);
        check_len!(ms1_coelution_score);
        check_len!(ms2_cosine_ref_sim);
        check_len!(ms2_coelution_score);
        check_len!(ms2_lazyscore);
        check_len!(ms2_lazyscore_vs_baseline);
        // check_len!(split_lazyscore);
        check_len!(ms2_corr_v_gauss);
        check_len!(ms1_corr_v_gauss);

        Ok(())
    }

    fn clear(&mut self) {
        let TimeResolvedScores {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore,
            ms2_lazyscore_vs_baseline,
            // split_lazyscore,
            ms2_corr_v_gauss,
            ms1_corr_v_gauss,
            ms2_lazyscore_vs_baseline_std: _,
        } = self;

        ms1_cosine_ref_sim.clear();
        ms1_coelution_score.clear();
        ms2_cosine_ref_sim.clear();
        ms2_coelution_score.clear();
        ms2_lazyscore.clear();
        ms2_lazyscore_vs_baseline.clear();
        // split_lazyscore.clear();
        ms2_corr_v_gauss.clear();
        ms1_corr_v_gauss.clear();
    }

    pub fn new_with_capacity(capacity: usize) -> Self {
        let ms1_cosine_ref_sim = Vec::with_capacity(capacity);
        let ms1_coelution_score = Vec::with_capacity(capacity);
        let ms2_cosine_ref_sim = Vec::with_capacity(capacity);
        let ms2_coelution_score = Vec::with_capacity(capacity);
        let ms2_lazyscore = Vec::with_capacity(capacity);
        let ms2_lazyscore_vs_baseline = Vec::with_capacity(capacity);
        // let split_lazyscore = Vec::with_capacity(capacity);
        let ms2_corr_v_gauss = Vec::with_capacity(capacity);
        let ms1_corr_v_gauss = Vec::with_capacity(capacity);
        Self {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore,
            ms2_lazyscore_vs_baseline,
            // split_lazyscore,
            ms2_lazyscore_vs_baseline_std: 1.0,
            ms2_corr_v_gauss,
            ms1_corr_v_gauss,
        }
    }

    fn find_apex_candidates(&self) -> [ScoreInTime; 20] {
        let mut candidate_groups: [TopNArray<2, ScoreInTime>; 21] = [TopNArray::new(); 21];
        let five_pct_index = self.ms1_corr_v_gauss.len() * 5 / 100;
        for (i, score) in self.main_score_iter().enumerate() {
            if score.is_nan() {
                continue;
            }
            let sit = ScoreInTime { score, index: i };
            candidate_groups[i / five_pct_index].push(sit);
        }
        let mut candidates: TopNArray<20, ScoreInTime> = TopNArray::new();
        for c in candidate_groups.iter() {
            for val in c.get_values().iter() {
                candidates.push(*val);
            }
        }
        candidates.get_values()
    }

    pub fn main_score_iter(&self) -> impl '_ + Iterator<Item = f32> {
        (0..self.ms1_corr_v_gauss.len()).map(|x| self.main_score_at(x))
    }

    fn main_score_at(&self, idx: usize) -> f32 {
        // The 0.75 - 0.25 means we are downscaling the main lazyscore up to 0.75
        // since the similarity is in the 0-1 range; even if the precursor has
        // similarity of 0, we still have a scoring value.
        const MS1_SCALING: f32 = 0.25;
        const MS1_OFFSET: f32 = 0.75;

        let ms1_cos_score =
            MS1_SCALING + (MS1_OFFSET * self.ms1_cosine_ref_sim[idx].max(1e-3).powi(2));
        let ms1_gauss_score =
            MS1_SCALING + (MS1_OFFSET * self.ms1_corr_v_gauss[idx].max(1e-3).powi(2));
        let mut loc_score = self.ms2_lazyscore[idx];
        loc_score *= ms1_cos_score;
        loc_score *= ms1_gauss_score;
        loc_score *= self.ms2_cosine_ref_sim[idx].max(1e-3).powi(2);
        loc_score *= self.ms2_coelution_score[idx].max(1e-3).powi(2);
        loc_score *= self.ms2_corr_v_gauss[idx].max(1e-3).powi(2);
        loc_score
    }
}

// TODO: Move this to somewhere that makes sense ...
/// If a negative step is passed, it essentially means we are
/// going left.
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

#[derive(Debug, Clone, Copy, PartialEq)]
struct ScoreInTime {
    score: f32,
    index: usize,
}

impl Default for ScoreInTime {
    fn default() -> Self {
        Self {
            score: f32::NAN,
            index: 0,
        }
    }
}

/// The primary entry point for scoring peptide candidates.
///
/// `PeptideScorer` encapsulates the complex scoring logic and manages reusable memory buffers
/// for efficiency. It is designed to be created once and then used to score multiple peptide
/// candidates (`PreScore` objects) against the raw data from a single mass spectrometry run.
///
/// # Example
///
/// ```ignore
/// // 1. Create a scorer instance.
/// let mut scorer = PeptideScorer::new()?;
///
/// // 2. Create a PreScore object with peptide info and raw data.
/// let prescore = PreScore { ... };
///
/// // 3. Score the peptide.
/// let main_score = scorer.score(&prescore)?;
///
/// // 4. Use the results.
/// println!("Score for peptide: {}", main_score.score);
/// ```
pub struct PeptideScorer {
    intensity_arrays: IntensityArrays,
    time_resolved_scores: TimeResolvedScores,
}

impl PeptideScorer {
    /// Creates a new `PeptideScorer`.
    ///
    /// This initializes the scorer and pre-allocates internal buffers. For performance,
    /// you should create one `PeptideScorer` per analysis run and reuse it for all peptides.
    pub fn new(
        n_ms1: usize,
        n_ms2: usize,
        rt_ms_arc: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        // TODO: The initialization could be more robust, perhaps by taking the maximum
        // number of time points expected as an argument to pre-allocate buffers more accurately.
        let time_resolved_scores = TimeResolvedScores::new_with_capacity(rt_ms_arc.len());
        let intensity_arrays = IntensityArrays::new_empty(n_ms1, n_ms2, rt_ms_arc)?;
        Ok(Self {
            intensity_arrays,
            time_resolved_scores,
        })
    }

    /// Scores a peptide candidate against the experimental data.
    ///
    /// This is the main method of the scorer. It takes a `PreScore` object, which contains
    /// the theoretical peptide information and the observed chromatogram data, and returns a
    /// `MainScore` struct summarizing the quality of the match.
    ///
    /// # Arguments
    ///
    /// * `prescore` - A reference to a `PreScore` object containing the data to be scored.
    ///
    /// # Returns
    ///
    /// A `Result` containing the `MainScore` on success, or a `DataProcessingError` if
    /// the scoring fails (e.g., due to insufficient data).
    pub fn score(&mut self, prescore: &PreScore) -> Result<MainScore, DataProcessingError> {
        self.intensity_arrays
            .reset_with(&prescore.query_values, &prescore.expected_intensities)?;
        let main_score = self.calculate_scores_with_intensities(prescore)?;

        // Final check to ensure the score is valid before returning.
        if main_score.score.is_nan() {
            warn!("Main score is NaN");
            return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
        }
        Ok(main_score)
    }

    /// Returns a reference to the time-resolved scores from the last scoring operation.
    ///
    /// This is useful for debugging or when detailed, time-series data is needed for analysis.
    pub fn last_time_resolved_scores(&self) -> &TimeResolvedScores {
        &self.time_resolved_scores
    }

    fn calculate_scores_with_intensities(
        &mut self,
        prescore: &PreScore,
    ) -> Result<MainScore, DataProcessingError> {
        // I prbably should pass an allocator here ...
        // Since there is so much stuff happening ...
        self.time_resolved_scores
            .try_reset_with(&self.intensity_arrays)?;

        let apex_candidates = self.time_resolved_scores.find_apex_candidates();
        let norm_lazy_std =
            calculate_centered_std(&self.time_resolved_scores.ms2_lazyscore_vs_baseline);
        let max_val = apex_candidates[0].score;
        let max_loc = apex_candidates[0].index;

        if max_val == 0.0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("No non-0 main score".into()),
            });
        }

        // This is a delta next with the constraint that it has to be more than 5% of the max
        // index apart from the max.
        let ten_pct_index = prescore.query_values.fragments.rts_ms.len() / 20;
        let max_window =
            max_loc.saturating_sub(ten_pct_index)..max_loc.saturating_add(ten_pct_index);
        let next = apex_candidates
            .iter()
            .find(|x| !max_window.contains(&x.index));
        let second_next = match next {
            Some(next) => {
                let next_window = next.index.saturating_sub(ten_pct_index)
                    ..next.index.saturating_add(ten_pct_index);
                apex_candidates
                    .iter()
                    .find(|x| !max_window.contains(&x.index) && !next_window.contains(&x.index))
            }
            None => None,
        };

        let delta_next = match next {
            Some(next) => max_val - next.score,
            None => f32::NAN,
        };
        let delta_second_next = match second_next {
            Some(next) => max_val - next.score,
            None => f32::NAN,
        };

        // FOR NOW I will leave this assert and if it holds this is an assumption I can make and
        // I will remove the dead code above.
        assert_eq!(
            prescore.query_values.fragments.rts_ms,
            prescore.query_values.precursors.rts_ms
        );
        let (ms1_loc, ms2_loc) = (max_loc, max_loc);
        let ref_time_ms = prescore.query_values.precursors.rts_ms[max_loc];

        let summed_ms1_int: f32 = match self.intensity_arrays.ms1_rtmajor.arr.get_row(ms1_loc) {
            Some(row) => row.iter().sum(),
            None => 0.0,
        };
        let summed_ms2_int: f32 = match self.intensity_arrays.ms2_rtmajor.arr.get_row(ms2_loc) {
            Some(row) => row.iter().sum(),
            None => 0.0,
        };
        let npeak_ms2 = match self.intensity_arrays.ms2_rtmajor.arr.get_row(ms2_loc) {
            Some(row) => {
                let mut count = 0;
                for x in row {
                    if *x > 1.0f32 {
                        count += 1;
                    }
                }
                count
            }
            None => 0,
        };

        let lazyscore_z = self.time_resolved_scores.ms2_lazyscore[max_loc] / norm_lazy_std;
        if lazyscore_z.is_nan() {
            let tmp = format!(
                "Lazy score is NaN {} and {}",
                self.time_resolved_scores.ms2_lazyscore[max_loc], norm_lazy_std
            );
            return Err(DataProcessingError::ExpectedFiniteNonNanData { context: tmp });
        }
        let raising_cycles = count_falling_steps(
            max_loc,
            -1,
            self.time_resolved_scores.ms2_lazyscore.as_slice(),
        );
        let falling_cycles = count_falling_steps(
            max_loc,
            1,
            self.time_resolved_scores.ms2_lazyscore.as_slice(),
        );

        Ok(MainScore {
            score: max_val,
            delta_next,
            delta_second_next,
            retention_time_ms: ref_time_ms,

            ms2_cosine_ref_sim: self.time_resolved_scores.ms2_cosine_ref_sim[max_loc],
            ms2_coelution_score: self.time_resolved_scores.ms2_coelution_score[max_loc],
            ms2_summed_intensity: summed_ms2_int,
            npeaks: npeak_ms2 as u8,
            lazyscore: self.time_resolved_scores.ms2_lazyscore[max_loc],
            lazyscore_vs_baseline: self.time_resolved_scores.ms2_lazyscore_vs_baseline[max_loc],
            lazyscore_z,
            ms2_corr_v_gauss: self.time_resolved_scores.ms2_corr_v_gauss[max_loc],
            ms1_corr_v_gauss: self.time_resolved_scores.ms1_corr_v_gauss[max_loc],
            // ms1_ms2_correlation: self.time_resolved_scores.ms1_ms2_correlation[max_loc],
            ms1_cosine_ref_sim: self.time_resolved_scores.ms1_cosine_ref_sim[max_loc],
            ms1_coelution_score: self.time_resolved_scores.ms1_coelution_score[max_loc],
            ms1_summed_intensity: summed_ms1_int,

            raising_cycles,
            falling_cycles,
        })
    }
}

impl PartialOrd for ScoreInTime {
    fn partial_cmp(&self, other: &ScoreInTime) -> Option<std::cmp::Ordering> {
        if self.score.is_nan() {
            return Some(std::cmp::Ordering::Less);
        }
        self.score.partial_cmp(&other.score)
    }
}

impl PreScore {}

/// Contains the final scores and features for a peptide candidate at its detected apex.
///
/// This struct is the output of the [`PeptideScorer::score`] method. It provides a comprehensive
/// set of metrics to evaluate the quality of a peptide-spectrum match. The primary metric is
/// `score`, with other fields providing detailed diagnostic information.
#[derive(Debug, Clone, Copy)]
pub struct MainScore {
    /// The main composite score, representing the overall quality of the match. Higher is better.
    pub score: f32,
    /// The difference between the main score and the score of the second-best candidate peak.
    /// A large delta suggests a more confident identification.
    pub delta_next: f32,
    /// The difference between the main score and the score of the third-best candidate peak.
    pub delta_second_next: f32,
    /// The retention time (in milliseconds) at the detected apex of the elution peak.
    pub retention_time_ms: u32,

    // --- MS2-level features ---
    /// Cosine similarity between observed and theoretical MS2 fragment intensities.
    pub ms2_cosine_ref_sim: f32,
    /// Score representing the co-elution of MS2 fragments.
    pub ms2_coelution_score: f32,
    /// Total summed intensity of all MS2 fragment ions at the apex.
    pub ms2_summed_intensity: f32,
    /// The number of MS2 fragment ions detected at the apex.
    pub npeaks: u8,
    /// A signal-to-noise-like score for MS2 fragments.
    pub lazyscore: f32,
    /// The `lazyscore` relative to a local baseline.
    pub lazyscore_vs_baseline: f32,
    /// The Z-score of the `lazyscore`, indicating its statistical significance.
    pub lazyscore_z: f32,
    /// Correlation of the MS2 signal with a Gaussian peak shape.
    pub ms2_corr_v_gauss: f32,

    // --- MS1-level features ---
    /// Correlation of the MS1 signal with a Gaussian peak shape.
    pub ms1_corr_v_gauss: f32,
    /// Cosine similarity between observed and theoretical MS1 isotope intensities.
    pub ms1_cosine_ref_sim: f32,
    /// Score representing the co-elution of MS1 isotopes.
    pub ms1_coelution_score: f32,
    /// Total summed intensity of all MS1 isotope ions at the apex.
    pub ms1_summed_intensity: f32,

    // --- Peak shape features ---
    /// The number of consecutive time points on the leading edge of the peak where the score is increasing.
    pub raising_cycles: u8,
    /// The number of consecutive time points on the trailing edge of the peak where the score is decreasing.
    pub falling_cycles: u8,
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
            if *k >= 0i8 && weight > 0.0 {
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
