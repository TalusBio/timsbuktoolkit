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

/// The PreScore is meant to be in essence a bundle of a query
/// with the collected intensities from the raw data for a specific
/// file.
#[derive(Debug)]
pub struct PreScore {
    pub digest: DigestSlice,
    pub charge: u8,
    pub expected_intensities: ExpectedIntensities,
    pub query_values: ChromatogramCollector<IonAnnot, f32>,
}

#[derive(Debug, Serialize, Clone)]
pub struct LongitudinalMainScoreElements {
    pub ms1_cosine_ref_sim: Vec<f32>,
    pub ms1_coelution_score: Vec<f32>,
    pub ms1_corr_v_gauss: Vec<f32>,
    pub ms2_cosine_ref_sim: Vec<f32>,
    pub ms2_coelution_score: Vec<f32>,
    pub ms2_lazyscore: Vec<f32>,
    pub ms2_lazyscore_vs_baseline: Vec<f32>,
    pub ms2_corr_v_gauss: Vec<f32>,
    pub split_lazyscore: Vec<f32>,

    /// END
    ms2_lazyscore_vs_baseline_std: f32,
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
    pub fn new(
        query_values: &ChromatogramCollector<IonAnnot, f32>,
        expected_intensities: &ExpectedIntensities,
    ) -> Result<Self, DataProcessingError> {
        let (_, ms1_mzmajor_arr, ms2_mzmajor_arr) = query_values.clone().unpack();
        let ms1_rtmajor_arr = ms1_mzmajor_arr.transpose_clone();
        let ms2_rtmajor_arr = ms2_mzmajor_arr.transpose_clone();
        let ms2_ref_vec: Vec<_> = ms2_mzmajor_arr
            .mz_order
            .iter()
            .map(|&(k, _)| {
                *expected_intensities
                    .fragment_intensities
                    .get(&k)
                    .expect("Failed to find expected intensity for fragment")
            })
            .collect();

        Ok(Self {
            ms1_rtmajor: ms1_rtmajor_arr,
            ms1_mzmajor: ms1_mzmajor_arr,
            ms2_rtmajor: ms2_rtmajor_arr,
            ms2_mzmajor: ms2_mzmajor_arr,
            ms1_expected_intensities: expected_intensities.precursor_intensities.clone(),
            ms2_expected_intensities: ms2_ref_vec,
        })
    }

    pub fn new_empty(
        num_ms1: usize,
        num_ms2: usize,
        ref_time_ms: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        // This is mainly used to pre-allocated the memory that will be used later by several.
        // runs.
        let ms1_order: Arc<[(i8, f64)]> = (0..=num_ms1).map(|o| (o as i8, 867.8309)).collect();
        let ms2_order: Arc<[(IonAnnot, f64)]> = (0..=num_ms2)
            .map(|o| {
                (
                    IonAnnot::new('b', Some((o + 1) as u8), 1, 0).unwrap(),
                    867.8309,
                )
            })
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

    pub fn get_top_ions<const N: usize>(&self) {
        todo!()
    }
}

fn gaussblur(x: &mut [f32], buffer: &mut Vec<f32>) {
    // Temp implementation ... shoudl make something nicer in the future.
    let len = x.len();
    if len < 3 {
        return;
    }

    // Using fixed kernel weights [0.5, 1.0, 0.5]
    // Note: These weights are already normalized (sum = 2)
    buffer.clear();
    buffer.resize(len, 0.0);

    // Handle first element
    buffer[0] = (x[0] * 1.5 + x[1] * 0.5) / 2.0;

    // Main convolution loop
    for i in 1..len - 1 {
        buffer[i] = (x[i - 1] * 0.5 + x[i] * 1.0 + x[i + 1] * 0.5) / 2.0;
    }

    // Handle last element
    buffer[len - 1] = (x[len - 1] * 1.5 + x[len - 2] * 0.5) / 2.0;

    // Copy results back
    x.copy_from_slice(buffer);
}

impl LongitudinalMainScoreElements {
    pub fn try_new(intensity_arrays: &IntensityArrays) -> Result<Self, DataProcessingError> {
        let mut out = Self::new_with_capacity(intensity_arrays.ms1_rtmajor.rts_ms.len());
        out.try_reset_with(intensity_arrays)?;
        Ok(out)
    }

    pub fn try_reset_with(
        &mut self,
        intensity_arrays: &IntensityArrays,
    ) -> Result<(), DataProcessingError> {
        self.clear();
        self.ms2_lazyscore
            .extend(hyperscore::lazyscore(&intensity_arrays.ms2_rtmajor));
        let max_lzs = self
            .ms2_lazyscore
            .iter()
            .max_by(|x, y| {
                if x.is_nan() {
                    return std::cmp::Ordering::Less;
                }
                if y.is_nan() {
                    return std::cmp::Ordering::Greater;
                }
                match x.partial_cmp(y) {
                    Some(x) => x,
                    None => panic!("Failed to compare lazyscore values {} vs {}", x, y),
                }
            })
            .unwrap_or(&0.0f32);
        if max_lzs <= &1e-3 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("No non-0 lazyscore".into()),
            });
        }

        // TODO: Reimplement as in-place
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
        // Fill missing
        self.ms1_cosine_ref_sim.iter_mut().for_each(|x| {
            if x.is_nan() {
                *x = 1e-3
            }
        });
        self.ms2_cosine_ref_sim.iter_mut().for_each(|x| {
            if x.is_nan() {
                *x = 1e-3
            }
        });

        let rt_len = intensity_arrays.ms1_rtmajor.rts_ms.len();

        // TODO: reimplement this to make it inplace ...
        // START: Burning hot code ...
        self.ms2_coelution_score = coelution_score::coelution_score::<10, IonAnnot>(
            &intensity_arrays.ms2_mzmajor,
            COELUTION_WINDOW_WIDTH,
        )?;
        let tmp = coelution_score::coelution_score_filter::<6, i8>(
            &intensity_arrays.ms1_mzmajor,
            7,
            &Some(|x: &i8| *x >= 0i8),
        );
        self.ms1_coelution_score = match tmp {
            Ok(scores) => scores,
            Err(_) => vec![0.0; rt_len],
        };
        // END: Burning hot code ...

        self.split_lazyscore.extend(hyperscore::split_ion_lazyscore(
            &intensity_arrays.ms2_rtmajor,
        ));

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

        // This gaussian blur makes it ~ 25% slower ... but it makes the scores
        // significantly "better" (+5% as measured by number of ID's at 1% FDR)
        // and they "look" better.
        let mut buffer = Vec::with_capacity(self.ms2_lazyscore.len());
        gaussblur(&mut self.ms2_lazyscore, &mut buffer);
        gaussblur(&mut self.ms1_coelution_score, &mut buffer);
        gaussblur(&mut self.ms2_coelution_score, &mut buffer);
        gaussblur(&mut self.ms2_cosine_ref_sim, &mut buffer);
        gaussblur(&mut self.ms1_cosine_ref_sim, &mut buffer);
        gaussblur(&mut self.ms2_corr_v_gauss, &mut buffer);
        gaussblur(&mut self.ms1_corr_v_gauss, &mut buffer);
        gaussblur(&mut self.split_lazyscore, &mut buffer);

        let five_pct_index = rt_len * 5 / 100;
        let half_five_pct_idnex = five_pct_index / 2;

        if five_pct_index > 100 {
            eprintln!("{:?}", intensity_arrays.ms1_rtmajor.rts_ms);
            eprintln!(
                "Got a very high five percent index ... {} down from {}",
                five_pct_index, rt_len
            );
            // TODO: handle gradefully cases where
            // long gradients/ultra fast cycle times go down.
            panic!();
        }

        // TODO trim here if the max size is too large
        calculate_value_vs_baseline_into(
            &self.ms2_lazyscore,
            five_pct_index,
            &mut self.ms2_lazyscore_vs_baseline,
        );
        let mut lzb_std = calculate_centered_std(
            &self.ms2_lazyscore_vs_baseline[(half_five_pct_idnex)
                ..(self.ms2_lazyscore_vs_baseline.len() - half_five_pct_idnex)],
        );

        // This feels incredibly dirty ...
        // TODO: Find an elegant way to set this threshold ...
        if lzb_std < 1.0 {
            lzb_std = 1.0;
        }

        self.ms2_lazyscore_vs_baseline_std = lzb_std;

        self.check_lengths()
            .expect("All vecs to be the same length");
        Ok(())
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
            split_lazyscore,
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
        check_len!(split_lazyscore);
        check_len!(ms2_corr_v_gauss);
        check_len!(ms1_corr_v_gauss);

        Ok(())
    }

    fn clear(&mut self) {
        let LongitudinalMainScoreElements {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore,
            ms2_lazyscore_vs_baseline,
            split_lazyscore,
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
        split_lazyscore.clear();
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
        let split_lazyscore = Vec::with_capacity(capacity);
        let ms2_corr_v_gauss = Vec::with_capacity(capacity);
        let ms1_corr_v_gauss = Vec::with_capacity(capacity);
        Self {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore,
            ms2_lazyscore_vs_baseline,
            split_lazyscore,
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

impl PartialOrd for ScoreInTime {
    fn partial_cmp(&self, other: &ScoreInTime) -> Option<std::cmp::Ordering> {
        if self.score.is_nan() {
            return Some(std::cmp::Ordering::Less);
        }
        self.score.partial_cmp(&other.score)
    }
}

impl PreScore {
    fn calc_with_intensities(
        &self,
        intensity_arrays: &IntensityArrays,
        buffer: &mut LongitudinalMainScoreElements,
    ) -> Result<MainScore, DataProcessingError> {
        // I prbably should pass an allocator here ...
        // Since there is so much stuff happening ...
        buffer.try_reset_with(intensity_arrays)?;
        let longitudinal_main_score_elements = buffer;

        let apex_candidates = longitudinal_main_score_elements.find_apex_candidates();
        let norm_lazy_std =
            calculate_centered_std(&longitudinal_main_score_elements.ms2_lazyscore_vs_baseline);
        let max_val = apex_candidates[0].score;
        let max_loc = apex_candidates[0].index;

        if max_val == 0.0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("No non-0 main score".into()),
            });
        }

        // This is a delta next with the constraint that it has to be more than 5% of the max
        // index apart from the max.
        let ten_pct_index = self.query_values.fragments.rts_ms.len() / 20;
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
            self.query_values.fragments.rts_ms,
            self.query_values.precursors.rts_ms
        );
        let (ms1_loc, ms2_loc) = (max_loc, max_loc);
        let ref_time_ms = self.query_values.precursors.rts_ms[max_loc];

        let summed_ms1_int: f32 = match intensity_arrays.ms1_rtmajor.arr.get_row(ms1_loc) {
            Some(row) => row.iter().sum(),
            None => 0.0,
        };
        let summed_ms2_int: f32 = match intensity_arrays.ms2_rtmajor.arr.get_row(ms2_loc) {
            Some(row) => row.iter().sum(),
            None => 0.0,
        };
        let npeak_ms2 = match intensity_arrays.ms2_rtmajor.arr.get_row(ms2_loc) {
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

        let lazyscore_z = longitudinal_main_score_elements.ms2_lazyscore[max_loc] / norm_lazy_std;
        if lazyscore_z.is_nan() {
            let tmp = format!(
                "Lazy score is NaN {} and {}",
                longitudinal_main_score_elements.ms2_lazyscore[max_loc], norm_lazy_std
            );
            return Err(DataProcessingError::ExpectedFiniteNonNanData { context: tmp });
        }
        let raising_cycles = count_falling_steps(
            max_loc,
            -1,
            longitudinal_main_score_elements.ms2_lazyscore.as_slice(),
        );
        let falling_cycles = count_falling_steps(
            max_loc,
            1,
            longitudinal_main_score_elements.ms2_lazyscore.as_slice(),
        );

        Ok(MainScore {
            score: max_val,
            delta_next,
            delta_second_next,
            retention_time_ms: ref_time_ms,

            ms2_cosine_ref_sim: longitudinal_main_score_elements.ms2_cosine_ref_sim[max_loc],
            ms2_coelution_score: longitudinal_main_score_elements.ms2_coelution_score[max_loc],
            ms2_summed_intensity: summed_ms2_int,
            npeaks: npeak_ms2 as u8,
            lazyscore: longitudinal_main_score_elements.ms2_lazyscore[max_loc],
            lazyscore_vs_baseline: longitudinal_main_score_elements.ms2_lazyscore_vs_baseline
                [max_loc],
            lazyscore_z,
            ms2_corr_v_gauss: longitudinal_main_score_elements.ms2_corr_v_gauss[max_loc],
            ms1_corr_v_gauss: longitudinal_main_score_elements.ms1_corr_v_gauss[max_loc],
            // ms1_ms2_correlation: longitudinal_main_score_elements.ms1_ms2_correlation[max_loc],
            ms1_cosine_ref_sim: longitudinal_main_score_elements.ms1_cosine_ref_sim[max_loc],
            ms1_coelution_score: longitudinal_main_score_elements.ms1_coelution_score[max_loc],
            ms1_summed_intensity: summed_ms1_int,

            raising_cycles,
            falling_cycles,
        })
    }

    fn calc_main_score(&self) -> Result<MainScore, DataProcessingError> {
        let intensity_arrays =
            IntensityArrays::new(&self.query_values, &self.expected_intensities)?;
        let mut buffer = LongitudinalMainScoreElements::new_with_capacity(
            self.query_values.precursors.rts_ms.len(),
        );
        self.calc_with_intensities(&intensity_arrays, &mut buffer)
    }

    fn calc_with_inten_buffer(
        &self,
        intensity_arrays: &mut IntensityArrays,
        longit_buffer: &mut LongitudinalMainScoreElements,
    ) -> Result<MainScore, DataProcessingError> {
        intensity_arrays.reset_with(&self.query_values, &self.expected_intensities)?;
        self.calc_with_intensities(intensity_arrays, longit_buffer)
    }

    pub fn localize(self) -> Result<MainScore, DataProcessingError> {
        let main_score = self.calc_main_score()?;
        if main_score.score.is_nan() {
            // TODO find a way to nicely log the reason why some are nan.
            return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
        }
        Ok(main_score)
    }

    pub fn localize_with_buffer(
        &self,
        intensity_arrays: &mut IntensityArrays,
        longit_buffer: &mut LongitudinalMainScoreElements,
    ) -> Result<MainScore, DataProcessingError> {
        let main_score = self.calc_with_inten_buffer(intensity_arrays, longit_buffer)?;
        if main_score.score.is_nan() {
            // TODO find a way to nicely log the reason why some are nan.
            warn!("Main score is NaN");
            return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
        }
        Ok(main_score)
    }
}

/// The main score is meant to be a single number that is used to
/// identify the location of the 'best' candidate peptide in the
/// chromatogram.
#[derive(Debug, Clone, Copy)]
pub struct MainScore {
    pub score: f32,
    pub delta_next: f32,
    pub delta_second_next: f32,
    pub retention_time_ms: u32,

    pub ms2_cosine_ref_sim: f32,
    pub ms2_coelution_score: f32,
    pub ms2_summed_intensity: f32,
    pub npeaks: u8,
    pub lazyscore: f32,
    pub lazyscore_vs_baseline: f32,
    pub lazyscore_z: f32,
    pub ms2_corr_v_gauss: f32,
    pub ms1_corr_v_gauss: f32,

    pub ms1_cosine_ref_sim: f32,
    pub ms1_coelution_score: f32,
    pub ms1_summed_intensity: f32,

    /// Top N ions sorted in decreasing order of intensity
    /// at the apex of the chromatogram.
    // pub top_ions: [Option<IonAnnot>; 6],
    pub raising_cycles: u8,
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
