use super::coelution::coelution_score;
use super::{
    corr_v_ref,
    hyperscore,
};
use crate::data_sources::speclib::ExpectedIntensities;
use crate::errors::{
    DataProcessingError,
    Result,
};
use crate::fragment_mass::IonAnnot;
use crate::models::{
    DigestSlice,
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
use crate::utils::rolling_calculators::{
    calculate_centered_std,
    calculate_value_vs_baseline,
};
use crate::utils::top_n_array::TopNArray;
use core::f32;
use serde::Serialize;
use std::sync::Arc;
use timsquery::ElutionGroup;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::NaturalFinalizedMultiCMGArrays;

#[derive(Debug)]
pub struct PreScore {
    pub digest: DigestSlice,
    pub charge: u8,
    pub reference: ElutionGroup<IonAnnot>,
    pub expected_intensities: ExpectedIntensities,
    pub query_values: NaturalFinalizedMultiCMGArrays<IonAnnot>,
    pub ref_time_ms: Arc<[u32]>,
}

#[derive(Debug, Serialize, Clone)]
pub struct LongitudinalMainScoreElements {
    pub ms1_cosine_ref_sim: Vec<f32>,
    pub ms1_coelution_score: Vec<f32>,
    pub ms2_cosine_ref_sim: Vec<f32>,
    pub ms2_coelution_score: Vec<f32>,
    pub ms2_lazyscore: Vec<f32>,
    pub ms2_lazyscore_vs_baseline: Vec<f32>,
    // TODO: REMOVE
    pub hyperscore: Vec<f32>,
    pub split_lazyscore: Vec<f32>,

    /// END
    pub ref_time_ms: Arc<[u32]>,
    ms2_lazyscore_vs_baseline_std: f32,
}

#[derive(Debug)]
pub struct IntensityArrays {
    pub ms1_rtmajor: RTMajorIntensityArray<usize>,
    pub ms1_mzmajor: MzMajorIntensityArray<usize>,
    pub ms2_rtmajor: RTMajorIntensityArray<IonAnnot>,
    pub ms2_mzmajor: MzMajorIntensityArray<IonAnnot>,
    pub ms1_expected_intensities: Vec<f32>,
    pub ms2_expected_intensities: Vec<f32>,
}

impl IntensityArrays {
    pub fn new(
        query_values: &NaturalFinalizedMultiCMGArrays<IonAnnot>,
        expected_intensities: &ExpectedIntensities,
    ) -> Result<Self> {
        let ms1_order: Arc<[usize]> = expected_intensities
            .precursor_intensities
            .iter()
            .enumerate()
            .map(|x| x.0)
            .collect();
        let (ms2_order, ms2_ref_vec): (Vec<IonAnnot>, Vec<f32>) = expected_intensities
            .fragment_intensities
            .iter()
            .map(|(pos, intensity)| (*pos, { *intensity }))
            .unzip();
        let ms2_order: Arc<[IonAnnot]> = ms2_order.into();

        let ms1_rtmajor_arr =
            RTMajorIntensityArray::new(&query_values.ms1_arrays, ms1_order.clone());
        let ms1_mzmajor_arr =
            MzMajorIntensityArray::new(&query_values.ms1_arrays, ms1_order.clone());
        let ms2_rtmajor_arr =
            RTMajorIntensityArray::new(&query_values.ms2_arrays, ms2_order.clone());
        let ms2_mzmajor_arr =
            MzMajorIntensityArray::new(&query_values.ms2_arrays, ms2_order.clone());

        match (
            ms1_rtmajor_arr,
            ms1_mzmajor_arr,
            ms2_rtmajor_arr,
            ms2_mzmajor_arr,
        ) {
            (
                Ok(ms1_rtmajor_arr),
                Ok(ms1_mzmajor_arr),
                Ok(ms2_rtmajor_arr),
                Ok(ms2_mzmajor_arr),
            ) => Ok(Self {
                ms1_rtmajor: ms1_rtmajor_arr,
                ms1_mzmajor: ms1_mzmajor_arr,
                ms2_rtmajor: ms2_rtmajor_arr,
                ms2_mzmajor: ms2_mzmajor_arr,
                ms1_expected_intensities: expected_intensities.precursor_intensities.clone(),
                ms2_expected_intensities: ms2_ref_vec,
            }),
            (Err(e), _, _, _) | (_, Err(e), _, _) => {
                let e = e.append_to_context("Failed to create IntensityArrays for MS1, ");
                Err(e.into())
            }
            (_, _, Err(e), _) | (_, _, _, Err(e)) => {
                let e = e.append_to_context("Failed to create IntensityArrays for MS2, ");
                Err(e.into())
            }
        }
    }
}

fn gaussblur(x: &mut [f32]) {
    // Temp implementation ... shoudl make something nicer in the future.
    let len = x.len();
    if len < 3 {
        return;
    }

    // Using fixed kernel weights [0.5, 1.0, 0.5]
    // Note: These weights are already normalized (sum = 2)
    let mut temp = vec![0.0; len];

    // Handle first element
    temp[0] = (x[0] * 1.5 + x[1] * 0.5) / 2.0;

    // Main convolution loop
    for i in 1..len - 1 {
        temp[i] = (x[i - 1] * 0.5 + x[i] * 1.0 + x[i + 1] * 0.5) / 2.0;
    }

    // Handle last element
    temp[len - 1] = (x[len - 1] * 1.5 + x[len - 2] * 0.5) / 2.0;

    // Copy results back
    x.copy_from_slice(&temp);
}

impl LongitudinalMainScoreElements {
    pub fn new(intensity_arrays: &IntensityArrays, ref_time_ms: Arc<[u32]>) -> Result<Self> {
        let mut ms1_cosine_ref_sim = corr_v_ref::calculate_cosine_with_ref(
            &intensity_arrays.ms1_rtmajor,
            &intensity_arrays.ms1_expected_intensities,
        )
        .unwrap();
        let mut ms2_cosine_ref_sim = corr_v_ref::calculate_cosine_with_ref(
            &intensity_arrays.ms2_rtmajor,
            &intensity_arrays.ms2_expected_intensities,
        )
        .unwrap();

        // In this section a "insufficient data error" will be returned if not enough
        // data exists to calculate a score.
        let ms2_coe_scores =
            coelution_score::coelution_score::<10, IonAnnot>(&intensity_arrays.ms2_mzmajor, 7);
        if let Err(e) = ms2_coe_scores {
            let e = e.append_to_context("Failed to calculate coelution score for MS2, ");
            return Err(crate::errors::TimsSeekError::DataProcessingError(e));
        }
        let mut ms2_coelution_score = ms2_coe_scores.unwrap();

        let tmp = coelution_score::coelution_score::<6, usize>(&intensity_arrays.ms1_mzmajor, 7);
        let mut ms1_coelution_score = match tmp {
            Ok(scores) => scores,
            Err(_) => vec![0.0; ref_time_ms.len()],
        };
        let mut lazyscore = hyperscore::lazyscore(&intensity_arrays.ms2_rtmajor);

        let mut hyperscore = hyperscore::hyperscore(&intensity_arrays.ms2_rtmajor);
        let mut split_lazyscore = hyperscore::split_ion_lazyscore(&intensity_arrays.ms2_rtmajor);

        gaussblur(&mut lazyscore);
        gaussblur(&mut ms1_coelution_score);
        gaussblur(&mut ms2_coelution_score);
        gaussblur(&mut ms2_cosine_ref_sim);
        gaussblur(&mut ms1_cosine_ref_sim);

        let five_pct_index = ref_time_ms.len() * 5 / 100;
        let half_five_pct_idnex = five_pct_index / 2;
        let lazyscore_vs_baseline = calculate_value_vs_baseline(&lazyscore, five_pct_index);
        let lzb_std = calculate_centered_std(
            &lazyscore_vs_baseline
                [(half_five_pct_idnex)..(&lazyscore_vs_baseline.len() - half_five_pct_idnex)],
        );

        Ok(Self {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore: lazyscore,
            ms2_lazyscore_vs_baseline: lazyscore_vs_baseline,
            split_lazyscore,
            hyperscore,
            ref_time_ms,
            ms2_lazyscore_vs_baseline_std: lzb_std,
        })
    }

    fn find_apex_candidates(&self) -> [ScoreInTime; 20] {
        let mut candidate_groups: [TopNArray<2, ScoreInTime>; 21] = [TopNArray::new(); 21];
        let five_pct_index = self.ref_time_ms.len() * 5 / 100;
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
        (0..self.ref_time_ms.len()).map(|i| {
            // The 0.75 - 0.25 means we are downscaling the main lazyscore up to 0.75
            // since the similarity is in the 0-1 range; even if the precursor has
            // similarity of 0, we still have a scoring value.

            let ms1_cos_score = 0.75 + (0.25 * self.ms1_cosine_ref_sim[i].powi(2));
            // let ms1_cos_score = self.ms1_cosine_ref_sim[i].powi(2);
            // let mut loc_score = self.ms2_lazyscore[i];
            let mut loc_score = self.ms2_lazyscore_vs_baseline[i];
            // let mut loc_score =  self.ms2_lazyscore_vs_baseline[i] / self.ms2_lazyscore_vs_baseline_std;
            loc_score *= ms1_cos_score;
            loc_score *= self.ms2_cosine_ref_sim[i].powi(2);
            loc_score *= self.ms2_coelution_score[i].powi(2);
            loc_score
        })
    }
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
    fn calc_main_score(&self) -> Result<MainScore> {
        let intensity_arrays =
            IntensityArrays::new(&self.query_values, &self.expected_intensities)?;
        let longitudinal_main_score_elements =
            LongitudinalMainScoreElements::new(&intensity_arrays, self.ref_time_ms.clone())?;

        let apex_candidates = longitudinal_main_score_elements.find_apex_candidates();
        let norm_lazy_std =
            calculate_centered_std(&longitudinal_main_score_elements.ms2_lazyscore_vs_baseline);
        let max_val = apex_candidates[0].score;
        let max_loc = apex_candidates[0].index;

        // TODO: Branch here if the results are empty.

        // This is a delta next with the constraint that it has to be more than 5% of the max
        // index apart from the max.
        let ten_pct_index = self.ref_time_ms.len() / 20;
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

        let ms2_loc = self
            .query_values
            .ms2_arrays
            .retention_time_miliseconds
            .partition_point(|&x| x <= self.ref_time_ms[max_loc]);

        let ms1_loc = self
            .query_values
            .ms1_arrays
            .retention_time_miliseconds
            .partition_point(|&x| x <= self.ref_time_ms[max_loc]);

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
                    if *x > 10.0f32 {
                        count += 1;
                    }
                }
                count
            }
            None => 0,
        };

        let ims_ms2 = self
            .query_values
            .ms2_arrays
            .weighted_ims_mean
            .get(ms2_loc)
            .copied()
            .unwrap_or(f64::NAN);

        let ims_ms1 = self
            .query_values
            .ms1_arrays
            .weighted_ims_mean
            .get(ms1_loc)
            .copied()
            .unwrap_or(f64::NAN);

        let ims = match (ims_ms1, ims_ms2) {
            (x, y) if !x.is_nan() && y.is_nan() => x,
            (x, y) if x.is_nan() && !y.is_nan() => y,
            (x, y) if x.is_nan() && y.is_nan() => f64::NAN,
            (x, y) => (x + y) / 2.0,
        };

        Ok(MainScore {
            score: max_val,
            delta_next,
            delta_second_next,
            observed_mobility: ims as f32,
            observed_mobility_ms1: ims_ms1 as f32,
            observed_mobility_ms2: ims_ms2 as f32,
            retention_time_ms: self.ref_time_ms[max_loc],

            ms2_cosine_ref_sim: longitudinal_main_score_elements.ms2_cosine_ref_sim[max_loc],
            ms2_coelution_score: longitudinal_main_score_elements.ms2_coelution_score[max_loc],
            ms2_summed_intensity: summed_ms2_int,
            npeaks: npeak_ms2 as u8,
            lazyscore: longitudinal_main_score_elements.ms2_lazyscore[max_loc],
            lazyscore_vs_baseline: longitudinal_main_score_elements.ms2_lazyscore_vs_baseline
                [max_loc],
            lazyscore_z: longitudinal_main_score_elements.ms2_lazyscore[max_loc] / norm_lazy_std,
            // ms1_ms2_correlation: longitudinal_main_score_elements.ms1_ms2_correlation[max_loc],
            ms1_cosine_ref_sim: longitudinal_main_score_elements.ms1_cosine_ref_sim[max_loc],
            ms1_coelution_score: longitudinal_main_score_elements.ms1_coelution_score[max_loc],
            ms1_summed_intensity: summed_ms1_int,

            ref_ms1_idx: ms1_loc,
            ref_ms2_idx: ms2_loc,
        })
    }

    pub fn localize(self) -> Result<LocalizedPreScore> {
        let main_score = self.calc_main_score()?;
        if main_score.score.is_nan() {
            // TODO find a way to nicely log the reason why some are nan.
            return Err(DataProcessingError::ExpectedNonEmptyData { context: None }.into());
        }
        Ok(LocalizedPreScore::new(self, main_score))
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
    pub observed_mobility: f32,
    pub observed_mobility_ms1: f32,
    pub observed_mobility_ms2: f32,
    pub retention_time_ms: u32,

    pub ms2_cosine_ref_sim: f32,
    pub ms2_coelution_score: f32,
    pub ms2_summed_intensity: f32,
    pub npeaks: u8,
    pub lazyscore: f32,
    pub lazyscore_vs_baseline: f32,
    pub lazyscore_z: f32,

    pub ms1_cosine_ref_sim: f32,
    pub ms1_coelution_score: f32,
    pub ms1_summed_intensity: f32,

    ref_ms1_idx: usize,
    ref_ms2_idx: usize,
}

#[derive(Debug)]
pub struct LocalizedPreScore {
    pub pre_score: PreScore,
    pub main_score: MainScore,
    ints_at_apex: SortedIntElemAtIndex,
}

impl LocalizedPreScore {
    pub fn inten_sorted_errors(&self) -> SortedErrors {
        sorted_err_at_idx(&self.ints_at_apex)
    }

    pub fn relative_intensities(&self) -> RelativeIntensities {
        RelativeIntensities::new(&self.ints_at_apex)
    }
}

impl LocalizedPreScore {
    fn new(pre_score: PreScore, main_score: MainScore) -> Self {
        let ints_at_apex = SortedIntElemAtIndex::new(
            main_score.ref_ms1_idx,
            main_score.ref_ms2_idx,
            &pre_score.query_values,
            &pre_score.reference,
        );
        Self {
            pre_score,
            main_score,
            ints_at_apex,
        }
    }
}

#[derive(Debug)]
pub struct SortedErrors {
    pub ms1_mz_errors: [f32; 3],
    pub ms2_mz_errors: [f32; 7],
    pub ms1_mobility_errors: [f32; 3],
    pub ms2_mobility_errors: [f32; 7],
}

#[derive(Debug)]
pub struct SortedIntElemAtIndex {
    ms1: [SortableError; 3],
    ms2: [SortableError; 7],
}

impl SortedIntElemAtIndex {
    fn new(
        ms1_idx: usize,
        ms2_idx: usize,
        cmgs: &NaturalFinalizedMultiCMGArrays<IonAnnot>,
        elution_group: &ElutionGroup<IonAnnot>,
    ) -> Self {
        // Get the elements at every index and sort them by intensity.
        // Once sorted calculate the pairwise diff.
        let mut ms1_elems: TopNArray<3, SortableError> = TopNArray::new();
        let mut ms2_elems: TopNArray<7, SortableError> = TopNArray::new();
        let ref_ims = elution_group.mobility;

        if ms1_idx < cmgs.ms1_arrays.retention_time_miliseconds.len() {
            for i in 0..3 {
                let expect_mz = elution_group.precursor_mzs.get(i);
                let mz_err = if let Some(mz) = expect_mz {
                    (mz - cmgs.ms1_arrays.mz_means[&i][ms1_idx]) as f32
                } else {
                    continue;
                };
                let ims_err = ref_ims - (cmgs.ms1_arrays.ims_means[&i][ms1_idx]) as f32;

                let tmp = SortableError {
                    intensity: cmgs.ms1_arrays.intensities[&i][ms1_idx],
                    mz_err,
                    ims_err,
                };
                ms1_elems.push(tmp);
            }
        }

        // TODO: make an impl to get the length on the RT axis.
        if ms2_idx < cmgs.ms2_arrays.retention_time_miliseconds.len() {
            for (k, v) in cmgs.ms2_arrays.intensities.iter() {
                let expect_mz = elution_group.fragment_mzs.get(k);
                let mz_err = if let Some(mz) = expect_mz {
                    (mz - cmgs.ms2_arrays.mz_means[k][ms2_idx]) as f32
                } else {
                    continue;
                };

                // TODO: figure out why this happens
                if mz_err.abs() > 1.0 {
                    println!("Large mz diff for fragment {} is {}", k, mz_err.abs());
                    println!("Expected mz: {}", expect_mz.unwrap());
                    println!("Actual mz: {}", cmgs.ms2_arrays.mz_means[k][ms2_idx]);
                    println!("EG: {:#?}", elution_group);
                    println!("CMGS: {:#?}", cmgs);
                    panic!();
                }
                let ims_err = ref_ims - (cmgs.ms2_arrays.ims_means[k][ms2_idx]) as f32;

                let tmp = SortableError {
                    intensity: v[ms2_idx],
                    mz_err,
                    ims_err,
                };
                ms2_elems.push(tmp);
            }
        }

        SortedIntElemAtIndex {
            ms1: ms1_elems.get_values(),
            ms2: ms2_elems.get_values(),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct SortableError {
    intensity: u64,
    mz_err: f32,
    ims_err: f32,
}

impl PartialOrd for SortableError {
    fn partial_cmp(&self, other: &SortableError) -> Option<std::cmp::Ordering> {
        self.intensity.partial_cmp(&other.intensity)
    }
}

impl Default for SortableError {
    fn default() -> Self {
        Self {
            intensity: 0,
            mz_err: f32::NAN,
            ims_err: f32::NAN,
        }
    }
}

fn sorted_err_at_idx(ints_at_apex: &SortedIntElemAtIndex) -> SortedErrors {
    let mut ms1_out_mz_err = [0.0; 3];
    let mut ms2_out_mz_err = [0.0; 7];
    let mut ms1_out_ims_err = [0.0; 3];
    let mut ms2_out_ims_err = [0.0; 7];

    let ms1_elems = ints_at_apex.ms1;
    let ms2_elems = ints_at_apex.ms2;

    for (i, val) in ms1_elems.iter().enumerate() {
        if i == 0 {
            ms1_out_mz_err[i] = val.mz_err;
            ms1_out_ims_err[i] = val.ims_err;
        } else {
            ms1_out_mz_err[i] = val.mz_err - ms1_elems[i - 1].mz_err;
            ms1_out_ims_err[i] = val.ims_err - ms1_elems[i - 1].ims_err;
        }
    }

    for (i, val) in ms2_elems.iter().enumerate() {
        if i == 0 {
            ms2_out_mz_err[i] = val.mz_err;
            ms2_out_ims_err[i] = val.ims_err;
        } else {
            ms2_out_mz_err[i] = val.mz_err - ms2_elems[i - 1].mz_err;
            ms2_out_ims_err[i] = val.ims_err - ms2_elems[i - 1].ims_err;
        }
    }

    SortedErrors {
        ms1_mz_errors: ms1_out_mz_err,
        ms2_mz_errors: ms2_out_mz_err,
        ms1_mobility_errors: ms1_out_ims_err,
        ms2_mobility_errors: ms2_out_ims_err,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RelativeIntensities {
    pub ms1: [f32; 3],
    pub ms2: [f32; 7],
}

impl RelativeIntensities {
    pub fn new(ints_at_apex: &SortedIntElemAtIndex) -> Self {
        let tot_ms1: u64 = ints_at_apex.ms1.iter().map(|x| x.intensity).sum();
        let tot_ms2: u64 = ints_at_apex.ms2.iter().map(|x| x.intensity).sum();
        let mut ms1: [f32; 3] = [0.0; 3];
        let mut ms2: [f32; 7] = [0.0; 7];
        ints_at_apex
            .ms1
            .iter()
            .enumerate()
            .for_each(|(i, x)| ms1[i] = (x.intensity as f32).ln_1p() - (tot_ms1 as f32).ln_1p());
        ints_at_apex
            .ms2
            .iter()
            .enumerate()
            .for_each(|(i, x)| ms2[i] = (x.intensity as f32).ln_1p() - (tot_ms2 as f32).ln_1p());

        Self { ms1, ms2 }
    }
}
