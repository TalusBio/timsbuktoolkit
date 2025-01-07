use super::coelution::coelution_score;
use super::{
    corr_v_ref,
    hyperscore,
};
use crate::data_sources::speclib::{
    ExpectedIntensities,
    ReferenceEG,
};
use crate::errors::Result;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::models::{
    DecoyMarking,
    DigestSlice,
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
use crate::utils::aligning::snap_to_reference;
use crate::utils::top_n_array::TopNArray;
use serde::Serialize;
use std::hash::Hash;
use std::sync::Arc;
use timsquery::ElutionGroup;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::NaturalFinalizedMultiCMGArrays;

#[derive(Debug)]
pub struct PreScore<'a> {
    pub digest: &'a DigestSlice,
    pub charge: u8,
    pub reference: &'a ElutionGroup<SafePosition>,
    pub expected_intensities: &'a ExpectedIntensities,
    pub query_values: &'a NaturalFinalizedMultiCMGArrays<SafePosition>,
    pub ref_time_ms: Arc<[u32]>,
}

#[derive(Debug, Serialize, Clone)]
pub struct LongitudinalMainScoreElements {
    pub ms1_cosine_ref_sim: Vec<f32>,
    pub ms1_coelution_score: Vec<f32>,
    pub ms2_cosine_ref_sim: Vec<f32>,
    pub ms2_coelution_score: Vec<f32>,
    pub ms2_lazyscore: Vec<f32>,
    pub ref_time_ms: Arc<[u32]>,
}

#[derive(Debug)]
pub struct IntensityArrays {
    pub ms1_rtmajor: RTMajorIntensityArray,
    pub ms1_mzmajor: MzMajorIntensityArray,
    pub ms2_rtmajor: RTMajorIntensityArray,
    pub ms2_mzmajor: MzMajorIntensityArray,
    pub ms1_expected_intensities: Vec<f32>,
    pub ms2_expected_intensities: Vec<f32>,
}

impl IntensityArrays {
    pub fn new(query_values: &NaturalFinalizedMultiCMGArrays<SafePosition>, expected_intensities: &ExpectedIntensities) -> Self {
        let ms1_order: Vec<usize> = expected_intensities
            .precursor_intensities
            .iter()
            .enumerate()
            .map(|x| x.0)
            .collect();
        let (ms2_order, ms2_ref_vec): (Vec<SafePosition>, Vec<f32>) = expected_intensities
            .fragment_intensities
            .iter()
            .map(|(pos, intensity)| (*pos, *intensity as f32))
            .unzip();

        let ms1_rtmajor_arr =
            RTMajorIntensityArray::new(&query_values.ms1_arrays, Some(&ms1_order));
        let ms1_mzmajor_arr =
            MzMajorIntensityArray::new(&query_values.ms1_arrays, Some(&ms1_order));
        let ms2_rtmajor_arr =
            RTMajorIntensityArray::new(&query_values.ms2_arrays, Some(&ms2_order));
        let ms2_mzmajor_arr =
            MzMajorIntensityArray::new(&query_values.ms2_arrays, Some(&ms2_order));
        
        Self {
            ms1_rtmajor: ms1_rtmajor_arr,
            ms1_mzmajor: ms1_mzmajor_arr,
            ms2_rtmajor: ms2_rtmajor_arr,
            ms2_mzmajor: ms2_mzmajor_arr,
            ms1_expected_intensities: expected_intensities.precursor_intensities.clone(),
            ms2_expected_intensities: ms2_ref_vec,
        }
    }
}

impl LongitudinalMainScoreElements {
    pub fn new(intensity_arrays: &IntensityArrays, ref_time_ms: Arc<[u32]>, ms1_rts: &[u32], ms2_rts: &[u32]) -> Self {
        let ms1_cosine_ref_sim = snap_to_reference(
            &corr_v_ref::calculate_cosine_with_ref(
                &intensity_arrays.ms1_rtmajor,
                &intensity_arrays.ms1_expected_intensities,
            )
            .unwrap(),
            &ms1_rts,
            &ref_time_ms,
        )
        .unwrap();
        let ms2_cosine_ref_sim = snap_to_reference(
            &corr_v_ref::calculate_cosine_with_ref(&intensity_arrays.ms2_rtmajor, &intensity_arrays.ms2_expected_intensities).unwrap(),
            &ms2_rts,
            &ref_time_ms,
        )
        .unwrap();

        let ms2_coelution_score = snap_to_reference(
            &coelution_score::coelution_score(&intensity_arrays.ms2_mzmajor, 7),
            &ms2_rts,
            &ref_time_ms,
        )
        .unwrap();
        let ms1_coelution_score = snap_to_reference(
            &coelution_score::coelution_score(&intensity_arrays.ms1_mzmajor, 7),
            &ms1_rts,
            &ref_time_ms,
        )
        .unwrap();
        let lazyscore = snap_to_reference(
            &hyperscore::lazyscore(&intensity_arrays.ms2_rtmajor),
            &ms2_rts,
            &ref_time_ms,
        )
        .unwrap();

        Self {
            ms1_cosine_ref_sim,
            ms1_coelution_score,
            ms2_cosine_ref_sim,
            ms2_coelution_score,
            ms2_lazyscore: lazyscore,
            ref_time_ms,
        }
    }

    fn find_apex(&self) -> (usize, f32) {
        // IN THEORY all scores could be iterator and snapped on the fly,
        // Will do that later ... feels like a lot of work.
        // I would love to have this abstracted where multiple accumulators can be used
        // so I could have a top N accumulator that preserves the best locations, but could
        // also have an accumulator that keeps the score across the the full chromatogram
        let mut max_loc = 0;
        let mut max_val: f32 = 0.0;
        for i in 0..self.ref_time_ms.len() {
            let ms1_cos_score = 0.5 + (0.5 * self.ms1_cosine_ref_sim[i]);
            let loc_score =
                ms1_cos_score * self.ms2_cosine_ref_sim[i] * self.ms2_coelution_score[i] * self.ms2_lazyscore[i];
            if loc_score > max_val {
                max_val = loc_score;
                max_loc = i;
            }
        }

        (max_loc, max_val)
    }
}

impl<'a> PreScore<'a> {
    fn calc_main_score(&self) -> MainScore {
        let intensity_arrays = IntensityArrays::new(self.query_values, self.expected_intensities);
        let longitudinal_main_score_elements = LongitudinalMainScoreElements::new(
            &intensity_arrays,
            self.ref_time_ms.clone(),
            &self.query_values.ms1_arrays.retention_time_miliseconds,
            &self.query_values.ms2_arrays.retention_time_miliseconds,
        );

        let (max_loc, max_val) = longitudinal_main_score_elements.find_apex();
        
        let ms2_loc = self
            .query_values
            .ms2_arrays
            .retention_time_miliseconds
            .partition_point(|&x| x < self.ref_time_ms[max_loc]);

        let ms1_loc = self
            .query_values
            .ms1_arrays
            .retention_time_miliseconds
            .partition_point(|&x| x < self.ref_time_ms[max_loc]);

        let summed_ms1_int: f32 = intensity_arrays.ms1_rtmajor.arr.get_row(ms1_loc).iter().sum();
        let summed_ms2_int: f32 = intensity_arrays.ms2_rtmajor.arr.get_row(ms2_loc).iter().sum();
        let npeak_ms2 = intensity_arrays.ms2_rtmajor.arr.get_row(ms2_loc).iter().filter(|&x| *x > 10.0).count();

        let ims = self.query_values.ms2_arrays.weighted_ims_mean[ms2_loc];
        MainScore {
            score: max_val,
            observed_mobility: ims as f32,
            retention_time_ms: self.ref_time_ms[max_loc],

            ms2_cosine_ref_sim: longitudinal_main_score_elements.ms2_cosine_ref_sim[max_loc],
            ms2_coelution_score: longitudinal_main_score_elements.ms2_coelution_score[max_loc],
            ms2_summed_intensity: summed_ms2_int,
            npeaks: npeak_ms2 as u8,
            lazyscore: longitudinal_main_score_elements.ms2_lazyscore[max_loc],

            ms1_cosine_ref_sim: longitudinal_main_score_elements.ms1_cosine_ref_sim[max_loc],
            ms1_coelution_score: longitudinal_main_score_elements.ms1_coelution_score[max_loc],
            ms1_summed_intensity: summed_ms1_int,

            ref_ms1_idx: ms1_loc,
            ref_ms2_idx: ms2_loc,
        }
    }

    pub fn localize(self) -> LocalizedPreScore<'a> {
        let main_score = self.calc_main_score();
        LocalizedPreScore {
            pre_score: self,
            main_score: main_score,
        }
    }
}

/// The main score is meant to be a single number that is used to
/// identify the location of the 'best' candidate peptide in the
/// chromatogram.
#[derive(Debug, Clone, Copy)]
pub struct MainScore {
    pub score: f32,
    pub observed_mobility: f32,
    pub retention_time_ms: u32,

    pub ms2_cosine_ref_sim: f32,
    pub ms2_coelution_score: f32,
    pub ms2_summed_intensity: f32,
    pub npeaks: u8,
    pub lazyscore: f32,

    pub ms1_cosine_ref_sim: f32,
    pub ms1_coelution_score: f32,
    pub ms1_summed_intensity: f32,

    ref_ms1_idx: usize,
    ref_ms2_idx: usize,
}

#[derive(Debug)]
pub struct LocalizedPreScore<'a> {
    pub pre_score: PreScore<'a>,
    pub main_score: MainScore,
}

impl<'a> LocalizedPreScore<'a> {
    pub fn inten_sorted_errors(&self) -> SortedErrors {
        sorted_err_at_idx(
            self.main_score.ref_ms1_idx,
            self.main_score.ref_ms2_idx,
            &self.pre_score.query_values,
            &self.pre_score.reference,
        )
    }
}

#[derive(Debug)]
pub struct SortedErrors {
    pub ms1_mz_errors: [f32; 3],
    pub ms2_mz_errors: [f32; 7],
    pub ms1_mobility_errors: [f32; 3],
    pub ms2_mobility_errors: [f32; 7],
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

fn sorted_err_at_idx(
    ms1_idx: usize,
    ms2_idx: usize,
    cmgs: &NaturalFinalizedMultiCMGArrays<SafePosition>,
    elution_group: &ElutionGroup<SafePosition>,
) -> SortedErrors {
    // Get the elements at every index and sort them by intensity.
    // Once sorted calculate the pairwise diff.
    let mut ms1_elems: TopNArray<3, SortableError> = TopNArray::new();
    let mut ms2_elems: TopNArray<7, SortableError> = TopNArray::new();
    let ref_ims = elution_group.mobility;
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

    for (k, v) in cmgs.ms2_arrays.intensities.iter() {
        let expect_mz = elution_group.fragment_mzs.get(k);
        let mz_err = if let Some(mz) = expect_mz {
            (mz - cmgs.ms2_arrays.mz_means[&k][ms2_idx]) as f32
        } else {
            continue;
        };

        // TODO: figure out why this happens
        if mz_err.abs() > 1.0 {
            println!("Large mz diff for fragment {} is {}", k, mz_err.abs());
            println!("Expected mz: {}", expect_mz.unwrap());
            println!("Actual mz: {}", cmgs.ms2_arrays.mz_means[&k][ms2_idx]);
            println!("EG: {:#?}", elution_group);
            println!("CMGS: {:#?}", cmgs);
            panic!();
        }
        let ims_err = ref_ims - (cmgs.ms2_arrays.ims_means[&k][ms2_idx]) as f32;

        let tmp = SortableError {
            intensity: v[ms2_idx],
            mz_err,
            ims_err,
        };
        ms2_elems.push(tmp);
    }

    // Sort them by decreasing intensity
    let ms1_elems = ms1_elems.get_values();
    let ms2_elems = ms2_elems.get_values();

    let mut ms1_out_mz_err = [0.0; 3];
    let mut ms2_out_mz_err = [0.0; 7];
    let mut ms1_out_ims_err = [0.0; 3];
    let mut ms2_out_ims_err = [0.0; 7];

    for (i, val) in ms1_elems.iter().enumerate() {
        if i == 0 {
            ms1_out_mz_err[i] = val.mz_err;
            ms1_out_ims_err[i] = val.ims_err;
        } else {
            ms1_out_mz_err[i] = (val.mz_err - ms1_elems[i - 1].mz_err) as f32;
            ms1_out_ims_err[i] = (val.ims_err - ms1_elems[i - 1].ims_err) as f32;
        }
    }

    for (i, val) in ms2_elems.iter().enumerate() {
        if i == 0 {
            ms2_out_mz_err[i] = val.mz_err;
            ms2_out_ims_err[i] = val.ims_err;
        } else {
            ms2_out_mz_err[i] = (val.mz_err - ms2_elems[i - 1].mz_err) as f32;
            ms2_out_ims_err[i] = (val.ims_err - ms2_elems[i - 1].ims_err) as f32;
        }
    }

    SortedErrors {
        ms1_mz_errors: ms1_out_mz_err,
        ms2_mz_errors: ms2_out_mz_err,
        ms1_mobility_errors: ms1_out_ims_err,
        ms2_mobility_errors: ms2_out_ims_err,
    }
}
