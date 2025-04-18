use std::cmp::Ordering;

use crate::fragment_mass::{
    IonAnnot,
    IonSeriesTerminality,
};
use timsquery::models::RTMajorIntensityArray;
use crate::utils::math::{
    lnfact,
    lnfact_f32,
};
use timsquery::traits::key_like::KeyLike;

pub fn peak_count<FH: KeyLike>(
    slices: &RTMajorIntensityArray<FH>,
    count_threshold: f32,
) -> Vec<u8> {
    slices
        .arr
        .row_apply(|slice| {
            let mut count: u8 = 0;
            for intensity in slice {
                if *intensity > count_threshold {
                    count += 1;
                }
            }
            count
        })
        .collect()
}

/// From: PMC5409104
///
/// `log(Nb! * Ny! * (Sum Intensity b) * (Sum Intensity y))`
///
/// Example:
/// ```
/// use timsseek::scoring::hyperscore::single_hyperscore;
/// let slice = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
/// let grouping = vec![0, 1, 2, 0, 1, 2];
/// let count_threshold = 0.0;
/// let score = single_hyperscore(&slice, Some(&grouping), count_threshold);
/// assert_eq!(score, 8.747712863479158);
/// ```
pub fn single_hyperscore(slice: &[f32], grouping: Option<&[u8]>, count_threshold: f32) -> f32 {
    let mut score = 0.0;
    let min_grouping: u8 = grouping.map_or(0, |x| *x.iter().min().unwrap_or(&0));
    let max_grouping: u8 = grouping.map_or(0, |x| *x.iter().max().unwrap_or(&0));
    for local_group in min_grouping..=max_grouping {
        let mut local_sum = 0.0;
        let mut local_count: u16 = 0;
        for (i, intensity) in slice.iter().enumerate() {
            if (local_group == grouping.map_or(0, |x| x[i])) && *intensity > count_threshold {
                local_sum += intensity;
                local_count += 1;
            }
        }

        // Addition of logs is the same as log of products
        score += lnfact(local_count) as f32;
        score += local_sum.ln();
    }

    score
}

pub fn single_hyperscore_labeled(slice: &[f32], labels: &[IonAnnot], count_threshold: f32) -> f32 {
    let mut nt_count = 0;
    let mut ct_count = 0;

    let mut nt_sum = 0.0;
    let mut ct_sum = 0.0;

    slice.iter().zip(labels).for_each(|(inten, lab)| {
        // Not smaller is diff than bigger because of Nans
        match inten.partial_cmp(&count_threshold) {
            Some(Ordering::Less) => return,
            Some(Ordering::Equal) => return,
            Some(Ordering::Greater) => {}
            None => return,
        }

        match lab.terminality() {
            IonSeriesTerminality::CTerm => {
                ct_count += 1;
                ct_sum += inten;
            }
            IonSeriesTerminality::NTerm => {
                nt_count += 1;
                nt_sum += inten;
            }
            // Unfragmented precursors dont get added to the score
            IonSeriesTerminality::None => {}
        }
    });

    let score =
        (nt_sum + ct_sum).ln_1p() + lnfact(ct_count as u16) as f32 + lnfact(nt_count as u16) as f32;
    if score.is_finite() { score } else { 255.0 }
}

pub fn hyperscore(slices: &RTMajorIntensityArray<IonAnnot>) -> Vec<f32> {
    slices
        .arr
        .row_apply(|slice| single_hyperscore_labeled(slice, &slices.order, 10.0))
        .collect()
}

pub fn calculate_lazy_hyperscore(npeaks: &[u8], summed_intensity: &[u64]) -> Vec<f32> {
    let mut scores = vec![0.0; npeaks.len()];
    for i in 0..npeaks.len() {
        let npeaks_i = npeaks[i];
        let summed_intensity_i = summed_intensity[i];
        let log1p_intensities_i = (summed_intensity_i as f32 + 1.0).ln();
        scores[i] = lnfact(npeaks_i as u16) as f32 + (2.0 * log1p_intensities_i);
    }
    scores
}

fn single_lazyscore(slc: &[f32]) -> f32 {
    lnfact_f32(slc.iter().map(|&x| x.max(1.0).ln()).sum())
}

fn single_split_ion_lazyscore(slc: &[f32], labels: &[IonAnnot]) -> f32 {
    let mut ct_lnsum = 0.0;
    let mut nt_lnsum = 0.0;

    for (i, label) in labels.iter().enumerate() {
        match label.terminality() {
            IonSeriesTerminality::CTerm => {
                ct_lnsum += slc[i].max(1.0).ln();
            }
            IonSeriesTerminality::NTerm => {
                nt_lnsum += slc[i].max(1.0).ln();
            }
            // Unfragmented precursors dont get added to the score
            IonSeriesTerminality::None => {}
        }
    }

    lnfact_f32(ct_lnsum) + lnfact_f32(nt_lnsum)
}

pub fn lazyscore<K: Clone>(slices: &RTMajorIntensityArray<K>) -> Vec<f32> {
    slices.arr.row_apply(single_lazyscore).collect()
}

pub fn split_ion_lazyscore(slices: &RTMajorIntensityArray<IonAnnot>) -> Vec<f32> {
    slices
        .arr
        .row_apply(|slc| single_split_ion_lazyscore(slc, &slices.order))
        .collect()
}
