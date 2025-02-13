use crate::models::RTMajorIntensityArray;
use crate::utils::math::{
    lnfact,
    lnfact_f32,
};

pub fn peak_count(slices: &RTMajorIntensityArray, count_threshold: f32) -> Vec<u8> {
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

pub fn hyperscore(
    slices: &RTMajorIntensityArray,
    grouping: Option<&[u8]>,
    peak_count_cache: Option<&[u8]>,
) -> Vec<f32> {
    if peak_count_cache.is_some() {
        todo!("Implement peak count cache usage.")
    }
    slices
        .arr
        .row_apply(|slice| single_hyperscore(slice, grouping, 10.0))
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

pub fn lazyscore(slices: &RTMajorIntensityArray) -> Vec<f32> {
    slices.arr.row_apply(single_lazyscore).collect()
}
