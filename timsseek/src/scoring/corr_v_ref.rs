use crate::errors::Result;
use crate::models::RTMajorIntensityArray;
use crate::utils::correlation::cosine_similarity;

pub fn calculate_cosine_with_ref(
    slices: &RTMajorIntensityArray,
    ref_slice: &[f32],
) -> Result<Vec<f32>> {
    slices
        .arr
        .row_apply(|slice| cosine_similarity(slice, ref_slice))
        .into_iter()
        .collect::<Result<Vec<f32>>>()
}
