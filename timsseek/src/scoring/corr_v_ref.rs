use crate::errors::Result;
use crate::models::RTMajorIntensityArray;
use crate::utils::correlation::cosine_similarity;
use serde::Serialize;
use std::hash::Hash;

pub fn calculate_cosine_with_ref<FH: Clone + Eq + Serialize + Hash + Send + Sync>(
    slices: &RTMajorIntensityArray<FH>,
    ref_slice: &[f32],
) -> Result<Vec<f32>> {
    slices
        .arr
        .row_apply(|slice| cosine_similarity(slice, ref_slice))
        .into_iter()
        .collect::<Result<Vec<f32>>>()
}
