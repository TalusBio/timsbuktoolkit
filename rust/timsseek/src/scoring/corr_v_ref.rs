use crate::errors::Result;
use crate::utils::correlation::cosine_similarity;
use timsquery::models::{
    MzMajorIntensityArray,
    RTMajorIntensityArray,
};
use timsquery::traits::key_like::KeyLike;

pub fn calculate_cosine_with_ref<FH: KeyLike>(
    slices: &RTMajorIntensityArray<FH>,
    ref_slice: &[f32],
) -> Result<Vec<f32>> {
    slices
        .arr
        .row_apply(|slice| cosine_similarity(slice, ref_slice))
        .collect()
}

// From https://doi.org/10.1101/2024.11.19.624419
const REF_GAUSSIAN: [f32; 7] = [0.0044, 0.054, 0.242, 0.399, 0.242, 0.054, 0.0044];
const REF_GAUSS_OFFSET: usize = 4;

fn slide_cosine_v_gaussian(slice: &[f32]) -> impl Iterator<Item = Result<f32>> + '_ {
    slice
        .windows(7)
        .map(|window| cosine_similarity(window, &REF_GAUSSIAN))
}

pub fn calculate_cosine_with_ref_gaussian<FH: KeyLike>(
    slices: &MzMajorIntensityArray<FH>,
) -> Result<Vec<f32>> {
    let mut result = vec![0.0; slices.arr.ncols()];
    slices
        .arr
        .row_apply(|row| {
            for (i, v) in slide_cosine_v_gaussian(row).enumerate() {
                match v {
                    Err(e) => return Err(e),
                    Ok(v) => {
                        if v.is_nan() {
                            continue;
                        }
                        result[i + REF_GAUSS_OFFSET] += v.max(0.0);
                    }
                }
            }
            Ok(())
        })
        .collect::<Result<()>>()?;

    let nrows = slices.arr.nrows();
    result.iter_mut().for_each(|v| *v /= nrows as f32);

    Ok(result)
}

#[test]
fn test_calculate_cosine_with_ref_gaussian() {
    let test_vec = (0..3)
        .into_iter()
        .map(|_| REF_GAUSSIAN)
        .flatten()
        .collect::<Vec<f32>>();
    let out = slide_cosine_v_gaussian(&test_vec).collect::<Result<Vec<f32>>>();
    let expect_out = [
        1.0000001,
        0.7786916,
        0.36945745,
        0.122937046,
        0.122937046,
        0.3694575,
        0.77869165,
        1.0000001,
        0.7786916,
        0.36945745,
        0.122937046,
        0.122937046,
        0.3694575,
        0.77869165,
        1.0000001,
    ];
    assert!(out.is_ok());

    let out = out.unwrap();
    println!("{:?}", out);
    for (a, b) in out.iter().zip(expect_out.iter()) {
        assert!((a - b).abs() < 1e-6);
    }
}
