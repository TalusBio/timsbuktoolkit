use crate::errors::DataProcessingError;
use crate::utils::correlation::cosine_similarity;
use timsquery::models::MzMajorIntensityArray;
use timsquery::traits::key_like::KeyLike;

// From https://doi.org/10.1101/2024.11.19.624419
const REF_GAUSSIAN: [f32; 7] = [0.0044, 0.054, 0.242, 0.399, 0.242, 0.054, 0.0044];
const REF_GAUSS_OFFSET: usize = 4;

// Note that padding is not hadled here, so the output will be smaller than the input.
fn slide_cosine_v_gaussian(
    slice: &[f32],
) -> impl Iterator<Item = Result<f32, DataProcessingError>> + '_ {
    slice
        .windows(7)
        .map(|window| cosine_similarity(window, &REF_GAUSSIAN))
}

pub fn calculate_cosine_with_ref_gaussian_into<FH: KeyLike>(
    slices: &MzMajorIntensityArray<FH, f32>,
    filter: impl Fn(&FH) -> bool,
    buffer: &mut Vec<f32>,
) -> Result<(), DataProcessingError> {
    buffer.clear();
    buffer.resize(slices.arr.ncols(), 0.0);

    let ratio = 1f32 / slices.arr.nrows() as f32;
    slices.iter_mzs().try_for_each(|((k, _mz), slc)| {
        if !filter(k) {
            return Ok(());
        }
        for (i, v) in slide_cosine_v_gaussian(slc).enumerate() {
            match v {
                Err(e) => return Err(e),
                Ok(v) => {
                    if v.is_nan() {
                        continue;
                    }
                    buffer[i + REF_GAUSS_OFFSET] += v.max(0.0) * ratio;
                }
            }
        }
        Ok(())
    })
}

// pub fn calculate_cosine_with_ref_gaussian<FH: KeyLike>(
//     slices: &MzMajorIntensityArray<FH, f32>,
//     filter: impl Fn(&FH) -> bool,
// ) -> Result<Vec<f32>, DataProcessingError> {
//     let mut result = vec![0.0; slices.arr.ncols()];
//     calculate_cosine_with_ref_gaussian_into(slices, filter, &mut result)?;
//
//     Ok(result)
// }

#[test]
fn test_calculate_cosine_with_ref_gaussian() {
    let test_vec = (0..3).flat_map(|_| REF_GAUSSIAN).collect::<Vec<f32>>();
    let out = slide_cosine_v_gaussian(&test_vec).collect::<Result<Vec<f32>, DataProcessingError>>();
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
