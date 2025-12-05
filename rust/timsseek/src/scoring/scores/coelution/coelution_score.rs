use crate::errors::DataProcessingError;
use crate::utils::correlation::rolling_cosine_similarity;
use timsquery::models::{
    Array2D,
    MzMajorIntensityArray,
};
use tracing::trace;

// /// Calculates the coelution score of a set of chromatograms.
// ///
// /// # Example
// ///
// /// ```
// /// use timsquery::Array2D;
// /// use timsseek::scoring::coelution;
// ///
// /// let slices = Array2D::new(
// ///     vec![[0., 1., 3., 22., 5.],
// ///          [0., 2., 4., 20., 5.],
// ///          [0., 1., 2., 19., 2.]]).unwrap();
// /// let window = 3;
// /// // Note that the generic type parameter is the top N of scores that will
// /// // be averaged to report the coelution.
// /// let scores = coelution::coelution_score_arr::<3>(&slices, window).unwrap();
// /// assert_eq!(scores, [0.0, 0.9866667, 0.9939658, 0.9849558, 0.0]);
// /// ```
// pub fn coelution_score_arr<'a, const TOP_N: usize>(
//     slices: &'a Array2D<f32>,
//     window_size: usize,
// ) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
//     let filter: Option<fn(usize) -> bool> = None;
//     coelution_score_iter_filter::<TOP_N>(slices, window_size, filter)
// }

// Assuming these types are in the current scope:
// use crate::{Array2D, DataProcessingError, TopNArray, rolling_cosine_similarity};

/// Calculates the coelution score of a set of chromatograms, returning a lazy iterator.
fn coelution_vref_score_filter_onto(
    slices: &Array2D<f32>,
    ref_slice: &[f32],
    window_size: usize,
    filter: impl Fn(usize) -> bool,
    buffer: &mut Vec<f32>,
) -> Result<(), DataProcessingError> {
    if slices.ncols() < window_size {
        trace!("Not enough data to calculate coelution score");
        return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
    }

    let num_elems = (0..slices.nrows()).filter(|&i| filter(i)).count();
    let norm_factor = 1f32 / (num_elems as f32).max(1.0);
    if norm_factor == 1.0 {
        trace!("No valid slices after filtering");
        return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
    }
    if num_elems > 50 {
        trace!(
            "There are too many valid slices after filtering, probably an mz-major and an rt-major array got mixed up"
        );
        return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
    }
    buffer.clear();
    buffer.resize(slices.ncols(), 0.0);

    // Collect all rolling similarity calculations v the reference.
    // This still computes the similarities upfront, but the aggregation into the final
    // score is done lazily, one time-point at a time.
    let res: Result<(), DataProcessingError> = (0..slices.nrows())
        .filter(|&i| filter(i))
        .try_for_each(|i| {
            let slice1 = slices.get_row(i).expect("Row index i is within bounds");
            let iter = rolling_cosine_similarity(slice1, ref_slice, window_size)?;
            for (i, v) in iter.enumerate() {
                if v.is_nan() {
                    continue;
                }
                buffer[i] += v.max(0.0)
            }
            Ok(())
        });
    res?;

    for x in buffer.iter_mut() {
        *x *= norm_factor;
    }
    Ok(())
}

/// Calculates the coelution score for a set of chromatograms against a reference slice.
///
/// This function is a variant of `coelution_vref_score_filter_onto` that works with
/// an `MzMajorIntensityArray`. It uses the m/z order within the `slices` to filter
/// which chromatograms to include in the score calculation.
///
/// # Arguments
///
/// * `slices` - An `MzMajorIntensityArray` containing the intensity data. The rows of this
///   array correspond to different m/z values, and the columns correspond to
///   different time points or cycles. The `mz_order` field of this struct provides
///   the mapping from row index to m/z value.
/// * `ref_slice` - A slice representing the reference chromatogram.
/// * `window` - The size of the rolling window for the cosine similarity calculation.
/// * `filter` - A closure that takes a key (of type `K`) from the `mz_order` and returns
///   `true` if the corresponding chromatogram should be included in the calculation.
/// * `buffer` - A mutable buffer to store the resulting coelution scores.
pub fn coelution_vref_score_filter_into<'a, K: Clone + Ord>(
    slices: &'a MzMajorIntensityArray<K, f32>,
    ref_slice: &'a [f32],
    window: usize,
    filter: &'a impl Fn(&K) -> bool,
    buffer: &mut Vec<f32>,
) -> Result<(), DataProcessingError> {
    let inner_filter = |i| slices.mz_order.get(i).is_some_and(|(k, _mz)| filter(k));
    coelution_vref_score_filter_onto(&slices.arr, ref_slice, window, inner_filter, buffer)?;
    Ok(())
}

// /// See the docs for [`coelution_score_arr`].
// pub fn coelution_score<'a, const TOP_N: usize, K: Clone + Ord>(
//     slices: &'a MzMajorIntensityArray<K, f32>,
//     window: usize,
// ) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
//     let filter: Option<fn(usize) -> bool> = None;
//     coelution_score_iter_filter::<TOP_N>(&slices.arr, window, filter)
// }

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close_enough(a: &[f32], b: &[f32], tol: f32) {
        assert_eq!(a.len(), b.len());
        for (i, (aa, bb)) in a.iter().zip(b.iter()).enumerate() {
            assert!(
                (aa - bb).abs() < tol,
                "Failed at index {}; Expected {:?}, got {:?}, within left: {:?}, right: {:?}",
                i,
                aa,
                bb,
                a,
                b
            );
        }
    }

    #[test]
    fn test_coelution_vref_score_filter_into() {
        // 1. Setup the MzMajorIntensityArray
        let mz_order = vec![(1, 100.0), (2, 200.0)];
        let n_cycles = 4;
        let cycle_offset = 0;
        let mut slices =
            MzMajorIntensityArray::try_new_empty(mz_order.into(), n_cycles, cycle_offset).unwrap();
        slices.arr = Array2D::new(vec![
            [1.0, 1.0, 0.0, 0.0], // id = 1
            [0.0, 0.0, 1.0, 1.0], // id = 2
        ])
        .unwrap();

        // 2. Setup other parameters
        let ref_slice = vec![1.0, 1.0, 0.0, 0.0];
        let window = 3;
        let mut buffer = Vec::new();

        // 3. Define a filter and call the function
        let filter = |k: &i32| *k == 1 || *k == 2; // Include both
        coelution_vref_score_filter_into(&slices, &ref_slice, window, &filter, &mut buffer)
            .unwrap();

        // 4. Define expected results and assert
        // For row 1 vs ref_slice (self): rolling cos sim is roughly [0, 1, 1, 0]
        // For row 2 vs ref_slice (orthogonal): rolling cos sim is roughly [0, 0, 0, 0]
        // The sum of similarities is [0.0, 1.0, 1.0, 0.0]
        // The number of elements is 2, so the norm_factor is 0.5.
        // Expected buffer = sum * norm_factor
        let expected = vec![0.0, 0.5, 0.5, 0.0];
        assert_close_enough(&buffer, &expected, 1e-7);

        // 5. Test the filter that returns an error (only one item selected)
        let filter_one = |k: &i32| *k == 1;
        let result =
            coelution_vref_score_filter_into(&slices, &ref_slice, window, &filter_one, &mut buffer);
        assert!(result.is_err());
    }
}
