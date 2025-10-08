use crate::errors::DataProcessingError;
use crate::utils::correlation::rolling_cosine_similarity;
use timsquery::models::{Array2D, MzMajorIntensityArray};
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
    (0..slices.nrows())
        .filter(|&i| filter(i))
        .try_for_each(|i| {
            let slice1 = slices.get_row(i).expect("Row index i is within bounds");
            let iter = rolling_cosine_similarity(slice1, ref_slice, window_size)?;
            for (i, v) in iter.enumerate() {
                if v.is_nan() {
                    continue;
                }
                buffer[i] += v.max(0.0) * norm_factor;
            }
            Ok(())
        })
}

/// See the docs for [`coelution_score_arr`].
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

    // #[test]
    // fn test_coelution_score() {
    //     let slices =
    //         Array2D::new(vec![[0., 1., 1., 3., 200., 5.], [1., 2., 1., 4., 200., 6.]]).unwrap();
    //     let window = 3;
    //     let scores = coelution_score_arr::<1>(&slices, window)
    //         .unwrap()
    //         .collect::<Vec<_>>();
    //     let expected = vec![0.0, 0.75, 0.974026, 0.99997497, 0.99995005, 0.0];
    //     assert_close_enough(&scores, &expected, 1e-2);
    // }

    // #[test]
    // fn test_coelution_score_eq() {
    //     let arr = vec![0., 1., 1., 3., 200., 5.];
    //     let slices = Array2D::new(vec![arr.clone(), arr.clone()]).unwrap();
    //     let window = 3;
    //     let scores = coelution_score_arr::<2>(&slices, window)
    //         .unwrap()
    //         .collect::<Vec<_>>();
    //     let expected = vec![0.0, 0.5, 0.5, 0.5, 0.5, 0.0];
    //     assert_close_enough(&scores, &expected, 1e-7);
    // }

    // #[test]
    // fn test_coelution_score_no_overlap() {
    //     let slices = Array2D::new(vec![[0., 1.], [1., 2.], [2., 3.]]).unwrap();
    //     // Since the window is larger than the size of the slices, we should get an error.
    //     let window = 3;
    //     let scores = coelution_score_arr::<2>(&slices, window);
    //     assert!(scores.is_err());
    // }
}
