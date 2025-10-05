use crate::errors::DataProcessingError;
use crate::utils::correlation::rolling_cosine_similarity;
use crate::utils::top_n_array::TopNArray;
use timsquery::models::{
    Array2D,
    MzMajorIntensityArray,
};
use tracing::trace;

/// Calculates the coelution score of a set of chromatograms.
///
/// # Example
///
/// ```
/// use timsquery::Array2D;
/// use timsseek::scoring::coelution;
///
/// let slices = Array2D::new(
///     vec![[0., 1., 3., 22., 5.],
///          [0., 2., 4., 20., 5.],
///          [0., 1., 2., 19., 2.]]).unwrap();
/// let window = 3;
/// // Note that the generic type parameter is the top N of scores that will
/// // be averaged to report the coelution.
/// let scores = coelution::coelution_score_arr::<3>(&slices, window).unwrap();
/// assert_eq!(scores, [0.0, 0.9866667, 0.9939658, 0.9849558, 0.0]);
/// ```
pub fn coelution_score_arr<'a, const TOP_N: usize>(
    slices: &'a Array2D<f32>,
    window_size: usize,
) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
    let filter: Option<fn(usize) -> bool> = None;
    coelution_score_iter_filter::<TOP_N>(slices, window_size, filter)
}

// Assuming these types are in the current scope:
// use crate::{Array2D, DataProcessingError, TopNArray, rolling_cosine_similarity};

struct TopMultiIter<I: Iterator<Item = f32>, const TOP_N: usize> {
    iters: Vec<I>,
}

impl<I: Iterator<Item = f32>, const TOP_N: usize> Iterator for TopMultiIter<I, TOP_N> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        let mut top_n = TopNArray::<TOP_N, f32>::new();
        let mut all_none = true;
        let mut any_none = false;
        for iter in &mut self.iters {
            if let Some(score) = iter.next() {
                all_none = false;
                if !score.is_nan() && score > 0.0 {
                    top_n.push(score.powi(2));
                }
            } else {
                any_none = true;
            }
        }
        if all_none {
            return None;
        }
        if any_none {
            panic!("Iterators should all have the same length");
        }

        // Calculate the final score for the current time point.
        let curr_vals = top_n.get_values();
        let out = curr_vals.iter().sum::<f32>() / top_n.capacity() as f32;
        // let out = curr_vals.iter().sum::<f32>();

        // This is required bc the main score assumes scores are between 0 and 1
        // For now ... and catching it earlier is good! (instead of later when the
        // main score is getting calculated)
        assert!(
            out <= 1.0,
            "Coelution score greater than 1.0, vals: {:?}, out: {}",
            curr_vals,
            out
        );
        Some(out)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.iters.is_empty() {
            return (0, Some(0));
        }
        let (mut min, mut max) = self.iters[0].size_hint();
        for it in &self.iters[1..] {
            let (it_min, it_max) = it.size_hint();
            min = min.min(it_min);
            max = match (max, it_max) {
                (Some(a), Some(b)) => Some(a.min(b)),
                _ => None,
            };
        }
        (min, max)
    }
}

// This alternative calculate K^2 pairwise comparisons, which is more expensive
// and doesn't seem to improve the score in practice.
//
// /// Calculates the coelution score of a set of chromatograms, returning a lazy iterator.
// fn coelution_score_iter_filter<'a, const TOP_N: usize>(
//     slices: &'a Array2D<f32>,
//     window_size: usize,
//     filter: Option<impl Fn(usize) -> bool>,
// ) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
//     if slices.ncols() < window_size {
//         trace!("Not enough data to calculate coelution score");
//         return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
//     }
// 
//     // Collect all pairwise rolling similarity calculations.
//     // This still computes the similarities upfront, but the aggregation into the final
//     // score is done lazily, one time-point at a time.
//     let similarity_iters: Vec<_> = (0..slices.nrows())
//         .filter(|&i| filter.as_ref().map_or(true, |f| f(i)))
//         .tuple_combinations()
//         .map(|(i, j)| {
//             let slice1 = slices.get_row(i).expect("Row index i is within bounds");
//             let slice2 = slices.get_row(j).expect("Row index j is within bounds");
//             rolling_cosine_similarity(slice1, slice2, window_size)
//         })
//         .collect::<Result<Vec<_>, _>>()? // Handle potential errors from `rolling_cosine_similarity`
//         .into_iter()
//         .map(|v| {
//             v.into_iter()
//                 .map(|x| if x.is_nan() { 0.0 } else { x.max(0.0) })
//         })
//         .collect();
// 
//     let out = TopMultiIter::<_, TOP_N> {
//         iters: similarity_iters,
//     };
// 
//     Ok(out)
// }

/// Calculates the coelution score of a set of chromatograms, returning a lazy iterator.
fn coelution_vref_score_iter_filter<'a, const TOP_N: usize>(
    slices: &'a Array2D<f32>,
    ref_slice: &'a [f32],
    window_size: usize,
    filter: Option<impl Fn(usize) -> bool>,
) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
    if slices.ncols() < window_size {
        trace!("Not enough data to calculate coelution score");
        return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
    }

    // Collect all rolling similarity calculations v the reference.
    // This still computes the similarities upfront, but the aggregation into the final
    // score is done lazily, one time-point at a time.
    let similarity_iters: Vec<_> = (0..slices.nrows())
        .filter(|&i| filter.as_ref().map_or(true, |f| f(i)))
        .map(|i| {
            let slice1 = slices.get_row(i).expect("Row index i is within bounds");
            rolling_cosine_similarity(slice1, ref_slice, window_size)
        })
        .collect::<Result<Vec<_>, _>>()? // Handle potential errors from `rolling_cosine_similarity`
        .into_iter()
        .map(|v| {
            v.into_iter()
                .map(|x| if x.is_nan() { 0.0 } else { x.max(0.0) })
        })
        .collect();

    assert!(
        similarity_iters.len() > 0,
        "No valid slices after filtering"
    );
    assert!(
        similarity_iters.len() < 50,
        "There are too many valid slices after filtering, probably an mz-major and an rt-major array got mixed up"
    );

    let out = TopMultiIter::<_, TOP_N> {
        iters: similarity_iters,
    };

    Ok(out)
}

/// See the docs for [`coelution_score_arr`].
pub fn coelution_score_filter<'a, const TOP_N: usize, K: Clone + Ord>(
    slices: &'a MzMajorIntensityArray<K, f32>,
    window: usize,
    filter: &'a Option<impl Fn(&K) -> bool>,
) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
    let inner_filter = filter
        .as_ref()
        .map(|f| |i| slices.mz_order.get(i).is_some_and(|(k, _mz)| f(k)));
    coelution_score_iter_filter::<TOP_N>(&slices.arr, window, inner_filter)
}

/// See the docs for [`coelution_score_arr`].
pub fn coelution_vref_score_filter<'a, const TOP_N: usize, K: Clone + Ord>(
    slices: &'a MzMajorIntensityArray<K, f32>,
    ref_slice: &'a [f32],
    window: usize,
    filter: &'a Option<impl Fn(&K) -> bool>,
) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
    let inner_filter = filter
        .as_ref()
        .map(|f| |i| slices.mz_order.get(i).is_some_and(|(k, _mz)| f(k)));
    coelution_vref_score_iter_filter::<TOP_N>(&slices.arr, ref_slice, window, inner_filter)
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
    fn test_coelution_score() {
        let slices =
            Array2D::new(vec![[0., 1., 1., 3., 200., 5.], [1., 2., 1., 4., 200., 6.]]).unwrap();
        let window = 3;
        let scores = coelution_score_arr::<1>(&slices, window)
            .unwrap()
            .collect::<Vec<_>>();
        let expected = vec![0.0, 0.75, 0.974026, 0.99997497, 0.99995005, 0.0];
        assert_close_enough(&scores, &expected, 1e-2);
    }

    #[test]
    fn test_coelution_score_eq() {
        let arr = vec![0., 1., 1., 3., 200., 5.];
        let slices = Array2D::new(vec![arr.clone(), arr.clone()]).unwrap();
        let window = 3;
        let scores = coelution_score_arr::<2>(&slices, window)
            .unwrap()
            .collect::<Vec<_>>();
        let expected = vec![0.0, 0.5, 0.5, 0.5, 0.5, 0.0];
        assert_close_enough(&scores, &expected, 1e-7);
    }

    #[test]
    fn test_coelution_score_no_overlap() {
        let slices = Array2D::new(vec![[0., 1.], [1., 2.], [2., 3.]]).unwrap();
        // Since the window is larger than the size of the slices, we should get an error.
        let window = 3;
        let scores = coelution_score_arr::<2>(&slices, window);
        assert!(scores.is_err());
    }
}
