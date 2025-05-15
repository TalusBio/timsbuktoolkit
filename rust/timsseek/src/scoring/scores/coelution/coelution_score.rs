
use crate::errors::DataProcessingError;
use crate::utils::correlation::rolling_cosine_similarity;
use crate::utils::top_n_array::TopNArray;
use timsquery::models::{
    Array2D,
    MzMajorIntensityArray,
};
use tracing::debug;

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
/// assert_eq!(scores, [0.0, 0.9866667, 0.9939657, 0.9849558, 0.0]);
/// ```
pub fn coelution_score_arr<const TOP_N: usize>(
    slices: &Array2D<f32>,
    window_size: usize,
) -> Result<Vec<f32>, DataProcessingError> {
    let filter: Option<fn(usize) -> bool> = None;
    coelution_score_arr_filter::<TOP_N>(slices, window_size, &filter)
}

/// Calculates the coelution score of a set of chromatograms.
/// with a filter ...
pub fn coelution_score_arr_filter<const TOP_N: usize>(
    slices: &Array2D<f32>,
    window_size: usize,
    filter: &Option<impl Fn(usize) -> bool>,
) -> Result<Vec<f32>, DataProcessingError> {
    if slices.ncols() < window_size {
        debug!("Not enough data to calculate coelution score");
        return Err(DataProcessingError::ExpectedNonEmptyData { context: None });
    }
    let mut scores = vec![TopNArray::<TOP_N, f32>::new(); slices.ncols()];
    for i in 0..slices.nrows() {
        // I also think this is a very ugly chunk of code ...
        // that I need to find a better way to make more elegant.
        match filter {
            Some(f) => {
                if !f(i) {
                    continue;
                }
            }
            None => {}
        };
        let slice1 = slices.get_row(i).expect("Using nrows to check length");
        for j in 0..slices.nrows() {
            if j >= i {
                continue;
            }
            match filter {
                Some(f) => {
                    if !f(j) {
                        continue;
                    }
                }
                None => {}
            };
            let slice2 = slices.get_row(j).expect("Using nrows to check length");
            let cosine_similarities = rolling_cosine_similarity(slice1, slice2, window_size)
                .expect("The passed array should already be checked for length and non-emptyness");

            for (si, score) in cosine_similarities.into_iter().enumerate() {
                // Q: Should this raise an error or warn?
                // Q: SHould I weight on the expected intensity?
                if !score.is_nan() && score > 0.0 {
                    scores[si].push(score.powi(2));
                }
            }
        }
    }

    // TODO: Make this a single pass
    // In theory I could make the calculation of the rolling similarities an iterator
    // and then do a single pass over the iterators.
    let mut out = scores
        .iter()
        .map(|x| {
            let curr_vals = x.get_values();
            let out = curr_vals.iter().sum::<f32>() / x.capacity() as f32;
            assert!(
                out <= 1.0,
                "Coelution score greater than 1.0, vals: {:?}, out: {}",
                curr_vals,
                out
            );
            out
        })
        .collect::<Vec<f32>>();

    let last_ind = out.len() - 1;
    for i in 0..(window_size / 2) {
        out[i] = 0.0;
        out[last_ind - i] = 0.0;
    }

    Ok(out)
}

/// See the docs for [`coelution_score_arr`].
pub fn coelution_score_filter<const TOP_N: usize, K: Clone + Ord>(
    slices: &MzMajorIntensityArray<K, f32>,
    window: usize,
    filter: &Option<impl Fn(&K) -> bool>,
) -> Result<Vec<f32>, DataProcessingError> {
    let inner_filter = match filter {
        Some(f) => Some(|i| slices.mz_order.get(i).map_or(false, |(k, _mz)| f(k))),
        None => None,
    };
    coelution_score_arr_filter::<TOP_N>(&slices.arr, window, &inner_filter)
}

/// See the docs for [`coelution_score_arr`].
pub fn coelution_score<const TOP_N: usize, K: Clone + Ord>(
    slices: &MzMajorIntensityArray<K, f32>,
    window: usize,
) -> Result<Vec<f32>, DataProcessingError> {
    let filter: Option<fn(usize) -> bool> = None;
    coelution_score_arr_filter::<TOP_N>(&slices.arr, window, &filter)
}

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
        let scores = coelution_score_arr::<1>(&slices, window);
        let expected = vec![0.0, 0.75, 0.974026, 0.99997497, 0.99995005, 0.0];
        assert_close_enough(&scores.unwrap(), &expected, 1e-2);
    }

    #[test]
    fn test_coelution_score_eq() {
        let arr = vec![0., 1., 1., 3., 200., 5.];
        let slices = Array2D::new(vec![arr.clone(), arr.clone()]).unwrap();
        let window = 3;
        let scores = coelution_score_arr::<2>(&slices, window);
        let expected = vec![0.0, 0.5, 0.5, 0.5, 0.5, 0.0];
        assert_close_enough(&scores.unwrap(), &expected, 1e-7);
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
