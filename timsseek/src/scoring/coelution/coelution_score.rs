use crate::models::{
    Array2D,
    MzMajorIntensityArray,
};
use crate::utils::correlation::rolling_cosine_similarity;
use crate::utils::top_n_array::TopNArray;
use tracing::warn;

/// Calculates the coelution score of a set of chromatograms.
///
/// # Example
///
/// ```
/// use timsseek::models::Array2D;
/// use timsseek::scoring::coelution::coelution_score;
///
/// let slices = Array2D::new(
///     vec![[0., 1., 3., 22., 5.],
///          [0., 2., 4., 20., 5.],
///          [0., 1., 2., 19., 2.]]).unwrap();
/// let window = 3;
/// let scores = coelution_score::coelution_score_arr(&slices, window);
/// assert_eq!(scores, [0.0, 0.9932996624407776, 1.0, 0.9977801004359416, 0.0]);
/// ```
fn coelution_score_arr(slices: &Array2D<f32>, window_size: usize) -> Vec<f32> {
    const TOP_N: usize = 6;
    if slices.ncols() < window_size {
        warn!("Not enough data to calculate coelution score");
        return vec![0.0; slices.ncols()];
    }
    let mut scores = vec![TopNArray::<TOP_N, f32>::new(); slices.ncols()];
    for i in 0..slices.nrows() {
        let slice1 = slices.get_row(i);
        for j in 0..slices.nrows() {
            if j >= i {
                continue;
            }
            let slice2 = slices.get_row(j);
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
            assert!(out <= 1.0, "Coelution score greater than 1.0, vals: {:?}, out: {}", curr_vals, out);
            out
        })
        .collect::<Vec<f32>>();

    let last_ind = out.len() - 1;
    for i in 0..(window_size / 2) {
        out[i] = 0.0;
        out[last_ind - i] = 0.0;
    }

    out
}

/// See the docs for [`coelution_score_arr`].
pub fn coelution_score(slices: &MzMajorIntensityArray, window: usize) -> Vec<f32> {
    coelution_score_arr(&slices.arr, window)
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
        let scores = coelution_score_arr(&slices, window);
        let expected = vec![0.0, 0.866, 0.816, 0.968, 0.999, 0.0];
        assert_close_enough(&scores, &expected, 1e-3);
    }

    #[test]
    fn test_coelution_score_eq() {
        let arr = vec![0., 1., 1., 3., 200., 5.];
        let slices = Array2D::new(vec![arr.clone(), arr.clone()]).unwrap();
        let window = 3;
        let scores = coelution_score_arr(&slices, window);
        let expected = vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0];
        assert_close_enough(&scores, &expected, 1e-8);
    }

    #[test]
    fn test_coelution_score_no_overlap() {
        let slices = Array2D::new(vec![[0., 1.], [1., 2.], [2., 3.]]).unwrap();
        let window = 3;
        let scores = coelution_score_arr(&slices, window);
        let expected = vec![0.0, 0.0];
        assert_close_enough(&scores, &expected, 1e-8);
    }
}
