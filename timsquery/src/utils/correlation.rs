use std::f64;

use crate::errors::{
    DataProcessingError,
    Result,
};

/// Calculates the cosine similarity between two vectors of the same size.
///
/// # Example
///
/// ```
/// use timsquery::utils::correlation::cosine_similarity;
///
/// let a = vec![1.0, 2.0, 3.0];
/// let b = vec![4.0, 5.0, 6.0];
/// let result = cosine_similarity(&a, &b).unwrap();
/// assert_eq!(result, 0.9746318461970762);
/// ```
pub fn cosine_similarity(a: &[f64], b: &[f64]) -> Result<f64> {
    // Check if vectors have the same length and are not empty
    if a.len() != b.len() || a.is_empty() {
        return Err(DataProcessingError::ExpectedVectorSameLength(a.len(), b.len()).into());
    }

    // Calculate dot product (numerator)
    let dot_product: f64 = a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum();

    // Calculate magnitudes (denominator)
    let magnitude_a: f64 = a.iter().map(|&x| x * x).sum::<f64>().sqrt();

    let magnitude_b: f64 = b.iter().map(|&x| x * x).sum::<f64>().sqrt();

    // Avoid division by zero
    if magnitude_a == 0.0 || magnitude_b == 0.0 {
        return Ok(f64::NAN);
    }

    // Calculate cosine similarity
    Ok(dot_product / (magnitude_a * magnitude_b))
}

/// Calculates the cosine similarity at every window for two vectors.
///
/// This means that for every window of size `window_size` the cosine similarity
/// is calculated between the two vectors. (and is padded with f64::NAN)
/// The output is a vector of the same length as the input.
///
/// # Example
///
/// ```
/// use timsquery::utils::correlation::rolling_cosine_similarity;
///
/// let a = vec![1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
/// let b = vec![4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0];
/// let expect_res: [f64; 9] = [
///     f64::NAN,
///     0.97463184, // Note that this is the same as above
///     0.97823,
///     0.95065,
///     0.97463184, // And this one
///     0.97823,
///     0.95065,
///     0.97463184,
///     f64::NAN,
/// ];
/// let results = rolling_cosine_similarity(&a, &b, 3).unwrap();
/// assert_eq!(results.len(), expect_res.len());
/// ```
pub fn rolling_cosine_similarity(a: &[f64], b: &[f64], window_size: usize) -> Result<Vec<f64>> {
    // Check if vectors have the same length and are long enough for the window
    if a.len() != b.len() {
        return Err(DataProcessingError::ExpectedVectorSameLength(a.len(), b.len()).into());
    }
    if a.len() < window_size {
        return Err(DataProcessingError::InsufficientData {
            real: a.len(),
            expected: window_size,
        }
        .into());
    }

    let offset = window_size / 2;
    let mut results = vec![f64::NAN; a.len()];

    // Initialize the first window
    let mut dot_product = 0.0;
    let mut sum_a_squared = 0.0;
    let mut sum_b_squared = 0.0;

    // Calculate initial values for first window
    for i in 0..window_size {
        dot_product += a[i] * b[i];
        sum_a_squared += a[i] * a[i];
        sum_b_squared += b[i] * b[i];
    }

    // Calculate similarity for first window
    results[offset] = calculate_similarity(dot_product, sum_a_squared, sum_b_squared);

    // Roll the window
    for i in (offset + 1)..(a.len() - offset) {
        // Remove contribution of element leaving the window
        let leaving_idx = i - offset;
        dot_product -= a[leaving_idx] * b[leaving_idx];
        sum_a_squared -= a[leaving_idx] * a[leaving_idx];
        sum_b_squared -= b[leaving_idx] * b[leaving_idx];

        // Add contribution of new element
        dot_product += a[i] * b[i];
        sum_a_squared += a[i] * a[i];
        sum_b_squared += b[i] * b[i];

        // Calculate similarity for current window
        results[i] = calculate_similarity(dot_product, sum_a_squared, sum_b_squared);
    }

    Ok(results)
}

#[inline(always)]
fn calculate_similarity(dot_product: f64, sum_a_squared: f64, sum_b_squared: f64) -> f64 {
    let magnitude_a = sum_a_squared.sqrt();
    let magnitude_b = sum_b_squared.sqrt();

    if magnitude_a == 0.0 || magnitude_b == 0.0 {
        0.0
    } else {
        dot_product / (magnitude_a * magnitude_b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cosine_similarity() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let result = cosine_similarity(&a, &b).unwrap();
        assert!((result - 0.974631846).abs() < 1e-8);
    }

    #[test]
    fn test_identical_vectors() {
        let a = vec![1.0, 1.0, 1.0];
        let result = cosine_similarity(&a, &a).unwrap();
        assert!((result - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_empty_vectors() {
        let a: Vec<f64> = vec![];
        let b: Vec<f64> = vec![];
        assert!(cosine_similarity(&a, &b).is_err());
    }

    #[test]
    fn test_different_lengths() {
        let a = vec![1.0, 2.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(cosine_similarity(&a, &b).is_err());
    }

    #[test]
    fn test_zero_vector() {
        let a = vec![0.0, 0.0, 0.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(cosine_similarity(&a, &b).unwrap().is_nan());
    }

    #[test]
    fn test_rolling_cosine_similarity() {
        let a = vec![1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0];
        let expect_res: [f64; 9] = [
            f64::NAN,
            0.97463184, // Note that this is the same as above
            0.97823,
            0.95065,
            0.97463184, // And this one
            0.97823,
            0.95065,
            0.97463184,
            f64::NAN,
        ];
        let results = rolling_cosine_similarity(&a, &b, 3).unwrap();
        assert_eq!(results.len(), expect_res.len());

        for (result, expect) in results.iter().zip(expect_res.iter()) {
            if result.is_nan() {
                assert!(expect.is_nan());
            } else {
                assert!(
                    (result - expect).abs() < 1e-3,
                    "Expected {:?}, got {:?}",
                    expect_res,
                    results,
                );
            }
        }
    }

    #[test]
    fn test_rolling_basic() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let expect_res: [f64; 4] = [f64::NAN, 1.0, 1.0, f64::NAN];
        let results = rolling_cosine_similarity(&a, &b, 2).unwrap();
        assert_eq!(results.len(), 4);
        for result in results.iter().zip(expect_res.iter()) {
            if result.0.is_nan() {
                assert!(result.1.is_nan());
            } else {
                assert!(
                    (result.0 - result.1).abs() < 1e-8,
                    "Expected {:?}, got {:?}",
                    results,
                    expect_res
                );
            }
        }
    }

    #[test]
    fn test_rolling_window_too_large() {
        let a = vec![1.0, 2.0];
        let b = vec![1.0, 2.0];
        let results = rolling_cosine_similarity(&a, &b, 3);
        assert!(results.is_err());
    }

    #[test]
    fn test_rolling_different_lengths() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0];
        let results = rolling_cosine_similarity(&a, &b, 2);
        assert!(results.is_err());
    }

    #[test]
    fn test_rolling_zero_window() {
        let a = vec![0.0, 0.0, 1.0, 1.0];
        let b = vec![0.0, 0.0, 1.0, 1.0];
        let results = rolling_cosine_similarity(&a, &b, 2).unwrap();
        assert_eq!(results.len(), 4);
        assert!(results[0].is_nan());
        assert!(results[3].is_nan());
        assert_eq!(results[1], 0.0, "{:?}", results);
        assert_eq!(results[2], 1.0, "{:?}", results);
    }
}
