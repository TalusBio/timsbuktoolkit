use std::f32;

use crate::errors::DataProcessingError;

// TODO: benchmark how much faster is f32

/// Calculates the cosine similarity between two vectors of the same size.
///
/// # Example
///
/// ```
/// use timsseek::utils::correlation::cosine_similarity;
///
/// let a = vec![1.0, 2.0, 3.0];
/// let b = vec![4.0, 5.0, 6.0];
/// let result = cosine_similarity(&a, &b).unwrap();
/// assert_eq!(result, 0.9746318);
/// ```
pub fn cosine_similarity(a: &[f32], b: &[f32]) -> Result<f32, DataProcessingError> {
    // Check if vectors have the same length and are not empty
    if a.len() != b.len() || a.is_empty() {
        return Err(DataProcessingError::ExpectedSlicesSameLength {
            expected: a.len(),
            other: b.len(),
            context: "cosine_similarity".to_string(),
        });
    }

    // Calculate dot product (numerator)
    let dot_product: f32 = a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum();

    // Calculate magnitudes (denominator)
    let magnitude_a: f32 = a.iter().map(|&x| x * x).sum::<f32>().sqrt();

    let magnitude_b: f32 = b.iter().map(|&x| x * x).sum::<f32>().sqrt();

    // Avoid division by zero
    if magnitude_a == 0.0 || magnitude_b == 0.0 {
        return Ok(f32::NAN);
    }

    // Calculate cosine similarity
    Ok(dot_product / (magnitude_a * magnitude_b))
}

#[derive(Debug, Copy, Clone)]
struct RollingElem {
    a_sq: u64,
    b_sq: u64,
    prod: u64,
}

const MAX_CAPACITY: usize = 20;

#[derive(Debug)]
struct CosineSimilarityCircularBuffer {
    a_sum_sq: u64,
    b_sum_sq: u64,
    dot_product: u64,
    buffer: [RollingElem; MAX_CAPACITY],
    size: usize,
    curr_index: usize,
}

impl CosineSimilarityCircularBuffer {
    fn new(a: &[f32], b: &[f32]) -> Self {
        if a.len() != b.len() {
            // I am panicking here bc I dont really want this to be
            // a recoverable error, if different lengths are passed
            // the program is already degenerate enough that it should
            // panic
            panic!("Vectors must be of the same length");
        }
        if a.len() > MAX_CAPACITY {
            // I am panicking here bc I dont really want this to be
            // a recoverable error, if different lengths are passed
            // the program is already degenerate enough that it should
            // panic
            panic!("Vectors must be of length <= {}", MAX_CAPACITY);
        }
        let buffer = [RollingElem {
            a_sq: 0,
            b_sq: 0,
            prod: 0,
        }; MAX_CAPACITY];

        let mut out = Self {
            a_sum_sq: 0,
            b_sum_sq: 0,
            dot_product: 0,
            buffer,
            size: a.len(),
            curr_index: 0,
        };
        for (&a, &b) in a.iter().zip(b.iter()) {
            out.update(a, b);
        }
        // Make sure the current index is at 0
        assert!(out.curr_index == 0);
        assert!(out.size == a.len());

        out
    }

    fn calculate_similarity(&mut self) -> f32 {
        let mag_a = (self.a_sum_sq as f64).sqrt();
        let mag_b = (self.b_sum_sq as f64).sqrt();
        let denom = mag_a * mag_b;
        if denom == 0.0 {
            return 0.0;
        }

        let out = self.dot_product as f32 / denom as f32;
        assert!(out >= 0.0);
        assert!(
            out <= 1.0001,
            "Expected <=1 got {} (1.0 + 1e-4 is acceptable)",
            out
        );

        out
    }

    fn update(&mut self, a: f32, b: f32) {
        let a: u64 = a as u64;
        let b: u64 = b as u64;
        let prod = a * b;

        let curr_elem = RollingElem {
            a_sq: a * a,
            b_sq: b * b,
            prod,
        };

        self.a_sum_sq -= self.buffer[self.curr_index].a_sq;
        self.a_sum_sq += curr_elem.a_sq;
        self.b_sum_sq -= self.buffer[self.curr_index].b_sq;
        self.b_sum_sq += curr_elem.b_sq;
        self.dot_product -= self.buffer[self.curr_index].prod;
        self.dot_product += curr_elem.prod;

        // Recalculate every 5x size ... this helps with numerical errors
        // Ideally we would never recalculate but I have not been able to
        // manage the stability of the sums.
        self.buffer[self.curr_index] = curr_elem;
        self.curr_index = (self.curr_index + 1) % self.size;
    }
}

/// Calculates the cosine similarity at every window for two vectors.
///
/// This means that for every window of size `window_size` the cosine similarity
/// is calculated between the two vectors. (and is padded with f32::NAN)
/// The output is a vector of the same length as the input.
///
/// # Example
///
/// ```
/// use timsseek::utils::correlation::rolling_cosine_similarity;
///
/// let a = vec![1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
/// let b = vec![4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0];
/// let expect_res: [f32; 9] = [
///     f32::NAN,
///     0.9746319, 0.9746319, 0.9746319, 0.9746319, 0.9746319, 0.9746319, 0.9746319,
///     f32::NAN,
/// ];
/// let results = rolling_cosine_similarity(&a, &b, 3).unwrap();
/// assert_eq!(results.len(), expect_res.len());
/// for (result, expect) in results.iter().zip(expect_res.iter()) {
///     if result.is_nan() {
///         assert!(expect.is_nan());
///     } else {
///         assert!((result - expect).abs() < 1e-2, "{:?} vs {:?} at {:?}", results, expect_res, result);
///     }
/// }
/// ```
pub fn rolling_cosine_similarity(
    a: &[f32],
    b: &[f32],
    window_size: usize,
) -> Result<Vec<f32>, DataProcessingError> {
    let mut results = vec![f32::NAN; a.len()];
    rolling_cosine_similarity_into(a, b, window_size, &mut results)?;
    Ok(results)
}

pub fn rolling_cosine_similarity_into(
    a: &[f32],
    b: &[f32],
    window_size: usize,
    result_vec: &mut Vec<f32>,
) -> Result<(), DataProcessingError> {
    // let mut results = vec![f32::NAN; a.len()];
    result_vec.clear();
    result_vec.resize(a.len(), f32::NAN);
    // Check if vectors have the same length and are long enough for the window
    if a.len() != b.len() {
        return Err(DataProcessingError::ExpectedSlicesSameLength {
            expected: a.len(),
            other: b.len(),
            context: "cosine_similarity".to_string(),
        });
    }
    if a.len() < window_size {
        return Err(DataProcessingError::ExpectedNonEmptyData {
            context: Some("cosine_similarity".to_string()),
        });
    }

    let offset = window_size / 2;

    // Initialize the first window
    let mut cosine_sim =
        CosineSimilarityCircularBuffer::new(&a[0..window_size], &b[0..window_size]);

    // Calculate similarity for first window
    result_vec[offset] = cosine_sim.calculate_similarity();

    // Roll the window
    // for i in (offset + 1)..(a.len() - offset) {
    for i in 0..(a.len() - window_size) {
        cosine_sim.update(a[i + window_size], b[i + window_size]);
        result_vec[i + offset + 1] = cosine_sim.calculate_similarity();
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cosine_similarity() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let result = cosine_similarity(&a, &b).unwrap();
        assert!((result - 0.974_631_85).abs() < 1e-4);
    }

    #[test]
    fn test_identical_vectors() {
        let a = vec![1.0, 1.0, 1.0];
        let result = cosine_similarity(&a, &a).unwrap();
        assert!((result - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_empty_vectors() {
        let a: Vec<f32> = vec![];
        let b: Vec<f32> = vec![];
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
        let a: Vec<f32> = vec![1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
        let b: Vec<f32> = vec![4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0];
        let expect_res: [f32; 9] = [
            f32::NAN,
            0.97463185,
            0.97463185,
            0.97463185,
            0.97463185,
            0.97463185,
            0.97463185,
            0.97463185,
            f32::NAN,
        ];
        let results = rolling_cosine_similarity(&a, &b, 3).unwrap();
        assert_eq!(results.len(), expect_res.len());

        for (result, expect) in results.iter().zip(expect_res.iter()) {
            if result.is_nan() {
                assert!(expect.is_nan());
            } else {
                assert!(
                    (result - expect).abs() < 1e-2,
                    "Expected {:?} in {:?}, got {:?} in {:?}",
                    expect,
                    expect_res,
                    result,
                    results,
                );
            }
        }
    }

    #[test]
    fn test_rolling_basic() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.0, 2.0, 3.0, 4.0];
        // Weird things happen when even number windows are used.
        let expect_res: [f32; 4] = [f32::NAN, 1.0, 1.0, 1.0];
        let results = rolling_cosine_similarity(&a, &b, 2).unwrap();
        assert_eq!(results.len(), 4);
        for result in results.iter().zip(expect_res.iter()) {
            if result.0.is_nan() {
                assert!(result.1.is_nan());
            } else {
                assert!(
                    (result.0 - result.1).abs() < 1e-4,
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
        let results = rolling_cosine_similarity(&a, &b, 3).unwrap();
        assert_eq!(results.len(), 4);
        assert!(results[0].is_nan());
        assert!(results[3].is_nan());
        assert_eq!(results[1], 1.0, "{:?}", results);
        assert_eq!(results[2], 1.0, "{:?}", results);
    }

    #[test]
    fn test_rolling_zeros_window() {
        let a = vec![0.0; 20];
        let b = vec![0.0; 20];
        let results = rolling_cosine_similarity(&a, &b, 5).unwrap();
        assert_eq!(results.len(), 20);
        let expect = [
            f32::NAN,
            f32::NAN,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            f32::NAN,
            f32::NAN,
        ];
        assert_eq!(results.len(), expect.len());
        for (result, expect) in results.iter().zip(expect.iter()) {
            if result.is_nan() {
                assert!(expect.is_nan());
            } else {
                assert!((result - expect).abs() < 1e-4, "{:?}", results);
            }
        }
    }
}
