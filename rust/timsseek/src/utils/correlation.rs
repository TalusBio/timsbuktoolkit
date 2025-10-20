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
    window_size: usize,
    curr_index: usize,
}

impl CosineSimilarityCircularBuffer {
    fn new(window_size: usize) -> Self {
        if window_size > MAX_CAPACITY {
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

        Self {
            a_sum_sq: 0,
            b_sum_sq: 0,
            dot_product: 0,
            buffer,
            window_size,
            curr_index: 0,
        }
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

    // This is too hot, cannot instrument for performance reasons
    // #[cfg_attr(
    //     feature = "instrumentation",
    //     tracing::instrument(skip_all, level = "trace")
    // )]
    fn update(&mut self, a: f32, b: f32) {
        let a: u64 = a as u64;
        let b: u64 = b as u64;
        let prod = a * b;

        // Get mutable reference to the current element in the buffer
        let elem = &mut self.buffer[self.curr_index];

        // Subtract the old values from the sums
        self.a_sum_sq -= elem.a_sq;
        self.b_sum_sq -= elem.b_sq;
        self.dot_product -= elem.prod;

        // Update the element in place
        elem.a_sq = a * a;
        elem.b_sq = b * b;
        elem.prod = prod;

        // Add the new values to the sums
        self.a_sum_sq += elem.a_sq;
        self.b_sum_sq += elem.b_sq;
        self.dot_product += elem.prod;

        // Move to the next index
        self.curr_index = (self.curr_index + 1) % self.window_size;
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
/// let results = rolling_cosine_similarity(&a, &b, 3).unwrap().collect::<Vec<_>>();
/// assert_eq!(results.len(), expect_res.len());
/// for (result, expect) in results.iter().zip(expect_res.iter()) {
///     if result.is_nan() {
///         assert!(expect.is_nan());
///     } else {
///         assert!((result - expect).abs() < 1e-2, "{:?} vs {:?} at {:?}", results, expect_res, result);
///     }
/// }
/// ```
pub fn rolling_cosine_similarity<'a>(
    a: &'a [f32],
    b: &'a [f32],
    window_size: usize,
) -> Result<impl Iterator<Item = f32> + 'a, DataProcessingError> {
    // As of RN ... Oct 5 2025, this is the hottest function in the codebase.

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

    // Initialize the first window
    let cosine_sim = CosineSimilarityCircularBuffer::new(window_size);
    Ok(RollingCosineIterator {
        a,
        b,
        buffer: cosine_sim,
        window_size,
        state: 0,
    })

    // // This alternative is only slightly slower ...
    // let left_padding = window_size / 2;
    // let right_padding = window_size - left_padding - 1;
    // let tmp = a
    //     .windows(window_size)
    //     .zip(b.windows(window_size))
    //     .map(|(a_win, b_win)| cosine_similarity(a_win, b_win).unwrap());
    // // Pad the beginning and end with NaNs
    // let iter = std::iter::repeat(f32::NAN)
    //     .take(left_padding)
    //     .chain(tmp)
    //     .chain(std::iter::repeat(f32::NAN).take(right_padding));
    // Ok(iter)
}

struct RollingCosineIterator<'a> {
    a: &'a [f32],
    b: &'a [f32],
    buffer: CosineSimilarityCircularBuffer,
    window_size: usize,
    state: usize,
}

impl Iterator for RollingCosineIterator<'_> {
    type Item = f32;

    // This function is too hot, cannot instrument for performance reasons
    // #[cfg_attr(
    //     feature = "instrumentation",
    //     tracing::instrument(skip_all, level = "trace")
    // )]
    fn next(&mut self) -> Option<Self::Item> {
        if self.state >= self.a.len() {
            return None;
        }
        self.buffer.update(self.a[self.state], self.b[self.state]);
        self.state += 1;
        // So if the window size is 3, we want to pad the first and last elements
        let left_padding = self.window_size / 2;
        let right_padding = self.window_size - left_padding - 1;
        let in_init_padding = self.state <= left_padding;
        let in_end_padding = self.state > (self.a.len() - right_padding);
        if in_init_padding || in_end_padding {
            Some(f32::NAN)
        } else {
            Some(self.buffer.calculate_similarity())
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.a.len() - self.state;
        (len, Some(len))
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
        let results: Vec<_> = rolling_cosine_similarity(&a, &b, 3).unwrap().collect();
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
        let expect_res: [f32; 4] = [f32::NAN, 1.0, 1.0, f32::NAN];
        // Weird things happen when even number windows are used.
        let results: Vec<_> = rolling_cosine_similarity(&a, &b, 3).unwrap().collect();
        assert_eq!(results.len(), 4);
        for result in results.iter().zip(expect_res.iter()) {
            if result.0.is_nan() {
                assert!(
                    result.1.is_nan(),
                    "Nan at position within {:?} expected {:?}",
                    results,
                    expect_res,
                );
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
        let results: Vec<_> = rolling_cosine_similarity(&a, &b, 3).unwrap().collect();
        assert_eq!(results.len(), 4);
        assert!(results[0].is_nan());
        assert!(results[3].is_nan());
        assert_eq!(results[1], 0.0, "{:?}", results);
        assert_eq!(results[2], 1.0, "{:?}", results);
    }

    #[test]
    fn test_rolling_zeros_window() {
        let a = vec![0.0; 20];
        let b = vec![0.0; 20];
        let results: Vec<_> = rolling_cosine_similarity(&a, &b, 5).unwrap().collect();
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
