use std::f32;

use crate::errors::{
    DataProcessingError,
    Result,
};

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
pub fn cosine_similarity(a: &[f32], b: &[f32]) -> Result<f32> {
    // Check if vectors have the same length and are not empty
    if a.len() != b.len() || a.is_empty() {
        return Err(DataProcessingError::ExpectedSlicesSameLength {
            expected: a.len(),
            other: b.len(),
            context: "cosine_similarity".to_string(),
        }
        .into());
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
    a_sq: f64,
    b_sq: f64,
    prod: f64,
}

const MAX_CAPACITY: usize = 20;

#[derive(Debug)]
struct CosineSimilarityCircularBuffer {
    a_sum_sq: f64,
    b_sum_sq: f64,
    dot_product: f64,
    buffer: [RollingElem; MAX_CAPACITY],
    size: usize,
    curr_index: usize,
    abs_index: usize,
    recalculations: usize,
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
        let mut buffer = [RollingElem {
            a_sq: 0.0,
            b_sq: 0.0,
            prod: 0.0,
        }; MAX_CAPACITY];
        a.iter()
            .zip(b.iter())
            .map(|(&a, &b)| RollingElem {
                a_sq: (a * a) as f64,
                b_sq: (b * b) as f64,
                prod: (a * b) as f64,
            })
            .enumerate()
            .for_each(|(i, elem)| {
                buffer[i] = elem;
            });
        Self {
            a_sum_sq: buffer[..a.len()].iter().map(|x| x.a_sq).sum(),
            b_sum_sq: buffer[..a.len()].iter().map(|x| x.b_sq).sum(),
            dot_product: buffer[..a.len()].iter().map(|x| x.prod).sum(),
            buffer,
            size: a.len(),
            abs_index: a.len(),
            curr_index: 0,
            recalculations: 0,
        }
    }

    fn recalculate(&mut self) {
        self.a_sum_sq = self.buffer[..self.size].iter().map(|x| x.a_sq).sum();
        self.b_sum_sq = self.buffer[..self.size].iter().map(|x| x.b_sq).sum();
        self.dot_product = self.buffer[..self.size].iter().map(|x| x.prod).sum();
        self.recalculations += 1;
    }

    fn calculate_similarity(&mut self) -> f32 {
        let mut mag_a = self.a_sum_sq.sqrt();
        let mut mag_b = self.b_sum_sq.sqrt();
        let mut denom = mag_a * mag_b;
        if denom == 0.0 {
            return 0.0;
        }
        if self.dot_product > denom {
            self.recalculate();
            mag_a = self.a_sum_sq.sqrt();
            mag_b = self.b_sum_sq.sqrt();
            denom = mag_a * mag_b;
        }
        let out = self.dot_product / denom;

        if out > 1.1 {
            // TODO: Once I am happy with the numeric stability errors I should make this panic
            if out > 2.0 {
                println!(
                    "Cosine similarity out of bounds, got {}, state: {:#?}",
                    out, self
                );
            }
            1.0
        } else {
            // Due to numerical errors, the cosine similarity can be slightly above 1
            // And we dont really care about similarities < 1
            out.clamp(0.0, 1.0) as f32
        }
    }

    fn update(&mut self, a: f32, b: f32) {
        self.abs_index += 1;
        let a: f64 = a as f64;
        let b: f64 = b as f64;
        let prod = a * b;

        let curr_elem = RollingElem {
            a_sq: a * a,
            b_sq: b * b,
            prod,
        };

        self.a_sum_sq += curr_elem.a_sq - self.buffer[self.curr_index].a_sq;
        self.b_sum_sq += curr_elem.b_sq - self.buffer[self.curr_index].b_sq;
        self.dot_product += curr_elem.prod - self.buffer[self.curr_index].prod;

        // Recalculate every 5x size ... this helps with numerical errors
        // Ideally we would never recalculate but I have not been able to
        // manage the stability of the sums.
        if self.abs_index % (5 * self.size) == 0 {
            self.recalculate();
        }
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
///     0.97463184, // Note that this is the same as above
///     0.97823,
///     0.95065,
///     0.97463184, // And this one
///     0.97823,
///     0.95065,
///     0.97463184,
///     f32::NAN,
/// ];
/// let results = rolling_cosine_similarity(&a, &b, 3).unwrap();
/// assert_eq!(results.len(), expect_res.len());
/// ```
pub fn rolling_cosine_similarity(a: &[f32], b: &[f32], window_size: usize) -> Result<Vec<f32>> {
    // Check if vectors have the same length and are long enough for the window
    if a.len() != b.len() {
        return Err(DataProcessingError::ExpectedSlicesSameLength {
            expected: a.len(),
            other: b.len(),
            context: "cosine_similarity".to_string(),
        }
        .into());
    }
    if a.len() < window_size {
        return Err(DataProcessingError::ExpectedNonEmptyData {
            context: Some("cosine_similarity".to_string()),
        }
        .into());
    }

    let offset = window_size / 2;
    let mut results = vec![f32::NAN; a.len()];

    // Initialize the first window
    let mut cosine_sim =
        CosineSimilarityCircularBuffer::new(&a[0..window_size], &b[0..window_size]);

    // Calculate similarity for first window
    results[offset] = cosine_sim.calculate_similarity();

    // Roll the window
    for i in (offset + 1)..(a.len() - offset) {
        cosine_sim.update(a[i + offset], b[i + offset]);
        results[i] = cosine_sim.calculate_similarity();
    }

    Ok(results)
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
        let expect_res: [f32; 4] = [f32::NAN, 1.0, 1.0, f32::NAN];
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
        let results = rolling_cosine_similarity(&a, &b, 2).unwrap();
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
                assert!((result - expect).abs() < 1e-4);
            }
        }
    }
}
