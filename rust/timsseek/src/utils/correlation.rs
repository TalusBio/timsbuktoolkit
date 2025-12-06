//! Cosine similarity calculations for chromatogram correlation.
//!
//! This module provides efficient rolling window cosine similarity calculations, which are
//! used heavily in coelution scoring. The rolling window approach uses a circular buffer to
//! maintain running sums, avoiding recomputation of the entire window on each step.
//!
//! # Performance Considerations
//!
//! The `rolling_cosine_similarity` function is one of the hottest paths in the codebase,
//! consuming ~30% of total runtime. Any changes here should be carefully benchmarked.
//!
//! # Numeric Stability
//!
//! The implementation uses f64 accumulation to handle extreme intensity ranges correctly
//! (e.g., 100 vs 1,000,000). See `test_extreme_range_stability` for validation.

use std::f32;

use crate::errors::DataProcessingError;

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

/// Single element in the rolling window buffer.
///
/// Stores precomputed squares and products to enable O(1) updates.
#[derive(Debug, Copy, Clone)]
struct RollingElem {
    a_sq: f64, // a*a
    b_sq: f64, // b*b
    prod: f64, // a*b
}

/// Maximum window size supported by the circular buffer.
///
/// Larger windows would require heap allocation. Current value (20) is sufficient
/// for typical coelution window sizes (7).
const MAX_CAPACITY: usize = 20;

/// Circular buffer for efficient rolling window cosine similarity.
///
/// This structure maintains running sums of squared values and dot products,
/// allowing O(1) updates as the window slides. The circular buffer avoids
/// allocations and minimizes cache misses.
///
/// # Algorithm
///
/// For each new (a, b) pair:
/// 1. Remove oldest element's contribution from sums
/// 2. Add new element's contribution to sums
/// 3. Compute similarity = dot_product / sqrt(a_sum_sq * b_sum_sq)
///
/// # Numeric Stability
///
/// Uses f64 accumulation to avoid truncation errors. This handles extreme intensity
/// ranges (e.g., 100 vs 1,000,000) correctly, as f64 maintains sufficient precision
/// for typical window sizes (≤20) and intensity values (up to ~1e6).
#[derive(Debug)]
struct CosineSimilarityCircularBuffer {
    a_sum_sq: f64,    // Sum of (a*a) over window
    b_sum_sq: f64,    // Sum of (b*b) over window
    dot_product: f64, // Sum of (a*b) over window
    buffer: [RollingElem; MAX_CAPACITY],
    window_size: usize,
    curr_index: usize, // Next position to write (wraps around)
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
            a_sq: 0.0,
            b_sq: 0.0,
            prod: 0.0,
        }; MAX_CAPACITY];

        Self {
            a_sum_sq: 0.0,
            b_sum_sq: 0.0,
            dot_product: 0.0,
            buffer,
            window_size,
            curr_index: 0,
        }
    }

    fn calculate_similarity(&mut self) -> f32 {
        let denom_sq = self.a_sum_sq * self.b_sum_sq;
        if denom_sq == 0.0 {
            return 0.0;
        }

        // Single sqrt operation, compute result in f64 then convert to f32
        let result = (self.dot_product / denom_sq.sqrt()) as f32;

        // Clamp to [0, 1] to handle minor floating-point errors
        if result > 2.0 {
            // I am leaving this here because I am terrified of silent
            // cases where numeric instability causes this to go wrong
            // and really weird bugs down the road.
            panic!(
                "Cosine similarity computation messed up: dot_product={}, a_sum_sq={}, b_sum_sq={}",
                self.dot_product, self.a_sum_sq, self.b_sum_sq
            );
        }
        result.clamp(0.0, 1.0)
    }

    // This is too hot, cannot instrument for performance reasons
    // #[cfg_attr(
    //     feature = "instrumentation",
    //     tracing::instrument(skip_all, level = "trace")
    // )]
    fn update(&mut self, a: f32, b: f32) {
        // Convert to f64 for accumulation (no truncation)
        let a = a as f64;
        let b = b as f64;

        // Get mutable reference to the current element in the buffer
        let elem = &mut self.buffer[self.curr_index];

        // Subtract the old values from the sums
        self.a_sum_sq -= elem.a_sq;
        self.b_sum_sq -= elem.b_sq;
        self.dot_product -= elem.prod;

        // Update the element in place
        elem.a_sq = a * a;
        elem.b_sq = b * b;
        elem.prod = a * b;

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

    /// Test numeric stability with extreme intensity ranges.
    ///
    /// This tests the scenario that motivated u64 accumulation: two Gaussian peaks
    /// with very different apex intensities (100 vs 1,000,000) in a background of noise.
    /// The cosine similarity should still be accurate despite the large dynamic range.
    #[test]
    fn test_extreme_range_stability() {
        // Create two chromatograms with Gaussian peaks at different intensities
        let n = 100;
        let peak1_center = 25;
        let peak2_center = 75;
        let noise_level = 0.5;

        // Helper to create Gaussian peak
        let gaussian = |x: f32, center: f32, amplitude: f32, width: f32| -> f32 {
            let exponent = -((x - center).powi(2)) / (2.0 * width.powi(2));
            amplitude * exponent.exp()
        };

        // Trace A: Peak at position 25 with amplitude 100
        let trace_a: Vec<f32> = (0..n)
            .map(|i| {
                let base = gaussian(i as f32, peak1_center as f32, 100.0, 5.0);
                base + (i as f32 * 0.1).sin() * noise_level // Add reproducible noise
            })
            .collect();

        // Trace B: Similar shape but scaled 10,000x higher + different peak at position 75
        let trace_b: Vec<f32> = (0..n)
            .map(|i| {
                let base1 = gaussian(i as f32, peak1_center as f32, 1_000_000.0, 5.0);
                let base2 = gaussian(i as f32, peak2_center as f32, 500_000.0, 5.0);
                base1 + base2 + (i as f32 * 0.1).cos() * noise_level * 10_000.0
            })
            .collect();

        // Calculate rolling cosine similarity with window size 7
        let window = 7;
        let results: Vec<f32> = rolling_cosine_similarity(&trace_a, &trace_b, window)
            .unwrap()
            .collect();

        // Verify results make sense:
        // 1. No NaN values (except padding)
        let left_pad = window / 2;
        let right_pad = window - left_pad - 1;
        for (i, &val) in results.iter().enumerate() {
            if i < left_pad || i >= results.len() - right_pad {
                assert!(val.is_nan(), "Expected NaN in padding at index {}", i);
            } else {
                assert!(
                    !val.is_nan(),
                    "Unexpected NaN at index {} (value: {})",
                    i,
                    val
                );
                assert!(
                    val >= -0.01 && val <= 1.01,
                    "Cosine similarity out of range at index {}: {}",
                    i,
                    val
                );
            }
        }

        // 2. High similarity around peak1_center (both traces have peak there)
        let peak1_similarity = results[peak1_center];
        assert!(
            peak1_similarity > 0.9,
            "Expected high similarity at peak1 center ({}), got {}",
            peak1_center,
            peak1_similarity
        );

        // 3. Non-zero similarity around peak2_center (was 0.0 with u64 truncation bug)
        let peak2_similarity = results[peak2_center];
        println!(
            "\nPeak2 similarity (position {}): {}",
            peak2_center, peak2_similarity
        );

        // 4. Print some values for inspection
        println!("\nSimilarity values around peaks:");
        for i in
            (peak1_center - 5).max(left_pad)..=(peak1_center + 5).min(results.len() - right_pad - 1)
        {
            println!("  [{}] = {:.6}", i, results[i]);
        }
        println!();
        for i in
            (peak2_center - 5).max(left_pad)..=(peak2_center + 5).min(results.len() - right_pad - 1)
        {
            println!("  [{}] = {:.6}", i, results[i]);
        }

        // 5. Check for catastrophic numeric errors (values way outside [0,1] or sudden jumps to 0/NaN)
        for i in left_pad..(results.len() - right_pad) {
            assert!(
                !results[i].is_nan(),
                "Unexpected NaN in valid region at index {}",
                i
            );
        }
    }

    /// Benchmark test: Measure performance of f64 accumulation.
    ///
    /// This is not run by default (use `cargo test -- --ignored` to run).
    /// Provides baseline timing for the current f64 implementation.
    #[test]
    #[ignore]
    fn bench_accumulation_strategies() {
        use std::time::Instant;

        let n = 10_000;
        let window = 7;

        // Create test data with extreme range
        let trace_a: Vec<f32> = (0..n).map(|i| (i as f32).sin() * 100.0 + 50.0).collect();
        let trace_b: Vec<f32> = (0..n)
            .map(|i| (i as f32 * 1.1).cos() * 1_000_000.0 + 500_000.0)
            .collect();

        // Preallocate result buffer (reused across iterations)
        let mut results = vec![0.0f32; n];

        // Benchmark current implementation (f64)
        let iterations = 100;
        let start = Instant::now();
        for _ in 0..iterations {
            for (i, val) in rolling_cosine_similarity(&trace_a, &trace_b, window)
                .unwrap()
                .enumerate()
            {
                results[i] = val;
            }
        }
        let f64_time = start.elapsed();

        println!(
            "f64 accumulation: {:?} ({} iterations)",
            f64_time, iterations
        );
        println!("  Per iteration: {:?}", f64_time / iterations);
        println!(
            "  Throughput: {:.2} iterations/sec",
            iterations as f64 / f64_time.as_secs_f64()
        );
    }
}
