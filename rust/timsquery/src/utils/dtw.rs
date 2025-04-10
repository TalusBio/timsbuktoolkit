use crate::errors::{
    DataProcessingError,
    Result,
};
use std::ops::{
    Add,
    Sub,
};

/// Defines how values should be combined when warping
pub trait WarpingStrategy {
    /// Combines two values according to the warping strategy
    fn combine(&self, current: f64, new_value: f64) -> f64;

    /// Provides the initial value for the warping operation
    fn initial_value(&self) -> f64;
}

pub trait WarpingX<
    T: PartialOrd + Copy + Clone + std::ops::Sub<Output = T> + std::ops::Add<Output = T>,
>
{
    fn maxval() -> T;
}

impl WarpingX<f64> for f64 {
    fn maxval() -> f64 {
        f64::MAX
    }
}

impl WarpingX<u64> for u64 {
    fn maxval() -> u64 {
        u64::MAX
    }
}

impl WarpingX<u32> for u32 {
    fn maxval() -> u32 {
        u32::MAX
    }
}

/// Strategy that keeps the maximum value when warping
#[derive(Debug, Clone, Copy)]
pub struct MaxWarping;

impl WarpingStrategy for MaxWarping {
    fn combine(&self, current: f64, new_value: f64) -> f64 {
        current.max(new_value)
    }

    fn initial_value(&self) -> f64 {
        f64::NEG_INFINITY
    }
}

/// Strategy that adds values when warping
#[derive(Debug, Clone, Copy)]
pub struct AddWarping;

impl WarpingStrategy for AddWarping {
    fn combine(&self, current: f64, new_value: f64) -> f64 {
        if current == self.initial_value() {
            new_value
        } else {
            current + new_value
        }
    }

    fn initial_value(&self) -> f64 {
        0.0
    }
}

/// Generic Dynamic Time Warping implementation that uses a specified warping strategy
pub fn dtw_with_strategy<
    S: WarpingStrategy,
    T: WarpingX<T> + PartialOrd + Copy + Clone + Sub<Output = T> + Add<Output = T>,
>(
    ref_x: &[T],
    target_x: &[T],
    target_vals: &[f64],
    strategy: &S,
) -> Result<Vec<f64>> {
    // Check if vectors have the same length and are not empty
    if target_vals.len() != target_x.len() {
        return Err(DataProcessingError::ExpectedVectorSameLength(
            target_vals.len(),
            target_x.len(),
        )
        .into());
    }

    // Initialize output with strategy's initial value
    let mut output = vec![strategy.initial_value(); ref_x.len()];
    let defval = T::maxval();

    let mut target_i = 0;
    for (ref_i, ref_x_val) in ref_x.iter().enumerate() {
        let next_ref_val = ref_x.get(ref_i + 1).unwrap_or(&defval);
        // 1. If out target is lower than the current reference add.
        // 2. If our value is between the current and next reference add to the closest one.
        // 3. If our value is larget than the next reference, skip.
        // Once we exhaust the target values add everything to the last one.

        while target_i < target_x.len() && target_x[target_i] < *ref_x_val {
            output[ref_i] = strategy.combine(output[ref_i], target_vals[target_i]);
            target_i += 1;
        }

        while target_i < target_x.len() && target_x[target_i] < *next_ref_val {
            let left_diff = target_x[target_i] - *ref_x_val;
            let right_diff = *next_ref_val - target_x[target_i];
            if left_diff < right_diff {
                output[ref_i] = strategy.combine(output[ref_i], target_vals[target_i]);
            } else {
                output[ref_i + 1] = strategy.combine(output[ref_i + 1], target_vals[target_i]);
            }
            target_i += 1;
        }
    }

    Ok(output)
}

/// Dynamic Time Warping to max
///
/// Attempts to align the values in `target_vals` which correspond
/// to the positions in `target_x` to the positions in `ref_x`.
///
/// This version preseves the maximum value of the `target_vals`
/// in each position of the `target_x` array.
///
/// # Note
/// This assumes that `ref_x` and `target_x` are sorted (monotonically increasing).
///
/// Thus .... If the reference is shorter:
///
/// If ref_x    =       [1, 2, 3, 4, 5, 6, 7, 8, 10]
/// If target_x =       [1, 2,    4, 5, 6, 7,    10]
/// If target_vals =    [5, 3,    2, 1, 2, 4,     9]
///
/// Then the output would be:
///                     [5, 3, 3, 2, 1, 2, 4, 4,  9] // Length = ref_x.len()
///
/// Thus .... If the reference is longer:
///
/// If ref_x    =       [1, 2,    4, 5, 6, 7,    10]
/// If target_x =       [1, 2, 3, 4, 5, 6, 7, 8, 10]
/// If target_vals =    [5, 3, 1, 2, 1, 2, 4, 4,  9]
///
/// Then the output would be:
///                     [5, 3,    2, 1, 2, 4,     9] // Length = ref_x.len()
pub fn dtw_max<T: WarpingX<T> + PartialOrd + Copy + Clone + Sub<Output = T> + Add<Output = T>>(
    ref_x: &[T],
    target_x: &[T],
    target_vals: &[f64],
) -> Result<Vec<f64>> {
    dtw_with_strategy(ref_x, target_x, target_vals, &MaxWarping)
}

/// Dynamic Time Warping to Sum
///
/// Attempts to align the values in `target_vals` which correspond
/// to the positions in `target_x` to the positions in `ref_x`.
///
/// This version preseves the maximum value of the `target_vals`
/// in each position of the `target_x` array.
///
/// # Note
/// This assumes that `ref_x` and `target_x` are sorted (monotonically increasing).
///
/// Thus .... If the reference is shorter:
///
/// If ref_x    =       [1, 2, 3, 4, 5, 6, 7, 8, 10]
/// If target_x =       [1, 2,    4, 5, 6, 7,    10]
/// If target_vals =    [5, 3,    2, 1, 2, 4,     9]
///
/// Then the output would be:
///                     [5, 3, 3, 2, 1, 2, 4, 4,  9] // Length = ref_x.len()
///
/// Thus .... If the reference is longer:
///
/// If ref_x    =       [1, 2,    4, 5, 6, 7,    10]
/// If target_x =       [1, 2, 3, 4, 5, 6, 7, 8, 10]
/// If target_vals =    [5, 3, 1, 2, 1, 2, 4, 4,  9]
///
/// Then the output would be:
///                     [5, 4,    2, 1, 2, 4,     13] // Length = ref_x.len()
pub fn dtw_add<T: WarpingX<T> + PartialOrd + Copy + Clone + Sub<Output = T> + Add<Output = T>>(
    ref_x: &[T],
    target_x: &[T],
    target_vals: &[f64],
) -> Result<Vec<f64>> {
    dtw_with_strategy(ref_x, target_x, target_vals, &AddWarping)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[rustfmt::skip]
    fn test_dtw_max_reference_shorter() {
        let ref_x =       vec![1.0, 2.0,      4.0, 5.0, 6.0, 7.0,      10.0];
        let target_x =    vec![1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 10.1];
        let target_vals = vec![5.0, 3.0, 1.0, 2.0, 1.0, 2.0, 4.0, 11.0,9.0];
        //                               ->                       <-
        let expected    = vec![5.0, 3.0,      2.0, 1.0, 2.0, 11.0,     9.0];

        let result = dtw_max(&ref_x, &target_x, &target_vals).unwrap();
        assert_eq!(ref_x.len(), result.len());
        assert_eq!(result, expected, "Expected {:?} got {:?}", result, expected);
    }

    #[test]
    #[rustfmt::skip]
    fn test_dtw_max_reference_longer() {
        let ref_x       = vec![1.0, 2.0,      4.0, 5.0, 6.0, 7.0,      10.0];
        let target_x    = vec![1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 10.1];
        let target_vals = vec![5.0, 1.0, 3.0, 2.0, 1.0, 2.0, 7.0, 11.0, 9.0];
        //                               ->                       <-
        let expected    = vec![5.0, 1.0,      3.0, 1.0, 2.0, 11.0,     9.0];

        let result = dtw_max(&ref_x, &target_x, &target_vals).unwrap();
        assert_eq!(ref_x.len(), result.len());
        assert_eq!(result, expected, "Expected {:?} got {:?}", result, expected);
    }

    #[test]
    #[rustfmt::skip]
    fn test_dtw_add_reference_shorter() {
        let ref_x =        vec![1.0, 2.0,      4.0, 5.0, 6.0, 7.0,      10.0];
        let target_x =     vec![1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 10.1];
        let target_vals =  vec![5.0, 3.0, 1.0, 2.0, 1.0, 2.0, 4.0, 4.0, 9.0];
        //                                 ->                       <-
        let expect_out =   vec![5.0, 3.0,      3.0, 1.0, 2.0, 8.0,      9.0];
        // Here 3 and 4 should add to the same value (of 4), same for 7 and 8.

        let result = dtw_add(&ref_x, &target_x, &target_vals).unwrap();
        assert_eq!(ref_x.len(), result.len());
        // Values get added when they map to the same reference point
        for (l, r) in result.iter().zip(expect_out.iter()) {
            assert!(
                (l - r).abs() < 1e-6,
                "Expected {:?} got {:?}",
                result,
                expect_out
            );
        }
    }

    #[test]
    fn test_custom_strategy() {
        // Example of a custom strategy that multiplies values
        struct MultiplyWarping;
        impl WarpingStrategy for MultiplyWarping {
            fn combine(&self, current: f64, new_value: f64) -> f64 {
                if current == self.initial_value() {
                    new_value
                } else {
                    current * new_value
                }
            }

            fn initial_value(&self) -> f64 {
                1.0
            }
        }

        let ref_x = vec![1.0, 2.0, 3.0];
        let target_x = vec![1.0, 1.1, 2.0, 2.1, 3.0];
        let target_vals = vec![2.0, 2.0, 3.0, 10.0, 4.0];

        let result = dtw_with_strategy(&ref_x, &target_x, &target_vals, &MultiplyWarping).unwrap();
        let expected_out = vec![4.0, 30.0, 4.0];
        for (l, r) in result.iter().zip(expected_out.iter()) {
            assert!(
                (l - r).abs() < 1e-6,
                "Expected {:?} got {:?}",
                result,
                expected_out
            );
        }
    }
}
