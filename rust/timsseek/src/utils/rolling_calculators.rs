use std::ops::Neg;

use arrayvec::ArrayVec;
use tracing::warn;

const MAX_WINDOW_SIZE: usize = 100;

pub struct RollingMedianCalculator<T: PartialOrd + Copy + Clone> {
    window_size: usize,
    data: ArrayVec<(T, usize), MAX_WINDOW_SIZE>,
    index: usize,
}

impl<T: PartialOrd + Copy + Clone> RollingMedianCalculator<T> {
    pub fn new(window_size: usize) -> Self {
        let mut window_size_use = window_size;
        if window_size > MAX_WINDOW_SIZE {
            warn!(
                "Window size {} is larger than max size {}. Clamping to max size.",
                window_size, MAX_WINDOW_SIZE
            );
            window_size_use = MAX_WINDOW_SIZE;
        }
        Self {
            window_size: window_size_use,
            data: ArrayVec::new(),
            index: 0,
        }
    }

    pub fn add(&mut self, value: T) {
        if self.data.len() > (self.window_size - 1) {
            let min_index_keep = (self.index - self.window_size) + 1;
            self.data.retain(|x| x.1 >= min_index_keep);
            assert!(
                self.data.len() < self.window_size,
                "Expected to trim to less than {}, but got {}",
                self.window_size,
                self.data.len(),
            );
        }
        self.insert_in_position((value, self.index));
        // self.push((value, self.index));
        self.index += 1;
        // self.reorder();
    }

    // fn reorder(&mut self) {
    //     // self.data
    //     //     .sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    //     // Since only the last element can be in the wrong position, we can do
    //     // an insertion sort step.
    // }

    fn insert_in_position(&mut self, last: (T, usize)) {
        if self.data.len() < 2 {
            self.data.push(last);
            return;
        }
        let mut pos = self.data.len();
        while pos > 0 && last.0 < self.data[pos - 1].0 {
            pos -= 1;
        }
        self.data.insert(pos, last);
    }

    pub fn median(&self) -> Option<T> {
        if self.data.len() < self.window_size {
            None
        } else {
            let median_index = self.data.len() / 2;
            Some(self.data[median_index].0)
        }
    }
}

pub fn rolling_median_into<T: PartialOrd + Copy + Clone>(
    values: &[T],
    window_size: usize,
    pad_value: T,
    out: &mut Vec<T>,
) {
    out.clear();
    out.resize(values.len(), pad_value);
    let mut rolling = RollingMedianCalculator::new(window_size);
    let offset = window_size / 2;
    for (i, value) in values.iter().enumerate() {
        rolling.add(*value);
        if i >= (window_size - 1) {
            out[i - offset] = rolling.median().unwrap();
        }
    }
}

pub fn calculate_value_vs_baseline_into(
    vals: &[f32],
    baseline_window_size: usize,
    out: &mut Vec<f32>,
) {
    rolling_median_into(vals, baseline_window_size, f32::NAN, out);
    for (x, bx) in vals.iter().zip(out.iter_mut()) {
        *bx = bx.neg();
        *bx += x;
    }
}

/// Calculate the centered standard deviation
///
/// This is a standard deviation for a slice, where we can assume
/// the mean is 0.
pub fn calculate_centered_std(vals: &[f32]) -> f32 {
    let sqsum = vals
        .iter()
        .filter(|x| !x.is_nan())
        .map(|x| x.powi(2))
        .sum::<f32>();
    let variance = sqsum / vals.len() as f32;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rolling_median_calculator() {
        let mut calc = RollingMedianCalculator::new(3);
        calc.add(10.0);
        calc.add(20.0);
        calc.add(30.0);
        assert_eq!(calc.median(), Some(20.0));
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        assert_eq!(calc.median(), Some(1.0));
    }

    #[test]
    fn test_rolling_median() {
        let input = vec![1.0, 2.0, 30.0, 4.0, 5.0, 60.0, 7.0, 8.0, 9.0];
        let mut out = Vec::new();
        rolling_median_into(&input, 3, f64::NAN, &mut out);
        let expect_out = vec![f64::NAN, 2.0, 4.0, 5.0, 5.0, 7.0, 8.0, 8.0, f64::NAN];

        // assert_eq!(
        //     out,
        //     expect_out,
        // );
        //
        for i in 0..out.len() {
            if expect_out[i].is_nan() {
                assert!(out[i].is_nan());
            } else {
                assert!(
                    (out[i] - expect_out[i]).abs() < 1e-6,
                    "Expected {:?}, got {:?}",
                    expect_out,
                    out
                );
            }
        }
    }

    #[test]
    fn test_value_vs_baseline() {
        let vals: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let baseline_window_size = 3;
        let mut _baseline = Vec::new();
        rolling_median_into(&vals, baseline_window_size, f32::NAN, &mut _baseline);
        let mut out = Vec::new();
        calculate_value_vs_baseline_into(&vals, baseline_window_size, &mut out);
        let expect_val = vec![f32::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, f32::NAN];
        let all_close = out
            .iter()
            .zip(expect_val.iter())
            .filter(|(a, b)| ((!a.is_nan()) && (!b.is_nan())))
            .all(|(a, b)| (a - b).abs() < 1e-6);

        let all_match_nan = out
            .iter()
            .zip(expect_val.iter())
            .filter(|(a, b)| ((a.is_nan()) || (b.is_nan())))
            .all(|(a, b)| a.is_nan() && b.is_nan());

        assert!(all_close, "Expected {:?}, got {:?}", expect_val, out);
        assert!(all_match_nan, "Expected {:?}, got {:?}", expect_val, out);
    }
}
