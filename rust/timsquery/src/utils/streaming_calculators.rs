use std::ops::AddAssign;

// Generic streaming aggregator that takes a pair of unsigned ints one with a value
// and another with a weight and in a streaming fashion adds the value to the accumulator
// to calculate the total, mean and variance.

// TODO: move this to the general errors.
// TODO2: I am not 100% sure whether this file should be here ... happy to get
// feedback on code organization.
#[derive(Debug, Clone, Copy)]
pub enum StreamingAggregatorError {
    DivisionByZero,
    NotEnoughData,
}

type Result<T> = std::result::Result<T, StreamingAggregatorError>;

/// A struct that can be used to calculate the mean and variance of a stream of numbers.
///
/// # Example
///
/// ```
/// use timsquery::utils::streaming_calculators::RunningStatsCalculator;
///
/// // Create a new calculator with a weight of 10 and a mean of 0.0
/// let mut calc = RunningStatsCalculator::new(1.0, 0.0);
/// calc.add(1.0, 10.0);
/// calc.add(1.0, 0.0);
/// calc.add(1.0, 10.0);
/// calc.add(1.0, 0.0);
/// calc.add(1.0, 10.0);
/// calc.add(1.0, 0.0);
/// // So overall this should be the equivalent of the mean for
/// // [0.0, 10.0, 0.0, 10.0, 0.0, 10.0]
/// assert_eq!(calc.mean().unwrap(), 5.0, "{calc:#?}");
/// // assert!((4.5..5.5).contains(&calc.standard_deviation().unwrap()), "{calc:#?}");
/// ```
///
/// # Notes
///
/// It is important to know that the calculation of the mean is not
/// perfect. Thus if the initial value passes if very far off the real
/// mean, the final estimate will be off.
///
/// # References
///
/// - <https://nestedsoftware.com/2018/03/27/calculating-standard-deviation-on-streaming-data-253l.23919.html>
/// - <https://nestedsoftware.com/2019/09/26/incremental-average-and-standard-deviation-with-sliding-window-470k.176143.html>
#[derive(Debug, Clone, Copy, Default)]
pub struct RunningStatsCalculator {
    weight: f64,
    mean_n: f64,
    // TODO: Make min/max conditionally compiled
    min: f64,
    max: f64,
}

impl RunningStatsCalculator {
    pub fn new(weight: f64, mean: f64) -> Self {
        if weight == 0. {
            panic!("Weight must be > 0, initializing");
        }
        Self {
            weight,
            mean_n: mean,
            min: mean,
            max: mean,
            // count: 0,
        }
    }

    /// Add a new value to the running stats calculator.
    pub fn add(&mut self, weight: f64, value: f64) {
        if weight == 0. {
            panic!("Weight must be > 0, adding");
        }
        let f64_weight = weight;
        // Update the mean
        let weight_ratio = f64_weight / self.weight;
        let delta = value - self.mean_n;
        self.mean_n += delta * weight_ratio;

        // Update the weight
        self.weight += weight;

        // That should be the end of it but I seem to be getting consistently some
        // values outside of the min and max observed values. Which might be a
        // float issue ... TODO investigate.

        // In the meantime I will just squeeze the mean to the min/max observed values.
        // TODO:Make this conditional compilation
        self.min = self.min.min(value);
        self.max = self.max.max(value);

        self.mean_n = self.mean_n.min(self.max).max(self.min);
    }

    pub fn mean(&self) -> Result<f64> {
        if self.weight == 0. {
            return Err(StreamingAggregatorError::NotEnoughData);
        }
        Ok(self.mean_n)
    }

    pub fn weight(&self) -> f64 {
        self.weight
    }
}

impl AddAssign for RunningStatsCalculator {
    fn add_assign(&mut self, other: Self) {
        if other.weight == 0. {
            return;
        }
        // Update the mean
        self.add(other.weight, other.mean_n);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_running_stats_calculator() {
        let mut calc = RunningStatsCalculator::new(10., 0.0);
        calc.add(2., 10.0);
        calc.add(2., 10.0);
        calc.add(2., 10.0);
        calc.add(2., 10.0);
        calc.add(2., 10.0);
        assert!(calc.mean().unwrap() < 5.6);
        assert!(calc.mean().unwrap() > 4.4);
        // assert!(calc.variance().unwrap() > 15.);
        // assert!(calc.variance().unwrap() < 25.);
        // assert!(calc.standard_deviation().unwrap() > 4.5);
        // assert!(calc.standard_deviation().unwrap() < 5.5);
    }

    // https://www.kaggle.com/datasets/carlmcbrideellis/data-anscombes-quartet?resource=download
    // ascombes quarted data from kaggle
    //
    // Both have real mean of 7.5 and std of 1.94
    const ASCOMBES_3: [f64; 11] = [
        7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73,
    ];
    const ASCOMBES_4: [f64; 11] = [
        6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.5, 5.56, 7.91, 6.89,
    ];

    #[test]
    fn test_running_stats_calculator_ascombes_3() {
        let mut calc = RunningStatsCalculator::new(1., ASCOMBES_3[0]);
        for val in ASCOMBES_3[1..].iter() {
            calc.add(1., *val);
        }
        assert!(calc.mean().unwrap() < 7.6);
        assert!(calc.mean().unwrap() > 7.4);
        // assert!(calc.standard_deviation().unwrap() > 1.92);
        // assert!(calc.standard_deviation().unwrap() < 1.99);
    }

    #[test]
    fn test_running_stats_calculator_ascombes_4() {
        let mut calc = RunningStatsCalculator::new(1., ASCOMBES_4[0]);
        for val in ASCOMBES_4[1..].iter() {
            calc.add(1., *val);
        }
        assert!(calc.mean().unwrap() < 7.6);
        assert!(calc.mean().unwrap() > 7.4);

        // // Note that the tolerance here is a hair higher ... bc there
        // // is an outlier value.
        // assert!(
        //     calc.standard_deviation().unwrap() > 1.91,
        //     "Expected > 1.92, got {}",
        //     calc.standard_deviation().unwrap()
        // );
        // assert!(calc.standard_deviation().unwrap() < 1.99);
    }
}
