//! Per-batch convergence metrics: the series a future early-stopping rule
//! would threshold on.
//!
//! These are computed at every batch, including ones the UI skips past, so the
//! series has no holes. Each is five floats, so keeping all of them costs
//! nothing next to the frame slab.

use crate::CalibrantPoint;
use calibrt::{
    CalibrationCurve,
    LibraryRT,
};
use std::collections::HashSet;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BatchMetrics {
    pub chunk: usize,
    pub n_points: usize,
    pub wrmse: f64,
    pub max_delta: f64,
    pub mean_delta: f64,
    pub path_nodes: usize,
    pub ridge_half_width: f64,
    pub admitted: usize,
    pub evicted: usize,
}

/// `(max, mean)` of `|cur.predict(x) - prev.predict(x)|` over `samples` evenly
/// spaced points in `x_range`.
///
/// Samples where either curve is out of bounds are skipped rather than counted
/// as zero — counting them would dilute the mean toward zero exactly when the
/// curves disagree most about their domain, which is the opposite of the signal
/// wanted here.
pub fn curve_delta(
    prev: &CalibrationCurve,
    cur: &CalibrationCurve,
    x_range: (f64, f64),
    samples: usize,
) -> (f64, f64) {
    let samples = samples.max(2);
    let span = x_range.1 - x_range.0;
    let mut max = 0.0f64;
    let mut sum = 0.0f64;
    let mut n = 0usize;
    for k in 0..samples {
        let x = x_range.0 + span * k as f64 / (samples - 1) as f64;
        let (Ok(a), Ok(b)) = (prev.predict(LibraryRT(x)), cur.predict(LibraryRT(x))) else {
            continue;
        };
        let d = (b.0 - a.0).abs();
        max = max.max(d);
        sum += d;
        n += 1;
    }
    if n == 0 {
        (f64::NAN, f64::NAN)
    } else {
        (max, sum / n as f64)
    }
}

/// `(admitted, evicted)` between two heap snapshots, by `speclib_index`.
pub fn churn(prev: &[CalibrantPoint], cur: &[CalibrantPoint]) -> (usize, usize) {
    let prev_set: HashSet<usize> = prev.iter().map(|p| p.speclib_index).collect();
    let cur_set: HashSet<usize> = cur.iter().map(|p| p.speclib_index).collect();
    let admitted = cur_set.difference(&prev_set).count();
    let evicted = prev_set.difference(&cur_set).count();
    (admitted, evicted)
}

#[cfg(test)]
mod tests {
    use super::*;
    use calibrt::{
        LibraryRT,
        ObservedRTSeconds,
    };

    fn pt(idx: usize) -> CalibrantPoint {
        CalibrantPoint {
            library_rt: 1.0,
            observed_rt: 1.0,
            speclib_index: idx,
        }
    }

    /// Builds a `CalibrationCurve` whose path is exactly `points`, by running a
    /// tiny real fit (`CalibrationCurve::new` is `pub(crate)` in `calibrt`).
    ///
    /// Preconditions: every weight at or above 1.0, and pairs strictly
    /// increasing in both coordinates, one per grid cell.
    fn curve_from(points: &[(f64, f64)]) -> CalibrationCurve {
        let n = points.len();
        let grid_size = (n * 4).max(8);
        let x_min = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let x_max = points.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
        let y_min = points.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let y_max = points.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);

        let mut state =
            calibrt::CalibrationState::new(grid_size, (x_min, x_max), (y_min, y_max), n)
                .expect("fixture ranges are nonzero");
        state
            .update(
                points
                    .iter()
                    .map(|&(x, y)| (LibraryRT(x), ObservedRTSeconds(y), 2.0)),
            )
            .expect("fixture points are finite");
        state.fit();
        let curve = state
            .curve()
            .expect("fixture points are monotone and well-separated, so a path always exists")
            .clone();
        // Guards the one-point-per-cell precondition: two points in one cell
        // would be centroid-merged into a shorter curve.
        assert_eq!(curve.points().len(), points.len());
        curve
    }

    #[test]
    fn churn_counts_admissions_and_evictions_by_index() {
        let prev = [pt(1), pt(2), pt(3)];
        let cur = [pt(2), pt(3), pt(4), pt(5)];
        assert_eq!(churn(&prev, &cur), (2, 1), "4 and 5 admitted, 1 evicted");
    }

    /// Both signs of the same offset. The delta is a *magnitude*, so a curve
    /// that moved down by 2 must report the same 2 as one that moved up by 2 —
    /// without the `.abs()`, the downward direction reports `max = 0` (the
    /// running max never rises above its 0.0 seed) and a negative mean.
    #[test]
    fn a_constant_offset_shows_up_as_that_offset_in_either_direction() {
        let low = curve_from(&[(0.0, 0.0), (10.0, 10.0)]);
        let high = curve_from(&[(0.0, 2.0), (10.0, 12.0)]);

        let (max, mean) = curve_delta(&low, &high, (0.0, 10.0), 11);
        assert!((max - 2.0).abs() < 1e-9, "moved up by 2, got max {max}");
        assert!((mean - 2.0).abs() < 1e-9, "moved up by 2, got mean {mean}");

        let (max, mean) = curve_delta(&high, &low, (0.0, 10.0), 11);
        assert!((max - 2.0).abs() < 1e-9, "moved down by 2, got max {max}");
        assert!(
            (mean - 2.0).abs() < 1e-9,
            "moved down by 2, got mean {mean}"
        );
    }

    /// `samples` is clamped to 2 because the sample positions divide by
    /// `samples - 1`. No production caller passes below `CURVE_DELTA_SAMPLES`,
    /// so nothing else would notice the clamp going missing — at 1 the position
    /// arithmetic becomes `0.0 / 0.0` and every sample lands on NaN, and at 0
    /// the loop never runs at all; both silently degrade to `(NaN, NaN)`.
    #[test]
    fn too_few_samples_are_clamped_to_the_two_endpoints() {
        let low = curve_from(&[(0.0, 0.0), (10.0, 10.0)]);
        let high = curve_from(&[(0.0, 2.0), (10.0, 12.0)]);
        for samples in [0, 1, 2] {
            let (max, mean) = curve_delta(&low, &high, (0.0, 10.0), samples);
            assert!(
                (max - 2.0).abs() < 1e-9 && (mean - 2.0).abs() < 1e-9,
                "{samples} samples must still measure the constant offset over \
                 the two endpoints, got ({max}, {mean})"
            );
        }
    }

    #[test]
    fn max_delta_exceeds_mean_when_the_change_is_local() {
        let a = curve_from(&[(0.0, 0.0), (5.0, 5.0), (10.0, 10.0)]);
        let b = curve_from(&[(0.0, 0.0), (5.0, 8.0), (10.0, 10.0)]);
        let (max, mean) = curve_delta(&a, &b, (0.0, 10.0), 21);
        assert!(max > mean, "max {max} should exceed mean {mean}");
        assert!((max - 3.0).abs() < 1e-6);
    }

    #[test]
    fn samples_outside_both_curves_are_skipped_not_counted_as_zero() {
        // b covers a narrower x range; samples beyond it must not dilute the mean.
        let a = curve_from(&[(0.0, 0.0), (10.0, 10.0)]);
        let b = curve_from(&[(0.0, 3.0), (5.0, 8.0)]);
        let (max, mean) = curve_delta(&a, &b, (0.0, 10.0), 11);
        assert!((max - 3.0).abs() < 1e-6);
        assert!(
            (mean - 3.0).abs() < 1e-6,
            "only the overlapping half contributes"
        );
    }

    #[test]
    fn disjoint_domains_never_overlap_and_delta_is_nan() {
        // a and b's x ranges don't intersect at all, so every sample is
        // out-of-bounds for at least one curve: n stays 0 and the result must
        // be NaN, not 0.0 (which would read as "the curves agree").
        let a = curve_from(&[(0.0, 0.0), (5.0, 5.0)]);
        let b = curve_from(&[(20.0, 20.0), (25.0, 25.0)]);
        let (max, mean) = curve_delta(&a, &b, (0.0, 5.0), 11);
        assert!(max.is_nan(), "max should be NaN, got {max}");
        assert!(mean.is_nan(), "mean should be NaN, got {mean}");
    }
}
