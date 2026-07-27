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
    RidgeMeasurement,
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

/// Ridge half-width, weighted by each column's in-ridge weight. NaN when there
/// are no measurements — heavier columns should carry more authority, and an
/// unweighted mean over an empty set is not a number.
pub fn weighted_ridge_half_width(widths: &[RidgeMeasurement]) -> f64 {
    let total: f64 = widths.iter().map(|m| m.ridge_weight).sum();
    if widths.is_empty() || total <= 0.0 {
        return f64::NAN;
    }
    widths
        .iter()
        .map(|m| m.half_width * m.ridge_weight)
        .sum::<f64>()
        / total
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
            score: 1.0,
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

    #[test]
    fn churn_on_an_unchanged_set_is_zero() {
        let prev = [pt(1), pt(2)];
        assert_eq!(churn(&prev, &prev), (0, 0));
    }

    #[test]
    fn churn_from_empty_is_all_admissions() {
        assert_eq!(churn(&[], &[pt(1), pt(2)]), (2, 0));
    }

    #[test]
    fn identical_curves_have_zero_delta() {
        let c = curve_from(&[(0.0, 0.0), (10.0, 10.0)]);
        let (max, mean) = curve_delta(&c, &c, (0.0, 10.0), 11);
        assert_eq!((max, mean), (0.0, 0.0));
    }

    #[test]
    fn a_constant_offset_shows_up_as_that_offset() {
        let a = curve_from(&[(0.0, 0.0), (10.0, 10.0)]);
        let b = curve_from(&[(0.0, 2.0), (10.0, 12.0)]);
        let (max, mean) = curve_delta(&a, &b, (0.0, 10.0), 11);
        assert!((max - 2.0).abs() < 1e-9);
        assert!((mean - 2.0).abs() < 1e-9);
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

    #[test]
    fn weighted_ridge_half_width_weights_by_ridge_weight() {
        let widths = vec![
            RidgeMeasurement {
                library: LibraryRT(1.0),
                half_width: 10.0,
                ridge_weight: 1.0,
                column_weight: 1.0,
            },
            RidgeMeasurement {
                library: LibraryRT(2.0),
                half_width: 20.0,
                ridge_weight: 3.0,
                column_weight: 3.0,
            },
        ];
        // (10*1 + 20*3) / 4 = 17.5
        assert!((weighted_ridge_half_width(&widths) - 17.5).abs() < 1e-9);
    }

    #[test]
    fn weighted_ridge_half_width_of_nothing_is_nan() {
        assert!(weighted_ridge_half_width(&[]).is_nan());
    }
}
