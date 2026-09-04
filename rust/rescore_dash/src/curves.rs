//! The FDR and calibration curves: cumulative counts over a grid, built once in
//! `Dashboard::build`.

/// Targets passing as a function of the q-value threshold: one curve per entry
/// of `zooms`, curve `i` spanning `n_points` evenly spaced thresholds in
/// `(0, zooms[i]]`.
///
/// Computed together because the sort is shared: every zoom is answered by
/// `partition_point` over the same sorted array.
pub(crate) fn qvalue_curves(
    qvalue: &[f32],
    is_target: &[bool],
    n_points: usize,
    zooms: &[f64],
) -> Vec<Vec<(f64, f64)>> {
    let mut targets: Vec<f32> = qvalue
        .iter()
        .zip(is_target)
        .filter(|&(q, &t)| t && !q.is_nan())
        .map(|(&q, _)| q)
        .collect();
    targets.sort_unstable_by(f32::total_cmp);

    zooms
        .iter()
        .map(|&zoom| {
            (1..=n_points)
                .map(|k| {
                    let thresh = zoom * k as f64 / n_points as f64;
                    let n = targets.partition_point(|&q| (q as f64) <= thresh);
                    (thresh, n as f64)
                })
                .collect()
        })
        .collect()
}

/// Empirical decoy CDF against empirical target CDF at a shared grid of score
/// thresholds spanning the lowest to the highest finite score seen; the last
/// point is always `(1,1)`. The y=x reference is drawn by the caller.
///
/// When every finite score is identical, `lo == hi` and the grid still runs
/// over `n_points` steps (the span falls back to `1.0`), but every threshold
/// already covers both classes, so the curve is `(1,1)` at every point rather
/// than starting at `(0,0)`.
///
/// Empty when either class is absent -- a PP plot of one class says nothing.
pub(crate) fn pp_curve(score: &[f32], is_target: &[bool], n_points: usize) -> Vec<(f64, f64)> {
    let mut targets: Vec<f32> = Vec::new();
    let mut decoys: Vec<f32> = Vec::new();
    for (&s, &t) in score.iter().zip(is_target) {
        if !s.is_finite() {
            continue;
        }
        if t { targets.push(s) } else { decoys.push(s) }
    }
    if targets.is_empty() || decoys.is_empty() {
        return Vec::new();
    }
    targets.sort_unstable_by(f32::total_cmp);
    decoys.sort_unstable_by(f32::total_cmp);

    let lo = targets[0].min(decoys[0]) as f64;
    let hi = (*targets.last().unwrap()).max(*decoys.last().unwrap()) as f64;
    let span = if hi > lo { hi - lo } else { 1.0 };

    (0..=n_points)
        .map(|k| {
            let thresh = (lo + span * k as f64 / n_points as f64) as f32;
            let dc = decoys.partition_point(|&s| s <= thresh) as f64 / decoys.len() as f64;
            let tc = targets.partition_point(|&s| s <= thresh) as f64 / targets.len() as f64;
            (dc, tc)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qvalue_curve_is_monotone_nondecreasing() {
        let q: Vec<f32> = (0..100).map(|i| i as f32 / 100.0).collect();
        let t = vec![true; 100];
        let curve = qvalue_curves(&q, &t, 20, &[1.0]).remove(0);
        assert_eq!(curve.len(), 20);
        assert_eq!(curve.last().unwrap().1, 100.0);

        // Decoys alone contribute nothing: the curve counts targets.
        let decoys_only = qvalue_curves(&[0.5, 0.6], &[false, false], 5, &[1.0]).remove(0);
        assert!(decoys_only.iter().all(|p| p.1 == 0.0));
    }

    /// A zoomed curve must spend all of its points inside its own range and
    /// agree with the wide curve wherever they overlap -- a zoom that merely
    /// rescaled the axis without re-gridding would be the bug this guards.
    #[test]
    fn zoomed_curves_resolve_the_low_q_region() {
        // 1000 targets, all under q = 0.02: the whole population lands in the
        // first two points of the (0, 1] curve and nowhere else.
        let q: Vec<f32> = (0..1000).map(|i| i as f32 * 2e-5).collect();
        let t = vec![true; 1000];
        // 100 points puts q = 0.01 exactly on all three grids, so the three
        // views can be compared at a shared threshold below.
        let curves = qvalue_curves(&q, &t, 100, &[1.0, 0.05, 0.01]);

        let wide = &curves[0];
        assert!(
            wide.iter().filter(|p| p.1 > 0.0 && p.1 < 1000.0).count() <= 1,
            "the wide view should be saturated, that is why zooms exist"
        );

        for (curve, zoom) in curves.iter().zip([1.0, 0.05, 0.01]) {
            assert_eq!(curve.len(), 100);
            assert!(
                curve.iter().all(|&(x, _)| x > 0.0 && x <= zoom + 1e-12),
                "points must stay inside the zoom range"
            );
            assert!(
                (curve.last().unwrap().0 - zoom).abs() < 1e-12,
                "the last point must sit on the range's edge"
            );
        }

        // q = 0.01 is a grid point of every zoom here, so all three must
        // report the same count there.
        let at = |c: &Vec<(f64, f64)>| {
            c.iter()
                .find(|(x, _)| (x - 0.01).abs() < 1e-9)
                .map(|&(_, y)| y)
                .expect("0.01 is on every grid")
        };
        assert_eq!(at(&curves[0]), at(&curves[1]));
        assert_eq!(at(&curves[1]), at(&curves[2]));

        // And the zoom actually resolves: the 0.01 view climbs gradually
        // instead of jumping from nothing to everything in one step.
        assert!(
            curves[2]
                .iter()
                .filter(|p| p.1 > 0.0 && p.1 < 500.0)
                .count()
                > 10,
            "zoomed curve should have many intermediate points"
        );
    }

    #[test]
    fn pp_curve_is_diagonal_when_classes_are_identical() {
        let score: Vec<f32> = (0..200).map(|i| (i / 2) as f32).collect();
        let t: Vec<bool> = (0..200).map(|i| i % 2 == 0).collect();
        for (d, tt) in pp_curve(&score, &t, 25) {
            assert!((d - tt).abs() < 0.05, "decoy {d} vs target {tt}");
        }
    }

    /// Axes are `(decoy CDF, target CDF)`; the panel draws the `y = x`
    /// reference line, so a curve point where `d > tt` sits BELOW that
    /// diagonal. When targets score higher, the decoy CDF races ahead of the
    /// target CDF in the shared middle range -- the curve dips below y = x,
    /// which is what good target/decoy separation looks like on this plot.
    #[test]
    fn pp_curve_dips_below_the_diagonal_when_targets_score_higher() {
        let mut score = Vec::new();
        let mut t = Vec::new();
        for i in 0..100 {
            score.push(i as f32);
            t.push(false);
        }
        for i in 100..200 {
            score.push(i as f32);
            t.push(true);
        }
        let curve = pp_curve(&score, &t, 25);
        // Somewhere in the middle the decoy CDF is far ahead of the target CDF
        // (d > tt), i.e. the curve point sits well below the diagonal.
        assert!(
            curve.iter().any(|(d, tt)| d - tt > 0.4),
            "expected separation below the diagonal, got {curve:?}"
        );
    }

    #[test]
    fn pp_curve_is_empty_without_both_classes() {
        assert!(pp_curve(&[1.0, 2.0], &[true, true], 10).is_empty());
    }

    #[test]
    fn pp_curve_handles_identical_scores_without_panicking() {
        let score = vec![5.0f32; 10];
        let t: Vec<bool> = (0..10).map(|i| i % 2 == 0).collect();
        let curve = pp_curve(&score, &t, 10);
        assert!(
            curve.iter().all(|&(d, tt)| d == 1.0 && tt == 1.0),
            "degenerate range: every threshold already covers both classes, got {curve:?}"
        );
    }
}
