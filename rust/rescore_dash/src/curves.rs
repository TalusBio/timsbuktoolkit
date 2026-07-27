//! The FDR and calibration curves.
//!
//! Both are cumulative counts over a grid, built once in `Dashboard::build`.

/// How many rows pass at one q-value cutoff.
///
/// `total == targets + decoys` always: an unlabeled row is not counted at all
/// rather than counted into the total alone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThresholdRow {
    pub q: f32,
    pub total: usize,
    pub targets: usize,
    pub decoys: usize,
}

/// [`ThresholdRow`] for each threshold, in the order given.
///
/// One pass over the rows for all thresholds, not one pass each: every row is
/// counted into each threshold it clears, so the cost is independent of how
/// many cutoffs are asked for.
///
/// Deliberately a re-derivation of the same quantity as
/// `timsseek::ml::qvalues::report_qvalues_at_thresholds` rather than a call to
/// it — this crate takes no pipeline dependency (see [`crate::view`]). The two
/// agree because they are both "rows at or below the cutoff, by class".
///
/// `is_target` is row-aligned with `qvalue`; a row with no label (past the end
/// of `is_target`) is skipped entirely.
pub fn threshold_table(
    qvalue: &[f32],
    is_target: &[bool],
    thresholds: &[f32],
) -> Vec<ThresholdRow> {
    let mut counts = vec![(0usize, 0usize); thresholds.len()];
    for (i, &row_q) in qvalue.iter().enumerate() {
        let Some(&is_t) = is_target.get(i) else {
            continue;
        };
        if row_q.is_nan() {
            continue;
        }
        for (&q, c) in thresholds.iter().zip(counts.iter_mut()) {
            if row_q <= q {
                if is_t {
                    c.0 += 1;
                } else {
                    c.1 += 1;
                }
            }
        }
    }
    thresholds
        .iter()
        .zip(counts)
        .map(|(&q, (targets, decoys))| ThresholdRow {
            q,
            total: targets + decoys,
            targets,
            decoys,
        })
        .collect()
}

/// Targets passing as a function of the q-value threshold: one curve per entry
/// of `zooms`, curve `i` spanning `n_points` evenly spaced thresholds in
/// `(0, zooms[i]]`.
///
/// Several curves rather than one because the full `(0, 1]` view answers the
/// least: on any run worth looking at, essentially every target has arrived by
/// q = 0.05, so the shape that matters — how many IDs a tighter cutoff costs —
/// is squeezed into the leftmost few percent of the panel and reads as a
/// vertical line. Sub-sampling a `(0, 1]` curve down to `q <= 0.01` would leave
/// one point in a hundred, so each zoom needs its own grid.
///
/// They are computed together because the sort is the expensive part and it is
/// shared: every zoom is answered by `partition_point` over the same sorted
/// array, so the extra views cost a scan and no extra pass over the run.
pub fn qvalue_curves(
    qvalue: &[f32],
    is_target: &[bool],
    n_points: usize,
    zooms: &[f64],
) -> Vec<Vec<(f64, f64)>> {
    let n_points = n_points.max(2);
    let mut targets: Vec<f32> = qvalue
        .iter()
        .enumerate()
        .filter(|(i, q)| !q.is_nan() && is_target.get(*i).copied().unwrap_or(false))
        .map(|(_, q)| *q)
        .collect();
    targets.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    zooms
        .iter()
        .map(|&zoom| {
            (1..=n_points)
                .map(|k| {
                    let thresh = zoom * k as f64 / n_points as f64;
                    // Sorted: partition_point is the count at or below the
                    // threshold.
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
/// Empty when either class is absent — a PP plot of one class says nothing.
pub fn pp_curve(score: &[f32], is_target: &[bool], n_points: usize) -> Vec<(f64, f64)> {
    let n_points = n_points.max(2);
    let mut targets: Vec<f32> = Vec::new();
    let mut decoys: Vec<f32> = Vec::new();
    for (i, &s) in score.iter().enumerate() {
        if !s.is_finite() {
            continue;
        }
        match is_target.get(i) {
            Some(true) => targets.push(s),
            Some(false) => decoys.push(s),
            None => {}
        }
    }
    if targets.is_empty() || decoys.is_empty() {
        return Vec::new();
    }
    let cmp = |a: &f32, b: &f32| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal);
    targets.sort_by(cmp);
    decoys.sort_by(cmp);

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

    fn row(q: f32, total: usize, targets: usize, decoys: usize) -> ThresholdRow {
        ThresholdRow {
            q,
            total,
            targets,
            decoys,
        }
    }

    #[test]
    fn threshold_table_matches_hand_counts() {
        let q = vec![0.001f32, 0.02, 0.2, 0.9];
        let t = vec![true, true, false, false];
        let got = threshold_table(&q, &t, &[0.01, 0.05, 1.0]);
        assert_eq!(got[0], row(0.01, 1, 1, 0));
        assert_eq!(got[1], row(0.05, 2, 2, 0));
        assert_eq!(got[2], row(1.0, 4, 2, 2));
    }

    #[test]
    fn threshold_table_counts_labeled_rows_only() {
        // is_target is shorter than qvalue: rows 3 and 4 have no label and
        // must be skipped entirely, not counted into `total` alone.
        let q = vec![0.001f32, 0.02, 0.2, 0.9, 0.5];
        let t = vec![true, true, false];
        let got = threshold_table(&q, &t, &[0.01, 0.05, 1.0]);
        assert_eq!(got[2], row(1.0, 3, 2, 1), "unlabeled rows never counted");
    }

    #[test]
    fn qvalue_curve_is_monotone_nondecreasing() {
        let q: Vec<f32> = (0..100).map(|i| i as f32 / 100.0).collect();
        let t = vec![true; 100];
        let curve = qvalue_curves(&q, &t, 20, &[1.0]).remove(0);
        assert_eq!(curve.len(), 20);
        for w in curve.windows(2) {
            assert!(w[1].0 >= w[0].0, "q must ascend");
            assert!(
                w[1].1 >= w[0].1,
                "target count must not shrink as q loosens"
            );
        }
        assert_eq!(curve.last().unwrap().1, 100.0);
    }

    #[test]
    fn qvalue_curve_handles_no_targets() {
        let q = vec![0.5f32, 0.6];
        let t = vec![false, false];
        let curve = qvalue_curves(&q, &t, 5, &[1.0]).remove(0);
        assert!(curve.iter().all(|p| p.1 == 0.0));
    }

    /// A zoomed curve must spend all of its points inside its own range and
    /// agree with the wide curve wherever they overlap — a zoom that merely
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
    /// target CDF in the shared middle range — the curve dips below y = x,
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
