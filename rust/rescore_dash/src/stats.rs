//! Per-column statistics for the histogram and feature panels.
//!
//! Table statistics come from a fixed-width histogram: one pass to find the
//! range, one to bin, both O(rows) with no per-row allocation. Exact
//! rank-based AUC is O(rows log rows) and is used only for the selected
//! feature.
//!
//! Every function here treats non-finite values (NaN, +-Inf) as missing and
//! reports how many it saw, rather than propagating them into a mean.

/// Bin count for every histogram. Wide enough that a 1e6-row column keeps
/// shape, narrow enough to fit the AUC accumulation in cache.
pub const N_BINS: usize = 512;

#[derive(Debug, Clone, PartialEq)]
pub struct Hist {
    pub lo: f64,
    pub hi: f64,
    pub target: Vec<u32>,
    pub decoy: Vec<u32>,
    /// Non-finite values, excluded from the bins.
    pub n_nan: usize,
    /// Finite values outside `[lo, hi]`, excluded from the bins.
    pub n_out: usize,
}

impl Hist {
    /// Center of bin `i` on the value axis.
    pub fn bin_center(&self, i: usize) -> f64 {
        let w = (self.hi - self.lo) / self.target.len() as f64;
        self.lo + w * (i as f64 + 0.5)
    }

    pub fn n_target(&self) -> u32 {
        self.target.iter().sum()
    }

    pub fn n_decoy(&self) -> u32 {
        self.decoy.iter().sum()
    }
}

/// Bin `values` into `nbins` over `[lo, hi]`, split by label.
///
/// `is_target` is row-aligned with `values`; values beyond its length are
/// ignored. A degenerate range (`hi <= lo`) yields a single occupied bin.
pub fn histogram(
    values: impl Iterator<Item = f64>,
    is_target: &[bool],
    lo: f64,
    hi: f64,
    nbins: usize,
) -> Hist {
    let nbins = nbins.max(1);
    let mut h = Hist {
        lo,
        hi,
        target: vec![0; nbins],
        decoy: vec![0; nbins],
        n_nan: 0,
        n_out: 0,
    };
    let span = hi - lo;
    for (i, v) in values.enumerate() {
        let Some(&is_t) = is_target.get(i) else { break };
        if !v.is_finite() {
            h.n_nan += 1;
            continue;
        }
        if v < lo || v > hi {
            h.n_out += 1;
            continue;
        }
        let bin = if span > 0.0 {
            (((v - lo) / span) * nbins as f64) as usize
        } else {
            0
        };
        let bin = bin.min(nbins - 1);
        if is_t {
            h.target[bin] += 1;
        } else {
            h.decoy[bin] += 1;
        }
    }
    h
}

/// `(min, max)` over the finite values, or `None` if nothing is finite.
///
/// O(n), no allocation. Use this to pick a histogram range for a whole table
/// of columns — `summarize` calls it once per feature, and at ~1e6 rows x 131
/// features a sort per column (as `percentile_range` does) would mean 131
/// full sorts just to find the endpoints. `percentile_range` is for the
/// single selected-feature panel, where trimming outlier tails is wanted and
/// a sorted copy is affordable.
pub fn finite_range(values: impl Iterator<Item = f64>) -> Option<(f64, f64)> {
    values
        .filter(|v| v.is_finite())
        .fold(None, |acc, v| match acc {
            None => Some((v, v)),
            Some((lo, hi)) => Some((lo.min(v), hi.max(v))),
        })
}

/// `(lo, hi)` at the given percentiles over the finite values, or `None` if
/// nothing is finite. Allocates a sorted copy — call it once per panel, not
/// per frame.
pub fn percentile_range(
    values: impl Iterator<Item = f64>,
    lo_p: f64,
    hi_p: f64,
) -> Option<(f64, f64)> {
    let mut v: Vec<f64> = values.filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let pick = |p: f64| -> f64 {
        let idx = ((p / 100.0) * (v.len() - 1) as f64).round() as usize;
        v[idx.min(v.len() - 1)]
    };
    let (lo, hi) = (pick(lo_p), pick(hi_p));
    if hi > lo {
        Some((lo, hi))
    } else {
        Some((lo, lo + 1.0))
    }
}

/// AUC = P(target > decoy) + 0.5 * P(tie), read off the bins.
///
/// Ties within a bin are counted as half, which is also how the exact version
/// handles equal values — the two agree to bin resolution.
///
/// Returns NaN when either class is empty (no comparison is defined).
pub fn auc_from_hist(h: &Hist) -> f64 {
    let (nt, nd) = (h.n_target() as f64, h.n_decoy() as f64);
    if nt == 0.0 || nd == 0.0 {
        return f64::NAN;
    }
    let mut decoys_below = 0.0f64;
    let mut acc = 0.0f64;
    for i in 0..h.target.len() {
        let d = h.decoy[i] as f64;
        acc += h.target[i] as f64 * (decoys_below + 0.5 * d);
        decoys_below += d;
    }
    acc / (nt * nd)
}

/// Rank-based (Mann-Whitney) AUC over the finite values. Exact, including
/// ties, at the cost of a sort. NaN when either class is empty.
pub fn auc_exact(values: impl Iterator<Item = f64>, is_target: &[bool]) -> f64 {
    let mut pairs: Vec<(f64, bool)> = values
        .enumerate()
        .filter_map(|(i, v)| {
            let &is_t = is_target.get(i)?;
            v.is_finite().then_some((v, is_t))
        })
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let nt = pairs.iter().filter(|p| p.1).count() as f64;
    let nd = pairs.len() as f64 - nt;
    if nt == 0.0 || nd == 0.0 {
        return f64::NAN;
    }

    // Mid-ranks so tied values contribute 0.5 each.
    let mut rank_sum_targets = 0.0f64;
    let mut i = 0usize;
    while i < pairs.len() {
        let mut j = i;
        while j + 1 < pairs.len() && pairs[j + 1].0 == pairs[i].0 {
            j += 1;
        }
        let mid_rank = (i + j) as f64 / 2.0 + 1.0;
        for p in &pairs[i..=j] {
            if p.1 {
                rank_sum_targets += mid_rank;
            }
        }
        i = j + 1;
    }
    (rank_sum_targets - nt * (nt + 1.0) / 2.0) / (nt * nd)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeatureSummary {
    pub target_mean: f64,
    pub decoy_mean: f64,
    /// Histogram-resolution AUC.
    pub auc: f64,
    /// |Cohen's d| — pooled-SD standardized mean difference.
    pub cohens_d: f64,
    /// Fraction of rows whose value was non-finite.
    pub nan_frac: f64,
}

/// One column's table row. Two passes, both O(rows): exact per-class
/// mean/variance via Welford's algorithm, then a histogram-binned AUC over
/// the full finite value range (min/max via `finite_range`, not a sort).
pub fn summarize(values: impl Iterator<Item = f64> + Clone, is_target: &[bool]) -> FeatureSummary {
    // Welford per class.
    let (mut nt, mut nd) = (0u64, 0u64);
    let (mut mt, mut md) = (0.0f64, 0.0f64);
    let (mut m2t, mut m2d) = (0.0f64, 0.0f64);
    let mut n_nan = 0usize;
    let mut n_seen = 0usize;

    for (i, v) in values.clone().enumerate() {
        let Some(&is_t) = is_target.get(i) else { break };
        n_seen += 1;
        if !v.is_finite() {
            n_nan += 1;
            continue;
        }
        if is_t {
            nt += 1;
            let delta = v - mt;
            mt += delta / nt as f64;
            m2t += delta * (v - mt);
        } else {
            nd += 1;
            let delta = v - md;
            md += delta / nd as f64;
            m2d += delta * (v - md);
        }
    }

    let nan_frac = if n_seen == 0 {
        0.0
    } else {
        n_nan as f64 / n_seen as f64
    };
    if nt == 0 || nd == 0 {
        return FeatureSummary {
            target_mean: if nt == 0 { f64::NAN } else { mt },
            decoy_mean: if nd == 0 { f64::NAN } else { md },
            auc: f64::NAN,
            cohens_d: f64::NAN,
            nan_frac,
        };
    }

    let vt = if nt > 1 { m2t / (nt - 1) as f64 } else { 0.0 };
    let vd = if nd > 1 { m2d / (nd - 1) as f64 } else { 0.0 };
    let pooled = ((vt + vd) / 2.0).sqrt();
    let diff = mt - md;
    // Zero pooled SD means a constant column: no difference to standardize.
    let cohens_d = if pooled > 0.0 {
        (diff / pooled).abs()
    } else {
        0.0
    };

    let auc = match finite_range(values.clone()) {
        Some((lo, hi)) => auc_from_hist(&histogram(values, is_target, lo, hi, N_BINS)),
        None => f64::NAN,
    };

    FeatureSummary {
        target_mean: mt,
        decoy_mean: md,
        auc,
        cohens_d,
        nan_frac,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two classes that do not overlap: AUC must be exactly 1.
    fn separated() -> (Vec<f64>, Vec<bool>) {
        let mut v = Vec::new();
        let mut t = Vec::new();
        for i in 0..100 {
            v.push(i as f64);
            t.push(false);
        }
        for i in 100..200 {
            v.push(i as f64);
            t.push(true);
        }
        (v, t)
    }

    #[test]
    fn histogram_counts_per_class_and_clamps_outside() {
        let v = vec![0.0, 1.0, 2.0, 9.0];
        let t = vec![true, true, false, false];
        let h = histogram(v.into_iter(), &t, 0.0, 3.0, 3);
        assert_eq!(h.target.iter().sum::<u32>(), 2);
        assert_eq!(h.decoy.iter().sum::<u32>(), 1);
        assert_eq!(h.n_out, 1, "9.0 is outside [0, 3]");
        assert_eq!(h.n_nan, 0);
    }

    #[test]
    fn histogram_counts_nans_separately() {
        let v = vec![1.0, f64::NAN, f64::INFINITY];
        let t = vec![true, true, false];
        let h = histogram(v.into_iter(), &t, 0.0, 2.0, 4);
        assert_eq!(h.n_nan, 2, "NaN and non-finite both count as missing");
        assert_eq!(h.target.iter().sum::<u32>(), 1);
    }

    #[test]
    fn auc_is_one_for_separated_classes() {
        let (v, t) = separated();
        assert!((auc_exact(v.iter().copied(), &t) - 1.0).abs() < 1e-9);
        let h = histogram(v.into_iter(), &t, 0.0, 199.0, N_BINS);
        assert!((auc_from_hist(&h) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn auc_is_half_for_identical_classes() {
        let v = vec![1.0; 100];
        let t: Vec<bool> = (0..100).map(|i| i % 2 == 0).collect();
        assert!(
            (auc_exact(v.iter().copied(), &t) - 0.5).abs() < 1e-9,
            "all ties"
        );
    }

    #[test]
    fn auc_is_zero_when_decoys_score_higher() {
        let (v, mut t) = separated();
        for b in &mut t {
            *b = !*b;
        }
        assert!(auc_exact(v.into_iter(), &t) < 1e-9);
    }

    #[test]
    fn summarize_reports_means_and_effect_size() {
        let (v, t) = separated();
        let s = summarize(v.into_iter(), &t);
        assert!((s.target_mean - 149.5).abs() < 1e-6);
        assert!((s.decoy_mean - 49.5).abs() < 1e-6);
        assert!(s.cohens_d > 3.0, "well-separated, got {}", s.cohens_d);
        assert_eq!(s.nan_frac, 0.0);
    }

    #[test]
    fn summarize_survives_an_all_nan_feature() {
        let v = vec![f64::NAN; 10];
        let t: Vec<bool> = (0..10).map(|i| i % 2 == 0).collect();
        let s = summarize(v.into_iter(), &t);
        assert_eq!(s.nan_frac, 1.0);
        assert!(s.target_mean.is_nan());
        assert!(s.auc.is_nan());
        assert!(s.cohens_d.is_nan());
    }

    #[test]
    fn summarize_survives_a_constant_feature() {
        let v = vec![3.0; 10];
        let t: Vec<bool> = (0..10).map(|i| i % 2 == 0).collect();
        let s = summarize(v.into_iter(), &t);
        assert_eq!(s.target_mean, 3.0);
        assert_eq!(s.decoy_mean, 3.0);
        assert_eq!(s.cohens_d, 0.0, "zero variance, zero difference");
        assert!((s.auc - 0.5).abs() < 1e-9);
    }

    #[test]
    fn summarize_survives_an_empty_decoy_set() {
        let v = vec![1.0, 2.0, 3.0];
        let t = vec![true, true, true];
        let s = summarize(v.into_iter(), &t);
        assert!(s.decoy_mean.is_nan());
        assert!(s.auc.is_nan());
    }

    #[test]
    fn percentile_range_clips_the_tails() {
        let v: Vec<f64> = (0..=1000).map(|i| i as f64).collect();
        let (lo, hi) = percentile_range(v.into_iter(), 0.5, 99.5).unwrap();
        assert!((1.0..=10.0).contains(&lo), "lo = {lo}");
        assert!((990.0..=999.0).contains(&hi), "hi = {hi}");
    }

    #[test]
    fn percentile_range_is_none_when_everything_is_missing() {
        let v = vec![f64::NAN, f64::NEG_INFINITY];
        assert!(percentile_range(v.into_iter(), 0.5, 99.5).is_none());
    }

    #[test]
    fn finite_range_spans_the_finite_values() {
        let v = vec![3.0, -1.0, f64::NAN, 7.0, f64::INFINITY];
        assert_eq!(finite_range(v.into_iter()), Some((-1.0, 7.0)));
    }

    #[test]
    fn finite_range_is_none_when_everything_is_missing() {
        let v = vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
        assert!(finite_range(v.into_iter()).is_none());
    }

    #[test]
    fn finite_range_of_a_single_value_is_that_value_twice() {
        let v = vec![42.0];
        assert_eq!(finite_range(v.into_iter()), Some((42.0, 42.0)));
    }

    #[test]
    fn histogram_degenerate_range_puts_everything_in_one_bin() {
        let v = vec![5.0, 5.0, 5.0];
        let t = vec![true, false, true];
        let h = histogram(v.into_iter(), &t, 5.0, 5.0, 8);
        assert_eq!(h.target.len(), 8);
        assert_eq!(h.decoy.len(), 8);
        assert_eq!(h.target[0], 2);
        assert_eq!(h.decoy[0], 1);
        assert_eq!(
            h.n_out, 0,
            "every value equals lo == hi, so none is out of range"
        );
    }

    #[test]
    fn histogram_zero_nbins_is_clamped_to_one_bin() {
        let v = vec![1.0, 2.0, 3.0];
        let t = vec![true, true, false];
        let h = histogram(v.into_iter(), &t, 0.0, 3.0, 0);
        assert_eq!(h.target.len(), 1, "nbins = 0 is clamped up to 1");
        assert_eq!(h.decoy.len(), 1);
        assert_eq!(h.target.iter().sum::<u32>(), 2);
        assert_eq!(h.decoy.iter().sum::<u32>(), 1);
    }
}
