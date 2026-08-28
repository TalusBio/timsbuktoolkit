//! Per-column statistics and rank statistics.
//!
//! Everything here scans data and runs at init.
//!
//! Non-finite values (NaN, +-Inf) are missing everywhere: they are counted, not
//! propagated into a mean, a range or a rank.

/// Exact, all-rows statistics for one column.
///
/// Doubles as the accumulator for the sweep that produces it: filled by
/// [`ColumnStats::push`], combined by [`ColumnStats::merge`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ColumnStats {
    pub(crate) n_target: u64,
    pub(crate) n_decoy: u64,
    /// Class means, read through [`Self::mean`].
    target_mean: f64,
    decoy_mean: f64,
    /// Welford's sum of squared deviations, per class. Read through
    /// [`Self::var`].
    m2_target: f64,
    m2_decoy: f64,
    /// Finite min/max over both classes. `INFINITY` / `NEG_INFINITY` when the
    /// column holds no finite value at all, which is the sentinel every
    /// consumer checks with `hi >= lo`.
    pub(crate) lo: f64,
    pub(crate) hi: f64,
    /// Exact lower bound of `log10`'s domain over all rows, or `INFINITY` if
    /// the column is wholly non-positive.
    pub(crate) min_pos: f64,
    /// Smallest `|v|` over the finite values, or `INFINITY` if there are none.
    /// The exact lower bound of `square`'s output over all rows, which is not
    /// `min(|lo|, |hi|)` for a column that straddles zero.
    pub(crate) min_abs: f64,
    /// Exact zeros, which `min_pos` deliberately skips.
    pub(crate) n_zero: u64,
    pub(crate) n_nan: u64,
}

impl ColumnStats {
    pub(crate) const IDENTITY: Self = Self {
        n_target: 0,
        n_decoy: 0,
        target_mean: 0.0,
        decoy_mean: 0.0,
        m2_target: 0.0,
        m2_decoy: 0.0,
        lo: f64::INFINITY,
        hi: f64::NEG_INFINITY,
        min_pos: f64::INFINITY,
        min_abs: f64::INFINITY,
        n_zero: 0,
        n_nan: 0,
    };

    /// Fold one value in.
    #[inline]
    pub(crate) fn push(&mut self, v: f64, is_target: bool) {
        if !v.is_finite() {
            self.n_nan += 1;
            return;
        }
        self.lo = self.lo.min(v);
        self.hi = self.hi.max(v);
        self.min_abs = self.min_abs.min(v.abs());
        self.n_zero += u64::from(v == 0.0);
        let positive = if v > 0.0 { v } else { f64::INFINITY };
        self.min_pos = self.min_pos.min(positive);
        if is_target {
            self.n_target += 1;
            let d = v - self.target_mean;
            self.target_mean += d / self.n_target as f64;
            self.m2_target += d * (v - self.target_mean);
        } else {
            self.n_decoy += 1;
            let d = v - self.decoy_mean;
            self.decoy_mean += d / self.n_decoy as f64;
            self.m2_decoy += d * (v - self.decoy_mean);
        }
    }

    /// Combine a partial sweep into this one. Chan et al.'s parallel variance
    /// combination, per class: plain summation of `m2` would drop the term for
    /// the shift between the two partial means.
    pub(crate) fn merge(&mut self, o: &Self) {
        merge_class(
            &mut self.n_target,
            &mut self.target_mean,
            &mut self.m2_target,
            o.n_target,
            o.target_mean,
            o.m2_target,
        );
        merge_class(
            &mut self.n_decoy,
            &mut self.decoy_mean,
            &mut self.m2_decoy,
            o.n_decoy,
            o.decoy_mean,
            o.m2_decoy,
        );
        self.lo = self.lo.min(o.lo);
        self.hi = self.hi.max(o.hi);
        self.min_pos = self.min_pos.min(o.min_pos);
        self.min_abs = self.min_abs.min(o.min_abs);
        self.n_zero += o.n_zero;
        self.n_nan += o.n_nan;
    }

    /// Sample variance of one class, or `0.0` for a class with fewer than two
    /// values.
    pub(crate) fn var(&self, is_target: bool) -> f64 {
        let (n, m2) = if is_target {
            (self.n_target, self.m2_target)
        } else {
            (self.n_decoy, self.m2_decoy)
        };
        if n > 1 { m2 / (n - 1) as f64 } else { 0.0 }
    }

    /// Whether any finite value was seen.
    pub(crate) fn has_finite(&self) -> bool {
        self.hi >= self.lo
    }

    /// `|Cohen's d|` -- pooled-SD standardized mean difference. NaN when either
    /// class contributed no finite value; `0.0` for a constant column, where
    /// there is no difference to standardize.
    pub(crate) fn cohens_d(&self) -> f64 {
        if self.n_target == 0 || self.n_decoy == 0 {
            return f64::NAN;
        }
        let pooled = ((self.var(true) + self.var(false)) / 2.0).sqrt();
        if pooled > 0.0 {
            ((self.target_mean - self.decoy_mean) / pooled).abs()
        } else {
            0.0
        }
    }

    /// Fraction of rows whose value was non-finite.
    pub(crate) fn nan_frac(&self) -> f64 {
        let seen = self.n_target + self.n_decoy + self.n_nan;
        if seen == 0 {
            0.0
        } else {
            self.n_nan as f64 / seen as f64
        }
    }

    /// Class mean, or NaN when that class contributed no finite value.
    pub(crate) fn mean(&self, is_target: bool) -> f64 {
        let (n, m) = if is_target {
            (self.n_target, self.target_mean)
        } else {
            (self.n_decoy, self.decoy_mean)
        };
        if n == 0 { f64::NAN } else { m }
    }
}

fn merge_class(n_a: &mut u64, mean_a: &mut f64, m2_a: &mut f64, n_b: u64, mean_b: f64, m2_b: f64) {
    if n_b == 0 {
        return;
    }
    if *n_a == 0 {
        (*n_a, *mean_a, *m2_a) = (n_b, mean_b, m2_b);
        return;
    }
    let (na, nb) = (*n_a as f64, n_b as f64);
    let delta = mean_b - *mean_a;
    let n = na + nb;
    *mean_a += delta * nb / n;
    *m2_a += m2_b + delta * delta * na * nb / n;
    *n_a += n_b;
}

/// Runs of equal values in an ascending slice, as inclusive `(first, last)`
/// index pairs.
///
/// Both consumers below are tie-aware and this is the scan they share: a run's
/// members get its mid-rank, whether that is being turned into a percentile or
/// summed into a rank sum.
fn tie_runs(sorted: &[(f64, bool)]) -> impl Iterator<Item = (usize, usize)> + '_ {
    let mut i = 0usize;
    std::iter::from_fn(move || {
        let first = i;
        if first >= sorted.len() {
            return None;
        }
        while i + 1 < sorted.len() && sorted[i + 1].0 == sorted[first].0 {
            i += 1;
        }
        let last = i;
        i += 1;
        Some((first, last))
    })
}

/// Mid-rank percentiles (`0..=100`) for an ascending slice, written into `out`
/// at the same positions.
///
/// Ties share the mid-rank of their run, so a wholly constant column maps
/// entirely to 50. This is the `RankPercentile` transform: because the input is
/// already sorted, the transform is a single linear pass rather than the sort
/// it costs on unordered data.
pub(crate) fn mid_rank_percentiles(sorted: &[(f64, bool)], out: &mut Vec<f64>) {
    out.clear();
    out.resize(sorted.len(), 0.0);
    if sorted.len() < 2 {
        out.fill(50.0);
        return;
    }
    let denom = (sorted.len() - 1) as f64;
    for (i, j) in tie_runs(sorted) {
        out[i..=j].fill(100.0 * ((i + j) as f64 / 2.0) / denom);
    }
}

/// Mann-Whitney AUC = P(target > decoy) + 0.5 * P(tie) over pairs already
/// sorted ascending by value. Exact, including ties. NaN when either class is
/// empty, since no comparison is defined.
pub(crate) fn auc_from_sorted(sorted: &[(f64, bool)]) -> f64 {
    let nt = sorted.iter().filter(|p| p.1).count() as f64;
    let nd = sorted.len() as f64 - nt;
    if nt == 0.0 || nd == 0.0 {
        return f64::NAN;
    }
    let mut rank_sum_targets = 0.0f64;
    for (i, j) in tie_runs(sorted) {
        let mid_rank = (i + j) as f64 / 2.0 + 1.0;
        rank_sum_targets += mid_rank * sorted[i..=j].iter().filter(|p| p.1).count() as f64;
    }
    (rank_sum_targets - nt * (nt + 1.0) / 2.0) / (nt * nd)
}

/// [`auc_from_sorted`] over an unsorted column, paying the sort. Used for the
/// discriminant score, which is one column and is therefore reported exactly
/// over all rows rather than from the sample.
pub(crate) fn auc_exact(values: impl Iterator<Item = f64>, is_target: &[bool]) -> f64 {
    let mut pairs: Vec<(f64, bool)> = values
        .zip(is_target)
        .filter(|(v, _)| v.is_finite())
        .map(|(v, &t)| (v, t))
        .collect();
    pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    auc_from_sorted(&pairs)
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

    fn sweep(values: &[f64], is_target: &[bool]) -> ColumnStats {
        let mut c = ColumnStats::IDENTITY;
        for (i, &v) in values.iter().enumerate() {
            c.push(v, is_target[i]);
        }
        c
    }

    #[test]
    fn stats_report_means_and_effect_size() {
        let (v, t) = separated();
        let c = sweep(&v, &t);
        assert!((c.mean(true) - 149.5).abs() < 1e-6);
        assert!((c.mean(false) - 49.5).abs() < 1e-6);
        assert!(c.cohens_d() > 3.0, "well separated, got {}", c.cohens_d());
        assert_eq!(c.nan_frac(), 0.0);
        assert_eq!((c.lo, c.hi), (0.0, 199.0));
    }

    #[test]
    fn stats_survive_an_all_nan_column() {
        let v = vec![f64::NAN; 10];
        let t: Vec<bool> = (0..10).map(|i| i % 2 == 0).collect();
        let c = sweep(&v, &t);
        assert_eq!(c.nan_frac(), 1.0);
        assert!(c.mean(true).is_nan());
        assert!(c.cohens_d().is_nan());
        assert!(!c.has_finite(), "the lo/hi sentinel must survive intact");
    }

    #[test]
    fn stats_survive_a_constant_column() {
        let v = vec![3.0; 10];
        let t: Vec<bool> = (0..10).map(|i| i % 2 == 0).collect();
        let c = sweep(&v, &t);
        assert_eq!(c.mean(true), 3.0);
        assert_eq!(c.mean(false), 3.0);
        assert_eq!(c.cohens_d(), 0.0, "zero variance, zero difference");
    }

    /// Neither floor is `min(|lo|, |hi|)` for a column that straddles zero.
    #[test]
    fn stats_track_the_positive_and_absolute_floors_across_zero() {
        let v = vec![-5.0, -0.25, 0.0, 0.5, 8.0];
        let t = vec![true, false, true, false, true];
        let c = sweep(&v, &t);
        assert_eq!(c.lo, -5.0);
        assert_eq!(c.hi, 8.0);
        assert_eq!(c.min_pos, 0.5, "zero and negatives are not positive");
        assert_eq!(c.min_abs, 0.0, "the exact zero has the smallest magnitude");
        assert_eq!(c.n_zero, 1);
        assert_ne!(
            c.min_abs,
            c.lo.abs().min(c.hi.abs()),
            "min_abs must not be derivable from the endpoints here"
        );

        // A wholly non-positive column has no `log10` domain at all, which is
        // the `INFINITY` sentinel `exact_range` checks for.
        let c = sweep(&[-3.0, -1.0, -7.0], &[true, false, true]);
        assert_eq!(c.min_pos, f64::INFINITY, "no positive value exists");
        assert_eq!(c.min_abs, 1.0);
        assert_eq!(c.n_zero, 0);
    }

    /// The parallel sweep splits rows across threads and merges. Chan's
    /// combination has to reproduce the sequential variance exactly, not
    /// approximately, or `cohens_d` drifts with the thread count.
    #[test]
    fn merging_partial_sweeps_matches_one_sequential_sweep() {
        let v: Vec<f64> = (0..1000).map(|i| (i as f64 * 0.37).sin() * 100.0).collect();
        let t: Vec<bool> = (0..1000).map(|i| i % 3 == 0).collect();
        let whole = sweep(&v, &t);

        let mut merged = ColumnStats::IDENTITY;
        for chunk in 0..7 {
            let (a, b) = (chunk * 1000 / 7, (chunk + 1) * 1000 / 7);
            merged.merge(&sweep(&v[a..b], &t[a..b]));
        }
        assert_eq!(merged.n_target, whole.n_target);
        assert!((merged.mean(true) - whole.mean(true)).abs() < 1e-9);
        assert!((merged.var(true) - whole.var(true)).abs() < 1e-9);
        assert!((merged.var(false) - whole.var(false)).abs() < 1e-9);
        assert_eq!((merged.lo, merged.hi), (whole.lo, whole.hi));
        assert_eq!(merged.min_pos, whole.min_pos);
        assert_eq!(merged.min_abs, whole.min_abs);

        // IDENTITY's `INFINITY`/`NEG_INFINITY` sentinels must be neutral under
        // the min/max folds, not poison them.
        let mut with_empty = whole;
        with_empty.merge(&ColumnStats::IDENTITY);
        assert_eq!(with_empty, whole);
    }

    fn sorted_pairs(values: &[f64], is_target: &[bool]) -> Vec<(f64, bool)> {
        let mut p: Vec<(f64, bool)> = values
            .iter()
            .zip(is_target)
            .filter(|(v, _)| v.is_finite())
            .map(|(&v, &t)| (v, t))
            .collect();
        p.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        p
    }

    /// The three fixed points of the statistic: perfectly separated, reversed,
    /// and wholly tied.
    #[test]
    fn auc_hits_one_zero_and_a_half() {
        let (v, t) = separated();
        assert!((auc_from_sorted(&sorted_pairs(&v, &t)) - 1.0).abs() < 1e-9);

        let flipped: Vec<bool> = t.iter().map(|b| !b).collect();
        assert!(auc_from_sorted(&sorted_pairs(&v, &flipped)) < 1e-9);

        let tied = vec![1.0; 100];
        let alternating: Vec<bool> = (0..100).map(|i| i % 2 == 0).collect();
        assert!((auc_from_sorted(&sorted_pairs(&tied, &alternating)) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn auc_is_nan_without_both_classes() {
        let v = vec![1.0, 2.0, 3.0];
        let t = vec![true, true, true];
        assert!(auc_from_sorted(&sorted_pairs(&v, &t)).is_nan());
        assert!(auc_exact(v.into_iter(), &t).is_nan());
    }

    #[test]
    fn mid_rank_percentiles_span_zero_to_one_hundred() {
        let s = sorted_pairs(&[10.0, 30.0, 50.0], &[true, true, true]);
        let mut pct = Vec::new();
        mid_rank_percentiles(&s, &mut pct);
        assert_eq!(pct[0], 0.0);
        assert!((pct[1] - 50.0).abs() < 1e-9);
        assert_eq!(pct[2], 100.0);
    }

    /// Ties share the mid-rank of their run, so an all-tied column collapses to
    /// 50 regardless of length -- not to 0, and not spread across the axis.
    #[test]
    fn mid_rank_percentiles_map_an_all_tied_column_entirely_to_fifty() {
        for len in [1, 2, 3, 4, 10, 11] {
            let s = sorted_pairs(&vec![7.0; len], &vec![true; len]);
            let mut pct = Vec::new();
            mid_rank_percentiles(&s, &mut pct);
            assert_eq!(pct.len(), len);
            assert!(
                pct.iter().all(|&v| (v - 50.0).abs() < 1e-9),
                "len {len}: got {pct:?}"
            );
        }
    }

    #[test]
    fn mid_rank_percentiles_handle_an_empty_column() {
        let mut pct = vec![1.0, 2.0];
        mid_rank_percentiles(&[], &mut pct);
        assert!(pct.is_empty(), "stale contents must be cleared");
    }
}
