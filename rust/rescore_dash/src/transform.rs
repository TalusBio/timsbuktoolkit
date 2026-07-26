//! Axis transforms for the histogram panels.
//!
//! A transform never errors. Values it cannot map (non-positives under a log,
//! negatives under a square root, anything non-finite) are dropped, and the
//! caller reports the drop count in the panel title — a silently shrinking
//! histogram would misread as "these rows do not exist".

/// Largest exponent `exp` is allowed to see, guarding against overflow to
/// `+inf`. Only the upper side needs guarding: `exp` overflows above ~709.78,
/// but very negative inputs already underflow to exactly `0.0` on their own
/// well before that. Deliberately unclamped below — do not "restore the
/// symmetry" by adding a lower bound: `exp(-EXP_CLAMP)` is `9.86e-305`, a
/// manufactured nonzero value where the true answer is zero to any precision
/// the display cares about.
const EXP_CLAMP: f64 = 700.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XTransform {
    Linear,
    Log10,
    SignedLog1p,
    Sqrt,
    Square,
    Exp,
    RankPercentile,
}

impl XTransform {
    const CYCLE: [XTransform; 7] = [
        Self::Linear,
        Self::Log10,
        Self::SignedLog1p,
        Self::Sqrt,
        Self::Square,
        Self::Exp,
        Self::RankPercentile,
    ];

    pub fn next(self) -> Self {
        let i = Self::CYCLE.iter().position(|t| *t == self).unwrap_or(0);
        Self::CYCLE[(i + 1) % Self::CYCLE.len()]
    }

    pub fn prev(self) -> Self {
        let i = Self::CYCLE.iter().position(|t| *t == self).unwrap_or(0);
        Self::CYCLE[(i + Self::CYCLE.len() - 1) % Self::CYCLE.len()]
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Log10 => "log10",
            Self::SignedLog1p => "signed-log1p",
            Self::Sqrt => "sqrt",
            Self::Square => "square",
            Self::Exp => "exp",
            Self::RankPercentile => "rank-pct",
        }
    }

    /// Pointwise map. `None` means the value is outside this transform's
    /// domain and must be dropped. `RankPercentile` is not pointwise — it
    /// returns `Some(v)` here and is handled in [`transform_column`].
    pub fn apply(self, v: f64) -> Option<f64> {
        if !v.is_finite() {
            return None;
        }
        let out = match self {
            Self::Linear | Self::RankPercentile => v,
            Self::Log10 => {
                if v <= 0.0 {
                    return None;
                }
                v.log10()
            }
            Self::SignedLog1p => v.signum() * v.abs().ln_1p(),
            Self::Sqrt => {
                if v < 0.0 {
                    return None;
                }
                v.sqrt()
            }
            Self::Square => v * v,
            Self::Exp => v.min(EXP_CLAMP).exp(),
        };
        out.is_finite().then_some(out)
    }
}

/// Transform a whole column, returning the surviving values (input order) and
/// the number dropped.
pub fn transform_column(t: XTransform, values: &[f64]) -> (Vec<f64>, usize) {
    if t == XTransform::RankPercentile {
        return rank_percentile(values);
    }
    let mut out = Vec::with_capacity(values.len());
    let mut dropped = 0;
    for &v in values {
        match t.apply(v) {
            Some(x) => out.push(x),
            None => dropped += 1,
        }
    }
    (out, dropped)
}

/// Map finite values onto 0..=100 by rank, preserving input order. Ties share
/// the mid-rank, so a constant column maps entirely to 50.
fn rank_percentile(values: &[f64]) -> (Vec<f64>, usize) {
    let mut idx: Vec<usize> = (0..values.len())
        .filter(|&i| values[i].is_finite())
        .collect();
    let dropped = values.len() - idx.len();
    if idx.is_empty() {
        return (Vec::new(), dropped);
    }
    if idx.len() == 1 {
        return (vec![50.0], dropped);
    }
    idx.sort_by(|&a, &b| {
        values[a]
            .partial_cmp(&values[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let denom = (idx.len() - 1) as f64;
    let mut pct = vec![0.0f64; values.len()];
    let mut i = 0usize;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && values[idx[j + 1]] == values[idx[i]] {
            j += 1;
        }
        let mid = (i + j) as f64 / 2.0;
        let p = 100.0 * mid / denom;
        for &k in &idx[i..=j] {
            pct[k] = p;
        }
        i = j + 1;
    }

    let out = (0..values.len())
        .filter(|&i| values[i].is_finite())
        .map(|i| pct[i])
        .collect();
    (out, dropped)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum YTransform {
    /// Fraction of the class's own total. The default: target and decoy counts
    /// differ by orders of magnitude, and raw counts flatten the decoy shape
    /// into the axis.
    Density,
    Count,
    Log10Count,
}

impl YTransform {
    const CYCLE: [YTransform; 3] = [Self::Density, Self::Count, Self::Log10Count];

    pub fn next(self) -> Self {
        let i = Self::CYCLE.iter().position(|t| *t == self).unwrap_or(0);
        Self::CYCLE[(i + 1) % Self::CYCLE.len()]
    }

    pub fn prev(self) -> Self {
        let i = Self::CYCLE.iter().position(|t| *t == self).unwrap_or(0);
        Self::CYCLE[(i + Self::CYCLE.len() - 1) % Self::CYCLE.len()]
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Density => "density",
            Self::Count => "count",
            Self::Log10Count => "log10 count",
        }
    }

    pub fn apply(self, count: u32, class_total: u32) -> f64 {
        match self {
            Self::Density => {
                if class_total == 0 {
                    0.0
                } else {
                    count as f64 / class_total as f64
                }
            }
            Self::Count => count as f64,
            Self::Log10Count => (count as f64 + 1.0).log10(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x_cycles_forward_and_back() {
        let t = XTransform::Linear;
        assert_eq!(t.next().prev(), t);
        // The cycle covers every variant and returns to the start.
        let mut seen = vec![t];
        let mut cur = t;
        for _ in 0..6 {
            cur = cur.next();
            seen.push(cur);
        }
        assert_eq!(cur.next(), XTransform::Linear);
        assert_eq!(seen.len(), 7, "seven transforms in the cycle");
    }

    #[test]
    fn log10_drops_non_positives() {
        let (out, dropped) = transform_column(XTransform::Log10, &[100.0, 0.0, -1.0, 10.0]);
        assert_eq!(dropped, 2);
        assert_eq!(out, vec![2.0, 1.0]);
    }

    #[test]
    fn sqrt_drops_negatives_and_keeps_zero() {
        let (out, dropped) = transform_column(XTransform::Sqrt, &[4.0, 0.0, -9.0]);
        assert_eq!(dropped, 1);
        assert_eq!(out, vec![2.0, 0.0]);
    }

    #[test]
    fn signed_log1p_is_total_and_odd() {
        let (out, dropped) = transform_column(XTransform::SignedLog1p, &[-9.0, 0.0, 9.0]);
        assert_eq!(dropped, 0);
        assert!((out[0] + out[2]).abs() < 1e-12, "odd about zero");
        assert_eq!(out[1], 0.0);
    }

    #[test]
    fn square_is_total() {
        let (out, dropped) = transform_column(XTransform::Square, &[-3.0, 0.0, 2.0]);
        assert_eq!(dropped, 0);
        assert_eq!(out, vec![9.0, 0.0, 4.0]);
    }

    #[test]
    fn exp_clamps_instead_of_overflowing() {
        let (out, dropped) = transform_column(XTransform::Exp, &[0.0, 1e9, -1e9]);
        assert_eq!(dropped, 0, "clamped, not dropped");
        assert_eq!(out[0], 1.0);
        assert!(out[1].is_finite(), "1e9 must clamp, got {}", out[1]);
        assert_eq!(out[2], 0.0);
    }

    #[test]
    fn non_finite_input_is_always_dropped() {
        for t in [
            XTransform::Linear,
            XTransform::Square,
            XTransform::SignedLog1p,
            XTransform::RankPercentile,
        ] {
            let (out, dropped) = transform_column(t, &[1.0, f64::NAN, f64::INFINITY]);
            assert_eq!(dropped, 2, "{t:?} must drop non-finite values");
            assert_eq!(out.len(), 1);
        }
    }

    #[test]
    fn rank_percentile_maps_onto_zero_to_one_hundred() {
        let (out, dropped) = transform_column(XTransform::RankPercentile, &[50.0, 10.0, 30.0]);
        assert_eq!(dropped, 0);
        // Output stays in input order: 50 is the largest, 10 the smallest.
        assert_eq!(out[1], 0.0);
        assert_eq!(out[0], 100.0);
        assert!((out[2] - 50.0).abs() < 1e-9);
    }

    /// Documented on `rank_percentile`: "a constant column maps entirely to
    /// 50". Ties share the mid-rank, and for an all-tied column every value
    /// shares the single mid-rank of the whole (0-indexed) run, which is
    /// exactly 50% of the way from 0 to 100 regardless of length — so this
    /// must hold for any length >= 2, not just the length-3 case above.
    #[test]
    fn rank_percentile_maps_an_all_tied_column_entirely_to_fifty() {
        for len in [2, 3, 4, 10, 11] {
            let values = vec![7.0; len];
            let (out, dropped) = transform_column(XTransform::RankPercentile, &values);
            assert_eq!(dropped, 0, "len {len}: nothing should be dropped");
            assert!(
                out.iter().all(|&v| (v - 50.0).abs() < 1e-9),
                "len {len}: an all-tied column must map entirely to 50, got {out:?}"
            );
        }
    }

    #[test]
    fn y_density_normalizes_per_class() {
        assert_eq!(YTransform::Density.apply(5, 10), 0.5);
        assert_eq!(YTransform::Count.apply(5, 10), 5.0);
        assert_eq!(YTransform::Log10Count.apply(9, 10), 1.0);
        assert_eq!(YTransform::Log10Count.apply(0, 10), 0.0, "log(0+1) = 0");
        assert_eq!(
            YTransform::Density.apply(5, 0),
            0.0,
            "empty class, no divide by zero"
        );
    }
}
