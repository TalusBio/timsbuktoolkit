//! Axis transforms for the histogram panels.
//!
//! A transform never errors. Values it cannot map (non-positives under a log,
//! negatives under a square root, anything non-finite) are dropped, and the
//! caller reports the drop count in the panel subtitle — a silently shrinking
//! histogram would misread as "these rows do not exist".
//!
//! [`XTransform`] is applied once, at init, over a sorted sample. (`YTransform`
//! is not: it maps stored counts and runs per bin per frame.) That order is
//! what [`XTransform::accepts`] and [`XTransform::is_monotone`] exist for: the
//! survivors of a domain restriction are a *suffix* of a sorted column, so a
//! `partition_point` finds them, and for a monotone map the p-th percentile of
//! `T(x)` is `T` of the p-th percentile of `x`.

use crate::cycle;

/// Largest exponent `exp` is allowed to see, guarding against overflow to
/// `+inf` above ~709.78.
const EXP_CLAMP: f64 = 700.0;

/// What a histogram's x axis shows.
///
/// `RankPercentile` is the row's sorted position, not a function of its value,
/// so it has no pointwise map and no domain to accept or refuse. Keeping it out
/// of [`XTransform`] is what lets [`XTransform::apply`] be total.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// A pointwise map of the value.
    Value(XTransform),
    /// Sorted position as a percentile in `0..=100`.
    RankPercentile,
}

impl Axis {
    /// Every axis, in cycle order. Stored histograms are indexed by position
    /// in this array.
    pub const ALL: [Axis; 7] = [
        Self::Value(XTransform::Linear),
        Self::Value(XTransform::Log10),
        Self::Value(XTransform::SignedLog1p),
        Self::Value(XTransform::Sqrt),
        Self::Value(XTransform::Square),
        Self::Value(XTransform::Exp),
        Self::RankPercentile,
    ];

    pub fn next(self) -> Self {
        cycle::step(&Self::ALL, self, 1)
    }

    pub fn prev(self) -> Self {
        cycle::step(&Self::ALL, self, -1)
    }

    /// Index into [`Self::ALL`], which is how stored histograms are addressed.
    pub fn index(self) -> usize {
        cycle::index_of(&Self::ALL, self)
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Value(t) => t.label(),
            Self::RankPercentile => "rank-pct",
        }
    }

    /// Whether a *finite* `v` can be plotted on this axis. Every value has a
    /// rank, so only a value transform can refuse one.
    pub fn accepts(self, v: f64) -> bool {
        match self {
            Self::Value(t) => t.accepts(v),
            Self::RankPercentile => true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XTransform {
    Linear,
    Log10,
    SignedLog1p,
    Sqrt,
    Square,
    Exp,
}

impl XTransform {
    pub fn label(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Log10 => "log10",
            Self::SignedLog1p => "signed-log1p",
            Self::Sqrt => "sqrt",
            Self::Square => "square",
            Self::Exp => "exp",
        }
    }

    /// Whether a *finite* `v` is inside this transform's domain — exactly the
    /// condition [`Self::apply`] tests.
    pub fn accepts(self, v: f64) -> bool {
        match self {
            Self::Log10 => v > 0.0,
            Self::Sqrt => v >= 0.0,
            _ => true,
        }
    }

    /// Whether this map is non-decreasing over the values it accepts. `Square`
    /// is the one that is not: it decreases below zero.
    pub fn is_monotone(self) -> bool {
        !matches!(self, Self::Square)
    }

    /// Pointwise map. `None` means the value is outside this transform's
    /// domain (or maps outside the finite range) and must be dropped.
    pub fn apply(self, v: f64) -> Option<f64> {
        if !v.is_finite() {
            return None;
        }
        let out = match self {
            Self::Linear => v,
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
    pub const ALL: [YTransform; 3] = [Self::Density, Self::Count, Self::Log10Count];

    pub fn next(self) -> Self {
        cycle::step(&Self::ALL, self, 1)
    }

    pub fn prev(self) -> Self {
        cycle::step(&Self::ALL, self, -1)
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

    /// The value transforms, derived from [`Axis::ALL`] rather than listed again.
    fn value_transforms() -> Vec<XTransform> {
        Axis::ALL
            .into_iter()
            .filter_map(|a| match a {
                Axis::Value(t) => Some(t),
                Axis::RankPercentile => None,
            })
            .collect()
    }

    /// Stored histograms are addressed by `index()`, so two axes sharing one
    /// index would alias their slots. Cycling itself is `cycle`'s test.
    #[test]
    fn axis_index_is_injective() {
        let mut seen = Axis::ALL.to_vec();
        seen.sort_by_key(|a| a.index());
        seen.dedup();
        assert_eq!(seen.len(), Axis::ALL.len(), "every axis, once");
    }

    /// Each map on a value it accepts, and on the values just outside its
    /// domain.
    #[test]
    fn each_transform_maps_its_domain_and_refuses_the_rest() {
        let cases = [
            (XTransform::Linear, 3.0, Some(3.0)),
            (XTransform::Log10, 100.0, Some(2.0)),
            (XTransform::Log10, 0.0, None),
            (XTransform::Log10, -1.0, None),
            (XTransform::Sqrt, 4.0, Some(2.0)),
            (XTransform::Sqrt, 0.0, Some(0.0)),
            (XTransform::Sqrt, -9.0, None),
            (XTransform::Square, -3.0, Some(9.0)),
            (XTransform::Square, 2.0, Some(4.0)),
        ];
        for (t, v, want) in cases {
            assert_eq!(t.apply(v), want, "{t:?} on {v}");
        }
    }

    #[test]
    fn signed_log1p_is_total_and_odd() {
        let neg = XTransform::SignedLog1p.apply(-9.0).unwrap();
        let pos = XTransform::SignedLog1p.apply(9.0).unwrap();
        assert!((neg + pos).abs() < 1e-12, "odd about zero");
        assert_eq!(XTransform::SignedLog1p.apply(0.0), Some(0.0));
    }

    #[test]
    fn exp_clamps_instead_of_overflowing() {
        assert_eq!(XTransform::Exp.apply(0.0), Some(1.0));
        let big = XTransform::Exp
            .apply(1e9)
            .expect("1e9 must clamp, not drop");
        assert!(big.is_finite(), "got {big}");
        assert_eq!(XTransform::Exp.apply(-1e9), Some(0.0));
    }

    #[test]
    fn non_finite_input_is_always_rejected() {
        for t in value_transforms() {
            assert_eq!(t.apply(f64::NAN), None, "{t:?} must reject NaN");
            assert_eq!(t.apply(f64::INFINITY), None, "{t:?} must reject +inf");
            assert_eq!(t.apply(f64::NEG_INFINITY), None, "{t:?} must reject -inf");
        }
    }

    /// `accepts` drives the sorted-suffix search that finds each transform's
    /// survivors, while `apply` decides what is actually plotted. If they
    /// disagree on any finite value, the clip range is taken over a different
    /// set than the one being binned.
    #[test]
    fn accepts_agrees_with_apply_on_every_finite_value() {
        let probes = [
            -1e6,
            -1.0,
            -1e-9,
            -f64::MIN_POSITIVE,
            0.0,
            f64::MIN_POSITIVE,
            1e-9,
            0.5,
            1.0,
            2.0,
            1e6,
        ];
        for t in value_transforms() {
            for v in probes {
                assert_eq!(
                    t.accepts(v),
                    t.apply(v).is_some(),
                    "{t:?} disagrees on {v:e}"
                );
            }
        }
    }

    /// The accepted set must be an up-set: everything at or above the smallest
    /// accepted value is accepted too. That is what makes the survivors a
    /// suffix of a sorted column rather than a scattered subset.
    #[test]
    fn accepted_values_form_a_suffix_of_the_value_axis() {
        let ascending = [-1e6, -1.0, -1e-9, 0.0, 1e-9, 0.5, 1.0, 1e6];
        for t in value_transforms() {
            let mut seen_accepted = false;
            for v in ascending {
                let ok = t.accepts(v);
                assert!(
                    ok || !seen_accepted,
                    "{t:?} rejects {v:e} after accepting a smaller value"
                );
                seen_accepted |= ok;
            }
        }
    }

    /// A monotone claim is load-bearing: the clip range is computed as
    /// `T(percentile of x)` rather than `percentile of T(x)`, which is only the
    /// same thing when `T` is non-decreasing over its accepted values.
    #[test]
    fn monotone_transforms_really_are_non_decreasing() {
        let ascending = [-1e3, -1.0, -0.5, 0.0, 1e-9, 0.5, 1.0, 2.0, 1e3];
        for t in value_transforms().into_iter().filter(|t| t.is_monotone()) {
            let mut prev = f64::NEG_INFINITY;
            for v in ascending {
                let Some(y) = t.apply(v) else { continue };
                assert!(y >= prev, "{t:?} decreased at {v:e}: {y} < {prev}");
                prev = y;
            }
        }
        assert!(
            !XTransform::Square.is_monotone(),
            "square decreases below zero"
        );
    }

    #[test]
    fn y_density_normalizes_per_class() {
        assert_eq!(YTransform::Log10Count.apply(0, 10), 0.0, "log(0+1) = 0");
        assert_eq!(
            YTransform::Density.apply(5, 0),
            0.0,
            "empty class, no divide by zero"
        );
    }
}
