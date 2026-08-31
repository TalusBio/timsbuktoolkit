//! Neutral losses, keyed by atomic composition rather than by spelling.
//!
//! Libraries may write the same chemical loss in different ways:
//!
//! | written | library | composition |
//! |---|---|---|
//! | `-CH3SOH` | NIST | C1H4O1S1 |
//! | `-CH4OS` | SpectraST | C1H4O1S1 |
//! | `-NH2-CO-CH2SH` | NIST | C2H5N1O1S1 |
//! | `-C2H5NOS` | SpectraST | C2H5N1O1S1 |
//!
//! Parsing goes `text -> composition -> discriminant`; the packed annotation
//! stores only the discriminant.

use std::fmt::Display;

use crate::IonParsingError;

/// Atom counts for the elements that appear in peptide neutral losses.
///
/// These losses use only C/H/N/O/S/P.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct Composition {
    c: u8,
    h: u8,
    n: u8,
    o: u8,
    s: u8,
    p: u8,
}

/// A [`Composition`] naming only the elements it contains: `C!(h: 2, o: 1)`.
macro_rules! C {
    ($($field:ident: $count:expr),+ $(,)?) => {
        Composition { $($field: $count,)+ ..Composition::ZERO }
    };
}

impl Composition {
    /// All-zero, for [`TABLE`] rows to fill in only the elements they contain.
    const ZERO: Self = Self {
        c: 0,
        h: 0,
        n: 0,
        o: 0,
        s: 0,
        p: 0,
    };

    /// The count for one element symbol, or `None` if this crate does not
    /// represent that element.
    fn count_mut(&mut self, symbol: u8) -> Option<&mut u8> {
        Some(match symbol {
            b'C' => &mut self.c,
            b'H' => &mut self.h,
            b'N' => &mut self.n,
            b'O' => &mut self.o,
            b'S' => &mut self.s,
            b'P' => &mut self.p,
            _ => return None,
        })
    }

    /// The six counts, for tests that need to range over them.
    #[cfg(test)]
    fn counts(self) -> [u8; 6] {
        [self.c, self.h, self.n, self.o, self.s, self.p]
    }

    /// Combine element-wise, saturating each count.
    fn zip(self, other: Self, f: impl Fn(u8, u8) -> u8) -> Self {
        Self {
            c: f(self.c, other.c),
            h: f(self.h, other.h),
            n: f(self.n, other.n),
            o: f(self.o, other.o),
            s: f(self.s, other.s),
            p: f(self.p, other.p),
        }
    }

    /// Multiply every count. Used for the `2H2O` multiplier form.
    fn scaled(self, k: u8) -> Self {
        self.zip(Self::ZERO, |n, _| n.saturating_mul(k))
    }

    fn plus(self, other: Self) -> Self {
        self.zip(other, u8::saturating_add)
    }

    /// Parse a bare formula like `H2O`, `CH4OS`, `C2H5NOS`.
    ///
    /// Only single-letter C/H/N/O/S/P elements are recognized.
    fn parse_formula(s: &str) -> Result<Self, IonParsingError> {
        if s.is_empty() {
            return Err(IonParsingError::parse(s, "Empty neutral-loss formula"));
        }
        let mut out = Composition::default();
        let b = s.as_bytes();
        let mut i = 0;
        while i < b.len() {
            let symbol = b[i];
            i += 1;
            let start = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            let count: u8 = if start == i {
                1
            } else {
                s[start..i].parse().map_err(|_| {
                    IonParsingError::parse(s, "Neutral-loss atom count out of range")
                })?
            };
            let slot = out
                .count_mut(symbol)
                .ok_or_else(|| IonParsingError::parse(s, "Unsupported element in neutral loss"))?;
            *slot = slot.saturating_add(count);
        }
        Ok(out)
    }

    /// Parse a full loss expression: `-` separated terms, each optionally
    /// prefixed by a repeat count. `2H2O`, `H2O-NH3`, `NH2-CO-CH2SH`.
    ///
    /// Terms are summed, so ordering and multiplier spelling are normalized.
    pub(crate) fn parse_expression(s: &str) -> Result<Self, IonParsingError> {
        let mut total = Composition::default();
        for term in s.split('-') {
            let term = term.trim();
            if term.is_empty() {
                return Err(IonParsingError::parse(
                    s,
                    "Empty term in neutral-loss expression",
                ));
            }
            // Leading digits are a repeat count for the whole term.
            let digits = term
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(term.len());
            let (mult, formula) = if digits > 0 {
                let m: u8 = term[..digits].parse().map_err(|_| {
                    IonParsingError::parse(s, "Neutral-loss repeat count out of range")
                })?;
                (m, &term[digits..])
            } else {
                (1, term)
            };
            total = total.plus(Composition::parse_formula(formula)?.scaled(mult));
        }
        Ok(total)
    }
}

/// The neutral losses this crate can represent, as a packed discriminant.
///
/// Inputs outside this table are reported as unrepresentable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[repr(u8)]
pub enum NeutralLoss {
    #[default]
    None = 0,
    /// H2O
    Water = 1,
    /// NH3
    Ammonia = 2,
    /// CO
    CarbonMonoxide = 3,
    /// CO2
    CarbonDioxide = 4,
    /// 2 H2O
    WaterX2 = 5,
    /// 2 NH3
    AmmoniaX2 = 6,
    /// H2O + NH3
    WaterAmmonia = 7,
    /// CH4OS -- methanesulfenic acid, off oxidized Met. Also spelled `CH3SOH`.
    Methanesulfenic = 8,
    /// C2H5NOS -- also spelled `NH2-CO-CH2SH`.
    Carbamidomethylthiol = 9,
    /// H3PO4 -- phospho-Ser/Thr.
    PhosphoricAcid = 10,
    /// HPO3 -- phospho-Tyr, and phospho-Ser/Thr.
    Metaphosphoric = 11,
    /// H3PO4 + H2O
    PhosphoricAcidWater = 12,
}

/// `(composition, discriminant, canonical spelling)`, indexed by discriminant.
/// Entries are ordered by discriminant; canonical spelling is emitted by
/// [`Display`].
const TABLE: &[(Composition, NeutralLoss, &str)] = &[
    (C!(h: 2, o: 1), NeutralLoss::Water, "H2O"),
    (C!(h: 3, n: 1), NeutralLoss::Ammonia, "NH3"),
    (C!(c: 1, o: 1), NeutralLoss::CarbonMonoxide, "CO"),
    (C!(c: 1, o: 2), NeutralLoss::CarbonDioxide, "CO2"),
    (C!(h: 4, o: 2), NeutralLoss::WaterX2, "2H2O"),
    (C!(h: 6, n: 2), NeutralLoss::AmmoniaX2, "2NH3"),
    (C!(h: 5, n: 1, o: 1), NeutralLoss::WaterAmmonia, "H2O-NH3"),
    (
        C!(c: 1, h: 4, o: 1, s: 1),
        NeutralLoss::Methanesulfenic,
        "CH4OS",
    ),
    (
        C!(c: 2, h: 5, n: 1, o: 1, s: 1),
        NeutralLoss::Carbamidomethylthiol,
        "C2H5NOS",
    ),
    (C!(h: 3, o: 4, p: 1), NeutralLoss::PhosphoricAcid, "H3PO4"),
    (C!(h: 1, o: 3, p: 1), NeutralLoss::Metaphosphoric, "HPO3"),
    (
        C!(h: 5, o: 5, p: 1),
        NeutralLoss::PhosphoricAcidWater,
        "H3PO4-H2O",
    ),
];

impl NeutralLoss {
    /// Number of loss entries in the table.
    pub(crate) const COUNT: u8 = TABLE.len() as u8;

    /// Decode a packed loss discriminant.
    pub(crate) fn from_discriminant(d: u8) -> Self {
        Self::at(d).map_or(Self::None, |(_, loss, _)| *loss)
    }

    /// The [`TABLE`] row a discriminant names, or `None` for [`Self::None`] and
    /// for the reserved values above the table.
    fn at(d: u8) -> Option<&'static (Composition, NeutralLoss, &'static str)> {
        TABLE.get((d as usize).checked_sub(1)?)
    }

    /// Resolve a loss expression (without the leading `-`) to a discriminant.
    ///
    /// `None` means the composition is valid but not in the supported table.
    ///
    /// The lookup is by composition, so `CH3SOH` and `CH4OS` resolve to the
    /// same loss. Callers holding a loss as a formula string, rather than as an
    /// mzPAF annotation, come in here.
    pub fn from_expression(s: &str) -> Result<Option<Self>, IonParsingError> {
        let comp = Composition::parse_expression(s)?;
        Ok(TABLE
            .iter()
            .find(|(c, _, _)| *c == comp)
            .map(|(_, l, _)| *l))
    }

    /// Canonical spelling, without the leading `-`. Empty for [`Self::None`].
    pub(crate) fn canonical(self) -> &'static str {
        Self::at(self as u8).map_or("", |(_, _, spelling)| *spelling)
    }
}

impl Display for NeutralLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if *self == NeutralLoss::None {
            return Ok(());
        }
        write!(f, "-{}", self.canonical())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_counts_parse_and_reject_out_of_range() {
        assert_eq!(
            Composition::parse_expression("H2O").unwrap(),
            C!(h: 2, o: 1),
            "an element with no digits is one atom"
        );
        assert_eq!(
            Composition::parse_expression("C10H12").unwrap(),
            C!(c: 10, h: 12),
            "counts are multi-digit, not one digit per element"
        );
        assert_eq!(
            Composition::parse_expression("C255").unwrap(),
            C!(c: 255),
            "the u8 slot is full at 255"
        );
        assert!(
            Composition::parse_expression("C256").is_err(),
            "one past the slot must be an error, not a wrap to 0"
        );
    }

    #[test]
    fn different_spellings_resolve_to_one_loss() {
        assert_eq!(
            NeutralLoss::from_expression("CH3SOH").unwrap(),
            Some(NeutralLoss::Methanesulfenic)
        );
        assert_eq!(
            NeutralLoss::from_expression("CH4OS").unwrap(),
            Some(NeutralLoss::Methanesulfenic)
        );
        assert_eq!(
            NeutralLoss::from_expression("NH2-CO-CH2SH").unwrap(),
            Some(NeutralLoss::Carbamidomethylthiol)
        );
        assert_eq!(
            NeutralLoss::from_expression("C2H5NOS").unwrap(),
            Some(NeutralLoss::Carbamidomethylthiol)
        );
    }

    #[test]
    fn ordering_and_multipliers_normalize() {
        assert_eq!(
            NeutralLoss::from_expression("H2O-NH3").unwrap(),
            NeutralLoss::from_expression("NH3-H2O").unwrap()
        );
        assert_eq!(
            NeutralLoss::from_expression("2H2O").unwrap(),
            NeutralLoss::from_expression("H2O-H2O").unwrap()
        );
        assert_eq!(
            NeutralLoss::from_expression("2H2O").unwrap(),
            Some(NeutralLoss::WaterX2)
        );
    }

    #[test]
    fn unrepresentable_is_distinct_from_malformed() {
        assert_eq!(NeutralLoss::from_expression("HCOOH").unwrap(), None);
        assert!(NeutralLoss::from_expression("Xe2").is_err());
        assert!(NeutralLoss::from_expression("").is_err());
        assert!(NeutralLoss::from_expression("H2O-").is_err());
    }

    #[test]
    fn table_is_indexed_by_discriminant() {
        for (i, (_comp, loss, canon)) in TABLE.iter().enumerate() {
            assert_eq!(
                *loss as usize,
                i + 1,
                "{canon} sits at offset {i} but its discriminant is {}",
                *loss as u8
            );
        }
        // Nothing outside the table decodes to a loss, in either direction.
        assert_eq!(NeutralLoss::from_discriminant(0), NeutralLoss::None);
        assert_eq!(
            NeutralLoss::from_discriminant(TABLE.len() as u8 + 1),
            NeutralLoss::None
        );
        assert_eq!(NeutralLoss::None.canonical(), "");
    }

    #[test]
    fn table_round_trips_through_canonical_spelling() {
        for (_comp, loss, canon) in TABLE {
            assert_eq!(
                NeutralLoss::from_expression(canon).unwrap(),
                Some(*loss),
                "canonical spelling {canon} must resolve to its own loss"
            );
            assert_eq!(loss.canonical(), *canon);
            assert_eq!(
                NeutralLoss::from_discriminant(*loss as u8),
                *loss,
                "{canon} does not survive the discriminant round trip"
            );
        }
    }

    #[test]
    fn table_compositions_are_unique() {
        for (i, (a, _, sa)) in TABLE.iter().enumerate() {
            for (b, _, sb) in TABLE.iter().skip(i + 1) {
                assert_ne!(a, b, "{sa} and {sb} share a composition");
            }
            assert!(
                a.counts().iter().all(|&n| n < u8::MAX),
                "{sa} holds a saturated atom count"
            );
        }
    }

    #[test]
    fn display_adds_the_prefix_and_none_renders_empty() {
        assert_eq!(NeutralLoss::None.to_string(), "");
        assert_eq!(NeutralLoss::Water.to_string(), "-H2O");
    }
}
