//! Neutral losses, keyed by atomic composition rather than by spelling.
//!
//! Libraries write the same chemical loss different ways. Two real examples
//! from the HUPO-PSI mzSpecLib corpus:
//!
//! | written | library | composition |
//! |---|---|---|
//! | `-CH3SOH` | NIST | C1H4O1S1 |
//! | `-CH4OS` | SpectraST | C1H4O1S1 |
//! | `-NH2-CO-CH2SH` | NIST | C2H5N1O1S1 |
//! | `-C2H5NOS` | SpectraST | C2H5N1O1S1 |
//!
//! Keying on the string would make `y5-CH4OS` and `y5-CH3SOH` distinct labels
//! for one ion. Since fragment labels must be unique within a precursor (see
//! `ExpectedIntensities::try_from_pairs` in timsseek), that is exactly the
//! wrong direction: it hides a genuine duplicate behind two spellings.
//!
//! So parsing goes `text -> composition -> discriminant`, and the composition
//! is a *parse-time* concept only. What gets stored on an `IonAnnot` is the
//! discriminant, so there is no per-annotation cost at runtime.

use std::fmt::Display;

use crate::IonParsingError;

/// Slots in [`Composition`], in the order [`Composition::new`] takes them.
const C: usize = 0;
const H: usize = 1;
const N: usize = 2;
const O: usize = 3;
const S: usize = 4;
const P: usize = 5;

/// Atom counts for the elements that appear in peptide neutral losses.
///
/// Deliberately not a general chemical formula: these losses only ever draw
/// from C/H/N/O/S/P, and keeping it to six `u8`s makes equality a single
/// 6-byte compare during the parse-time table lookup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct Composition([u8; 6]);

impl Composition {
    pub(crate) const fn new(c: u8, h: u8, n: u8, o: u8, s: u8, p: u8) -> Self {
        Self([c, h, n, o, s, p])
    }

    /// Multiply every count, saturating. Used for the `2H2O` multiplier form.
    fn scaled(self, k: u8) -> Self {
        Self(self.0.map(|n| n.saturating_mul(k)))
    }

    fn plus(self, other: Self) -> Self {
        Self(std::array::from_fn(|i| {
            self.0[i].saturating_add(other.0[i])
        }))
    }

    /// Parse a bare formula like `H2O`, `CH4OS`, `C2H5NOS`.
    ///
    /// Only single-letter C/H/N/O/S/P are recognized; anything else is an
    /// error rather than a silent skip, so an unsupported loss surfaces as
    /// "not representable" instead of being mistaken for a smaller one.
    fn parse_formula(s: &str) -> Result<Self, IonParsingError> {
        if s.is_empty() {
            return Err(IonParsingError::ParsingError {
                error: s.to_string(),
                context: Some("Empty neutral-loss formula"),
            });
        }
        let mut out = Composition::default();
        let b = s.as_bytes();
        let mut i = 0;
        while i < b.len() {
            let elem = b[i];
            i += 1;
            let start = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            let count: u8 = if start == i {
                1
            } else {
                s[start..i]
                    .parse()
                    .map_err(|_| IonParsingError::ParsingError {
                        error: s.to_string(),
                        context: Some("Neutral-loss atom count out of range"),
                    })?
            };
            let slot = match elem {
                b'C' => C,
                b'H' => H,
                b'N' => N,
                b'O' => O,
                b'S' => S,
                b'P' => P,
                _ => {
                    return Err(IonParsingError::ParsingError {
                        error: s.to_string(),
                        context: Some("Unsupported element in neutral loss"),
                    });
                }
            };
            out.0[slot] = out.0[slot].saturating_add(count);
        }
        Ok(out)
    }

    /// Parse a full loss expression: `-` separated terms, each optionally
    /// prefixed by a repeat count. `2H2O`, `H2O-NH3`, `NH2-CO-CH2SH`.
    ///
    /// Because terms are summed, ordering and multiplier spelling collapse for
    /// free: `H2O-NH3` == `NH3-H2O`, and `2H2O` == `H2O-H2O`.
    pub(crate) fn parse_expression(s: &str) -> Result<Self, IonParsingError> {
        let mut total = Composition::default();
        for term in s.split('-') {
            let term = term.trim();
            if term.is_empty() {
                return Err(IonParsingError::ParsingError {
                    error: s.to_string(),
                    context: Some("Empty term in neutral-loss expression"),
                });
            }
            // Leading digits are a repeat count for the whole term.
            let digits = term
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(term.len());
            let (mult, formula) = if digits > 0 {
                let m: u8 = term[..digits]
                    .parse()
                    .map_err(|_| IonParsingError::ParsingError {
                        error: s.to_string(),
                        context: Some("Neutral-loss repeat count out of range"),
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
/// Scoped deliberately: DIA-NN emits none, Spectronaut two (`-H2O`, `-NH3`),
/// NIST eight, plus the phospho losses that no non-phospho corpus can show.
/// SpectraST's wider combinatorics are NOT here -- an unlisted loss parses to a
/// composition that misses the table and is reported as unrepresentable, which
/// routes the peak to an unknown label rather than silently mislabelling it.
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

/// `(composition, discriminant, canonical spelling)`.
///
/// The canonical spelling is what `Display` emits, so a non-canonical input
/// (`-CH3SOH`) round-trips to the canonical form (`-CH4OS`). Round-trip tests
/// must therefore compare parsed values, not bytes.
const TABLE: &[(Composition, NeutralLoss, &str)] = &[
    (
        Composition::new(0, 2, 0, 1, 0, 0),
        NeutralLoss::Water,
        "H2O",
    ),
    (
        Composition::new(0, 3, 1, 0, 0, 0),
        NeutralLoss::Ammonia,
        "NH3",
    ),
    (
        Composition::new(1, 0, 0, 1, 0, 0),
        NeutralLoss::CarbonMonoxide,
        "CO",
    ),
    (
        Composition::new(1, 0, 0, 2, 0, 0),
        NeutralLoss::CarbonDioxide,
        "CO2",
    ),
    (
        Composition::new(0, 4, 0, 2, 0, 0),
        NeutralLoss::WaterX2,
        "2H2O",
    ),
    (
        Composition::new(0, 6, 2, 0, 0, 0),
        NeutralLoss::AmmoniaX2,
        "2NH3",
    ),
    (
        Composition::new(0, 5, 1, 1, 0, 0),
        NeutralLoss::WaterAmmonia,
        "H2O-NH3",
    ),
    (
        Composition::new(1, 4, 0, 1, 1, 0),
        NeutralLoss::Methanesulfenic,
        "CH4OS",
    ),
    (
        Composition::new(2, 5, 1, 1, 1, 0),
        NeutralLoss::Carbamidomethylthiol,
        "C2H5NOS",
    ),
    (
        Composition::new(0, 3, 0, 4, 0, 1),
        NeutralLoss::PhosphoricAcid,
        "H3PO4",
    ),
    (
        Composition::new(0, 1, 0, 3, 0, 1),
        NeutralLoss::Metaphosphoric,
        "HPO3",
    ),
    (
        Composition::new(0, 5, 0, 5, 0, 1),
        NeutralLoss::PhosphoricAcidWater,
        "H3PO4-H2O",
    ),
];

impl NeutralLoss {
    /// Inverse of the `#[repr(u8)]` discriminant, for unpacking out of a bit
    /// field. Lives next to the enum so the two cannot drift apart.
    ///
    /// An unrecognized value maps to [`Self::None`]: the only way to produce
    /// one is a reserved discriminant, which no constructor emits.
    pub(crate) fn from_discriminant(d: u8) -> Self {
        match d {
            1 => Self::Water,
            2 => Self::Ammonia,
            3 => Self::CarbonMonoxide,
            4 => Self::CarbonDioxide,
            5 => Self::WaterX2,
            6 => Self::AmmoniaX2,
            7 => Self::WaterAmmonia,
            8 => Self::Methanesulfenic,
            9 => Self::Carbamidomethylthiol,
            10 => Self::PhosphoricAcid,
            11 => Self::Metaphosphoric,
            12 => Self::PhosphoricAcidWater,
            _ => Self::None,
        }
    }

    /// Resolve a loss expression (without the leading `-`) to a discriminant.
    ///
    /// `None` means "parsed as a valid composition, but not one we represent" --
    /// distinct from `Err`, which means the text was not a loss expression at
    /// all. Callers route the former to an unknown label and the latter to a
    /// parse failure.
    pub(crate) fn from_expression(s: &str) -> Result<Option<Self>, IonParsingError> {
        let comp = Composition::parse_expression(s)?;
        Ok(TABLE
            .iter()
            .find(|(c, _, _)| *c == comp)
            .map(|(_, l, _)| *l))
    }

    /// Canonical spelling, without the leading `-`. Empty for [`Self::None`].
    pub(crate) fn canonical(self) -> &'static str {
        if self == NeutralLoss::None {
            return "";
        }
        TABLE
            .iter()
            .find(|(_, l, _)| *l == self)
            .map(|(_, _, s)| *s)
            .unwrap_or("")
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

    /// The count parser: implicit 1, multi-digit counts, and a count past the
    /// `u8` slot. Single-digit cases are covered by every TABLE row.
    #[test]
    fn atom_counts_parse_and_reject_out_of_range() {
        assert_eq!(
            Composition::parse_expression("H2O").unwrap(),
            Composition::new(0, 2, 0, 1, 0, 0),
            "an element with no digits is one atom"
        );
        assert_eq!(
            Composition::parse_expression("C10H12").unwrap(),
            Composition::new(10, 12, 0, 0, 0, 0),
            "counts are multi-digit, not one digit per element"
        );
        assert_eq!(
            Composition::parse_expression("C255").unwrap(),
            Composition::new(255, 0, 0, 0, 0, 0),
            "the u8 slot is full at 255"
        );
        assert!(
            Composition::parse_expression("C256").is_err(),
            "one past the slot must be an error, not a wrap to 0"
        );
    }

    /// The two cross-library spelling collisions this module exists for.
    #[test]
    fn different_spellings_resolve_to_one_loss() {
        // NIST vs SpectraST, methanesulfenic acid.
        assert_eq!(
            NeutralLoss::from_expression("CH3SOH").unwrap(),
            Some(NeutralLoss::Methanesulfenic)
        );
        assert_eq!(
            NeutralLoss::from_expression("CH4OS").unwrap(),
            Some(NeutralLoss::Methanesulfenic)
        );
        // NIST structural notation vs SpectraST molecular notation.
        assert_eq!(
            NeutralLoss::from_expression("NH2-CO-CH2SH").unwrap(),
            Some(NeutralLoss::Carbamidomethylthiol)
        );
        assert_eq!(
            NeutralLoss::from_expression("C2H5NOS").unwrap(),
            Some(NeutralLoss::Carbamidomethylthiol)
        );
    }

    /// Summing terms collapses ordering and multiplier spelling for free.
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

    /// A well-formed composition outside the table is `Ok(None)` -- "valid but
    /// not representable" -- while malformed text is `Err`. Callers need to
    /// tell those apart to route one to an unknown label and the other to a
    /// parse failure.
    #[test]
    fn unrepresentable_is_distinct_from_malformed() {
        assert_eq!(NeutralLoss::from_expression("HCOOH").unwrap(), None);
        assert!(NeutralLoss::from_expression("Xe2").is_err());
        assert!(NeutralLoss::from_expression("").is_err());
        assert!(NeutralLoss::from_expression("H2O-").is_err());
    }

    /// Every table entry must survive canonical -> composition -> discriminant,
    /// and back out through the bit field. `from_discriminant` hand-mirrors the
    /// `#[repr(u8)]` values, so adding a loss without updating it would
    /// silently decode as [`NeutralLoss::None`] -- this is what catches that.
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
        assert_eq!(NeutralLoss::from_discriminant(0), NeutralLoss::None);
        // Every non-None variant must be in TABLE, or it has no spelling and no
        // composition and could never be produced by parsing.
        assert_eq!(
            TABLE.len(),
            (1..=u8::MAX)
                .filter(|d| NeutralLoss::from_discriminant(*d) != NeutralLoss::None)
                .count(),
            "a variant is decodable but missing from TABLE"
        );
    }

    /// Compositions must be unique: two entries sharing one would make the
    /// table lookup order-dependent.
    #[test]
    fn table_compositions_are_unique() {
        for (i, (a, _, sa)) in TABLE.iter().enumerate() {
            for (b, _, sb) in TABLE.iter().skip(i + 1) {
                assert_ne!(a, b, "{sa} and {sb} share a composition");
            }
            // `scaled` and `plus` saturate, so an absurd input like `200H2O`
            // lands on 255 in some slot. That is only safe to report as
            // unrepresentable while no table entry holds a saturated count --
            // otherwise the saturation would alias onto a real loss.
            assert!(
                a.0.iter().all(|&n| n < u8::MAX),
                "{sa} holds a saturated atom count"
            );
        }
    }

    /// The two things the TABLE round trip cannot see: the `-` prefix Display
    /// adds, and that `None` renders as nothing at all rather than "-".
    #[test]
    fn display_adds_the_prefix_and_none_renders_empty() {
        assert_eq!(NeutralLoss::None.to_string(), "");
        assert_eq!(NeutralLoss::Water.to_string(), "-H2O");
    }
}
