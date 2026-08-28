//! What kind of ion an annotation names, and how that fits in a packed payload.
//!
//! This module owns three related details:
//!
//! 1. the mzPAF **spelling** of each kind (`Display` and `IonSeriesOrdinal::parse`,
//!    which are exact inverses),
//! 2. the **payload layout** each kind packs into (`IonSeriesOrdinal::to_parts`
//!    and `from_parts`, also inverses),
//! 3. the **bounds** that layout imposes, enforced by the constructors.
//!
//! [`IonAnnot`](crate::IonAnnot) owns the surrounding word -- charge, isotope,
//! loss, and their positions.

use std::fmt::Display;

use crate::{
    IonParsingError,
    mask,
};

/// Width of the `kind` discriminant in the packed word.
pub(crate) const KIND_BITS: u32 = 4;
/// Width of the `payload` this module encodes into.
pub(crate) const PAYLOAD_BITS: u32 = 12;
/// Width of each internal-fragment endpoint inside `payload`. Two share it.
const INTERNAL_POS_BITS: u32 = 6;
/// Width of the immonium residue index inside `payload`.
const IMMONIUM_BITS: u32 = 5;

/// Widest residue position an internal fragment endpoint holds.
///
/// Narrower than a backbone ordinal because two endpoints share one payload.
/// Internal fragments are bounded by peptide length, so 63 is well past tryptic.
pub const INTERNAL_POS_MAX: u8 = mask(INTERNAL_POS_BITS) as u8;

const _: () = assert!(
    2 * INTERNAL_POS_BITS <= PAYLOAD_BITS,
    "two internal-fragment endpoints must fit in one payload"
);
const _: () = assert!(
    IMMONIUM_BITS <= PAYLOAD_BITS && mask(IMMONIUM_BITS) >= (b'Z' - b'A') as u32,
    "the immonium field must fit the payload and hold every uppercase residue"
);

/// `kind` discriminants. The numbering lives only here; [`IonSeriesOrdinal::to_parts`]
/// and `from_parts` are its only readers.
///
/// Backbone ions take 1..=9 from [`Series`] itself, so a series can never drift
/// out of step with its letter. `unknown` is 0 so that the all-zero word decodes
/// to a real annotation. 13..=15 are unassigned and decode as `unknown`.
const KIND_UNKNOWN: u32 = 0;
const KIND_PRECURSOR: u32 = 10;
const KIND_INTERNAL: u32 = 11;
const KIND_IMMONIUM: u32 = 12;

const _: () = assert!(
    KIND_IMMONIUM <= mask(KIND_BITS),
    "a kind discriminant overflows KIND_BITS"
);

/// A backbone fragment ion series.
///
/// The nine mzPAF backbone series differ only by their letter, so the
/// letter-to-discriminant pairing lives in exactly one private table and every
/// match over them is a single arm.
#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
#[allow(non_camel_case_types)]
#[repr(u8)]
pub enum Series {
    a = 1,
    b,
    c,
    d,
    v,
    w,
    x,
    y,
    z,
}

impl Series {
    /// Every series, in discriminant order.
    pub const ALL: [Self; 9] = [
        Self::a,
        Self::b,
        Self::c,
        Self::d,
        Self::v,
        Self::w,
        Self::x,
        Self::y,
        Self::z,
    ];
    /// The mzPAF letters, in discriminant order.
    const CHARS: &'static [u8; 9] = b"abcdvwxyz";

    /// The mzPAF letter for this series.
    pub const fn as_char(self) -> char {
        Self::CHARS[self as usize - 1] as char
    }

    /// The series for an mzPAF letter, or `None` if it names no backbone series.
    fn from_char(c: char) -> Option<Self> {
        let idx = Self::CHARS.iter().position(|&b| b as char == c)?;
        Some(Self::ALL[idx])
    }
}

impl Display for Series {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

/// What kind of ion an annotation names, with whatever that kind carries.
///
/// A *view*: [`IonAnnot`](crate::IonAnnot) stores a packed word and reconstructs
/// this on demand, so building one allocates no annotation. Fields are public
/// because reading them is the point; the validating constructors
/// ([`Self::try_internal`], [`Self::try_immonium`]) are how you get one that is known to
/// fit the payload.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Copy)]
#[allow(non_camel_case_types)]
pub enum IonSeriesOrdinal {
    /// One of the nine backbone series, at position `ordinal` in the ladder.
    backbone {
        series: Series,
        ordinal: u8,
    },
    /// An unannotated peak. `ordinal` is a uniqueness counter, not a position.
    unknown {
        ordinal: u8,
    },
    precursor,
    /// An internal fragment spanning residues `start..=end`.
    internal {
        start: u8,
        end: u8,
    },
    /// A bare immonium ion for an uppercase residue code.
    immonium {
        residue: char,
    },
}

impl IonSeriesOrdinal {
    /// An internal fragment spanning residues `start..=end`.
    ///
    /// Endpoints are bounded by [`INTERNAL_POS_MAX`], which is this module's
    /// payload width and is checked here rather than by the caller.
    pub fn try_internal(start: u8, end: u8) -> Result<Self, IonParsingError> {
        if start > INTERNAL_POS_MAX || end > INTERNAL_POS_MAX {
            return Err(IonParsingError::OrdinalOutOfRange {
                ordinal: start.max(end),
                series: 'm',
            });
        }
        Ok(Self::internal { start, end })
    }

    /// A bare immonium ion for an uppercase residue code.
    ///
    /// Modified immonium (`IC[Carbamidomethyl]`) carries an arbitrary mod string
    /// that no fixed-width field holds, so only a bare uppercase residue is
    /// accepted.
    pub fn try_immonium(residue: char) -> Result<Self, IonParsingError> {
        if !residue.is_ascii_uppercase() {
            return Err(IonParsingError::UnsupportedFragmentType {
                fragment_type: residue,
            });
        }
        Ok(Self::immonium { residue })
    }

    /// A backbone / precursor / unknown ion from its mzPAF letter.
    ///
    /// `p` is the only letter that takes no ordinal; every other one requires
    /// one. Internal fragments and immonium ions are spelled differently and
    /// have their own constructors.
    pub fn from_series_char(c: char, ordinal: Option<u8>) -> Result<Self, IonParsingError> {
        if c == 'p' {
            return match ordinal {
                None => Ok(Self::precursor),
                Some(ordinal) => Err(IonParsingError::UnexpectedOrdinal { series: c, ordinal }),
            };
        }
        let ordinal = ordinal.ok_or(IonParsingError::MissingOrdinal { series: c })?;
        if c == '?' {
            return Ok(Self::unknown { ordinal });
        }
        Series::from_char(c)
            .map(|series| Self::backbone { series, ordinal })
            .ok_or(IonParsingError::UnsupportedFragmentType { fragment_type: c })
    }

    /// Parse the ion-kind part of an annotation: everything left of the charge,
    /// isotope and loss suffixes.
    ///
    /// The inverse of this type's [`Display`].
    pub(crate) fn parse(core: &str) -> Result<Self, IonParsingError> {
        // Internal fragment: `m<start>:<end>`. Checked before the single-letter
        // forms because `m` is not a backbone letter, so there is no ambiguity.
        if let Some(spans) = core.strip_prefix('m')
            && let Some((a, b)) = spans.split_once(':')
        {
            let start = a
                .parse::<u8>()
                .map_err(|_| IonParsingError::parse(core, "internal-fragment start is not a u8"))?;
            let end = b
                .parse::<u8>()
                .map_err(|_| IonParsingError::parse(core, "internal-fragment end is not a u8"))?;
            return Self::try_internal(start, end);
        }

        // Immonium: `I<residue>`, bare only.
        if let Some(res) = core.strip_prefix('I') {
            let mut ch = res.chars();
            return match (ch.next(), ch.next()) {
                (Some(r), None) => Self::try_immonium(r),
                (None, _) => Err(IonParsingError::parse(
                    core,
                    "immonium ion names no residue",
                )),
                // Carries a mod string no fixed-width field can hold.
                _ => Err(IonParsingError::UnsupportedModifiedImmonium {
                    annotation: core.to_string(),
                }),
            };
        }

        // Backbone / precursor / unknown: a series letter then an ordinal.
        let mut chars = core.chars();
        let series = chars
            .next()
            .ok_or(IonParsingError::parse(core, "empty annotation"))?;
        let rest = chars.as_str();
        let ordinal =
            if rest.is_empty() {
                None
            } else {
                Some(rest.parse::<u8>().map_err(|_| {
                    IonParsingError::parse(core, "ordinal is not a number in 0..=255")
                })?)
            };
        Self::from_series_char(series, ordinal)
    }

    /// Split into the `kind` discriminant and its `payload`, the two fields
    /// [`IonAnnot`](crate::IonAnnot) packs.
    ///
    /// Split into the fields used by the packed word.
    pub(crate) const fn to_parts(self) -> (u32, u32) {
        match self {
            Self::backbone { series, ordinal } => (series as u32, ordinal as u32),
            Self::unknown { ordinal } => (KIND_UNKNOWN, ordinal as u32),
            Self::precursor => (KIND_PRECURSOR, 0),
            Self::internal { start, end } => (
                KIND_INTERNAL,
                ((start as u32) & mask(INTERNAL_POS_BITS))
                    | (((end as u32) & mask(INTERNAL_POS_BITS)) << INTERNAL_POS_BITS),
            ),
            // `saturating_sub` rather than `-`: the constructor rejects a
            // non-uppercase residue, but the variant is public and this must not
            // underflow for one built directly.
            Self::immonium { residue } => (
                KIND_IMMONIUM,
                (residue as u32).saturating_sub(b'A' as u32) & mask(IMMONIUM_BITS),
            ),
        }
    }

    /// Inverse of [`Self::to_parts`]. Unknown discriminants decode as unknown
    /// ions.
    pub(crate) const fn from_parts(kind: u32, payload: u32) -> Self {
        let ordinal = payload as u8;
        match kind {
            1..=9 => Self::backbone {
                // `kind` is in range, so this is the inverse of `series as u32`.
                series: Series::ALL[kind as usize - 1],
                ordinal,
            },
            KIND_PRECURSOR => Self::precursor,
            KIND_INTERNAL => Self::internal {
                start: (payload & mask(INTERNAL_POS_BITS)) as u8,
                end: ((payload >> INTERNAL_POS_BITS) & mask(INTERNAL_POS_BITS)) as u8,
            },
            KIND_IMMONIUM => Self::immonium {
                residue: (b'A' + (payload & mask(IMMONIUM_BITS)) as u8) as char,
            },
            _ => Self::unknown { ordinal },
        }
    }
}

impl Display for IonSeriesOrdinal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::backbone { series, ordinal } => write!(f, "{}{}", series.as_char(), ordinal),
            Self::unknown { ordinal } => write!(f, "?{}", ordinal),
            Self::precursor => write!(f, "p"),
            Self::internal { start, end } => write!(f, "m{}:{}", start, end),
            Self::immonium { residue } => write!(f, "I{}", residue),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_series() -> Vec<IonSeriesOrdinal> {
        let mut out: Vec<IonSeriesOrdinal> = Series::ALL
            .iter()
            .enumerate()
            .map(|(i, &series)| IonSeriesOrdinal::backbone {
                series,
                ordinal: i as u8 + 1,
            })
            .collect();
        out.extend([
            IonSeriesOrdinal::precursor,
            IonSeriesOrdinal::unknown { ordinal: 10 },
            IonSeriesOrdinal::internal { start: 2, end: 11 },
            IonSeriesOrdinal::internal {
                start: INTERNAL_POS_MAX,
                end: INTERNAL_POS_MAX,
            },
            IonSeriesOrdinal::immonium { residue: 'A' },
            IonSeriesOrdinal::immonium { residue: 'Z' },
        ]);
        out
    }

    #[test]
    fn every_variant_round_trips_through_the_packed_parts() {
        for series in all_series() {
            let (kind, payload) = series.to_parts();
            assert_eq!(IonSeriesOrdinal::from_parts(kind, payload), series);
            assert!(kind <= mask(KIND_BITS), "{series:?} kind overflows");
            assert!(
                payload <= mask(PAYLOAD_BITS),
                "{series:?} payload overflows"
            );
        }

        let mut kinds: Vec<u32> = Series::ALL
            .iter()
            .map(|&series| {
                IonSeriesOrdinal::backbone { series, ordinal: 1 }
                    .to_parts()
                    .0
            })
            .collect();
        kinds.sort_unstable();
        kinds.dedup();
        assert_eq!(kinds.len(), Series::ALL.len(), "two series share a kind");
    }

    #[test]
    fn every_kind_bit_pattern_decodes_without_panicking() {
        for kind in 0..=mask(KIND_BITS) {
            for payload in [0, 1, mask(PAYLOAD_BITS)] {
                let _ = IonSeriesOrdinal::from_parts(kind, payload).to_string();
            }
        }
    }

    #[test]
    fn every_variant_round_trips_through_its_mzpaf_spelling() {
        for series in all_series() {
            let text = series.to_string();
            assert_eq!(
                IonSeriesOrdinal::parse(&text).expect("own spelling parses"),
                series,
                "{text}"
            );
        }
    }

    #[test]
    fn series_letters_and_discriminants_agree() {
        assert_eq!(Series::ALL.len(), Series::CHARS.len());
        for &series in &Series::ALL {
            assert_eq!(Series::from_char(series.as_char()), Some(series));
        }
        for c in ['p', '?', 'm', 'I', 'q', 'A'] {
            assert_eq!(Series::from_char(c), None, "{c} is not a backbone series");
        }
    }

    #[test]
    fn constructors_reject_what_the_payload_cannot_hold() {
        assert!(IonSeriesOrdinal::try_internal(INTERNAL_POS_MAX, INTERNAL_POS_MAX).is_ok());
        assert!(matches!(
            IonSeriesOrdinal::try_internal(INTERNAL_POS_MAX + 1, 1),
            Err(IonParsingError::OrdinalOutOfRange { .. })
        ));
        assert!(IonSeriesOrdinal::try_immonium('A').is_ok());
        assert!(IonSeriesOrdinal::try_immonium('Z').is_ok());
        assert!(IonSeriesOrdinal::try_immonium('a').is_err());

        // A bare `I` names no residue; a modified one carries a mod string.
        assert!(IonSeriesOrdinal::parse("I").is_err());
        assert!(matches!(
            IonSeriesOrdinal::parse("IC[Carbamidomethyl]"),
            Err(IonParsingError::UnsupportedModifiedImmonium { .. })
        ));
        // `p` is the only letter that takes no ordinal, and requires none.
        assert!(matches!(
            IonSeriesOrdinal::parse("p1"),
            Err(IonParsingError::UnexpectedOrdinal { .. })
        ));
        assert!(matches!(
            IonSeriesOrdinal::parse("y"),
            Err(IonParsingError::MissingOrdinal { .. })
        ));
        assert!(matches!(
            IonSeriesOrdinal::parse("q1"),
            Err(IonParsingError::UnsupportedFragmentType { .. })
        ));
    }
}
