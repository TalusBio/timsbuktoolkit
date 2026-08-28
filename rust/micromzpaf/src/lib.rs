//! Compact representation of fragment ion annotations for mass spectrometry.
//!
//! A spectral library carries one annotation per fragment, so this type is
//! replicated millions of times in a loaded arena. It is a packed `u32` so that
//! the annotations an mzSpecLib-shaped library needs -- neutral losses,
//! internal fragments, immonium ions -- fit *without* growing the type.
//!
//! Bolting a loss field onto a struct of fields would have cost three bytes,
//! not one: a 5-byte key paired with an `f32` pads to a 12-byte tuple, growing
//! the inline `TinyVec` storage in timsseek's `ExpectedIntensities` from 104 to
//! 156 bytes. Packed, the tuple stays 8 bytes.
//!
//! # Bit layout
//!
//! ```text
//! bit: 31 30  29            18 17      12 11     8 7      4 3     0
//!     ┌───────┬────────────────┬──────────┬────────┬────────┬───────┐
//!     │ spare │    payload     │   loss   │isotope │ charge │ kind  │
//!     │  2b   │      12b       │    6b    │   4b   │ 4b zz  │  4b   │
//!     └───────┴────────────────┴──────────┴────────┴────────┴───────┘
//! ```
//!
//! The widths and their shifts are checked by `const` assertions next to the
//! constants, so this diagram cannot drift away from the code.
//!
//! `payload` is reinterpreted per `kind` -- a tagged union inside the word:
//!
//! | kind | payload |
//! |---|---|
//! | backbone (a/b/c/d/v/w/x/y/z), unknown | ordinal, 8b (0..=255) |
//! | internal | start 6b │ end 6b (peptides to 63 residues) |
//! | immonium | residue index 5b |
//! | precursor | unused |
//!
//! `unknown` is discriminant 0 and `charge` is stored biased by one, so the
//! all-zero word is the valid annotation `?0` at charge 1. That matters
//! because `IonAnnot: Default` is forced by `tinyvec::Array` and a default can
//! reach any serde path.
//!
//! `charge` is zigzag-encoded to stay signed in 4 bits; `isotope` is unsigned
//! (see [`ISOTOPE_MIN`]). A bit field truncates rather than wrapping loudly, so
//! every constructor range-checks -- see [`IonAnnot::try_new`].
//!
//! # The mzPAF subset
//!
//! Parses: the a/b/c/d/v/w/x/y/z series, precursor (`p`), unknown (`?`),
//! internal fragments (`m<start>:<end>`), bare immonium (`IA`), charge (`^N`),
//! positive isotopes (`+Ni`), neutral losses from [`NeutralLoss`], and the
//! mass-error suffix (`/-0.0003`, `/1.2ppm`).
//!
//! Rejects, rather than coercing onto a nearby representable ion: negative
//! isotope offsets, modified immonium (`IC[Carbamidomethyl]` carries an
//! arbitrary mod string), and losses outside the [`NeutralLoss`] table.
//!
//! # Examples
//!
//! ```
//! use micromzpaf::{IonAnnot, NeutralLoss, split_mass_error};
//!
//! // A backbone ion, with charge and isotope suffixes.
//! let ion: IonAnnot = "b12+i^3".try_into().unwrap();
//! assert_eq!(ion.try_get_ordinal(), Some(12));
//! assert_eq!((ion.get_charge(), ion.get_isotope()), (3, 1));
//!
//! // A neutral loss. Two spellings of one chemical loss are one annotation,
//! // and render canonically.
//! let ion: IonAnnot = "y5-CH3SOH".try_into().unwrap();
//! assert_eq!(ion.loss(), NeutralLoss::Methanesulfenic);
//! assert_eq!(ion.to_string(), "y5-CH4OS");
//!
//! // An internal fragment spanning residues 2..=11, and a bare immonium ion.
//! // Neither sits on a ladder, so neither has an ordinal.
//! assert_eq!(IonAnnot::try_from("m2:11").unwrap().try_get_ordinal(), None);
//! assert_eq!(IonAnnot::try_from("IA").unwrap().to_string(), "IA");
//!
//! // mzSpecLib peaks carry observed m/z; the suffix recovers theoretical.
//! let (ion, err) = split_mass_error("y1/-0.0005").unwrap();
//! assert_eq!(ion, "y1");
//! assert!((err.unwrap().theoretical_from_observed(175.1184) - 175.1189).abs() < 1e-4);
//!
//! // Unrepresentable annotations fail rather than losing the detail.
//! assert!(IonAnnot::try_from("y1-HCOOH").is_err());
//! assert!(IonAnnot::try_from("IC[Carbamidomethyl]").is_err());
//! ```

pub mod loss;

pub use loss::NeutralLoss;
use serde::{
    Deserialize,
    Serialize,
};
use std::fmt::Display;
use std::hash::Hash;
use thiserror::Error;

// ── Bit layout ───────────────────────────────────────────────────────────────

const KIND_SHIFT: u32 = 0;
const KIND_BITS: u32 = 4;
const CHARGE_SHIFT: u32 = 4;
const CHARGE_BITS: u32 = 4;
const ISOTOPE_SHIFT: u32 = 8;
const ISOTOPE_BITS: u32 = 4;
const LOSS_SHIFT: u32 = 12;
const LOSS_BITS: u32 = 6;
const PAYLOAD_SHIFT: u32 = 18;
const PAYLOAD_BITS: u32 = 12;

/// Widest charge the 4-bit zigzag field holds. Observed maximum is 3.
///
/// The bias applied when packing makes `-7..=8` representable, but the range is
/// capped symmetrically: nothing has ever needed charge 8, and an asymmetric
/// public bound invites the reader to check the arithmetic rather than trust it.
pub const CHARGE_MIN: i8 = -7;
pub const CHARGE_MAX: i8 = 7;
/// Widest isotope offset the 4-bit field holds. Observed maximum is 3.
///
/// Isotope offsets are unsigned. mzPAF spells a negative offset `-Ni`, which
/// this crate does not parse, so allowing one in the field would let `Display`
/// emit `y5+-2i` -- text no other mzPAF reader accepts, through a type whose
/// wire format *is* that text. Rejecting it here keeps the invalid spelling
/// unreachable and spends the whole 4 bits on offsets that can round-trip.
pub const ISOTOPE_MIN: i8 = 0;
pub const ISOTOPE_MAX: i8 = mask(ISOTOPE_BITS) as i8;
/// Width of each internal-fragment endpoint inside `payload`.
const INTERNAL_POS_BITS: u32 = 6;
/// Width of the immonium residue index inside `payload`.
const IMMONIUM_BITS: u32 = 5;
/// Widest residue index an internal fragment endpoint holds.
pub const INTERNAL_POS_MAX: u8 = mask(INTERNAL_POS_BITS) as u8;

// The layout is otherwise enforced by prose and a diagram. These make it
// self-checking, so a field width cannot be widened without the build failing.
const _: () = assert!(
    KIND_BITS + CHARGE_BITS + ISOTOPE_BITS + LOSS_BITS + PAYLOAD_BITS <= 32,
    "the fields overflow the word"
);
const _: () = assert!(KIND_SHIFT == 0);
const _: () = assert!(CHARGE_SHIFT == KIND_SHIFT + KIND_BITS, "fields must abut");
const _: () = assert!(
    ISOTOPE_SHIFT == CHARGE_SHIFT + CHARGE_BITS,
    "fields must abut"
);
const _: () = assert!(
    LOSS_SHIFT == ISOTOPE_SHIFT + ISOTOPE_BITS,
    "fields must abut"
);
const _: () = assert!(PAYLOAD_SHIFT == LOSS_SHIFT + LOSS_BITS, "fields must abut");
const _: () = assert!(
    2 * INTERNAL_POS_BITS <= PAYLOAD_BITS,
    "two internal-fragment endpoints must fit in one payload"
);
const _: () = assert!(
    mask(IMMONIUM_BITS) >= (b'Z' - b'A') as u32,
    "the immonium field must hold every uppercase residue"
);
// The zigzag bounds are hand-derived; these are what make them checked.
const _: () = assert!(
    zigzag_charge(CHARGE_MIN) <= mask(CHARGE_BITS)
        && zigzag_charge(CHARGE_MAX) <= mask(CHARGE_BITS),
    "the charge range does not survive zigzag inside CHARGE_BITS"
);
const _: () = assert!(
    ISOTOPE_MIN >= 0 && ISOTOPE_MAX as u32 <= mask(ISOTOPE_BITS),
    "the isotope range does not fit ISOTOPE_BITS"
);
// `pack` masks the loss discriminant, so a table grown past the field would
// truncate into a *different* loss rather than fail.
const _: () = assert!(
    NeutralLoss::COUNT as u32 <= mask(LOSS_BITS),
    "the loss table has outgrown LOSS_BITS"
);

#[inline]
const fn mask(bits: u32) -> u32 {
    (1u32 << bits) - 1
}
/// Zigzag: map a small signed value onto an unsigned one without losing the
/// sign bit to the field width.
#[inline]
const fn zigzag(v: i8) -> u32 {
    (((v as i32) << 1) ^ ((v as i32) >> 31)) as u32
}
#[inline]
const fn unzigzag(u: u32) -> i8 {
    (((u >> 1) as i32) ^ -((u & 1) as i32)) as i8
}

/// Charge is stored biased by one, so the zero field decodes to charge 1.
///
/// Charge 0 is rejected by every constructor, so it is not a value the field
/// needs to represent -- and spending the zero word on it would make
/// `IonAnnot::default()` render an annotation that cannot be parsed back.
/// `CHARGE_MIN..=CHARGE_MAX` minus one still zigzags inside 4 bits.
#[inline]
const fn zigzag_charge(charge: i8) -> u32 {
    zigzag(charge - 1)
}
#[inline]
const fn unzigzag_charge(u: u32) -> i8 {
    unzigzag(u) + 1
}

/// Compact representation of fragment annotations.
///
/// A packed `u32`; see the crate docs for the bit layout.
///
/// Deliberately not `Ord`. Ordering the packed word sorts by ordinal first,
/// then loss, then isotope, then charge, then series -- an order nobody means.
/// Nothing in the workspace sorts annotations, and `KeyLike` does not require
/// it, so the trait is not offered rather than offered and meaningless.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct IonAnnot(u32);

impl Serialize for IonAnnot {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&format!("{}", self))
    }
}

/// Deserializes simple annotations for fragments.
///
/// b12^3 -> b12 charge 3 (implicit 0 isotope)
/// b12+i^3 -> b12 charge 3 isotope 1
/// b12+3i^3 -> b12 charge 3 isotope 2
/// b13 -> b13 (implicit charge 1 and isotope 0)
///
/// The wire format is the mzPAF string, not the packed word, so the bit layout
/// can change without breaking existing files.
impl<'de> Deserialize<'de> for IonAnnot {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        Self::try_from(s.as_str()).map_err(serde::de::Error::custom)
    }
}

impl IonAnnot {
    /// Build a backbone / precursor / unknown annotation.
    ///
    /// Errors when `charge` is zero or either of `charge`/`isotope` falls
    /// outside the packed field's range. That range check is load-bearing: the
    /// bit field truncates silently, so an unchecked value would corrupt the
    /// annotation rather than fail.
    pub fn try_new(
        ion_type: char,
        ordinal: Option<u8>,
        charge: i8,
        isotope: i8,
    ) -> Result<Self, IonParsingError> {
        Self::try_new_with_loss(ion_type, ordinal, charge, isotope, NeutralLoss::None)
    }

    /// As [`Self::try_new`], carrying a neutral loss.
    pub fn try_new_with_loss(
        ion_type: char,
        ordinal: Option<u8>,
        charge: i8,
        isotope: i8,
        loss: NeutralLoss,
    ) -> Result<Self, IonParsingError> {
        Self::pack(
            IonSeriesOrdinal::from_series_char(ion_type, ordinal)?,
            loss,
            charge,
            isotope,
        )
    }

    /// Build an internal fragment spanning residues `start..=end`.
    ///
    /// Endpoints are capped at [`INTERNAL_POS_MAX`], narrower than a backbone
    /// ordinal, because two of them share the 12-bit payload. Internal
    /// fragments are bounded by peptide length, so 63 is well past tryptic.
    pub fn try_new_internal(
        start: u8,
        end: u8,
        charge: i8,
        isotope: i8,
        loss: NeutralLoss,
    ) -> Result<Self, IonParsingError> {
        if start > INTERNAL_POS_MAX || end > INTERNAL_POS_MAX {
            return Err(IonParsingError::OrdinalOutOfRange {
                ordinal: start.max(end),
                series: 'm',
            });
        }
        Self::pack(
            IonSeriesOrdinal::internal { start, end },
            loss,
            charge,
            isotope,
        )
    }

    /// Build a bare immonium ion for an uppercase residue code.
    ///
    /// Modified immonium (`IC[Carbamidomethyl]`) is not representable at any
    /// field width and is rejected by the parser rather than silently losing
    /// the modification.
    pub fn try_new_immonium(
        residue: char,
        charge: i8,
        isotope: i8,
        loss: NeutralLoss,
    ) -> Result<Self, IonParsingError> {
        if !residue.is_ascii_uppercase() {
            return Err(IonParsingError::UnsupportedFragmentType {
                fragment_type: residue,
            });
        }
        Self::pack(
            IonSeriesOrdinal::immonium { residue },
            loss,
            charge,
            isotope,
        )
    }

    fn pack(
        series: IonSeriesOrdinal,
        loss: NeutralLoss,
        charge: i8,
        isotope: i8,
    ) -> Result<Self, IonParsingError> {
        if charge == 0 {
            return Err(IonParsingError::ChargeCannotBeZero);
        }
        if !(CHARGE_MIN..=CHARGE_MAX).contains(&charge) {
            return Err(IonParsingError::ChargeOutOfRange { charge });
        }
        if !(ISOTOPE_MIN..=ISOTOPE_MAX).contains(&isotope) {
            return Err(IonParsingError::IsotopeOutOfRange { isotope });
        }
        let (kind, payload) = series.to_parts();
        debug_assert!(payload <= mask(PAYLOAD_BITS), "payload overflows its field");
        Ok(IonAnnot(
            (kind << KIND_SHIFT)
                | ((zigzag_charge(charge) & mask(CHARGE_BITS)) << CHARGE_SHIFT)
                | (((isotope as u32) & mask(ISOTOPE_BITS)) << ISOTOPE_SHIFT)
                | ((loss as u32 & mask(LOSS_BITS)) << LOSS_SHIFT)
                | ((payload & mask(PAYLOAD_BITS)) << PAYLOAD_SHIFT),
        ))
    }

    #[inline]
    fn payload(self) -> u32 {
        (self.0 >> PAYLOAD_SHIFT) & mask(PAYLOAD_BITS)
    }

    #[inline]
    pub fn get_charge(&self) -> i8 {
        unzigzag_charge((self.0 >> CHARGE_SHIFT) & mask(CHARGE_BITS))
    }

    #[inline]
    pub fn get_isotope(&self) -> i8 {
        ((self.0 >> ISOTOPE_SHIFT) & mask(ISOTOPE_BITS)) as i8
    }

    /// The neutral loss this ion carries, [`NeutralLoss::None`] if it carries
    /// none.
    #[inline]
    pub fn loss(&self) -> NeutralLoss {
        NeutralLoss::from_discriminant(((self.0 >> LOSS_SHIFT) & mask(LOSS_BITS)) as u8)
    }

    /// Shift the isotope by `offset_neutrons`.
    ///
    /// Errors when the result leaves [`ISOTOPE_MIN`]..=[`ISOTOPE_MAX`]. That
    /// bound is an order of magnitude past any observed isotope offset.
    pub fn try_with_offset_neutrons(&self, offset_neutrons: i8) -> Result<Self, IonParsingError> {
        // Saturating rather than checked: the field is far narrower than `i8`,
        // so the range check below is the one that matters and it reports the
        // bound that actually applies.
        let new_isotope = self.get_isotope().saturating_add(offset_neutrons);
        if !(ISOTOPE_MIN..=ISOTOPE_MAX).contains(&new_isotope) {
            return Err(IonParsingError::IsotopeOutOfRange {
                isotope: new_isotope,
            });
        }
        Ok(IonAnnot(
            (self.0 & !(mask(ISOTOPE_BITS) << ISOTOPE_SHIFT))
                | (((new_isotope as u32) & mask(ISOTOPE_BITS)) << ISOTOPE_SHIFT),
        ))
    }

    /// The series ordinal, for the kinds that have one.
    ///
    /// `None` for precursor, unknown, internal and immonium: `?1` is a
    /// uniqueness counter, not a position in a ladder.
    pub fn try_get_ordinal(&self) -> Option<u8> {
        use IonSeriesOrdinal as S;
        match self.series_ordinal() {
            S::backbone { ordinal, .. } => Some(ordinal),
            S::unknown { .. } | S::precursor | S::internal { .. } | S::immonium { .. } => None,
        }
    }

    /// The backbone series this annotation belongs to, if any.
    ///
    /// `None` for precursor, unknown, internal and immonium ions, none of
    /// which sit on a backbone ladder.
    pub fn try_get_series(&self) -> Option<Series> {
        match self.series_ordinal() {
            IonSeriesOrdinal::backbone { series, .. } => Some(series),
            _ => None,
        }
    }

    /// The logical series-and-payload view of this annotation.
    pub fn series_ordinal(&self) -> IonSeriesOrdinal {
        IonSeriesOrdinal::from_parts((self.0 >> KIND_SHIFT) & mask(KIND_BITS), self.payload())
    }
}

/// A mass-error suffix on an mzPAF annotation: observed minus theoretical.
///
/// mzSpecLib peak lists carry *observed* m/z, so the theoretical value a
/// library reader wants is `observed - error`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MassError {
    /// Absolute, in daltons.
    Da(f64),
    /// Relative, in parts per million.
    Ppm(f64),
}

impl MassError {
    /// Recover the theoretical m/z from the observed one.
    pub fn theoretical_from_observed(&self, observed: f64) -> f64 {
        match self {
            MassError::Da(d) => observed - d,
            MassError::Ppm(p) => observed / (1.0 + p * 1e-6),
        }
    }
}

/// Split the trailing `/<error>[ppm]` off an annotation, if present.
///
/// Public because a library reader needs the error even when the ion itself is
/// unrepresentable: the error is what recovers theoretical m/z, so a peak that
/// ends up with an unknown label still gets a correct mass.
pub fn split_mass_error(s: &str) -> Result<(&str, Option<MassError>), IonParsingError> {
    let Some((head, tail)) = s.rsplit_once('/') else {
        return Ok((s, None));
    };
    let (num, is_ppm) = match tail.strip_suffix("ppm") {
        Some(n) => (n, true),
        None => (tail, false),
    };
    let v: f64 = num
        .parse()
        .map_err(|_| IonParsingError::parse(s, "Unable to parse the mass-error suffix"))?;
    Ok((
        head,
        Some(if is_ppm {
            MassError::Ppm(v)
        } else {
            MassError::Da(v)
        }),
    ))
}

impl IonAnnot {
    /// Parse the ion itself, with no mass-error suffix.
    ///
    /// A comma-separated list of alternatives is NOT handled here -- that is
    /// ambiguity, and resolving it needs the caller's policy. Split on `,` and
    /// parse each alternative.
    fn parse_ion(value: &str) -> Result<Self, IonParsingError> {
        // charge: trailing ^N
        let (rest, charge) = match value.split_once('^') {
            Some((rest, charge)) => {
                let charge = charge.parse::<i8>().map_err(|_| {
                    IonParsingError::parse(value, "Unable to parse the charge number")
                })?;
                (rest, charge)
            }
            None => (value, 1),
        };

        // isotope: +Ni. Negative isotope offsets are not supported.
        let (rest, isotope) = match rest.split_once('+') {
            Some((rest, adducts)) => {
                let adducts = adducts.strip_suffix('i').ok_or(IonParsingError::parse(
                    value,
                    "Only the isotope adduct '+Ni' is supported",
                ))?;
                let isotope = if adducts.is_empty() {
                    1
                } else {
                    adducts.parse::<i8>().map_err(|_| {
                        IonParsingError::parse(value, "Unable to parse the isotope number")
                    })?
                };
                (rest, isotope)
            }
            None => (rest, 0),
        };

        // neutral loss: everything from the first '-'. Internal fragments use
        // 'm<start>:<end>' and never contain '-' before the loss, so splitting
        // at the first '-' is unambiguous.
        let (core, loss) = match rest.split_once('-') {
            Some((core, loss_expr)) => {
                let loss = NeutralLoss::from_expression(loss_expr)?.ok_or_else(|| {
                    IonParsingError::UnsupportedNeutralLoss {
                        loss: loss_expr.to_string(),
                    }
                })?;
                (core, loss)
            }
            None => (rest, NeutralLoss::None),
        };

        // internal fragment: m<start>:<end>
        if let Some(spans) = core.strip_prefix('m')
            && let Some((a, b)) = spans.split_once(':')
        {
            let start = a.parse::<u8>().map_err(|_| {
                IonParsingError::parse(value, "Unable to parse internal-fragment start")
            })?;
            let end = b.parse::<u8>().map_err(|_| {
                IonParsingError::parse(value, "Unable to parse internal-fragment end")
            })?;
            return Self::try_new_internal(start, end, charge, isotope, loss);
        }

        // immonium: I<residue>, bare only.
        if let Some(res) = core.strip_prefix('I') {
            let mut ch = res.chars();
            return match (ch.next(), ch.next()) {
                (Some(r), None) => Self::try_new_immonium(r, charge, isotope, loss),
                // A bare `I` names no residue at all; `IC[Carbamidomethyl]` and
                // friends carry a mod string no fixed-width field can hold.
                // Both are unrepresentable, and the message says which it is.
                (None, _) => Err(IonParsingError::parse(
                    value,
                    "Immonium ion names no residue",
                )),
                _ => Err(IonParsingError::UnsupportedModifiedImmonium {
                    annotation: value.to_string(),
                }),
            };
        }

        // Backbone / precursor / unknown: a series char then an ordinal.
        let mut chars = core.chars();
        let series = chars
            .next()
            .ok_or(IonParsingError::parse(value, "Empty string"))?;
        let rest = chars.as_str();
        let ordinal =
            if rest.is_empty() {
                None
            } else {
                Some(rest.parse::<u8>().map_err(|_| {
                    IonParsingError::parse(value, "Ordinal is not a number in 0..=255")
                })?)
            };
        Self::try_new_with_loss(series, ordinal, charge, isotope, loss)
    }
}

impl TryFrom<&str> for IonAnnot {
    type Error = IonParsingError;

    /// Parses an annotation, discarding any mass-error suffix. Use
    /// [`split_mass_error`] first to keep it.
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Self::parse_ion(split_mass_error(value)?.0)
    }
}

impl Display for IonAnnot {
    /// Renders the canonical mzPAF spelling. Not byte-inverse to parsing: a
    /// non-canonical loss spelling (`-CH3SOH`) renders canonically (`-CH4OS`).
    /// The mass-error suffix is not part of the annotation and is not rendered.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.series_ordinal())?;
        write!(f, "{}", self.loss())?;

        match self.get_isotope() {
            0 => {}
            1 => write!(f, "+i")?,
            i => write!(f, "+{}i", i)?,
        }

        let charge = self.get_charge();
        if charge != 1 {
            write!(f, "^{}", charge)?;
        }

        Ok(())
    }
}

/// Why an annotation could not be parsed or represented.
///
/// Never matched outside this crate today, so the variants exist for the
/// message a user sees when a library fails to load. Each one names the
/// offending value, because the readers that surface these discard the variant
/// and show only the rendered string.
#[derive(Debug, Error)]
pub enum IonParsingError {
    #[error("Ordinal {ordinal} out of range for series '{series}'")]
    OrdinalOutOfRange { ordinal: u8, series: char },
    #[error("Series '{series}' requires an ordinal")]
    MissingOrdinal { series: char },
    #[error("Series '{series}' takes no ordinal, got {ordinal}")]
    UnexpectedOrdinal { series: char, ordinal: u8 },
    #[error("Unsupported fragment type: '{fragment_type}'")]
    UnsupportedFragmentType { fragment_type: char },
    #[error("Charge cannot be 0")]
    ChargeCannotBeZero,
    #[error("Charge {charge} outside the representable range {CHARGE_MIN}..={CHARGE_MAX}")]
    ChargeOutOfRange { charge: i8 },
    #[error(
        "Isotope offset {isotope} outside the representable range {ISOTOPE_MIN}..={ISOTOPE_MAX}"
    )]
    IsotopeOutOfRange { isotope: i8 },
    #[error("Neutral loss '{loss}' is not representable")]
    UnsupportedNeutralLoss { loss: String },
    #[error("Immonium ions must be a bare uppercase residue, got '{annotation}'")]
    UnsupportedModifiedImmonium { annotation: String },
    #[error("Ran out of distinct unknown-ion labels: the 8-bit ordinal is exhausted")]
    UnknownIonsExhausted,
    #[error("Could not parse '{annotation}': {context}")]
    ParsingError {
        /// The whole annotation, not the fragment of it that failed: a reader
        /// reports this with no row index, so the full text is the only handle
        /// the user gets on which row broke.
        annotation: String,
        context: &'static str,
    },
}

impl IonParsingError {
    /// Shorthand for [`Self::ParsingError`], which is built at thirteen sites.
    fn parse(annotation: &str, context: &'static str) -> Self {
        Self::ParsingError {
            annotation: annotation.to_string(),
            context,
        }
    }
}

/// Hands out `?1`, `?2`, ... for peaks whose annotation this crate cannot
/// represent.
///
/// Fragment labels must be unique within a precursor -- lookup is by first
/// match, so a repeated label makes every later peak carrying it unreachable.
/// A monotonic counter makes that uniqueness structural, and returning an
/// error once the 8-bit ordinal is spent keeps the overflow from being
/// something each reader has to remember to check.
///
/// Deliberately neither `Copy` nor `Clone`: a duplicated counter forks, and
/// each fork reissues labels the other already handed out -- exactly the bug
/// this type exists to prevent. Pass it by `&mut`.
#[derive(Debug, Default)]
pub struct UnknownIonCounter(u8);

impl UnknownIonCounter {
    /// The next unused unknown label at `charge`.
    pub fn next_unknown(&mut self, charge: i8) -> Result<IonAnnot, IonParsingError> {
        let ordinal = self
            .0
            .checked_add(1)
            .ok_or(IonParsingError::UnknownIonsExhausted)?;
        let annot = IonAnnot::try_new('?', Some(ordinal), charge, 0)?;
        self.0 = ordinal;
        Ok(annot)
    }
}

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
    /// Every series, in discriminant order, and parallel to the private letter
    /// table -- `series_letters_and_discriminants_agree` pins that pairing.
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
    /// The mzPAF letters, in discriminant order. The single place the
    /// letter↔discriminant pairing lives.
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

/// The logical series-and-payload view of an [`IonAnnot`].
///
/// This is a *view*: `IonAnnot` stores a packed word and reconstructs this on
/// demand. Constructing one directly does not allocate an annotation.
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
    /// Split into the `kind` discriminant and its `payload`, the two fields
    /// [`IonAnnot`] packs.
    ///
    /// This and [`Self::from_parts`] are the only place the numbering lives.
    /// The nine backbone series share one arm, so their discriminants come
    /// from [`Series`] itself and cannot drift out of step with their letters.
    const fn to_parts(self) -> (u32, u32) {
        match self {
            Self::backbone { series, ordinal } => (series as u32, ordinal as u32),
            // Discriminant 0, so the all-zero word -- `IonAnnot::default()` --
            // is the unknown ion `?0` rather than an undecodable value.
            Self::unknown { ordinal } => (0, ordinal as u32),
            Self::precursor => (10, 0),
            Self::internal { start, end } => {
                (11, (start as u32) | ((end as u32) << INTERNAL_POS_BITS))
            }
            Self::immonium { residue } => (12, (residue as u8 - b'A') as u32),
        }
    }

    /// Inverse of [`Self::to_parts`]. Total by construction: `unknown` is the
    /// catch-all discriminant, so a value this build does not recognise -- only
    /// reachable from a corrupted word -- decodes as an unknown ion instead of
    /// panicking on a path `Display` (and therefore `Serialize`) reaches.
    const fn from_parts(kind: u32, payload: u32) -> Self {
        let ordinal = payload as u8;
        match kind {
            1..=9 => Self::backbone {
                // `kind` is in range, so this is the inverse of `series as u32`.
                series: Series::ALL[kind as usize - 1],
                ordinal,
            },
            10 => Self::precursor,
            11 => Self::internal {
                start: (payload & mask(INTERNAL_POS_BITS)) as u8,
                end: ((payload >> INTERNAL_POS_BITS) & mask(INTERNAL_POS_BITS)) as u8,
            },
            12 => Self::immonium {
                residue: (b'A' + (payload & mask(IMMONIUM_BITS)) as u8) as char,
            },
            _ => Self::unknown { ordinal },
        }
    }

    /// Build the series view from an mzPAF series letter and its ordinal.
    ///
    /// `p` is the only letter that takes no ordinal; every other one requires
    /// one. Internal fragments and immonium ions are spelled differently and
    /// have their own constructors.
    fn from_series_char(c: char, ordinal: Option<u8>) -> Result<Self, IonParsingError> {
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

    fn ion(s: &str) -> IonAnnot {
        IonAnnot::try_from(s).unwrap_or_else(|e| panic!("{s:?} must parse: {e}"))
    }

    /// The size claim the design rests on: a loss, an internal-fragment span
    /// and an immonium residue all ride along without growing the tuple that
    /// timsseek stores inline.
    ///
    /// The predecessor was also 4 bytes, so `size_of` alone pins nothing --
    /// what this asserts is that the *added* fields cost nothing, which is only
    /// meaningful together with the round-trip tests proving they are really in
    /// there.
    #[test]
    fn the_added_fields_cost_no_space() {
        assert_eq!(size_of::<IonAnnot>(), 4);
        // Paired with an intensity on the scoring hot path; padding here would
        // grow the inline TinyVec storage in `ExpectedIntensities` from 104 to
        // 156 bytes at the current inline capacity of 13.
        assert_eq!(size_of::<(IonAnnot, f32)>(), 8);

        // All four of these are new capacity, and none of them widened the word.
        let loaded = IonAnnot::try_new_internal(2, 11, 3, 2, NeutralLoss::PhosphoricAcidWater)
            .expect("every field at once");
        assert_eq!(size_of_val(&loaded), 4);
        assert_eq!(loaded.loss(), NeutralLoss::PhosphoricAcidWater);
        assert_eq!(
            loaded.series_ordinal(),
            IonSeriesOrdinal::internal { start: 2, end: 11 }
        );
        assert_eq!((loaded.get_charge(), loaded.get_isotope()), (3, 2));
    }

    /// `IonAnnot: Default` is not optional -- `tinyvec::Array` requires
    /// `Item: Default`, `TimsElutionGroup` stores labels in a `TinyVec`, and
    /// timsquery's `KeyLike` propagates the bound. So a default can reach any
    /// serde path, and the zero word has to mean something.
    ///
    /// It means `?0`: charge is stored as `zigzag(charge - 1)`, so the zero
    /// field is charge 1 rather than the impossible charge 0. The default is
    /// therefore a real annotation and `Serialize`/`Deserialize` are inverses
    /// on it, instead of rendering a value no constructor accepts.
    #[test]
    fn the_default_annotation_is_a_real_annotation() {
        let d = IonAnnot::default();
        assert_eq!(d.to_string(), "?0");
        assert_eq!(d.get_charge(), 1);
        assert_eq!(serde_json::to_string(&d).unwrap(), "\"?0\"");
        assert_eq!(IonAnnot::try_from("?0").expect("the default parses"), d);
    }

    #[test]
    fn test_deserialize() {
        let cases = [
            ("b12", 'b', Some(12u8), 1i8, 0i8),
            ("b12^3", 'b', Some(12), 3, 0),
            ("y12^3", 'y', Some(12), 3, 0),
            ("b12+i^3", 'b', Some(12), 3, 1),
            ("b12+3i^3", 'b', Some(12), 3, 3),
            ("b13", 'b', Some(13), 1, 0),
            ("p^2", 'p', None, 2, 0),
            ("p", 'p', None, 1, 0),
            ("?12^2", '?', Some(12), 2, 0),
        ];

        for (input, series, ordinal, charge, isotope) in cases {
            let annot = ion(input);
            let expected = IonAnnot::try_new(series, ordinal, charge, isotope).unwrap();
            assert_eq!(annot, expected, "{input}");
            assert_eq!(annot.get_charge(), charge, "{input}");
            assert_eq!(annot.get_isotope(), isotope, "{input}");
            // Round-trips byte-identically when no loss is involved.
            assert_eq!(format!("{}", annot), input);
        }
    }

    /// Every field must survive the pack/unpack at its extremes.
    #[test]
    fn every_field_round_trips_at_its_extremes() {
        for charge in CHARGE_MIN..=CHARGE_MAX {
            if charge == 0 {
                continue;
            }
            for isotope in ISOTOPE_MIN..=ISOTOPE_MAX {
                for ordinal in [1u8, 2, 127, 255] {
                    for loss in [
                        NeutralLoss::None,
                        NeutralLoss::Water,
                        NeutralLoss::PhosphoricAcidWater,
                    ] {
                        let a =
                            IonAnnot::try_new_with_loss('y', Some(ordinal), charge, isotope, loss)
                                .expect("in range");
                        assert_eq!(a.get_charge(), charge);
                        assert_eq!(a.get_isotope(), isotope);
                        assert_eq!(a.try_get_ordinal(), Some(ordinal));
                        assert_eq!(a.loss(), loss);
                    }
                }
            }
        }
    }

    #[test]
    fn out_of_range_is_rejected_not_truncated() {
        assert!(matches!(
            IonAnnot::try_new('y', Some(1), CHARGE_MAX + 1, 0),
            Err(IonParsingError::ChargeOutOfRange { .. })
        ));
        assert!(matches!(
            IonAnnot::try_new('y', Some(1), 1, ISOTOPE_MAX + 1),
            Err(IonParsingError::IsotopeOutOfRange { .. })
        ));
        assert!(IonAnnot::try_new('y', Some(1), 0, 0).is_err());
        assert!(matches!(
            IonAnnot::try_new_internal(INTERNAL_POS_MAX + 1, 1, 1, 0, NeutralLoss::None),
            Err(IonParsingError::OrdinalOutOfRange { .. })
        ));
    }

    /// The crate docs say negative isotope offsets are unsupported. When the
    /// field held them, `Display` emitted `y5+-2i` -- not mzPAF, and reparsed
    /// happily, so a library could round-trip through serde into text no other
    /// mzPAF reader accepts. The field is unsigned so that spelling is
    /// unreachable rather than merely undocumented.
    #[test]
    fn negative_isotopes_are_unrepresentable_not_just_unparsed() {
        assert!(matches!(
            IonAnnot::try_new('y', Some(5), 1, -1),
            Err(IonParsingError::IsotopeOutOfRange { isotope: -1 })
        ));
        assert!(ion("y5").try_with_offset_neutrons(-1).is_err());
        // mzPAF's own spelling for a negative offset is still not parsed.
        assert!(IonAnnot::try_from("y5-2i").is_err());
        // And the spelling that used to leak out is not accepted either.
        assert!(IonAnnot::try_from("y5+-2i").is_err());
        // Every representable isotope renders as something that parses back.
        for isotope in ISOTOPE_MIN..=ISOTOPE_MAX {
            let a = IonAnnot::try_new('y', Some(5), 1, isotope).expect("in range");
            let text = a.to_string();
            assert!(!text.contains("+-"), "{text} is not mzPAF");
            assert_eq!(ion(&text), a, "{text}");
        }
    }

    #[test]
    fn isotope_offset_respects_the_field_bound() {
        let a = ion("y5");
        assert_eq!(a.try_with_offset_neutrons(2).unwrap().get_isotope(), 2);
        assert!(a.try_with_offset_neutrons(ISOTOPE_MAX + 1).is_err());
    }

    #[test]
    fn parses_neutral_losses() {
        let a = ion("y5-H2O");
        assert_eq!(a.loss(), NeutralLoss::Water);
        assert_eq!(a.try_get_ordinal(), Some(5));
        assert_eq!(format!("{}", a), "y5-H2O");

        // Non-canonical spelling resolves to the same annotation, and renders
        // canonically -- so this pair is equal, which is the property that
        // upholds per-precursor label uniqueness.
        assert_eq!(ion("y5-CH3SOH"), ion("y5-CH4OS"));
        assert_eq!(format!("{}", ion("y5-CH3SOH")), "y5-CH4OS");

        // Loss combines with charge and isotope.
        let b = ion("y10-NH3+i^2");
        assert_eq!(b.loss(), NeutralLoss::Ammonia);
        assert_eq!(b.get_isotope(), 1);
        assert_eq!(b.get_charge(), 2);
    }

    #[test]
    fn parses_internal_fragments() {
        let a = ion("m2:11");
        assert_eq!(
            a.series_ordinal(),
            IonSeriesOrdinal::internal { start: 2, end: 11 }
        );
        // An internal fragment has no ladder position.
        assert_eq!(a.try_get_ordinal(), None);
        assert_eq!(format!("{}", a), "m2:11");

        let b = ion("m11:12-CO");
        assert_eq!(b.loss(), NeutralLoss::CarbonMonoxide);
        assert_eq!(format!("{}", b), "m11:12-CO");
    }

    #[test]
    fn parses_bare_immonium_and_rejects_modified() {
        let a = ion("IA");
        assert_eq!(
            a.series_ordinal(),
            IonSeriesOrdinal::immonium { residue: 'A' }
        );
        assert_eq!(format!("{}", a), "IA");

        // Carries an arbitrary mod string -- no field width represents it, so
        // it must fail rather than silently degrade to a bare immonium.
        assert!(matches!(
            IonAnnot::try_from("IC[Carbamidomethyl]"),
            Err(IonParsingError::UnsupportedModifiedImmonium { .. })
        ));
    }

    #[test]
    fn parses_and_applies_the_mass_error_suffix() {
        let (rest, err) = split_mass_error("y1/-0.0005").unwrap();
        assert_eq!(ion(rest), ion("y1"));
        assert_eq!(err, Some(MassError::Da(-0.0005)));
        // theoretical = observed - error; verified against a real SpectraST
        // peak: y1 for C-terminal R, observed 175.1184, theoretical 175.11895.
        let theo = err.unwrap().theoretical_from_observed(175.1184);
        assert!((theo - 175.1189).abs() < 1e-9, "got {theo}");

        let (_, err) = split_mass_error("y6/1.2ppm").unwrap();
        assert_eq!(err, Some(MassError::Ppm(1.2)));
        let theo = err.unwrap().theoretical_from_observed(700.0);
        assert!((theo - 699.99916).abs() < 1e-4, "got {theo}");

        // Absent suffix is not an error.
        assert_eq!(split_mass_error("y6").unwrap(), ("y6", None));
        // TryFrom discards it rather than failing.
        assert_eq!(IonAnnot::try_from("y1/-0.0005").unwrap(), ion("y1"));
    }

    /// One representative of every `IonSeriesOrdinal` case: all nine backbone
    /// series (each with a distinct ordinal, so a transposition in
    /// `to_parts`/`from_parts` cannot cancel out) plus the four others, at
    /// their field boundaries.
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
            // Both endpoints at their 6-bit ceiling: the widest payload the
            // 12 bits hold, and the case a narrowed field would silently clip.
            IonSeriesOrdinal::internal {
                start: INTERNAL_POS_MAX,
                end: INTERNAL_POS_MAX,
            },
            IonSeriesOrdinal::immonium { residue: 'A' },
            // Highest residue index the 5-bit immonium field must hold.
            IonSeriesOrdinal::immonium { residue: 'Z' },
        ]);
        out
    }

    /// `to_parts` and `from_parts` are hand-written inverses. Without this,
    /// swapping two arms (`v` encoding as `w`) mislabels a whole ion series
    /// and every other test still passes.
    #[test]
    fn every_series_variant_round_trips_through_the_packed_word() {
        for series in all_series() {
            let (kind, payload) = series.to_parts();
            assert_eq!(IonSeriesOrdinal::from_parts(kind, payload), series);
            assert!(kind <= mask(KIND_BITS), "{series:?} kind overflows");
            assert!(
                payload <= mask(PAYLOAD_BITS),
                "{series:?} payload overflows"
            );
        }

        // Every backbone series must land on its own discriminant; the other
        // four are singletons and are covered by the round trip above.
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

    /// `from_parts` claims to be total, and three of the sixteen `kind`
    /// values are unassigned. Renumbering the assigned ones is safe only
    /// while the unassigned ones stay decodable, since a corrupted word
    /// reaches this on a `Display` (and so `Serialize`) path.
    #[test]
    fn every_kind_bit_pattern_decodes_without_panicking() {
        for kind in 0..=mask(KIND_BITS) {
            for payload in [0, 1, mask(PAYLOAD_BITS)] {
                // Also exercises `Display`, which is where a partial decode
                // would have panicked.
                let _ = IonSeriesOrdinal::from_parts(kind, payload).to_string();
            }
        }
    }

    /// `Display`, `from_series_char` and the parser must agree with the
    /// packing for every case.
    #[test]
    fn every_series_variant_round_trips_through_its_mzpaf_spelling() {
        for series in all_series() {
            let annot = IonAnnot::pack(series, NeutralLoss::None, 1, 0).expect("valid");
            assert_eq!(annot.series_ordinal(), series);
            let text = annot.to_string();
            assert_eq!(ion(&text).series_ordinal(), series, "{text}");
        }
    }

    /// The letters are the mzPAF spelling of the discriminants, and
    /// `from_char`/`as_char` are hand-written inverses of each other.
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

    /// An unrepresentable loss must fail loudly. Parsing `y1-HCOOH` as plain
    /// `y1` would put a loss peak's m/z on the `y1` label and collide with the
    /// real `y1`.
    #[test]
    fn unrepresentable_loss_is_rejected_not_stripped() {
        assert!(matches!(
            IonAnnot::try_from("y1-HCOOH"),
            Err(IonParsingError::UnsupportedNeutralLoss { .. })
        ));
    }
}
