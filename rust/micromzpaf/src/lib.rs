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

pub mod error;
pub mod loss;
pub mod parse;
pub mod series;

pub use error::IonParsingError;
pub use loss::NeutralLoss;
pub use parse::{
    MassError,
    split_mass_error,
};
pub use series::{
    INTERNAL_POS_MAX,
    IonSeriesOrdinal,
    Series,
};

use serde::{
    Deserialize,
    Serialize,
};
use series::{
    KIND_BITS,
    PAYLOAD_BITS,
};
use std::fmt::Display;
use std::hash::Hash;

// ── The word ─────────────────────────────────────────────────────────────────
//
// `kind` and `payload` are `series`'s to define; this module only decides where
// they sit and what surrounds them.

const KIND_SHIFT: u32 = 0;
const CHARGE_SHIFT: u32 = 4;
const CHARGE_BITS: u32 = 4;
const ISOTOPE_SHIFT: u32 = 8;
const ISOTOPE_BITS: u32 = 4;
const LOSS_SHIFT: u32 = 12;
const LOSS_BITS: u32 = 6;
const PAYLOAD_SHIFT: u32 = 18;

/// Charge bounds. Observed maximum in the HUPO-PSI corpus is 3.
///
/// Asymmetric because the field stores `charge - 1` zigzagged: that shifts the
/// whole window up by one, so `8` fits and `-8` does not. The `const` assertion
/// below is what holds these to the field rather than to this comment --
/// `CHARGE_MIN = -8` would zigzag to 17 and truncate to 1, silently decoding as
/// charge 0.
pub const CHARGE_MIN: i8 = -7;
pub const CHARGE_MAX: i8 = 8;
/// Widest isotope offset the 4-bit field holds. Observed maximum is 3.
///
/// Isotope offsets are unsigned. mzPAF spells a negative offset `-Ni`, which
/// this crate does not parse, so allowing one in the field would let `Display`
/// emit `y5+-2i` -- text no other mzPAF reader accepts, through a type whose
/// wire format *is* that text. Rejecting it here keeps the invalid spelling
/// unreachable and spends the whole 4 bits on offsets that can round-trip.
pub const ISOTOPE_MIN: i8 = 0;
pub const ISOTOPE_MAX: i8 = mask(ISOTOPE_BITS) as i8;

// The layout is otherwise enforced by prose and a diagram. These make it
// self-checking, so a field cannot be widened or moved without the build
// failing. `series` asserts its own payload bounds.
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
pub(crate) const fn mask(bits: u32) -> u32 {
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

/// Deserializes an annotation from its mzPAF string.
///
/// The wire format is the string, not the packed word, so the bit layout can
/// change without breaking existing files.
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
    /// Build an annotation from its parts.
    ///
    /// The one constructor. `series` has already been validated against the
    /// payload width by [`IonSeriesOrdinal`]'s own constructors, so all this
    /// checks is `charge` and `isotope` -- and it must, because a bit field
    /// truncates silently, so an unchecked value would corrupt the annotation
    /// rather than fail.
    pub fn new(
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
        Ok(IonAnnot(
            (kind << KIND_SHIFT)
                | ((zigzag_charge(charge) & mask(CHARGE_BITS)) << CHARGE_SHIFT)
                | (((isotope as u32) & mask(ISOTOPE_BITS)) << ISOTOPE_SHIFT)
                | ((loss as u32 & mask(LOSS_BITS)) << LOSS_SHIFT)
                | ((payload & mask(PAYLOAD_BITS)) << PAYLOAD_SHIFT),
        ))
    }

    /// A backbone / precursor / unknown annotation from its mzPAF letter.
    ///
    /// The shape a file reader has: a series character straight out of a column.
    /// Internal fragments and immonium ions are spelled differently and go
    /// through [`IonSeriesOrdinal::try_internal`] / [`IonSeriesOrdinal::try_immonium`]
    /// and [`Self::new`].
    pub fn try_new(
        ion_type: char,
        ordinal: Option<u8>,
        charge: i8,
        isotope: i8,
    ) -> Result<Self, IonParsingError> {
        Self::new(
            IonSeriesOrdinal::from_series_char(ion_type, ordinal)?,
            NeutralLoss::None,
            charge,
            isotope,
        )
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
    /// Errors when the result leaves [`ISOTOPE_MIN`]..=[`ISOTOPE_MAX`], naming
    /// that range, since the caller cannot see the field width.
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

    /// What kind of ion this is, and whatever that kind carries.
    pub fn series_ordinal(&self) -> IonSeriesOrdinal {
        IonSeriesOrdinal::from_parts((self.0 >> KIND_SHIFT) & mask(KIND_BITS), self.payload())
    }
}

/// Guards the one-shot warning in [`IonAnnot::try_from`].
static MASS_ERROR_DISCARDED: std::sync::Once = std::sync::Once::new();

impl TryFrom<&str> for IonAnnot {
    type Error = IonParsingError;

    /// Parses an annotation, **discarding any mass-error suffix**, and warns
    /// once per process the first time it discards a real one.
    ///
    /// The suffix is a property of one observed peak, not of the ion: two peaks
    /// annotated `b12` in different spectra carry different errors, so keeping
    /// it on the annotation would make two `b12`s unequal and break the
    /// per-precursor label uniqueness the whole crate keys on.
    ///
    /// A caller that needs the error wants it for the *m/z*, not the label:
    /// call [`split_mass_error`] first, recover the theoretical m/z with
    /// [`MassError::theoretical_from_observed`], and store that. Nothing is
    /// lost that way -- only the residual, which belongs to the measurement.
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let (ion, mass_error) = split_mass_error(value)?;
        if mass_error.is_some() {
            MASS_ERROR_DISCARDED.call_once(|| {
                tracing::warn!(
                    "annotation {value:?} carries a mass-error suffix, which \
                     `IonAnnot` does not store; the m/z used is whatever the \
                     caller supplied. Use `split_mass_error` to recover the \
                     theoretical m/z. Warning once per process."
                );
            });
        }
        let (series, suffixes) = parse::parse_ion(ion)?;
        Self::new(series, suffixes.loss, suffixes.charge, suffixes.isotope)
    }
}

impl Display for IonAnnot {
    /// Renders the canonical mzPAF spelling, the exact inverse of parsing except
    /// that a non-canonical loss spelling (`-CH3SOH`) renders canonically
    /// (`-CH4OS`). The mass-error suffix is not part of the annotation.
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
        let loaded = IonAnnot::new(
            IonSeriesOrdinal::try_internal(2, 11).expect("in range"),
            NeutralLoss::PhosphoricAcidWater,
            3,
            2,
        )
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
                        let series = IonSeriesOrdinal::backbone {
                            series: Series::y,
                            ordinal,
                        };
                        let a = IonAnnot::new(series, loss, charge, isotope).expect("in range");
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
        // The payload bounds are `series`'s to enforce; see
        // `constructors_reject_what_the_payload_cannot_hold` there. This pins
        // that an out-of-range span cannot reach a word through the parser.
        assert!(matches!(
            IonAnnot::try_from("m64:1"),
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
