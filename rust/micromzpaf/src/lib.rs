//! Compact representation of fragment ion annotations for mass spectrometry.
//!
//! A packed `u32` representation of the fragment annotations used by spectral
//! libraries, including neutral losses, internal fragments and immonium ions.
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
//! The widths and shifts are checked by `const` assertions.
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
//! all-zero word is the valid annotation `?0` at charge 1.
//!
//! `charge` is zigzag-encoded to stay signed in 4 bits; `isotope` is unsigned.
//! Constructors validate both before packing.
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

const KIND_SHIFT: u32 = 0;
const CHARGE_SHIFT: u32 = 4;
const CHARGE_BITS: u32 = 4;
const ISOTOPE_SHIFT: u32 = 8;
const ISOTOPE_BITS: u32 = 4;
const LOSS_SHIFT: u32 = 12;
const LOSS_BITS: u32 = 6;
const PAYLOAD_SHIFT: u32 = 18;

/// Charge bounds representable by the packed field.
pub const CHARGE_MIN: i8 = -7;
pub const CHARGE_MAX: i8 = 8;
/// Isotope offset bounds representable by the packed field.
pub const ISOTOPE_MIN: i8 = 0;
pub const ISOTOPE_MAX: i8 = mask(ISOTOPE_BITS) as i8;

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
    zigzag_charge(CHARGE_MIN) <= mask(CHARGE_BITS)
        && zigzag_charge(CHARGE_MAX) <= mask(CHARGE_BITS),
    "the charge range does not survive zigzag inside CHARGE_BITS"
);
const _: () = assert!(
    ISOTOPE_MIN >= 0 && ISOTOPE_MAX as u32 <= mask(ISOTOPE_BITS),
    "the isotope range does not fit ISOTOPE_BITS"
);
const _: () = assert!(
    NeutralLoss::COUNT as u32 <= mask(LOSS_BITS),
    "the loss table has outgrown LOSS_BITS"
);

#[inline]
pub(crate) const fn mask(bits: u32) -> u32 {
    (1u32 << bits) - 1
}

/// Map a small signed value to an unsigned value without losing its sign.
#[inline]
const fn zigzag(v: i8) -> u32 {
    (((v as i32) << 1) ^ ((v as i32) >> 31)) as u32
}
#[inline]
const fn unzigzag(u: u32) -> i8 {
    (((u >> 1) as i32) ^ -((u & 1) as i32)) as i8
}

/// Encode charge with a one-unit bias so zero decodes to charge 1.
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
    /// Build an annotation from its parts, validating the packed fields.
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

    /// Build a backbone, precursor or unknown annotation from its series letter.
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

    /// Parse an annotation, ignoring any mass-error suffix.
    ///
    /// Use [`split_mass_error`] when the suffix is needed for m/z correction.
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

/// Generates unique unknown-ion labels within a precursor.
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

    /// Keep the annotation and annotation/intensity pair compact.
    #[test]
    fn the_added_fields_cost_no_space() {
        assert_eq!(size_of::<IonAnnot>(), 4);
        assert_eq!(size_of::<(IonAnnot, f32)>(), 8);

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

    /// The default packed word is a valid annotation and round-trips via serde.
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
        assert!(matches!(
            IonAnnot::try_from("m64:1"),
            Err(IonParsingError::OrdinalOutOfRange { .. })
        ));
    }

    #[test]
    fn negative_isotopes_are_unrepresentable_not_just_unparsed() {
        assert!(matches!(
            IonAnnot::try_new('y', Some(5), 1, -1),
            Err(IonParsingError::IsotopeOutOfRange { isotope: -1 })
        ));
        assert!(ion("y5").try_with_offset_neutrons(-1).is_err());
        assert!(IonAnnot::try_from("y5-2i").is_err());
        assert!(IonAnnot::try_from("y5+-2i").is_err());
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

        assert_eq!(ion("y5-CH3SOH"), ion("y5-CH4OS"));
        assert_eq!(format!("{}", ion("y5-CH3SOH")), "y5-CH4OS");

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
        let theo = err.unwrap().theoretical_from_observed(175.1184);
        assert!((theo - 175.1189).abs() < 1e-9, "got {theo}");

        let (_, err) = split_mass_error("y6/1.2ppm").unwrap();
        assert_eq!(err, Some(MassError::Ppm(1.2)));
        let theo = err.unwrap().theoretical_from_observed(700.0);
        assert!((theo - 699.99916).abs() < 1e-4, "got {theo}");

        assert_eq!(split_mass_error("y6").unwrap(), ("y6", None));
        assert_eq!(IonAnnot::try_from("y1/-0.0005").unwrap(), ion("y1"));
    }

    #[test]
    fn unrepresentable_loss_is_rejected_not_stripped() {
        assert!(matches!(
            IonAnnot::try_from("y1-HCOOH"),
            Err(IonParsingError::UnsupportedNeutralLoss { .. })
        ));
    }
}
