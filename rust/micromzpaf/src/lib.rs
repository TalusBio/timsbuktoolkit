//! Compact representation of fragment ion annotations for mass spectrometry.
//!
//! A spectral library carries one annotation per fragment, so this type is
//! replicated millions of times in a loaded arena and compared on the scoring
//! hot path. It is therefore a packed `u32` rather than a struct of fields:
//! equality is a single word compare, and `(IonAnnot, f32)` stays 8 bytes with
//! no padding.
//!
//! # Bit layout
//!
//! ```text
//! bit: 31   30 29           18 17      12 11     8 7      4 3     0
//!     ┌───────┬────────────────┬──────────┬────────┬────────┬───────┐
//!     │ spare │    payload     │   loss   │isotope │ charge │ kind  │
//!     │  2b   │      12b       │    6b    │ 4b zz  │ 4b zz  │  4b   │
//!     └───────┴────────────────┴──────────┴────────┴────────┴───────┘
//! ```
//!
//! `payload` is reinterpreted per `kind` — a tagged union inside the word:
//!
//! | kind | payload |
//! |---|---|
//! | backbone (a/b/c/d/v/w/x/y/z), unknown | ordinal, 8b (0..=255) |
//! | internal | start 6b │ end 6b (peptides to 63 residues) |
//! | immonium | residue index 5b |
//! | precursor | unused |
//!
//! `charge` and `isotope` are zigzag-encoded so they stay signed in 4 bits.
//! Their ranges (±7) are far wider than anything observed: the HUPO-PSI corpus
//! tops out at charge 3 and isotope 3, with no negative charges at all. Because
//! the field truncates rather than wrapping loudly, every constructor
//! range-checks — see [`IonAnnot::try_new`].
//!
//! # mzPAF compliance
//!
//! Supported: the a/b/c/d/v/w/x/y/z series, precursor (`p`), unknown (`?`),
//! internal fragments (`m<start>:<end>`), bare immonium (`IA`), charge (`^N`),
//! positive isotopes (`+Ni`), neutral losses from [`NeutralLoss`], and the
//! mass-error suffix (`/-0.0003`, `/1.2ppm`).
//!
//! Not supported: negative isotope offsets, modified immonium
//! (`IC[Carbamidomethyl]` carries an arbitrary mod string), and losses outside
//! the [`NeutralLoss`] table. These are reported as errors, never coerced into
//! a nearby representable ion.
//!
//! # Examples
//!
//! ```
//! use micromzpaf::IonAnnot;
//!
//! // Parse a simple b-ion annotation
//! let ion: IonAnnot = "b12".try_into().unwrap();
//! assert_eq!(format!("{}", ion), "b12");
//!
//! // Parse with charge and isotope
//! let ion: IonAnnot = "b12+i^3".try_into().unwrap();
//! assert_eq!(ion.get_charge(), 3);
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
pub const CHARGE_MIN: i8 = -7;
pub const CHARGE_MAX: i8 = 7;
/// Widest isotope offset the 4-bit zigzag field holds. Observed maximum is 3.
pub const ISOTOPE_MIN: i8 = -7;
pub const ISOTOPE_MAX: i8 = 7;
/// Widest residue index an internal fragment endpoint holds (6 bits).
pub(crate) const INTERNAL_POS_MAX: u8 = 63;

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

/// Compact representation of fragment annotations.
///
/// A packed `u32`; see the crate docs for the bit layout. Ordering is by the
/// packed word, not field-by-field.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
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
                | ((zigzag(charge) & mask(CHARGE_BITS)) << CHARGE_SHIFT)
                | ((zigzag(isotope) & mask(ISOTOPE_BITS)) << ISOTOPE_SHIFT)
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
        unzigzag((self.0 >> CHARGE_SHIFT) & mask(CHARGE_BITS))
    }

    #[inline]
    pub fn get_isotope(&self) -> i8 {
        unzigzag((self.0 >> ISOTOPE_SHIFT) & mask(ISOTOPE_BITS))
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
                | ((zigzag(new_isotope) & mask(ISOTOPE_BITS)) << ISOTOPE_SHIFT),
        ))
    }

    /// The series ordinal, for the kinds that have one.
    ///
    /// `None` for precursor, unknown, internal and immonium: `?1` is a
    /// uniqueness counter, not a position in a ladder.
    pub fn try_get_ordinal(&self) -> Option<u8> {
        use IonSeriesOrdinal as S;
        match self.series_ordinal() {
            S::a { ordinal }
            | S::b { ordinal }
            | S::c { ordinal }
            | S::d { ordinal }
            | S::v { ordinal }
            | S::w { ordinal }
            | S::x { ordinal }
            | S::y { ordinal }
            | S::z { ordinal } => Some(ordinal),
            S::unknown { .. }
            | S::precursor
            | S::internal { .. }
            | S::immonium { .. }
            | S::None => None,
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
    let v: f64 = num.parse().map_err(|_| IonParsingError::ParsingError {
        error: s.to_string(),
        context: Some("Unable to parse the mass-error suffix"),
    })?;
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
    /// A comma-separated list of alternatives is NOT handled here — that is
    /// ambiguity, and resolving it needs the caller's policy. Split on `,` and
    /// parse each alternative.
    fn parse_ion(value: &str) -> Result<Self, IonParsingError> {
        // charge: trailing ^N
        let (rest, charge) = match value.split_once('^') {
            Some((rest, charge)) => {
                let charge = charge
                    .parse::<i8>()
                    .map_err(|_| IonParsingError::ParsingError {
                        error: value.to_string(),
                        context: Some("Unable to parse the charge number"),
                    })?;
                (rest, charge)
            }
            None => (value, 1),
        };

        // isotope: +Ni. Negative isotope offsets are not supported.
        let (rest, isotope) = match rest.split_once('+') {
            Some((rest, adducts)) => {
                let adducts = adducts
                    .strip_suffix('i')
                    .ok_or(IonParsingError::ParsingError {
                        error: adducts.to_string(),
                        context: Some("Unsupported adduct found"),
                    })?;
                let isotope = if adducts.is_empty() {
                    1
                } else {
                    adducts
                        .parse::<i8>()
                        .map_err(|_| IonParsingError::ParsingError {
                            error: value.to_string(),
                            context: Some("Unable to parse the isotope number"),
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
            let start = a.parse::<u8>().map_err(|_| IonParsingError::ParsingError {
                error: value.to_string(),
                context: Some("Unable to parse internal-fragment start"),
            })?;
            let end = b.parse::<u8>().map_err(|_| IonParsingError::ParsingError {
                error: value.to_string(),
                context: Some("Unable to parse internal-fragment end"),
            })?;
            return Self::try_new_internal(start, end, charge, isotope, loss);
        }

        // immonium: I<residue>, bare only.
        if let Some(res) = core.strip_prefix('I') {
            let mut ch = res.chars();
            return match (ch.next(), ch.next()) {
                (Some(r), None) => Self::try_new_immonium(r, charge, isotope, loss),
                // `IC[Carbamidomethyl]` and friends carry a mod string that no
                // fixed-width field can hold.
                _ => Err(IonParsingError::UnsupportedModifiedImmonium {
                    annotation: value.to_string(),
                }),
            };
        }

        // Backbone / precursor / unknown: a series char then an ordinal.
        let mut chars = core.chars();
        let series = chars.next().ok_or(IonParsingError::ParsingError {
            error: value.to_string(),
            context: Some("Empty string"),
        })?;
        let rest = chars.as_str();
        let ordinal = if rest.is_empty() {
            None
        } else {
            Some(
                rest.parse::<u8>()
                    .map_err(|e| IonParsingError::ParsingError {
                        error: format!("{rest} -> {e:?}"),
                        context: Some("Unable to parse the ordinal number"),
                    })?,
            )
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

#[derive(Debug, Error)]
pub enum IonParsingError {
    #[error("Ordinal {ordinal} out of range for series '{series}'")]
    OrdinalOutOfRange { ordinal: u8, series: char },
    #[error("Series '{series}' requires an ordinal")]
    MissingOrdinal { series: char },
    #[error("Unsupported fragment type: '{fragment_type}'")]
    UnsupportedFragmentType { fragment_type: char },
    #[error("Charge cannot be 0")]
    ChargeCannotBeZero,
    #[error("Charge {charge} outside the representable range")]
    ChargeOutOfRange { charge: i8 },
    #[error("Isotope offset {isotope} outside the representable range")]
    IsotopeOutOfRange { isotope: i8 },
    #[error("Neutral loss '{loss}' is not representable")]
    UnsupportedNeutralLoss { loss: String },
    #[error("Modified immonium ions are not representable: '{annotation}'")]
    UnsupportedModifiedImmonium { annotation: String },
    #[error("Ran out of distinct unknown-ion labels: the 8-bit ordinal is exhausted")]
    UnknownIonsExhausted,
    #[error("Parsing error: {error}{}", .context.map(|c| format!(" ({})", c)).unwrap_or_default())]
    ParsingError {
        error: String,
        context: Option<&'static str>,
    },
}

/// Hands out `?1`, `?2`, ... for peaks whose annotation this crate cannot
/// represent.
///
/// Fragment labels must be unique within a precursor — lookup is by first
/// match, so a repeated label makes every later peak carrying it unreachable.
/// A monotonic counter makes that uniqueness structural, and returning an
/// error once the 8-bit ordinal is spent keeps the overflow from being
/// something each reader has to remember to check.
#[derive(Debug, Default, Clone, Copy)]
pub struct UnknownIonCounter(u8);

impl UnknownIonCounter {
    pub fn new() -> Self {
        Self::default()
    }

    /// The next unused unknown label at `charge`.
    pub fn next(&mut self, charge: i8) -> Result<IonAnnot, IonParsingError> {
        let ordinal = self
            .0
            .checked_add(1)
            .ok_or(IonParsingError::UnknownIonsExhausted)?;
        let annot = IonAnnot::try_new('?', Some(ordinal), charge, 0)?;
        self.0 = ordinal;
        Ok(annot)
    }
}

/// The logical series-and-payload view of an [`IonAnnot`].
///
/// This is a *view*: `IonAnnot` stores a packed word and reconstructs this on
/// demand. Constructing one directly does not allocate an annotation.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Copy, Default)]
#[allow(non_camel_case_types)]
pub enum IonSeriesOrdinal {
    a {
        ordinal: u8,
    },
    b {
        ordinal: u8,
    },
    c {
        ordinal: u8,
    },
    d {
        ordinal: u8,
    },
    v {
        ordinal: u8,
    },
    w {
        ordinal: u8,
    },
    x {
        ordinal: u8,
    },
    y {
        ordinal: u8,
    },
    z {
        ordinal: u8,
    },
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

    /// This variant should not be used directly ... its mainly added to satisfy trait constraints by TinyVec
    #[default]
    None,
}

impl IonSeriesOrdinal {
    /// Split into the `kind` discriminant and its `payload`, the two fields
    /// [`IonAnnot`] packs.
    ///
    /// This and [`Self::from_parts`] are the only place the numbering lives.
    /// Both are exhaustive over this enum, so adding a variant is a compile
    /// error here rather than a silently mislabelled ion series.
    const fn to_parts(self) -> (u32, u32) {
        match self {
            Self::a { ordinal } => (1, ordinal as u32),
            Self::b { ordinal } => (2, ordinal as u32),
            Self::c { ordinal } => (3, ordinal as u32),
            Self::d { ordinal } => (4, ordinal as u32),
            Self::v { ordinal } => (5, ordinal as u32),
            Self::w { ordinal } => (6, ordinal as u32),
            Self::x { ordinal } => (7, ordinal as u32),
            Self::y { ordinal } => (8, ordinal as u32),
            Self::z { ordinal } => (9, ordinal as u32),
            Self::precursor => (10, 0),
            Self::unknown { ordinal } => (11, ordinal as u32),
            Self::internal { start, end } => (12, (start as u32) | ((end as u32) << 6)),
            Self::immonium { residue } => (13, (residue as u8 - b'A') as u32),
            Self::None => (0, 0),
        }
    }

    /// Inverse of [`Self::to_parts`]. An unrecognised discriminant decodes to
    /// [`Self::None`] rather than panicking: it can only come from a
    /// corrupted word, and the render path must stay total.
    const fn from_parts(kind: u32, payload: u32) -> Self {
        let ordinal = payload as u8;
        match kind {
            1 => Self::a { ordinal },
            2 => Self::b { ordinal },
            3 => Self::c { ordinal },
            4 => Self::d { ordinal },
            5 => Self::v { ordinal },
            6 => Self::w { ordinal },
            7 => Self::x { ordinal },
            8 => Self::y { ordinal },
            9 => Self::z { ordinal },
            10 => Self::precursor,
            11 => Self::unknown { ordinal },
            12 => Self::internal {
                start: (payload & mask(6)) as u8,
                end: ((payload >> 6) & mask(6)) as u8,
            },
            13 => Self::immonium {
                residue: (b'A' + (payload & mask(5)) as u8) as char,
            },
            _ => Self::None,
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
                Some(ordinal) => Err(IonParsingError::OrdinalOutOfRange { ordinal, series: c }),
            };
        }
        let ordinal = ordinal.ok_or(IonParsingError::MissingOrdinal { series: c })?;
        Ok(match c {
            'a' => Self::a { ordinal },
            'b' => Self::b { ordinal },
            'c' => Self::c { ordinal },
            'd' => Self::d { ordinal },
            'v' => Self::v { ordinal },
            'w' => Self::w { ordinal },
            'x' => Self::x { ordinal },
            'y' => Self::y { ordinal },
            'z' => Self::z { ordinal },
            '?' => Self::unknown { ordinal },
            _ => {
                return Err(IonParsingError::UnsupportedFragmentType { fragment_type: c });
            }
        })
    }
}

impl Display for IonSeriesOrdinal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IonSeriesOrdinal::a { ordinal } => write!(f, "a{}", ordinal),
            IonSeriesOrdinal::b { ordinal } => write!(f, "b{}", ordinal),
            IonSeriesOrdinal::c { ordinal } => write!(f, "c{}", ordinal),
            IonSeriesOrdinal::d { ordinal } => write!(f, "d{}", ordinal),
            IonSeriesOrdinal::v { ordinal } => write!(f, "v{}", ordinal),
            IonSeriesOrdinal::w { ordinal } => write!(f, "w{}", ordinal),
            IonSeriesOrdinal::x { ordinal } => write!(f, "x{}", ordinal),
            IonSeriesOrdinal::y { ordinal } => write!(f, "y{}", ordinal),
            IonSeriesOrdinal::z { ordinal } => write!(f, "z{}", ordinal),
            IonSeriesOrdinal::unknown { ordinal } => write!(f, "?{}", ordinal),
            IonSeriesOrdinal::precursor => write!(f, "p"),
            IonSeriesOrdinal::internal { start, end } => write!(f, "m{}:{}", start, end),
            IonSeriesOrdinal::immonium { residue } => write!(f, "I{}", residue),
            // Reached only via `IonAnnot::default()`, which packs to zero.
            // That `Default` is not optional: `tinyvec::Array` requires
            // `Item: Default`, `TimsElutionGroup` stores labels in a
            // `TinyVec<[T; 13]>`, and timsquery's `KeyLike` propagates the
            // bound. Since `Serialize` renders through `format!`, panicking
            // here is reachable from any serde path — so render inertly.
            IonSeriesOrdinal::None => write!(f, "?0"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ion(s: &str) -> IonAnnot {
        IonAnnot::try_from(s).unwrap_or_else(|e| panic!("{s:?} must parse: {e}"))
    }

    /// The whole point of the packed representation.
    #[test]
    fn packs_into_one_word() {
        assert_eq!(size_of::<IonAnnot>(), 4);
        // Paired with an intensity on the scoring hot path; padding here would
        // grow the inline TinyVec storage in `ExpectedIntensities`.
        assert_eq!(size_of::<(IonAnnot, f32)>(), 8);
    }

    #[test]
    fn series_ordinal_is_a_view_of_the_packed_word() {
        assert_eq!(
            ion("b12").series_ordinal(),
            IonSeriesOrdinal::b { ordinal: 12 }
        );
    }

    /// `IonAnnot: Default` packs to zero, which decodes to `IonSeriesOrdinal::None` and a
    /// charge of 0. `Serialize` renders through `format!`, so a panicking
    /// `Display` arm is reachable from any serde path — `TinyVec` alone can
    /// hand out a default.
    ///
    /// The rendered value is degenerate and deliberately does NOT round-trip:
    /// charge 0 is rejected by every constructor. Not panicking is the
    /// property being pinned.
    #[test]
    fn default_annotation_renders_instead_of_panicking() {
        let d = IonAnnot::default();
        assert_eq!(d.to_string(), "?0^0");
        assert_eq!(serde_json::to_string(&d).unwrap(), "\"?0^0\"");
        assert!(
            IonAnnot::try_from("?0^0").is_err(),
            "the default is not a valid annotation, only a printable one"
        );
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

    /// Every field must survive the pack/unpack at its extremes. The bit
    /// fields truncate rather than wrap, so a missing range check corrupts
    /// silently — this is the test that catches it.
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
        // keeps per-precursor label uniqueness honest.
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

    /// One representative of every `IonSeriesOrdinal` variant, each with a
    /// distinct payload so a transposition in `to_parts`/`from_parts` cannot
    /// cancel out.
    const ALL_SERIES: &[IonSeriesOrdinal] = &[
        IonSeriesOrdinal::a { ordinal: 1 },
        IonSeriesOrdinal::b { ordinal: 2 },
        IonSeriesOrdinal::c { ordinal: 3 },
        IonSeriesOrdinal::d { ordinal: 4 },
        IonSeriesOrdinal::v { ordinal: 5 },
        IonSeriesOrdinal::w { ordinal: 6 },
        IonSeriesOrdinal::x { ordinal: 7 },
        IonSeriesOrdinal::y { ordinal: 8 },
        IonSeriesOrdinal::z { ordinal: 9 },
        IonSeriesOrdinal::precursor,
        IonSeriesOrdinal::unknown { ordinal: 10 },
        IonSeriesOrdinal::internal { start: 2, end: 11 },
        IonSeriesOrdinal::immonium { residue: 'W' },
        IonSeriesOrdinal::None,
    ];

    /// `to_parts` and `from_parts` are hand-written inverses. Without this,
    /// swapping two arms (`v` encoding as `w`) mislabels a whole ion series
    /// and every other test still passes.
    #[test]
    fn every_series_variant_round_trips_through_the_packed_word() {
        for &series in ALL_SERIES {
            let (kind, payload) = series.to_parts();
            assert_eq!(IonSeriesOrdinal::from_parts(kind, payload), series);
            assert!(kind <= mask(KIND_BITS), "{series:?} kind overflows");
            assert!(
                payload <= mask(PAYLOAD_BITS),
                "{series:?} payload overflows"
            );
        }

        let mut kinds: Vec<u32> = ALL_SERIES.iter().map(|s| s.to_parts().0).collect();
        kinds.sort_unstable();
        kinds.dedup();
        assert_eq!(
            kinds.len(),
            ALL_SERIES.len(),
            "two series share a discriminant"
        );
    }

    /// The other three tables — `Display`, `from_series_char` and the parser —
    /// must agree with the packing for every variant, not just the handful the
    /// other tests happen to spell out.
    #[test]
    fn every_series_variant_round_trips_through_its_mzpaf_spelling() {
        for &series in ALL_SERIES {
            // `None` is the `Default` filler; it renders inertly but is not a
            // real annotation, so it has no spelling to parse back.
            if series == IonSeriesOrdinal::None {
                continue;
            }
            let annot = IonAnnot::pack(series, NeutralLoss::None, 1, 0).expect("valid");
            assert_eq!(annot.series_ordinal(), series);
            let text = annot.to_string();
            assert_eq!(ion(&text).series_ordinal(), series, "{text}");
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
