//! Compact representation of fragment ion annotations for mass spectrometry.
//!
//! This crate provides a minimal, memory-efficient encoding of fragment ion annotations
//! that pack ion series, charge state, ordinal position, and isotope offset into a
//! compact structure suitable for high-throughput proteomics workflows.
//!
//! # mzPAF Format Compliance
//!
//! This implementation supports a subset of the mzPAF annotation format:
//! - Common ion series (a, b, c, d, v, w, x, y, z, precursor)
//! - Charge states (^N notation)
//! - Positive isotope offsets (+Ni notation)
//!
//! NOTABLY boes not support:
//! - Negative isotope offsets (not yet implemented)
//! - Complex neutral losses or modifications
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

use rustyms::fragment::FragmentType;
use serde::{
    Deserialize,
    Serialize,
};
use std::fmt::Display;
use std::hash::Hash;
use std::str::FromStr;
use thiserror::Error;

/// Compact representation of fragment annotations.
///
/// This is a very compressed representation of a fragment
/// ion annotation. Essentially we are packing in 32 bytes
/// the ion series (b, y, ...), charge (+1 / -1 ...),
/// ordinal (12 in the ion series) and isotope.
///
/// It is not meant to represent all possible ions but rather have
/// a very compact representation of the common ones.
///
/// # Invariants
///
/// - **charge**: Must be non-zero (±1 to ±127). Zero charge is invalid and rejected by constructors.
/// - **ordinal**: Limited to u8 range (1-255). Peptides with >255 residues cannot be represented.
/// - **isotope**: Isotope offset relative to monoisotopic peak (M+0), range -128 to +127.
///
/// # Memory Layout
///
/// The struct fits in 32 bytes with optimal packing:
/// - `IonSeriesOrdinal`: 2 bytes (enum discriminant + u8 ordinal)
/// - `charge`: 1 byte (i8)
/// - `isotope`: 1 byte (i8)
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct IonAnnot {
    series_ordinal: IonSeriesOrdinal,
    charge: i8,
    isotope: i8,
}

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
    pub fn try_new(
        ion_type: char,
        ordinal: Option<u8>,
        charge: i8,
        isotope: i8,
    ) -> Result<Self, IonParsingError> {
        Ok(Self {
            series_ordinal: IonSeriesOrdinal::try_new(ion_type, ordinal)?,
            charge,
            isotope,
        })
    }

    pub fn from_fragment(
        frag: FragmentType,
        charge: i8,
        isotope: i8,
    ) -> Result<Self, IonParsingError> {
        Ok(Self {
            series_ordinal: IonSeriesOrdinal::try_from(frag)?,
            charge,
            isotope,
        })
    }

    pub fn terminality(&self) -> IonSeriesTerminality {
        self.series_ordinal.terminality()
    }

    pub fn try_with_offset_neutrons(&self, offset_neutrons: i8) -> Result<Self, IonParsingError> {
        let new_isotope =
            self.isotope
                .checked_add(offset_neutrons)
                .ok_or_else(|| IonParsingError::Custom {
                    error: format!(
                        "Isotope offset overflow: {} + {} exceeds i8 range",
                        self.isotope, offset_neutrons
                    ),
                })?;

        Ok(Self {
            series_ordinal: self.series_ordinal,
            charge: self.charge,
            isotope: new_isotope,
        })
    }

    pub fn get_charge(&self) -> i8 {
        self.charge
    }
}

impl TryFrom<&str> for IonAnnot {
    type Error = IonParsingError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
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

        // TODO: Implement 'negative isotopes' parsing ... right now I dont use them
        // for serialization ...
        // Note that this is not 100% compliant with mzPAF
        let (rest, isotope) = match rest.split_once('+') {
            Some((rest, adducts)) => {
                // Make sure the adduct is only +{number}?i
                let adducts = adducts
                    .strip_suffix('i')
                    .ok_or(IonParsingError::ParsingError {
                        error: adducts.to_string(),
                        context: Some("Unsupported adduct found"),
                    })?;
                // If its empty its an implicit 1 isotope
                // Since we stripped the 'i' from '+i'
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
        let series_ord = IonSeriesOrdinal::from_str(rest)?;
        if charge == 0 {
            return Err(IonParsingError::ParsingError {
                error: value.to_string(),
                context: Some("Charge cannot be 0"),
            });
        }
        Ok(Self {
            series_ordinal: series_ord,
            charge,
            isotope,
        })
    }
}

impl Display for IonAnnot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.series_ordinal)?;

        match self.isotope {
            0 => {}
            1 => write!(f, "+i")?,
            i => write!(f, "+{}i", i)?,
        }

        if self.charge != 1 {
            write!(f, "^{}", self.charge)?;
        }

        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum IonParsingError {
    #[error("Ordinal {ordinal} out of range for series '{series:?}'")]
    OrdinalOutOfRange { ordinal: i32, series: Option<char> },
    #[error("Unsupported fragment type: '{fragment_type}'")]
    UnsupportedFragmentType { fragment_type: char },
    #[error("Parsing error: {error}{}", .context.map(|c| format!(" ({})", c)).unwrap_or_default())]
    ParsingError {
        error: String,
        context: Option<&'static str>,
    },
    #[error("{error}")]
    Custom { error: String },
}

/// Refers to what terminus of the original peptide retains the
/// charge after a fragmentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IonSeriesTerminality {
    NTerm,
    CTerm,
    None,
}

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

    /// This variant should not be used directly ... its mainly added to satisfy trait constraints by TinyVec
    #[default]
    None,
}

impl IonSeriesOrdinal {
    pub fn try_new(series: char, ordinal: Option<u8>) -> Result<Self, IonParsingError> {
        let tmp = match (series, ordinal) {
            ('a', Some(ordinal)) => Self::a { ordinal },
            ('b', Some(ordinal)) => Self::b { ordinal },
            ('c', Some(ordinal)) => Self::c { ordinal },
            ('d', Some(ordinal)) => Self::d { ordinal },
            ('v', Some(ordinal)) => Self::v { ordinal },
            ('w', Some(ordinal)) => Self::w { ordinal },
            ('x', Some(ordinal)) => Self::x { ordinal },
            ('y', Some(ordinal)) => Self::y { ordinal },
            ('z', Some(ordinal)) => Self::z { ordinal },
            ('?', Some(ordinal)) => Self::unknown { ordinal },
            ('p', None) => Self::precursor,
            ('p', Some(ordinal)) => {
                return Err(IonParsingError::OrdinalOutOfRange {
                    ordinal: ordinal as i32,
                    series: Some(series),
                });
            }
            _ => {
                return Err(IonParsingError::UnsupportedFragmentType {
                    fragment_type: series,
                });
            }
        };

        Ok(tmp)
    }

    pub fn terminality(&self) -> IonSeriesTerminality {
        match self {
            IonSeriesOrdinal::a { ordinal: _ } => IonSeriesTerminality::NTerm,
            IonSeriesOrdinal::b { ordinal: _ } => IonSeriesTerminality::NTerm,
            IonSeriesOrdinal::c { ordinal: _ } => IonSeriesTerminality::NTerm,
            IonSeriesOrdinal::d { ordinal: _ } => IonSeriesTerminality::NTerm,
            IonSeriesOrdinal::v { ordinal: _ } => IonSeriesTerminality::CTerm,
            IonSeriesOrdinal::w { ordinal: _ } => IonSeriesTerminality::CTerm,
            IonSeriesOrdinal::x { ordinal: _ } => IonSeriesTerminality::CTerm,
            IonSeriesOrdinal::y { ordinal: _ } => IonSeriesTerminality::CTerm,
            IonSeriesOrdinal::z { ordinal: _ } => IonSeriesTerminality::CTerm,
            IonSeriesOrdinal::unknown { ordinal: _ } => IonSeriesTerminality::None,
            IonSeriesOrdinal::precursor => IonSeriesTerminality::None,
            IonSeriesOrdinal::None => panic!("IonSeriesOrdinal::None should not be used directly"),
        }
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
            IonSeriesOrdinal::None => panic!("IonSeriesOrdinal::None should not be used directly"),
        }
    }
}

impl FromStr for IonSeriesOrdinal {
    type Err = IonParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // "b12" split into "b" and "12"
        let (series_chunk, ordinal_chunk) = s.split_at(1);
        let series_id = series_chunk.chars().next();
        let series_ordinal = ordinal_chunk.parse::<u8>();

        match (series_id, series_ordinal) {
            (None, _) => Err(IonParsingError::ParsingError {
                error: s.to_string(),
                context: Some("Empty string"),
            }),
            (Some(x), Ok(y)) => IonSeriesOrdinal::try_new(x, Some(y)),
            (Some('p'), Err(_)) => Ok(IonSeriesOrdinal::precursor),
            (Some(_), Err(err)) => Err(IonParsingError::ParsingError {
                error: format!("{ordinal_chunk} -> {err:?}"),
                context: Some("Unable to parse the ordinal number"),
            }),
        }
    }
}

impl TryFrom<FragmentType> for IonSeriesOrdinal {
    type Error = IonParsingError;

    fn try_from(value: FragmentType) -> Result<Self, Self::Error> {
        fn try_convert_ordinal(ordinal: usize, series: char) -> Result<u8, IonParsingError> {
            ordinal
                .try_into()
                .map_err(|_| IonParsingError::OrdinalOutOfRange {
                    ordinal: ordinal as i32,
                    series: Some(series),
                })
        }
        let tmp = match value {
            FragmentType::a(ordinal, _) => IonSeriesOrdinal::a {
                ordinal: try_convert_ordinal(ordinal.series_number, 'a')?,
            },
            FragmentType::b(ordinal, _) => IonSeriesOrdinal::b {
                ordinal: try_convert_ordinal(ordinal.series_number, 'b')?,
            },
            FragmentType::c(ordinal, _) => IonSeriesOrdinal::c {
                ordinal: try_convert_ordinal(ordinal.series_number, 'c')?,
            },
            FragmentType::d(ordinal, _, _, _, _) => IonSeriesOrdinal::d {
                ordinal: try_convert_ordinal(ordinal.series_number, 'd')?,
            },
            FragmentType::v(ordinal, _, _, _) => IonSeriesOrdinal::v {
                ordinal: try_convert_ordinal(ordinal.series_number, 'v')?,
            },
            FragmentType::w(ordinal, _, _, _, _) => IonSeriesOrdinal::w {
                ordinal: try_convert_ordinal(ordinal.series_number, 'w')?,
            },
            FragmentType::x(ordinal, _) => IonSeriesOrdinal::x {
                ordinal: try_convert_ordinal(ordinal.series_number, 'x')?,
            },
            FragmentType::y(ordinal, _) => IonSeriesOrdinal::y {
                ordinal: try_convert_ordinal(ordinal.series_number, 'y')?,
            },
            FragmentType::z(ordinal, _) => IonSeriesOrdinal::z {
                ordinal: try_convert_ordinal(ordinal.series_number, 'z')?,
            },
            FragmentType::Precursor => IonSeriesOrdinal::precursor,
            _ => {
                return Err(IonParsingError::Custom {
                    error: format!("Unsupported fragment type: {value:?}"),
                });
            }
        };
        Ok(tmp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ion_series_ord_from_str() {
        let ion: IonSeriesOrdinal = IonSeriesOrdinal::from_str("b12").unwrap();
        assert_eq!(ion, IonSeriesOrdinal::b { ordinal: 12 });
    }

    #[test]
    fn test_deserialize() {
        let serde_pairs = vec![
            (
                "b12",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::b { ordinal: 12 },
                    charge: 1,
                    isotope: 0,
                },
            ),
            (
                "b12^3",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::b { ordinal: 12 },
                    charge: 3,
                    isotope: 0,
                },
            ),
            (
                "y12^3",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::y { ordinal: 12 },
                    charge: 3,
                    isotope: 0,
                },
            ),
            (
                "b12+i^3",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::b { ordinal: 12 },
                    charge: 3,
                    isotope: 1,
                },
            ),
            (
                "b12+3i^3",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::b { ordinal: 12 },
                    charge: 3,
                    isotope: 3,
                },
            ),
            (
                "b13",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::b { ordinal: 13 },
                    charge: 1,
                    isotope: 0,
                },
            ),
            (
                "p^2",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::precursor,
                    charge: 2,
                    isotope: 0,
                },
            ),
            (
                "p",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::precursor,
                    charge: 1,
                    isotope: 0,
                },
            ),
            (
                "?12^2",
                IonAnnot {
                    series_ordinal: IonSeriesOrdinal::unknown { ordinal: 12 },
                    charge: 2,
                    isotope: 0,
                },
            ),
        ];

        for (input, expected) in serde_pairs {
            let annot = IonAnnot::try_from(input).unwrap();
            assert_eq!(annot, expected);

            // Re-serialize and check that its the same
            let serialized = format!("{}", annot);
            assert_eq!(serialized, input);
        }
    }
}
