//! Why an annotation could not be parsed or represented.

use thiserror::Error;

use crate::{
    CHARGE_MAX,
    CHARGE_MIN,
    ISOTOPE_MAX,
    ISOTOPE_MIN,
};

/// Why an annotation could not be parsed or represented.
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
        annotation: String,
        context: &'static str,
    },
}

impl IonParsingError {
    /// Create a parsing error with context.
    pub(crate) fn parse(annotation: &str, context: &'static str) -> Self {
        Self::ParsingError {
            annotation: annotation.to_string(),
            context,
        }
    }
}
