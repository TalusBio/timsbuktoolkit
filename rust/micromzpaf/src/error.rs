//! Why an annotation could not be parsed or represented.

use thiserror::Error;

use crate::{
    CHARGE_MAX,
    CHARGE_MIN,
    ISOTOPE_MAX,
    ISOTOPE_MIN,
};

/// Why an annotation could not be parsed or represented.
///
/// Never matched outside this crate today, so the variants exist for the message
/// a user sees when a library fails to load. Each one names the offending value
/// and, where there is one, the bound it missed -- a reader surfaces this with
/// no other context.
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
    /// Shorthand for [`Self::ParsingError`], which is built at a dozen sites.
    pub(crate) fn parse(annotation: &str, context: &'static str) -> Self {
        Self::ParsingError {
            annotation: annotation.to_string(),
            context,
        }
    }
}
