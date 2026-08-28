//! Peeling an mzPAF annotation apart, and the mass-error suffix.
//!
//! The suffixes come off right to left, and that order is not arbitrary: each
//! delimiter can only appear to the right of the previous one, so peeling in
//! this sequence is unambiguous while any other order is not.
//!
//! ```text
//! y5 - H2O + 2i ^ 3 / -0.0003
//! └core┘└loss┘└iso┘└ch┘└error┘
//!    4     3     2    1     0    <- order peeled
//! ```
//!
//! `0` is [`split_mass_error`], `1..=3` are `peel_suffixes`, and `4` is
//! `IonSeriesOrdinal::parse`, which owns the spelling of the ion itself.

use crate::series::IonSeriesOrdinal;
use crate::{
    IonParsingError,
    NeutralLoss,
};

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
        .map_err(|_| IonParsingError::parse(s, "mass-error suffix is not a number"))?;
    Ok((
        head,
        Some(if is_ppm {
            MassError::Ppm(v)
        } else {
            MassError::Da(v)
        }),
    ))
}

/// The three suffixes an annotation can carry after its ion kind.
pub(crate) struct Suffixes {
    pub(crate) charge: i8,
    pub(crate) isotope: i8,
    pub(crate) loss: NeutralLoss,
}

/// Peel `^N`, `+Ni` and `-<loss>` off `value`, returning the ion core and the
/// suffixes. `value` must already have had any mass-error suffix removed.
///
/// Defaults when a suffix is absent: charge 1, isotope 0, no loss.
pub(crate) fn peel_suffixes(value: &str) -> Result<(&str, Suffixes), IonParsingError> {
    // Charge: `^N`, rightmost of the three.
    let (rest, charge) = match value.split_once('^') {
        Some((rest, charge)) => (
            rest,
            charge
                .parse::<i8>()
                .map_err(|_| IonParsingError::parse(value, "charge is not a number"))?,
        ),
        None => (value, 1),
    };

    // Isotope: `+Ni`, or bare `+i` for one.
    let (rest, isotope) = match rest.split_once('+') {
        Some((rest, adducts)) => {
            let digits = adducts.strip_suffix('i').ok_or(IonParsingError::parse(
                value,
                "only the isotope adduct '+Ni' is supported",
            ))?;
            let isotope = if digits.is_empty() {
                1
            } else {
                digits
                    .parse::<i8>()
                    .map_err(|_| IonParsingError::parse(value, "isotope count is not a number"))?
            };
            (rest, isotope)
        }
        None => (rest, 0),
    };

    // Loss: everything from the first `-`. Internal fragments spell themselves
    // `m<start>:<end>` and so never contain a `-` left of the loss.
    let (core, loss) = match rest.split_once('-') {
        Some((core, expr)) => {
            let loss = NeutralLoss::from_expression(expr)?.ok_or_else(|| {
                IonParsingError::UnsupportedNeutralLoss {
                    loss: expr.to_string(),
                }
            })?;
            (core, loss)
        }
        None => (rest, NeutralLoss::None),
    };

    Ok((
        core,
        Suffixes {
            charge,
            isotope,
            loss,
        },
    ))
}

/// Parse a whole annotation, minus any mass-error suffix, into its parts.
///
/// A comma-separated list of alternatives is NOT handled here -- that is
/// ambiguity, and resolving it needs the caller's policy. Split on `,` and
/// parse each alternative.
pub(crate) fn parse_ion(value: &str) -> Result<(IonSeriesOrdinal, Suffixes), IonParsingError> {
    let (core, suffixes) = peel_suffixes(value)?;
    Ok((IonSeriesOrdinal::parse(core)?, suffixes))
}
