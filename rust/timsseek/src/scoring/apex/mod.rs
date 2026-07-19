//! Tunable knobs for the composite apex profile.
//!
//! The apex profile (`cos^cos_pow * I^i_exp * (s_ratio + s_norm)`, optionally
//! blurred) is built in [`crate::scoring::apex_finding`]. This module holds the
//! bench-validated shipping values for its knobs.

/// Tunable apex-profile knobs. `Default` holds the bench-validated shipping
/// values.
#[derive(Debug, Clone)]
pub struct ApexConfig {
    /// Cosine exponent in the base profile: `cos^cos_pow`.
    pub cos_pow: f32,
    /// Log-intensity exponent in the base profile: `I^i_exp`.
    pub i_exp: f32,
    /// Additive scribe floor: `apex = C * (s_ratio + s_norm)`.
    pub s_ratio: f32,
    /// Gaussian-blur passes applied to the base profile (0 = none).
    pub blur_passes: usize,
}

impl Default for ApexConfig {
    fn default() -> Self {
        // Optuna-tuned on the apex_sim canonical+broad+narrow recovery suites
        // (summed pass2%), rounded to robust non-degenerate values.
        Self {
            cos_pow: 0.25,
            i_exp: 1.0,
            s_ratio: 0.25,
            blur_passes: 2,
        }
    }
}
