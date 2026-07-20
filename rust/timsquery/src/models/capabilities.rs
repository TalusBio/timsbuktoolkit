#[derive(Debug, Clone, PartialEq)]
pub struct LibCapabilities {
    pub sequence_features: SeqFeatureState,
    pub isotopes: IsotopeStrategy,
    pub decoys: DecoyStrategy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFeatureState {
    Available,
    Unavailable,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IsotopeStrategy {
    /// Per-peptide: C/S countable -> composition envelope; else -> averagine.
    FromComposition { n_isotopes: u8 },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DecoyStrategy {
    LazyMassShift { offset: f64, n_decoys: u8 },
    Passthrough,
    None,
}

impl LibCapabilities {
    /// The default DIA-NN `.speclib` profile: sequence features assumed
    /// available (re-gated at load), 3-isotope composition envelopes, and
    /// lazily-generated ±CH2 decoys.
    pub fn default_diann() -> Self {
        Self {
            sequence_features: SeqFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::LazyMassShift {
                offset: 14.0,
                n_decoys: 2,
            },
        }
    }
}
