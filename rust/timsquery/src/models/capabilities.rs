#[derive(Debug, Clone, PartialEq)]
pub struct LibCapabilities {
    pub sequence_features: SeqFeatureState,
    pub fragment_features: FragmentFeatureState,
    pub isotopes: IsotopeStrategy,
    pub decoys: DecoyStrategy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFeatureState {
    Available,
    Unavailable,
}

/// Runtime reflection of whether this arena's label carries ion chemistry
/// (`FragmentLabel`). `IonAnnot` arenas => `Available`; string-labelled
/// arenas => `Unavailable`. Mirror of the `sequence_features` gate; consumed
/// later by annotation-dependent ML features.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragmentFeatureState {
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
            fragment_features: FragmentFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::LazyMassShift {
                offset: 14.0,
                n_decoys: 2,
            },
        }
    }

    /// The DIA-NN profile for EXTRACTION readers: identical to
    /// [`default_diann`](Self::default_diann) (3-isotope composition envelopes,
    /// sequence/fragment features available) but with `decoys =
    /// DecoyStrategy::None`.
    ///
    /// Decoy generation is a SCORING concern, decided by the consumer, not an
    /// assertion an extraction reader should bake into its output. timsseek
    /// stamps `caps.decoys` (via `map_decoy_strategy`) before sealing, so it is
    /// unaffected by this reader default; timsquery_cli does not override, so it
    /// correctly extracts target geometry only (no mass-shifted decoy variants).
    pub fn default_diann_no_decoys() -> Self {
        Self {
            decoys: DecoyStrategy::None,
            ..Self::default_diann()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_diann_declares_fragment_features_available() {
        assert_eq!(
            LibCapabilities::default_diann().fragment_features,
            FragmentFeatureState::Available
        );
    }
}
