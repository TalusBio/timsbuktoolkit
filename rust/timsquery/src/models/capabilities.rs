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

/// The unified ±CH2 mass-shift offset (Da) and variant count for lazily-generated
/// decoys. Single source of truth: `LibCapabilities::default_diann` and
/// timsseek's `map_decoy_strategy` both reference these, so an offset change
/// cannot drift between the reader default and the consumer mapping.
pub const DECOY_CH2_OFFSET_DA: f64 = 14.0;
pub const DECOY_N_DECOYS: u8 = 2;

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
                offset: DECOY_CH2_OFFSET_DA,
                n_decoys: DECOY_N_DECOYS,
            },
        }
    }

    /// Profile for string-labelled (unannotated) arenas: no ion chemistry, so
    /// sequence/fragment features are unavailable and no decoys are generated.
    /// Same isotope model as [`default_diann`](Self::default_diann).
    pub fn default_unlabeled() -> Self {
        Self {
            sequence_features: SeqFeatureState::Unavailable,
            fragment_features: FragmentFeatureState::Unavailable,
            ..Self::default_diann_no_decoys()
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
