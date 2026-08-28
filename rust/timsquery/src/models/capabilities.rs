#[derive(Debug, Clone, PartialEq)]
pub struct TargetCapabilities {
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
/// decoys. Single source of truth for timsseek's `map_decoy_strategy`, which
/// is the only production caller that turns decoys on.
pub const DECOY_CH2_OFFSET_DA: f64 = 14.0;
pub const DECOY_N_DECOYS: u8 = 2;

impl TargetCapabilities {
    /// The default DIA-NN `.speclib` profile: sequence features assumed
    /// available (re-gated at load), 3-isotope composition envelopes, and no
    /// decoys.
    ///
    /// Decoy generation is a scoring decision, so no constructor in this crate
    /// produces it -- `DecoyStrategy::LazyMassShift` has to be named outright,
    /// and timsseek is the only thing that names it (`map_decoy_strategy`,
    /// resolved once in `finalize_reference_library` before `seal`). This is
    /// not a style preference: readers used to default to decoys here, and
    /// timsquery_cli consequently emitted two mass-shifted variants per row
    /// into Carafe's results, which key by `id`. See the regression test in
    /// `timsquery_cli`'s `commands.rs`.
    pub fn default_diann() -> Self {
        Self {
            sequence_features: SeqFeatureState::Available,
            fragment_features: FragmentFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::None,
        }
    }

    /// Profile for string-labelled (unannotated) arenas: no ion chemistry, so
    /// sequence/fragment features are unavailable. Same isotope model as
    /// [`default_diann`](Self::default_diann).
    pub fn default_unlabeled() -> Self {
        Self {
            sequence_features: SeqFeatureState::Unavailable,
            fragment_features: FragmentFeatureState::Unavailable,
            ..Self::default_diann()
        }
    }

    /// The lazy ±CH2 decoy profile, for exercising the index transform in this
    /// crate's own tests. Not available to production code: turning decoys on
    /// belongs to the searcher.
    #[cfg(test)]
    pub(crate) fn test_lazy_decoys() -> Self {
        Self {
            decoys: DecoyStrategy::LazyMassShift {
                offset: DECOY_CH2_OFFSET_DA,
                n_decoys: DECOY_N_DECOYS,
            },
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
            TargetCapabilities::default_diann().fragment_features,
            FragmentFeatureState::Available
        );
    }
}
