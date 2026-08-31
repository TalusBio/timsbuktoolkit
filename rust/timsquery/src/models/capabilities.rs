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

/// How an arena expresses decoys: the two behaviours the index transform
/// implements, and no more.
///
/// There is no third "no decoys at all" variant. An arena that ships no decoys
/// and generates none is [`Stored`](Self::Stored) with
/// `n_stored_decoys() == 0`; that is a property of the rows, not of the
/// strategy, and reporting it belongs to whoever counts rows.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DecoyStrategy {
    /// Score the stored rows as they are. Any decoys the file shipped are the
    /// decoys; nothing is derived.
    Stored,
    /// Every stored row is a target, and the ± pair
    /// `QueryItem::variant_shift` derives from it are its decoys. Never
    /// materialized, so the arena holds targets only.
    MassShift { offset: f64 },
}

/// What a reader does with decoy rows the file ships.
///
/// Decided before reading rather than after, so a caller that is going to
/// generate its own decoys never pays to parse the file's. It also makes the
/// arena's own `MassShift -> Stored` downgrade an observation instead
/// of an override: with no shipped decoys in the arena there is nothing to
/// downgrade, so the strategy the caller asked for is the one that runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoyHandling {
    /// Read them. What a search wants unless it is replacing them.
    #[default]
    Keep,
    /// Drop them at the row level, before they reach the arena.
    Skip,
}

impl DecoyHandling {
    /// Whether a row with this decoy flag should be read.
    pub fn accepts(self, is_decoy: bool) -> bool {
        !(is_decoy && self == Self::Skip)
    }
}

/// The ±CH2 mass-shift offset (Da) for derived decoys. Every construction of
/// [`DecoyStrategy::MassShift`] uses it, including the test profiles, so a
/// change here reaches all of them.
///
/// Monoisotopic CH2: `12.0` (12-C, exact by definition) `+ 2 * 1.00782503207`
/// (1-H). Rounding it to `14.0` is 15.65 mDa light of the group it is named
/// after, or 22 ppm at m/z 700 charge 1 against a default 15 ppm tolerance,
/// which moves every decoy m/z by more than the tolerance allows.
pub const DECOY_CH2_OFFSET_DA: f64 = 14.01565006414;

/// Slots a [`DecoyStrategy::MassShift`] row expands into: the target plus the
/// `+offset`/`-offset` pair.
///
/// Not a knob. `QueryItem::variant_shift` maps variant 1 to `+offset`, variant
/// 2 to `-offset` and everything else to `0.0`, so a larger count would mint
/// slots at the target's own m/z that `is_target()` then reports as decoys --
/// target-identical scores on the decoy side of the FDR estimate.
pub const MASS_SHIFT_VARIANTS: usize = 3;

impl TargetCapabilities {
    /// The default DIA-NN `.speclib` profile: sequence features assumed
    /// available (re-gated at load), 3-isotope composition envelopes, and no
    /// decoys.
    ///
    /// Decoy generation is a scoring decision, so no constructor in this crate
    /// produces it: `DecoyStrategy::MassShift` has to be named outright, and
    /// timsseek is the only thing that names it, in `map_decoy_strategy`.
    ///
    /// A reader defaulting to decoys here makes `timsquery_cli` emit two
    /// mass-shifted variants per row into Carafe's results, which key by `id`.
    /// `test_extraction_does_not_expand_decoys` in `timsquery_cli`'s
    /// `commands.rs` is what holds that.
    pub fn default_diann() -> Self {
        Self {
            sequence_features: SeqFeatureState::Available,
            fragment_features: FragmentFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::Stored,
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

    /// The derived ±CH2 decoy profile, for exercising the index transform in
    /// this crate's own tests. Not available to production code: turning decoys
    /// on belongs to the searcher.
    #[cfg(test)]
    pub(crate) fn test_lazy_decoys() -> Self {
        Self {
            decoys: DecoyStrategy::MassShift {
                offset: DECOY_CH2_OFFSET_DA,
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
