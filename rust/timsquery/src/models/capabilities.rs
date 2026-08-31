use serde::{
    Deserialize,
    Serialize,
};

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

/// What the caller wants done about decoys: the whole decoy decision, stated
/// once, before the file is read.
///
/// Distinct from [`DecoyStrategy`], which is the mechanism the arena ends up
/// using. A policy plus the rows the file turned out to ship determines the
/// strategy ([`Self::strategy`]), and [`TargetColumns::seal`](crate::models::TargetColumns::seal) is the only place
/// that resolves one into the other -- so `caps.decoys` is written exactly once
/// per arena, by the seal, and nothing downstream can revise it.
///
/// Stated up front rather than after the load because [`Force`](Self::Force)
/// also means "do not read the file's decoys at all" ([`Self::accepts`]): a
/// caller replacing them never pays to parse them, and there is nothing left in
/// the arena for a later pass to have to undo.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DecoyPolicy {
    /// Derive ±CH2 decoys only if the file ships none. The default.
    #[default]
    IfMissing,
    /// Drop the file's decoys and derive ±CH2 decoys instead.
    Force,
    /// Derive nothing; score the file's own rows, decoys included.
    Never,
}

impl DecoyPolicy {
    /// Whether a row with this decoy flag should be read at all.
    ///
    /// Only [`Force`](Self::Force) rejects one, and rejecting it at the row
    /// level is what makes `Force` mean anything: with no shipped decoy in the
    /// arena, [`Self::strategy`] resolves to `MassShift` and the decoys the
    /// caller asked for are the ones scored.
    pub fn accepts(self, is_decoy: bool) -> bool {
        !(is_decoy && self == Self::Force)
    }

    /// The mechanism this policy comes to, given what the file shipped.
    ///
    /// `IfMissing` with shipped decoys and `Never` land on the same strategy,
    /// which is the point: once decoys are in the arena, "use them" and "derive
    /// nothing" are the same instruction.
    pub fn strategy(self, has_shipped_decoys: bool) -> DecoyStrategy {
        let shift = DecoyStrategy::MassShift {
            offset: DECOY_CH2_OFFSET_DA,
        };
        match self {
            Self::Force => shift,
            Self::IfMissing if !has_shipped_decoys => shift,
            Self::IfMissing | Self::Never => DecoyStrategy::Stored,
        }
    }
}

impl std::fmt::Display for DecoyPolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IfMissing => write!(f, "if-missing"),
            Self::Force => write!(f, "force"),
            Self::Never => write!(f, "never"),
        }
    }
}

impl std::str::FromStr for DecoyPolicy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "if-missing" | "ifmissing" | "if_missing" => Ok(Self::IfMissing),
            "force" => Ok(Self::Force),
            "never" | "none" => Ok(Self::Never),
            _ => Err(format!(
                "Invalid decoy policy: '{s}'. Valid options: if-missing, force, never"
            )),
        }
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
    /// `decoys` is a placeholder: [`TargetColumns::seal`](crate::models::TargetColumns::seal) overwrites it from the
    /// caller's [`DecoyPolicy`], so nothing reads it before the seal and no
    /// constructor here can decide it. `Stored` rather than `MassShift` so that
    /// a hypothetical unsealed read sees target geometry only -- a decoying
    /// default would make `timsquery_cli` emit two mass-shifted variants per row
    /// into Carafe's results, which key by `id`.
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
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `Force` is the only policy that drops what a file ships, and dropping is
    /// what makes it mean "replace": with no shipped decoy left, `strategy`
    /// resolves to `MassShift` and the derived decoys are the ones scored.
    #[test]
    fn only_force_drops_the_decoys_a_file_ships() {
        assert!(!DecoyPolicy::Force.accepts(true));
        assert!(DecoyPolicy::IfMissing.accepts(true));
        assert!(DecoyPolicy::Never.accepts(true));

        // A target survives every policy.
        for p in [
            DecoyPolicy::Force,
            DecoyPolicy::IfMissing,
            DecoyPolicy::Never,
        ] {
            assert!(p.accepts(false), "{p} rejected a target");
        }
    }

    #[test]
    fn force_is_mass_shift_whatever_the_file_shipped() {
        let shift = DecoyStrategy::MassShift {
            offset: DECOY_CH2_OFFSET_DA,
        };
        for shipped in [false, true] {
            assert_eq!(DecoyPolicy::Force.strategy(shipped), shift);
        }
    }

    /// The two policies that differ only when the file ships nothing.
    #[test]
    fn if_missing_derives_only_what_the_file_lacks() {
        assert_eq!(
            DecoyPolicy::IfMissing.strategy(false),
            DecoyStrategy::MassShift {
                offset: DECOY_CH2_OFFSET_DA,
            }
        );
        assert_eq!(DecoyPolicy::IfMissing.strategy(true), DecoyStrategy::Stored);
        for shipped in [false, true] {
            assert_eq!(DecoyPolicy::Never.strategy(shipped), DecoyStrategy::Stored);
        }
    }

    /// The cli round-trips this through both a flag string and a config file.
    #[test]
    fn the_cli_spellings_round_trip() {
        for p in [
            DecoyPolicy::IfMissing,
            DecoyPolicy::Force,
            DecoyPolicy::Never,
        ] {
            assert_eq!(p.to_string().parse::<DecoyPolicy>(), Ok(p));
        }
        assert_eq!("if_missing".parse(), Ok(DecoyPolicy::IfMissing));
        assert_eq!("none".parse(), Ok(DecoyPolicy::Never));
        assert!("sometimes".parse::<DecoyPolicy>().is_err());
    }

    #[test]
    fn default_diann_declares_fragment_features_available() {
        assert_eq!(
            TargetCapabilities::default_diann().fragment_features,
            FragmentFeatureState::Available
        );
    }
}
