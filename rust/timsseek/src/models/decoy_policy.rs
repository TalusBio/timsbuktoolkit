use serde::{
    Deserialize,
    Serialize,
};
use timsquery::models::capabilities::{
    DECOY_CH2_OFFSET_DA,
    DecoyHandling,
    DecoyStrategy,
};

/// CLI-facing decoy *policy*: what the user asks for, independent of how the
/// arena realizes it. Resolved by [`map_decoy_strategy`] into timsquery's
/// arena-side [`DecoyStrategy`] mechanism (`MassShift` / `Stored`).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
#[derive(Default)]
pub enum DecoyPolicy {
    /// Generate mass-shift decoys only if library has none (default)
    #[default]
    IfMissing,

    /// Force generation: drop library decoys and regenerate mass-shift decoys
    Force,

    /// Never generate decoys, use library as-is
    Never,
}

impl DecoyPolicy {
    /// What the reader should do with decoy rows the file ships.
    ///
    /// `Force` is the only policy that drops them, and dropping them is what
    /// makes it mean anything: `map_decoy_strategy` returns `MassShift` for it,
    /// but `seal()` rewrites that to `Stored` whenever the arena holds a shipped
    /// decoy. Filtering at the reader leaves nothing to rewrite, so the
    /// mass-shift decoys the user asked for are the ones that get scored.
    pub fn decoy_handling(self) -> DecoyHandling {
        match self {
            Self::Force => DecoyHandling::Skip,
            Self::IfMissing | Self::Never => DecoyHandling::Keep,
        }
    }
}

impl std::fmt::Display for DecoyPolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecoyPolicy::IfMissing => write!(f, "if-missing"),
            DecoyPolicy::Force => write!(f, "force"),
            DecoyPolicy::Never => write!(f, "never"),
        }
    }
}

impl std::str::FromStr for DecoyPolicy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "if-missing" | "ifmissing" | "if_missing" => Ok(DecoyPolicy::IfMissing),
            "force" => Ok(DecoyPolicy::Force),
            "never" | "none" => Ok(DecoyPolicy::Never),
            _ => Err(format!(
                "Invalid decoy policy: '{}'. Valid options: if-missing, force, never",
                s
            )),
        }
    }
}

/// Map the CLI-facing decoy policy to the arena's decoy strategy, given whether
/// the source library already ships its own decoys. The ±CH2 offset comes from
/// timsquery's [`DECOY_CH2_OFFSET_DA`], so `Force` and the reader default can
/// never drift apart.
///
/// - `Force`: always derive mass-shift decoys, ignoring any the file carries
///   (the reader drops those; see [`DecoyPolicy::decoy_handling`]).
/// - `IfMissing` + no file decoys: derive mass-shift decoys.
/// - `IfMissing` + file already has decoys: `Stored`.
/// - `Never`: `Stored`.
///
/// The last two landing on the same strategy is the point: with decoys already
/// in the arena, "use them" and "generate nothing" are the same instruction.
pub fn map_decoy_strategy(policy: DecoyPolicy, has_file_decoys: bool) -> DecoyStrategy {
    let shift = DecoyStrategy::MassShift {
        offset: DECOY_CH2_OFFSET_DA,
    };
    match policy {
        DecoyPolicy::Force => shift,
        DecoyPolicy::IfMissing if !has_file_decoys => shift,
        DecoyPolicy::IfMissing | DecoyPolicy::Never => DecoyStrategy::Stored,
    }
}

#[cfg(test)]
mod map_decoy_strategy_tests {
    use super::*;

    #[test]
    fn force_always_mass_shift_regardless_of_file_decoys() {
        for has_file_decoys in [false, true] {
            assert_eq!(
                map_decoy_strategy(DecoyPolicy::Force, has_file_decoys),
                DecoyStrategy::MassShift {
                    offset: DECOY_CH2_OFFSET_DA,
                }
            );
        }
    }

    /// The half that `map_decoy_strategy` alone cannot deliver.
    ///
    /// `Force` resolves to `MassShift`, and `seal()` rewrites that to `Stored`
    /// whenever the arena holds a shipped decoy -- so asserting the mapping
    /// proves nothing about what runs. What makes `Force` mean "replace them" is
    /// the reader dropping them, leaving `seal()` nothing to rewrite.
    #[test]
    fn only_force_drops_the_decoys_a_file_ships() {
        assert_eq!(DecoyPolicy::Force.decoy_handling(), DecoyHandling::Skip);
        assert_eq!(DecoyPolicy::IfMissing.decoy_handling(), DecoyHandling::Keep);
        assert_eq!(DecoyPolicy::Never.decoy_handling(), DecoyHandling::Keep);

        // `Skip` rejects only decoys, so targets survive either way.
        assert!(!DecoyHandling::Skip.accepts(true));
        assert!(DecoyHandling::Skip.accepts(false));
        assert!(DecoyHandling::Keep.accepts(true));
    }

    #[test]
    fn if_missing_without_file_decoys_is_mass_shift() {
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::IfMissing, false),
            DecoyStrategy::MassShift {
                offset: DECOY_CH2_OFFSET_DA,
            }
        );
    }

    #[test]
    fn if_missing_with_file_decoys_is_stored() {
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::IfMissing, true),
            DecoyStrategy::Stored
        );
    }

    #[test]
    fn never_is_stored_either_way() {
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::Never, false),
            DecoyStrategy::Stored
        );
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::Never, true),
            DecoyStrategy::Stored
        );
    }
}
