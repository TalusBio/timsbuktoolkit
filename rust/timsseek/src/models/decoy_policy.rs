use serde::{
    Deserialize,
    Serialize,
};
use timsquery::models::capabilities::DecoyStrategy;

/// CLI-facing decoy *policy*: what the user asks for, independent of how the
/// arena realizes it. Resolved by [`map_decoy_strategy`] into timsquery's
/// arena-side [`DecoyStrategy`] mechanism (LazyMassShift / Passthrough / None).
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
                "Invalid decoy strategy: '{}'. Valid options: if-missing, force, never",
                s
            )),
        }
    }
}

/// The unified mass-shift offset for lazily-generated decoys (±CH2, replacing
/// the old, inconsistent 12.0/14.0 constants split across the materialized
/// branches). Lives here because it is a property of the *mapping*, not of
/// `LibCapabilities::default_diann()` alone — `Force` always resolves to it
/// regardless of what capabilities a given caller starts from.
const UNIFIED_CH2_OFFSET_DA: f64 = 14.0;
const UNIFIED_N_DECOYS: u8 = 2;

/// Map the CLI-facing decoy policy to the timsquery arena's lazy decoy
/// strategy, given whether the source library already ships its own decoys.
///
/// - `Force`: always (re)generate lazy mass-shift decoys, ignoring any
///   decoys the file already carries.
/// - `IfMissing` + no file decoys: generate lazy mass-shift decoys.
/// - `IfMissing` + file already has decoys: `Passthrough` — use the file's
///   own decoy rows as-is (no arena-side generation).
/// - `Never`: no decoy generation; library rows are used as-is.
pub fn map_decoy_strategy(policy: DecoyPolicy, has_file_decoys: bool) -> DecoyStrategy {
    match policy {
        DecoyPolicy::Force => DecoyStrategy::LazyMassShift {
            offset: UNIFIED_CH2_OFFSET_DA,
            n_decoys: UNIFIED_N_DECOYS,
        },
        DecoyPolicy::IfMissing if !has_file_decoys => DecoyStrategy::LazyMassShift {
            offset: UNIFIED_CH2_OFFSET_DA,
            n_decoys: UNIFIED_N_DECOYS,
        },
        DecoyPolicy::IfMissing => DecoyStrategy::Passthrough,
        DecoyPolicy::Never => DecoyStrategy::None,
    }
}

#[cfg(test)]
mod map_decoy_strategy_tests {
    use super::*;

    #[test]
    fn force_always_lazy_mass_shift_regardless_of_file_decoys() {
        for has_file_decoys in [false, true] {
            assert_eq!(
                map_decoy_strategy(DecoyPolicy::Force, has_file_decoys),
                DecoyStrategy::LazyMassShift {
                    offset: UNIFIED_CH2_OFFSET_DA,
                    n_decoys: UNIFIED_N_DECOYS,
                }
            );
        }
    }

    #[test]
    fn if_missing_without_file_decoys_is_lazy_mass_shift() {
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::IfMissing, false),
            DecoyStrategy::LazyMassShift {
                offset: UNIFIED_CH2_OFFSET_DA,
                n_decoys: UNIFIED_N_DECOYS,
            }
        );
    }

    #[test]
    fn if_missing_with_file_decoys_is_passthrough() {
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::IfMissing, true),
            DecoyStrategy::Passthrough
        );
    }

    #[test]
    fn never_is_none() {
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::Never, false),
            DecoyStrategy::None
        );
        assert_eq!(
            map_decoy_strategy(DecoyPolicy::Never, true),
            DecoyStrategy::None
        );
    }
}
