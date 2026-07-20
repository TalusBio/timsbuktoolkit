use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
#[derive(Default)]
pub enum DecoyStrategy {
    /// Generate mass-shift decoys only if library has none (default)
    #[default]
    IfMissing,

    /// Force generation: drop library decoys and regenerate mass-shift decoys
    Force,

    /// Never generate decoys, use library as-is
    Never,
}

impl std::fmt::Display for DecoyStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecoyStrategy::IfMissing => write!(f, "if-missing"),
            DecoyStrategy::Force => write!(f, "force"),
            DecoyStrategy::Never => write!(f, "never"),
        }
    }
}

impl std::str::FromStr for DecoyStrategy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "if-missing" | "ifmissing" | "if_missing" => Ok(DecoyStrategy::IfMissing),
            "force" => Ok(DecoyStrategy::Force),
            "never" | "none" => Ok(DecoyStrategy::Never),
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

/// Map the CLI-facing decoy strategy to the timsquery arena's lazy decoy
/// strategy, given whether the source library already ships its own decoys.
///
/// - `Force`: always (re)generate lazy mass-shift decoys, ignoring any
///   decoys the file already carries.
/// - `IfMissing` + no file decoys: generate lazy mass-shift decoys.
/// - `IfMissing` + file already has decoys: `Passthrough` — use the file's
///   own decoy rows as-is (no arena-side generation).
/// - `Never`: no decoy generation; library rows are used as-is.
pub fn map_decoy_strategy(
    cli: DecoyStrategy,
    has_file_decoys: bool,
) -> timsquery::models::capabilities::DecoyStrategy {
    use timsquery::models::capabilities::DecoyStrategy as TqDecoyStrategy;
    match cli {
        DecoyStrategy::Force => TqDecoyStrategy::LazyMassShift {
            offset: UNIFIED_CH2_OFFSET_DA,
            n_decoys: UNIFIED_N_DECOYS,
        },
        DecoyStrategy::IfMissing if !has_file_decoys => TqDecoyStrategy::LazyMassShift {
            offset: UNIFIED_CH2_OFFSET_DA,
            n_decoys: UNIFIED_N_DECOYS,
        },
        DecoyStrategy::IfMissing => TqDecoyStrategy::Passthrough,
        DecoyStrategy::Never => TqDecoyStrategy::None,
    }
}

#[cfg(test)]
mod map_decoy_strategy_tests {
    use super::*;
    use timsquery::models::capabilities::DecoyStrategy as TqDecoyStrategy;

    #[test]
    fn force_always_lazy_mass_shift_regardless_of_file_decoys() {
        for has_file_decoys in [false, true] {
            assert_eq!(
                map_decoy_strategy(DecoyStrategy::Force, has_file_decoys),
                TqDecoyStrategy::LazyMassShift {
                    offset: UNIFIED_CH2_OFFSET_DA,
                    n_decoys: UNIFIED_N_DECOYS,
                }
            );
        }
    }

    #[test]
    fn if_missing_without_file_decoys_is_lazy_mass_shift() {
        assert_eq!(
            map_decoy_strategy(DecoyStrategy::IfMissing, false),
            TqDecoyStrategy::LazyMassShift {
                offset: UNIFIED_CH2_OFFSET_DA,
                n_decoys: UNIFIED_N_DECOYS,
            }
        );
    }

    #[test]
    fn if_missing_with_file_decoys_is_passthrough() {
        assert_eq!(
            map_decoy_strategy(DecoyStrategy::IfMissing, true),
            TqDecoyStrategy::Passthrough
        );
    }

    #[test]
    fn never_is_none() {
        assert_eq!(
            map_decoy_strategy(DecoyStrategy::Never, false),
            TqDecoyStrategy::None
        );
        assert_eq!(
            map_decoy_strategy(DecoyStrategy::Never, true),
            TqDecoyStrategy::None
        );
    }
}
