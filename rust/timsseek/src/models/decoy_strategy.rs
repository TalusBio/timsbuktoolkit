use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DecoyStrategy {
    /// Generate mass-shift decoys only if library has none (default)
    IfMissing,

    /// Force generation: drop library decoys and regenerate mass-shift decoys
    Force,

    /// Never generate decoys, use library as-is
    Never,
}

impl Default for DecoyStrategy {
    fn default() -> Self {
        DecoyStrategy::IfMissing
    }
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
