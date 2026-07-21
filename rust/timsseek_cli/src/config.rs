use serde::{
    Deserialize,
    Serialize,
};
use timsquery::Tolerance;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};
use timsseek::DecoyPolicy;
use timsseek::scoring::CalibrationConfig;

/// Hand-authored default configuration template. Kept in sync with
/// `Config::default_config()` by the `default_template_matches_default_config`
/// test (drift = CI failure).
pub const DEFAULT_CONFIG_TOML: &str = include_str!("../assets/default_config.toml");

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub input: Option<InputConfig>,
    pub analysis: AnalysisConfig,
    #[serde(default = "CalibrationConfig::default")]
    pub calibration: CalibrationConfig,
    pub output: Option<OutputConfig>,
    #[serde(default)]
    pub staging: Option<StagingConfig>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct StagingConfig {
    #[serde(default)]
    pub tempdir_root: Option<std::path::PathBuf>,
    #[serde(default = "default_max_prefix_keys")]
    pub max_prefix_keys: usize,
    #[serde(default)]
    pub save_sidecar: bool,
    #[serde(default = "default_sweep_age")]
    pub stale_sweep_age_hours: u64,
}

fn default_max_prefix_keys() -> usize {
    256
}
fn default_sweep_age() -> u64 {
    24
}

impl Default for StagingConfig {
    fn default() -> Self {
        Self {
            tempdir_root: None,
            max_prefix_keys: 256,
            save_sidecar: false,
            stale_sweep_age_hours: 24,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum InputConfig {
    #[serde(rename = "speclib")]
    Speclib {
        #[serde(alias = "path")]
        uri: String,
    },
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct AnalysisConfig {
    #[serde(alias = "dotd_files")]
    pub raw_inputs: Option<Vec<String>>,
    pub chunk_size: usize,
    pub tolerance: Tolerance,

    #[serde(default)]
    pub decoy_strategy: DecoyPolicy,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig {
    #[serde(alias = "directory")]
    pub uri: String,
}

impl Config {
    /// Creates a default configuration with sensible defaults:
    /// - MS tolerance: 15 ppm
    /// - Mobility tolerance: 5%
    /// - Quad tolerance: 0.1 absolute
    /// - RT tolerance: unrestricted
    /// - Chunk size: 20000
    /// - Decoy strategy: if-missing
    pub fn default_config() -> Self {
        Config {
            input: None,
            analysis: AnalysisConfig {
                raw_inputs: None,
                chunk_size: 20000,
                tolerance: Tolerance {
                    ms: MzTolerance::Ppm((15.0, 15.0)),
                    mobility: MobilityTolerance::Pct((5.0, 5.0)),
                    quad: QuadTolerance::Absolute((0.1, 0.1)),
                    rt: RtTolerance::Unrestricted,
                },
                decoy_strategy: DecoyPolicy::default(),
            },
            calibration: CalibrationConfig::default(),
            output: None,
            staging: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_TOML: &str = r#"
[analysis]
chunk_size = 20000

[analysis.tolerance]
ms = { Ppm = [15.0, 15.0] }
mobility = { Pct = [5.0, 5.0] }
quad = { Absolute = [0.1, 0.1] }
rt = "Unrestricted"
"#;

    #[test]
    fn parses_minimal_toml() {
        let c: Config = toml::from_str(MINIMAL_TOML).unwrap();
        assert_eq!(c.analysis.chunk_size, 20000);
    }

    #[test]
    fn rejects_unknown_top_level_field() {
        let bad = format!("bogus_field = 123\n{MINIMAL_TOML}");
        let err = toml::from_str::<Config>(&bad).unwrap_err().to_string();
        assert!(
            err.contains("bogus_field") || err.contains("unknown field"),
            "got: {err}"
        );
    }

    #[test]
    fn rejects_unknown_nested_field() {
        let bad = MINIMAL_TOML.replace("chunk_size = 20000", "chunk_size = 20000\nchunk_siez = 1");
        let err = toml::from_str::<Config>(&bad).unwrap_err().to_string();
        assert!(
            err.contains("chunk_siez") || err.contains("unknown field"),
            "got: {err}"
        );
    }

    #[test]
    fn default_config_roundtrips_through_toml() {
        let a = Config::default_config();
        let s = toml::to_string(&a).unwrap();
        let b: Config = toml::from_str(&s).unwrap();
        assert_eq!(b.analysis.chunk_size, a.analysis.chunk_size);
    }

    /// Drift guard: the embedded TOML template must deserialize to the same
    /// Config as `default_config()`. Comparison via JSON so we don't need
    /// PartialEq on every nested type.
    #[test]
    fn default_template_matches_default_config() {
        let from_template: Config =
            toml::from_str(DEFAULT_CONFIG_TOML).expect("default template must parse");
        let a = serde_json::to_string(&from_template).unwrap();
        let b = serde_json::to_string(&Config::default_config()).unwrap();
        assert_eq!(
            a, b,
            "default_config.toml drifted from Config::default_config()"
        );
    }
}
