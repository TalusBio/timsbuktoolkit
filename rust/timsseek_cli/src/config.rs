use serde::{
    Deserialize,
    Serialize,
};
use timsquery::Tolerance;
use timsseek::DecoyPolicy;
use timsseek::ml::RescoreModel;
use timsseek::scoring::CalibrationConfig;

/// Hand-authored default configuration template, and the SINGLE source of the
/// built-in defaults: `Config::default_config()` parses this string, and
/// `--print-default-config` / `--write-default-config` emit it verbatim.
/// Edit the defaults here, not in Rust.
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

    #[serde(default)]
    pub rescore_model: RescoreModel,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig {
    #[serde(alias = "directory")]
    pub uri: String,
}

impl Config {
    /// The built-in default configuration, parsed from the embedded template
    /// [`DEFAULT_CONFIG_TOML`] so there is exactly one place the defaults are
    /// written down. Optional sections (`[input]`, `[output]`, `[staging]`)
    /// are commented out in the template and therefore land as `None`.
    ///
    /// Panics only on a malformed template, which is a compile-time-embedded
    /// asset and is covered by `default_template_parses_with_expected_values`.
    pub fn default_config() -> Self {
        toml::from_str(DEFAULT_CONFIG_TOML).expect("embedded default template must parse")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::models::tolerance::{
        MzTolerance,
        RtTolerance,
    };

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

    /// Template sanity guard. `default_config()` now *is* a parse of the
    /// embedded template, so a template-vs-literal comparison would be
    /// tautological — but a malformed or silently-edited template would
    /// otherwise only blow up at runtime. This pins that it parses and that
    /// the load-bearing values survive.
    ///
    /// The `[calibration]` comparison is NOT tautological: the template spells
    /// every calibration field out explicitly, so this is the only thing
    /// keeping it tied to `CalibrationConfig::default()` (the `mz_sigma`
    /// 1.5→3.0 drift was exactly this failure). Compared via JSON so we don't
    /// need `PartialEq` on the nested type.
    #[test]
    fn default_template_parses_with_expected_values() {
        let c = Config::default_config();

        assert_eq!(c.analysis.chunk_size, 20000);
        assert_eq!(c.analysis.decoy_strategy, DecoyPolicy::IfMissing);
        assert_eq!(c.analysis.rescore_model, RescoreModel::Gbm);
        assert_eq!(c.analysis.tolerance.ms, MzTolerance::Ppm((15.0, 15.0)));
        assert_eq!(c.analysis.tolerance.rt, RtTolerance::Unrestricted);

        // Every optional section must stay commented out in the template, so
        // the no-config fallback leaves them unset for the CLI flags to fill.
        assert!(c.input.is_none(), "[input] must be commented out");
        assert!(c.output.is_none(), "[output] must be commented out");
        assert!(c.staging.is_none(), "[staging] must be commented out");
        assert!(c.analysis.raw_inputs.is_none());

        assert_eq!(
            serde_json::to_string(&c.calibration).unwrap(),
            serde_json::to_string(&CalibrationConfig::default()).unwrap(),
            "[calibration] in default_config.toml drifted from CalibrationConfig::default()"
        );
    }
}
