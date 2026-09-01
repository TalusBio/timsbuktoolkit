use serde::{
    Deserialize,
    Serialize,
};
use timsquery::Tolerance;
use timsseek::DecoyPolicy;
use timsseek::ml::RescoreModel;
use timsseek::scoring::CalibrationConfig;
use tracing::info;

use crate::errors;

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
    #[serde(default)]
    pub library: Option<LibraryConfig>,
    #[serde(default)]
    pub sequences: Option<SequencesConfig>,
}

/// The slice of a configuration file `build-library` reads.
///
/// A separate type from [`Config`], which requires `[analysis]`: a build has no
/// analysis, so a library-only file has to parse without one.
///
/// Unknown top-level keys are tolerated, which is what lets one file serve both
/// a build and the searches that read what it built. `[input]`, `[analysis]` and
/// the rest describe a search and are ignored. [`LibraryConfig`] still denies
/// unknown fields, so a typo *inside* `[library]` is an error.
#[derive(Debug, Deserialize, Clone, Default)]
pub struct BuildConfig {
    #[serde(default)]
    pub library: Option<LibraryConfig>,
}

/// How to predict a library.
///
/// Every field is optional and falls through to msspeculator's own default, so
/// the section can be omitted entirely and a build flag can beat any single
/// field without the others having to be spelled. This project owns the flag
/// surface and none of the digestion, modification or prediction logic.
#[derive(Debug, Serialize, Deserialize, Clone, Default, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct LibraryConfig {
    #[serde(default)]
    pub model: Option<String>,
    #[serde(default)]
    pub missed_cleavages: Option<usize>,
    #[serde(default)]
    pub min_length: Option<usize>,
    #[serde(default)]
    pub max_length: Option<usize>,
    #[serde(default)]
    pub min_charge: Option<i64>,
    #[serde(default)]
    pub max_charge: Option<i64>,
    #[serde(default)]
    pub fixed_mods: Option<Vec<String>>,
    #[serde(default)]
    pub variable_mods: Option<Vec<String>>,
    #[serde(default)]
    pub max_variable_mods: Option<usize>,
    #[serde(default)]
    pub min_intensity: Option<f64>,
    #[serde(default)]
    pub max_fragments: Option<usize>,
    #[serde(default)]
    pub decoys: Option<bool>,
}

/// The sequence database a run was given.
///
/// Its own section rather than a field of [`LibraryConfig`], which says how to
/// predict a library: a run that predicts nothing still wants the sequences, and
/// a library is not a substitute for the database it was digested from.
///
/// A path rather than a URI, because nothing stages a FASTA from a remote store.
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SequencesConfig {
    pub fasta: std::path::PathBuf,
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
        /// A second library, used only to fit the retention-time calibration.
        ///
        /// Lives here rather than staying a command-line-only flag because it
        /// changes what a run produces, and `config_used.json` is the artifact a
        /// run is reproduced from. A flag that never reached the config was a
        /// flag the artifact could not record.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        calib_uri: Option<String>,
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
    /// asset and is covered by `default_template_parses`.
    pub fn default_config() -> Self {
        toml::from_str(DEFAULT_CONFIG_TOML).expect("embedded default template must parse")
    }
}

/// Read a configuration file, or the built-in defaults when none was named.
///
/// TOML and JSON are both accepted, sniffed by extension.
pub(crate) fn load_config(path: Option<&std::path::Path>) -> Result<Config, errors::CliError> {
    let Some(config_path) = path else {
        info!("No config file provided, using default configuration");
        return Ok(Config::default_config());
    };
    parse_config_file(config_path)
}

/// Parse a configuration file into whichever view the caller needs.
///
/// TOML and JSON are both accepted, sniffed by extension. Generic over the
/// target so a build reads only `[library]` while a search reads the whole
/// thing; see [`BuildConfig`].
pub(crate) fn parse_config_file<T: serde::de::DeserializeOwned>(
    config_path: &std::path::Path,
) -> Result<T, errors::CliError> {
    let text = std::fs::read_to_string(config_path).map_err(|e| errors::CliError::Io {
        source: e.to_string(),
        path: Some(config_path.to_string_lossy().to_string()),
    })?;
    let is_toml = config_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("toml"))
        .unwrap_or(false);
    let parsed: Result<T, String> = if is_toml {
        toml::from_str(&text).map_err(|e| e.to_string())
    } else {
        serde_json::from_str(&text).map_err(|e| e.to_string())
    };
    parsed.map_err(|e| errors::CliError::ParseError {
        msg: format!(
            "Failed to parse config file {}: {e}\n\n\
             Run `timsseek --print-default-config` for a reference template, \
             or `--write-default-config <path>` to drop one to disk.\n\n\
             Reference default:\n```toml\n{}```\n",
            config_path.display(),
            DEFAULT_CONFIG_TOML,
        ),
    })
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
        assert_eq!(a.analysis.rescore_model, RescoreModel::Mlp);
        let s = toml::to_string(&a).unwrap();
        let b: Config = toml::from_str(&s).unwrap();
        assert_eq!(b.analysis.chunk_size, a.analysis.chunk_size);
    }

    /// Serde and clap derive names independently, so verify their spellings.
    #[test]
    fn rescore_model_spellings_agree_between_cli_and_toml() {
        use crate::cli::CliRescoreModel;
        use clap::ValueEnum;

        for cli_variant in CliRescoreModel::value_variants() {
            let model: RescoreModel = (*cli_variant).into();
            assert_eq!(CliRescoreModel::from(model), *cli_variant);

            let cli_name = cli_variant
                .to_possible_value()
                .expect("no rescore model variant is #[value(skip)]ed")
                .get_name()
                .to_string();

            let json = serde_json::to_string(&model).unwrap();
            let toml_name = json.trim_matches('"');

            assert_eq!(
                cli_name, toml_name,
                "`--rescore-model {cli_name}` and `rescore_model = \"{toml_name}\"` name \
                 {model:?} differently"
            );

            let src = MINIMAL_TOML.replace(
                "chunk_size = 20000",
                &format!("chunk_size = 20000\nrescore_model = \"{toml_name}\""),
            );
            let c: Config = toml::from_str(&src)
                .unwrap_or_else(|e| panic!("rescore_model = \"{toml_name}\" must parse: {e}"));
            assert_eq!(c.analysis.rescore_model, model);
        }
    }

    /// Template sanity guard. `default_config()` now *is* a parse of the
    /// embedded template, so a template-vs-literal comparison would be
    /// tautological -- but a malformed or silently-edited template would
    /// otherwise only blow up at runtime.
    ///
    /// No literal values are pinned here -- those live in the template. The
    /// `[calibration]` comparison is not a value check either: the template
    /// spells every calibration field out, so this is the only thing tying it
    /// to `CalibrationConfig::default()` (the `mz_sigma` 1.5→3.0 drift was
    /// exactly this failure). JSON so we don't need `PartialEq` on it.
    #[test]
    fn default_template_parses() {
        let c = Config::default_config();

        // Optional sections stay commented out, so the no-config fallback
        // leaves them for the CLI flags to fill.
        assert!(c.input.is_none(), "[input] must be commented out");
        assert!(c.output.is_none(), "[output] must be commented out");
        assert!(c.staging.is_none(), "[staging] must be commented out");
        assert!(c.library.is_none(), "[library] must be commented out");
        assert!(c.sequences.is_none(), "[sequences] must be commented out");
        assert!(c.analysis.raw_inputs.is_none());

        assert_eq!(
            serde_json::to_string(&c.calibration).unwrap(),
            serde_json::to_string(&CalibrationConfig::default()).unwrap(),
            "[calibration] in default_config.toml drifted from CalibrationConfig::default()"
        );
    }
}
