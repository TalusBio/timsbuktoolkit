use msspeculator_inference::MsContext;
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
    #[serde(default)]
    pub library: Option<LibraryConfig>,
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

/// What acquisition the model should predict for.
///
/// One type for the two spellings a user has, so "a named setup" and "a
/// collision energy" cannot both be in force at once. They were two independent
/// `Option` fields, resolved independently against the configuration file, and a
/// `--nce` given alongside a configured setup was discarded without a word.
///
/// [`FromStr`](std::str::FromStr) so a malformed spelling fails while the
/// command line is being parsed rather than minutes into a build.
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum Acquisition {
    /// A setup name fitted into the model artifact.
    Named(String),
    /// The four acquisition factors, spelled out.
    Factors {
        instrument: String,
        detector: String,
        fragmentation: String,
        energy: Option<f32>,
    },
    /// Collision energy alone, for an otherwise unknown setup.
    Nce(f32),
}

impl Acquisition {
    /// The shape msspeculator takes. `Nce` becomes `Factors` with the three
    /// unknown factors empty, which is how the artifact's fallback is selected.
    pub fn to_ms_context(&self) -> MsContext {
        match self {
            Self::Named(name) => MsContext::Named(name.clone()),
            Self::Factors {
                instrument,
                detector,
                fragmentation,
                energy,
            } => MsContext::Factors {
                instrument: instrument.clone(),
                detector: detector.clone(),
                fragmentation: fragmentation.clone(),
                energy: *energy,
            },
            Self::Nce(energy) => MsContext::Factors {
                instrument: String::new(),
                detector: String::new(),
                fragmentation: String::new(),
                energy: Some(*energy),
            },
        }
    }
}

impl std::str::FromStr for Acquisition {
    type Err = String;

    /// A bare name, or four `::`-separated factors. The separator is what tells
    /// them apart: a setup name is a label, and nothing about it parses as a
    /// factor list.
    fn from_str(spec: &str) -> Result<Self, Self::Err> {
        if !spec.contains("::") {
            if spec.trim().is_empty() {
                return Err(
                    "expected a setup name or INSTRUMENT::DETECTOR::FRAGMENTATION::ENERGY"
                        .to_string(),
                );
            }
            return Ok(Self::Named(spec.to_string()));
        }
        let parts: Vec<&str> = spec.split("::").collect();
        let [instrument, detector, fragmentation, energy] = parts.as_slice() else {
            return Err(format!(
                "expected four ::-separated factors \
                 (INSTRUMENT::DETECTOR::FRAGMENTATION::ENERGY), got {}",
                parts.len()
            ));
        };
        // An empty energy is legal and means "the artifact's own"; a non-numeric
        // one is a typo worth reporting.
        let energy = if energy.trim().is_empty() {
            None
        } else {
            Some(
                energy
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| format!("energy {energy:?} is not a number"))?,
            )
        };
        Ok(Self::Factors {
            instrument: instrument.to_string(),
            detector: detector.to_string(),
            fragmentation: fragmentation.to_string(),
            energy,
        })
    }
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
    /// One field, so the two ways of naming an acquisition cannot both be set.
    /// Two independent fields let a configured one and a flagged other resolve
    /// separately, and then one of them silently lost.
    #[serde(default)]
    pub acquisition: Option<Acquisition>,
    #[serde(default)]
    pub chrom_context: Option<String>,
    #[serde(default)]
    pub min_intensity: Option<f64>,
    #[serde(default)]
    pub max_fragments: Option<usize>,
    #[serde(default)]
    pub decoys: Option<bool>,
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
        assert!(c.analysis.raw_inputs.is_none());

        assert_eq!(
            serde_json::to_string(&c.calibration).unwrap(),
            serde_json::to_string(&CalibrationConfig::default()).unwrap(),
            "[calibration] in default_config.toml drifted from CalibrationConfig::default()"
        );
    }
}
