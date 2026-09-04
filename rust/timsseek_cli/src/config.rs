use serde::{
    Deserialize,
    Serialize,
};
use timscentroid::centroiding::{
    MzTolerance,
    WindowCap,
};
use timsquery::{
    CentroidingConfig,
    IndexingCentroidingConfig,
    Tolerance,
};
use timsseek::ml::RescoreModel;
use timsseek::scoring::CalibrationConfig;
use timsseek::{
    DecoyPolicy,
    UnannotatedPeaks,
};
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
    /// Optional library used only to fit retention-time calibration.
    /// Independent of whether the searched library is loaded or predicted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub calibration_library: Option<LibraryInputConfig>,
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
    /// Centroiding overrides for an index built from raw. A cached `.idx`
    /// ignores it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indexing: Option<IndexingConfig>,
    /// What produced the library a run predicted, as msspeculator reports it:
    /// the model and the digests of what went into it.
    ///
    /// A record rather than a setting. Nothing reads it back and no
    /// configuration file needs to carry it, but `config_used.json` is
    /// serialized `Config`, so a run whose library exists only in memory has
    /// nowhere else to leave its provenance -- and this type denies unknown
    /// fields, so the artifact would not parse again without the field being
    /// declared.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub library_provenance: Option<serde_json::Value>,
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

/// Per-MS-level centroiding overrides on top of
/// [`IndexingCentroidingConfig::default`]. Every field is optional.
#[derive(Debug, Serialize, Deserialize, Clone, Default, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct IndexingConfig {
    #[serde(default)]
    pub ms1: CentroidingOverride,
    #[serde(default)]
    pub ms2: CentroidingOverride,
}

#[derive(Debug, Serialize, Deserialize, Clone, Default, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CentroidingOverride {
    #[serde(default)]
    pub max_peaks: Option<usize>,
    /// `{ Bins = n }` or `{ Ppm = x }`.
    #[serde(default)]
    pub mz_tol: Option<MzTolerance>,
    #[serde(default)]
    pub im_pct_tol: Option<f64>,
    #[serde(default)]
    pub early_stop_iterations: Option<u32>,
    /// `{ max_peaks = N, window_da = W }`. `Some` replaces the default cap.
    #[serde(default)]
    pub window_cap: Option<WindowCap>,
}

impl CentroidingOverride {
    fn apply(&self, base: CentroidingConfig) -> CentroidingConfig {
        let Self {
            max_peaks,
            mz_tol,
            im_pct_tol,
            early_stop_iterations,
            window_cap,
        } = self;
        CentroidingConfig {
            max_peaks: max_peaks.unwrap_or(base.max_peaks),
            mz_tol: mz_tol.unwrap_or(base.mz_tol),
            im_pct_tol: im_pct_tol.unwrap_or(base.im_pct_tol),
            early_stop_iterations: early_stop_iterations.unwrap_or(base.early_stop_iterations),
            window_cap: window_cap.or(base.window_cap),
        }
    }
}

impl IndexingConfig {
    pub fn resolve(&self) -> IndexingCentroidingConfig {
        let base = IndexingCentroidingConfig::default();
        IndexingCentroidingConfig {
            ms1: self.ms1.apply(base.ms1),
            ms2: self.ms2.apply(base.ms2),
        }
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
    },
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct LibraryInputConfig {
    #[serde(alias = "path")]
    pub uri: String,
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
    pub unannotated_peaks: UnannotatedPeaks,

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
    fn indexing_override_touches_only_the_named_field() {
        let toml = format!(
            "{MINIMAL_TOML}\n[indexing.ms2]\nwindow_cap = {{ max_peaks = 1000, window_da = 100.0 }}\n"
        );
        let cfg: Config = toml::from_str(&toml).unwrap();
        let resolved = cfg.indexing.as_ref().unwrap().resolve();
        let default = IndexingCentroidingConfig::default();
        assert_eq!(
            resolved.ms2.window_cap,
            Some(WindowCap {
                max_peaks: 1000,
                window_da: 100.0
            })
        );
        assert_eq!(resolved.ms2.max_peaks, default.ms2.max_peaks);
        assert_eq!(resolved.ms1.window_cap, default.ms1.window_cap);
        assert_eq!(resolved.ms1.im_pct_tol, default.ms1.im_pct_tol);
    }

    #[test]
    fn indexing_mz_tol_takes_either_unit() {
        let toml = format!(
            "{MINIMAL_TOML}\n[indexing.ms1]\nmz_tol = {{ Ppm = 5.0 }}\n[indexing.ms2]\nmz_tol = {{ Bins = 1 }}\n"
        );
        let resolved = toml::from_str::<Config>(&toml)
            .unwrap()
            .indexing
            .unwrap()
            .resolve();
        assert_eq!(resolved.ms1.mz_tol, MzTolerance::Ppm(5.0));
        assert_eq!(resolved.ms2.mz_tol, MzTolerance::Bins(1));
    }

    #[test]
    fn indexing_rejects_unknown_field() {
        let toml = format!("{MINIMAL_TOML}\n[indexing.ms2]\nmax_peak = 20000\n");
        assert!(toml::from_str::<Config>(&toml).is_err());
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

    /// A run that predicted its library leaves the provenance in
    /// `config_used.json`; a run that loaded one from a file predicted nothing,
    /// and an empty key would say it had.
    #[test]
    fn the_provenance_is_serialized_when_a_run_has_one_and_the_key_is_absent_otherwise() {
        let mut c: Config = toml::from_str(MINIMAL_TOML).unwrap();
        let without = serde_json::to_string(&c).unwrap();
        assert!(!without.contains("library_provenance"), "got: {without}");

        c.library_provenance = Some(serde_json::json!({
            "generator": { "tool": "msspeculator", "model": "builtin:small-v0" },
        }));
        let with = serde_json::to_string(&c).unwrap();
        assert!(with.contains("builtin:small-v0"), "got: {with}");

        let reread: Config =
            serde_json::from_str(&with).expect("config_used.json has to parse as a Config");
        assert_eq!(reread.library_provenance, c.library_provenance);
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
        assert!(
            c.library_provenance.is_none(),
            "provenance is recorded by a run, so the template must not mention it"
        );
        assert!(c.analysis.raw_inputs.is_none());

        assert_eq!(
            serde_json::to_string(&c.calibration).unwrap(),
            serde_json::to_string(&CalibrationConfig::default()).unwrap(),
            "[calibration] in default_config.toml drifted from CalibrationConfig::default()"
        );
    }
}
