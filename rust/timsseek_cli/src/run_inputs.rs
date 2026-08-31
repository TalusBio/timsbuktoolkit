//! Where a run's inputs come from.
//!
//! A run names its raw files, its library and its output directory either on
//! the command line or in a configuration file. The command line wins, and the
//! merge happens once, before anything reads a value, so there is one resolved
//! value per input rather than two sources consulted at different points.

use crate::cli::Cli;
use crate::config::{
    Config,
    InputConfig,
    OutputConfig,
};
use crate::errors::CliError;
use tims_stage::expand_local_uri;

/// A run's inputs after the merge, with `~` expanded so validation, staging and
/// existence probes all see the same spelling.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedInputs {
    pub raw_inputs: Vec<String>,
    pub speclib_uri: String,
    pub calib_lib_uri: Option<String>,
    pub output_uri: String,
    pub overwrite: bool,
}

/// Merge the command line into `config`, then read the run's inputs back out of
/// it. Returns the merged configuration alongside the inputs, because the rest
/// of the run reads tolerances, chunk size and the decoy strategy from it and
/// writes it to `config_used.json`.
///
/// `--calib-lib` and `--overwrite` have no configuration counterpart and are
/// read straight from the command line. See the spec's Out of Scope for why
/// calibration input is not a configuration field yet.
pub fn resolve_run_inputs(
    args: &Cli,
    mut config: Config,
) -> Result<(Config, ResolvedInputs), CliError> {
    merge_cli_into_config(&mut config, args);
    let resolved = read_inputs(args, &config)?;
    Ok((config, resolved))
}

/// A flag beats the configuration file whenever the flag was given, whatever
/// its value.
fn merge_cli_into_config(config: &mut Config, args: &Cli) {
    if let Some(raw_inputs) = &args.raw_inputs {
        config.analysis.raw_inputs = Some(raw_inputs.clone());
    }
    if let Some(uri) = &args.speclib_uri {
        config.input = Some(InputConfig::Speclib { uri: uri.clone() });
    }
    if let Some(uri) = &args.output_uri {
        config.output = Some(OutputConfig { uri: uri.clone() });
    }
    if let Some(strategy) = args.decoy_strategy {
        config.analysis.decoy_strategy = strategy;
    }
    if let Some(model) = args.rescore_model {
        config.analysis.rescore_model = model.into();
    }
}

fn read_inputs(args: &Cli, config: &Config) -> Result<ResolvedInputs, CliError> {
    let raw_inputs = match &config.analysis.raw_inputs {
        Some(files) => files.iter().map(|f| expand_local_uri(f)).collect(),
        None => {
            return Err(CliError::Config {
                source: "No raw files provided, please provide raw_inputs in either the config file or with the --raw-inputs flag".to_string(),
            });
        }
    };

    let speclib_uri = match &config.input {
        Some(InputConfig::Speclib { uri }) => expand_local_uri(uri),
        None => {
            return Err(CliError::Config {
                source: "No input provided, please provide one in either the config file or with the --speclib-uri flag".to_string(),
            });
        }
    };

    let output_uri = match &config.output {
        Some(output) => expand_local_uri(&output.uri),
        None => {
            return Err(CliError::Config {
                source: "No output directory specified".to_string(),
            });
        }
    };

    Ok(ResolvedInputs {
        raw_inputs,
        speclib_uri,
        calib_lib_uri: args.calib_lib.as_deref().map(expand_local_uri),
        output_uri,
        overwrite: args.overwrite,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use timsseek::DecoyPolicy;

    fn cli(argv: &[&str]) -> Cli {
        Cli::parse_from(std::iter::once("timsseek").chain(argv.iter().copied()))
    }

    fn resolve(argv: &[&str], config: Config) -> Result<ResolvedInputs, CliError> {
        resolve_run_inputs(&cli(argv), config).map(|(_, resolved)| resolved)
    }

    fn config_with(extra: &str) -> Config {
        let base = r#"
[analysis]
chunk_size = 20000

[analysis.tolerance]
ms = { Ppm = [15.0, 15.0] }
mobility = { Pct = [5.0, 5.0] }
quad = { Absolute = [0.1, 0.1] }
rt = "Unrestricted"
"#;
        toml::from_str(&format!("{extra}{base}")).expect("test configuration must parse")
    }

    /// The invocation in the project's README, which spells all three inputs
    /// with their pre-rename aliases.
    #[test]
    fn readme_invocation_resolves_through_the_deprecated_aliases() {
        let resolved = resolve(
            &[
                "--speclib-file",
                "vimentin.ndjson",
                "--output-dir",
                "vimentin_search_results",
                "--dotd-files",
                "my_data.d",
            ],
            Config::default_config(),
        )
        .unwrap();

        assert_eq!(resolved.speclib_uri, "vimentin.ndjson");
        assert_eq!(resolved.output_uri, "vimentin_search_results");
        assert_eq!(resolved.raw_inputs, vec!["my_data.d".to_string()]);
    }

    #[test]
    fn canonical_spellings_resolve_to_the_same_inputs_as_the_aliases() {
        let aliased = resolve(
            &[
                "--speclib-file",
                "lib.ndjson",
                "--output-dir",
                "out",
                "--dotd-files",
                "a.d",
            ],
            Config::default_config(),
        )
        .unwrap();
        let canonical = resolve(
            &[
                "--speclib-uri",
                "lib.ndjson",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
            ],
            Config::default_config(),
        )
        .unwrap();
        assert_eq!(aliased, canonical);
    }

    #[test]
    fn a_flag_beats_the_configuration_file() {
        let config = config_with(
            r#"
[input]
type = "speclib"
uri = "from_config.ndjson"

[output]
uri = "config_results"
"#,
        );
        let resolved = resolve(
            &[
                "--speclib-uri",
                "from_flag.ndjson",
                "--output-uri",
                "flag_results",
                "--raw-inputs",
                "from_flag.d",
            ],
            config,
        )
        .unwrap();

        assert_eq!(resolved.speclib_uri, "from_flag.ndjson");
        assert_eq!(resolved.output_uri, "flag_results");
        assert_eq!(resolved.raw_inputs, vec!["from_flag.d".to_string()]);
    }

    #[test]
    fn the_configuration_file_is_used_where_no_flag_was_given() {
        let config = config_with(
            r#"
[input]
type = "speclib"
uri = "from_config.ndjson"

[output]
uri = "config_results"
"#,
        );
        let mut config = config;
        config.analysis.raw_inputs = Some(vec!["from_config.d".to_string()]);

        let resolved = resolve(&[], config).unwrap();

        assert_eq!(resolved.speclib_uri, "from_config.ndjson");
        assert_eq!(resolved.output_uri, "config_results");
        assert_eq!(resolved.raw_inputs, vec!["from_config.d".to_string()]);
    }

    /// `--decoy-strategy` overrides the same way the input flags do, and the
    /// merged configuration is what the rest of the run reads.
    #[test]
    fn a_flag_beats_the_configuration_file_for_non_input_settings() {
        let (merged, _) = resolve_run_inputs(
            &cli(&[
                "--speclib-uri",
                "lib.ndjson",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--decoy-strategy",
                "never",
            ]),
            config_with(""),
        )
        .unwrap();
        assert_eq!(merged.analysis.decoy_strategy, DecoyPolicy::Never);
    }

    #[test]
    fn a_run_with_no_library_names_the_flag_that_supplies_one() {
        let err = resolve(
            &["--raw-inputs", "a.d", "--output-uri", "out"],
            config_with(""),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("--speclib-uri"), "got: {err}");
    }

    #[test]
    fn calib_lib_is_read_from_the_command_line_and_never_from_the_configuration() {
        let resolved = resolve(
            &[
                "--speclib-uri",
                "lib.ndjson",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--calib-lib",
                "calib.ndjson",
            ],
            config_with(""),
        )
        .unwrap();
        assert_eq!(resolved.calib_lib_uri.as_deref(), Some("calib.ndjson"));
    }

    /// Every argument is a flag, so nothing can absorb a stray token. Adding a
    /// positional later would silently remove that guarantee.
    #[test]
    fn an_unrecognized_flag_is_an_error() {
        let err = Cli::try_parse_from(["timsseek", "--speclib-urk", "lib.ndjson"]).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
    }

    #[test]
    fn an_unrecognized_first_word_is_an_error() {
        let err = Cli::try_parse_from(["timsseek", "serach"]).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
    }
}
