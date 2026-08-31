//! Where a run's inputs come from.
//!
//! A run names its raw files, its library and its output directory either on
//! the command line or in a configuration file. The command line wins, and the
//! merge happens once, before anything reads a value, so there is one resolved
//! value per input rather than two sources consulted at different points.

use crate::cli::SearchArgs;
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
    args: &SearchArgs,
    mut config: Config,
) -> Result<(Config, ResolvedInputs), CliError> {
    merge_cli_into_config(&mut config, args);
    let resolved = read_inputs(args, &config)?;
    Ok((config, resolved))
}

/// A flag beats the configuration file whenever the flag was given, whatever
/// its value.
fn merge_cli_into_config(config: &mut Config, args: &SearchArgs) {
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

fn read_inputs(args: &SearchArgs, config: &Config) -> Result<ResolvedInputs, CliError> {
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
    use crate::cli::Cli;
    use clap::Parser;
    use timsseek::DecoyPolicy;

    /// The search arguments a command line resolves to, whether they were
    /// given bare or under `search`.
    fn search(argv: &[&str]) -> Cli {
        Cli::parse_from(std::iter::once("timsseek").chain(argv.iter().copied()))
    }

    fn resolve(argv: &[&str], config: Config) -> Result<ResolvedInputs, CliError> {
        let cli = search(argv);
        resolve_run_inputs(cli.search_args(), config).map(|(_, resolved)| resolved)
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

    /// The invocation the README documented before the rename, which spells all
    /// three inputs with their deprecated aliases. Existing scripts and
    /// container entrypoints still contain it.
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
            search(&[
                "--speclib-uri",
                "lib.ndjson",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--decoy-strategy",
                "never",
            ])
            .search_args(),
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

    /// Nothing is swallowed. Every argument is a flag, so no positional can
    /// absorb a stray token, and a near-miss spelling is an error rather than a
    /// prefix match. Adding a positional later would remove that guarantee
    /// without any other test noticing, which is why this one is here.
    #[test]
    fn nothing_is_swallowed() {
        use clap::error::ErrorKind;

        let unknown_flag = |argv: [&str; 3]| Cli::try_parse_from(argv).unwrap_err().kind();
        assert_eq!(
            unknown_flag(["timsseek", "--speclib-urk", "lib.ndjson"]),
            ErrorKind::UnknownArgument,
            "a misspelled flag"
        );
        assert_eq!(
            unknown_flag(["timsseek", "--speclib", "lib.ndjson"]),
            ErrorKind::UnknownArgument,
            "a prefix of a real flag is not that flag"
        );

        let stray_word = |argv: [&str; 2]| Cli::try_parse_from(argv).unwrap_err().kind();
        assert_eq!(
            stray_word(["timsseek", "serach"]),
            ErrorKind::InvalidSubcommand,
            "a typo of a real subcommand"
        );
        assert_eq!(
            stray_word(["timsseek", "lib.ndjson"]),
            ErrorKind::InvalidSubcommand,
            "a stray word that is not a subcommand at all"
        );
    }

    /// `search` names what a bare invocation already does, so the two have to
    /// resolve identically.
    #[test]
    fn the_search_subcommand_is_the_bare_invocation_named() {
        let bare = resolve(
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
        let named = resolve(
            &[
                "search",
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
        assert_eq!(bare, named);
    }

    /// Flags belong on one side of the subcommand or the other, never both.
    ///
    /// Two copies of `SearchArgs` exist, the flattened one and the subcommand's,
    /// and `search_args()` returns one of them. A flag before the subcommand
    /// lands in the copy it does not return, so anything short of a parse error
    /// here silently discards it.
    #[test]
    fn flags_cannot_be_split_across_a_subcommand() {
        for argv in [
            vec!["timsseek", "--speclib-uri", "lib.ndjson", "search"],
            vec!["timsseek", "--raw-inputs", "a.d", "search"],
            vec!["timsseek", "--overwrite", "build-library"],
        ] {
            let err = Cli::try_parse_from(&argv)
                .expect_err("a flag before a subcommand must not be silently dropped");
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::ArgumentConflict,
                "{argv:?}"
            );
        }
    }

    /// The deprecation notice is driven by the raw arguments, because clap
    /// reports only the variant it selected and not the word that chose it.
    #[test]
    fn the_deprecated_subcommand_is_still_accepted() {
        assert!(
            crate::cli::DEPRECATED_SUBCOMMANDS
                .iter()
                .any(|(spelling, _)| *spelling == "search"),
            "`search` is what the warning is keyed on"
        );
        assert!(
            Cli::try_parse_from(["timsseek", "search", "--speclib-uri", "lib.ndjson"]).is_ok(),
            "deprecated, not removed"
        );
    }
}
