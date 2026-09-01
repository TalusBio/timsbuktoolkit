//! Where a run's inputs come from.
//!
//! A run names its raw files, its library, its sequence database and its output
//! directory either on the command line or in a configuration file. The command
//! line wins, and the merge happens once, before anything reads a value, so
//! there is one resolved value per input rather than two sources consulted at
//! different points.

use std::path::{
    Path,
    PathBuf,
};

use crate::artifacts::built_library_path;
use crate::cli::SearchArgs;
use crate::config::{
    Config,
    InputConfig,
    OutputConfig,
    SequencesConfig,
};
use crate::errors::CliError;
use tims_stage::{
    expand_local_uri,
    is_remote_uri,
};

/// A run's inputs after the merge, with `~` expanded so validation, staging and
/// existence probes all see the same spelling.
///
/// Holds what the pipeline opens. An input a run records without opening is read
/// from the merged [`Config`], which is what `config_used.json` is written from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedInputs {
    pub raw_inputs: Vec<String>,
    pub library: LibrarySource,
    pub calib_lib_uri: Option<String>,
    pub output_uri: String,
    pub overwrite: bool,
}

/// Where the library a run scores against comes from.
///
/// Not an `Option<String>` with the FASTA read from elsewhere, because the two
/// are one decision: a run has exactly one library and either opens it or
/// predicts it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LibrarySource {
    /// A library file, named by `--speclib-uri` or `[input]`. May be remote.
    File(String),
    /// No library was named, so one is predicted from this sequence database and
    /// never written down.
    Fasta(PathBuf),
    /// The library to search is not on disk yet, so `--build-if-missing`
    /// predicts one from this sequence database, writes it to `out`, and opens
    /// what it wrote. Local only: writing is done through the filesystem.
    Build { out: PathBuf, fasta: PathBuf },
}

/// Whether the library a run named is already there.
///
/// Injected so the routing table can be exercised without touching a filesystem
/// or making a network round trip; the run itself passes
/// [`probe_uri_exists`](crate::output_sink::probe_uri_exists), which is the same
/// probe the collision check uses and answers for remote URIs too.
type ExistsFn<'a> = &'a dyn Fn(&str) -> Result<bool, CliError>;

/// Which of the three supplies the library.
///
/// A named library is used whether or not a FASTA was also given. Naming a
/// library is what asks for it, and the FASTA is an input in its own right --
/// the sequences a run reports against are not a substitute for its library --
/// so one arriving does not reinterpret the other.
///
/// `build_if_missing` is what turns a FASTA into a library on disk. Without it a
/// library that is not there is an error, because inferring "predict it" from a
/// missing file would turn a mistyped path into minutes of prediction; with it,
/// the library named is the library written, and a run that finds one already
/// there just opens it.
fn library_source(
    input: Option<&InputConfig>,
    sequences: Option<&SequencesConfig>,
    build_if_missing: bool,
    output_uri: &str,
    exists: ExistsFn<'_>,
) -> Result<LibrarySource, CliError> {
    let fasta = sequences.map(|SequencesConfig { fasta }| expand_local_path(fasta));
    match (input, fasta) {
        (Some(InputConfig::Speclib { uri, .. }), fasta) => {
            let uri = expand_local_uri(uri);
            // The probe is skipped without the flag: the library is opened
            // either way, and `validate_inputs` is where a missing one is
            // reported.
            if !build_if_missing || exists(&uri)? {
                return Ok(LibrarySource::File(uri));
            }
            let Some(fasta) = fasta else {
                return Err(CliError::Config {
                    source: format!(
                        "{uri} does not exist, and --build-if-missing has nothing to predict it \
                         from; name a sequence database with --fasta, in the config file or on \
                         the command line"
                    ),
                });
            };
            reject_remote_build(&uri)?;
            Ok(LibrarySource::Build {
                out: PathBuf::from(uri),
                fasta,
            })
        }
        // A derived library is an output of the run rather than an input it was
        // pointed at, so it is not probed here: the collision check reports one
        // an earlier run left behind, alongside every other artifact in the way.
        (None, Some(fasta)) if build_if_missing => {
            reject_remote_build(output_uri)?;
            Ok(LibrarySource::Build {
                out: built_library_path(output_uri, &fasta),
                fasta,
            })
        }
        (None, Some(fasta)) => Ok(LibrarySource::Fasta(fasta)),
        (None, None) => Err(CliError::Config {
            source: "No library and no sequences provided; either name a library with \
                     --speclib-uri or name a sequence database to predict one from with --fasta, \
                     in the config file or on the command line"
                .to_string(),
        }),
    }
}

/// A library is written through the filesystem, so a remote destination is
/// named rather than left to fail at open.
///
/// The same limitation `build-library` states for `--out`, reached the same way:
/// a search stages the remote inputs it *reads*, and staging a write back is a
/// separate piece of machinery.
fn reject_remote_build(destination: &str) -> Result<(), CliError> {
    if !is_remote_uri(destination) {
        return Ok(());
    }
    Err(CliError::Config {
        source: format!(
            "{destination} is a remote URI, and --build-if-missing writes a library through the \
             filesystem. Build it locally with `timsseek build-library`, upload it, and name it \
             with --speclib-uri."
        ),
    })
}

/// `~` expanded on a path, so a FASTA is spelled the way every URI is.
fn expand_local_path(path: &Path) -> PathBuf {
    PathBuf::from(expand_local_uri(&path.to_string_lossy()))
}

/// Merge the command line into `config`, then read the run's inputs back out of
/// it. Returns the merged configuration alongside the inputs, because the rest
/// of the run reads tolerances, chunk size and the decoy strategy from it and
/// writes it to `config_used.json`.
///
/// Every input has a configuration home, so `config_used.json` records the whole
/// run. `--overwrite` and `--build-if-missing` are the exceptions and are read
/// straight from argv: each is a decision about one invocation rather than about
/// how to search. `--overwrite` is recorded in `run_report.json` instead, and
/// what `--build-if-missing` builds is recorded by the library it writes.
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
    // The two libraries merge independently: naming one on the command line
    // must not silently drop the other's configured value.
    if let Some(uri) = &args.speclib_uri {
        match &mut config.input {
            Some(InputConfig::Speclib {
                uri: configured, ..
            }) => *configured = uri.clone(),
            None => {
                config.input = Some(InputConfig::Speclib {
                    uri: uri.clone(),
                    calib_uri: None,
                })
            }
        }
    }
    if let Some(calib) = &args.calib_lib {
        match &mut config.input {
            Some(InputConfig::Speclib { calib_uri, .. }) => *calib_uri = Some(calib.clone()),
            // A calibration library with no library to calibrate is rejected by
            // `read_inputs`, which reports the missing one by flag name.
            None => {}
        }
    }
    // The sequence database is merged independently of the library, because
    // `[sequences]` and `[input]` describe two inputs that are given together as
    // readily as they are given apart.
    if let Some(fasta) = &args.fasta {
        config.sequences = Some(SequencesConfig {
            fasta: fasta.clone(),
        });
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

    // Ahead of the library, which is named after the output directory when
    // `--build-if-missing` was given nowhere else to write it.
    let output_uri = match &config.output {
        Some(output) => expand_local_uri(&output.uri),
        None => {
            return Err(CliError::Config {
                source: "No output directory specified".to_string(),
            });
        }
    };

    let library = library_source(
        config.input.as_ref(),
        config.sequences.as_ref(),
        args.build_if_missing,
        &output_uri,
        &crate::output_sink::probe_uri_exists,
    )?;
    let calib_lib_uri = match &config.input {
        Some(InputConfig::Speclib { calib_uri, .. }) => calib_uri.as_deref().map(expand_local_uri),
        None => None,
    };

    Ok(ResolvedInputs {
        raw_inputs,
        library,
        calib_lib_uri,
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

    /// The library a run resolved to, for the runs that name a file.
    fn file(uri: &str) -> LibrarySource {
        LibrarySource::File(uri.to_string())
    }

    /// One row of the routing table: what was named, whether the named library
    /// is already there, and whether the flag that builds one was given.
    ///
    /// The output directory is a literal, because every row that reads it reads
    /// it the same way.
    fn route(
        uri: Option<&str>,
        fasta: Option<&str>,
        build_if_missing: bool,
        exists: bool,
    ) -> Result<LibrarySource, CliError> {
        route_into("out", uri, fasta, build_if_missing, exists)
    }

    /// A row of the table for a run whose output directory is worth naming.
    fn route_into(
        output_uri: &str,
        uri: Option<&str>,
        fasta: Option<&str>,
        build_if_missing: bool,
        exists: bool,
    ) -> Result<LibrarySource, CliError> {
        let input = uri.map(|uri| InputConfig::Speclib {
            uri: uri.to_string(),
            calib_uri: None,
        });
        let sequences = fasta.map(|fasta| SequencesConfig {
            fasta: PathBuf::from(fasta),
        });
        library_source(
            input.as_ref(),
            sequences.as_ref(),
            build_if_missing,
            output_uri,
            &move |_| Ok(exists),
        )
    }

    /// The library a run builds, for the rows that predict one and write it.
    fn build(out: &str, fasta: &str) -> LibrarySource {
        LibrarySource::Build {
            out: PathBuf::from(out),
            fasta: PathBuf::from(fasta),
        }
    }

    /// The search arguments a command line resolves to, whether they were
    /// given bare or under `search`.
    fn search(argv: &[&str]) -> Cli {
        Cli::parse_from(std::iter::once("timsseek").chain(argv.iter().copied()))
    }

    fn resolve(argv: &[&str], config: Config) -> Result<ResolvedInputs, CliError> {
        resolve_run_inputs(search(argv).search_args(), config).map(|(_, resolved)| resolved)
    }

    /// Both halves of the merge: the configuration as it would be written to
    /// `config_used.json`, and the inputs read back out of it.
    fn resolve_pair(argv: &[&str], config: Config) -> (Config, ResolvedInputs) {
        resolve_run_inputs(search(argv).search_args(), config).expect("fixture resolves")
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

    /// The deprecated aliases still resolve, and to exactly what the canonical
    /// spellings do.
    ///
    /// The values are the invocation the README documented before the rename,
    /// because that is the one still sitting in scripts and container
    /// entrypoints. Asserting the two spellings against each other proves the
    /// aliases are wired up; asserting one of them against literals proves both
    /// are wired to the right place, which agreement alone would not -- two
    /// spellings that resolved to the same wrong value would pass.
    #[test]
    fn the_deprecated_aliases_resolve_exactly_as_the_canonical_spellings_do() {
        let aliased = resolve(
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
        .expect("the README invocation resolves");
        let canonical = resolve(
            &[
                "--speclib-uri",
                "vimentin.ndjson",
                "--output-uri",
                "vimentin_search_results",
                "--raw-inputs",
                "my_data.d",
            ],
            Config::default_config(),
        )
        .expect("the canonical invocation resolves");

        assert_eq!(aliased, canonical);
        assert_eq!(aliased.library, file("vimentin.ndjson"));
        assert_eq!(aliased.output_uri, "vimentin_search_results");
        assert_eq!(aliased.raw_inputs, vec!["my_data.d".to_string()]);
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

        assert_eq!(resolved.library, file("from_flag.ndjson"));
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

        assert_eq!(resolved.library, file("from_config.ndjson"));
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

    /// Both flags, because a library file is no longer the only way to supply
    /// one: a FASTA to predict from is the other.
    #[test]
    fn a_run_with_no_library_names_both_flags_that_supply_one() {
        let err = resolve(
            &["--raw-inputs", "a.d", "--output-uri", "out"],
            config_with(""),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("--speclib-uri"), "got: {err}");
        assert!(err.contains("--fasta"), "got: {err}");
    }

    /// A named library is what a run searches, whether or not a FASTA came with
    /// it, and a lone FASTA is predicted from.
    #[test]
    fn a_named_library_wins_and_a_lone_fasta_is_predicted_from() {
        assert_eq!(
            route(Some("lib.mzspeclib.txt"), None, false, true).expect("a library alone resolves"),
            file("lib.mzspeclib.txt")
        );
        assert_eq!(
            route(
                Some("lib.mzspeclib.txt"),
                Some("proteome.fasta"),
                false,
                true
            )
            .expect("both resolve"),
            file("lib.mzspeclib.txt"),
            "a named library is used, so the FASTA predicts nothing"
        );
        assert_eq!(
            route(None, Some("proteome.fasta"), false, false).expect("a FASTA alone resolves"),
            LibrarySource::Fasta(PathBuf::from("proteome.fasta")),
            "a lone FASTA is predicted into memory and never written down"
        );
    }

    /// Naming neither is the one row no flag rescues, so the message has to name
    /// both of the things that would.
    #[test]
    fn a_run_that_names_no_library_and_no_sequences_says_what_to_name() {
        for build_if_missing in [false, true] {
            let err = route(None, None, build_if_missing, false)
                .expect_err("a run with no library and no sequences has nothing to score against")
                .to_string();
            assert!(err.contains("--speclib-uri"), "got: {err}");
            assert!(err.contains("--fasta"), "got: {err}");
        }
    }

    /// The point of the flag: a library that is already there is opened, so a
    /// script can carry the flag every run and pay for the prediction once.
    #[test]
    fn the_build_flag_is_inert_when_the_library_is_already_there() {
        for fasta in [None, Some("proteome.fasta")] {
            let without = route(Some("lib.mzspeclib.txt.gz"), fasta, false, true)
                .expect("an existing library resolves");
            let with = route(Some("lib.mzspeclib.txt.gz"), fasta, true, true)
                .expect("an existing library resolves with the flag too");
            assert_eq!(with, without, "fasta: {fasta:?}");
            assert_eq!(with, file("lib.mzspeclib.txt.gz"));
        }
    }

    /// The row the flag exists for: the library named is the library written,
    /// and the search reads back what it just predicted.
    #[test]
    fn a_missing_library_with_the_flag_is_built_where_it_was_named() {
        assert_eq!(
            route(
                Some("lib.mzspeclib.txt.gz"),
                Some("proteome.fasta"),
                true,
                false
            )
            .expect("a missing library with a FASTA and the flag is built"),
            build("lib.mzspeclib.txt.gz", "proteome.fasta")
        );
    }

    /// Without the flag a missing library stays a missing library, even with a
    /// FASTA sitting right there: `validate_inputs` is what reports it, and a
    /// mistyped path must not buy minutes of prediction.
    #[test]
    fn a_missing_library_without_the_flag_is_still_opened_rather_than_built() {
        assert_eq!(
            route(
                Some("lib.mzspeclib.txt.gz"),
                Some("proteome.fasta"),
                false,
                false
            )
            .expect("the missing file is reported at validation, not here"),
            file("lib.mzspeclib.txt.gz")
        );
    }

    /// The flag asks for a library to be built and names nothing to build it
    /// from, so the error is about the sequence database rather than about the
    /// file that is missing.
    #[test]
    fn a_missing_library_with_the_flag_and_no_fasta_names_the_flag_that_supplies_one() {
        let err = route(Some("lib.mzspeclib.txt.gz"), None, true, false)
            .expect_err("there is nothing to predict from")
            .to_string();
        assert!(err.contains("--fasta"), "got: {err}");
        assert!(err.contains("lib.mzspeclib.txt.gz"), "got: {err}");
    }

    /// With nowhere named to write it, the library lands in the output
    /// directory under the FASTA's own name -- which is what makes the second
    /// run of the same command find one.
    #[test]
    fn a_fasta_and_the_flag_alone_write_the_library_into_the_output_directory() {
        assert_eq!(
            route_into("results", None, Some("/seq/proteome.fasta"), true, false)
                .expect("a FASTA and the flag resolve"),
            build(
                &PathBuf::from("results")
                    .join("proteome.mzspeclib.txt.gz")
                    .to_string_lossy(),
                "/seq/proteome.fasta"
            )
        );
    }

    /// A build writes through the filesystem, so a remote destination is named
    /// up front rather than surfacing as a failure to open a bucket as a file.
    /// Both destinations a build can be pointed at get the same treatment.
    #[test]
    fn a_library_that_would_have_to_be_built_remotely_is_rejected_by_name() {
        let named = route(
            Some("s3://bkt/lib.mzspeclib.txt.gz"),
            Some("proteome.fasta"),
            true,
            false,
        )
        .expect_err("a remote library cannot be written")
        .to_string();
        assert!(named.contains("remote URI"), "got: {named}");
        assert!(named.contains("build-library"), "got: {named}");

        let derived = route_into("s3://bkt/out", None, Some("proteome.fasta"), true, false)
            .expect_err("a remote output directory cannot hold a library a search built")
            .to_string();
        assert!(derived.contains("remote URI"), "got: {derived}");

        assert_eq!(
            route(Some("s3://bkt/lib.mzspeclib.txt.gz"), None, true, true)
                .expect("a remote library that is already there is just opened"),
            file("s3://bkt/lib.mzspeclib.txt.gz"),
            "the flag stays inert against a library that exists, remote or not"
        );
    }

    /// A search from a FASTA and nothing else, which is the invocation that
    /// writes no library to disk.
    #[test]
    fn a_run_that_names_only_a_fasta_predicts_its_own_library() {
        let (_, resolved) = resolve_pair(
            &[
                "--fasta",
                "proteome.fasta",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
            ],
            config_with(""),
        );
        assert_eq!(
            resolved.library,
            LibrarySource::Fasta(PathBuf::from("proteome.fasta"))
        );
    }

    /// The flag is read off argv rather than out of the merged configuration,
    /// so the wiring from the command line to the routing is worth one run
    /// through `resolve_run_inputs` -- with a library path that is not there,
    /// which is the whole premise of the flag.
    #[test]
    fn the_build_flag_reaches_the_routing_from_the_command_line() {
        let (_, resolved) = resolve_pair(
            &[
                "--speclib-uri",
                "absent.mzspeclib.txt.gz",
                "--fasta",
                "proteome.fasta",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--build-if-missing",
            ],
            config_with(""),
        );
        assert_eq!(
            resolved.library,
            build("absent.mzspeclib.txt.gz", "proteome.fasta")
        );
    }

    #[test]
    fn the_calibration_library_comes_from_either_source() {
        let from_flag = resolve(
            &[
                "--speclib-uri",
                "lib.mzspeclib.txt",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--calib-lib",
                "calib.mzspeclib.txt",
            ],
            config_with(""),
        )
        .unwrap();
        assert_eq!(
            from_flag.calib_lib_uri.as_deref(),
            Some("calib.mzspeclib.txt")
        );

        let from_config = resolve(
            &["--output-uri", "out", "--raw-inputs", "a.d"],
            config_with(
                "[input]\ntype = \"speclib\"\nuri = \"lib.mzspeclib.txt\"\n\
                 calib_uri = \"calib.mzspeclib.txt\"\n",
            ),
        )
        .unwrap();
        assert_eq!(
            from_config.calib_lib_uri.as_deref(),
            Some("calib.mzspeclib.txt")
        );
    }

    /// `--speclib-uri` and `--calib-lib` reach the same config field, so naming
    /// one on the command line must not drop the other's configured value.
    #[test]
    fn naming_one_library_keeps_the_other() {
        let resolved = resolve(
            &[
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--speclib-uri",
                "flag.mzspeclib.txt",
            ],
            config_with(
                "[input]\ntype = \"speclib\"\nuri = \"config.mzspeclib.txt\"\n\
                 calib_uri = \"calib.mzspeclib.txt\"\n",
            ),
        )
        .unwrap();
        assert_eq!(resolved.library, file("flag.mzspeclib.txt"));
        assert_eq!(
            resolved.calib_lib_uri.as_deref(),
            Some("calib.mzspeclib.txt"),
            "the flag replaced the library, not the calibration library"
        );
    }

    /// A sequence database and a spectral library are independent inputs: one
    /// can be predicted from the other, and neither replaces it, so naming both
    /// is an ordinary invocation rather than a conflict -- and the named library
    /// is the one searched.
    #[test]
    fn a_fasta_and_a_library_are_accepted_together_and_the_library_is_searched() {
        let (merged, resolved) = resolve_pair(
            &[
                "--speclib-uri",
                "lib.mzspeclib.txt",
                "--fasta",
                "proteome.fasta",
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
            ],
            config_with(""),
        );
        assert_eq!(
            resolved.library,
            file("lib.mzspeclib.txt"),
            "the FASTA does not turn a named library into a prediction"
        );
        assert_eq!(
            merged.sequences.map(|sequences| sequences.fasta),
            Some(PathBuf::from("proteome.fasta"))
        );
    }

    /// `[sequences]` is a home for the FASTA that does not go through
    /// `[library]`, so a run that predicts nothing can still name one.
    #[test]
    fn a_configured_sequences_section_supplies_the_fasta() {
        let (merged, _) = resolve_pair(
            &["--output-uri", "out", "--raw-inputs", "a.d"],
            config_with(
                "[input]\ntype = \"speclib\"\nuri = \"lib.mzspeclib.txt\"\n\
                 [sequences]\nfasta = \"config.fasta\"\n",
            ),
        );
        assert_eq!(
            merged.sequences.map(|sequences| sequences.fasta),
            Some(PathBuf::from("config.fasta"))
        );
    }

    #[test]
    fn the_fasta_flag_beats_the_configured_sequences_section() {
        let (merged, _) = resolve_pair(
            &[
                "--output-uri",
                "out",
                "--raw-inputs",
                "a.d",
                "--speclib-uri",
                "lib.mzspeclib.txt",
                "--fasta",
                "flag.fasta",
            ],
            config_with("[sequences]\nfasta = \"config.fasta\"\n"),
        );
        assert_eq!(
            merged.sequences.map(|sequences| sequences.fasta),
            Some(PathBuf::from("flag.fasta"))
        );
    }

    /// `config_used.json` is the artifact a run is reproduced from, and it is
    /// serialized `Config`. So every input has to survive a round trip through
    /// it -- which is what a command-line-only `--calib-lib` could not do.
    #[test]
    fn every_input_survives_the_config_used_round_trip() {
        let argv = [
            "--speclib-uri",
            "lib.mzspeclib.txt",
            "--output-uri",
            "out",
            "--raw-inputs",
            "a.d",
            "--calib-lib",
            "calib.mzspeclib.txt",
            "--fasta",
            "proteome.fasta",
        ];
        let (config, first) = resolve_pair(&argv, config_with(""));

        let written = serde_json::to_string(&config).expect("config serializes");
        let reread: Config = serde_json::from_str(&written).expect("config_used.json is a Config");

        // Re-resolved with no flags at all: whatever the artifact holds is all
        // there is to go on the second time.
        let (rebuilt, second) = resolve_pair(&["--raw-inputs", "a.d"], reread);
        assert_eq!(first, second);
        assert_eq!(rebuilt.sequences, config.sequences);
        assert!(
            rebuilt.sequences.is_some(),
            "the FASTA reached the artifact"
        );
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
