use clap::{
    Parser,
    ValueEnum,
};
use std::path::PathBuf;
use timsseek::ml::RescoreModel;
use timsseek::{
    DecoyPolicy,
    UnannotatedPeaks,
};

/// Clap mirror of [`timsseek::ml::RescoreModel`]: `ValueEnum` is a foreign
/// trait, so it cannot be implemented for the lib-owned type from this crate.
///
/// `rename_all` matches `RescoreModel`'s serde spelling so model names are
/// identical on the command line and in configuration files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[clap(rename_all = "snake_case")]
pub enum CliRescoreModel {
    Gbm,
    Lda,
    Hybrid,
    Mlp,
}

impl From<CliRescoreModel> for RescoreModel {
    fn from(v: CliRescoreModel) -> Self {
        match v {
            CliRescoreModel::Gbm => RescoreModel::Gbm,
            CliRescoreModel::Lda => RescoreModel::Lda,
            CliRescoreModel::Hybrid => RescoreModel::Hybrid,
            CliRescoreModel::Mlp => RescoreModel::Mlp,
        }
    }
}

#[cfg(test)]
impl From<RescoreModel> for CliRescoreModel {
    fn from(v: RescoreModel) -> Self {
        match v {
            RescoreModel::Gbm => CliRescoreModel::Gbm,
            RescoreModel::Lda => CliRescoreModel::Lda,
            RescoreModel::Hybrid => CliRescoreModel::Hybrid,
            RescoreModel::Mlp => CliRescoreModel::Mlp,
        }
    }
}

/// `timsseek` with no subcommand is a search, so the search arguments are
/// flattened at the top level and the subcommand is optional. Every existing
/// invocation keeps parsing to the same values.
///
/// `args_conflicts_with_subcommands` is what makes the flattened arguments safe
/// alongside a subcommand: two copies of [`SearchArgs`] exist and only one is
/// read, so mixing the two has to be a parse error rather than a silent choice
/// between them.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(args_conflicts_with_subcommands = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    #[command(flatten)]
    pub search: SearchArgs,
}

#[derive(clap::Subcommand, Debug)]
pub enum Command {
    /// Deprecated spelling of a bare invocation; use `timsseek <SEARCH ARGS>`.
    /// Removed after 2026-12-31.
    Search(Box<SearchArgs>),
    /// Predict a spectral library from a FASTA and write it.
    BuildLibrary(Box<BuildLibraryArgs>),
}

/// Subcommand spelling put on a clock, warned about by
/// [`warn_deprecated_spellings`], paired with the instruction that replaces it.
pub const DEPRECATED_SUBCOMMANDS: &[(&str, &str)] =
    &[("search", "run timsseek with no subcommand")];

/// Flag spellings put on a clock, paired with the flag that replaces each.
///
/// Keyed with the leading `--` because these are matched against raw arguments,
/// where `dotd-files` on its own is a value rather than a flag.
pub const DEPRECATED_FLAGS: &[(&str, &str)] = &[
    ("--dotd-files", "--raw-inputs"),
    ("--speclib-file", "--speclib-uri"),
    ("--output-dir", "--output-uri"),
];

/// The date the deprecated spellings stop being accepted, in the first release
/// after it. The removal is scheduled work rather than a constant nobody reads:
/// see <https://github.com/TalusBio/timsbuktoolkit/issues/116>.
pub const DEPRECATION_REMOVAL_DATE: &str = "2026-12-31";

impl Cli {
    /// The search arguments this invocation resolved to, whether they were
    /// given bare or under the deprecated `search`.
    ///
    /// Reading one and discarding the other is only safe because
    /// `args_conflicts_with_subcommands` makes both-populated unrepresentable.
    pub fn search_args(&self) -> &SearchArgs {
        match &self.command {
            Some(Command::Search(args)) => args,
            _ => &self.search,
        }
    }
}

/// Warn once per deprecated spelling actually typed.
///
/// Called before the subscriber is installed, so this writes to stderr rather
/// than going through `tracing`.
pub fn warn_deprecated_spellings<I, S>(argv: I)
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    for warning in deprecation_warnings(argv) {
        eprintln!("{warning}");
    }
}

/// One warning per distinct deprecated spelling in `argv`, in the order typed.
///
/// clap resolves an alias to its canonical name and does not report which
/// spelling produced it, so the raw arguments are what has to be scanned. That
/// covers `--flag value` and `--flag=value` alike, and a subcommand name, which
/// clap likewise reports only as the variant it selected.
///
/// Separate from the printing so the wording is reachable from a test without
/// capturing stderr.
pub fn deprecation_warnings<I, S>(argv: I) -> Vec<String>
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let mut seen: Vec<&str> = Vec::new();
    let mut warnings = Vec::new();
    for arg in argv {
        // `--flag=value` and `--flag value` both have to match, so compare the
        // part left of any `=`.
        let arg = arg.as_ref();
        let Some((spelling, advice)) = replacement_advice(arg.split('=').next().unwrap_or(arg))
        else {
            continue;
        };
        if seen.contains(&spelling) {
            continue;
        }
        seen.push(spelling);
        warnings.push(format!(
            "warning: `{spelling}` is deprecated and will be removed after \
             {DEPRECATION_REMOVAL_DATE}; {advice}"
        ));
    }
    warnings
}

/// What to tell someone who typed `spelling`, or `None` if it is not deprecated.
///
/// The two tables are looked up alike but do not read alike: a subcommand is
/// replaced by a way of invoking timsseek, a flag by another flag.
fn replacement_advice(spelling: &str) -> Option<(&'static str, String)> {
    let found = |table: &'static [(&'static str, &'static str)]| {
        table
            .iter()
            .find(|(deprecated, _)| *deprecated == spelling)
            .copied()
    };
    if let Some((spelling, instruction)) = found(DEPRECATED_SUBCOMMANDS) {
        return Some((spelling, format!("{instruction} instead")));
    }
    let (spelling, flag) = found(DEPRECATED_FLAGS)?;
    Some((spelling, format!("use `{flag}` instead")))
}

#[derive(clap::Args, Debug)]
pub struct SearchArgs {
    /// Path to the log file.
    /// Defaults to {output_dir}/timsseek.log.
    /// Use "-" to send logs to stderr instead of a file.
    #[arg(long, value_name = "PATH")]
    pub log_path: Option<PathBuf>,

    /// Log level for the log file (default: info)
    #[arg(long, value_name = "LEVEL", default_value = "info")]
    pub log_level: String,

    /// Path to the configuration file, TOML or JSON (sniffed by extension). Optional; defaults are used if omitted. Use `--write-default-config` to scaffold one.
    #[arg(short, long)]
    pub config: Option<PathBuf>,

    /// Path to the .d file(s) (will over-write the config file)
    /// Can be specified multiple times for batch processing
    #[arg(long = "raw-inputs", alias = "dotd-files", value_name = "URI")]
    pub raw_inputs: Option<Vec<String>>,

    /// Path to the speclib file (will over-write the config file)
    #[arg(long = "speclib-uri", alias = "speclib-file")]
    pub speclib_uri: Option<String>,

    /// URI of a calibration library (optional).
    /// If provided, Phase 1 prescore uses this library for calibrant selection,
    /// while Phase 3 scoring uses the main speclib.
    /// If not provided, the main speclib is used for both phases.
    /// Accepts local paths or s3:// URIs.
    #[arg(long, value_name = "URI")]
    pub calib_lib: Option<String>,

    /// Sequence database for this run, overriding `[sequences]` in the
    /// configuration file. Local paths only.
    ///
    /// Independent of `--speclib-uri`: without a library this predicts one;
    /// with a library it compares the recorded build provenance against the
    /// requested prediction settings and warns on differences.
    #[arg(long, value_name = "PATH")]
    pub fasta: Option<PathBuf>,

    /// Predict the library named by `--speclib-uri` when it is not there yet,
    /// write it, and search it. With no library named, it is written to
    /// `<output>/<FASTA stem>.mzspeclib.txt.gz`. Requires `--fasta`.
    ///
    /// A library that is already there is opened and the flag does nothing, so
    /// a script can carry it every run and pay for the prediction once.
    ///
    /// No configuration counterpart, for the reason `--overwrite` has none:
    /// whether to spend minutes predicting a library is a decision about this
    /// invocation rather than about how to search.
    #[arg(long)]
    pub build_if_missing: bool,

    /// Path to the output directory
    #[arg(long = "output-uri", short = 'o', alias = "output-dir")]
    pub output_uri: Option<String>,

    /// Overwrite existing output directory if it exists
    #[arg(short = 'O', long)]
    pub overwrite: bool,

    /// Maximum q-value for output. Only results at or below this
    /// threshold are written to the Parquet file.
    /// Default: 0.5
    #[arg(long, default_value = "0.5")]
    pub max_qvalue: f32,

    /// Decoy generation strategy
    /// Options: if-missing (default), force, never
    /// - if-missing: Generate mass-shift decoys only if library has none
    /// - force: Drop library decoys and regenerate mass-shift decoys
    /// - never: Use library as-is without decoy generation
    #[arg(long, value_name = "STRATEGY")]
    pub decoy_strategy: Option<DecoyPolicy>,

    /// What to do with a library peak that carries no annotation.
    /// Options: keep (default), skip, fail, keep-all
    /// - keep: store it at the m/z the file measured, under a `?N` label
    /// - skip: drop it, leaving the entry only its annotated peaks
    /// - fail: refuse the library at the first one
    /// - keep-all: store every peak that way, annotated or not
    ///
    /// Read only for a library that annotates nothing at all, keep-all aside:
    /// a peak with no annotation beside peaks that have one is mzPAF's `?`,
    /// which is the writer calling it noise.
    #[arg(long, value_name = "POLICY")]
    pub unannotated_peaks: Option<UnannotatedPeaks>,

    /// Rescore model; overrides the config file.
    #[arg(long, value_enum)]
    pub rescore_model: Option<CliRescoreModel>,

    /// Print the default TOML configuration to stdout and exit.
    #[arg(long, conflicts_with = "write_default_config")]
    pub print_default_config: bool,

    /// Write the default TOML configuration to the given path and exit.
    #[arg(long, value_name = "PATH")]
    pub write_default_config: Option<PathBuf>,

    /// Skip writing the post-rescore sidecars (results.feature_stats.tsv and
    /// results.feature_importance.tsv).
    #[arg(long)]
    pub no_feature_stats: bool,
}

/// Prediction settings, mirroring msspeculator's own library command rather
/// than inventing a second vocabulary for the same settings.
///
/// Every field is `Option`, because "the flag was given" is what makes a flag
/// beat the configuration file. The fallbacks live in
/// [`LibraryConfig`](crate::config::LibraryConfig), which in turn falls back to
/// msspeculator's defaults.
#[derive(clap::Args, Debug)]
pub struct BuildLibraryArgs {
    /// Sequence database to digest.
    #[arg(long, value_name = "PATH")]
    pub fasta: PathBuf,

    /// Where to write the library. The suffix picks the format:
    /// `.mzspeclib.txt` writes mzSpecLib, anything else writes DIA-NN TSV, and
    /// a trailing `.gz` compresses either. `.mzspeclib.txt.gz` is the
    /// recommended spelling.
    #[arg(long, short = 'o', value_name = "PATH")]
    pub out: PathBuf,

    /// Path to the configuration file, TOML or JSON. Its `[library]` section
    /// supplies anything not given as a flag.
    #[arg(long, short = 'c', value_name = "PATH")]
    pub config: Option<PathBuf>,

    /// A `.safetensors` artifact, or `builtin:NAME` for one compiled in.
    #[arg(long, value_name = "MODEL")]
    pub model: Option<String>,

    #[arg(long, value_name = "N")]
    pub missed_cleavages: Option<usize>,
    #[arg(long, value_name = "N")]
    pub min_length: Option<usize>,
    #[arg(long, value_name = "N")]
    pub max_length: Option<usize>,
    #[arg(long, value_name = "Z")]
    pub min_charge: Option<i64>,
    #[arg(long, value_name = "Z")]
    pub max_charge: Option<i64>,

    /// Fixed modification rule, repeatable.
    #[arg(long = "fixed-mod", value_name = "TARGETS[MOD]")]
    fixed_mods: Option<Vec<String>>,
    /// Predict with no fixed modifications at all.
    ///
    /// Needed because the default is carbamidomethyl and an empty
    /// `--fixed-mod` list cannot be spelled: repeating a flag zero times is
    /// indistinguishable from not passing it.
    #[arg(long, conflicts_with = "fixed_mods")]
    pub no_fixed_mods: bool,
    /// Variable modification rule, repeatable.
    #[arg(long = "variable-mod", value_name = "TARGETS[MOD]")]
    pub variable_mods: Option<Vec<String>>,
    #[arg(long, value_name = "N")]
    pub max_variable_mods: Option<usize>,

    // No acquisition or chromatography context. The model artifact has its own
    // default for both, and selecting a different one is a decision about the
    // model rather than about the library this command writes -- so it belongs
    // wherever the model is chosen, not on a second set of flags here.
    /// Drop fragments below this base-peak-relative intensity.
    #[arg(long, value_name = "FRACTION")]
    pub min_intensity: Option<f64>,
    /// Keep at most this many of the strongest fragments per precursor, applied
    /// after `--min-intensity`.
    #[arg(long, value_name = "N")]
    pub max_fragments: Option<usize>,

    /// Add pseudo-reversed decoy precursors, on by default. `--no-decoys` turns
    /// them off, leaving a search to derive weaker mass-shift decoys instead.
    ///
    /// Two flags on one field, so "not given" stays distinct from "off": a bare
    /// `bool` could not disable a configured `decoys = true`, which left the
    /// configuration file with the last word on a flag.
    #[arg(long, overrides_with = "no_decoys")]
    decoys: bool,
    #[arg(long, overrides_with = "decoys")]
    no_decoys: bool,

    /// Replace the library if it already exists.
    ///
    /// Predicting a proteome takes minutes, so overwriting one by accident is
    /// expensive in both directions: the old library is gone and the new one
    /// costs the same wait again.
    #[arg(short = 'O', long)]
    pub overwrite: bool,

    /// Where to write the resolved-configuration sidecar. Defaults to
    /// `<out>.config.json`.
    #[arg(long, value_name = "PATH", conflicts_with = "no_config_out")]
    pub config_out: Option<PathBuf>,
    /// Skip the resolved-configuration sidecar.
    #[arg(long)]
    pub no_config_out: bool,
}

impl BuildLibraryArgs {
    /// `Some(true)`/`Some(false)` when either decoy flag was given, `None` when
    /// neither was -- which is what lets the configuration file supply the value
    /// only in the absence of a flag.
    ///
    /// `overrides_with` makes the later flag win rather than erroring, so
    /// `--decoys --no-decoys` is `Some(false)` and both being set at once is
    /// unreachable.
    pub fn decoys(&self) -> Option<bool> {
        match (self.decoys, self.no_decoys) {
            (false, false) => None,
            (decoys, _) => Some(decoys),
        }
    }

    /// The fixed-modification rules asked for, where an empty list is a real
    /// answer and `None` means "not given".
    pub fn fixed_mods(&self) -> Option<Vec<String>> {
        if self.no_fixed_mods {
            return Some(Vec::new());
        }
        self.fixed_mods.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every deprecated flag names its replacement and the day it stops
    /// working, in whichever of the two spellings a script happens to use.
    #[test]
    fn a_deprecated_flag_warns_and_names_its_replacement() {
        for (deprecated, replacement) in DEPRECATED_FLAGS {
            for argv in [
                vec![deprecated.to_string(), "value".to_string()],
                vec![format!("{deprecated}=value")],
            ] {
                let warnings = deprecation_warnings(&argv);
                assert_eq!(warnings.len(), 1, "{argv:?}");
                let warning = &warnings[0];
                assert!(warning.contains(deprecated), "{warning}");
                assert!(warning.contains(replacement), "{warning}");
                assert!(warning.contains(DEPRECATION_REMOVAL_DATE), "{warning}");
            }
        }
    }

    /// Nobody using the current spellings should be told to change anything.
    #[test]
    fn the_canonical_spellings_warn_about_nothing() {
        let argv = [
            "--raw-inputs",
            "my_data.d",
            "--speclib-uri=lib.ndjson",
            "--output-uri",
            "out",
        ];
        assert!(deprecation_warnings(argv).is_empty());
    }

    /// `--dotd-files` is repeatable, and a warning per `.d` file would bury the
    /// one line that matters.
    #[test]
    fn a_repeated_deprecated_flag_warns_once() {
        let argv = ["--dotd-files", "a.d", "--dotd-files", "b.d"];
        assert_eq!(deprecation_warnings(argv).len(), 1);
    }

    /// A value that happens to spell a deprecated flag without the dashes is a
    /// path, not a flag.
    #[test]
    fn a_positional_value_spelled_like_a_flag_warns_about_nothing() {
        assert!(deprecation_warnings(["--raw-inputs", "dotd-files"]).is_empty());
    }

    /// The subcommand's replacement is a way of invoking timsseek, so its
    /// warning must not read as if `search` were a flag with a flag to swap in.
    #[test]
    fn the_deprecated_subcommand_warns_as_a_subcommand() {
        let warnings = deprecation_warnings(["search", "--speclib-uri", "lib.ndjson"]);
        assert_eq!(warnings.len(), 1);
        let warning = &warnings[0];
        assert!(
            warning.contains("run timsseek with no subcommand"),
            "{warning}"
        );
        assert!(warning.contains(DEPRECATION_REMOVAL_DATE), "{warning}");
        assert!(!warning.contains("use `--"), "{warning}");
    }

    /// Deprecated, not removed: each alias still has to land in the field its
    /// replacement fills, which is what a later cleanup of the aliases would
    /// break.
    #[test]
    fn the_deprecated_flags_still_parse_to_their_canonical_fields() {
        let parse = |argv: [&str; 3]| {
            Cli::try_parse_from(argv)
                .expect("a deprecated spelling is still accepted")
                .search
        };
        assert_eq!(
            parse(["timsseek", "--dotd-files", "my_data.d"]).raw_inputs,
            Some(vec!["my_data.d".to_string()])
        );
        assert_eq!(
            parse(["timsseek", "--speclib-file", "lib.ndjson"]).speclib_uri,
            Some("lib.ndjson".to_string())
        );
        assert_eq!(
            parse(["timsseek", "--output-dir", "out"]).output_uri,
            Some("out".to_string())
        );
    }

    /// The tables are what the warnings are keyed on, so a flag listed there
    /// without a matching `alias` in clap would warn about a spelling that does
    /// not parse.
    #[test]
    fn every_deprecated_flag_is_a_clap_alias() {
        for (deprecated, _) in DEPRECATED_FLAGS {
            assert!(
                Cli::try_parse_from(["timsseek", deprecated, "value"]).is_ok(),
                "{deprecated} is warned about but not accepted"
            );
        }
    }

    /// The direction that bites: an alias added to clap and not to the table is
    /// accepted forever without ever saying so. clap knows its own aliases, so
    /// the table can be checked against them rather than against a memory.
    #[test]
    fn every_clap_alias_is_a_deprecated_flag() {
        use clap::CommandFactory;
        let command = Cli::command();
        let arguments = command.get_arguments().chain(
            command
                .get_subcommands()
                .flat_map(clap::Command::get_arguments),
        );
        for argument in arguments {
            for alias in argument.get_all_aliases().unwrap_or_default() {
                let spelling = format!("--{alias}");
                assert!(
                    DEPRECATED_FLAGS
                        .iter()
                        .any(|(deprecated, _)| *deprecated == spelling),
                    "{spelling} is an accepted alias of `--{}` but nothing warns about it",
                    argument.get_long().unwrap_or_default()
                );
            }
        }
    }
}
