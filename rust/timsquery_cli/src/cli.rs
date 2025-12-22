use clap::{
    Parser,
    Subcommand,
};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Query the index.
    QueryIndex(QueryIndexArgs),
    /// Write template configuration files.
    WriteTemplate(WriteTemplateArgs),
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum PossibleAggregator {
    PointIntensityAggregator,
    #[default]
    ChromatogramAggregator,
    SpectrumAggregator,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum SerializationFormat {
    Json,
    #[default]
    PrettyJson,
    Ndjson,
}

#[derive(Parser, Debug, Clone)]
pub struct QueryIndexArgs {
    /// The path to the raw file to query.
    #[arg(short, long)]
    pub raw_file_path: PathBuf,

    /// The path to the json file with the tolerance settings.
    #[arg(short, long)]
    pub tolerance_settings_path: PathBuf,

    /// The path to the json file with the elution groups.
    #[arg(short, long)]
    pub elution_groups_path: PathBuf,

    /// The path to the output files.
    #[arg(short, long)]
    pub output_path: PathBuf,

    /// The format to use for the output
    #[arg(short, long, default_value_t, value_enum)]
    pub format: SerializationFormat,

    /// The aggregator to use.
    #[arg(short, long, default_value_t, value_enum)]
    pub aggregator: PossibleAggregator,

    /// Batch size for streaming serialization (default: 500)
    #[arg(short, long, default_value_t = 500)]
    pub batch_size: usize,
}

#[derive(Parser, Debug)]
pub struct WriteTemplateArgs {
    /// The path to the output files.
    #[arg(short, long)]
    pub output_path: PathBuf,
}
