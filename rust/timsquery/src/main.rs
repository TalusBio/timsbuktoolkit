use clap::{
    Parser,
    Subcommand,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::sync::Arc;
use std::time::Instant;

use timscentroid::{
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsquery::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use timsquery::models::elution_group::ElutionGroup;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
    Tolerance,
};
use tracing::instrument;
use tracing::subscriber::set_global_default;
use tracing_bunyan_formatter::{
    BunyanFormattingLayer,
    JsonStorageLayer,
};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

// Require for the trait.
use timsquery::QueriableData;

// mimalloc seems to work better for windows
// ... more accurately ... not using it causes everyting to
// be extremely slow on windows...
#[cfg(target_os = "windows")]
use mimalloc::MiMalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    let formatting_layer = BunyanFormattingLayer::new("timsquery".into(), std::io::stdout);
    let subscriber = Registry::default()
        .with(env_filter)
        .with(JsonStorageLayer)
        .with(formatting_layer);

    set_global_default(subscriber).expect("Setting default subscriber failed");
    let args = Args::parse();

    match args.command {
        Some(Commands::QueryIndex(args)) => main_query_index(args),
        Some(Commands::WriteTemplate(args)) => main_write_template(args),
        None => {
            println!("No command provided");
        }
    }
}

fn main_write_template(args: WriteTemplateArgs) {
    let output_path = args.output_path;
    let num_elution_groups = args.num_elution_groups;
    let egs = template_elution_groups(num_elution_groups);
    let tolerance = template_tolerance_settings();

    // Serialize both and write as files in the output path.
    // Do pretty serialization.
    let egs_json = serde_json::to_string_pretty(&egs).unwrap();
    let tolerance_json = serde_json::to_string_pretty(&tolerance).unwrap();
    let tolerance_json_narrow = serde_json::to_string_pretty(
        &tolerance.with_rt_tolerance(RtTolerance::Minutes((5.0, 5.0))),
    )
    .unwrap();

    let put_path = std::path::Path::new(&output_path);
    std::fs::create_dir_all(put_path).unwrap();
    println!("Writing to {}", put_path.display());
    let egs_json_path = put_path.join("elution_groups.json");
    let tolerance_json_path = put_path.join("tolerance_settings.json");
    let tolerance_json_narrow_path = put_path.join("tolerance_settings_narrow.json");

    println!("\n>>> Example elution groups: \n{}\n", &egs_json);
    std::fs::write(egs_json_path.clone(), egs_json).unwrap();
    println!("\n>>> Example tolerances: \n{}\n", &tolerance_json);
    std::fs::write(tolerance_json_path.clone(), tolerance_json).unwrap();
    println!(
        "\n>>> Example tolerances (narrow): \n{}\n",
        &tolerance_json_narrow
    );
    std::fs::write(tolerance_json_narrow_path.clone(), tolerance_json_narrow).unwrap();
    println!(
        "use as `timsquery query-index --output-path '.' --raw-file-path 'your_file.d' --tolerance-settings-path {:#?} --elution-groups-path {:#?}`",
        tolerance_json_path, egs_json_path,
    );
}

fn template_elution_groups(num: usize) -> Vec<ElutionGroup<usize>> {
    let mut egs = Vec::with_capacity(num);

    let min_mz = 400.0;
    let max_mz = 900.0;
    let max_mobility = 0.8;
    let min_mobility = 0.6;
    let min_rt = 66.0;
    let max_rt = 11.0 * 60.0;

    let rt_step = (max_rt - min_rt) / (num as f32);
    let mobility_step = (max_mobility - min_mobility) / (num as f32);
    let mz_step = (max_mz - min_mz) / (num as f64);

    for i in 1..num {
        let rt = min_rt + (i as f32 * rt_step);
        let mobility = min_mobility + (i as f32 * mobility_step);
        let mz = min_mz + (i as f64 * mz_step);
        let fragment_mzs = (0..10).map(|x| (x as usize, mz + x as f64));
        let fragments = Arc::from(Vec::from_iter(fragment_mzs));
        egs.push(ElutionGroup {
            id: i as u64,
            rt_seconds: rt,
            mobility,
            precursors: Arc::from(vec![(0, mz), (1, mz + 1.0)]),
            fragments,
        });
    }
    egs
}

fn get_ms1_rts_as_millis(file: &TimsTofPath) -> Arc<[u32]> {
    let reader = file.load_frame_reader().unwrap();
    let mut rts: Vec<_> = reader
        .frame_metas
        .iter()
        .map(|f| (f.rt_in_seconds * 1000.0).round() as u32)
        .collect();
    rts.sort_unstable();
    rts.dedup();
    rts.into()
}

#[instrument]
fn main_query_index(args: QueryIndexArgs) {
    let raw_file_path = args.raw_file_path;
    let tolerance_settings_path = args.tolerance_settings_path;
    let elution_groups_path = args.elution_groups_path;
    let aggregator_use = args.aggregator;

    let tolerance_settings: Tolerance =
        serde_json::from_str(&std::fs::read_to_string(&tolerance_settings_path).unwrap()).unwrap();
    let elution_groups: Vec<ElutionGroup<String>> =
        serde_json::from_str(&std::fs::read_to_string(&elution_groups_path).unwrap()).unwrap();

    let file = TimsTofPath::new(&raw_file_path).unwrap();
    let centroiding_config = timscentroid::CentroidingConfig {
        max_peaks: 50_000,
        mz_ppm_tol: 5.0,
        im_pct_tol: 3.0,
        early_stop_iterations: 200,
    };
    let (index, building_stats) = IndexedTimstofPeaks::from_timstof_file(&file, centroiding_config);
    let rts = get_ms1_rts_as_millis(&file);
    println!("Indexing Stats: {:#?}", building_stats);
    let mut queries = AggregatorContainer::new(elution_groups, aggregator_use, rts);

    let output_path = args.output_path;
    std::fs::create_dir_all(&output_path).unwrap();
    let serialization_format = args.format;
    queries.add_query(&index, &tolerance_settings);
    queries.serialize_write(serialization_format, &output_path);
}

fn template_tolerance_settings() -> Tolerance {
    Tolerance {
        ms: MzTolerance::Ppm((15.0, 15.0)),
        // rt: RtTolerance::Absolute((120.0, 120.0)),
        rt: RtTolerance::Unrestricted,
        mobility: MobilityTolerance::Pct((10.0, 10.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum PossibleAggregator {
    #[default]
    PointIntensityAggregator,
    ChromatogramAggregator,
    SpectrumAggregator,
}

pub enum AggregatorContainer {
    Point(Vec<PointIntensityAggregator<String>>),
    Chromatogram(Vec<ChromatogramCollector<String, f32>>),
    Spectrum(Vec<SpectralCollector<String, f32>>),
}

impl AggregatorContainer {
    fn new(
        queries: Vec<ElutionGroup<String>>,
        aggregator: PossibleAggregator,
        ref_rts: Arc<[u32]>,
    ) -> Self {
        match aggregator {
            PossibleAggregator::PointIntensityAggregator => AggregatorContainer::Point(
                queries
                    .into_iter()
                    .map(|x| PointIntensityAggregator::new_with_elution_group(x.into()))
                    .collect(),
            ),
            PossibleAggregator::ChromatogramAggregator => AggregatorContainer::Chromatogram(
                queries
                    .into_iter()
                    .map(|x| ChromatogramCollector::new(x, ref_rts.clone()).unwrap())
                    .collect(),
            ),
            PossibleAggregator::SpectrumAggregator => AggregatorContainer::Spectrum(
                queries
                    .into_iter()
                    .map(|x| SpectralCollector::new(x))
                    .collect(),
            ),
        }
    }

    fn add_query(&mut self, index: &IndexedTimstofPeaks, tolerance: &Tolerance) {
        match self {
            AggregatorContainer::Point(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            }
            AggregatorContainer::Spectrum(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            }
        }
    }

    fn serialize_inner(&self, format: SerializationFormat) -> String {
        // There has to be a better way of doing this ...
        macro_rules! serialize_format {
            ($format:ident, $aggregators:ident) => {
                match format {
                    SerializationFormat::PrettyJson => {
                        println!("Pretty printing enabled");
                        serde_json::to_string_pretty(&$aggregators).unwrap()
                    }
                    SerializationFormat::Json => serde_json::to_string(&$aggregators).unwrap(),
                }
            };
        }
        match self {
            AggregatorContainer::Point(aggregators) => {
                serialize_format!(format, aggregators)
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                serialize_format!(format, aggregators)
            }
            AggregatorContainer::Spectrum(aggregators) => {
                serialize_format!(format, aggregators)
            }
        }
    }

    fn serialize_write(&self, serialization_format: SerializationFormat, output_path: &str) {
        let serialization_start = Instant::now();
        let put_path = std::path::Path::new(&output_path).join("results.json");
        let serialized = self.serialize_inner(serialization_format);
        std::fs::write(put_path.clone(), serialized).unwrap();
        println!("Wrote to {}", put_path.display());
        let serialization_elapsed = serialization_start.elapsed();
        println!("Serialization took {:#?}", serialization_elapsed);
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum PossibleIndex {
    #[default]
    ExpandedRawFrameIndex,
    TransposedQuadIndex,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum SerializationFormat {
    #[default]
    PrettyJson,
    Json,
}

#[derive(Parser, Debug, Clone)]
pub struct QueryIndexArgs {
    /// The path to the raw file to query.
    #[arg(short, long)]
    raw_file_path: String,

    /// The path to the json file with the tolerance settings.
    #[arg(short, long)]
    tolerance_settings_path: String,

    /// The path to the json file with the elution groups.
    #[arg(short, long)]
    elution_groups_path: String,

    /// The path to the output files.
    #[arg(short, long)]
    output_path: String,

    /// The format to use for the output
    #[arg(short, long, default_value_t, value_enum)]
    format: SerializationFormat,

    // The aggregator to use.
    #[arg(short, long, default_value_t, value_enum)]
    aggregator: PossibleAggregator,
}

#[derive(Parser, Debug)]
struct WriteTemplateArgs {
    /// The path to the output files.
    #[arg(short, long)]
    output_path: String,

    /// The number of elution groups to generate.
    #[arg(short, long, default_value_t = 10)]
    num_elution_groups: usize,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Query the index.
    QueryIndex(QueryIndexArgs),
    WriteTemplate(WriteTemplateArgs),
}

#[derive(Debug, Serialize, Deserialize)]
struct ElutionGroupResults<T: Serialize> {
    elution_group: ElutionGroup<String>,
    result: T,
}
