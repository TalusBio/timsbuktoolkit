use clap::{
    Parser,
    Subcommand,
};
use serde::{
    Deserialize,
    Serialize,
};
use timsquery::{QueriableData, GenerallyQueriable};
use std::collections::HashMap;
use std::time::Instant;
use timsquery::models::aggregators::{
    EGSAggregator,
    PointIntensityAggregator,
    EGCAggregator,
};
use timsquery::models::elution_group::ElutionGroup;
use timsquery::models::indices::{
    ExpandedRawFrameIndex,
    QuadSplittedTransposedIndex,
};
use timsquery::models::tolerance::{
    Tolerance,
    MobilityTolerance,
    MzToleramce,
    QuadTolerance,
    RtTolerance,
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
use std::sync::Arc;

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

    let put_path = std::path::Path::new(&output_path);
    std::fs::create_dir_all(put_path).unwrap();
    println!("Writing to {}", put_path.display());
    let egs_json_path = put_path.join("elution_groups.json");
    let tolerance_json_path = put_path.join("tolerance_settings.json");
    std::fs::write(egs_json_path.clone(), egs_json).unwrap();
    std::fs::write(tolerance_json_path.clone(), tolerance_json).unwrap();
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

#[instrument]
fn main_query_index(args: QueryIndexArgs) {
    let args_clone = args.clone();

    let raw_file_path = args.raw_file_path;
    let tolerance_settings_path = args.tolerance_settings_path;
    let elution_groups_path = args.elution_groups_path;
    let index_use = args.index;
    let aggregator_use = args.aggregator;

    let tolerance_settings: Tolerance =
        serde_json::from_str(&std::fs::read_to_string(&tolerance_settings_path).unwrap()).unwrap();
    let elution_groups: Vec<Arc<ElutionGroup<String>>> =
        serde_json::from_str(&std::fs::read_to_string(&elution_groups_path).unwrap()).unwrap();

    let (index, rts) = index_use.build_index(&raw_file_path);
    let mut queries = AggregatorContainer::new(elution_groups.clone(), aggregator_use, rts);

    execute_query(
        index_use,
        aggregator_use,
        raw_file_path,
        tolerance_settings,
        elution_groups,
        args_clone,
    );
}

fn template_tolerance_settings() -> Tolerance {
    Tolerance {
        ms: MzToleramce::Ppm((15.0, 15.0)),
        // rt: RtTolerance::Absolute((120.0, 120.0)),
        rt: RtTolerance::None,
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
    Chromatogram(Vec<EGSAggregator<String>>),
    Spectrum(Vec<EGSAggregator<String>>),
}

impl AggregatorContainer {
    fn new(queries: Vec<Arc<ElutionGroup<String>>>, aggregator: PossibleAggregator, ref_rts: Arc<u32>) -> Self {
        match aggregator {
            PossibleAggregator::PointIntensityAggregator => {
                AggregatorContainer::Point(queries.iter().map(|x| PointIntensityAggregator::new_with_elution_group(x.clone())).collect())
            }
            PossibleAggregator::ChromatogramAggregator => {
                AggregatorContainer::Chromatogram(queries.iter().map(|x| EGCAggregator::new(x.clone(), ref_rts.clone())).collect())
            }
            PossibleAggregator::SpectrumAggregator => {
                AggregatorContainer::Spectrum(queries.iter().map(|x| EGSAggregator::new(x.clone())).collect())
            }
        }
    }

    fn add_query(&mut self, index: &dyn GenerallyQueriable<String>, tolerance: &Tolerance) {
        match self {
            AggregatorContainer::Point(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            },
            AggregatorContainer::Chromatogram(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            },
            AggregatorContainer::Spectrum(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            },
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum PossibleIndex {
    #[default]
    ExpandedRawFrameIndex,
    TransposedQuadIndex,
}

impl PossibleIndex {
    pub fn build_index(&self, raw_file_path: &str) -> (Box<dyn GenerallyQueriable<String>>, Arc<[u32]>) {
        match self {
            PossibleIndex::ExpandedRawFrameIndex => {
                let tmp = ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
                let rts = tmp.cycle_rt_ms.clone();
                (Box::new(tmp), rts)
            }
            PossibleIndex::TransposedQuadIndex => {
                let tmp =  QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap();
                let rts = tmp.cycle_rt_ms.clone();
                (Box::new(tmp), rts)
            },
        }
    }
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

    // The index to use.
    #[arg(short, long, value_enum)]
    index: PossibleIndex,
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

pub fn execute_query(
    index: PossibleIndex,
    aggregator: PossibleAggregator,
    raw_file_path: String,
    tolerance: Tolerance,
    elution_groups: Vec<Arc<ElutionGroup<String>>>,
    args: QueryIndexArgs,
) {
    let output_path = args.output_path;
    let serialization_format = args.format;

    macro_rules! execute_query_inner {
        ($index:expr, $agg:expr) => {
            let mut tmp: Vec<_> = elution_groups.iter().map(|x| $agg(x.clone())).collect();
            let tmp = &$index.par_add_query_multi(&tolerance, &elution_groups, &$agg);

            std::fs::create_dir_all(output_path.clone()).unwrap();

            let serialization_start = Instant::now();
            // TODO make this a macro ...
            let put_path = match serialization_format {
                SerializationFormat::PrettyJson => {
                    println!("Pretty printing enabled");
                    let put_path = std::path::Path::new(&output_path).join("results.json");
                    let serialized = serde_json::to_string_pretty(&out).unwrap();
                    std::fs::write(put_path.clone(), serialized).unwrap();
                    put_path
                }
                SerializationFormat::Json => {
                    let serialized = serde_json::to_string(&out).unwrap();
                    let put_path = std::path::Path::new(&output_path).join("results.json");
                    std::fs::write(put_path.clone(), serialized).unwrap();
                    put_path
                }
            };
            println!("Wrote to {}", put_path.display());
            let serialization_elapsed = serialization_start.elapsed();
            println!("Serialization took {:#?}", serialization_elapsed);
        };
    }

    match (index, aggregator) {
        (PossibleIndex::TransposedQuadIndex, aggregator) => {
            let index = QuadSplittedTransposedIndex::from_path_centroided(&(raw_file_path.clone()))
                .unwrap();
            match aggregator {
                PossibleAggregator::RawPeakIntensityAggregator => {
                    let aggregator = PointIntensityAggregator::new_with_elution_group;
                    execute_query_inner!(index, aggregator);
                }
                PossibleAggregator::RawPeakVectorAggregator => {
                    let aggregator = RawPeakVectorAggregator::new_with_elution_group;
                    execute_query_inner!(index, aggregator);
                }
                PossibleAggregator::MultiCMGStats => {
                    let factory = MultiCMGStatsFactory {
                        converters: (index.mz_converter, index.im_converter),
                        _phantom: std::marker::PhantomData::<String>,
                    };
                    execute_query_inner!(index, |x| factory.build_with_elution_group(x, None));
                }
            }
        }
        (PossibleIndex::ExpandedRawFrameIndex, aggregator) => {
            let index = ExpandedRawFrameIndex::from_path(&(raw_file_path.clone())).unwrap();
            match aggregator {
                PossibleAggregator::RawPeakIntensityAggregator => {
                    let aggregator = PointIntensityAggregator::new_with_elution_group;
                    execute_query_inner!(index, aggregator);
                }
                PossibleAggregator::RawPeakVectorAggregator => {
                    let aggregator = RawPeakVectorAggregator::new_with_elution_group;
                    execute_query_inner!(index, aggregator);
                }
                PossibleAggregator::MultiCMGStats => {
                    let factory = MultiCMGStatsFactory {
                        converters: (index.mz_converter, index.im_converter),
                        _phantom: std::marker::PhantomData::<String>,
                    };
                    execute_query_inner!(index, |x| factory.build_with_elution_group(x, None));
                }
            }
        }

        (PossibleIndex::Unspecified, _) => {
            panic!("Should have been handled");
        }
    }
}
