use clap::{
    Parser,
    Subcommand,
};
use rayon::iter::{
    ParallelDrainRange,
    ParallelIterator,
};
use serde::ser::{
    SerializeSeq,
    Serializer,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};
use timsquery::serde::load_index_caching;
use timsrust::MSLevel;
use tracing_subscriber::fmt::format::FmtSpan;

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
use tracing_subscriber::EnvFilter;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

// Require for the trait.
use timsquery::QueriableData;
use tracing::{
    error,
    info,
    warn,
};

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
    let subscriber = Registry::default()
        .with(env_filter)
        .with(tracing_subscriber::fmt::layer().with_span_events(FmtSpan::CLOSE));

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

#[instrument]
fn main_query_index(args: QueryIndexArgs) {
    let raw_file_path = args.raw_file_path;
    let tolerance_settings_path = args.tolerance_settings_path;
    let elution_groups_path = args.elution_groups_path;
    let aggregator_use = args.aggregator;

    let tolerance_settings: Tolerance =
        serde_json::from_str(&std::fs::read_to_string(&tolerance_settings_path).unwrap()).unwrap();
    info!("Using tolerance settings: {:#?}", tolerance_settings);
    info!("Loading elution groups from {}", elution_groups_path);
    let elution_groups: Vec<ElutionGroup<usize>> = read_query_elution_groups(&elution_groups_path);
    info!("Loaded {} elution groups", elution_groups.len());

    let index = load_index_caching(&raw_file_path).unwrap();
    let rts = get_ms1_rts_as_millis(&raw_file_path);

    let output_path = args.output_path;
    let serialization_format = args.format;
    let batch_size = args.batch_size;

    std::fs::create_dir_all(&output_path).unwrap();
    let put_path = std::path::Path::new(&output_path).join("results.json");

    stream_process_batches(
        elution_groups,
        aggregator_use,
        rts,
        &index,
        &tolerance_settings,
        serialization_format,
        &put_path,
        batch_size,
    );
}

fn read_query_elution_groups(path: &impl AsRef<std::path::Path>) -> Vec<ElutionGroup<usize>> {
    // We attempt to read as from elution group intpus, if we fail we try to use the direct elution group format.
    let tmp_inner = std::fs::read_to_string(&path.as_ref()).unwrap();
    let try_input: Result<Vec<ElutionGroupInput>, _> = serde_json::from_str(&tmp_inner);
    if let Ok(eg_inputs) = try_input {
        let out: Vec<ElutionGroup<usize>> = eg_inputs.into_iter().map(|x| x.into()).collect();
        return out;
    }
    // Fallback to direct elution group format. and if that does not work try to convert to Vec<ElutionGroup<String>>
    // warn that the ids will be dropped and converted to usize indices.
    let out = serde_json::from_str(&tmp_inner);
    let out = match out {
        Ok(egs) => egs,
        Err(e) => {
            warn!(
                "Failed to read elution groups as ElutionGroupInput or ElutionGroup<usize>: {}",
                e
            );
            warn!("Attempting to read as ElutionGroup<String> and convert ids to usize indices.");
            warn!("This will drop the original ids and replace them with usize indices.");
            let egs_string: Vec<ElutionGroup<String>> = serde_json::from_str(&tmp_inner).unwrap();
            let mut out: Vec<ElutionGroup<usize>> = Vec::with_capacity(egs_string.len());
            for (i, eg) in egs_string.into_iter().enumerate() {
                let eg_usize = ElutionGroup {
                    id: i as u64,
                    mobility: eg.mobility,
                    rt_seconds: eg.rt_seconds,
                    precursors: eg.precursors,
                    fragments: Arc::from(
                        eg.fragments
                            .iter()
                            .enumerate()
                            .map(|(i, (_lab, mz))| (i, *mz))
                            .collect::<Vec<(usize, f64)>>(),
                    ),
                };
                out.push(eg_usize);
            }
            out
        }
    };

    out
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

fn template_elution_groups(num: usize) -> Vec<ElutionGroupInput> {
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
        let fragments = (0..10).map(|x| mz + x as f64).collect();
        let precursors = (0..2).map(|x| mz + x as f64).collect();
        egs.push(ElutionGroupInput {
            id: i as u64,
            rt_seconds: rt,
            mobility,
            fragments,
            precursors,
        });
    }

    egs
}

fn get_ms1_rts_as_millis(file: &impl AsRef<Path>) -> Arc<[u32]> {
    let ttp = TimsTofPath::new(file.as_ref()).unwrap();
    let reader = ttp.load_frame_reader().unwrap();
    let mut rts: Vec<_> = reader
        .frame_metas
        .iter()
        .filter(|x| x.ms_level == MSLevel::MS1)
        .map(|f| (f.rt_in_seconds * 1000.0).round() as u32)
        .collect();
    rts.sort_unstable();
    rts.dedup();
    rts.into()
}

#[instrument(skip_all)]
fn stream_process_batches(
    elution_groups: Vec<ElutionGroup<usize>>,
    aggregator_use: PossibleAggregator,
    rts: Arc<[u32]>,
    index: &IndexedTimstofPeaks,
    tolerance: &Tolerance,
    serialization_format: SerializationFormat,
    output_path: &Path,
    batch_size: usize,
) {
    let total_groups = elution_groups.len();
    let total_batches = (total_groups + batch_size - 1) / batch_size;

    info!(
        "Processing {} elution groups in {} batches of up to {}",
        total_groups, total_batches, batch_size
    );

    let serialization_start = Instant::now();
    let file = File::create(output_path).unwrap();
    let writer = BufWriter::new(file);

    match serialization_format {
        SerializationFormat::PrettyJson => {
            let mut ser = serde_json::Serializer::pretty(writer);
            process_and_serialize(
                elution_groups,
                aggregator_use,
                rts,
                index,
                tolerance,
                &mut ser,
                total_groups,
                total_batches,
                batch_size,
                serialization_start,
                output_path,
            );
        }
        SerializationFormat::Json => {
            let mut ser = serde_json::Serializer::new(writer);
            process_and_serialize(
                elution_groups,
                aggregator_use,
                rts,
                index,
                tolerance,
                &mut ser,
                total_groups,
                total_batches,
                batch_size,
                serialization_start,
                output_path,
            );
        }
    }
}

#[instrument(skip_all)]
fn process_and_serialize<S>(
    elution_groups: Vec<ElutionGroup<usize>>,
    aggregator_use: PossibleAggregator,
    rts: Arc<[u32]>,
    index: &IndexedTimstofPeaks,
    tolerance: &Tolerance,
    ser: S,
    total_groups: usize,
    total_batches: usize,
    batch_size: usize,
    serialization_start: Instant,
    output_path: &Path,
) where
    S: Serializer,
{
    let mut seq = ser.serialize_seq(Some(total_groups)).unwrap();

    // Report max once per 2s progress
    let mut last_progress = Instant::now();
    let progress_interval = Duration::from_secs(2);

    for (batch_idx, chunk) in elution_groups.chunks(batch_size).enumerate() {
        let mut batch_start = None;
        if last_progress.elapsed() >= progress_interval {
            info!(
                "Processing batch {}/{} ({} groups)",
                batch_idx + 1,
                total_batches,
                chunk.len(),
            );
            batch_start = Some(Instant::now());
            last_progress = Instant::now();
        }

        // TODO: reset container instead of recreating it every time.
        let mut container = AggregatorContainer::new(chunk.to_vec(), aggregator_use, rts.clone());

        container.add_query(index, tolerance);

        // Serialize each result in the batch
        container.serialize_to_seq(&mut seq);

        if let Some(batch_start) = batch_start {
            let batch_elapsed = batch_start.elapsed();
            info!(
                "Batch {}/{} completed in {:.2?}",
                batch_idx + 1,
                total_batches,
                batch_elapsed,
            );
        }
    }

    seq.end().unwrap();

    let serialization_elapsed = serialization_start.elapsed();
    println!("Wrote to {}", output_path.display());
    println!(
        "Total processing and serialization took {:#?}",
        serialization_elapsed
    );
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ElutionGroupInput {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursors: Vec<f64>,
    pub fragments: Vec<f64>,
}

impl Into<ElutionGroup<usize>> for ElutionGroupInput {
    fn into(self) -> ElutionGroup<usize> {
        let precursors: Arc<[(i8, f64)]> = Arc::from(
            self.precursors
                .into_iter()
                .enumerate()
                .map(|(i, mz)| (i as i8, mz))
                .collect::<Vec<(i8, f64)>>(),
        );
        let fragments: Arc<[(usize, f64)]> = Arc::from(
            self.fragments
                .into_iter()
                .enumerate()
                .map(|(i, mz)| (i, mz))
                .collect::<Vec<(usize, f64)>>(),
        );
        ElutionGroup {
            id: self.id,
            mobility: self.mobility,
            rt_seconds: self.rt_seconds,
            precursors,
            fragments,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct SpectrumOutput {
    id: u64,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mzs: Vec<f64>,
    fragment_mzs: Vec<f64>,
    precursor_intensities: Vec<f32>,
    fragment_intensities: Vec<f32>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ChromatogramOutput {
    id: u64,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mzs: Vec<f64>,
    fragment_mzs: Vec<f64>,
    precursor_intensities: Vec<Vec<f32>>,
    fragment_intensities: Vec<Vec<f32>>,
    retention_time_results_seconds: Vec<f32>,
}

impl TryInto<ChromatogramOutput> for ChromatogramCollector<usize, f32> {
    type Error = ();

    fn try_into(mut self) -> Result<ChromatogramOutput, Self::Error> {
        let mut non_zero_min_idx = self.ref_rt_ms.len();
        let mut non_zero_max_idx = 0usize;

        self.iter_mut_precursors().for_each(|((_idx, _mz), cmg)| {
            let slc = cmg.as_slice();
            for (i, &inten) in slc.iter().enumerate() {
                if inten > 0.0 {
                    let abs_idx = i;
                    if abs_idx < non_zero_min_idx {
                        non_zero_min_idx = abs_idx;
                    }
                    if abs_idx > non_zero_max_idx {
                        non_zero_max_idx = abs_idx;
                    }
                }
            }
        });

        self.iter_mut_fragments().for_each(|((_idx, _mz), cmg)| {
            let slc = cmg.as_slice();
            for (i, &inten) in slc.iter().enumerate() {
                if inten > 0.0 {
                    let abs_idx = i;
                    if abs_idx < non_zero_min_idx {
                        non_zero_min_idx = abs_idx;
                    }
                    if abs_idx > non_zero_max_idx {
                        non_zero_max_idx = abs_idx;
                    }
                }
            }
        });

        if non_zero_min_idx > non_zero_max_idx {
            // No non zero intensities found
            return Err(());
        }

        let (precursor_mzs, precursor_intensities): (Vec<f64>, Vec<Vec<f32>>) = self
            .iter_mut_precursors()
            .map(|(&(idx, mz), cmg)| {
                let out_vec = match cmg.try_get_slice(non_zero_min_idx, non_zero_max_idx + 1) {
                    Some(slc) => slc.to_vec(),
                    None => {
                        error!(
                            "Failed to get slice for precursor mz {} in chromatogram id {}",
                            mz, idx,
                        );
                        panic!(
                            "Failed to get slice for precursor mz {} in chromatogram id {}, len = {}, min_idx = {}, max_idx = {}",
                            mz, idx,
                            cmg.as_slice().len(),
                            non_zero_min_idx,
                            non_zero_max_idx,
                        );
                    }
                };
                if out_vec.iter().all(|&x| x == 0.0) {
                    return None;
                }
                Some((mz, out_vec))
            })
            .flatten()
            .unzip();
        let (fragment_mzs, fragment_intensities): (Vec<f64>, Vec<Vec<f32>>) = self
            .iter_mut_fragments()
            .map(|(&(idx, mz), cmg)| {
                let out_vec = match cmg.try_get_slice(non_zero_min_idx, non_zero_max_idx + 1) {
                    Some(slc) => slc.to_vec(),
                    None => {
                        error!(
                            "Failed to get slice for fragment mz {} in chromatogram id {}",
                            mz, idx,
                        );
                        panic!(
                            "Failed to get slice for fragment mz {} in chromatogram id {}, len = {}, min_idx = {}, max_idx = {}",
                            mz, idx,
                            cmg.as_slice().len(),
                            non_zero_min_idx,
                            non_zero_max_idx,
                        );
                    }
                };
                if out_vec.iter().all(|&x| x == 0.0) {
                    return None;
                }

                Some((
                    mz,
                    out_vec,
                ))
            })
            .flatten()
            .unzip();

        Ok(ChromatogramOutput {
            id: self.eg.id,
            mobility_ook0: self.eg.mobility,
            rt_seconds: self.eg.rt_seconds,
            precursor_mzs,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
            retention_time_results_seconds: self.ref_rt_ms[non_zero_min_idx..=non_zero_max_idx]
                .iter()
                .map(|&x| x as f32 / 1000.0)
                .collect(),
        })
    }
}

impl From<&SpectralCollector<usize, f32>> for SpectrumOutput {
    fn from(agg: &SpectralCollector<usize, f32>) -> Self {
        let (precursor_mzs, precursor_intensities): (Vec<f64>, Vec<f32>) = agg
            .iter_precursors()
            .map(|((_idx, mz), inten)| (mz, inten))
            .unzip();
        let (fragment_mzs, fragment_intensities) = agg
            .iter_fragments()
            .map(|((_idx, mz), inten)| (mz, inten))
            .unzip();

        SpectrumOutput {
            id: agg.eg.id,
            mobility_ook0: agg.eg.mobility,
            rt_seconds: agg.eg.rt_seconds,
            precursor_mzs,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
        }
    }
}

fn template_tolerance_settings() -> Tolerance {
    Tolerance {
        ms: MzTolerance::Ppm((15.0, 15.0)),
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
    Point(Vec<PointIntensityAggregator<usize>>),
    Chromatogram(Vec<ChromatogramCollector<usize, f32>>),
    Spectrum(Vec<SpectralCollector<usize, f32>>),
}

impl AggregatorContainer {
    fn new(
        queries: Vec<ElutionGroup<usize>>,
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

    #[instrument(skip_all, level = "debug")]
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

    fn serialize_to_seq<S>(&mut self, seq: &mut S)
    where
        S: SerializeSeq,
    {
        match self {
            AggregatorContainer::Point(aggregators) => {
                for agg in aggregators {
                    seq.serialize_element(agg).unwrap();
                }
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                let converted: Vec<_> = aggregators
                    .par_drain(..)
                    .map(|agg| {
                        let agg_id = agg.eg.id;
                        let ser_agg: Result<ChromatogramOutput, _> = agg.try_into();
                        match ser_agg {
                            Ok(ser_agg) => Some(ser_agg),
                            Err(_) => {
                                warn!(
                                    "Skipping empty chromatogram for elution group id {}",
                                    agg_id,
                                );
                                // Skip empty chromatograms
                                None
                            }
                        }
                    })
                    .flatten()
                    .collect();

                for ser_agg in converted.iter() {
                    seq.serialize_element(&ser_agg).unwrap();
                }
            }
            AggregatorContainer::Spectrum(aggregators) => {
                for agg in aggregators.iter() {
                    let ser_agg = SpectrumOutput::from(agg);
                    seq.serialize_element(&ser_agg).unwrap();
                }
            }
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum SerializationFormat {
    Json,
    #[default]
    PrettyJson,
    // Ndjson,
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

    /// The aggregator to use.
    #[arg(short, long, default_value_t, value_enum)]
    aggregator: PossibleAggregator,

    /// Batch size for streaming serialization (default: 500)
    #[arg(short, long, default_value_t = 500)]
    batch_size: usize,
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
