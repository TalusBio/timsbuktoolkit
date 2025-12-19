use std::fs::File;
use std::io::BufWriter;
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};

use serde::ser::{
    SerializeSeq,
    Serializer,
};
use timscentroid::{
    IndexedTimstofPeaks,
    TimsTofPath,
};
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
    Tolerance,
};
use timsquery::serde::load_index_caching;
use timsrust::MSLevel;
use tracing::{
    info,
    instrument,
    warn,
};

use crate::cli::{
    PossibleAggregator,
    QueryIndexArgs,
    SerializationFormat,
    WriteTemplateArgs,
};
use crate::error::CliError;
use crate::processing::AggregatorContainer;

/// Main function for the 'query-index' subcommand.
#[instrument]
pub fn main_query_index(args: QueryIndexArgs) -> Result<(), CliError> {
    let raw_file_path = args.raw_file_path;
    let tolerance_settings_path = args.tolerance_settings_path;
    let elution_groups_path = args.elution_groups_path;
    let aggregator_use = args.aggregator;

    let tolerance_settings: Tolerance =
        serde_json::from_str(&std::fs::read_to_string(&tolerance_settings_path)?)?;
    info!("Using tolerance settings: {:#?}", tolerance_settings);
    info!(
        "Loading elution groups from {}",
        elution_groups_path.display()
    );
    let elution_groups: Vec<TimsElutionGroup<u8>> =
        read_query_elution_groups(&elution_groups_path)?;
    info!("Loaded {} elution groups", elution_groups.len());

    let index = load_index_caching(&raw_file_path)
        .map_err(|e| CliError::DataReading(format!("{:?}", e)))?;
    let rts = get_ms1_rts_as_millis(&raw_file_path)?;

    let output_path = args.output_path;
    let serialization_format = args.format;
    let batch_size = args.batch_size;

    std::fs::create_dir_all(&output_path)?;
    let put_path = output_path.join("results.json");

    stream_process_batches(
        elution_groups,
        aggregator_use,
        rts,
        &index,
        &tolerance_settings,
        serialization_format,
        &put_path,
        batch_size,
    )?;
    Ok(())
}

/// Reads elution groups from a given path, attempting to parse them in several formats.
pub fn read_query_elution_groups(path: &PathBuf) -> Result<Vec<TimsElutionGroup<u8>>, CliError> {
    match timsquery::serde::read_library_file(path) {
        Ok(timsquery::serde::ElutionGroupCollection::TinyIntLabels(egs, _)) => Ok(egs),
        Ok(other) => Err(CliError::DataReading(format!(
            "Expected elution groups with u8 labels, but got different label type: {:?}",
            other
        ))),
        Err(e) => Err(CliError::DataReading(format!(
            "Failed to read elution groups from {}: {:?}",
            path.display(),
            e
        ))),
    }
}

/// Main function for the 'write-template' subcommand.
pub fn main_write_template(args: WriteTemplateArgs) -> Result<(), CliError> {
    todo!("Ooops, this feature is not yet implemented.");
}

/// Retrieves MS1 retention times from a TIMS-TOF file, sorted and deduped.
pub fn get_ms1_rts_as_millis(file: &PathBuf) -> Result<Arc<[u32]>, CliError> {
    let ttp = TimsTofPath::new(file).map_err(|e| CliError::TimsFileLoad {
        path: file.clone(),
        source: e,
    })?;
    let reader = ttp.load_frame_reader()?;
    let mut rts: Vec<_> = reader
        .frame_metas
        .iter()
        .filter(|x| x.ms_level == MSLevel::MS1)
        .map(|f| (f.rt_in_seconds * 1000.0).round() as u32)
        .collect();
    rts.sort_unstable();
    rts.dedup();
    Ok(rts.into())
}

/// Streams and processes elution groups in batches, then serializes the results.
#[instrument(skip_all)]
pub fn stream_process_batches(
    elution_groups: Vec<TimsElutionGroup<u8>>,
    aggregator_use: PossibleAggregator,
    rts: Arc<[u32]>,
    index: &IndexedTimstofPeaks,
    tolerance: &Tolerance,
    serialization_format: SerializationFormat,
    output_path: &Path,
    batch_size: usize,
) -> Result<(), CliError> {
    let total_groups = elution_groups.len();
    let total_batches = total_groups.div_ceil(batch_size);

    info!(
        "Processing {} elution groups in {} batches of up to {}",
        total_groups, total_batches, batch_size
    );

    let serialization_start = Instant::now();
    let file = File::create(output_path)?;
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
            )?;
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
            )?;
        }
    }
    Ok(())
}

/// Processes batches of elution groups and serializes the aggregated results.
#[instrument(skip_all)]
pub fn process_and_serialize<S>(
    elution_groups: Vec<TimsElutionGroup<u8>>,
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
) -> Result<(), CliError>
where
    S: Serializer,
{
    let mut seq = ser.serialize_seq(Some(total_groups)).unwrap();

    let mut last_progress = Instant::now();
    let progress_interval = Duration::from_secs(2);

    for (batch_idx, chunk) in elution_groups.chunks(batch_size).enumerate() {
        if last_progress.elapsed() >= progress_interval {
            info!(
                "Processing batch {}/{} ({} groups)",
                batch_idx + 1,
                total_batches,
                chunk.len(),
            );
            last_progress = Instant::now();
        }

        let mut container =
            AggregatorContainer::new(chunk.to_vec(), aggregator_use, rts.clone(), tolerance)?;

        container.add_query(index, tolerance);

        container.serialize_to_seq(&mut seq, &rts)?;
    }

    seq.end().unwrap();

    let serialization_elapsed = serialization_start.elapsed();
    println!("Wrote to {}", output_path.display());
    println!(
        "Total processing and serialization took {:#?}",
        serialization_elapsed
    );
    Ok(())
}

/// Generates a default set of tolerance settings.
pub fn template_tolerance_settings() -> Tolerance {
    Tolerance {
        ms: MzTolerance::Ppm((15.0, 15.0)),
        rt: RtTolerance::Unrestricted,
        mobility: MobilityTolerance::Pct((10.0, 10.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    }
}
