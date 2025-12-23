use serde::Serialize;
use std::fmt::Display;
use std::fs::File;
use std::io::{
    self,
    BufWriter,
    Write,
};
use std::path::{
    Path,
    PathBuf,
};
use std::time::{
    Duration,
    Instant,
};

use timscentroid::IndexedTimstofPeaks;
use timsquery::KeyLike;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::load_index_auto;
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
use timsquery::serde::ElutionGroupCollection;

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
    let elution_groups: ElutionGroupCollection = read_query_elution_groups(&elution_groups_path)?;
    info!("Loaded {} elution groups", elution_groups.len());

    let index = load_index_auto(
        raw_file_path.to_str().ok_or_else(|| {
            CliError::DataReading("Invalid path encoding".to_string())
        })?,
        None,  // Use default config
    )
    .map_err(|e| CliError::DataReading(format!("{:?}", e)))?
    .into_eager()
    .map_err(|e| CliError::DataReading(format!("{:?}", e)))?;

    let output_path = args.output_path;
    let serialization_format = args.format;
    let batch_size = args.batch_size;

    std::fs::create_dir_all(&output_path)?;
    let put_path = output_path.join("results.json");

    match elution_groups {
        ElutionGroupCollection::StringLabels(egs, _) => stream_process_batches(
            egs,
            aggregator_use,
            &index,
            &tolerance_settings,
            serialization_format,
            &put_path,
            batch_size,
        ),
        ElutionGroupCollection::MzpafLabels(egs, _) => stream_process_batches(
            egs,
            aggregator_use,
            &index,
            &tolerance_settings,
            serialization_format,
            &put_path,
            batch_size,
        ),
        ElutionGroupCollection::TinyIntLabels(egs, _) => stream_process_batches(
            egs,
            aggregator_use,
            &index,
            &tolerance_settings,
            serialization_format,
            &put_path,
            batch_size,
        ),
        ElutionGroupCollection::IntLabels(egs, _) => stream_process_batches(
            egs,
            aggregator_use,
            &index,
            &tolerance_settings,
            serialization_format,
            &put_path,
            batch_size,
        ),
    }?;
    Ok(())
}

/// Reads elution groups from a given path, attempting to parse them in several formats.
pub fn read_query_elution_groups(path: &PathBuf) -> Result<ElutionGroupCollection, CliError> {
    match timsquery::serde::read_library_file(path) {
        Ok(egs) => Ok(egs),
        Err(e) => Err(CliError::DataReading(format!(
            "Failed to read elution groups from {}: {:?}",
            path.display(),
            e
        ))),
    }
}

const NARROW_TOLERANCE_TEMPLATE: &str = r#"{
  "ms": { "ppm": [10.0, 10.0] },
  "rt": { "minutes": [0.5, 0.5] },
  "mobility": { "pct": [5.0, 5.0] },
  "quad": { "da": [1.05, 1.05] }
}"#;

const WIDE_TOLERANCE_TEMPLATE: &str = r#"{
  "ms": { "da": [0.04, 0.04] },
  "rt": "Unrestricted",
  "mobility": { "pct": [20.0, 20.0] },
  "quad": { "da": [0.2, 0.2] }
}"#;

const ELUTION_GROUP_TEMPLATE: &str = r#"[
    {
		"fragment_labels":[ "y1", "y1^2", "y2", "y2^2", "b1", "b1^2", "p^2" ],
		"fragments":[ 147.1128, 74.06004, 248.1604, 124.58387, 347.22889, 174.11808, 418.26601 ],
		"id":0,
		"mobility":0.9851410984992981,
		"precursor_mz":723.844601280237,
		"precursor_charge":2,
		"precursor_isotopes":[0,1,2],
		"rt_seconds":302.2712
	},
    {
		"fragment_labels":[ "y1", "y1^2", "y2", "y2^2", "b1", "b1^2", "p^2" ],
		"fragments":[ 147.1128, 74.06004, 248.1604, 124.58387, 347.22889, 174.11808, 418.26601 ],
		"id":1,
		"mobility":0.9851410984992981,
		"precursor_mz":723.844601280237,
		"precursor_charge":1,
		"precursor_isotopes":[0],
		"rt_seconds":354.2712
	}
]"#;

/// Main function for the 'write-template' subcommand.
pub fn main_write_template(args: WriteTemplateArgs) -> Result<(), CliError> {
    let target_dir = args.output_path;
    std::fs::create_dir_all(&target_dir)?;

    let narrow_path = target_dir.join("narrow_tolerance_template.json");
    let wide_path = target_dir.join("wide_tolerance_template.json");
    std::fs::write(&narrow_path, NARROW_TOLERANCE_TEMPLATE)?;
    std::fs::write(&wide_path, WIDE_TOLERANCE_TEMPLATE)?;
    println!(
        "Wrote tolerance templates to:\n- {}\n- {}",
        narrow_path.display(),
        wide_path.display()
    );

    // Elution group template
    let elution_group_path = target_dir.join("elution_group_template.json");
    std::fs::write(&elution_group_path, ELUTION_GROUP_TEMPLATE)?;
    println!(
        "Wrote elution group template to: {}",
        elution_group_path.display()
    );
    Ok(())
}

pub struct JsonStreamSerializer<W: Write> {
    writer: W,
    format: SerializationFormat,
    is_first: bool,
}

impl<W: Write> JsonStreamSerializer<W> {
    pub fn new(writer: W, format: SerializationFormat) -> Self {
        Self {
            writer,
            format,
            is_first: true,
        }
    }

    /// Serializes an item based on the selected format.
    pub fn serialize<T: Serialize>(&mut self, item: &T) -> io::Result<()> {
        match self.format {
            SerializationFormat::Ndjson => {
                serde_json::to_writer(&mut self.writer, item).map_err(io::Error::other)?;
                self.writer.write_all(b"\n")?;
            }
            SerializationFormat::Json | SerializationFormat::PrettyJson => {
                if self.is_first {
                    self.writer.write_all(b"[")?;
                    self.is_first = false;
                } else {
                    self.writer.write_all(b",")?;
                }

                if matches!(self.format, SerializationFormat::PrettyJson) {
                    serde_json::to_writer_pretty(&mut self.writer, item)
                } else {
                    serde_json::to_writer(&mut self.writer, item)
                }
                .map_err(io::Error::other)?;
            }
        }
        Ok(())
    }

    /// Finalizes the output (crucial for closing JSON arrays).
    pub fn finish(mut self) -> io::Result<()> {
        match self.format {
            SerializationFormat::Json | SerializationFormat::PrettyJson => {
                if self.is_first {
                    // Handle case where no items were ever serialized
                    self.writer.write_all(b"[]")?;
                } else {
                    self.writer.write_all(b"]")?;
                }
            }
            SerializationFormat::Ndjson => {} // No closing tag needed
        }
        self.writer.flush()
    }
}

/// Streams and processes elution groups in batches, then serializes the results.
#[instrument(skip_all)]
pub fn stream_process_batches<T: KeyLike + Display>(
    elution_groups: Vec<TimsElutionGroup<T>>,
    aggregator_use: PossibleAggregator,
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
            let ser = JsonStreamSerializer::new(writer, SerializationFormat::PrettyJson);
            process_and_serialize(
                elution_groups,
                aggregator_use,
                index,
                tolerance,
                ser,
                total_batches,
                batch_size,
                serialization_start,
                output_path,
            )?;
        }
        SerializationFormat::Json => {
            let ser = JsonStreamSerializer::new(writer, SerializationFormat::Json);
            process_and_serialize(
                elution_groups,
                aggregator_use,
                index,
                tolerance,
                ser,
                total_batches,
                batch_size,
                serialization_start,
                output_path,
            )?;
        }
        SerializationFormat::Ndjson => {
            let ser = JsonStreamSerializer::new(writer, SerializationFormat::Ndjson);
            process_and_serialize(
                elution_groups,
                aggregator_use,
                index,
                tolerance,
                ser,
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
pub fn process_and_serialize<T: KeyLike + Display>(
    elution_groups: Vec<TimsElutionGroup<T>>,
    aggregator_use: PossibleAggregator,
    index: &IndexedTimstofPeaks,
    tolerance: &Tolerance,
    mut ser: JsonStreamSerializer<impl Write>,
    total_batches: usize,
    batch_size: usize,
    serialization_start: Instant,
    output_path: &Path,
) -> Result<(), CliError> {
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

        let mut container = AggregatorContainer::new(
            chunk.to_vec(),
            aggregator_use,
            index.ms1_cycle_mapping(),
            tolerance,
        )?;

        container.add_query(index, tolerance);
        container.serialize_to_seq(&mut ser, index.ms1_cycle_mapping())?;
    }

    ser.finish()?;
    let serialization_elapsed = serialization_start.elapsed();
    println!("Wrote to {}", output_path.display());
    println!(
        "Total processing and serialization took {:#?}",
        serialization_elapsed
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::models::tolerance::{
        MobilityTolerance,
        MzTolerance,
        QuadTolerance,
        RtTolerance,
        Tolerance,
    };

    #[test]
    fn test_we_can_read_data_contract() {
        let manifest_path = env!("CARGO_MANIFEST_DIR");
        let elution_groups_path =
            PathBuf::from(manifest_path).join("data_contracts/single_elution_group.json");

        let elution_groups = read_query_elution_groups(&elution_groups_path).unwrap();
        // Do not change that file, it means its a data contract test
        // AKA, we promised we would be compatible with that file.
        // IF you do want to change it contact the Carafe developers first.
        assert!(elution_groups.len() == 1);
    }

    #[test]
    fn test_templates_tolerance_deserializable() {
        let narrow: Tolerance = serde_json::from_str(NARROW_TOLERANCE_TEMPLATE).unwrap();
        let wide: Tolerance = serde_json::from_str(WIDE_TOLERANCE_TEMPLATE).unwrap();

        assert!(matches!(narrow.ms, MzTolerance::Ppm(_)));
        assert!(matches!(narrow.rt, RtTolerance::Minutes(_)));
        assert!(matches!(narrow.mobility, MobilityTolerance::Pct(_)));
        assert!(matches!(narrow.quad, QuadTolerance::Absolute(_)));

        assert!(matches!(wide.ms, MzTolerance::Absolute(_)));
        assert!(matches!(wide.rt, RtTolerance::Unrestricted));
        assert!(matches!(wide.mobility, MobilityTolerance::Pct(_)));
        assert!(matches!(wide.quad, QuadTolerance::Absolute(_)));
    }

    #[test]
    fn test_elution_group_template_deserializable() {
        use timsquery::IonAnnot;
        let elution_groups =
            serde_json::from_str::<Vec<TimsElutionGroup<IonAnnot>>>(ELUTION_GROUP_TEMPLATE)
                .unwrap();
        assert!(elution_groups.len() == 2);
        // Write to a temp file ... while I implement direct reading api
        let tmp_file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp_file.path(), ELUTION_GROUP_TEMPLATE).unwrap();
        let elution_groups = read_query_elution_groups(&tmp_file.path().to_path_buf()).unwrap();
        assert!(elution_groups.len() == 2);
    }
}
