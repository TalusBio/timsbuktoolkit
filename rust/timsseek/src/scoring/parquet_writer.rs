use arrow::array::*;
use arrow::datatypes::*;
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use tracing::debug;

use super::results::{FinalResult, ScoringFields};

// ---------------------------------------------------------------------------
// Macro: build typed columns from accessor closures
// ---------------------------------------------------------------------------

macro_rules! columns {
    ($results:expr; $( $name:expr => $dtype:expr, $array_type:ident ( $accessor:expr ) );* $(;)? ) => {{
        let mut fields: Vec<Field> = Vec::new();
        let mut arrays: Vec<Arc<dyn Array>> = Vec::new();
        $(
            fields.push(Field::new($name, $dtype, true));
            arrays.push(Arc::new(<$array_type>::from_iter(
                $results.iter().map($accessor)
            )));
        )*
        (fields, arrays)
    }};
}

// ---------------------------------------------------------------------------
// Helper: expand fixed-size [f32; N] arrays into numbered columns
// ---------------------------------------------------------------------------

fn expand_f32_array_columns<const N: usize>(
    results: &[FinalResult],
    prefix: &str,
    accessor: impl Fn(&FinalResult) -> &[f32; N],
    fields: &mut Vec<Field>,
    arrays: &mut Vec<Arc<dyn Array>>,
) {
    for i in 0..N {
        fields.push(Field::new(
            format!("{}_{}", prefix, i),
            DataType::Float32,
            true,
        ));
        arrays.push(Arc::new(Float32Array::from_iter(
            results.iter().map(|r| Some(accessor(r)[i])),
        )));
    }
}

// ---------------------------------------------------------------------------
// Build a RecordBatch from a slice of FinalResult
// ---------------------------------------------------------------------------

/// Build a `RecordBatch` from a slice of `FinalResult`.
///
/// **COMPILE-TIME SAFETY:** The exhaustive destructure at the top ensures that
/// adding a field to `ScoringFields` or `FinalResult` without updating this
/// function causes a compile error.
pub fn build_record_batch(results: &[FinalResult]) -> RecordBatch {
    // -----------------------------------------------------------------------
    // Exhaustive destructure for compile-time completeness check.
    // Every field must be listed -- no `..` allowed.
    // The bindings are prefixed with _ to suppress unused warnings;
    // the actual column building uses accessor closures below.
    // -----------------------------------------------------------------------
    if let Some(r) = results.first() {
        let FinalResult {
            scoring:
                ScoringFields {
                    sequence: _,
                    library_id: _,
                    decoy_group_id: _,
                    precursor_mz: _,
                    precursor_charge: _,
                    precursor_mobility: _,
                    is_target: _,
                    library_rt: _,
                    calibrated_rt_seconds: _,
                    obs_rt_seconds: _,
                    calibrated_sq_delta_rt: _,
                    obs_mobility: _,
                    delta_ms1_ms2_mobility: _,
                    sq_delta_ms1_ms2_mobility: _,
                    main_score: _,
                    delta_next: _,
                    delta_second_next: _,
                    apex_lazyscore: _,
                    ms2_lazyscore: _,
                    ms2_isotope_lazyscore: _,
                    ms2_isotope_lazyscore_ratio: _,
                    lazyscore_z: _,
                    lazyscore_vs_baseline: _,
                    split_product_score: _,
                    cosine_au: _,
                    scribe_au: _,
                    cosine_cg: _,
                    scribe_cg: _,
                    cosine_weighted_coelution: _,
                    cosine_gradient_consistency: _,
                    scribe_weighted_coelution: _,
                    scribe_gradient_consistency: _,
                    peak_shape: _,
                    ratio_cv: _,
                    centered_apex: _,
                    precursor_coelution: _,
                    fragment_coverage: _,
                    precursor_apex_match: _,
                    xic_quality: _,
                    fragment_apex_agreement: _,
                    isotope_correlation: _,
                    gaussian_correlation: _,
                    per_frag_gaussian_corr: _,
                    rising_cycles: _,
                    falling_cycles: _,
                    npeaks: _,
                    n_scored_fragments: _,
                    ms2_summed_intensity: _,
                    ms1_summed_intensity: _,
                    ms2_mz_errors: _,
                    ms2_mobility_errors: _,
                    ms1_mz_errors: _,
                    ms1_mobility_errors: _,
                    ms2_intensity_ratios: _,
                    ms1_intensity_ratios: _,
                },
            delta_group: _,
            delta_group_ratio: _,
            discriminant_score: _,
            qvalue: _,
        } = r;
    }

    // -----------------------------------------------------------------------
    // Scalar columns via the `columns!` macro
    // -----------------------------------------------------------------------
    let (mut fields, mut arrays) = columns!(results;
        // Identity
        "sequence"         => DataType::Utf8,    StringArray(|r: &FinalResult| Some(r.scoring.sequence.as_str()));
        "library_id"       => DataType::UInt32,  UInt32Array(|r: &FinalResult| Some(r.scoring.library_id));
        "decoy_group_id"   => DataType::UInt32,  UInt32Array(|r: &FinalResult| Some(r.scoring.decoy_group_id));
        "precursor_mz"     => DataType::Float64, Float64Array(|r: &FinalResult| Some(r.scoring.precursor_mz));
        "precursor_charge" => DataType::UInt8,   UInt8Array(|r: &FinalResult| Some(r.scoring.precursor_charge));
        "precursor_mobility" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.precursor_mobility));
        "is_target"        => DataType::Boolean, BooleanArray(|r: &FinalResult| Some(r.scoring.is_target));

        // RT
        "library_rt"             => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.library_rt));
        "calibrated_rt_seconds"  => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.calibrated_rt_seconds));
        "obs_rt_seconds"         => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.obs_rt_seconds));
        "calibrated_sq_delta_rt" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.calibrated_sq_delta_rt));

        // Mobility
        "obs_mobility"              => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.obs_mobility));
        "delta_ms1_ms2_mobility"    => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.delta_ms1_ms2_mobility));
        "sq_delta_ms1_ms2_mobility" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.sq_delta_ms1_ms2_mobility));

        // Primary scores
        "main_score"          => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.main_score));
        "delta_next"          => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.delta_next));
        "delta_second_next"   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.delta_second_next));

        // Lazyscores
        "apex_lazyscore"              => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.apex_lazyscore));
        "ms2_lazyscore"               => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.ms2_lazyscore));
        "ms2_isotope_lazyscore"       => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.ms2_isotope_lazyscore));
        "ms2_isotope_lazyscore_ratio" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.ms2_isotope_lazyscore_ratio));
        "lazyscore_z"                 => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.lazyscore_z));
        "lazyscore_vs_baseline"       => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.lazyscore_vs_baseline));

        // Split product (9 components)
        "split_product_score"         => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.split_product_score));
        "cosine_au"                   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.cosine_au));
        "scribe_au"                   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.scribe_au));
        "cosine_cg"                   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.cosine_cg));
        "scribe_cg"                   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.scribe_cg));
        "cosine_weighted_coelution"   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.cosine_weighted_coelution));
        "cosine_gradient_consistency" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.cosine_gradient_consistency));
        "scribe_weighted_coelution"   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.scribe_weighted_coelution));
        "scribe_gradient_consistency" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.scribe_gradient_consistency));

        // 11 apex features
        "peak_shape"               => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.peak_shape));
        "ratio_cv"                 => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.ratio_cv));
        "centered_apex"            => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.centered_apex));
        "precursor_coelution"      => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.precursor_coelution));
        "fragment_coverage"        => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.fragment_coverage));
        "precursor_apex_match"     => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.precursor_apex_match));
        "xic_quality"              => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.xic_quality));
        "fragment_apex_agreement"  => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.fragment_apex_agreement));
        "isotope_correlation"      => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.isotope_correlation));
        "gaussian_correlation"     => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.gaussian_correlation));
        "per_frag_gaussian_corr"   => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.per_frag_gaussian_corr));

        // Peak shape
        "rising_cycles"  => DataType::UInt8, UInt8Array(|r: &FinalResult| Some(r.scoring.rising_cycles));
        "falling_cycles" => DataType::UInt8, UInt8Array(|r: &FinalResult| Some(r.scoring.falling_cycles));

        // Counts
        "npeaks"             => DataType::UInt8, UInt8Array(|r: &FinalResult| Some(r.scoring.npeaks));
        "n_scored_fragments" => DataType::UInt8, UInt8Array(|r: &FinalResult| Some(r.scoring.n_scored_fragments));

        // Intensities
        "ms2_summed_intensity" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.ms2_summed_intensity));
        "ms1_summed_intensity" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.scoring.ms1_summed_intensity));

        // FinalResult-level fields
        "delta_group"       => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.delta_group));
        "delta_group_ratio" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.delta_group_ratio));
        "discriminant_score" => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.discriminant_score));
        "qvalue"            => DataType::Float32, Float32Array(|r: &FinalResult| Some(r.qvalue))
    );

    // -----------------------------------------------------------------------
    // Array columns: expand [f32; N] into numbered columns
    // -----------------------------------------------------------------------
    expand_f32_array_columns(results, "ms2_mz_error",       |r| &r.scoring.ms2_mz_errors,       &mut fields, &mut arrays);
    expand_f32_array_columns(results, "ms2_mobility_error",  |r| &r.scoring.ms2_mobility_errors,  &mut fields, &mut arrays);
    expand_f32_array_columns(results, "ms1_mz_error",       |r| &r.scoring.ms1_mz_errors,       &mut fields, &mut arrays);
    expand_f32_array_columns(results, "ms1_mobility_error",  |r| &r.scoring.ms1_mobility_errors,  &mut fields, &mut arrays);
    expand_f32_array_columns(results, "ms2_intensity_ratio", |r| &r.scoring.ms2_intensity_ratios, &mut fields, &mut arrays);
    expand_f32_array_columns(results, "ms1_intensity_ratio", |r| &r.scoring.ms1_intensity_ratios, &mut fields, &mut arrays);

    // -----------------------------------------------------------------------
    // Build the RecordBatch
    // -----------------------------------------------------------------------
    let schema = Arc::new(Schema::new(fields));
    RecordBatch::try_new(schema, arrays).expect("schema/array mismatch in build_record_batch")
}

// ---------------------------------------------------------------------------
// Buffered Parquet writer
// ---------------------------------------------------------------------------

pub struct ResultParquetWriter {
    writer: ArrowWriter<File>,
    buffer: Vec<FinalResult>,
    row_group_size: usize,
}

impl ResultParquetWriter {
    pub fn new(path: impl AsRef<Path>, row_group_size: usize) -> std::io::Result<Self> {
        let file = match File::create_new(path.as_ref()) {
            Ok(f) => f,
            Err(err) => {
                tracing::error!(
                    "Failed to create file {:?}: {}",
                    path.as_ref(),
                    err
                );
                return Err(err);
            }
        };

        // Build schema from a zero-row batch
        let empty_batch = build_record_batch(&[]);
        let schema = empty_batch.schema();

        let props = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .build();

        let writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        Ok(Self {
            writer,
            buffer: Vec::with_capacity(row_group_size),
            row_group_size,
        })
    }

    pub fn add(&mut self, result: FinalResult) {
        self.buffer.push(result);
        if self.buffer.len() >= self.row_group_size {
            self.flush();
        }
    }

    fn flush(&mut self) {
        if self.buffer.is_empty() {
            return;
        }
        debug!("Flushing {} results to parquet", self.buffer.len());
        let batch = build_record_batch(&self.buffer);
        self.writer.write(&batch).expect("parquet write failed");
        self.buffer.clear();
    }

    pub fn close(mut self) {
        self.flush();
        self.writer.close().expect("parquet close failed");
    }
}
