//! High-performance serialization for indexed timsTOF peaks.
//!
//! The only struct here that actually matters is [`IndexedTimstofPeaks`]
//!
//! # Examples
//!
//! ```no_run
//! use timscentroid::{IndexedTimstofPeaks, CentroidingConfig, StorageLocation};
//! use timscentroid::serialization::SerializationConfig;
//! use timsrust::TimsTofPath;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Index peaks from raw data
//! let file = TimsTofPath::new("data.d")?;
//! let config = CentroidingConfig::default();
//! let (index, _) = IndexedTimstofPeaks::from_timstof_file(&file, config);
//!
//! // Save with default settings (balanced speed/size)
//! index.save_to_directory("indexed_peaks/")?;
//!
//! // Load from disk (much faster than re-indexing)
//! let location = StorageLocation::from_path("indexed_peaks/");
//! let loaded = IndexedTimstofPeaks::load_from_storage(location)?;
//!
//! # Ok(())
//! # }
//! ```

use crate::geometry::QuadrupoleIsolationScheme;
use crate::indexing::{
    IndexedPeak,
    IndexedPeakGroup,
    IndexedTimstofPeaks,
};
use crate::rt_mapping::{
    CycleToRTMapping,
    MS1CycleIndex,
    RTIndex,
    WindowCycleIndex,
};
use crate::storage::{
    StorageLocation,
    StorageProvider,
};
use arrow::array::{
    Array,
    Float16Array,
    Float32Array,
    UInt32Array,
};
use arrow::datatypes::{
    DataType,
    Field,
    Schema,
};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::{
    Compression,
    Encoding,
    // ZstdLevel, // Sometimes I want to use zstd level 1.
};
use parquet::file::properties::{
    EnabledStatistics,
    WriterProperties,
};
use parquet::schema::types::ColumnPath;
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Errors that can occur during serialization/deserialization
#[derive(Debug)]
pub enum SerializationError {
    Io(std::io::Error),
    Parquet(parquet::errors::ParquetError),
    Arrow(arrow::error::ArrowError),
    Json(serde_json::Error),
    ObjectStore(object_store::Error),
    UrlParse(url::ParseError),
    MissingColumn {
        name: String,
        available: Vec<String>,
    },
    WrongColumnType {
        column: String,
        expected: &'static str,
        got: String,
    },
    SchemaVersionMismatch {
        expected: &'static str,
        found: String,
    },
}

impl fmt::Display for SerializationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {}", e),
            Self::Parquet(e) => write!(f, "Parquet error: {}", e),
            Self::Arrow(e) => write!(f, "Arrow error: {}", e),
            Self::Json(e) => write!(f, "JSON error: {}", e),
            Self::ObjectStore(e) => write!(f, "Object store error: {}", e),
            Self::UrlParse(e) => write!(f, "URL parse error: {}", e),
            Self::MissingColumn { name, available } => {
                write!(
                    f,
                    "Missing column '{}'. Available columns: {:?}",
                    name, available
                )
            }
            Self::WrongColumnType {
                column,
                expected,
                got,
            } => {
                write!(
                    f,
                    "Column '{}' has wrong type. Expected {}, got {}",
                    column, expected, got
                )
            }
            Self::SchemaVersionMismatch { expected, found } => {
                write!(
                    f,
                    "Schema version mismatch. Expected {}, found {}",
                    expected, found
                )
            }
        }
    }
}

impl std::error::Error for SerializationError {}

impl From<std::io::Error> for SerializationError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<parquet::errors::ParquetError> for SerializationError {
    fn from(e: parquet::errors::ParquetError) -> Self {
        Self::Parquet(e)
    }
}

impl From<arrow::error::ArrowError> for SerializationError {
    fn from(e: arrow::error::ArrowError) -> Self {
        Self::Arrow(e)
    }
}

impl From<serde_json::Error> for SerializationError {
    fn from(e: serde_json::Error) -> Self {
        Self::Json(e)
    }
}

impl From<object_store::Error> for SerializationError {
    fn from(e: object_store::Error) -> Self {
        Self::ObjectStore(e)
    }
}

impl From<url::ParseError> for SerializationError {
    fn from(e: url::ParseError) -> Self {
        Self::UrlParse(e)
    }
}

const SCHEMA_VERSION: &str = "1.0";

/// Configuration for parquet serialization
#[derive(Debug, Clone, Copy)]
pub struct SerializationConfig {
    /// Compression algorithm and level
    pub compression: Compression,
    // Number of rows per row group (affects memory usage and parallelism)
    pub row_group_size: usize,
    /// Write batch size for internal buffering, the square of this is the row group size
    pub write_batch_size: usize,
}

impl Default for SerializationConfig {
    fn default() -> Self {
        // Optimized through extensive benchmarking across compression algorithms,
        // levels, and row group sizes. This configuration provides the best balance
        // of read/write speed and disk space for typical MS data (100M+ peaks).
        Self {
            // compression: Compression::ZSTD(ZstdLevel::try_new(1).expect("ZSTD level 1 is valid")),
            compression: Compression::SNAPPY,
            row_group_size: 100_000,
            write_batch_size: 8192,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct TimscentroidMetadata {
    pub version: String,
    pub created_at: String,
    pub ms1_peaks: PeakGroupMetadata<MS1CycleIndex>,
    pub ms2_window_groups: Vec<Ms2GroupMetadata>,
}

#[derive(Serialize, Deserialize)]
pub struct PeakGroupMetadata<T: RTIndex> {
    pub relative_path: PathBuf,
    pub cycle_to_rt_ms: CycleToRTMapping<T>,
    pub bucket_size: usize,
}

#[derive(Serialize, Deserialize)]
pub struct Ms2GroupMetadata {
    pub id: usize,
    pub quadrupole_isolation: QuadrupoleIsolationScheme,
    pub group_info: PeakGroupMetadata<WindowCycleIndex>,
}

/// Expected schema for IndexedPeak parquet files
pub(crate) struct PeakSchema {
    mz_idx: usize,
    intensity_idx: usize,
    mobility_idx: usize,
    cycle_idx: usize,
}

impl PeakSchema {
    pub(crate) fn mz_idx(&self) -> usize {
        self.mz_idx
    }

    #[allow(dead_code)]
    pub(crate) fn intensity_idx(&self) -> usize {
        self.intensity_idx
    }

    pub(crate) fn mobility_idx(&self) -> usize {
        self.mobility_idx
    }

    #[allow(dead_code)]
    pub(crate) fn cycle_idx(&self) -> usize {
        self.cycle_idx
    }

    /// Create the canonical Arrow schema for IndexedPeak data
    pub(crate) fn canonical() -> Schema {
        Schema::new(vec![
            Field::new("mz", DataType::Float32, false),
            Field::new("intensity", DataType::Float32, false),
            Field::new("mobility_ook0", DataType::Float16, false),
            Field::new("cycle_index", DataType::UInt32, false),
        ])
    }

    /// Validate that a schema matches our expected format and return column indices
    pub(crate) fn validate(schema: &Schema) -> Result<Self, SerializationError> {
        let field_names: Vec<String> = schema.fields().iter().map(|f| f.name().clone()).collect();

        let find_col =
            |name: &str, expected_type: &DataType| -> Result<usize, SerializationError> {
                let idx = schema
                    .fields()
                    .iter()
                    .position(|f| f.name() == name)
                    .ok_or_else(|| SerializationError::MissingColumn {
                        name: name.to_string(),
                        available: field_names.clone(),
                    })?;

                let field = &schema.fields()[idx];
                if field.data_type() != expected_type {
                    return Err(SerializationError::WrongColumnType {
                        column: name.to_string(),
                        expected: format!("{:?}", expected_type).leak(),
                        got: format!("{:?}", field.data_type()),
                    });
                }

                Ok(idx)
            };

        Ok(Self {
            mz_idx: find_col("mz", &DataType::Float32)?,
            intensity_idx: find_col("intensity", &DataType::Float32)?,
            mobility_idx: find_col("mobility_ook0", &DataType::Float16)?,
            cycle_idx: find_col("cycle_index", &DataType::UInt32)?,
        })
    }
}

impl IndexedTimstofPeaks {
    /// Save indexed peaks to a directory (async version)
    pub async fn save_to_directory_async(
        &self,
        directory: impl AsRef<Path>,
    ) -> Result<(), SerializationError> {
        let location = StorageLocation::from_path(directory);
        self.save_to_storage_async(location, SerializationConfig::default())
            .await
    }

    /// Save indexed peaks to a directory (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `save_to_directory_async`.
    pub fn save_to_directory(&self, directory: impl AsRef<Path>) -> Result<(), SerializationError> {
        let location = StorageLocation::from_path(directory);
        self.save_to_storage(location, SerializationConfig::default())
    }

    /// Save indexed peaks with custom config (async version)
    pub async fn save_to_directory_with_config_async(
        &self,
        directory: impl AsRef<Path>,
        config: SerializationConfig,
    ) -> Result<(), SerializationError> {
        let location = StorageLocation::from_path(directory);
        self.save_to_storage_async(location, config).await
    }

    /// Save indexed peaks with custom config (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `save_to_directory_with_config_async`.
    pub fn save_to_directory_with_config(
        &self,
        directory: impl AsRef<Path>,
        config: SerializationConfig,
    ) -> Result<(), SerializationError> {
        let location = StorageLocation::from_path(directory);
        self.save_to_storage(location, config)
    }

    /// Save to cloud storage URL (async version)
    pub async fn save_to_url_async(&self, url: impl AsRef<str>) -> Result<(), SerializationError> {
        let location = StorageLocation::from_url(url)?;
        self.save_to_storage_async(location, SerializationConfig::default())
            .await
    }

    /// Save to cloud storage URL (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `save_to_url_async`.
    pub fn save_to_url(&self, url: impl AsRef<str>) -> Result<(), SerializationError> {
        let location = StorageLocation::from_url(url)?;
        self.save_to_storage(location, SerializationConfig::default())
    }

    /// Save to any storage location (async version)
    pub async fn save_to_storage_async(
        &self,
        location: StorageLocation,
        config: SerializationConfig,
    ) -> Result<(), SerializationError> {
        // Note: Currently this is not truly async - it uses blocking I/O internally.
        // This wrapper exists for API consistency and future optimization.
        self.save_to_storage_impl(location, config)
    }

    /// Save to any storage location (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `save_to_storage_async`.
    pub fn save_to_storage(
        &self,
        location: StorageLocation,
        config: SerializationConfig,
    ) -> Result<(), SerializationError> {
        self.save_to_storage_impl(location, config)
    }

    /// Internal implementation of save_to_storage
    fn save_to_storage_impl(
        &self,
        location: StorageLocation,
        config: SerializationConfig,
    ) -> Result<(), SerializationError> {
        let storage = StorageProvider::new(location)?;

        // Ensure directories exist (no-op for cloud, creates dirs for local)
        storage.ensure_directory("")?;
        storage.ensure_directory("ms2")?;

        // Parallel write: MS1 and all MS2 groups written concurrently
        let (ms1_result, ms2_results): (
            Result<_, SerializationError>,
            Vec<Result<_, SerializationError>>,
        ) = rayon::join(
            // Write MS1 in parallel thread
            || {
                let ms1_bytes = write_peaks_to_parquet_bytes(&self.ms1_peaks.peaks, config)?;
                storage.write_bytes("ms1.parquet", ms1_bytes)?;

                Ok(PeakGroupMetadata {
                    relative_path: PathBuf::from("ms1.parquet"),
                    cycle_to_rt_ms: self.ms1_peaks.cycle_to_rt_ms.clone(),
                    bucket_size: self.ms1_peaks.bucket_size,
                })
            },
            // Write all MS2 groups in parallel
            || {
                self.ms2_window_groups
                    .par_iter()
                    .enumerate()
                    .map(|(i, (quad, group))| {
                        let filename = format!("group_{}.parquet", i);
                        let path = format!("ms2/{}", filename);

                        let bytes = write_peaks_to_parquet_bytes(&group.peaks, config)?;
                        storage.write_bytes(&path, bytes)?;

                        Ok(Ms2GroupMetadata {
                            id: i,
                            quadrupole_isolation: quad.clone(),
                            group_info: PeakGroupMetadata {
                                relative_path: PathBuf::from("ms2").join(filename),
                                cycle_to_rt_ms: group.cycle_to_rt_ms.clone(),
                                bucket_size: group.bucket_size,
                            },
                        })
                    })
                    .collect()
            },
        );

        let ms1_meta = ms1_result?;
        let ms2_metas: Vec<_> = ms2_results.into_iter().collect::<Result<_, _>>()?;

        // Write metadata.json
        let meta = TimscentroidMetadata {
            version: SCHEMA_VERSION.to_string(),
            created_at: chrono::Utc::now().to_rfc3339(),
            ms1_peaks: ms1_meta,
            ms2_window_groups: ms2_metas,
        };

        let metadata_json = serde_json::to_string_pretty(&meta)?;
        storage.write_bytes("metadata.json", metadata_json.into_bytes())?;

        Ok(())
    }

    /// Load indexed peaks from cloud storage (async version)
    pub async fn load_from_storage_async(
        location: StorageLocation,
    ) -> Result<Self, SerializationError> {
        let storage = StorageProvider::new(location)?;

        // Read metadata.json
        let metadata_json = storage.read_to_string_async("metadata.json").await?;
        let meta: TimscentroidMetadata = serde_json::from_str(&metadata_json)?;

        // Validate schema version
        if meta.version != SCHEMA_VERSION {
            return Err(SerializationError::SchemaVersionMismatch {
                expected: SCHEMA_VERSION,
                found: meta.version,
            });
        }

        // Use storage-aware loading strategy
        Self::load_from_storage_impl(storage, meta)
    }

    /// Load indexed peaks from cloud storage (blocking version)
    pub fn load_from_storage(location: StorageLocation) -> Result<Self, SerializationError> {
        let storage = StorageProvider::new(location)?;

        // Read metadata.json
        let metadata_json = storage.read_to_string("metadata.json")?;
        let meta: TimscentroidMetadata = serde_json::from_str(&metadata_json)?;

        // Validate schema version
        if meta.version != SCHEMA_VERSION {
            return Err(SerializationError::SchemaVersionMismatch {
                expected: SCHEMA_VERSION,
                found: meta.version,
            });
        }

        // Use storage-aware loading strategy (serial for cloud, parallel for local)
        Self::load_from_storage_impl(storage, meta)
    }

    fn load_from_storage_impl(
        storage: StorageProvider,
        meta: TimscentroidMetadata,
    ) -> Result<Self, SerializationError> {
        // Choose loading strategy based on storage type:
        // - Local: CPU bottleneck → parallel loading benefits from multi-core
        // - Cloud: Network I/O bottleneck → serial loading avoids connection overhead
        let use_parallel = storage.is_local();

        if use_parallel {
            Self::load_from_storage_parallel(storage, meta)
        } else {
            Self::load_from_storage_serial(storage, meta)
        }
    }

    /// Parallel loading: optimized for local storage where CPU is the bottleneck
    fn load_from_storage_parallel(
        storage: StorageProvider,
        meta: TimscentroidMetadata,
    ) -> Result<Self, SerializationError> {
        let (ms1_result, ms2_results): (
            Result<_, SerializationError>,
            Vec<Result<_, SerializationError>>,
        ) = rayon::join(
            // Load MS1 in parallel thread
            || {
                let ms1_peaks_vec = storage
                    .read_parquet_peaks(meta.ms1_peaks.relative_path.to_str().unwrap())?;

                let (ms1_peaks, _stats) = IndexedPeakGroup::new(
                    ms1_peaks_vec,
                    meta.ms1_peaks.cycle_to_rt_ms.clone(),
                    meta.ms1_peaks.bucket_size,
                );
                Ok(ms1_peaks)
            },
            // Load all MS2 groups in parallel
            || {
                meta.ms2_window_groups
                    .par_iter()
                    .map(|group_meta| {
                        let peaks_vec = storage
                            .read_parquet_peaks(group_meta.group_info.relative_path.to_str().unwrap())?;

                        let (group, _stats) = IndexedPeakGroup::new(
                            peaks_vec,
                            group_meta.group_info.cycle_to_rt_ms.clone(),
                            group_meta.group_info.bucket_size,
                        );
                        Ok((group_meta.quadrupole_isolation.clone(), group))
                    })
                    .collect()
            },
        );

        let ms1_peaks = ms1_result?;
        let ms2_window_groups: Vec<_> = ms2_results.into_iter().collect::<Result<_, _>>()?;

        Ok(Self {
            ms1_peaks,
            ms2_window_groups,
        })
    }

    /// Serial loading: optimized for cloud storage where network I/O is the bottleneck
    ///
    /// Avoids creating multiple concurrent connections which can cause:
    /// - Rate limiting from cloud providers
    /// - Connection overhead
    /// - Inefficient use of network bandwidth
    fn load_from_storage_serial(
        storage: StorageProvider,
        meta: TimscentroidMetadata,
    ) -> Result<Self, SerializationError> {
        // Load MS1 first
        let ms1_peaks_vec = storage
            .read_parquet_peaks(meta.ms1_peaks.relative_path.to_str().unwrap())?;
        let (ms1_peaks, _stats) = IndexedPeakGroup::new(
            ms1_peaks_vec,
            meta.ms1_peaks.cycle_to_rt_ms.clone(),
            meta.ms1_peaks.bucket_size,
        );

        // Load MS2 groups sequentially
        let ms2_window_groups: Result<Vec<_>, SerializationError> = meta
            .ms2_window_groups
            .iter()
            .map(|group_meta| {
                let peaks_vec = storage
                    .read_parquet_peaks(group_meta.group_info.relative_path.to_str().unwrap())?;

                let (group, _stats) = IndexedPeakGroup::new(
                    peaks_vec,
                    group_meta.group_info.cycle_to_rt_ms.clone(),
                    group_meta.group_info.bucket_size,
                );
                Ok((group_meta.quadrupole_isolation.clone(), group))
            })
            .collect();

        Ok(Self {
            ms1_peaks,
            ms2_window_groups: ms2_window_groups?,
        })
    }

}

// Helper function to write parquet to bytes (used by cloud storage)
fn write_peaks_to_parquet_bytes<T: RTIndex>(
    peaks: &[IndexedPeak<T>],
    config: SerializationConfig,
) -> Result<Vec<u8>, SerializationError> {
    let mut buffer = Vec::new();
    let cursor = std::io::Cursor::new(&mut buffer);

    // Build writer properties with provided configuration
    let props = WriterProperties::builder()
        .set_compression(config.compression)
        .set_max_row_group_size(config.row_group_size)
        .set_statistics_enabled(EnabledStatistics::Page)
        .set_column_encoding(ColumnPath::from("mz"), Encoding::RLE)
        .set_column_encoding(ColumnPath::from("intensity"), Encoding::BYTE_STREAM_SPLIT)
        .set_column_encoding(
            ColumnPath::from("cycle_index"),
            Encoding::DELTA_BINARY_PACKED,
        )
        .build();

    let schema = Arc::new(PeakSchema::canonical());
    let mut writer = ArrowWriter::try_new(cursor, schema.clone(), Some(props))?;

    // Write peaks in chunks to control memory usage and row group size
    for chunk in peaks.chunks(config.row_group_size) {
        let mz_array = Float32Array::from_iter_values(chunk.iter().map(|p| p.mz));
        let intensity_array = Float32Array::from_iter_values(chunk.iter().map(|p| p.intensity));
        let mobility_array = Float16Array::from_iter_values(chunk.iter().map(|p| p.mobility_ook0));
        let cycle_array =
            UInt32Array::from_iter_values(chunk.iter().map(|p| p.cycle_index.as_u32()));

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(mz_array),
                Arc::new(intensity_array),
                Arc::new(mobility_array),
                Arc::new(cycle_array),
            ],
        )?;

        writer.write(&batch)?;
    }

    writer.close()?;
    Ok(buffer)
}

/// Convert a RecordBatch to peaks
///
/// This helper extracts peaks from an Arrow RecordBatch, validating
/// column types and indices. Used by both file-based and cloud-based
/// parquet readers.
pub(crate) fn batch_to_peaks<T: RTIndex>(
    batch: &RecordBatch,
) -> Result<Vec<IndexedPeak<T>>, SerializationError> {
    // Validate schema and get column indices
    let peak_schema = PeakSchema::validate(batch.schema().as_ref())?;

    // Use validated column indices - no string lookups in hot loop
    let mz = batch
        .column(peak_schema.mz_idx)
        .as_any()
        .downcast_ref::<Float32Array>()
        .ok_or_else(|| SerializationError::WrongColumnType {
            column: "mz".to_string(),
            expected: "Float32",
            got: format!("{:?}", batch.column(peak_schema.mz_idx).data_type()),
        })?;

    let intensity = batch
        .column(peak_schema.intensity_idx)
        .as_any()
        .downcast_ref::<Float32Array>()
        .ok_or_else(|| SerializationError::WrongColumnType {
            column: "intensity".to_string(),
            expected: "Float32",
            got: format!("{:?}", batch.column(peak_schema.intensity_idx).data_type()),
        })?;

    // Read as f16 directly - no conversion!
    let mobility = batch
        .column(peak_schema.mobility_idx)
        .as_any()
        .downcast_ref::<Float16Array>()
        .ok_or_else(|| SerializationError::WrongColumnType {
            column: "mobility_ook0".to_string(),
            expected: "Float16",
            got: format!("{:?}", batch.column(peak_schema.mobility_idx).data_type()),
        })?;

    let cycle = batch
        .column(peak_schema.cycle_idx)
        .as_any()
        .downcast_ref::<UInt32Array>()
        .ok_or_else(|| SerializationError::WrongColumnType {
            column: "cycle_index".to_string(),
            expected: "UInt32",
            got: format!("{:?}", batch.column(peak_schema.cycle_idx).data_type()),
        })?;

    // Pre-allocate for this batch
    let mut peaks = Vec::with_capacity(batch.num_rows());

    // Use iterator-based approach for better optimization
    // Arrow guarantees all columns have the same length in a RecordBatch
    for i in 0..batch.num_rows() {
        peaks.push(IndexedPeak {
            mz: mz.value(i),
            intensity: intensity.value(i),
            mobility_ook0: mobility.value(i), // No conversion needed!
            cycle_index: T::new(cycle.value(i)),
        });
    }

    Ok(peaks)
}


