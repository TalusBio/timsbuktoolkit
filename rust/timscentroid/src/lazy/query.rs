use std::ops::Range;
use std::sync::Arc;

use crate::storage::{
    RUNTIME,
    StorageProvider,
};
use arrow::array::{
    AsArray,
    BooleanArray,
    Float16Array,
    Float32Array,
};
use arrow::compute::kernels::cmp as array_cmp;
use arrow::compute::{
    and,
    concat_batches,
};
use arrow::datatypes::{
    Float16Type,
    Float32Type,
};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use futures::stream::StreamExt;
use half::f16;
use object_store::path::Path as ObjectPath;
use parquet::arrow::arrow_reader::{
    ArrowPredicate,
    ArrowReaderMetadata,
    RowFilter,
};
use parquet::arrow::async_reader::{
    ParquetObjectReader,
    ParquetRecordBatchStreamBuilder,
};
use parquet::file::metadata::{
    ParquetMetaData,
    RowGroupMetaData,
};
use parquet::file::statistics::Statistics;

use crate::serialization::PeakSchema;

pub struct ParquetQuerier {
    storage: StorageProvider,
    relative_path: String,
    peak_schema: PeakSchema,
    file_metadata: ArrowReaderMetadata,
}

struct FilterPredicate {
    mz_range: Range<f32>,
    ims_range: Option<Range<f16>>,
    projection_mask: parquet::arrow::ProjectionMask,
}

impl FilterPredicate {
    fn new(
        mz_range: Range<f32>,
        ims_range: Option<Range<f16>>,
        index_mz: usize,
        index_ims: usize,
        file_metadata: Arc<ParquetMetaData>,
    ) -> Self {
        // Define the projection mask to only include the necessary columns
        let projection_mask = parquet::arrow::ProjectionMask::leaves(
            file_metadata.file_metadata().schema_descr(),
            vec![index_mz, index_ims],
        );

        Self {
            mz_range,
            ims_range,
            projection_mask,
        }
    }
}

impl ArrowPredicate for FilterPredicate {
    fn evaluate(&mut self, batch: RecordBatch) -> Result<BooleanArray, ArrowError> {
        // 1. Filter Column A
        // Note we are referencing by index 0 because of the projection order
        // tht we defined earlier.
        let col_mz = batch
            .column_by_name("mz")
            .unwrap()
            .as_primitive::<Float32Type>();

        // Logic: A >= start && A < end
        // Note: Arrow compute kernels handle nulls automatically (usually propagating them or treating as false)
        let a_gte = array_cmp::gt_eq(&col_mz, &Float32Array::new_scalar(self.mz_range.start))?;
        let a_lt = array_cmp::lt_eq(&col_mz, &Float32Array::new_scalar(self.mz_range.end))?;
        let mask_a = and(&a_gte, &a_lt)?;

        if let Some(ims_range) = &self.ims_range {
            // 2. Filter Column B
            let col_ims = batch
                .column_by_name("mobility_ook0")
                .unwrap()
                .as_primitive::<Float16Type>();
            // Logic: B >= start && B < end
            let b_gte = array_cmp::gt_eq(&col_ims, &Float16Array::new_scalar(ims_range.start))?;
            let b_lt = array_cmp::lt_eq(&col_ims, &Float16Array::new_scalar(ims_range.end))?;
            let mask_b = and(&b_gte, &b_lt)?;

            // 3. Combine: Mask A AND Mask B
            and(&mask_a, &mask_b)
        } else {
            Ok(mask_a)
        }
    }

    fn projection(&self) -> &parquet::arrow::ProjectionMask {
        &self.projection_mask
    }
}

impl ParquetQuerier {
    /// Initialize with storage provider and relative path, load metadata, and validate schema (async version).
    ///
    /// This async method fetches the parquet file metadata once and caches it for reuse in subsequent queries.
    /// This eliminates the need to refetch metadata on every query, significantly improving performance
    /// for cloud storage where each metadata fetch can add 50-100ms of latency.
    pub async fn new_async(
        storage: StorageProvider,
        relative_path: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let object_store = storage.as_object_store();
        let full_path = storage.build_path(relative_path);
        let object_path = ObjectPath::from(full_path.as_str());

        // Fetch metadata (just the footer - lightweight operation)
        let mut reader = ParquetObjectReader::new(object_store, object_path.clone());

        // Load metadata using ArrowReaderMetadata for caching
        let file_metadata =
            ArrowReaderMetadata::load_async(&mut reader, Default::default()).await?;

        let schema = file_metadata.schema();
        let peak_schema =
            PeakSchema::validate(schema).map_err(|e| format!("Schema validation error: {}", e))?;

        Ok(Self {
            storage,
            relative_path: relative_path.to_string(),
            peak_schema,
            file_metadata,
        })
    }

    /// Initialize with storage provider and relative path, load metadata, and validate schema (blocking version).
    ///
    /// This is a blocking wrapper around `new_async()`. For better performance in async contexts,
    /// prefer using `new_async()` directly.
    pub fn new(
        storage: StorageProvider,
        relative_path: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        RUNTIME.block_on(Self::new_async(storage, relative_path))
    }

    /// Query parquet file with predicates (async version).
    ///
    /// This async method reuses the cached metadata from `new_async()`, avoiding the overhead
    /// of refetching metadata on every query. This is critical for cloud storage performance
    /// where metadata fetches can add significant latency.
    ///
    /// The query uses predicate pushdown to:
    /// 1. Skip row groups based on statistics
    /// 2. Fetch only needed row groups via concurrent HTTP range requests
    /// 3. Apply predicates during decoding
    pub async fn query_async(
        &self,
        mz_range: Range<f32>,
        ims_range: Option<Range<f16>>,
    ) -> Result<RecordBatch, Box<dyn std::error::Error>> {
        let object_store = self.storage.as_object_store();
        let full_path = self.storage.build_path(&self.relative_path);
        let object_path = ObjectPath::from(full_path.as_str());

        // Create async reader
        let reader = ParquetObjectReader::new(object_store, object_path);

        // 1. Identify which Row Groups to read based on Statistics
        let parquet_metadata = self.file_metadata.metadata();
        let mz_col_idx = self.peak_schema.mz_idx();
        let ims_col_idx = self.peak_schema.mobility_idx();

        // Filter the row group indices
        let row_groups_to_fetch: Vec<usize> = parquet_metadata
            .row_groups()
            .iter()
            .enumerate()
            .filter(|(_idx, rg_meta)| {
                // Check MZ Range
                if !overlap_check_f32(rg_meta, mz_col_idx, &mz_range) {
                    return false;
                }

                // Check IMS Range (if exists)
                if let Some(r) = &ims_range
                    && !overlap_check_f16(rg_meta, ims_col_idx, r)
                {
                    return false;
                }

                true
            })
            .map(|(idx, _)| idx)
            .collect();

        // If no row groups match, return empty early
        if row_groups_to_fetch.is_empty() {
            let schema = self.file_metadata.schema();
            return Ok(RecordBatch::new_empty(schema.clone()));
        }

        // 2. Initialize Builder with Cached Metadata
        let mut builder =
            ParquetRecordBatchStreamBuilder::new_with_metadata(reader, self.file_metadata.clone());

        // 3. APPLY THE PRUNING (This is the missing link)
        builder = builder.with_row_groups(row_groups_to_fetch);

        // 4. Set up the RowFilter (Keep this! It filters specific rows within the kept groups)
        let pred = FilterPredicate::new(
            mz_range,
            ims_range,
            mz_col_idx,
            ims_col_idx,
            parquet_metadata.clone(),
        );
        let filter = RowFilter::new(vec![Box::new(pred)]);
        builder = builder.with_row_filter(filter);

        let mut stream = builder.build()?;
        let schema = stream.schema().clone();
        let mut batches = Vec::new();

        while let Some(batch) = stream.next().await {
            batches.push(batch?);
        }

        if batches.is_empty() {
            return Ok(RecordBatch::new_empty(schema));
        }

        let combined = concat_batches(&schema, &batches)?;
        Ok(combined)
    }

    /// Query parquet file with predicates (blocking version).
    ///
    /// This is a blocking wrapper around `query_async()`. For better performance in async contexts,
    /// prefer using `query_async()` directly.
    pub fn query(
        &self,
        mz_range: Range<f32>,
        ims_range: Option<Range<f16>>,
    ) -> Result<RecordBatch, Box<dyn std::error::Error>> {
        RUNTIME.block_on(self.query_async(mz_range, ims_range))
    }
}

// --- Helper Functions for Statistic Checking ---

fn overlap_check_f32(rg: &RowGroupMetaData, col_idx: usize, range: &Range<f32>) -> bool {
    let col_meta = rg.column(col_idx);

    if let Some(stats) = col_meta.statistics() {
        match stats {
            Statistics::Float(value_stats) => {
                match (value_stats.min_opt(), value_stats.max_opt()) {
                    (Some(min), Some(max)) => {
                        // Check for overlap: !(Max < Start || Min >= End)
                        // Simplified: Max >= Start && Min < End
                        return *max >= range.start && *min < range.end;
                    }
                    _ => return true, // Missing min/max stats? Safety fallback: read it.
                }
            }
            _ => return true, // Stats exist but wrong type? Safety fallback: read it.
        }
    }
    // If stats are missing, we must read the group to be safe
    true
}

fn overlap_check_f16(rg: &RowGroupMetaData, col_idx: usize, range: &Range<f16>) -> bool {
    let col_meta = rg.column(col_idx);

    // 1. Check if stats exist
    let stats = match col_meta.statistics() {
        Some(s) => s,
        None => return true, // No stats, must read to be safe
    };

    // 2. Match strict Physical Types
    match stats {
        // Float16 is ALWAYS stored as FixedLenByteArray in Parquet
        Statistics::FixedLenByteArray(stats) => {
            let (min_bytes, max_bytes) = match (stats.min_opt(), stats.max_opt()) {
                (Some(min), Some(max)) => (min, max),
                _ => return true, // Missing min/max stats? Safety fallback: read it.
            };

            // Safety: Float16 must be exactly 2 bytes
            if min_bytes.len() != 2 || max_bytes.len() != 2 {
                return true;
            }

            // Copy to array for conversion
            let min_arr: [u8; 2] = min_bytes.as_ref().try_into().unwrap();
            let max_arr: [u8; 2] = max_bytes.as_ref().try_into().unwrap();

            // Decode Little Endian (Parquet Standard)
            let min_val = f16::from_le_bytes(min_arr);
            let max_val = f16::from_le_bytes(max_arr);

            // Check overlap: Max >= Start && Min < End
            max_val >= range.start && min_val < range.end
        }
        // If the writer did something weird (like storing it as Int32),
        // we can't prune safely, so we default to reading the row group.
        _ => true,
    }
}

// Tests are in the integration test file
