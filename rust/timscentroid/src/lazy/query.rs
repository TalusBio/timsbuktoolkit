use std::ops::Range;
use std::sync::Arc;

use crate::storage::{
    RUNTIME,
    StorageProvider,
};
use futures::stream::StreamExt;
use object_store::path::Path as ObjectPath;
use parquet::arrow::async_reader::{
    ParquetObjectReader,
    ParquetRecordBatchStreamBuilder,
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
use half::f16;
use parquet::arrow::arrow_reader::{
    ArrowPredicate,
    RowFilter,
};

use crate::serialization::PeakSchema;
use parquet::file::metadata::ParquetMetaData;

pub struct ParquetQuerier {
    storage: StorageProvider,
    relative_path: String,
    peak_schema: PeakSchema,
    file_metadata: Arc<ParquetMetaData>,
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
    /// Initialize with storage provider and relative path, load metadata, and validate schema.
    pub fn new(
        storage: StorageProvider,
        relative_path: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let object_store = storage.as_object_store();
        let object_path = ObjectPath::from(relative_path);

        // Fetch metadata (just the footer - lightweight operation)
        let (file_metadata, peak_schema) = RUNTIME.block_on(async {
            let reader = ParquetObjectReader::new(object_store, object_path.clone());
            let builder = ParquetRecordBatchStreamBuilder::new(reader).await?;

            let schema = builder.schema();
            let peak_schema = PeakSchema::validate(schema)
                .map_err(|e| format!("Schema validation error: {}", e))?;

            let file_metadata = builder.metadata().clone();

            Ok::<_, Box<dyn std::error::Error>>((file_metadata, peak_schema))
        })?;

        Ok(Self {
            storage,
            relative_path: relative_path.to_string(),
            peak_schema,
            file_metadata,
        })
    }

    /// Query parquet file with predicates
    pub fn query(
        &self,
        mz_range: Range<f32>,
        ims_range: Option<Range<f16>>,
    ) -> Result<RecordBatch, Box<dyn std::error::Error>> {
        let object_store = self.storage.as_object_store();
        let object_path = ObjectPath::from(self.relative_path.as_str());

        RUNTIME.block_on(async {
            // Create async reader
            let reader = ParquetObjectReader::new(object_store, object_path);
            let mut builder = ParquetRecordBatchStreamBuilder::new(reader).await?;

            // Build and apply filter predicate
            let pred = FilterPredicate::new(
                mz_range,
                ims_range,
                self.peak_schema.mz_idx(),
                self.peak_schema.mobility_idx(),
                self.file_metadata.clone(),
            );
            let filter = RowFilter::new(vec![Box::new(pred)]);
            builder = builder.with_row_filter(filter);

            // Build stream and collect batches
            // Parquet will automatically:
            // 1. Use row group statistics to skip irrelevant row groups
            // 2. Fetch only needed row groups via concurrent HTTP range requests
            // 3. Apply predicates during decoding
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
        })
    }
}

// Tests are in the integration test file
