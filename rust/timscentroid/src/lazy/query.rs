use arrow::array::RecordBatchReader;
use std::fs::File;
use std::ops::Range;
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;

use arrow::array::{
    Array,
    AsArray,
    BooleanArray,
    Float16Array,
    Float32Array,
    Int32Array,
};
use arrow::compute::kernels::cmp as array_cmp;
use arrow::compute::{
    and,
    concat_batches,
};
use arrow::datatypes::{
    DataType,
    Float16Type,
    Float32Type,
    Schema,
};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use half::f16;
// use parquet::arrow::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_reader::{
    ArrowPredicate,
    RowFilter,
};
use parquet::file::metadata::FileMetaData;

use crate::serialization::{
    PeakSchema,
    SerializationError,
};
use parquet::arrow::arrow_reader::{
    ArrowPredicateFn,
    ParquetRecordBatchReaderBuilder,
};

use parquet::file::metadata::ParquetMetaData;

pub struct ParquetQuerier {
    path: PathBuf,
    peak_schema: PeakSchema,
    file_metadata: Arc<ParquetMetaData>,
}

struct FilterPredicate {
    mz_range: Range<f32>,
    ims_range: Option<Range<f16>>,
    index_mz: usize,
    index_ims: usize,
    file_metadata: Arc<ParquetMetaData>,
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
            index_mz,
            index_ims,
            file_metadata,
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
    /// 1. Initialize with path, load metadata, and validate hard-coded schema.
    pub fn new(path: impl Into<PathBuf>) -> Result<Self, Box<dyn std::error::Error>> {
        let path = path.into();
        let file = File::open(&path)?;

        // Create a builder just to inspect the schema from the file metadata
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let schema = builder.schema();
        let peak_schema =
            PeakSchema::validate(&schema).map_err(|e| format!("Schema validation error: {}", e))?;

        let file_metadata = builder.metadata().clone();

        Ok(Self {
            path,
            peak_schema,
            file_metadata,
        })
    }

    /// Single exposed query method
    pub fn query(
        &self,
        mz_range: Range<f32>,
        ims_range: Option<Range<f16>>,
    ) -> Result<RecordBatch, Box<dyn std::error::Error>> {
        let file = File::open(&self.path)?;
        let mut builder = ParquetRecordBatchReaderBuilder::try_new(file)?;

        // 1. Build the predicate
        let pred = FilterPredicate::new(
            mz_range,
            ims_range,
            self.peak_schema.mz_idx(),
            self.peak_schema.mobility_idx(),
            self.file_metadata.clone(),
        );
        let filter = RowFilter::new(vec![Box::new(pred)]);

        // 2. Apply the filter to the builder
        builder = builder.with_row_filter(filter);

        // 3. Create the reader
        let reader = builder.build()?;
        let schema = reader.schema();

        // 4. Read all batches
        let batches: Vec<RecordBatch> = reader.collect::<Result<Vec<_>, _>>()?;

        if batches.is_empty() {
            // Return an empty batch with the correct schema
            return Ok(RecordBatch::new_empty(schema));
        }

        // 5. Concatenate into a single RecordBatch
        let combined_batch = concat_batches(&schema, &batches)?;

        Ok(combined_batch)
    }
}

// --- Example Usage Test ---
#[cfg(test)]
mod tests {
    use super::*;
    // This test assumes a file exists. It's mostly to show calling syntax.
    #[test]
    fn test_signature() {
        let path = PathBuf::from("dummy.parquet");
        // We won't actually run this as the file doesn't exist
        if path.exists() {
            let querier = ParquetQuerier::new(path).unwrap();

            let range_a = 0.0f32..100.0f32;
            let range_b = f16::from_f32(0.0)..f16::from_f32(1.0);

            let result = querier.query(range_a, Some(range_b));
            match result {
                Ok(batch) => println!("Rows: {}", batch.num_rows()),
                Err(e) => println!("Error: {}", e),
            }
        }
    }
}
