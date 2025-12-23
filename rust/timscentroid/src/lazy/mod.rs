//! Lazy loading for indexed timsTOF peaks.
//!
//! # Example
//!
//! ```no_run
//! use timscentroid::lazy::LazyIndexedTimstofPeaks;
//! use timscentroid::utils::{TupleRange, OptionallyRestricted::*};
//! use timscentroid::StorageLocation;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Fast initialization: loads only metadata.json
//! let location = StorageLocation::from_path("indexed_peaks/");
//! let lazy_index = LazyIndexedTimstofPeaks::load_from_storage(location)?;
//!
//! // Query peaks (loads needed row groups on first access)
//! let peaks: Vec<_> = lazy_index.query_peaks_ms1(
//!     TupleRange::try_new(400.0, 500.0).unwrap(),  // m/z range
//!     Unrestricted,                    // all retention times
//!     Unrestricted,                    // all ion mobilities
//! ).collect();
//! # Ok(())
//! # }
//! ```
//!
//! # Future Optimizations
//!
//! The current implementation uses the existing parquet file layout where row groups
//! (~2M peaks) don't align with buckets (4096 peaks). This causes some over-reading.
//!
//! Future optimization: Reorganize parquet files to align bucket boundaries with row
//! group boundaries. This would enable more precise loading without per-bucket statistics,
//! at the cost of more row groups (smaller row group size).
mod query;

use crate::geometry::QuadrupoleIsolationScheme;
use crate::indexing::IndexedPeak;
use crate::lazy::query::ParquetQuerier;
use crate::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
    WindowCycleIndex,
};
use crate::serialization::{
    PeakGroupMetadata,
    SerializationError,
    TimscentroidMetadata,
};
use crate::storage::{
    StorageLocation,
    StorageProvider,
    RUNTIME,
};
use crate::utils::{
    OptionallyRestricted,
    TupleRange,
};
use half::f16;

/// Lazy-loading indexed peaks with on-demand row group loading
pub struct LazyIndexedTimstofPeaks {
    storage: StorageProvider,
    ms1_metadata: PeakGroupMetadata<MS1CycleIndex>,
    ms2_metadata: Vec<(
        QuadrupoleIsolationScheme,
        PeakGroupMetadata<WindowCycleIndex>,
    )>,
}

impl LazyIndexedTimstofPeaks {

    /// Load from cloud storage URL (async version)
    ///
    /// Supports s3://, gs://, az:// URLs depending on enabled feature flags.
    pub async fn load_from_url_async(url: impl AsRef<str>) -> Result<Self, SerializationError> {
        let location = StorageLocation::from_url(url)?;
        Self::load_from_storage_async(location).await
    }

    /// Load from cloud storage URL (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `load_from_url_async`.
    pub fn load_from_url(url: impl AsRef<str>) -> Result<Self, SerializationError> {
        RUNTIME.block_on(Self::load_from_url_async(url))
    }

    /// Load from any storage location (async version)
    pub async fn load_from_storage_async(
        location: StorageLocation,
    ) -> Result<Self, SerializationError> {
        let start = std::time::Instant::now();

        let storage = StorageProvider::new(location)?;

        // Read metadata.json (using async version)
        let metadata_json = storage.read_to_string_async("metadata.json").await?;
        let meta: TimscentroidMetadata = serde_json::from_str(&metadata_json)?;

        // Validate version
        if meta.version != "1.0" {
            // This might be over-engineered for now but something tells me we might want to
            // support multiple versions in the future ...
            return Err(SerializationError::SchemaVersionMismatch {
                expected: "1.0",
                found: meta.version,
            });
        }

        let ms2_metadata = meta
            .ms2_window_groups
            .into_iter()
            .map(|m| (m.quadrupole_isolation, m.group_info))
            .collect();

        let elapsed = start.elapsed();
        eprintln!("Lazy loading initialization: {:?}", elapsed);

        Ok(Self {
            storage,
            ms1_metadata: meta.ms1_peaks,
            ms2_metadata,
        })
    }

    /// Load from any storage location (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `load_from_storage_async`.
    pub fn load_from_storage(location: StorageLocation) -> Result<Self, SerializationError> {
        RUNTIME.block_on(Self::load_from_storage_async(location))
    }
}

impl LazyIndexedTimstofPeaks {
    /// Query MS1 peaks with lazy loading
    ///
    /// This method performs query planning to identify relevant row groups,
    /// loads them on-demand (with caching), and returns an iterator over matching peaks.
    pub fn query_peaks_ms1(
        &self,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = IndexedPeak<MS1CycleIndex>> {
        let relative_path = self.ms1_metadata.relative_path.to_str().unwrap();
        self.query_peaks_file(relative_path, mz_range, cycle_range, im_range)
    }

    fn query_peaks_file<T: RTIndex>(
        &self,
        relative_path: &str,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = IndexedPeak<T>> {
        // Create querier with storage provider
        let querier = match ParquetQuerier::new(self.storage.clone(), relative_path) {
            Ok(q) => q,
            Err(e) => {
                eprintln!("Error initializing ParquetQuerier: {}", e);
                return vec![].into_iter();
            }
        };

        let mz_range_range = mz_range.start()..mz_range.end();
        let im_range_opt = match im_range {
            OptionallyRestricted::Restricted(r) => Some(r.start()..r.end()),
            OptionallyRestricted::Unrestricted => None,
        };
        let record_batch = match querier.query(mz_range_range, im_range_opt) {
            Ok(batch) => batch,
            Err(e) => {
                eprintln!("Error querying Parquet: {}", e);
                // return vec![].into_iter();
                eprintln!("Error querying Parquet file: {}", e);
                todo!("Sebastian has been to lazy to make this a good error ... fix me!");
            }
        };

        // Convert RecordBatch to iterator of IndexedPeak
        let peaks: Vec<_> = record_batch
            .column_by_name("mz")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow::array::Float32Array>()
            .unwrap()
            .iter()
            .zip(
                record_batch
                    .column_by_name("intensity")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::Float32Array>()
                    .unwrap()
                    .iter(),
            )
            .zip(
                record_batch
                    .column_by_name("mobility_ook0")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::Float16Array>()
                    .unwrap()
                    .iter(),
            )
            .zip(
                record_batch
                    .column_by_name("cycle_index")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<arrow::array::UInt32Array>()
                    .unwrap()
                    .iter(),
            )
            .filter_map(|(((mz, intensity), mobility), cycle)| {
                if let (Some(mz), Some(intensity), Some(mobility), Some(cycle)) =
                    (mz, intensity, mobility, cycle)
                {
                    let mobility_f16: f16 = mobility;
                    let cycle_u32: u32 = cycle;

                    // Apply cycle range filter if specified
                    if let OptionallyRestricted::Restricted(cycle_range) = &cycle_range
                        && !cycle_range.contains(cycle_u32)
                    {
                        return None;
                    }

                    Some(IndexedPeak {
                        mz,
                        intensity,
                        mobility_ook0: mobility_f16,
                        cycle_index: T::new(cycle_u32),
                    })
                } else {
                    None
                }
            })
            .collect();

        peaks.into_iter()
    }

    /// Query MS2 peaks with lazy loading
    ///
    /// Returns an iterator over (isolation_scheme, peaks_iterator) pairs for
    /// window groups that match the precursor range.
    pub fn query_peaks_ms2(
        &self,
        precursor_range_mz: TupleRange<f32>,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> Vec<(
        QuadrupoleIsolationScheme,
        Vec<IndexedPeak<WindowCycleIndex>>,
    )> {
        // TODO: most of this logic should be re-implemented as an iterator, and we can just return
        // the opaque iterator type.

        // We could expose a version of this that takes an input mutable reference to a vector
        // of peaks to avoid allocations in tight loops. (bring your own mutable state to work ...)
        let mut results = Vec::new();

        for (isolation_scheme, group_metadata) in &self.ms2_metadata {
            // Filter by precursor range using geometric intersection
            let f64_mz_range = (
                precursor_range_mz.start() as f64,
                precursor_range_mz.end() as f64,
            );
            let f64_im_range = im_range.map(|r| (r.start().to_f64(), r.end().to_f64()));

            if !isolation_scheme
                .intersects(f64_mz_range, f64_im_range.unwrap_or((f64::MIN, f64::MAX)))
            {
                continue;
            }
            let local_im_range = im_range.map(|r| {
                isolation_scheme
                    .intersects_ranges(
                        (
                            precursor_range_mz.start() as f64,
                            precursor_range_mz.end() as f64,
                        ),
                        (r.start().to_f64(), r.end().to_f64()),
                    )
                    .map(|(im_start, im_end)| {
                        TupleRange::try_new(f16::from_f64(im_start), f16::from_f64(im_end))
                            .expect("Intersect should always be valid")
                    })
                    .expect("Since we filtered before the precursor range. It should intersect")
            });

            let relative_path = group_metadata.relative_path.to_str().unwrap();
            let peaks_iter =
                self.query_peaks_file(relative_path, mz_range, cycle_range, local_im_range);
            results.push((isolation_scheme.clone(), peaks_iter.collect()));
        }

        results
    }
}
