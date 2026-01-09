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
pub mod query;

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
    RUNTIME,
    StorageLocation,
    StorageProvider,
};
use crate::utils::{
    OptionallyRestricted,
    TupleRange,
};
use half::f16;
use std::fmt::Debug;
use std::sync::Arc;

/// Lazy-loading indexed peaks with on-demand row group loading
#[derive(Clone)]
pub struct LazyIndexedTimstofPeaks {
    storage: StorageProvider,
    ms1_metadata: PeakGroupMetadata<MS1CycleIndex>,
    ms2_metadata: Vec<(
        QuadrupoleIsolationScheme,
        PeakGroupMetadata<WindowCycleIndex>,
    )>,
    // OPTIMIZATION: Pre-initialized queriers for reuse across queries
    // This eliminates the need to create new queriers (and fetch metadata) on every query
    ms1_querier: Arc<ParquetQuerier>,
    ms2_queriers: Vec<Arc<ParquetQuerier>>,
}

impl Debug for LazyIndexedTimstofPeaks {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LazyIndexedTimstofPeaks")
            .field("storage", &self.storage)
            .field("ms1_metadata", &self.ms1_metadata)
            .field("ms2_metadata", &self.ms2_metadata)
            .finish()
    }
}

impl LazyIndexedTimstofPeaks {
    /// Enable instrumentation for this lazy index
    ///
    /// This wraps the storage provider with metrics tracking and recreates all queriers
    /// with the instrumented storage. Returns a new instance with instrumentation enabled.
    pub fn with_instrumentation(mut self, label: impl Into<String>) -> Self {
        // Wrap the storage with instrumentation
        self.storage = self.storage.with_instrumentation(label);

        // CRITICAL: Recreate all queriers with the instrumented storage
        // The old queriers hold references to the unwrapped storage and would bypass instrumentation

        // Recreate MS1 querier
        let ms1_relative_path = self.ms1_metadata.relative_path.to_str().unwrap();
        self.ms1_querier = Arc::new(
            RUNTIME
                .block_on(ParquetQuerier::new_async(
                    self.storage.clone(),
                    ms1_relative_path,
                ))
                .expect("Failed to recreate MS1 querier with instrumentation"),
        );

        // Recreate MS2 queriers
        self.ms2_queriers.clear();
        for (_, group_metadata) in &self.ms2_metadata {
            let relative_path = group_metadata.relative_path.to_str().unwrap();
            let querier = Arc::new(
                RUNTIME
                    .block_on(ParquetQuerier::new_async(
                        self.storage.clone(),
                        relative_path,
                    ))
                    .expect("Failed to recreate MS2 querier with instrumentation"),
            );
            self.ms2_queriers.push(querier);
        }

        self
    }

    /// Print storage metrics if instrumentation is enabled
    pub fn print_metrics(&self, label: &str) {
        self.storage.print_metrics(label);
    }

    /// Get storage metrics if instrumentation is enabled
    pub fn metrics(&self) -> Option<&crate::instrumentation::StorageMetrics> {
        self.storage.metrics()
    }

    /// Add fake latency to simulate network delays (for testing cloud storage performance)
    ///
    /// This must be called AFTER `with_instrumentation()`. It simulates network latency
    /// by adding a sleep before each storage operation.
    pub fn with_fake_latency(mut self, latency: std::time::Duration) -> Self {
        // Wrap the storage with fake latency
        self.storage = self.storage.with_fake_latency(latency);

        // CRITICAL: Recreate all queriers with the latency-wrapped storage
        // The old queriers hold references to the unwrapped storage and would bypass latency

        // Recreate MS1 querier
        self.ms1_querier = Arc::new(
            ParquetQuerier::new(
                self.storage.clone(),
                self.ms1_metadata.relative_path.to_str().unwrap(),
            )
            .expect("Failed to recreate MS1 querier"),
        );

        // Recreate MS2 queriers
        self.ms2_queriers = self
            .ms2_metadata
            .iter()
            .map(|(_, meta)| {
                Arc::new(
                    ParquetQuerier::new(self.storage.clone(), meta.relative_path.to_str().unwrap())
                        .expect("Failed to recreate MS2 querier"),
                )
            })
            .collect();

        self
    }
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

        let ms2_metadata: Vec<(
            QuadrupoleIsolationScheme,
            PeakGroupMetadata<WindowCycleIndex>,
        )> = meta
            .ms2_window_groups
            .into_iter()
            .map(|m| (m.quadrupole_isolation, m.group_info))
            .collect();

        // OPTIMIZATION: Create queriers during initialization for reuse
        // This eliminates the need to create queriers (and fetch metadata) on every query

        // Create MS1 querier
        let ms1_relative_path = meta.ms1_peaks.relative_path.to_str().unwrap();
        let ms1_querier = Arc::new(
            ParquetQuerier::new_async(storage.clone(), ms1_relative_path)
                .await
                .map_err(|e| {
                    SerializationError::Io(std::io::Error::other(format!(
                        "Error creating MS1 querier: {}",
                        e
                    )))
                })?,
        );

        // Create MS2 queriers for each window group
        let mut ms2_queriers = Vec::with_capacity(ms2_metadata.len());
        for (_, group_metadata) in &ms2_metadata {
            let relative_path = group_metadata.relative_path.to_str().unwrap();
            let querier = Arc::new(
                ParquetQuerier::new_async(storage.clone(), relative_path)
                    .await
                    .map_err(|e| {
                        SerializationError::Io(std::io::Error::other(format!(
                            "Error creating MS2 querier: {}",
                            e
                        )))
                    })?,
            );
            ms2_queriers.push(querier);
        }

        let elapsed = start.elapsed();
        eprintln!("Lazy loading initialization: {:?}", elapsed);

        Ok(Self {
            storage,
            ms1_metadata: meta.ms1_peaks,
            ms2_metadata,
            ms1_querier,
            ms2_queriers,
        })
    }

    /// Load from any storage location (blocking version)
    ///
    /// This blocks the current thread. If you're in an async context,
    /// prefer `load_from_storage_async`.
    pub fn load_from_storage(location: StorageLocation) -> Result<Self, SerializationError> {
        RUNTIME.block_on(Self::load_from_storage_async(location))
    }

    pub fn ms1_metadata(&self) -> &PeakGroupMetadata<MS1CycleIndex> {
        &self.ms1_metadata
    }

    pub fn ms1_cycle_mapping(&self) -> &crate::rt_mapping::CycleToRTMapping<MS1CycleIndex> {
        &self.ms1_metadata.cycle_to_rt_ms
    }

    pub fn rt_ms_to_cycle_index(&self, rt_ms: u32) -> MS1CycleIndex {
        self.ms1_metadata.cycle_to_rt_ms.ms_to_closest_index(rt_ms)
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

    /// Query peaks from a parquet file asynchronously (async helper)
    ///
    /// OPTIMIZATION: Takes a pre-created querier instead of creating one,
    /// eliminating metadata fetching on every query
    async fn query_peaks_file_with_querier_async<T: RTIndex>(
        &self,
        querier: &ParquetQuerier,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> Result<Vec<IndexedPeak<T>>, SerializationError> {
        // Use the provided querier (no metadata fetch!)

        let mz_range_range = mz_range.start()..mz_range.end();
        let im_range_opt = match im_range {
            OptionallyRestricted::Restricted(r) => Some(r.start()..r.end()),
            OptionallyRestricted::Unrestricted => None,
        };

        // Query asynchronously
        let record_batch = querier
            .query_async(mz_range_range, im_range_opt)
            .await
            .map_err(|e| {
                SerializationError::Io(std::io::Error::other(format!(
                    "Error querying Parquet: {}",
                    e
                )))
            })?;

        // Convert RecordBatch to Vec of IndexedPeak
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

        Ok(peaks)
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

    /// Query MS2 peaks with lazy loading (async version with concurrent execution)
    ///
    /// This async method queries multiple MS2 groups concurrently, significantly improving
    /// performance for cloud storage by overlapping network I/O operations.
    ///
    /// Returns a vector of (isolation_scheme, peaks) pairs for window groups that match
    /// the precursor range.
    pub async fn query_peaks_ms2_async(
        &self,
        precursor_range_mz: TupleRange<f32>,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> Result<
        Vec<(
            QuadrupoleIsolationScheme,
            Vec<IndexedPeak<WindowCycleIndex>>,
        )>,
        SerializationError,
    > {
        // Use futures::stream to handle async iteration properly
        use futures::stream::{
            FuturesUnordered,
            StreamExt,
        };

        let mut futures = FuturesUnordered::new();

        // OPTIMIZATION: Use pre-created queriers instead of creating new ones
        // Iterate through metadata and queriers together
        for ((isolation_scheme, _group_metadata), querier) in
            self.ms2_metadata.iter().zip(self.ms2_queriers.iter())
        {
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

            let isolation_scheme = isolation_scheme.clone();
            let querier = querier.clone(); // Arc clone is cheap

            // Push the future into the unordered set - this will execute concurrently
            futures.push(async move {
                let peaks = self
                    .query_peaks_file_with_querier_async::<WindowCycleIndex>(
                        &querier,
                        mz_range,
                        cycle_range,
                        local_im_range,
                    )
                    .await?;
                Ok::<_, SerializationError>((isolation_scheme, peaks))
            });
        }

        // Collect all results concurrently
        let mut results = Vec::new();
        while let Some(result) = futures.next().await {
            results.push(result?);
        }

        Ok(results)
    }

    /// Query MS2 peaks with lazy loading (blocking version)
    ///
    /// This is a blocking wrapper around `query_peaks_ms2_async()`. For better performance,
    /// prefer using `query_peaks_ms2_async()` in async contexts.
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
        RUNTIME
            .block_on(self.query_peaks_ms2_async(
                precursor_range_mz,
                mz_range,
                cycle_range,
                im_range,
            ))
            .unwrap_or_else(|e| {
                eprintln!("Error in query_peaks_ms2: {}", e);
                Vec::new()
            })
    }
}
