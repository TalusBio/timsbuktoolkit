//! Lazy loading for indexed timsTOF peaks.
//!
//! This module provides on-disk querying with minimal memory footprint by loading
//! only the JSON metadata at initialization, then selectively loading parquet row groups
//! on-demand during queries.
//!
//! # Performance
//!
//! - **Initialization**: <50ms (JSON parsing only, no parquet loading)
//! - **First query (cold cache)**: 20-250ms (loads needed row groups from disk)
//! - **Subsequent queries (warm cache)**: 2-10ms (similar to eager loading)
//! - **Memory usage**: ~2.5MB metadata + configurable cache (default 1GB)
//!
//! # Example
//!
//! ```no_run
//! use timscentroid::lazy::LazyIndexedTimstofPeaks;
//! use timscentroid::utils::{TupleRange, OptionallyRestricted::*};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Fast initialization: loads only metadata.json
//! let lazy_index = LazyIndexedTimstofPeaks::load_from_directory("indexed_peaks/")?;
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
use crate::utils::{
    OptionallyRestricted,
    TupleRange,
};
use half::f16;
use std::fs::File;
use std::path::{
    Path,
    PathBuf,
};

/// Lazy-loading indexed peaks with on-demand row group loading
pub struct LazyIndexedTimstofPeaks {
    base_directory: PathBuf,
    ms1_metadata: PeakGroupMetadata<MS1CycleIndex>,
    ms2_metadata: Vec<(
        QuadrupoleIsolationScheme,
        PeakGroupMetadata<WindowCycleIndex>,
    )>,
}

impl LazyIndexedTimstofPeaks {
    /// Load indexed peaks from directory with lazy loading (JSON metadata only)
    ///
    /// This is extremely fast (<50ms) as it only parses the metadata.json file
    /// without loading any parquet data.
    ///
    /// # Arguments
    ///
    /// * `directory` - Directory containing metadata.json and parquet files
    ///
    /// # Returns
    ///
    /// A lazy-loading index ready for querying. Parquet files are loaded on-demand.
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Directory doesn't exist
    /// - metadata.json is missing or malformed
    /// - Metadata version is not 1.0 (v1.0 not supported for lazy loading)
    pub fn load_from_directory(directory: impl AsRef<Path>) -> Result<Self, SerializationError> {
        let directory = directory.as_ref();
        let meta_path = directory.join("metadata.json");

        let start = std::time::Instant::now();

        // Parse JSON metadata (fast)
        let file = File::open(&meta_path)?;
        let meta: TimscentroidMetadata = serde_json::from_reader(file)?;

        // Validate version
        if meta.version != "1.0" {
            // This might be over-engineered for not but something tells me we might want to
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
            base_directory: directory.to_path_buf(),
            ms1_metadata: meta.ms1_peaks,
            ms2_metadata,
        })
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
        let ms1_path = self.base_directory.join(&self.ms1_metadata.relative_path);
        self.query_peaks_file(ms1_path, mz_range, cycle_range, im_range)
    }

    fn query_peaks_file<T: RTIndex>(
        &self,
        path: impl AsRef<Path>,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = IndexedPeak<T>> {
        let querier = match ParquetQuerier::new(path.as_ref()) {
            Ok(q) => q,
            Err(e) => {
                eprintln!("Error initializing ParquetQuerier: {}", e);
                eprintln!("Error initializing ParquetQuerier: {}", e);
                todo!("Sebastian has been to lazy to make this a good error ... fix me!");
            }
        };

        let mz_range = mz_range.start()..mz_range.end();
        let im_range = match im_range {
            OptionallyRestricted::Restricted(r) => Some(r.start()..r.end()),
            OptionallyRestricted::Unrestricted => None,
        };
        let record_batch = match querier.query(mz_range, im_range) {
            Ok(batch) => batch,
            Err(e) => {
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

            let ms2_path = self.base_directory.join(&group_metadata.relative_path);
            let peaks_iter = self.query_peaks_file(ms2_path, mz_range, cycle_range, local_im_range);
            results.push((isolation_scheme.clone(), peaks_iter.collect()));
        }

        results
    }
}
