use half::f16;
use rayon::prelude::*;
use timsrust::converters::{
    ConvertableDomain,
    Tof2MzConverter2,
};
use timsrust::readers::{
    FrameReader,
    FrameReaderError,
    TdfBlob,
};
use timsrust::{
    FrameMeta,
    FramePeaks,
    Metadata,
};

use crate::rt_mapping::{
    CycleToRTMapping,
    MS1CycleIndex,
    RTIndex,
    WindowCycleIndex,
};
use crate::utils::OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use crate::utils::{
    OptionallyRestricted,
    TupleRange,
};
pub use timsrust::TimsTofPath;

use crate::centroiding::{
    AggregatedClusteringSummary,
    CentroidingConfig,
    PeakCentroider,
};
use crate::geometry::QuadrupoleIsolationScheme;

// Context numbers:
// 15_241 - number of frames in a 22 min run
// ~ 1694 per window group
// if post-centroid our max peaks/frame is 20k
// then we can pre-allocate 1694 * 20k = 33_880_000
// peaks per window group ... or calculate it based
// on the number of frames in each group!
// at 12 bytes per peak (4 + 4 + 2 + 2) = 406 MB per window group
// Which is not horrendous tbh ...

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct IndexedPeak<T: RTIndex> {
    pub mz: f32,
    pub intensity: f32,
    pub mobility_ook0: f16,
    pub cycle_index: T,
}

/// Main struct of the whole crate.
/// Holds the pre-indexed peaks after reading a TimsTof DIA file.
/// MS2 peaks are organized by window groups.
/// MS1 peaks are stored in a single block.
///
/// The main flow is:
/// 1. Read the TimsTof file and centroid the frames. (call [IndexedTimstofPeaks::from_timstof_file])
/// 2. Query the peaks using [IndexedTimstofPeaks::query_peaks_ms1] or [IndexedTimstofPeaks::query_peaks_ms2]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndexedTimstofPeaks {
    pub(crate) ms2_window_groups: Vec<(
        QuadrupoleIsolationScheme,
        IndexedPeakGroup<WindowCycleIndex>,
    )>,
    pub(crate) ms1_peaks: IndexedPeakGroup<MS1CycleIndex>,
}

/// Statistics about the indexing process.
///
/// Usually used for logging and performance monitoring.
#[derive(Debug, Clone)]
pub struct IndexBuildingStats {
    pub ms1_total_time: std::time::Duration,
    pub ms2_total_time: std::time::Duration,
    pub ms1_stats: IndexedPeakGroupBuildingStats,
    pub ms2_stats: IndexedPeakGroupBuildingStats,
}

impl std::fmt::Display for IndexBuildingStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Index Building Stats:")?;
        writeln!(f, "MS1 Indexing Time: {:.2?}", self.ms1_total_time)?;
        writeln!(f, "MS1 Indexing Stats: {}", self.ms1_stats)?;
        writeln!(f, "MS2 Indexing Time: {:.2?}", self.ms2_total_time)?;
        writeln!(f, "MS2 Indexing Stats: {}", self.ms2_stats)?;
        Ok(())
    }
}

impl IndexedTimstofPeaks {
    /// Create an IndexedTimstofPeaks from pre-built components.
    ///
    /// This is useful for testing - you can build the MS1 and MS2 groups
    /// using `IndexedPeakGroup::new()` with mock data, then combine them here.
    ///
    /// # Example
    /// ```ignore
    /// let ms1_peaks = vec![/* ... */];
    /// let (ms1_group, _) = IndexedPeakGroup::new(ms1_peaks, cycle_to_rt, 4096);
    ///
    /// let ms2_groups = vec![
    ///     (quad_geometry1, ms2_group1),
    ///     (quad_geometry2, ms2_group2),
    /// ];
    ///
    /// let index = IndexedTimstofPeaks::from_parts(ms1_group, ms2_groups);
    /// ```
    pub fn from_parts(
        ms1_peaks: IndexedPeakGroup<MS1CycleIndex>,
        ms2_window_groups: Vec<(
            QuadrupoleIsolationScheme,
            IndexedPeakGroup<WindowCycleIndex>,
        )>,
    ) -> Self {
        Self {
            ms1_peaks,
            ms2_window_groups,
        }
    }

    pub fn from_timstof_file(
        file: &TimsTofPath,
        centroiding_config: CentroidingConfig,
    ) -> (Self, IndexBuildingStats) {
        let frame_reader = file.load_frame_reader().unwrap();
        let metadata = file.load_metadata().unwrap();

        // Read MS1 peaks
        let st = std::time::Instant::now();
        let (ms1_peaks, ms1_summ) = Self::read_ms1(&frame_reader, &metadata, centroiding_config);
        let read_time_ms1 = st.elapsed();

        // Read MS2 peaks organized by window groups
        let st = std::time::Instant::now();
        let (ms2_window_groups, ms2_summ) =
            Self::read_ms2_window_groups(&frame_reader, &metadata, centroiding_config).unwrap();
        let read_time_ms2 = st.elapsed();

        let out = Self {
            ms2_window_groups,
            ms1_peaks,
        };
        // out.print_glimpse();
        (
            out,
            IndexBuildingStats {
                ms1_total_time: read_time_ms1,
                ms2_total_time: read_time_ms2,
                ms1_stats: ms1_summ,
                ms2_stats: ms2_summ,
            },
        )
    }

    pub fn fragmented_range(&self) -> TupleRange<f64> {
        let mut min_mz = f64::MAX;
        let mut max_mz = f64::MIN;
        for (wg, _peaks) in self.ms2_window_groups.iter() {
            let (start, end) = wg
                .fragmented_range()
                .expect("Window group should have a valid range");
            if start < min_mz {
                min_mz = start;
            }
            if end > max_mz {
                max_mz = end;
            }
        }
        if min_mz <= max_mz {
            TupleRange::try_new(min_mz, max_mz).unwrap()
        } else {
            panic!("No window groups found when calculating fragmented range");
        }
    }

    /// Print a summary of the indexed peaks.
    /// Meant to be a human-readable summary.
    pub fn print_glimpse(&self) {
        println!("Indexed {} MS2 window groups", self.ms2_window_groups.len());
        println!("MS1 Peaks:");
        self.ms1_peaks.print_glimpse();
        for (i, (_wg, peaks)) in self.ms2_window_groups.iter().enumerate() {
            println!("Window Group {}: ", i);
            peaks.print_glimpse();
        }
    }

    /// Query MS1 peaks based on m/z, rt, and im ranges.
    /// All ranges are inclusive.
    pub fn query_peaks_ms1(
        &self,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<MS1CycleIndex>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = &IndexedPeak<MS1CycleIndex>> {
        self.ms1_peaks.query_peaks(mz_range, cycle_range, im_range)
    }

    /// Query MS2 peaks based on precursor m/z range, fragment m/z range, rt range, and im range.
    /// All ranges are inclusive.
    ///
    /// In essence this is equivalent to sequentially:
    /// 1. Filtering the window groups that intersect with the precursor m/z range
    ///    (and im range if provided) (call [IndexedTimstofPeaks::filter_precursor_ranges])
    /// 2. For each matching window group:
    ///     1. Narrow down the im range for each window group to the intersection of the
    ///        provided im range and the window group's im range (if im range is provided)
    ///        (call [QuadrupoleIsolationScheme::intersects_ranges])
    ///     2. query the peaks using [IndexedPeakGroup::query_peaks] with the fragment
    ///        m/z range, rt range, and the narrowed im range (if provided)
    pub fn query_peaks_ms2(
        &self,
        precursor_range_mz: TupleRange<f32>,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<WindowCycleIndex>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<
        Item = (
            &QuadrupoleIsolationScheme,
            impl Iterator<Item = &IndexedPeak<WindowCycleIndex>>,
        ),
    > {
        self.filter_precursor_ranges(precursor_range_mz, im_range)
            .map(move |(wg_info, peak_group)| {
                let local_im_range = im_range.map(|r| {
                    wg_info
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

                (
                    wg_info,
                    peak_group.query_peaks(mz_range, cycle_range, local_im_range),
                )
            })
    }

    pub fn filter_precursor_ranges(
        &self,
        precursor_range_mz: TupleRange<f32>,
        ion_mobility_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<
        Item = &(
            QuadrupoleIsolationScheme,
            IndexedPeakGroup<WindowCycleIndex>,
        ),
    > {
        let f64_mz_range = (
            precursor_range_mz.start() as f64,
            precursor_range_mz.end() as f64,
        );
        let f64_im_range = ion_mobility_range.map(|r| (r.start().to_f64(), r.end().to_f64()));

        self.ms2_window_groups.iter().filter(move |(wg, _peaks)| {
            wg.intersects(f64_mz_range, f64_im_range.unwrap_or((f64::MIN, f64::MAX)))
        })
    }

    pub fn ms1_cycle_mapping(&self) -> &CycleToRTMapping<MS1CycleIndex> {
        &self.ms1_peaks.cycle_to_rt_ms
    }

    pub fn rt_ms_to_cycle_index(&self, rt_ms: u32) -> MS1CycleIndex {
        self.ms1_peaks.cycle_to_rt_ms.ms_to_closest_index(rt_ms)
    }

    /// Read MS1 frames and return an IndexedPeakGroup along with building stats.
    ///
    /// To use this outside of the indexing use []
    fn read_ms1(
        frame_reader: &FrameReader,
        metadata: &Metadata,
        centroiding_config: CentroidingConfig,
    ) -> (
        IndexedPeakGroup<MS1CycleIndex>,
        IndexedPeakGroupBuildingStats,
    ) {
        IndexedPeakGroup::read_with_filter(
            frame_reader,
            metadata,
            |meta| matches!(meta.ms_level, timsrust::MSLevel::MS1),
            centroiding_config,
        )
    }

    fn read_ms2_window_groups(
        frame_reader: &FrameReader,
        metadata: &Metadata,
        centroiding_config: CentroidingConfig,
    ) -> Result<
        (
            Vec<(
                QuadrupoleIsolationScheme,
                IndexedPeakGroup<WindowCycleIndex>,
            )>,
            IndexedPeakGroupBuildingStats,
        ),
        (),
    > {
        let windows = frame_reader.dia_windows.as_ref().ok_or(())?;
        let mut out = Vec::with_capacity(windows.len());
        let mut out_stats: Option<IndexedPeakGroupBuildingStats> = None;

        windows.iter().for_each(|quad_query| {
            let filter = |meta: &FrameMeta| match (meta.ms_level, &meta.window_group) {
                (timsrust::MSLevel::MS1, _) => false,
                (timsrust::MSLevel::MS2, Some(wg)) => &wg.quadrupole_settings == quad_query,
                (timsrust::MSLevel::MS2, None) => unreachable!(),
                (timsrust::MSLevel::Unknown, _) => unreachable!(),
            };
            let (indexed_peaks, building_stats) = IndexedPeakGroup::read_with_filter(
                frame_reader,
                metadata,
                filter,
                centroiding_config,
            );

            let (_mz_calibrations, ims_calibrations) = metadata.get_calibration().unwrap();
            let out_wg = QuadrupoleIsolationScheme::from_quad(quad_query, |x| {
                ims_calibrations
                    .get_by_id(1)
                    .unwrap()
                    .get_conversion_function()(x)
            });

            if let Some(stats) = out_stats.as_mut() {
                *stats = stats.clone().combine(building_stats);
            } else {
                out_stats = Some(building_stats);
            }
            out.push((out_wg, indexed_peaks))
        });

        Ok((out, out_stats.unwrap()))
    }
}

/// Represents a group of indexed peaks, organized into buckets based on m/z ranges.
/// Each bucket internally sorted by retention time (rt).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndexedPeakGroup<T: RTIndex> {
    pub(crate) peaks: Vec<IndexedPeak<T>>,
    pub(crate) bucket_mz_ranges: Vec<TupleRange<f32>>,
    pub(crate) bucket_size: usize,
    pub(crate) cycle_to_rt_ms: CycleToRTMapping<T>,
}

#[derive(Debug, Clone, Copy)]
struct PeakBucket<'a, T: RTIndex> {
    inner: &'a [IndexedPeak<T>],
}

impl<'a, T: RTIndex> PeakBucket<'a, T> {
    // Returns a range such that all peaks within self.inner[range]
    // have cycle_index within the provided cycle_range.
    fn find_cycle_range(&'a self, cycle_range: TupleRange<T>) -> std::ops::Range<usize> {
        let start_idx = self
            .inner
            .partition_point(|x| x.cycle_index < cycle_range.start());
        let end_idx = start_idx
            + self.inner[start_idx..].partition_point(|x| x.cycle_index <= cycle_range.end());
        start_idx..end_idx
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

impl<'a, T: RTIndex> From<&'a [IndexedPeak<T>]> for PeakBucket<'a, T> {
    fn from(value: &'a [IndexedPeak<T>]) -> Self {
        assert!(value.first().unwrap().cycle_index <= value.last().unwrap().cycle_index);
        Self { inner: value }
    }
}

/// Statistics about an IndexedPeakGroup.
#[derive(Debug, Clone)]
pub struct IndexedPeakGroupStats {
    pub num_peaks: usize,
    pub num_buckets: usize,
    pub memory_usage_bytes: usize,
    pub sorting_time: std::time::Duration,
    pub bucketing_time: std::time::Duration,
}

impl IndexedPeakGroupStats {
    fn combine(&self, other: &Self) -> Self {
        Self {
            num_peaks: self.num_peaks + other.num_peaks,
            num_buckets: self.num_buckets + other.num_buckets,
            memory_usage_bytes: self.memory_usage_bytes + other.memory_usage_bytes,
            sorting_time: self.sorting_time + other.sorting_time,
            bucketing_time: self.bucketing_time + other.bucketing_time,
        }
    }
}

/// Statistics about the building of an IndexedPeakGroup.
#[derive(Debug, Clone)]
pub struct IndexedPeakGroupBuildingStats {
    pub indexing_stats: IndexedPeakGroupStats,
    pub clustering_stats: Option<AggregatedClusteringSummary>,
    pub read_time: std::time::Duration,
}

impl std::fmt::Display for IndexedPeakGroupBuildingStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Indexed Peak Group Building Stats:")?;
        writeln!(f, "Read Time (wall): {:.2?}", self.read_time)?;
        writeln!(f, "Indexing Stats: {:#?}", self.indexing_stats)?;
        match &self.clustering_stats {
            None => writeln!(f, "No Clustering Stats")?,
            Some(s) => {
                let avg_read_time = self.read_time / s.frames_processed as u32;
                writeln!(
                    f,
                    "Average Read Time (wall) per Frame: {:#.2?}",
                    avg_read_time
                )?;
                writeln!(f, "Clustering Stats: {}", s)?;
            }
        }
        Ok(())
    }
}

impl IndexedPeakGroupBuildingStats {
    fn combine(mut self, other: Self) -> Self {
        self.indexing_stats = self.indexing_stats.combine(&other.indexing_stats);
        self.clustering_stats = match (self.clustering_stats, other.clustering_stats) {
            (None, None) => None,
            (Some(s), None) => Some(s),
            (None, Some(s)) => Some(s),
            (Some(s1), Some(s2)) => Some(s1.combine(&s2)),
        };
        self.read_time += other.read_time;
        self
    }
}

impl<T: RTIndex> IndexedPeakGroup<T> {
    /// Query peaks based on m/z, rt, and im ranges.
    pub fn query_peaks(
        &self,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<T>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = &IndexedPeak<T>> {
        QueryPeaksIterator::new(self, mz_range, cycle_range, im_range)
    }

    /// Create a new IndexedPeakGroup from a vector of peaks.
    ///
    /// This is the **canonical** way to build an IndexedPeakGroup - all bucketing
    /// logic lives here to ensure consistency.
    ///
    /// NOTE: This internally uses `par_sort_unstable` to sort the peaks
    /// so in theory it should not be called within a parallel loop.
    pub(crate) fn new(
        mut peaks: Vec<IndexedPeak<T>>,
        cycle_to_rt_ms: CycleToRTMapping<T>,
        bucket_size: usize,
    ) -> (Self, IndexedPeakGroupStats) {
        let st = std::time::Instant::now();
        peaks.par_sort_unstable_by(|x, y| x.mz.partial_cmp(&y.mz).unwrap());
        let sort_time = st.elapsed();
        assert!(peaks.first().unwrap().mz <= peaks.last().unwrap().mz);
        assert!(peaks.iter().all(|x| x.intensity >= 0.0));
        let max_cycle = T::new(cycle_to_rt_ms.len() as u32 - 1);
        assert!(peaks.iter().all(|x| x.cycle_index <= max_cycle));

        let st = std::time::Instant::now();
        let bucket_mz_ranges: Vec<_> = peaks
            .par_chunks_mut(bucket_size)
            .map(|chunk| {
                let start = chunk.first().unwrap().mz;
                let end = chunk.last().unwrap().mz;
                chunk.sort_unstable_by(|x, y| {
                    x.cycle_index
                        .partial_cmp(&y.cycle_index)
                        .unwrap()
                        .then(x.mobility_ook0.partial_cmp(&y.mobility_ook0).unwrap())
                });
                TupleRange::try_new(start, end).expect("Incoming vec should have been sorted")
            })
            .collect();

        let bucket_time = st.elapsed();

        let tmp = Self {
            peaks,
            bucket_size,
            bucket_mz_ranges,
            cycle_to_rt_ms,
        };
        let stats = IndexedPeakGroupStats {
            num_peaks: tmp.peaks.len(),
            num_buckets: tmp.bucket_mz_ranges.len(),
            memory_usage_bytes: tmp.aproximate_memory_usage(),
            sorting_time: sort_time,
            bucketing_time: bucket_time,
        };
        (tmp, stats)
    }

    pub fn unpack(self) -> (Vec<IndexedPeak<T>>, Vec<u32>, usize, Vec<TupleRange<f32>>) {
        let Self {
            peaks,
            bucket_mz_ranges,
            bucket_size,
            cycle_to_rt_ms,
        } = self;
        (
            peaks,
            cycle_to_rt_ms.unpack(),
            bucket_size,
            bucket_mz_ranges,
        )
    }

    /// Create a new IndexedPeakGroup for testing purposes.
    ///
    /// I only make this public for testing - in real use cases
    /// You should never be used
    #[doc(hidden)]
    pub fn testing_new(
        peaks: Vec<IndexedPeak<T>>,
        cycle_to_rt_ms: CycleToRTMapping<T>,
        bucket_size: usize,
    ) -> (Self, IndexedPeakGroupStats) {
        if cfg!(test) {
            panic!("Not intended to run in production code")
        }
        Self::new(peaks, cycle_to_rt_ms, bucket_size)
    }

    fn aproximate_memory_usage(&self) -> usize {
        let self_mem = std::mem::size_of::<Self>();
        let bucket_mem =
            std::mem::size_of::<TupleRange<f32>>() * (self.bucket_mz_ranges.capacity());
        let peak_mem = self.peaks.capacity() * std::mem::size_of::<IndexedPeak<T>>();
        self_mem + bucket_mem + peak_mem
    }

    /// Return a vec of tuples, where each tuple is (frame_index, cycle_index, rt_ms)
    ///
    /// We define each cycle as an occurrence of an MS1 frame.
    /// So if we have              [MS2, MS1, MS2, MS2, MS1, MS2]
    /// the cycle indices would be [  0,   1,   1,   1,   2,   2]
    /// TODO: consider if we want to assign the index to the "closest" MS1 frame
    ///      instead of the previous one.
    fn get_frame_indices_matching(
        frame_reader: &FrameReader,
        filter: impl Fn(&FrameMeta) -> bool + Sync,
    ) -> Vec<(usize, u32, u32)> {
        let mut out = Vec::new();
        let mut cycle_index = 0;

        for (i, meta) in frame_reader.frame_metas.iter().enumerate() {
            if matches!(meta.ms_level, timsrust::MSLevel::MS1) {
                cycle_index += 1;
            }
            if filter(meta) {
                let rt_ms = (meta.rt_in_seconds * 1000.0).round() as u32;
                out.push((i, cycle_index, rt_ms));
            }
        }
        // TODO: check that cycle indices are continuous AND without replicates
        // or gaps ...
        out.as_slice().windows(2).for_each(|w| {
            // Note: These are assertions instead of errors because
            // they are invariants I am checking, not recoverable errors, none of the logic
            // after this point would make sense if this is not true.
            assert!(
                w[0].1 < w[1].1,
                "Cycle indices should be strictly increasing"
            );
            assert!(w[0].2 <= w[1].2, "Retention times should be non-decreasing");
            assert!(
                w[1].1 - w[0].1 == 1,
                "Cycle indices should be continuous without gaps, found gap between {} and {}",
                w[0].1,
                w[1].1,
            );
        });

        out
    }

    /// Read frames from a FrameReader that match a given filter function.
    #[tracing::instrument(level = "debug", skip_all)]
    fn read_with_filter(
        frame_reader: &FrameReader,
        metadata: &Metadata,
        filter: impl Fn(&FrameMeta) -> bool + Sync,
        centroiding_config: CentroidingConfig,
    ) -> (Self, IndexedPeakGroupBuildingStats) {
        // I dont like this allocation but its not that big of a deal RN ...
        let indices = Self::get_frame_indices_matching(frame_reader, filter);

        // I am still not sure what I prefer here ...
        // On the greater schme of things allocating a vec for every frame is not thaaaat bad
        // And according to benchmarks its actually faster than the lock + trim overhead
        //
        // use std::sync::{Arc, Mutex};
        // let total_peaks_estimate = indices.len() * 20_000; // assuming max 20k peaks per frame post-centroid
        // let all_peaks = Arc::new(Mutex::new(Vec::with_capacity(total_peaks_estimate)));

        let (mz_calibrations, ims_calibrations) = metadata.get_calibration().unwrap();
        let st = std::time::Instant::now();
        let _x = indices
            .par_iter()
            .with_min_len(200)
            .map_init(
                || {
                    (
                        PeakCentroider::with_capacity(
                            500_000,
                            centroiding_config,
                            metadata.mz_converter,
                            metadata.im_converter,
                        ),
                        TdfBlob::with_capacity(500_000),
                        timsrust::Frame {
                            meta: Default::default(),
                            peaks: FramePeaks::with_capacity(1000, 500_000),
                        },
                    )
                },
                |(centroider, blob_buffer, frame_buffer), idx_tpl| {
                    let (idx, cycle_index, _rt_ms) = *idx_tpl;
                    match frame_reader.get_buffered(idx, frame_buffer, blob_buffer) {
                        Ok(_) => {}
                        Err(FrameReaderError::CorruptFrame) => {
                            eprintln!("Corrupt frame found at index {}", idx);
                            return Err(FrameReaderError::CorruptFrame);
                        }
                        Err(e) => panic!("Unhandled Error reading frame at index {}: {:?}", idx, e),
                    };

                    let calibration = mz_calibrations
                        .get_by_id(frame_buffer.meta.calibration.calibration_id)
                        .unwrap();
                    let mz_converter = Tof2MzConverter2::try_from_calibration(
                        calibration,
                        frame_buffer.meta.calibration.t1,
                        frame_buffer.meta.calibration.t2,
                    )
                    .unwrap(); // TODO: make this an error instead of an option...
                    let im_converter = ims_calibrations
                        .get_by_id(frame_buffer.meta.calibration.calibration_id)
                        .unwrap()
                        .get_conversion_function();

                    let (reason, peaks_iter) = centroider.centroid_frame(frame_buffer);
                    let tmp = peaks_iter.map(|peak| {
                        let mz = mz_converter.convert(peak.tof_index as f64) as f32;
                        let mobility_ook0 = f16::from_f64(im_converter(peak.scan_index as f64));
                        let intensity = peak.corrected_intensity as f32;
                        IndexedPeak {
                            mz,
                            intensity,
                            mobility_ook0,
                            cycle_index: T::new(cycle_index),
                        }
                    });
                    // let res_arc = all_peaks.clone();
                    // let mut res = res_arc.lock().unwrap();
                    // res.extend(tmp);
                    Ok((tmp.collect::<Vec<IndexedPeak<T>>>(), reason))
                },
            )
            .map(|x| match x {
                Ok((x, summ)) => (1, 0, x, AggregatedClusteringSummary::from(summ)),
                Err(_) => (0, 1, Vec::new(), AggregatedClusteringSummary::new()),
            })
            .reduce(
                || (0, 0, Vec::new(), AggregatedClusteringSummary::new()),
                |mut a, mut b| {
                    a.2.append(&mut b.2);
                    let out_summ = a.3.combine(&b.3);
                    (a.0 + b.0, a.1 + b.1, a.2, out_summ)
                },
            );

        let read_time = st.elapsed();

        // let unwrapped_mutex = Arc::try_unwrap(all_peaks)
        //     .map_err(|_x| Err::<(), ()>(()))
        //     .unwrap();
        // let mut inner_vec = unwrapped_mutex.into_inner().unwrap();
        // inner_vec.shrink_to_fit();
        let mut inner_vec = _x.2;
        inner_vec.shrink_to_fit();
        let mut cycle_to_rt_ms: Vec<_> = Vec::with_capacity(indices.len());
        for (_i, cycle_idx, rt_ms) in indices {
            if cycle_to_rt_ms.len() <= cycle_idx as usize {
                cycle_to_rt_ms.resize(cycle_idx as usize + 1, 0);
            }
            cycle_to_rt_ms[cycle_idx as usize] = rt_ms;
        }
        // In theory there should be no zeros ...

        // 2**12 = 4096 peaks per bucket
        let cycle_mapping = CycleToRTMapping::new(cycle_to_rt_ms);
        let (out, stats) = IndexedPeakGroup::new(inner_vec, cycle_mapping, 2usize.pow(12));
        (
            out,
            IndexedPeakGroupBuildingStats {
                indexing_stats: stats,
                clustering_stats: Some(_x.3),
                read_time,
            },
        )
    }

    /// Query the bucket indices that overlap with the given m/z range.
    fn query_bucket_range(&self, mz_range: TupleRange<f32>) -> std::ops::Range<usize> {
        let start_idx = self
            .bucket_mz_ranges
            .partition_point(|x| x.end() <= mz_range.start());
        let end_idx = start_idx
            + self.bucket_mz_ranges[start_idx..].partition_point(|x| x.start() <= mz_range.end());
        // TODO: test since multiple buckets can have the same start/end values.
        start_idx..end_idx
    }

    /// Get a specific bucket by its index.
    /// Returns None if the index is out of bounds.
    fn get_bucket(&self, bucket_idx: usize) -> Option<PeakBucket<'_, T>> {
        self.get_bucket_range(bucket_idx)
            .map(|r| PeakBucket::from(&self.peaks[r]))
    }

    /// Get the range of peak indices for a specific bucket index.
    fn get_bucket_range(&self, bucket_idx: usize) -> Option<std::ops::Range<usize>> {
        let start = bucket_idx * self.bucket_size;
        if start >= self.peaks.len() {
            return None;
        }
        let end = ((bucket_idx + 1) * self.bucket_size).min(self.peaks.len());
        Some(start..end)
    }

    fn print_glimpse(&self) {
        let num_buckets = self.bucket_mz_ranges.len();
        let num_peaks = self.peaks.len();
        let mem_usage = num_peaks * std::mem::size_of::<IndexedPeak<T>>();
        println!("IndexedPeakGroup Glimpse:");
        println!("  Number of peaks: {}", num_peaks);
        println!("  Number of buckets: {}", num_buckets);
        println!("  Bucket size: {}", self.bucket_size);
        println!(
            "  Estimated memory usage: {:.2} MB",
            mem_usage as f64 / (1024.0 * 1024.0)
        );
    }
}

/// This whole struct is just a to solve lifetime issues that arise
/// from attempting a flat map over multiple iterators that borrow from self
/// when querying the peaks.
#[derive(Debug)]
struct QueryPeaksIterator<'a, T: RTIndex> {
    indexed_window_group: &'a IndexedPeakGroup<T>,
    mz_range: TupleRange<f32>,
    cycle_range: OptionallyRestricted<TupleRange<T>>,
    im_range: OptionallyRestricted<TupleRange<f16>>,
    bucket_idx: usize,
    bucket_end: usize,
    position_in_bucket: usize,
    end_of_current_bucket: usize,
    current_bucket: Option<PeakBucket<'a, T>>,
}

impl<'a, T: RTIndex> QueryPeaksIterator<'a, T> {
    pub fn new(
        indexed_window_group: &'a IndexedPeakGroup<T>,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<T>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> Self {
        let bucket_range = indexed_window_group.query_bucket_range(mz_range);

        Self {
            indexed_window_group,
            mz_range,
            cycle_range,
            im_range,
            bucket_idx: bucket_range.start,
            bucket_end: bucket_range.end,
            position_in_bucket: 0,
            end_of_current_bucket: 0,
            current_bucket: None,
        }
    }
}

impl<'a, T: RTIndex> QueryPeaksIterator<'a, T> {
    // Advance to the next bucket that matches the m/z range.
    // None if there are no more buckets.
    // Returns the index of the advanced bucket.
    fn advance_bucket(&mut self) -> Option<usize> {
        if self.bucket_idx >= self.bucket_end {
            return None;
        }
        let curr_bucket = self
            .indexed_window_group
            .get_bucket(self.bucket_idx)
            .expect("Bucket index should be valid since we calculated the bucket range in new()");
        self.bucket_idx += 1;
        match self.cycle_range.as_ref() {
            Restricted(cycle_range) => {
                let rt_idx_range = curr_bucket.find_cycle_range(*cycle_range);
                self.position_in_bucket = rt_idx_range.start;
                self.end_of_current_bucket = rt_idx_range.end;
            }
            Unrestricted => {
                self.position_in_bucket = 0;
                self.end_of_current_bucket = curr_bucket.len();
            }
        }
        self.current_bucket = Some(curr_bucket);
        Some(self.bucket_idx - 1)
    }

    fn next_in_current_bucket(&mut self) -> Option<&'a IndexedPeak<T>> {
        // Use a loop instead of recursion to avoid stack overflow when many
        // consecutive peaks don't match the filter criteria (common in dense spectra).
        while let Some(bucket) = self.current_bucket.as_ref() {
            if self.position_in_bucket >= self.end_of_current_bucket {
                return None;
            }

            let peak = &bucket.inner[self.position_in_bucket];
            self.position_in_bucket += 1;

            if self.mz_range.contains(peak.mz)
                && self
                    .im_range
                    .as_ref()
                    .is_unrestricted_or(|r| r.contains(peak.mobility_ook0))
            {
                return Some(peak);
            }
            // Continue loop to check next peak
        }
        None
    }
}

impl<'a, T: RTIndex> Iterator for QueryPeaksIterator<'a, T> {
    type Item = &'a IndexedPeak<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let peak = self.next_in_current_bucket();
            if peak.is_some() {
                return peak;
            }
            match self.advance_bucket() {
                Some(_) => continue,
                None => return None,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use half::f16;

    fn tuples_to_peaks<T: RTIndex>(data: &[(f32, f32, f32, u32)]) -> Vec<IndexedPeak<T>> {
        data.iter()
            .map(|&(mz, intensity, im, cycle_index)| IndexedPeak::<T> {
                mz,
                intensity,
                mobility_ook0: f16::from_f32(im),
                cycle_index: T::new(cycle_index),
            })
            .collect()
    }

    #[test]
    fn test_peak_bucket_query_rt() {
        let test_data = vec![
            (100.0, 200.0, 1.0, 0u32),
            (100.0, 250.0, 1.0, 1u32),
            (100.0, 250.0, 1.0, 1u32),
            (100.0, 250.0, 1.0, 1u32),
            (100.0, 300.0, 1.0, 2u32),
            (100.0, 350.0, 1.0, 3u32),
            (100.0, 400.0, 1.0, 4u32),
        ];

        let peaks = tuples_to_peaks::<MS1CycleIndex>(&test_data);
        let bucket = PeakBucket::<MS1CycleIndex>::from(&peaks[..]);
        // let rt_range: std::ops::Range<u32> = 1900..4100;
        let cycle_range = (MS1CycleIndex::new(1), MS1CycleIndex::new(3));
        let rt_idx_range = bucket.find_cycle_range(cycle_range.try_into().unwrap());
        let out: Vec<_> = bucket.inner[rt_idx_range].to_vec();
        let expected = tuples_to_peaks::<MS1CycleIndex>(&[
            (100.0, 250.0, 1.0, 1),
            (100.0, 250.0, 1.0, 1),
            (100.0, 250.0, 1.0, 1),
            (100.0, 300.0, 1.0, 2),
            (100.0, 350.0, 1.0, 3),
        ]);

        assert_eq!(out, expected);
    }
}
