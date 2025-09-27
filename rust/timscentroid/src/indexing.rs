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
pub struct IndexedPeak {
    pub mz: f32,
    pub intensity: f32,
    pub mobility_ook0: f16,
    pub rt_ms: u32,
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
    ms2_window_groups: Vec<(QuadrupoleIsolationScheme, IndexedPeakGroup)>,
    ms1_peaks: IndexedPeakGroup,
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
        rt_milliseconds_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = &IndexedPeak> {
        self.ms1_peaks
            .query_peaks(mz_range, rt_milliseconds_range, im_range)
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
        rt_milliseconds_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<
        Item = (
            &QuadrupoleIsolationScheme,
            impl Iterator<Item = &IndexedPeak>,
        ),
    > {
        self.filter_precursor_ranges(precursor_range_mz, im_range)
            .map(move |(wg_info, peak_group)| {
                let im_range = im_range.map(|r| {
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
                    peak_group.query_peaks(mz_range, rt_milliseconds_range, im_range),
                )
            })
    }

    pub fn filter_precursor_ranges(
        &self,
        precursor_range_mz: TupleRange<f32>,
        ion_mobility_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = &(QuadrupoleIsolationScheme, IndexedPeakGroup)> {
        let f64_mz_range = (
            precursor_range_mz.start() as f64,
            precursor_range_mz.end() as f64,
        );
        let f64_im_range = ion_mobility_range.map(|r| (r.start().to_f64(), r.end().to_f64()));
        let matching_iter = self.ms2_window_groups.iter().filter(move |(wg, _peaks)| {
            wg.intersects(f64_mz_range, f64_im_range.unwrap_or((f64::MIN, f64::MAX)))
        });

        matching_iter
    }

    /// Read MS1 frames and return an IndexedPeakGroup along with building stats.
    ///
    /// To use this outside of the indexing use []
    fn read_ms1(
        frame_reader: &FrameReader,
        metadata: &Metadata,
        centroiding_config: CentroidingConfig,
    ) -> (IndexedPeakGroup, IndexedPeakGroupBuildingStats) {
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
            Vec<(QuadrupoleIsolationScheme, IndexedPeakGroup)>,
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
pub struct IndexedPeakGroup {
    // TODO: Implement a way to represent the quad settings as polygons
    peaks: Vec<IndexedPeak>,
    bucket_mz_ranges: Vec<TupleRange<f32>>,
    bucket_size: usize,
}

#[derive(Debug, Clone, Copy)]
struct PeakBucket<'a> {
    inner: &'a [IndexedPeak],
}

impl<'a> PeakBucket<'a> {
    fn find_rt_range(&'a self, rt_miliseconds_range: TupleRange<u32>) -> std::ops::Range<usize> {
        let start_idx = self
            .inner
            .partition_point(|x| x.rt_ms < rt_miliseconds_range.start());
        let end_idx = start_idx
            + self.inner[start_idx..].partition_point(|x| x.rt_ms <= rt_miliseconds_range.end());
        start_idx..end_idx
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

impl<'a> From<&'a [IndexedPeak]> for PeakBucket<'a> {
    fn from(value: &'a [IndexedPeak]) -> Self {
        assert!(value.first().unwrap().rt_ms <= value.last().unwrap().rt_ms);
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
    pub clustering_stats: AggregatedClusteringSummary,
    pub read_time: std::time::Duration,
}

impl std::fmt::Display for IndexedPeakGroupBuildingStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Indexed Peak Group Building Stats:")?;
        writeln!(f, "Read Time (wall): {:.2?}", self.read_time)?;
        let avg_read_time = self.read_time / self.clustering_stats.frames_processed as u32;
        writeln!(
            f,
            "Average Read Time (wall) per Frame: {:#.2?}",
            avg_read_time
        )?;
        writeln!(f, "Indexing Stats: {:#?}", self.indexing_stats)?;
        writeln!(f, "Clustering Stats: {}", self.clustering_stats)?;
        Ok(())
    }
}

impl IndexedPeakGroupBuildingStats {
    fn combine(mut self, other: Self) -> Self {
        self.indexing_stats = self.indexing_stats.combine(&other.indexing_stats);
        self.clustering_stats = self.clustering_stats.combine(&other.clustering_stats);
        self.read_time += other.read_time;
        self
    }
}

impl IndexedPeakGroup {
    /// Query peaks based on m/z, rt, and im ranges.
    pub fn query_peaks(
        &self,
        mz_range: TupleRange<f32>,
        rt_milliseconds_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = &IndexedPeak> {
        QueryPeaksIterator::new(self, mz_range, rt_milliseconds_range, im_range)
    }

    /// Create a new IndexedPeakGroup from a vector of peaks.
    ///
    /// NOTE: This internally uses `par_sort_unstable` to sort the peaks
    /// so in theory it should not be called within a parallel loop.
    fn new(mut peaks: Vec<IndexedPeak>, bucket_size: usize) -> (Self, IndexedPeakGroupStats) {
        let st = std::time::Instant::now();
        peaks.par_sort_unstable_by(|x, y| {
            x.mz.partial_cmp(&y.mz)
                .unwrap()
                .then(x.rt_ms.partial_cmp(&y.rt_ms).unwrap())
        });
        let sort_time = st.elapsed();
        assert!(peaks.first().unwrap().mz <= peaks.last().unwrap().mz);

        let st = std::time::Instant::now();
        let bucket_mz_ranges: Vec<_> = peaks
            .par_chunks_mut(bucket_size)
            .map(|chunk| {
                let start = chunk.first().unwrap().mz;
                let end = chunk.last().unwrap().mz;
                chunk.sort_unstable_by(|x, y| {
                    x.rt_ms
                        .partial_cmp(&y.rt_ms)
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

    fn aproximate_memory_usage(&self) -> usize {
        let self_mem = std::mem::size_of::<Self>();
        let bucket_mem =
            std::mem::size_of::<TupleRange<f32>>() * (self.bucket_mz_ranges.capacity());
        let peak_mem = self.peaks.capacity() * std::mem::size_of::<IndexedPeak>();
        self_mem + bucket_mem + peak_mem
    }

    /// Read frames from a FrameReader that match a given filter function.
    fn read_with_filter(
        frame_reader: &FrameReader,
        metadata: &Metadata,
        filter: impl Fn(&FrameMeta) -> bool + Sync,
        centroiding_config: CentroidingConfig,
    ) -> (Self, IndexedPeakGroupBuildingStats) {
        // I dont like this allocation but its not that big of a deal RN ...
        let indices: Vec<usize> = frame_reader
            .frame_metas
            .iter()
            .enumerate()
            .filter_map(|(i, meta)| if filter(meta) { Some(i) } else { None })
            .collect();

        // I am still not sure what I prefer here ...
        // On the greater schme of things allocating a vec for every frame is not thaaaat bad
        //
        // use std::sync::{Arc, Mutex};
        // let total_peaks_estimate = indices.len() * 20_000; // assuming max 20k peaks per frame post-centroid
        // let all_peaks = Arc::new(Mutex::new(Vec::with_capacity(total_peaks_estimate)));

        let (mz_calibrations, ims_calibrations) = metadata.get_calibration().unwrap();
        let st = std::time::Instant::now();
        let _x = indices
            .into_par_iter()
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
                |(centroider, blob_buffer, frame_buffer), idx| {
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

                    let rt_milliseconds = (frame_buffer.meta.rt_in_seconds * 1000.0).floor() as u32;
                    let (reason, peaks_iter) = centroider.centroid_frame(frame_buffer);
                    let tmp = peaks_iter.map(|peak| {
                        let mz = mz_converter.convert(peak.tof_index as f64) as f32;
                        let im = f16::from_f64(im_converter(peak.scan_index as f64)); // TODO: convert to actual IM using calibration
                        let intensity = peak.corrected_intensity as f32;
                        IndexedPeak {
                            mz,
                            intensity,
                            mobility_ook0: im,
                            rt_ms: rt_milliseconds,
                        }
                    });
                    // let res_arc = all_peaks.clone();
                    // let mut res = res_arc.lock().unwrap();
                    // res.extend(tmp);
                    Ok((tmp.collect::<Vec<_>>(), reason))
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

        // 2**12 = 4096 peaks per bucket
        let (out, stats) = IndexedPeakGroup::new(inner_vec, 2usize.pow(12));
        (
            out,
            IndexedPeakGroupBuildingStats {
                indexing_stats: stats,
                clustering_stats: _x.3,
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
    fn get_bucket(&self, bucket_idx: usize) -> Option<PeakBucket<'_>> {
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
        let mem_usage = num_peaks * std::mem::size_of::<IndexedPeak>();
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
struct QueryPeaksIterator<'a> {
    indexed_window_group: &'a IndexedPeakGroup,
    mz_range: TupleRange<f32>,
    rt_milliseconds_range: OptionallyRestricted<TupleRange<u32>>,
    im_range: OptionallyRestricted<TupleRange<f16>>,
    bucket_idx: usize,
    bucket_end: usize,
    position_in_bucket: usize,
    end_of_current_bucket: usize,
    current_bucket: Option<PeakBucket<'a>>,
}

impl<'a> QueryPeaksIterator<'a> {
    pub fn new(
        indexed_window_group: &'a IndexedPeakGroup,
        mz_range: TupleRange<f32>,
        rt_milliseconds_range: OptionallyRestricted<TupleRange<u32>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> Self {
        let bucket_range = indexed_window_group.query_bucket_range(mz_range);
        Self {
            indexed_window_group,
            mz_range,
            rt_milliseconds_range,
            im_range,
            bucket_idx: bucket_range.start,
            bucket_end: bucket_range.end,
            position_in_bucket: 0,
            end_of_current_bucket: 0,
            current_bucket: None,
        }
    }
}

impl<'a> QueryPeaksIterator<'a> {
    fn advance_bucket(&mut self) -> Option<usize> {
        if self.bucket_idx >= self.bucket_end {
            return None;
        }
        let curr_bucket = self
            .indexed_window_group
            .get_bucket(self.bucket_idx)
            .expect("Bucket index should be valid since we calculated the bucket range in new()");
        self.bucket_idx += 1;
        match self.rt_milliseconds_range.as_ref() {
            Restricted(rt_seconds_range) => {
                let rt_idx_range = curr_bucket.find_rt_range(*rt_seconds_range);
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

    fn next_in_current_bucket(&mut self) -> Option<&'a IndexedPeak> {
        if let Some(bucket) = self.current_bucket.as_ref() {
            if self.position_in_bucket < self.end_of_current_bucket {
                let peak = &bucket.inner[self.position_in_bucket];
                self.position_in_bucket += 1;
                if self.mz_range.contains(peak.mz)
                    && self
                        .im_range
                        .as_ref()
                        .is_unrestricted_or(|r| r.contains(peak.mobility_ook0))
                {
                    return Some(peak);
                } else {
                    return self.next_in_current_bucket();
                }
            }
        }
        None
    }
}

impl<'a> Iterator for QueryPeaksIterator<'a> {
    type Item = &'a IndexedPeak;

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

    fn tuples_to_peaks(data: &[(f32, f32, f32, u32)]) -> Vec<IndexedPeak> {
        data.iter()
            .map(|&(mz, intensity, im, rt_milliseconds)| IndexedPeak {
                mz,
                intensity,
                mobility_ook0: f16::from_f32(im),
                rt_ms: rt_milliseconds,
            })
            .collect()
    }

    #[test]
    fn test_peak_bucket_query_rt() {
        let test_data = vec![
            (100.0, 200.0, 1.0, 1_000u32),
            (100.0, 250.0, 1.0, 2_000u32),
            (100.0, 250.0, 1.0, 2_000u32),
            (100.0, 250.0, 1.0, 2_000u32),
            (100.0, 300.0, 1.0, 3_000u32),
            (100.0, 350.0, 1.0, 4_000u32),
            (100.0, 400.0, 1.0, 5_000u32),
        ];

        let peaks = tuples_to_peaks(&test_data);
        let bucket = PeakBucket::from(&peaks[..]);
        let rt_range: std::ops::Range<u32> = 1900..4100;
        let rt_idx_range = bucket.find_rt_range(rt_range.try_into().unwrap());
        let out: Vec<_> = bucket.inner[rt_idx_range].to_vec();
        let expected = tuples_to_peaks(&[
            (100.0, 250.0, 1.0, 2_000),
            (100.0, 250.0, 1.0, 2_000),
            (100.0, 250.0, 1.0, 2_000),
            (100.0, 300.0, 1.0, 3_000),
            (100.0, 350.0, 1.0, 4_000),
        ]);

        assert_eq!(out, expected);
    }
}
