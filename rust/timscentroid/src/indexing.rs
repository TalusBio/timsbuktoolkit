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
use tracing::{
    instrument,
    warn,
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
    MobInt,
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

    /// Rebuild every peak group (MS1 + all MS2 window groups) at a new
    /// bucket size. See `IndexedPeakGroup::rebucket` for the tradeoff.
    pub fn rebucket(self, new_bucket_size: usize) -> Self {
        let ms1_peaks = self.ms1_peaks.rebucket(new_bucket_size);
        let ms2_window_groups = self
            .ms2_window_groups
            .into_iter()
            .map(|(quad, pg)| (quad, pg.rebucket(new_bucket_size)))
            .collect();
        Self {
            ms1_peaks,
            ms2_window_groups,
        }
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
    ) -> impl Iterator<Item = IndexedPeak<MS1CycleIndex>> {
        self.ms1_peaks.query_peaks(mz_range, cycle_range, im_range)
    }

    /// Callback-style MS1 peak scan (see `IndexedPeakGroup::for_each_peak`).
    #[inline]
    pub fn for_each_ms1_peak<F>(
        &self,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<MS1CycleIndex>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
        f: F,
    ) where
        F: FnMut(&IndexedPeak<MS1CycleIndex>),
    {
        self.ms1_peaks
            .for_each_peak(mz_range, cycle_range, im_range, f)
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
            impl Iterator<Item = IndexedPeak<WindowCycleIndex>>,
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

/// Owned SoA column buffers for a peak set.
///
/// Four length-aligned columns carried together as a single value. Fields
/// are private so the "all columns same length" invariant can't be broken
/// by a direct push to one column — callers go through `push`, which
/// advances all four in lockstep.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
pub struct PeakColumns<T: RTIndex> {
    mz: Vec<f32>,
    intensity: Vec<f32>,
    mobility: Vec<MobInt>,
    cycle_index: Vec<T>,
}

impl<T: RTIndex> PeakColumns<T> {
    pub fn with_capacity(n: usize) -> Self {
        Self {
            mz: Vec::with_capacity(n),
            intensity: Vec::with_capacity(n),
            mobility: Vec::with_capacity(n),
            cycle_index: Vec::with_capacity(n),
        }
    }

    pub fn len(&self) -> usize {
        self.mz.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mz.is_empty()
    }

    /// Reserve `n` more slots across all four columns.
    pub fn reserve(&mut self, n: usize) {
        self.mz.reserve(n);
        self.intensity.reserve(n);
        self.mobility.reserve(n);
        self.cycle_index.reserve(n);
    }

    /// Atomic row append — all four columns advance together.
    #[inline]
    pub fn push(&mut self, mz: f32, intensity: f32, mobility: MobInt, cycle_index: T) {
        self.mz.push(mz);
        self.intensity.push(intensity);
        self.mobility.push(mobility);
        self.cycle_index.push(cycle_index);
    }

    /// Borrow all four columns as aligned slices (whole-group view).
    pub fn view(&self) -> PeakColumnsView<'_, T> {
        PeakColumnsView {
            mz: &self.mz,
            intensity: &self.intensity,
            mobility: &self.mobility,
            cycle_index: &self.cycle_index,
        }
    }

    /// Consume into raw column Vecs. Private — drops the bundled invariant;
    /// used only by the AoS sort round-trip inside this module.
    fn into_parts(self) -> (Vec<f32>, Vec<f32>, Vec<MobInt>, Vec<T>) {
        (self.mz, self.intensity, self.mobility, self.cycle_index)
    }
}

/// Borrowed view over four aligned SoA columns. Used for whole-group access
/// and, when narrowed via `slice(...)`, for bucket-scoped access.
///
/// Any sub-view obtained by slicing a cycle-sorted range (i.e. a bucket slice
/// from the group) preserves that ordering, so `find_cycle_range` is valid on
/// bucket-scoped views.
#[derive(Debug, Clone, Copy)]
pub struct PeakColumnsView<'a, T: RTIndex> {
    mz: &'a [f32],
    intensity: &'a [f32],
    mobility: &'a [MobInt],
    cycle_index: &'a [T],
}

impl<'a, T: RTIndex> PeakColumnsView<'a, T> {
    pub fn mz(&self) -> &'a [f32] {
        self.mz
    }

    pub fn intensity(&self) -> &'a [f32] {
        self.intensity
    }

    pub fn mobility(&self) -> &'a [MobInt] {
        self.mobility
    }

    pub fn cycle_index(&self) -> &'a [T] {
        self.cycle_index
    }

    pub fn len(&self) -> usize {
        self.cycle_index.len()
    }

    pub fn is_empty(&self) -> bool {
        self.cycle_index.is_empty()
    }

    /// Materialize a single peak from column index `i`.
    ///
    /// LLVM SROA elides the temporary when the consumer reads one field.
    #[inline(always)]
    pub fn materialize(&self, i: usize) -> IndexedPeak<T>
    where
        T: Copy,
    {
        IndexedPeak {
            mz: self.mz[i],
            intensity: self.intensity[i],
            mobility_ook0: self.mobility[i].to_f16(),
            cycle_index: self.cycle_index[i],
        }
    }

    /// Narrow to a sub-range across all four columns in lockstep.
    #[inline(always)]
    pub fn slice(&self, range: std::ops::Range<usize>) -> Self {
        Self {
            mz: &self.mz[range.clone()],
            intensity: &self.intensity[range.clone()],
            mobility: &self.mobility[range.clone()],
            cycle_index: &self.cycle_index[range],
        }
    }

    /// Binary-search the cycle_index column for rows with cycle in range.
    ///
    /// Caller's responsibility: this view's `cycle_index` slice is sorted
    /// ascending. True for bucket-scoped sub-views of an `IndexedPeakGroup`
    /// (the group's bucketing invariant), undefined otherwise.
    pub fn find_cycle_range(&self, cycle_range: TupleRange<T>) -> std::ops::Range<usize> {
        let start_idx = self
            .cycle_index
            .partition_point(|x| *x < cycle_range.start());
        let end_idx =
            start_idx + self.cycle_index[start_idx..].partition_point(|x| *x <= cycle_range.end());
        start_idx..end_idx
    }

    /// Split into an iterator of fixed-N chunks plus a sub-view over the
    /// `< N` remainder. Mirrors `slice::as_chunks` across all four columns.
    #[inline(always)]
    pub fn as_chunks<const N: usize>(
        &self,
    ) -> (
        impl Iterator<Item = PeakColumnsChunk<'a, T, N>> + use<'a, T, N>,
        PeakColumnsView<'a, T>,
    )
    where
        T: Copy,
    {
        let (mz_c, mz_t) = self.mz.as_chunks::<N>();
        let (mob_c, mob_t) = self.mobility.as_chunks::<N>();
        let (int_c, int_t) = self.intensity.as_chunks::<N>();
        let (cyc_c, cyc_t) = self.cycle_index.as_chunks::<N>();
        let iter = mz_c
            .iter()
            .zip(mob_c.iter())
            .zip(int_c.iter())
            .zip(cyc_c.iter())
            .map(
                |(((mz, mobility), intensity), cycle_index)| PeakColumnsChunk {
                    mz,
                    mobility,
                    intensity,
                    cycle_index,
                },
            );
        let tail = PeakColumnsView {
            mz: mz_t,
            mobility: mob_t,
            intensity: int_t,
            cycle_index: cyc_t,
        };
        (iter, tail)
    }
}

/// Borrowed fixed-N chunk across the four SoA columns. Yielded by
/// `PeakColumnsView::as_chunks`.
#[derive(Debug, Clone, Copy)]
pub struct PeakColumnsChunk<'a, T: RTIndex, const N: usize> {
    mz: &'a [f32; N],
    mobility: &'a [MobInt; N],
    intensity: &'a [f32; N],
    cycle_index: &'a [T; N],
}

impl<'a, T: RTIndex + Copy, const N: usize> PeakColumnsChunk<'a, T, N> {
    pub fn mz(&self) -> &'a [f32; N] {
        self.mz
    }

    pub fn mobility(&self) -> &'a [MobInt; N] {
        self.mobility
    }

    /// Materialize a single peak from chunk-local index `i`.
    #[inline(always)]
    pub fn materialize(&self, i: usize) -> IndexedPeak<T> {
        IndexedPeak {
            mz: self.mz[i],
            intensity: self.intensity[i],
            mobility_ook0: self.mobility[i].to_f16(),
            cycle_index: self.cycle_index[i],
        }
    }
}

/// Represents a group of indexed peaks, organized into buckets based on m/z ranges.
/// Each bucket internally sorted by retention time (rt).
///
/// Storage is SoA via `PeakColumns<T>`. `bucket_size` + `bucket_mz_ranges`
/// describe bucket boundaries implicitly via column offsets
/// `[k*bucket_size, (k+1)*bucket_size)`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndexedPeakGroup<T: RTIndex> {
    pub(crate) cols: PeakColumns<T>,
    pub(crate) bucket_mz_ranges: Vec<TupleRange<f32>>,
    pub(crate) bucket_size: usize,
    pub(crate) cycle_to_rt_ms: CycleToRTMapping<T>,
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

// Filter-funnel counters for `for_each_peak`, gated behind the
// `query-instr` feature. Inside `for_each_peak` we accumulate into a
// `FepLocal` (stack-local u32 fields) and flush once per call via a
// handful of atomic adds. With the feature off, `FepLocal` becomes a
// zero-sized struct and every bump / flush compiles to nothing.
#[cfg(feature = "query-instr")]
use std::sync::atomic::AtomicU64;
#[cfg(feature = "query-instr")]
pub static FEP_CALLS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "query-instr")]
pub static FEP_BUCKETS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "query-instr")]
pub static FEP_BUCKETS_MZ_CONTAINED: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "query-instr")]
pub static FEP_PEAKS_IN_CYCLE_RANGE: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "query-instr")]
pub static FEP_PEAKS_AFTER_MZ: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "query-instr")]
pub static FEP_PEAKS_AFTER_IM: AtomicU64 = AtomicU64::new(0);

/// Per-call filter-funnel counters. Fields only exist when the
/// `query-instr` feature is enabled; the struct is zero-sized
/// otherwise and every method body compiles to nothing.
#[derive(Default)]
struct FepLocal {
    #[cfg(feature = "query-instr")]
    buckets: u32,
    #[cfg(feature = "query-instr")]
    buckets_mz_contained: u32,
    #[cfg(feature = "query-instr")]
    peaks_in_cycle_range: u32,
    #[cfg(feature = "query-instr")]
    peaks_after_mz: u32,
    #[cfg(feature = "query-instr")]
    peaks_after_im: u32,
}

impl FepLocal {
    #[inline(always)]
    #[cfg_attr(not(feature = "query-instr"), allow(unused_variables))]
    fn bump_bucket(&mut self, mz_contained: bool) {
        #[cfg(feature = "query-instr")]
        {
            self.buckets += 1;
            if mz_contained {
                self.buckets_mz_contained += 1;
            }
        }
    }

    #[inline(always)]
    #[cfg_attr(not(feature = "query-instr"), allow(unused_variables))]
    fn bump_in_cycle(&mut self, n: u32) {
        #[cfg(feature = "query-instr")]
        {
            self.peaks_in_cycle_range += n;
        }
    }

    #[inline(always)]
    #[cfg_attr(not(feature = "query-instr"), allow(unused_variables))]
    fn bump_after_mz(&mut self, n: u32) {
        #[cfg(feature = "query-instr")]
        {
            self.peaks_after_mz += n;
        }
    }

    #[inline(always)]
    fn bump_after_im(&mut self) {
        #[cfg(feature = "query-instr")]
        {
            self.peaks_after_im += 1;
        }
    }

    /// Count survivors in a chunk mask and add them to `peaks_after_mz`.
    /// The count loop is compiled out when `query-instr` is disabled.
    #[inline(always)]
    #[cfg_attr(not(feature = "query-instr"), allow(unused_variables))]
    fn count_after_mz_mask<const N: usize>(&mut self, mask: &[bool; N]) {
        #[cfg(feature = "query-instr")]
        {
            let mut c: u32 = 0;
            for i in 0..N {
                c += mask[i] as u32;
            }
            self.peaks_after_mz += c;
        }
    }

    /// Count survivors in a chunk mask and add them to `peaks_after_im`.
    /// The count loop is compiled out when `query-instr` is disabled.
    #[inline(always)]
    #[cfg_attr(not(feature = "query-instr"), allow(unused_variables))]
    fn count_after_im_mask<const N: usize>(&mut self, mask: &[bool; N]) {
        #[cfg(feature = "query-instr")]
        {
            let mut c: u32 = 0;
            for i in 0..N {
                c += mask[i] as u32;
            }
            self.peaks_after_im += c;
        }
    }

    /// Push local counts into the global atomics. One fetch_add per
    /// counter, once per `for_each_peak` call. Destructured so adding
    /// a new field to `FepLocal` fails to compile here unless
    /// matched + flushed.
    #[inline(always)]
    fn flush(&self) {
        #[cfg(feature = "query-instr")]
        {
            use std::sync::atomic::Ordering::Relaxed;
            let Self {
                buckets,
                buckets_mz_contained,
                peaks_in_cycle_range,
                peaks_after_mz,
                peaks_after_im,
            } = *self;
            FEP_CALLS.fetch_add(1, Relaxed);
            FEP_BUCKETS.fetch_add(buckets as u64, Relaxed);
            FEP_BUCKETS_MZ_CONTAINED.fetch_add(buckets_mz_contained as u64, Relaxed);
            FEP_PEAKS_IN_CYCLE_RANGE.fetch_add(peaks_in_cycle_range as u64, Relaxed);
            FEP_PEAKS_AFTER_MZ.fetch_add(peaks_after_mz as u64, Relaxed);
            FEP_PEAKS_AFTER_IM.fetch_add(peaks_after_im as u64, Relaxed);
        }
    }
}

#[cfg_attr(not(feature = "query-instr"), allow(unused_variables))]
pub fn dump_for_each_peak_funnel(label: &str) {
    #[cfg(feature = "query-instr")]
    {
        use std::sync::atomic::Ordering::Relaxed;
        let calls = FEP_CALLS.load(Relaxed);
        let buckets = FEP_BUCKETS.load(Relaxed);
        let buckets_contained = FEP_BUCKETS_MZ_CONTAINED.load(Relaxed);
        let cyc = FEP_PEAKS_IN_CYCLE_RANGE.load(Relaxed);
        let mz_ok = FEP_PEAKS_AFTER_MZ.load(Relaxed);
        let im_ok = FEP_PEAKS_AFTER_IM.load(Relaxed);
        let denom = cyc.max(1);
        eprintln!(
            "[fep-funnel:{label}] calls={calls} buckets={buckets} \
             buckets_mz_contained={buckets_contained} ({:.1}% of visited) \
             in_cycle_range={cyc} after_mz={mz_ok} ({:.1}% pass) \
             after_im={im_ok} ({:.1}% pass, {:.1}% of total) \
             mz_killed={} ({:.1}%) im_killed={} ({:.1}%)",
            100.0 * buckets_contained as f64 / buckets.max(1) as f64,
            100.0 * mz_ok as f64 / denom as f64,
            100.0 * im_ok as f64 / mz_ok.max(1) as f64,
            100.0 * im_ok as f64 / denom as f64,
            cyc - mz_ok,
            100.0 * (cyc - mz_ok) as f64 / denom as f64,
            mz_ok - im_ok,
            100.0 * (mz_ok - im_ok) as f64 / denom as f64,
        );
        if buckets > 0 {
            eprintln!(
                "[fep-funnel:{label}] avg peaks-in-cycle per bucket = {:.1}, avg buckets per call = {:.2}",
                cyc as f64 / buckets as f64,
                buckets as f64 / calls.max(1) as f64,
            );
        }
    }
}

/// Zero the funnel counters — call between phases so each phase's dump
/// shows its own contribution, not cumulative. No-op when `query-instr`
/// is disabled.
pub fn reset_for_each_peak_funnel() {
    #[cfg(feature = "query-instr")]
    {
        use std::sync::atomic::Ordering::Relaxed;
        FEP_CALLS.store(0, Relaxed);
        FEP_BUCKETS.store(0, Relaxed);
        FEP_BUCKETS_MZ_CONTAINED.store(0, Relaxed);
        FEP_PEAKS_IN_CYCLE_RANGE.store(0, Relaxed);
        FEP_PEAKS_AFTER_MZ.store(0, Relaxed);
        FEP_PEAKS_AFTER_IM.store(0, Relaxed);
    }
}

/// AoS scratch → SoA.
///
/// Mobility is re-validated at flush time. In-memory values are known-good
/// here because they either came from a validated parquet column (via
/// `pack_soa_to_aos`) or from the centroider's `f16::from_f64` on a
/// non-negative converter. The unwrap is a debug-style guard.
fn flush_aos_to_soa<T: RTIndex + Copy>(peaks: Vec<IndexedPeak<T>>) -> PeakColumns<T> {
    let mut cols = PeakColumns::with_capacity(peaks.len());
    for p in peaks {
        let mob = MobInt::from_f16(p.mobility_ook0).expect("mobility must be non-neg non-NaN");
        cols.push(p.mz, p.intensity, mob, p.cycle_index);
    }
    cols
}

/// SoA → AoS scratch (inverse of `flush_aos_to_soa`).
fn pack_soa_to_aos<T: RTIndex + Copy>(
    mz: Vec<f32>,
    intensity: Vec<f32>,
    mobility: Vec<MobInt>,
    cycle_index: Vec<T>,
) -> Vec<IndexedPeak<T>> {
    let n = mz.len();
    let mut out = Vec::with_capacity(n);
    out.extend(
        mz.into_iter()
            .zip(intensity)
            .zip(mobility)
            .zip(cycle_index)
            .map(|(((mz, intensity), mobility), cycle_index)| IndexedPeak {
                mz,
                intensity,
                mobility_ook0: mobility.to_f16(),
                cycle_index,
            }),
    );
    assert_eq!(out.len(), n, "SoA columns must be length-aligned");
    out
}

/// Pre-convert an `f16` mobility query range to validated `MobInt` bounds.
/// Hoisted out of the inner scan so the filter loop does a raw u16 compare.
#[inline]
fn im_range_to_mob_bounds(
    im_range: OptionallyRestricted<TupleRange<f16>>,
) -> OptionallyRestricted<(MobInt, MobInt)> {
    im_range.map(|r| {
        let lo = MobInt::from_f16(r.start()).expect("im lo must be non-neg non-NaN");
        let hi = MobInt::from_f16(r.end()).expect("im hi must be non-neg non-NaN");
        (lo, hi)
    })
}

/// Heuristic check that a peak vec is already bucket-sorted (mz-monotonic
/// buckets, cycle-sorted within). Cheap enough to run on every `new()` to
/// skip the par-sort when the input came from a round-trip.
#[instrument(level = "info", skip_all, fields(num_peaks = peaks.len(), result))]
fn check_bucket_sorted_heuristic_aos<T: RTIndex>(
    peaks: &[IndexedPeak<T>],
    bucket_size: usize,
) -> bool {
    // 1. max mz of each bucket <= min mz of the next bucket
    // 2. each bucket is sorted by cycle_index (fully checked on bucket 0,
    //    first/last only on later buckets — cheap heuristic)
    let mut last_max = f32::MIN;
    let buckets_ordered = peaks.chunks(bucket_size).all(|bucket| {
        let curr_min = bucket
            .iter()
            .map(|x| x.mz)
            .fold(f32::INFINITY, |a, b| a.min(b));
        if curr_min < last_max {
            return false;
        }
        let curr_max = bucket
            .iter()
            .map(|x| x.mz)
            .fold(f32::NEG_INFINITY, |a, b| a.max(b));
        last_max = curr_max;
        true
    });
    if !buckets_ordered {
        tracing::Span::current().record("result", false);
        return false;
    }

    let internally_sorted_heuristic: bool = peaks.chunks(bucket_size).enumerate().all(|(i, x)| {
        if i == 0 {
            return x.windows(2).all(|w| w[0].cycle_index <= w[1].cycle_index);
        }
        let first_val = x.first().unwrap().cycle_index;
        let last_val = x.last().unwrap().cycle_index;
        last_val >= first_val
    });

    let res = internally_sorted_heuristic && buckets_ordered;
    tracing::Span::current().record("result", res);
    res
}

/// Per-bucket inner scan for `for_each_peak`. `mz_filter =
/// Unrestricted` means "the whole bucket is known to be inside the
/// query mz range, skip the per-peak mz compare". With
/// `#[inline(always)]` + the `OptionallyRestricted` split at each call
/// site, LLVM specializes two code paths — one that branches on mz
/// per peak, one that doesn't.
///
/// Filter reads only the column it needs (mz or mobility). Full
/// `IndexedPeak<T>` is materialized on survivors before calling `f`,
/// so the callback sees the same `&IndexedPeak<T>` API as the AoS
/// version. LLVM SROA elides the temporary struct when the callback
/// reads a single field.
/// Chunk size for the staged autovec path. 8 hits AVX2 (8×f32) and
/// unrolls NEON (4×f32 × 2). Remainder falls through to scalar tail.
const SCAN_CHUNK: usize = 8;

#[inline(always)]
fn apply_mz_mask<const N: usize>(
    mask: &mut [bool; N],
    mz_chunk: &[f32; N],
    range: TupleRange<f32>,
) {
    let lo = range.start();
    let hi = range.end();
    for i in 0..N {
        let v = mz_chunk[i];
        mask[i] = lo <= v && v <= hi;
    }
}

#[inline(always)]
fn apply_mob_mask<const N: usize>(
    mask: &mut [bool; N],
    mob_chunk: &[MobInt; N],
    lo: MobInt,
    hi: MobInt,
) {
    for i in 0..N {
        let v = mob_chunk[i];
        mask[i] &= lo <= v && v <= hi;
    }
}

#[inline(always)]
fn scan_bucket_slice<T, F>(
    view: PeakColumnsView<'_, T>,
    mz_filter: OptionallyRestricted<TupleRange<f32>>,
    im_range: OptionallyRestricted<TupleRange<f16>>,
    f: &mut F,
    local: &mut FepLocal,
) where
    T: RTIndex + Copy,
    F: FnMut(&IndexedPeak<T>),
{
    const N: usize = SCAN_CHUNK;
    let im_bounds = im_range_to_mob_bounds(im_range);
    let (chunks, tail) = view.as_chunks::<N>();

    // Bulk chunked path: fixed-N mask compute is autovec-friendly.
    // Unrestricted filters leave the mask at [true; N] and the counter
    // tallies N; `count_*_mask` methods are compiled out without the
    // `query-instr` feature.
    for chunk in chunks {
        let mut mask = [true; N];
        if let Restricted(r) = mz_filter {
            apply_mz_mask::<N>(&mut mask, chunk.mz(), r);
        }
        local.count_after_mz_mask::<N>(&mask);
        if let Restricted((lo, hi)) = im_bounds {
            apply_mob_mask::<N>(&mut mask, chunk.mobility(), lo, hi);
        }
        local.count_after_im_mask::<N>(&mask);
        for (i, &pass) in mask.iter().enumerate() {
            if pass {
                let peak = chunk.materialize(i);
                f(&peak);
            }
        }
    }

    // Tail scalar scan over remainder < N.
    let im_ok = |m: MobInt| match im_bounds {
        Restricted((lo, hi)) => lo <= m && m <= hi,
        Unrestricted => true,
    };
    let tail_mz = tail.mz();
    let tail_mob = tail.mobility();
    for i in 0..tail.len() {
        let mz_pass = match mz_filter {
            Restricted(r) => r.contains(tail_mz[i]),
            Unrestricted => true,
        };
        if !mz_pass {
            continue;
        }
        local.bump_after_mz(1);
        if !im_ok(tail_mob[i]) {
            continue;
        }
        local.bump_after_im();
        let peak = tail.materialize(i);
        f(&peak);
    }
}

impl<T: RTIndex> IndexedPeakGroup<T> {
    /// Query peaks based on m/z, rt, and im ranges.
    ///
    /// Returns an iterator of owned `IndexedPeak<T>` — peaks are materialized
    /// from the SoA columns on demand.
    pub fn query_peaks(
        &self,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<T>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> impl Iterator<Item = IndexedPeak<T>> {
        QueryPeaksIterator::new(self, mz_range, cycle_range, im_range)
    }

    /// Callback-style peak scan. Mirrors `query_peaks` but fuses the
    /// consumer body into the inner bucket loop — no
    /// `Iterator::next` call boundary per peak. Per the flamegraph,
    /// `QueryPeaksIterator::next` takes ~63% of wall at
    /// `RAYON_NUM_THREADS=1`; inlining the consumer via `#[inline]`
    /// lets LLVM schedule the mz/im filter + user body as a single
    /// basic block.
    #[inline]
    pub fn for_each_peak<F>(
        &self,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<T>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
        mut f: F,
    ) where
        F: FnMut(&IndexedPeak<T>),
    {
        let mut local = FepLocal::default();
        let bucket_range = self.query_bucket_range(mz_range);
        for bucket_idx in bucket_range {
            let Some(bucket) = self.get_bucket(bucket_idx) else {
                continue;
            };
            // Bucket fully inside query mz -> skip per-peak mz compare.
            // `scan_bucket_slice` has two specialized bodies the LLVM
            // inliner folds the `OptionallyRestricted` match down to.
            let bucket_mz_contained = mz_range.encloses(self.bucket_mz_ranges[bucket_idx]);
            local.bump_bucket(bucket_mz_contained);
            let (start, end) = match cycle_range.as_ref() {
                Restricted(cr) => {
                    let r = bucket.find_cycle_range(*cr);
                    (r.start, r.end)
                }
                Unrestricted => (0, bucket.len()),
            };
            local.bump_in_cycle((end - start) as u32);
            let sub = bucket.slice(start..end);
            if bucket_mz_contained {
                scan_bucket_slice(sub, Unrestricted, im_range, &mut f, &mut local);
            } else {
                scan_bucket_slice(sub, Restricted(mz_range), im_range, &mut f, &mut local);
            }
        }
        local.flush();
    }

    /// Create a new IndexedPeakGroup from a vector of peaks.
    ///
    /// This is the **canonical** way to build an IndexedPeakGroup - all bucketing
    /// logic lives here to ensure consistency.
    ///
    /// Sort + bucket happen in AoS form (`Vec<IndexedPeak<T>>`), then the four
    /// columns are unpacked into SoA storage. Callers that already have column
    /// data (e.g. parquet read) should use `new_from_soa` instead to skip the
    /// AoS round-trip.
    ///
    /// NOTE: This internally uses `par_sort_unstable` to sort the peaks
    /// so in theory it should not be called within a parallel loop.
    #[instrument(level = "info", skip_all, fields(num_peaks = peaks.len()))]
    pub(crate) fn new(
        mut peaks: Vec<IndexedPeak<T>>,
        cycle_to_rt_ms: CycleToRTMapping<T>,
        bucket_size: usize,
    ) -> (Self, IndexedPeakGroupStats) {
        let st = std::time::Instant::now();

        // Fast-pass check for already-sorted-and-bucketed inputs (e.g. freshly
        // deserialized). If it passes we skip the par_sort.
        let needs_sorting = !check_bucket_sorted_heuristic_aos(&peaks, bucket_size);
        if needs_sorting {
            peaks.par_sort_unstable_by(|x, y| x.mz.partial_cmp(&y.mz).unwrap());
            assert!(
                peaks.first().unwrap().mz <= peaks.last().unwrap().mz,
                "Peaks should be sorted by m/z after sorting step {:?} - {:?} [{};{}]",
                peaks.first().unwrap(),
                peaks.last().unwrap(),
                peaks.len(),
                bucket_size
            );
        }
        // Cheap but worth it — <40ms even on Hela.
        assert!(peaks.iter().all(|x| x.intensity >= 0.0));
        let max_cycle = T::new(cycle_to_rt_ms.len() as u32 - 1);
        assert!(peaks.iter().all(|x| x.cycle_index <= max_cycle));
        let sort_time = st.elapsed();

        let st = std::time::Instant::now();
        let bucket_mz_ranges: Vec<_> = peaks
            .par_chunks_mut(bucket_size)
            .map(|chunk| {
                let (start, end) = if needs_sorting {
                    let start = chunk.first().unwrap().mz;
                    let end = chunk.last().unwrap().mz;
                    chunk.sort_unstable_by(|x, y| {
                        x.cycle_index
                            .partial_cmp(&y.cycle_index)
                            .unwrap()
                            .then(x.mobility_ook0.partial_cmp(&y.mobility_ook0).unwrap())
                    });
                    (start, end)
                } else {
                    let mut start = f32::MAX;
                    let mut end = f32::MIN;
                    for peak in chunk.iter() {
                        if peak.mz < start {
                            start = peak.mz;
                        }
                        if peak.mz > end {
                            end = peak.mz;
                        }
                    }
                    (start, end)
                };
                TupleRange::try_new(start, end).expect("Incoming vec should have been sorted")
            })
            .collect();

        if needs_sorting {
            assert!(
                check_bucket_sorted_heuristic_aos(&peaks, bucket_size),
                "Peaks should be bucket sorted after bucketing step"
            );
        }

        let bucket_time = st.elapsed();

        let cols = flush_aos_to_soa(peaks);

        let tmp = Self {
            cols,
            bucket_size,
            bucket_mz_ranges,
            cycle_to_rt_ms,
        };
        let stats = IndexedPeakGroupStats {
            num_peaks: tmp.num_peaks(),
            num_buckets: tmp.bucket_mz_ranges.len(),
            memory_usage_bytes: tmp.aproximate_memory_usage(),
            sorting_time: sort_time,
            bucketing_time: bucket_time,
        };
        (tmp, stats)
    }

    /// Construct from owned SoA columns. Skips the AoS pack step when the
    /// caller already has column-oriented data (e.g. parquet reader).
    pub(crate) fn new_from_soa(
        cols: PeakColumns<T>,
        cycle_to_rt_ms: CycleToRTMapping<T>,
        bucket_size: usize,
    ) -> (Self, IndexedPeakGroupStats) {
        let (mz, intensity, mobility, cycle_index) = cols.into_parts();
        let peaks = pack_soa_to_aos(mz, intensity, mobility, cycle_index);
        Self::new(peaks, cycle_to_rt_ms, bucket_size)
    }

    /// Borrow the four SoA columns as aligned slices.
    pub(crate) fn columns(&self) -> PeakColumnsView<'_, T> {
        self.cols.view()
    }

    /// Rebuild the bucket layout at a different `new_bucket_size`.
    ///
    /// Consumes `self` and returns a new `IndexedPeakGroup` with the
    /// same peaks, re-chunked. Useful to adjust bucket granularity
    /// between query stages: narrow-tolerance queries (Phase 3) prefer
    /// small buckets, wide queries (Phase 1) tolerate larger ones.
    ///
    /// Cost: pack SoA→AoS scratch, par_sort by mz, par_chunks_mut cycle-sort,
    /// flush AoS→SoA. Dominated by the O(N log N) mz sort.
    pub fn rebucket(self, new_bucket_size: usize) -> Self {
        assert!(new_bucket_size > 0, "bucket_size must be > 0");
        if new_bucket_size == self.bucket_size {
            return self;
        }
        let Self {
            cols,
            cycle_to_rt_ms,
            ..
        } = self;
        let (mz, intensity, mobility, cycle_index) = cols.into_parts();
        let mut peaks = pack_soa_to_aos(mz, intensity, mobility, cycle_index);
        peaks.par_sort_unstable_by(|x, y| x.mz.total_cmp(&y.mz));
        let bucket_mz_ranges: Vec<TupleRange<f32>> = peaks
            .par_chunks_mut(new_bucket_size)
            .map(|chunk| {
                let start = chunk.first().unwrap().mz;
                let end = chunk.last().unwrap().mz;
                chunk.sort_unstable_by(|x, y| {
                    x.cycle_index
                        .partial_cmp(&y.cycle_index)
                        .unwrap()
                        .then(x.mobility_ook0.total_cmp(&y.mobility_ook0))
                });
                TupleRange::try_new(start, end).unwrap()
            })
            .collect();
        let cols = flush_aos_to_soa(peaks);
        Self {
            cols,
            bucket_mz_ranges,
            bucket_size: new_bucket_size,
            cycle_to_rt_ms,
        }
    }

    /// Decompose into raw SoA columns plus metadata. Inverse of `new_from_soa`
    /// in shape (no re-sort).
    pub fn unpack(self) -> (PeakColumns<T>, Vec<u32>, usize, Vec<TupleRange<f32>>) {
        let Self {
            cols,
            bucket_mz_ranges,
            bucket_size,
            cycle_to_rt_ms,
        } = self;
        (cols, cycle_to_rt_ms.unpack(), bucket_size, bucket_mz_ranges)
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
        // All four columns have equal length; use the row count × per-peak bytes.
        let per_peak_bytes = std::mem::size_of::<f32>() * 2
            + std::mem::size_of::<MobInt>()
            + std::mem::size_of::<T>();
        let col_mem = self.cols.len() * per_peak_bytes;
        self_mem + bucket_mem + col_mem
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
        let mut last_pushed_cycle = None;

        for (i, meta) in frame_reader.frame_metas.iter().enumerate() {
            if matches!(meta.ms_level, timsrust::MSLevel::MS1) {
                cycle_index += 1;
            }
            if filter(meta) {
                if let Some(last_cycle) = last_pushed_cycle
                    && last_cycle == cycle_index
                {
                    warn!(
                        "Found multiple frames in the same cycle,
                        skipping the non-first ones (cycle: {}, rt: {})
                        this might point to an error when collecting the data
                        (missing ms1 frame)
                        ",
                        last_cycle, meta.rt_in_seconds,
                    );
                    continue;
                }
                let rt_ms = (meta.rt_in_seconds * 1000.0).round() as u32;
                out.push((i, cycle_index, rt_ms));
                last_pushed_cycle = Some(cycle_index);
            }
        }
        // TODO: check that cycle indices are continuous AND without replicates
        // or gaps ...
        out.as_slice().windows(2).for_each(|w| {
            // Note: These are assertions instead of errors because
            // they are invariants I am checking, not recoverable errors, none of the logic
            // after this point would make sense if this is not true.
            //
            // This checks that whatever filter is used returns only a single frame within each
            // cycle ...
            assert!(
                w[0].1 < w[1].1,
                "Cycle indices should be strictly increasing, {:?} vs {:?}",
                w[0],
                w[1],
            );
            assert!(
                w[0].2 <= w[1].2,
                "Retention times should be non-decreasing, {:?} vs {:?}",
                w[0],
                w[1],
            );
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
    fn get_bucket(&self, bucket_idx: usize) -> Option<PeakColumnsView<'_, T>> {
        let r = self.get_bucket_range(bucket_idx)?;
        Some(self.cols.view().slice(r))
    }

    /// Total number of peaks.
    #[inline]
    pub fn num_peaks(&self) -> usize {
        self.cols.len()
    }

    /// Get the range of peak indices for a specific bucket index.
    fn get_bucket_range(&self, bucket_idx: usize) -> Option<std::ops::Range<usize>> {
        let start = bucket_idx * self.bucket_size;
        if start >= self.num_peaks() {
            return None;
        }
        let end = ((bucket_idx + 1) * self.bucket_size).min(self.num_peaks());
        Some(start..end)
    }

    fn print_glimpse(&self) {
        let num_buckets = self.bucket_mz_ranges.len();
        let num_peaks = self.num_peaks();
        // Four columns: f32 mz + f32 intensity + u16 mobility + sizeof<T> cycle_index.
        let per_peak_bytes = std::mem::size_of::<f32>() * 2
            + std::mem::size_of::<MobInt>()
            + std::mem::size_of::<T>();
        let mem_usage = num_peaks * per_peak_bytes;
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

#[derive(Debug)]
struct QueryPeaksIterator<'a, T: RTIndex> {
    indexed_window_group: &'a IndexedPeakGroup<T>,
    mz_range: TupleRange<f32>,
    cycle_range: OptionallyRestricted<TupleRange<T>>,
    im_bounds: OptionallyRestricted<(MobInt, MobInt)>,
    bucket_idx: usize,
    bucket_end: usize,
    position_in_bucket: usize,
    end_of_current_bucket: usize,
    current_bucket: Option<PeakColumnsView<'a, T>>,
}

impl<'a, T: RTIndex> QueryPeaksIterator<'a, T> {
    pub fn new(
        indexed_window_group: &'a IndexedPeakGroup<T>,
        mz_range: TupleRange<f32>,
        cycle_range: OptionallyRestricted<TupleRange<T>>,
        im_range: OptionallyRestricted<TupleRange<f16>>,
    ) -> Self {
        let bucket_range = indexed_window_group.query_bucket_range(mz_range);
        let im_bounds = im_range_to_mob_bounds(im_range);

        Self {
            indexed_window_group,
            mz_range,
            cycle_range,
            im_bounds,
            bucket_idx: bucket_range.start,
            bucket_end: bucket_range.end,
            position_in_bucket: 0,
            end_of_current_bucket: 0,
            current_bucket: None,
        }
    }
}

impl<'a, T: RTIndex + Copy> QueryPeaksIterator<'a, T> {
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

    fn next_in_current_bucket(&mut self) -> Option<IndexedPeak<T>> {
        // Loop (not recurse) so long filter-reject streaks in dense spectra
        // don't blow the stack.
        while let Some(bucket) = self.current_bucket.as_ref() {
            if self.position_in_bucket >= self.end_of_current_bucket {
                return None;
            }

            let i = self.position_in_bucket;
            self.position_in_bucket += 1;

            if !self.mz_range.contains(bucket.mz()[i]) {
                continue;
            }
            let mob = bucket.mobility()[i];
            let im_ok = match self.im_bounds {
                Restricted((lo, hi)) => lo <= mob && mob <= hi,
                Unrestricted => true,
            };
            if !im_ok {
                continue;
            }
            return Some(bucket.materialize(i));
        }
        None
    }
}

impl<'a, T: RTIndex + Copy> Iterator for QueryPeaksIterator<'a, T> {
    type Item = IndexedPeak<T>;

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
        let cols = flush_aos_to_soa(peaks);
        let bucket = cols.view();
        let cycle_range = (MS1CycleIndex::new(1), MS1CycleIndex::new(3));
        let rt_idx_range = bucket.find_cycle_range(cycle_range.try_into().unwrap());
        let out: Vec<IndexedPeak<MS1CycleIndex>> =
            rt_idx_range.map(|i| bucket.materialize(i)).collect();
        let expected = tuples_to_peaks::<MS1CycleIndex>(&[
            (100.0, 250.0, 1.0, 1),
            (100.0, 250.0, 1.0, 1),
            (100.0, 250.0, 1.0, 1),
            (100.0, 300.0, 1.0, 2),
            (100.0, 350.0, 1.0, 3),
        ]);

        assert_eq!(out, expected);
    }

    #[test]
    fn test_checking_bucketing() {
        let test_data = vec![
            // mz, intensity, im, cycle_index
            // So everything here is sorted by both mz and cycle_index
            (100.0, 200.0, 1.0, 0u32),
            (100.0, 250.0, 1.0, 1u32),
            (100.0, 250.0, 1.0, 1u32),
            (100.0, 250.0, 1.0, 1u32),
            (150.0, 300.0, 1.0, 2u32),
            (150.0, 350.0, 1.0, 3u32),
            (200.0, 400.0, 1.0, 4u32),
        ];
        let peaks = tuples_to_peaks::<MS1CycleIndex>(&test_data);
        assert!(check_bucket_sorted_heuristic_aos(&peaks, 4));

        let test_data = vec![
            // mz, intensity, im, cycle_index
            // all mz's from bucket N are less than (or equal to) all mz's from bucket N+1
            (200.0, 200.0, 1.0, 0u32),
            (100.0, 250.0, 1.0, 2u32),
            (50.0, 250.0, 1.0, 23u32),
            // End bucket 1 (internally sorted by cycle, not mz)
            (2000.0, 250.0, 1.0, 1u32),
            (1900.0, 300.0, 1.0, 2u32),
            (1800.0, 350.0, 1.0, 3u32),
            // End bucket 2
            (20_000.0, 400.0, 1.0, 4u32),
        ];
        let peaks = tuples_to_peaks::<MS1CycleIndex>(&test_data);
        assert!(check_bucket_sorted_heuristic_aos(&peaks, 3));
        // This should fail because in bucket 2 we have mz's less than in bucket 1
        assert!(!check_bucket_sorted_heuristic_aos(&peaks, 2));
        // Any bucket number other than 3 should fail
        for bucket_size in [1, 2, 4, 5, 10] {
            if bucket_size != 3 {
                assert!(
                    !check_bucket_sorted_heuristic_aos(&peaks, bucket_size),
                    "Bucket size {} should fail the check",
                    bucket_size
                );
            }
        }
    }
}
