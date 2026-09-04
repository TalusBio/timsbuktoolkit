use std::fmt::Display;

use timsrust::CorrectedFramePeak;
use timsrust::converters::ConvertableDomain;

/// Buffer that gets re-used on each thread to store the intermediates
/// of the centroiding for a single frame.
///
/// IN GENERAL, you should not be using this directly BUT it is exposed
pub struct PeakCentroider<T1: ConvertableDomain, T2: ConvertableDomain> {
    peaks: Vec<CorrectedFramePeak>,
    order: Vec<usize>,
    order_intensity: Vec<f64>,
    taken_buff: Vec<TakenState>,
    agg_buff: Vec<PeakAggregator>,
    neighbor_ranges: Vec<(usize, usize)>, // (start_idx, end_idx) for each peak
    ims_ranges: Vec<(u16, u16)>,          // (min_im, max_im) for scan index
    max_peaks: usize,
    mz_tol: MzTolerance,
    im_tol: ImTolerance,
    mz_converter: T1,
    im_converter: T2,
    early_stop_iterations: u32,
    window_cap: Option<WindowCap>,
    transitive: bool,
    /// Parents claimed so far in each `window_da`-wide m/z window of the
    /// current frame. Indexed by `mz / window_da`, grown on demand.
    window_counts: Vec<u32>,
}

/// Configuration for the centroiding algorithm
///
/// For details on the algorithm, check the documentation
/// of the [PeakCentroider::centroid_frame] method.
#[derive(Clone, Copy, Debug)]
pub struct CentroidingConfig {
    /// Maximum number of peaks to retain after centroiding
    /// This is a hard limit, if the number of peaks after
    /// centroiding is higher, the centroiding will stop
    /// early. A number ~20,000-50,000 seems to work well in practice.
    pub max_peaks: usize,
    /// m/z neighborhood for clustering, see [`MzTolerance`].
    pub mz_tol: MzTolerance,
    /// Mobility neighborhood for clustering, see [`ImTolerance`].
    pub im_tol: ImTolerance,
    /// Number of consecutive iterations with no new peaks clustered
    /// after which the centroiding will stop early (instead of going
    /// through all peaks, which will very likely be noise).
    /// A number ~200 seems to work well in practice. **`0` disables
    /// early-stop** -- the whole frame is clustered and only `max_peaks`
    /// can truncate it. Disabling is the right choice under a tight
    /// `mz_ppm_tol`, where most peaks are singletons and early-stop would
    /// otherwise clip real signal (see `IndexingCentroidingConfig`).
    pub early_stop_iterations: u32,
    /// Cap on centroids per m/z window, per frame. `None` leaves only
    /// `max_peaks` to bound a frame.
    ///
    /// `max_peaks` ranks the whole frame by intensity, so a dense region can
    /// spend the entire budget and leave a sparse one with nothing. A window
    /// cap spends the budget per region: the `max_peaks` most intense
    /// centroids in each `window_da` slice survive, whatever the rest of the
    /// frame looks like. Stored peaks per frame are then bounded by
    /// `windows * max_peaks` independent of the sample.
    pub window_cap: Option<WindowCap>,
    /// Whether a peak already absorbed into a cluster may recruit its own
    /// neighbors into that cluster (DBSCAN-style transitive growth).
    ///
    /// `true` is the historical behavior: a cluster can grow one tolerance
    /// step at a time as far as the peak density carries it, and its centroid
    /// is the weighted mean of everything it swallowed. `false` bounds a
    /// cluster to the parent's own tolerance box, so its radius is the
    /// tolerance and nothing more. With a sub-bin m/z tolerance the two are
    /// nearly the same, since chaining can only run along scans at one TOF
    /// index. With a one-bin tolerance chaining runs along TOF through
    /// crowded regions and pulls centroids off their ions.
    pub transitive: bool,
}

/// How far apart in m/z two raw peaks may be and still cluster.
///
/// The stored data is line spectra: the firmware has already picked one
/// integer TOF index per ion per scan, and a TOF bin is 2.4 to 7 ppm wide
/// across the m/z range (narrower in ppm at high m/z, since m/z grows with the
/// square of the index). The same ion on two scans lands in adjacent bins when
/// its true position sits near a bin edge, so the natural merge radius is a
/// whole number of bins, the same everywhere on the axis.
///
/// `Ppm` is converted to a bin delta per peak and rounded; below half a bin
/// it rounds to zero and only identical TOF indices merge. `Bins` skips the
/// conversion entirely.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum MzTolerance {
    /// Symmetric tolerance in parts per million of the peak's m/z.
    Ppm(f64),
    /// Symmetric tolerance in raw TOF bins.
    Bins(u32),
}

/// How far apart in mobility two raw peaks may be and still cluster.
///
/// Along the scan axis the data is a real profile: one ion appears on
/// consecutive scans as it leaves the TIMS tunnel, typically 5 to 20 of them.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum ImTolerance {
    /// Symmetric tolerance as a percentage of the peak's 1/K0.
    Pct(f64),
    /// Symmetric tolerance in scan indices.
    Scans(u16),
}

/// Per-window centroid quota; see [`CentroidingConfig::window_cap`].
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindowCap {
    /// Centroids kept per window per frame.
    pub max_peaks: u32,
    /// Window width in m/z. Windows tile from zero, so a 100 Da width puts
    /// the boundaries at 100, 200, 300, ...
    pub window_da: f32,
}

impl Default for CentroidingConfig {
    fn default() -> Self {
        Self {
            max_peaks: 20_000,
            mz_tol: MzTolerance::Ppm(5.0),
            im_tol: ImTolerance::Pct(3.0),
            early_stop_iterations: 200,
            window_cap: None,
            transitive: true,
        }
    }
}

/// Per-MS-level centroiding configuration for building an index.
///
/// MS1 and MS2 are centroided independently: MS1 favors precursor m/z
/// precision (peak counts are small, so a tight merge is cheap), while MS2
/// keeps a tight m/z merge with early-stop disabled and a moderate
/// `max_peaks` cap -- the setting that preserves fragment signal without the
/// peak-count explosion of a full raw pass.
///
/// Use [`IndexingCentroidingConfig::uniform`] to apply one config to both
/// levels (e.g. for non-TIMS sources or tests).
#[derive(Clone, Copy, Debug)]
pub struct IndexingCentroidingConfig {
    pub ms1: CentroidingConfig,
    pub ms2: CentroidingConfig,
    /// Build the index from this retention-time slice only, in seconds.
    ///
    /// A cycle is in the slice when its MS1 frame is; MS2 frames follow their
    /// cycle, so every window group covers the same cycles and cycle indices
    /// are renumbered from the first one in the slice. Not a centroiding
    /// parameter, but this is the one struct that reaches the raw build, and
    /// a two-minute slice turns a benchmark iteration from minutes into
    /// seconds. A sliced index is never written as an `.idx` sidecar.
    pub rt_range_seconds: Option<(f64, f64)>,
}

impl IndexingCentroidingConfig {
    /// Apply the same `CentroidingConfig` to both MS levels.
    pub fn uniform(cfg: CentroidingConfig) -> Self {
        Self {
            ms1: cfg,
            ms2: cfg,
            rt_range_seconds: None,
        }
    }
}

impl Default for IndexingCentroidingConfig {
    /// timsTOF defaults. MS1: identical-index m/z merge, 2% mobility, 20k
    /// per frame, early-stop on. MS2: identical-index merge, 3% mobility, 50k
    /// per frame as a backstop, early-stop off, and a 500-per-50-Da window
    /// cap, which is what bounds index memory. The tolerance and cap history
    /// is in `DATA_MODEL.md`.
    fn default() -> Self {
        Self {
            ms1: CentroidingConfig {
                max_peaks: 20_000,
                // Identical TOF index only. The sweep that tuned these found
                // "less merge, more IDs"; this is the explicit spelling of the
                // 0.5 / 1 ppm it landed on, which round to zero bins anyway.
                mz_tol: MzTolerance::Bins(0),
                im_tol: ImTolerance::Pct(2.0),
                early_stop_iterations: 200,
                window_cap: None,
                transitive: true,
            },
            rt_range_seconds: None,
            ms2: CentroidingConfig {
                // 500 centroids per 50 Da per frame. On a 21-min diaPASEF HeLa
                // run this took the index from 13.4 GB to 4.5 GB and raised
                // 1% IDs by ~10%: what it removes is dim clutter the scorer
                // was not using. MS1 is left uncapped by window; capping it
                // cost IDs. See timscentroid/DATA_MODEL.md.
                window_cap: Some(WindowCap {
                    max_peaks: 500,
                    window_da: 50.0,
                }),
                transitive: true,
                max_peaks: 50_000,
                mz_tol: MzTolerance::Bins(0),
                im_tol: ImTolerance::Pct(3.0),
                early_stop_iterations: 0, // disabled -- see field doc
            },
        }
    }
}

#[derive(Clone)]
struct PeakAggregator {
    weighed_tof_sum: u64,
    weighed_im_sum: u64,
    total_weight: u64,
    total_intensity: f64,
}

impl PeakAggregator {
    fn new(peak: &CorrectedFramePeak) -> Self {
        assert!(peak.corrected_intensity >= 0.0);

        let weight = peak.corrected_intensity as u64;
        Self {
            weighed_tof_sum: peak.tof_index as u64 * weight,
            weighed_im_sum: peak.scan_index as u64 * weight,
            total_weight: weight,
            total_intensity: peak.corrected_intensity,
        }
    }

    fn add_peak(&mut self, peak: &CorrectedFramePeak) {
        assert!(peak.corrected_intensity >= 0.0);

        let weight = peak.corrected_intensity as u64;

        // Use checked arithmetic to prevent silent overflow in weight accumulation.
        // With ~500k peaks at high intensities, overflow is theoretically possible.
        self.weighed_tof_sum = self
            .weighed_tof_sum
            .checked_add(peak.tof_index as u64 * weight)
            .expect("Weight sum overflow: tof_index * weight exceeded u64::MAX");

        self.weighed_im_sum = self
            .weighed_im_sum
            .checked_add(peak.scan_index as u64 * weight)
            .expect("Weight sum overflow: scan_index * weight exceeded u64::MAX");

        self.total_weight = self
            .total_weight
            .checked_add(weight)
            .expect("Total weight overflow: accumulated weight exceeded u64::MAX");

        self.total_intensity += peak.corrected_intensity;
    }

    fn finalize(&self) -> CorrectedFramePeak {
        CorrectedFramePeak {
            tof_index: self.calc_tof_index(),
            scan_index: (self.weighed_im_sum / self.total_weight) as u16,
            corrected_intensity: self.total_intensity,
        }
    }

    fn calc_tof_index(&self) -> u32 {
        (self.weighed_tof_sum / self.total_weight) as u32
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TakenState {
    Untaken,
    Taken {
        parent_idx: usize,
    },
    Parent,
    /// Would have been a parent, but its m/z window was at quota. Never
    /// absorbed by a later parent either: every later parent is less intense,
    /// and a centroid must not be pulled toward a peak that outranks it.
    Dropped,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StoppingReason {
    EarlyStop,
    MaxPeaks,
    AllTaken,
}

/// Summary of the clustering for a single frame
///
/// This will usually be returned alongside the centroided peaks
/// from the [PeakCentroider::centroid_frame] method.
#[derive(Debug, Clone, Copy)]
pub struct ClusteringSummary {
    pub initial_peaks: usize,
    pub aggregated_peaks: usize,
    pub final_peaks: usize,
    /// Peaks that would have become centroids but fell to the window cap.
    pub window_capped: usize,
    pub stopping_reason: StoppingReason,
    pub elapsed: std::time::Duration,
}

/// Summary of the clustering over multiple frames
#[derive(Debug, Clone)]
pub struct AggregatedClusteringSummary {
    pub total_initial_peaks: usize,
    pub total_aggregated_peaks: usize,
    pub total_final_peaks: usize,
    pub total_window_capped: usize,
    pub frames_processed: usize,
    pub stopping_reasons: [usize; 3],
    pub elapsed: std::time::Duration,
}

impl AggregatedClusteringSummary {
    pub(crate) fn combine(mut self, right: &Self) -> Self {
        self.total_final_peaks += right.total_final_peaks;
        self.total_aggregated_peaks += right.total_aggregated_peaks;
        self.total_initial_peaks += right.total_initial_peaks;
        self.total_window_capped += right.total_window_capped;
        self.frames_processed += right.frames_processed;
        self.elapsed += right.elapsed;
        for (i, v) in right.stopping_reasons.iter().enumerate() {
            self.stopping_reasons[i] += v;
        }
        self
    }

    fn reason_to_idx(reason: StoppingReason) -> usize {
        match reason {
            StoppingReason::EarlyStop => 0,
            StoppingReason::MaxPeaks => 1,
            StoppingReason::AllTaken => 2,
        }
    }

    pub(crate) fn fold_summary(mut left: Self, other: &ClusteringSummary) -> Self {
        left.total_initial_peaks += other.initial_peaks;
        left.total_aggregated_peaks += other.aggregated_peaks;
        left.total_final_peaks += other.final_peaks;
        left.total_window_capped += other.window_capped;
        left.frames_processed += 1;
        let reason_idx = Self::reason_to_idx(other.stopping_reason);
        left.stopping_reasons[reason_idx] += 1;
        left.elapsed += other.elapsed;
        left
    }

    pub(crate) fn new() -> Self {
        Self {
            total_initial_peaks: 0,
            total_aggregated_peaks: 0,
            total_final_peaks: 0,
            total_window_capped: 0,
            frames_processed: 0,
            stopping_reasons: [0; 3],
            elapsed: std::time::Duration::new(0, 0),
        }
    }
}

impl From<ClusteringSummary> for AggregatedClusteringSummary {
    fn from(value: ClusteringSummary) -> Self {
        let mut out = Self::new();
        out = Self::fold_summary(out, &value);
        out
    }
}

impl Display for AggregatedClusteringSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Frames processed: {}", self.frames_processed)?;
        writeln!(
            f,
            "Total elapsed time centroiding (cpu): {:#.2?}",
            self.elapsed
        )?;
        let avg_time = self.elapsed / self.frames_processed as u32;
        writeln!(f, "Average time per frame: {:#.2?}", avg_time)?;
        writeln!(f, "Total initial peaks: {}", self.total_initial_peaks)?;
        writeln!(f, "Total aggregated peaks: {}", self.total_aggregated_peaks)?;
        writeln!(f, "Total final peaks: {}", self.total_final_peaks)?;
        writeln!(f, "Total window-capped peaks: {}", self.total_window_capped)?;
        writeln!(
            f,
            "Average reduction: {:.2}x",
            self.total_initial_peaks as f64 / self.total_final_peaks as f64
        )?;
        writeln!(
            f,
            "Stopping reasons: EarlyStop: {}, MaxPeaks: {}, AllTaken: {}",
            self.stopping_reasons[Self::reason_to_idx(StoppingReason::EarlyStop)],
            self.stopping_reasons[Self::reason_to_idx(StoppingReason::MaxPeaks)],
            self.stopping_reasons[Self::reason_to_idx(StoppingReason::AllTaken)],
        )?;
        Ok(())
    }
}

impl<T1: ConvertableDomain, T2: ConvertableDomain> PeakCentroider<T1, T2> {
    pub fn with_capacity(
        capacity: usize,
        config: CentroidingConfig,
        mz_converter: T1,
        im_converter: T2,
    ) -> Self {
        Self {
            peaks: Vec::with_capacity(capacity),
            order: Vec::with_capacity(capacity),
            order_intensity: Vec::with_capacity(capacity),
            taken_buff: Vec::with_capacity(capacity),
            neighbor_ranges: Vec::with_capacity(capacity),
            agg_buff: Vec::with_capacity(config.max_peaks),
            // Note there are A LOT less scan indices than peaks
            // so we can cache all values in the 0-max scan index range
            // for each peak ... (and even extend only if needed)
            ims_ranges: Vec::new(),
            max_peaks: config.max_peaks,
            mz_tol: config.mz_tol,
            im_tol: config.im_tol,
            mz_converter,
            im_converter,
            early_stop_iterations: config.early_stop_iterations,
            window_cap: config.window_cap,
            transitive: config.transitive,
            window_counts: Vec::new(),
        }
    }

    /// Charges one parent slot in the m/z window holding peak `idx`. Returns
    /// `false` when that window is already at quota. Always `true` without a
    /// window cap.
    fn try_claim_window_slot(&mut self, idx: usize) -> bool {
        let Some(cap) = self.window_cap else {
            return true;
        };
        let mz = self.mz_converter.convert(self.peaks[idx].tof_index);
        let w = (mz / cap.window_da as f64).max(0.0) as usize;
        if w >= self.window_counts.len() {
            self.window_counts.resize(w + 1, 0);
        }
        if self.window_counts[w] >= cap.max_peaks {
            return false;
        }
        self.window_counts[w] += 1;
        true
    }

    /// Given a TOF index, returns the (inclusive) bounds of TOF indices
    /// that fall within the ppm tolerance of the mz corresponding to the TOF index.
    /// This is used to find neighboring peaks for clustering.
    fn tof_index_bounds(&self, tof_idex: u32) -> (u32, u32) {
        match self.mz_tol {
            MzTolerance::Bins(n) => (tof_idex.saturating_sub(n), tof_idex.saturating_add(n)),
            MzTolerance::Ppm(ppm) => {
                let mz = self.mz_converter.convert(tof_idex as f64);
                let delta_mz = mz * ppm * 1e-6;
                let left_mz = mz - delta_mz;
                let right_mz = mz + delta_mz;
                let left_tof = self.mz_converter.invert(left_mz).round() as u32;
                let right_tof = self.mz_converter.invert(right_mz).round() as u32;
                (left_tof, right_tof)
            }
        }
    }

    fn im_index_bounds(&self, im_index: u16) -> (u16, u16) {
        // Invariant: ims_ranges is populated by maybe_extend_ims_ranges() in with_frame()
        // before any lookups occur. If this panics, it indicates a logic bug in cache setup.
        self.ims_ranges[im_index as usize]
    }

    fn uncached_im_index_bounds(&self, im_index: u16) -> (u16, u16) {
        let pct = match self.im_tol {
            ImTolerance::Scans(n) => {
                return (im_index.saturating_sub(n), im_index.saturating_add(n));
            }
            ImTolerance::Pct(pct) => pct,
        };
        let im = self.im_converter.convert(im_index as f64);
        let delta_im = im * pct * 0.01;
        let left_im = im - delta_im;
        let right_im = im + delta_im;
        let left_im_index = self.im_converter.invert(left_im).round() as u16;
        let right_im_index = self.im_converter.invert(right_im).round() as u16;
        // Note: low ims have higher index values
        (
            left_im_index.min(right_im_index),
            right_im_index.max(left_im_index),
        )
    }

    fn maybe_extend_ims_ranges(&mut self) {
        let max_ims_index = self.peaks.iter().map(|x| x.scan_index).max().unwrap_or(0);
        let curr_len = self.ims_ranges.len();
        if max_ims_index as usize >= curr_len {
            self.ims_ranges.resize(max_ims_index as usize + 1, (0, 0));
            for idx in curr_len..=max_ims_index as usize {
                let bounds = self.uncached_im_index_bounds(idx as u16);
                self.ims_ranges[idx] = bounds;
            }
        }
    }

    /// Carries out the setup of the internal buffers with the frame
    /// to be centroided.
    fn with_frame(&mut self, frame: &timsrust::Frame) {
        self.clear();
        let expect_len = frame.peaks.len();
        self.expand_to_capacity(expect_len);
        self.peaks.extend(frame.iter_corrected_peaks());
        assert_eq!(self.peaks.len(), expect_len);
        self.maybe_extend_ims_ranges();

        // sort by mz ... bc binary searching on the mz space
        // for neighbors is the fastest way to find neighbors that I have tried.
        self.peaks.sort_unstable_by_key(|a| a.tof_index);
        self.compute_neighbor_ranges_and_intensity();
        // self.compute_neighbor_ranges();
        // self.compute_neighborhood_intensity();

        // The "order" is sorted by intensity
        // This will be used later during the centroiding (for details check that implementation)
        self.set_orders();
    }

    fn set_orders(&mut self) {
        self.order.extend(0..self.peaks.len());
        self.order
            .sort_unstable_by(|&a, &b| self.order_intensity[b].total_cmp(&self.order_intensity[a]));
    }

    fn clear(&mut self) {
        self.peaks.clear();
        self.order.clear();
        self.taken_buff.clear();
        self.agg_buff.clear();
        self.window_counts.fill(0);
    }

    /// Expands the internal buffers to the specified capacity
    /// if they are not already at that capacity.
    /// This is useful to avoid reallocations when processing
    /// multiple frames of similar size.
    ///
    /// IN GENERAL, you should not be using this directly.
    pub fn expand_to_capacity(&mut self, capacity: usize) {
        if capacity <= self.peaks.len() {
            return;
        }
        let additional = capacity - self.peaks.len();

        self.peaks.reserve(additional);
        self.order.reserve(additional);
        self.neighbor_ranges.reserve(additional);
        self.order_intensity.reserve(additional);
        self.taken_buff.reserve(additional);
        // self.agg_buff.reserve(capacity);
    }

    /// Centroiding of the IM-containing spectra
    ///
    /// It iterativelty goes over the peaks in decreasing order of the accumulated
    /// intensity of its neighbors.
    /// During each iteration, it accumulates the intensity of the peaks surrounding
    /// the peak. If the peak already has been "taken" by a parent, it gives its
    /// "neighbor" intensity to its parent (thus should be in essence equivalent
    /// to dbscan where each cluster is represented by the cumulative intensity
    /// of all points and the position of the most intense parent).
    ///
    /// The preserved mobility and mz are a weighted average of all the
    /// peaks in the cluster (with weight being the intensity of each peak, rounded
    /// to the nearest whole number).
    ///
    /// This dramatically reduces the number of peaks in the spectra
    /// which saves a ton of memory and time when doing LFQ, since we
    /// iterate over each peak.
    /// Returns the stopping reason, the number of peaks accounted for (taken,
    /// parent or dropped) and the number dropped by the window cap.
    fn itercentroid_frame(&mut self) -> (StoppingReason, usize, usize) {
        assert!(self.agg_buff.is_empty(), "agg_buff is not empty");

        self.taken_buff
            .resize(self.peaks.len(), TakenState::Untaken);

        // Make sure the array is mz sorted ... I should delete
        // this assertions once I am confident of the implementation.
        // but tbh, its not that slow and its simple.
        // More formally ... the check is O(n) and the sort is O(n log n)
        // and the clustering O(n^2) in the worst case (but usually much better
        // since we limit the search space with ppm and pct tolerances).
        debug_assert!(
            self.peaks
                .windows(2)
                .all(|x| x[0].tof_index <= x[1].tof_index),
            "mz_array is not sorted"
        );
        assert!(self.agg_buff.is_empty(), "agg_buff is not empty");

        let mut global_num_included = 0;
        let mut num_window_capped = 0;
        let mut early_stop_remaining = self.early_stop_iterations;
        let mut out = StoppingReason::AllTaken;

        // Indexed rather than `for &idx in &self.order`: claiming a window
        // slot needs `&mut self` inside the body.
        for oi in 0..self.order.len() {
            let idx = self.order[oi];
            if self.agg_buff.len() > self.max_peaks {
                out = StoppingReason::MaxPeaks;
                break;
            }

            let mut num_includable = 0;

            let (target_index, is_self_parent) = match self.taken_buff[idx] {
                TakenState::Parent => unreachable!("This should never happen"),
                TakenState::Dropped => unreachable!("a peak is dropped only at its own visit"),
                // Non-transitive: a taken peak recruits nothing. Everything in
                // the parent's box was swept when the parent was created, so
                // there is nothing left for it to add except neighbors outside
                // that box, which is exactly the growth being refused.
                TakenState::Taken { .. } if !self.transitive => continue,
                TakenState::Taken { parent_idx } => (parent_idx, false),
                TakenState::Untaken => {
                    if !self.try_claim_window_slot(idx) {
                        self.taken_buff[idx] = TakenState::Dropped;
                        global_num_included += 1;
                        num_window_capped += 1;
                        if global_num_included == self.peaks.len() {
                            out = StoppingReason::AllTaken;
                            break;
                        }
                        continue;
                    }
                    self.taken_buff[idx] = TakenState::Parent;
                    self.agg_buff.push(PeakAggregator::new(&self.peaks[idx]));
                    num_includable += 1;
                    (self.agg_buff.len() - 1, true)
                }
            };

            let (peak_idx_left, peak_idx_right) = self.neighbor_ranges[idx];
            let search_range = peak_idx_left..peak_idx_right;

            let im_index = self.peaks[idx].scan_index;
            let (left_im, right_im) = self.im_index_bounds(im_index);

            let curr_aggregator = &mut self.agg_buff[target_index];

            for i in search_range {
                match self.taken_buff[i] {
                    TakenState::Taken { .. } | TakenState::Parent | TakenState::Dropped => {
                        continue;
                    }
                    TakenState::Untaken => {}
                }
                let im_i = self.peaks[i].scan_index;
                if im_i >= left_im && im_i <= right_im {
                    // Since by definition we are iterating in decreasing intensity order
                    // we can only add intensities "uphill" to peaks of higher intensity.
                    curr_aggregator.add_peak(&self.peaks[i]);
                    assert!(i != idx);
                    self.taken_buff[i] = TakenState::Taken {
                        parent_idx: target_index,
                    };
                    num_includable += 1;
                }
            }

            global_num_included += num_includable;

            // `early_stop_iterations == 0` disables early-stop entirely: the
            // whole frame is clustered (only `max_peaks` can still truncate).
            // The guard also avoids the `0 - 1` u32 underflow that a zero
            // config would otherwise hit on the first lonely peak.
            if is_self_parent && self.early_stop_iterations != 0 {
                if num_includable < 3 {
                    // If we have only incliuded 'self' peaks for MAX_EARLY_STOP
                    // iterations consecutively, we can stop early.
                    early_stop_remaining -= 1;
                    if early_stop_remaining == 0 {
                        out = StoppingReason::EarlyStop;
                        break;
                    }
                } else {
                    early_stop_remaining = self.early_stop_iterations;
                }
            }

            if global_num_included == self.peaks.len() {
                out = StoppingReason::AllTaken;
                break;
            }
        }

        // println!("Centroiding: Start len: {}; end len: {};", arr_len, result.len());
        // Ultra data is usually start: 40k end 10k,
        // HT2 data is usually start 400k end 40k, limiting to 10k
        // rarely leaves peaks with intensity > 200 ... ive never seen
        // it happen. -JSP 2025-Jan
        (out, global_num_included, num_window_capped)
    }

    fn drain_aggregated_peaks(&mut self) -> impl Iterator<Item = CorrectedFramePeak> + '_ {
        self.agg_buff.drain(..).map(|agg| agg.finalize())
    }

    /// Centroiding of the IM-containing spectra
    ///
    /// It iterativelty goes over the peaks in decreasing order of the accumulated
    /// intensity of its neighbors.
    /// During each iteration, it accumulates the intensity of the peaks surrounding
    /// the peak. If the peak already has been "taken" by a parent, it gives its
    /// "neighbor" intensity to its parent (thus should be in essence equivalent
    /// to dbscan where each cluster is represented by the cumulative intensity
    /// of all points and the position of the most intense parent).
    ///
    /// The preserved mobility and mz are a weighted average of all the
    /// peaks in the cluster (with weight being the intensity of each peak, rounded
    /// to the nearest whole number).
    ///
    /// This dramatically reduces the number of peaks in the spectra
    /// which saves a ton of memory and time when doing LFQ, since we
    /// iterate over each peak.
    pub fn centroid_frame(
        &mut self,
        frame: &timsrust::Frame,
    ) -> (
        ClusteringSummary,
        impl Iterator<Item = CorrectedFramePeak> + '_,
    ) {
        let start = std::time::Instant::now();
        self.with_frame(frame);
        let (stop_cause, num_accumulated, window_capped) = self.itercentroid_frame();
        let elapsed = start.elapsed();
        let summary = ClusteringSummary {
            initial_peaks: self.peaks.len(),
            aggregated_peaks: num_accumulated,
            final_peaks: self.agg_buff.len(),
            window_capped,
            stopping_reason: stop_cause,
            elapsed,
        };
        (summary, self.drain_aggregated_peaks())
    }

    #[inline(always)]
    fn advance_pointers(&self, left_ptr: &mut usize, right_ptr: &mut usize, tof_index: u32) {
        let (left_tof, right_tof) = self.tof_index_bounds(tof_index);

        // Advance left pointer to first peak in range
        while *left_ptr < self.peaks.len() && self.peaks[*left_ptr].tof_index < left_tof {
            *left_ptr += 1;
        }

        // Advance right pointer to last peak in range
        // Start from max(right_ptr, left_ptr) to avoid going backwards
        *right_ptr = (*right_ptr).max(*left_ptr);
        while *right_ptr < self.peaks.len() && self.peaks[*right_ptr].tof_index <= right_tof {
            *right_ptr += 1;
        }
    }

    /// Pre-computes the neighbor ranges for each peak in the frame.
    /// This is used to speed up the centroiding process.
    /// It should be called after the peaks have been sorted by TOF index.
    /// and before the centroiding process.
    ///
    /// Stores them in the `neighbor_ranges` field.
    /// ... in otherwords ... for peak i, the neighbors are in
    /// `self.peaks[self.neighbor_ranges[i].0 .. self.neighbor_ranges[i].1]`
    /// also stores the total intensity of the neighbors
    /// in the `order_intensity` field.
    /// This is used to determine the order in which
    /// peaks are processed during centroiding.
    /// ... this seems to be critical to solve ties because a lot of peaks have the
    /// same intensity.
    fn compute_neighbor_ranges_and_intensity(&mut self) {
        self.neighbor_ranges.clear();
        self.order_intensity.clear();
        self.neighbor_ranges.reserve(self.peaks.len());
        self.order_intensity.reserve(self.peaks.len());

        let mut left_ptr = 0;
        let mut right_ptr = 0;

        let mut last_tof = u32::MAX; // impossible TOF index
        // ALSO we are startting with sorted tof indices, so the first should be a really low value
        // This just assures we compute the bounds for the first peak

        for idx in 0..self.peaks.len() {
            let peak = &self.peaks[idx];
            // If the last tof == current tof, we can reuse the left and right pointers
            if last_tof != peak.tof_index {
                // If the TOF has changed, we need to compute the new bounds
                self.advance_pointers(&mut left_ptr, &mut right_ptr, peak.tof_index);
                last_tof = peak.tof_index;
            }
            self.neighbor_ranges.push((left_ptr, right_ptr));

            let (left_im, right_im) = self.im_index_bounds(peak.scan_index);
            // Compute intensity directly here
            let mut summed_int = 0.0;
            for i in left_ptr..right_ptr {
                let scan_idx = self.peaks[i].scan_index;
                if scan_idx >= left_im && scan_idx <= right_im {
                    summed_int += self.peaks[i].corrected_intensity;
                }
            }
            self.order_intensity.push(summed_int);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// tof index is m/z, scan index is mobility.
    struct Identity;
    impl ConvertableDomain for Identity {
        fn convert<T: Into<f64> + Copy>(&self, v: T) -> f64 {
            v.into()
        }

        fn invert<T: Into<f64> + Copy>(&self, v: T) -> f64 {
            v.into()
        }
    }

    /// One scan holding `(tof, intensity)` peaks.
    fn frame(peaks: &[(u32, u32)]) -> timsrust::Frame {
        let mut fp = timsrust::FramePeaks::with_capacity(1, peaks.len());
        fp.scan_offsets = vec![0, peaks.len() as u32];
        for &(tof, i) in peaks {
            fp.tof_indices.push(tof);
            fp.intensities.push(i);
        }
        let meta = timsrust::FrameMeta {
            intensity_correction_factor: 1.0,
            ..Default::default()
        };
        timsrust::Frame { peaks: fp, meta }
    }

    fn config(window_cap: Option<WindowCap>) -> CentroidingConfig {
        CentroidingConfig {
            max_peaks: 1_000,
            mz_tol: MzTolerance::Ppm(5.0),
            im_tol: ImTolerance::Pct(3.0),
            // Every peak here is a singleton, which is exactly what early-stop
            // would clip.
            early_stop_iterations: 0,
            window_cap,
            transitive: true,
        }
    }

    /// Fifty peaks in 100-200, intensity rising with m/z; five in 300-400.
    fn two_window_frame() -> timsrust::Frame {
        let mut peaks: Vec<(u32, u32)> = (0..50u32).map(|k| (100 + 2 * k, 1_000 + k)).collect();
        peaks.extend((0..5u32).map(|k| (300 + 10 * k, 500)));
        frame(&peaks)
    }

    fn centroid(cfg: CentroidingConfig, f: &timsrust::Frame) -> (ClusteringSummary, Vec<u32>) {
        let mut c = PeakCentroider::with_capacity(64, cfg, Identity, Identity);
        let (summary, peaks) = c.centroid_frame(f);
        let mut tofs: Vec<u32> = peaks.map(|p| p.tof_index).collect();
        tofs.sort_unstable();
        (summary, tofs)
    }

    #[test]
    fn no_window_cap_keeps_every_singleton() {
        let (summary, tofs) = centroid(config(None), &two_window_frame());
        assert_eq!(tofs.len(), 55);
        assert_eq!(summary.window_capped, 0);
        assert_eq!(summary.stopping_reason, StoppingReason::AllTaken);
    }

    #[test]
    fn window_cap_keeps_the_most_intense_per_window_and_spares_sparse_windows() {
        let cap = WindowCap {
            max_peaks: 10,
            window_da: 100.0,
        };
        let (summary, tofs) = centroid(config(Some(cap)), &two_window_frame());

        let dense: Vec<u32> = tofs.iter().copied().filter(|t| *t < 200).collect();
        let sparse: Vec<u32> = tofs.iter().copied().filter(|t| *t >= 300).collect();
        // The ten most intense in 100-200 are k = 40..50, at tof 180..=198.
        assert_eq!(dense, (40..50u32).map(|k| 100 + 2 * k).collect::<Vec<_>>());
        // The sparse window is under quota and untouched.
        assert_eq!(sparse, vec![300, 310, 320, 330, 340]);
        assert_eq!(summary.window_capped, 40);
        assert_eq!(summary.final_peaks, 15);
        // Dropped peaks still count as accounted for, so the frame finishes.
        assert_eq!(summary.stopping_reason, StoppingReason::AllTaken);
    }

    #[test]
    fn window_counts_reset_between_frames() {
        let cap = WindowCap {
            max_peaks: 10,
            window_da: 100.0,
        };
        let mut c = PeakCentroider::with_capacity(64, config(Some(cap)), Identity, Identity);
        let f = two_window_frame();
        let (_, first) = c.centroid_frame(&f);
        let n_first = first.count();
        let (_, second) = c.centroid_frame(&f);
        assert_eq!(second.count(), n_first, "a full window must not stay full");
    }

    #[test]
    fn bins_tolerance_merges_adjacent_tof_indices_and_ppm_below_a_bin_does_not() {
        // The same ion quantized to adjacent bins on two scans. With the
        // identity converter one bin is one m/z unit, so at m/z 1000 a bin is
        // 1000 ppm and a 1 ppm tolerance rounds to zero bins.
        let f = frame(&[(1000, 100), (1001, 90)]);
        let ppm = CentroidingConfig {
            mz_tol: MzTolerance::Ppm(1.0),
            ..config(None)
        };
        let (_, tofs) = centroid(ppm, &f);
        assert_eq!(tofs.len(), 2, "sub-bin ppm keeps the two apart");

        let bins = CentroidingConfig {
            mz_tol: MzTolerance::Bins(1),
            ..config(None)
        };
        let (_, tofs) = centroid(bins, &f);
        assert_eq!(tofs.len(), 1, "one bin of tolerance merges them");
        // Intensity-weighted mean, truncated to the bin: 1000.47 -> 1000.
        assert_eq!(tofs, vec![1000]);
    }

    #[test]
    fn scans_tolerance_is_a_plain_scan_index_radius() {
        // Same TOF index, scans 10 and 13. A radius of 2 keeps them apart, 3
        // merges them. `frame` puts everything on one scan, so build by hand.
        let mut fp = timsrust::FramePeaks::with_capacity(14, 2);
        // scan_offsets: 14 scans; peaks at scan 10 and scan 13.
        fp.scan_offsets = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2];
        fp.tof_indices = vec![500, 500];
        fp.intensities = vec![100, 50];
        let meta = timsrust::FrameMeta {
            intensity_correction_factor: 1.0,
            ..Default::default()
        };
        let f = timsrust::Frame { peaks: fp, meta };

        for (radius, expect) in [(2u16, 2usize), (3, 1)] {
            let cfg = CentroidingConfig {
                im_tol: ImTolerance::Scans(radius),
                ..config(None)
            };
            let (_, tofs) = centroid(cfg, &f);
            assert_eq!(tofs.len(), expect, "radius {radius}");
        }
    }

    #[test]
    fn transitive_growth_chains_bin_to_bin_and_non_transitive_does_not() {
        // A run of five adjacent bins. Parents are picked by neighborhood
        // intensity, so 1001 (100+90+80) leads and takes 1000 and 1002.
        // Transitive with a 1-bin tolerance: 1002 later recruits 1003, which
        // recruits 1004 -- one cluster spanning five bins. Non-transitive:
        // 1003 becomes its own parent and takes 1004 -- two clusters, each no
        // wider than the tolerance around its parent.
        let f = frame(&[(1000, 100), (1001, 90), (1002, 80), (1003, 70), (1004, 60)]);
        let base = CentroidingConfig {
            mz_tol: MzTolerance::Bins(1),
            ..config(None)
        };
        let (_, chained) = centroid(base, &f);
        assert_eq!(chained.len(), 1, "transitive: the whole run is one cluster");

        let bounded = CentroidingConfig {
            transitive: false,
            ..base
        };
        let (_, tofs) = centroid(bounded, &f);
        assert_eq!(tofs.len(), 2, "non-transitive: radius is the tolerance");
    }

    #[test]
    fn a_dropped_peak_is_not_absorbed_by_a_weaker_neighbor() {
        // Two peaks within ppm tolerance in a window of quota 1, plus a third
        // strong peak so the quota is taken before either is visited. The
        // strong peak at 150 claims the slot; 100 is dropped; 100.0003 must
        // not become a parent with 100 folded into it.
        let cap = WindowCap {
            max_peaks: 1,
            window_da: 100.0,
        };
        // tof units are m/z here; 100 and 100 + 0.0003 are 3 ppm apart, but
        // tof indices are integers, so scale the whole frame by 1e4.
        let f = frame(&[(1_500_000, 900), (1_000_000, 500), (1_000_003, 100)]);
        let cfg = CentroidingConfig {
            window_cap: Some(WindowCap {
                window_da: 1_000_000.0,
                ..cap
            }),
            ..config(None)
        };
        let (summary, tofs) = centroid(cfg, &f);
        assert_eq!(tofs, vec![1_500_000]);
        assert_eq!(summary.window_capped, 2);
    }
}
