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
    mz_ppm_tol: f64,
    im_pct_tol: f64,
    mz_converter: T1,
    im_converter: T2,
    early_stop_iterations: u32,
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
    /// M/Z tolerance in ppm for clustering, 5.0 (5 ppm)
    /// seems to work well in practice.
    pub mz_ppm_tol: f64,
    /// IM tolerance in percentage for clustering
    /// 3.0 (3%) seems to work well in practice.
    pub im_pct_tol: f64,
    /// Number of consecutive iterations with no new peaks clustered
    /// after which the centroiding will stop early (instead of going
    /// through all peaks, which will very likely be noise).
    /// A number ~200 seems to work well in practice.
    pub early_stop_iterations: u32,
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
        self.weighed_tof_sum += peak.tof_index as u64 * weight;
        self.weighed_im_sum += peak.scan_index as u64 * weight;
        self.total_weight += weight;
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
    Taken { parent_idx: usize },
    Parent,
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
    pub stopping_reason: StoppingReason,
    pub elapsed: std::time::Duration,
}

/// Summary of the clustering over multiple frames
#[derive(Debug, Clone)]
pub struct AggregatedClusteringSummary {
    pub total_initial_peaks: usize,
    pub total_aggregated_peaks: usize,
    pub total_final_peaks: usize,
    pub frames_processed: usize,
    pub stopping_reasons: [usize; 3],
    pub elapsed: std::time::Duration,
}

impl AggregatedClusteringSummary {
    pub(crate) fn combine(mut self, right: &Self) -> Self {
        self.total_final_peaks += right.total_final_peaks;
        self.total_aggregated_peaks += right.total_aggregated_peaks;
        self.total_initial_peaks += right.total_initial_peaks;
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
            mz_ppm_tol: config.mz_ppm_tol,
            im_pct_tol: config.im_pct_tol,
            mz_converter,
            im_converter,
            early_stop_iterations: config.early_stop_iterations,
        }
    }

    /// Given a TOF index, returns the (inclusive) bounds of TOF indices
    /// that fall within the ppm tolerance of the mz corresponding to the TOF index.
    /// This is used to find neighboring peaks for clustering.
    fn tof_index_bounds(&self, tof_idex: u32) -> (u32, u32) {
        let mz = self.mz_converter.convert(tof_idex as f64);
        let delta_mz = mz * self.mz_ppm_tol * 1e-6;
        let left_mz = mz - delta_mz;
        let right_mz = mz + delta_mz;
        let left_tof = self.mz_converter.invert(left_mz).round() as u32;
        let right_tof = self.mz_converter.invert(right_mz).round() as u32;
        (left_tof, right_tof)
    }

    fn im_index_bounds(&self, im_index: u16) -> (u16, u16) {
        self.ims_ranges[im_index as usize]
    }

    fn uncached_im_index_bounds(&self, im_index: u16) -> (u16, u16) {
        let im = self.im_converter.convert(im_index as f64);
        let delta_im = im * self.im_pct_tol * 0.01;
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
        self.peaks
            .sort_unstable_by(|a, b| a.tof_index.cmp(&b.tof_index));
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
    fn itercentroid_frame(&mut self) -> (StoppingReason, usize) {
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
        let mut early_stop_remaining = self.early_stop_iterations;
        let mut out = StoppingReason::AllTaken;

        for &idx in &self.order {
            if self.agg_buff.len() > self.max_peaks {
                out = StoppingReason::MaxPeaks;
                break;
            }

            let mut num_includable = 0;

            let (target_index, is_self_parent) = match self.taken_buff[idx] {
                TakenState::Parent => unreachable!("This should never happen"),
                TakenState::Taken { parent_idx } => (parent_idx, false),
                TakenState::Untaken => {
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
                    TakenState::Taken { .. } => continue,
                    TakenState::Parent => continue,
                    _ => { /* continue processing */ }
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

            if is_self_parent {
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
        (out, global_num_included)
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
        let (stop_cause, num_accumulated) = self.itercentroid_frame();
        let elapsed = start.elapsed();
        let summary = ClusteringSummary {
            initial_peaks: self.peaks.len(),
            aggregated_peaks: num_accumulated,
            final_peaks: self.agg_buff.len(),
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

    #[test]
    fn it_works() {}
}
