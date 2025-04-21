use crate::sort_vecs_by_first;
use crate::utils::display::{
    GlimpseConfig,
    glimpse_vec,
};
use crate::utils::tolerance_ranges::IncludedRange;
use std::fmt::Display;

pub struct PeakInBucket {
    pub scan_index: u16,
    pub corrected_intensity: f32,
    pub retention_time_ms: u32,
}

#[derive(Debug)]
pub struct PeakBucket {
    corrected_intensities: Vec<f32>,
    retention_times_ms: Vec<u32>,
    scan_offsets: Vec<u16>,
    tof_index: u32,
    // TODO add this ...
    // mz: f64,
}

impl Display for PeakBucket {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PeakBucket tof={}: \n    len={},\n    retention_times={},\n    scan_offsets={},\n    intensities={}",
            self.tof_index,
            self.len(),
            glimpse_vec(
                &self.retention_times_ms,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
            glimpse_vec(
                &self.scan_offsets,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
            glimpse_vec(
                &self.corrected_intensities,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
        )
    }
}

#[derive(Debug)]
pub struct PeakBucketBuilder {
    corrected_intensities: Vec<f32>,
    retention_times_ms: Vec<u32>,
    scan_offsets: Vec<u16>,
    tof_index: u32,
}

impl PeakBucketBuilder {
    pub fn new(capacity: usize, tof_index: u32) -> Self {
        Self {
            corrected_intensities: Vec::with_capacity(capacity),
            retention_times_ms: Vec::with_capacity(capacity),
            scan_offsets: Vec::with_capacity(capacity),
            tof_index,
        }
    }

    pub fn len(&self) -> usize {
        self.corrected_intensities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.corrected_intensities.is_empty()
    }

    pub fn add_peak(&mut self, scan_index: u16, corrected_intensity: f32, retention_time_ms: u32) {
        self.corrected_intensities.push(corrected_intensity);
        self.retention_times_ms.push(retention_time_ms);
        self.scan_offsets.push(scan_index);
    }

    pub fn extend_peaks(
        &mut self,
        scan_indices: &[u16],
        corrected_intensities: &[f32],
        retention_time_ms: &[u32],
    ) {
        self.scan_offsets.extend(scan_indices);
        self.corrected_intensities.extend(corrected_intensities);
        self.retention_times_ms.extend(retention_time_ms);
    }

    pub fn build(self) -> PeakBucket {
        let sorted = sort_vecs_by_first!(
            &self.scan_offsets,
            &self.retention_times_ms,
            &self.corrected_intensities
        );
        let scan_offsets = sorted.0;
        let retention_times_ms = sorted.1;
        let corrected_intensities = sorted.2;

        let out = PeakBucket {
            corrected_intensities,
            retention_times_ms,
            scan_offsets,
            tof_index: self.tof_index,
        };

        debug_assert!(out.verify(), "PeakBucket::build failed at verify");
        out
    }
}

impl PeakBucket {
    pub fn len(&self) -> usize {
        self.corrected_intensities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.corrected_intensities.is_empty()
    }

    pub fn query_peaks(
        &self,
        scan_range: Option<IncludedRange<u16>>,
        rt_range_ms: Option<IncludedRange<u32>>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let (scan_min, scan_max) = match scan_range {
            Some(x) => x.into(),
            None => (
                0,
                *(self
                    .scan_offsets
                    .last()
                    .expect("It should not be possible to build an empty PeakBucket")),
            ),
        };

        // TODO do a binary search if the data is large-ish... maybe ...
        let mut start_min = 0;
        let mut end_max = self.len();
        while start_min < end_max && self.scan_offsets[start_min] < scan_min {
            start_min += 1;
        }
        while start_min < end_max && self.scan_offsets[end_max - 1] >= scan_max {
            end_max -= 1;
        }

        (start_min..end_max).filter_map(move |x| {
            let scan_index = self.scan_offsets[x];
            if scan_index < scan_min || scan_index > scan_max {
                // TODO search what cases make this branch happen ...
                return None;
            }

            let retention_time_ms = self.retention_times_ms[x];
            if let Some(x) = rt_range_ms {
                if !x.contains(retention_time_ms) {
                    return None;
                }
            }
            Some(PeakInBucket {
                scan_index,
                corrected_intensity: self.corrected_intensities[x],
                retention_time_ms,
            })
        })
    }

    fn verify(&self) -> bool {
        if self.corrected_intensities.len() != self.retention_times_ms.len() {
            println!("PeakBucket::verify failed at length check");
            return false;
        }
        if self.corrected_intensities.len() != self.scan_offsets.len() {
            println!("PeakBucket::verify failed at length check on sorted mode");
            return false;
        }
        for i in 1..self.scan_offsets.len() {
            if self.scan_offsets[i] < self.scan_offsets[i - 1] {
                println!("PeakBucket::verify failed at scan order check");
                return false;
            }
        }

        true
    }
}
