use crate::sort_vecs_by_first;
use crate::utils::compress_explode::compress_vec;
use crate::utils::display::{
    GlimpseConfig,
    glimpse_vec,
};
use crate::utils::tolerance_ranges::IncludedRange;
use std::fmt::Display;

pub struct PeakInBucket {
    pub scan_index: usize,
    pub intensity: u32,
    pub retention_time: f32,
}

#[derive(Debug, Clone, Copy)]
pub enum PeakBucketMode {
    Compressed,
    Sorted,
}

#[derive(Debug)]
pub struct PeakBucket {
    intensities: Vec<u32>,
    retention_times: Vec<f32>,
    scan_offsets: Vec<usize>,
    tof_index: u32,
    mode: PeakBucketMode,
}

impl Display for PeakBucket {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PeakBucket tof={}: \n    len={},\n    retention_times={},\n    scan_offsets={},\n    intensities={}, \n    mode={:?}",
            self.tof_index,
            self.len(),
            glimpse_vec(
                &self.retention_times,
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
                &self.intensities,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
            self.mode,
        )
    }
}

#[derive(Debug)]
pub struct PeakBucketBuilder {
    intensities: Vec<u32>,
    retention_times: Vec<f32>,
    scan_offsets: Vec<usize>,
    tof_index: u32,
}

impl PeakBucketBuilder {
    pub fn new(capacity: usize, tof_index: u32) -> Self {
        Self {
            intensities: Vec::with_capacity(capacity),
            retention_times: Vec::with_capacity(capacity),
            scan_offsets: Vec::with_capacity(capacity),
            tof_index,
        }
    }

    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intensities.is_empty()
    }

    pub fn add_peak(&mut self, scan_index: usize, intensity: u32, retention_time: f32) {
        self.intensities.push(intensity);
        self.retention_times.push(retention_time);
        self.scan_offsets.push(scan_index);
    }

    pub fn extend_peaks(
        &mut self,
        scan_indices: &[usize],
        intensities: &[u32],
        retention_times: &[f32],
    ) {
        self.scan_offsets.extend(scan_indices);
        self.intensities.extend(intensities);
        self.retention_times.extend(retention_times);
    }

    pub fn build(mut self) -> PeakBucket {
        let sorted =
            sort_vecs_by_first!(&self.scan_offsets, &self.retention_times, &self.intensities);
        self.scan_offsets = sorted.0;
        self.retention_times = sorted.1;
        self.intensities = sorted.2;

        // TODO consider if I really need to compress this.
        // Options:
        // 1. Change the compression of scans (I use the default in timsrust but im sure we can do
        //    better) ...
        //    - tuple of val-offset to prevent the long runs of non changing scans.
        //    - btree of offsets?
        // 2. Not compressing makes queries <10% slower, but building the index 10% faster ...
        //    compressing all makes building the index 3x slower.
        //
        // let out = if false {
        let out = if self.scan_offsets.len() > 500 {
            let compressed = compress_vec(&self.scan_offsets);
            PeakBucket {
                intensities: self.intensities,
                retention_times: self.retention_times,
                scan_offsets: compressed,
                mode: PeakBucketMode::Compressed,
                tof_index: self.tof_index,
            }
        } else {
            PeakBucket {
                intensities: self.intensities,
                retention_times: self.retention_times,
                scan_offsets: self.scan_offsets,
                mode: PeakBucketMode::Sorted,
                tof_index: self.tof_index,
            }
        };

        debug_assert!(out.verify(), "PeakBucket::build failed at verify");
        out
    }
}

impl PeakBucket {
    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intensities.is_empty()
    }

    pub fn query_peaks(
        &self,
        scan_range: Option<IncludedRange<usize>>,
        rt_range: Option<IncludedRange<f32>>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        match self.mode {
            PeakBucketMode::Compressed => self
                .query_peaks_compressed(scan_range, rt_range)
                .collect::<Vec<PeakInBucket>>()
                .into_iter(),
            PeakBucketMode::Sorted => self
                .query_peaks_sorted(scan_range, rt_range)
                .collect::<Vec<PeakInBucket>>()
                .into_iter(),
        }
    }

    pub fn query_peaks_sorted(
        &self,
        scan_range: Option<IncludedRange<usize>>,
        rt_range: Option<IncludedRange<f32>>,
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

            let retention_time = self.retention_times[x];
            if let Some(x) = rt_range {
                if !x.contains(retention_time) {
                    return None;
                }
            }
            Some(PeakInBucket {
                scan_index,
                intensity: self.intensities[x],
                retention_time,
            })
        })
    }

    pub fn query_peaks_compressed(
        &self,
        scan_range: Option<IncludedRange<usize>>,
        rt_range: Option<IncludedRange<f32>>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let scan_range = match scan_range {
            Some(x) => x.start()..x.end().min(self.scan_offsets.len() - 1),
            None => 0..(self.scan_offsets.len() - 1),
        };
        scan_range
            .flat_map(move |scan_index| {
                let peak_index_start = self.scan_offsets[scan_index];
                let peak_index_end = self.scan_offsets[scan_index + 1];

                (peak_index_start..peak_index_end).map(move |peak_index| {
                    let retention_time = self.retention_times[peak_index];

                    if let Some(x) = rt_range {
                        if !x.contains(retention_time) {
                            return None;
                        };
                    }

                    let intensity = self.intensities[peak_index];
                    Some(PeakInBucket {
                        scan_index,
                        intensity,
                        retention_time,
                    })
                })
            })
            .flatten()
    }

    fn verify(&self) -> bool {
        if self.intensities.len() != self.retention_times.len() {
            println!("PeakBucket::verify failed at length check");
            return false;
        }
        match self.mode {
            PeakBucketMode::Compressed => {
                match self.scan_offsets.last() {
                    Some(last) => {
                        if *last != self.retention_times.len() {
                            println!(
                                "PeakBucket::verify failed at last scan check, in compressed mode."
                            );
                            return false;
                        }
                    }
                    None => {
                        println!("PeakBucket::verify failed at last scan check, no scans");
                        return false;
                    }
                }

                for i in 1..self.scan_offsets.len() {
                    if self.scan_offsets[i] < self.scan_offsets[i - 1] {
                        println!("PeakBucket::verify failed at scan order check");
                        return false;
                    }
                }
            }
            PeakBucketMode::Sorted => {
                if self.intensities.len() != self.scan_offsets.len() {
                    println!("PeakBucket::verify failed at length check on sorted mode");
                    return false;
                }
                for i in 1..self.scan_offsets.len() {
                    if self.scan_offsets[i] < self.scan_offsets[i - 1] {
                        println!("PeakBucket::verify failed at scan order check");
                        return false;
                    }
                }
            }
        };

        true
    }
}
