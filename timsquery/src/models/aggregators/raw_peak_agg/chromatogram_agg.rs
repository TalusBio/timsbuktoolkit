use super::super::streaming_aggregator::RunningStatsCalculator;
use crate::sort_vecs_by_first;

use nohash_hasher::BuildNoHashHasher;
use serde::Serialize;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;

/// A struct that can be used to calculate the mean and variance
/// of a stream of weighted tof and scan numbers.
#[derive(Debug, Clone)]
pub struct ScanTofStatsCalculatorPair {
    pub scan: RunningStatsCalculator,
    pub tof: RunningStatsCalculator,
}

impl ScanTofStatsCalculatorPair {
    pub fn new(intensity: u64, scan_index: usize, tof_index: u32) -> Self {
        let scan_index = scan_index as f64;
        let tof_index = tof_index as f64;

        Self {
            scan: RunningStatsCalculator::new(intensity, scan_index),
            tof: RunningStatsCalculator::new(intensity, tof_index),
        }
    }

    pub fn add(&mut self, intensity: u64, scan_index: usize, tof_index: u32) {
        self.scan.add(scan_index as f64, intensity);
        self.tof.add(tof_index as f64, intensity);
    }

    pub fn weight(&self) -> u64 {
        self.scan.weight()
    }
}

#[derive(Debug, Clone)]
pub struct ChromatomobilogramStats {
    // TODO OPTIMIZE THIS ... as needed.
    // In theory we can optimize this to make a single aggregator struct
    // that shares the weight (intensity), since all will have the same weight
    // and retention times.
    pub scan_tof_mapping: HashMap<u32, ScanTofStatsCalculatorPair, BuildNoHashHasher<u32>>,
    pub id: u64,
    pub expected_tof_index: u32,
    pub expected_scan_index: usize,
}

impl ChromatomobilogramStats {
    pub fn new(id: u64, expected_tof_index: u32, expected_scan_index: usize) -> Self {
        Self {
            scan_tof_mapping: HashMap::with_hasher(BuildHasherDefault::default()),
            id,
            expected_tof_index,
            expected_scan_index,
        }
    }
}

/// A struct that calculates stores the mean/sd of an isolation window (transition)
/// across the chromatogram.
#[derive(Debug, Clone, Serialize, Default)]
pub struct ChromatomobilogramStatsArrays {
    pub retention_time_miliseconds: Vec<u32>,
    pub tof_index_means: Vec<f64>,
    pub tof_index_sds: Vec<f64>,
    pub scan_index_means: Vec<f64>,
    pub scan_index_sds: Vec<f64>,
    pub intensities: Vec<u64>,
    pub expected_intensity: Option<f64>,
    pub expected_tof_index: u32,
    pub expected_scan_index: usize,
}

impl ChromatomobilogramStatsArrays {
    // TODO use default instead of new everywhere ..
    pub fn new() -> Self {
        Self::default()
    }

    pub fn sort_by_rt(&mut self) {
        let x = sort_vecs_by_first!(
            &mut self.retention_time_miliseconds,
            &mut self.tof_index_means,
            &mut self.tof_index_sds,
            &mut self.scan_index_means,
            &mut self.scan_index_sds,
            &mut self.intensities
        );
        self.retention_time_miliseconds = x.0;
        self.tof_index_means = x.1;
        self.tof_index_sds = x.2;
        self.scan_index_means = x.3;
        self.scan_index_sds = x.4;
        self.intensities = x.5;
    }

    pub fn is_sorted(&self) -> bool {
        self.retention_time_miliseconds
            .windows(2)
            .all(|x| x[0] <= x[1])
    }

    pub fn len(&self) -> usize {
        self.retention_time_miliseconds.len()
    }

    pub fn is_empty(&self) -> bool {
        self.retention_time_miliseconds.is_empty()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromatomobilogramVectorArrayTuples {
    pub scan_indices: Vec<(usize, u64)>,
    pub tof_indices: Vec<(u32, u64)>,
    pub retention_times: Vec<(f32, u64)>,
}
