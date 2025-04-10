use crate::models::elution_group;
use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::{
    Aggregator,
    NoContext,
};
use serde::Serialize;

#[derive(Debug, Clone, Copy)]
pub struct RawPeakIntensityAggregator {
    pub id: u64,
    pub intensity: u64,
}

impl RawPeakIntensityAggregator {
    pub fn new(id: u64) -> Self {
        Self { id, intensity: 0 }
    }

    pub fn new_with_elution_group<
        T: Send + Sync + std::hash::Hash + Clone + Eq + PartialEq + Serialize + std::fmt::Debug,
    >(
        elution_group: &elution_group::ElutionGroup<T>,
    ) -> Self {
        Self {
            id: elution_group.id,
            intensity: 0,
        }
    }
}

impl Aggregator for RawPeakIntensityAggregator {
    type Context = NoContext;
    type Item = RawPeak;
    type Output = u64;

    fn add(&mut self, peak: impl Into<RawPeak>) {
        let peak = peak.into();
        self.intensity += peak.intensity as u64;
    }

    fn finalize(self) -> u64 {
        self.intensity
    }

    fn supports_context(&self) -> bool {
        false
    }
}

#[derive(Debug, Clone)]
pub struct RawPeakVectorAggregator {
    pub id: u64,
    pub peaks: RawPeakVectorArrays,
}

impl RawPeakVectorAggregator {
    pub fn new(id: u64) -> Self {
        Self {
            id,
            peaks: RawPeakVectorArrays::new(),
        }
    }

    pub fn new_with_elution_group<
        T: Send + Sync + std::hash::Hash + Clone + Eq + PartialEq + Serialize + std::fmt::Debug,
    >(
        elution_group: &elution_group::ElutionGroup<T>,
    ) -> Self {
        Self::new(elution_group.id)
    }
}

impl RawPeakVectorArrays {
    pub fn new() -> Self {
        Self {
            scans: Vec::new(),
            tofs: Vec::new(),
            intensities: Vec::new(),
            retention_times: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct RawPeakVectorArrays {
    pub scans: Vec<usize>,
    pub tofs: Vec<u32>,
    pub intensities: Vec<u32>,
    pub retention_times: Vec<f32>,
}

impl Aggregator for RawPeakVectorAggregator {
    type Context = NoContext;
    type Item = RawPeak;
    type Output = RawPeakVectorArrays;

    fn add(&mut self, peak: impl Into<RawPeak>) {
        let peak = peak.into();
        self.peaks.scans.push(peak.scan_index);
        self.peaks.tofs.push(peak.tof_index);
        self.peaks.intensities.push(peak.intensity);
        self.peaks.retention_times.push(peak.retention_time);
    }

    fn finalize(self) -> RawPeakVectorArrays {
        self.peaks
    }

    fn supports_context(&self) -> bool {
        false
    }
}
