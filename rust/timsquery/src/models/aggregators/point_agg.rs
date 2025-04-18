use crate::models::elution_group::ElutionGroup;
use crate::KeyLike;
use serde::Serialize;
use std::sync::Arc;

#[derive(Debug, Clone, Serialize)]
pub struct PointIntensityAggregator<T: KeyLike> {
    pub query: Arc<ElutionGroup<T>>,
    pub intensity: u64,
}

impl <T: KeyLike>PointIntensityAggregator<T> {
    pub fn new_with_elution_group(
        elution_group: Arc<ElutionGroup<T>>,
    ) -> Self {
        Self {
            query: elution_group,
            intensity: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RawPeakVectorAggregator<T: KeyLike> {
    pub query: Arc<ElutionGroup<T>>,
    pub peaks: RawPeakVectorArrays,
}

impl <T: KeyLike>RawPeakVectorAggregator<T> {
    pub fn new_with_elution_group(
        elution_group: Arc<ElutionGroup<T>>,
    ) -> Self {
        Self {
            query: elution_group,
            peaks: RawPeakVectorArrays::new(),
        }
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

