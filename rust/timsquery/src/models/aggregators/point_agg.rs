use crate::models::elution_group::ElutionGroup;
use crate::KeyLike;
use serde::Serialize;

#[derive(Debug, Clone, Copy)]
pub struct RawPeakIntensityAggregator<'a, T: KeyLike> {
    pub query: &'a ElutionGroup<T>,
    pub intensity: u64,
}

impl <'a, T: KeyLike>RawPeakIntensityAggregator<'_, T> {
    pub fn new_with_elution_group(
        elution_group: &ElutionGroup<T>,
    ) -> Self {
        Self {
            query: elution_group,
            intensity: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RawPeakVectorAggregator<'a,T: KeyLike> {
    pub query: &'a ElutionGroup<T>,
    pub peaks: RawPeakVectorArrays,
}

impl <'a, T: KeyLike>RawPeakVectorAggregator<'_ , T> {
    pub fn new_with_elution_group(
        elution_group: &ElutionGroup<T>,
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

