use crate::KeyLike;
use crate::models::elution_group::TimsElutionGroup;
use crate::traits::queriable_data::HasQueryData;
use serde::Serialize;
use std::sync::Arc;
use tinyvec::TinyVec;

const POINT_INLINE_CAP: usize = 13;

#[derive(Debug, Clone, Serialize)]
pub struct PointIntensityAggregator<T: KeyLike> {
    pub id: u64,
    pub mobility_ook0: f32,
    pub rt_seconds: f32,
    pub precursor_mono_mz: f64,
    pub precursor_charge: u8,
    pub precursor_mz_limits: (f64, f64),
    pub precursor_labels: TinyVec<[i8; POINT_INLINE_CAP]>,
    pub precursor_mzs: TinyVec<[f64; POINT_INLINE_CAP]>,
    pub fragment_labels: TinyVec<[T; POINT_INLINE_CAP]>,
    pub fragment_mzs: TinyVec<[f64; POINT_INLINE_CAP]>,
    pub intensity: f64,
}

impl<T: KeyLike> PointIntensityAggregator<T> {
    pub fn new_with_elution_group(elution_group: Arc<TimsElutionGroup<T>>) -> Self {
        Self::new(&elution_group)
    }

    pub fn new(eg: &TimsElutionGroup<T>) -> Self {
        let mut precursor_labels = TinyVec::new();
        let mut precursor_mzs = TinyVec::new();
        for (lbl, mz) in eg.iter_precursors() {
            precursor_labels.push(lbl);
            precursor_mzs.push(mz);
        }
        let mut fragment_labels = TinyVec::new();
        let mut fragment_mzs = TinyVec::new();
        for (lbl, mz) in eg.iter_fragments_refs() {
            fragment_labels.push(lbl.clone());
            fragment_mzs.push(*mz);
        }
        Self {
            id: eg.id(),
            mobility_ook0: eg.mobility_ook0(),
            rt_seconds: eg.rt_seconds(),
            precursor_mono_mz: eg.mono_precursor_mz(),
            precursor_charge: eg.precursor_charge(),
            precursor_mz_limits: eg.get_precursor_mz_limits(),
            precursor_labels,
            precursor_mzs,
            fragment_labels,
            fragment_mzs,
            intensity: 0.0,
        }
    }
}

impl<T: KeyLike> HasQueryData<T> for PointIntensityAggregator<T> {
    fn id(&self) -> u64 {
        self.id
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.precursor_mz_limits
    }

    fn mobility_ook0(&self) -> f32 {
        self.mobility_ook0
    }

    fn rt_seconds(&self) -> f32 {
        self.rt_seconds
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> + '_ {
        self.precursor_labels
            .iter()
            .copied()
            .zip(self.precursor_mzs.iter().copied())
    }

    fn iter_fragments<'a>(&'a self) -> impl Iterator<Item = (&'a T, f64)> + 'a
    where
        T: 'a,
    {
        self.fragment_labels
            .iter()
            .zip(self.fragment_mzs.iter().copied())
    }
}

#[derive(Debug, Clone)]
pub struct RawPeakVectorAggregator<T: KeyLike> {
    pub query: Arc<TimsElutionGroup<T>>,
    pub peaks: RawPeakVectorArrays,
}

impl<T: KeyLike> RawPeakVectorAggregator<T> {
    pub fn new_with_elution_group(elution_group: Arc<TimsElutionGroup<T>>) -> Self {
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
