use crate::models::target::OwnedTarget;
use crate::traits::Target;
use crate::traits::queriable_data::HasQueryData;
use crate::{
    ExtractionQuery,
    KeyLike,
    ResolvedRt,
};
use serde::Serialize;
use std::sync::Arc;
use tinyvec::TinyVec;

const POINT_INLINE_CAP: usize = 13;

fn serialize_rt_center<S: serde::Serializer>(
    rt: &ResolvedRt,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    rt.center().map_or(f32::NAN, |r| r.0).serialize(serializer)
}

#[derive(Debug, Clone, Serialize)]
pub struct PointIntensityAggregator<T: KeyLike> {
    pub id: crate::models::OwnedSourceId,
    pub mobility_ook0: f32,
    #[serde(rename = "rt_seconds", serialize_with = "serialize_rt_center")]
    pub rt: ResolvedRt,
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
    pub fn new(query: &ExtractionQuery<'_, impl Target<Label = T>>) -> Self {
        let eg = query.source();
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
            fragment_mzs.push(mz);
        }
        Self {
            id: eg.output_id().to_owned_id(),
            mobility_ook0: query.mobility_center(),
            rt: query.rt(),
            precursor_mono_mz: eg.mono_precursor_mz(),
            precursor_charge: eg.precursor_charge(),
            precursor_mz_limits: eg.precursor_mz_limits(),
            precursor_labels,
            precursor_mzs,
            fragment_labels,
            fragment_mzs,
            intensity: 0.0,
        }
    }
}

impl<T: KeyLike> HasQueryData<T> for PointIntensityAggregator<T> {
    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.precursor_mz_limits
    }

    fn mobility_ook0(&self) -> f32 {
        self.mobility_ook0
    }

    fn rt(&self) -> ResolvedRt {
        self.rt
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
    pub query: Arc<OwnedTarget<T>>,
    pub peaks: RawPeakVectorArrays,
}

impl<T: KeyLike> RawPeakVectorAggregator<T> {
    pub fn new_with_elution_group(elution_group: Arc<OwnedTarget<T>>) -> Self {
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
