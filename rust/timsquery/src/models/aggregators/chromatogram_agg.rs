use serde::Serialize;

use crate::errors::DataProcessingError;
use crate::models::base::{
    MutableChromatogram,
    MzMajorIntensityArray,
};
use crate::utils::tolerance_ranges::IncludedRange;
use crate::{
    ElutionGroup,
    KeyLike,
};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize)]
pub struct ChromatogramCollector<T: KeyLike> {
    pub eg: Arc<ElutionGroup<T>>,
    pub precursors: MzMajorIntensityArray<i8, f32>, // TODO: Move this to a generic ...
    pub fragments: MzMajorIntensityArray<T, f32>,
}

impl<T: KeyLike> ChromatogramCollector<T> {
    pub fn new(
        eg: Arc<ElutionGroup<T>>,
        ref_rt_mss: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        let precursors =
            MzMajorIntensityArray::try_new_empty(eg.precursors.clone(), ref_rt_mss.clone())?;
        let fragments = MzMajorIntensityArray::try_new_empty(eg.fragments.clone(), ref_rt_mss)?;
        Ok(Self {
            eg,
            precursors,
            fragments,
        })
    }

    pub fn iter_mut_precursors(
        &mut self,
    ) -> impl Iterator<Item = (&(i8, f64), MutableChromatogram<f32>)> {
        self.precursors.iter_mut_mzs()
    }

    pub fn iter_mut_fragments(
        &mut self,
    ) -> impl Iterator<Item = (&(T, f64), MutableChromatogram<f32>)> {
        self.fragments.iter_mut_mzs()
    }

    pub fn unpack(
        self,
    ) -> (
        Arc<ElutionGroup<T>>,
        MzMajorIntensityArray<i8, f32>,
        MzMajorIntensityArray<T, f32>,
    ) {
        (self.eg, self.precursors, self.fragments)
    }

    pub fn rt_range(&self) -> IncludedRange<u32> {
        let min = self
            .fragments
            .rts_ms
            .first()
            .unwrap()
            .min(self.precursors.rts_ms.first().unwrap());
        let max = self
            .fragments
            .rts_ms
            .last()
            .unwrap()
            .min(self.precursors.rts_ms.last().unwrap());
        (*min, *max).into()
    }
}
