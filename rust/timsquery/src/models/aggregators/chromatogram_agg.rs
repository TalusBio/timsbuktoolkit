use serde::Serialize;

use crate::errors::DataProcessingError;
use crate::models::base::{
    ArrayElement,
    MutableChromatogram,
    MzMajorIntensityArray,
};
use crate::{
    ElutionGroup,
    KeyLike,
    ValueLike,
};
use std::sync::Arc;
use timscentroid::utils::TupleRange;

#[derive(Debug, Clone, Serialize)]
pub struct ChromatogramCollector<T: KeyLike, V: ArrayElement + ValueLike> {
    pub eg: ElutionGroup<T>,
    pub precursors: MzMajorIntensityArray<i8, V>,
    pub fragments: MzMajorIntensityArray<T, V>,
    pub ref_rt_ms: Arc<[u32]>,
}

impl<T: KeyLike, V: ValueLike + ArrayElement> ChromatogramCollector<T, V> {
    pub fn new(eg: ElutionGroup<T>, ref_rt_ms: Arc<[u32]>) -> Result<Self, DataProcessingError> {
        let precursors =
            MzMajorIntensityArray::try_new_empty(eg.precursors.clone(), ref_rt_ms.len(), 0)?;
        let fragments =
            MzMajorIntensityArray::try_new_empty(eg.fragments.clone(), ref_rt_ms.len(), 0)?;
        Ok(Self {
            eg,
            precursors,
            fragments,
            ref_rt_ms,
        })
    }

    pub fn iter_mut_precursors(
        &mut self,
    ) -> impl Iterator<Item = (&(i8, f64), MutableChromatogram<'_, V>)> {
        self.precursors.iter_mut_mzs()
    }

    pub fn iter_mut_fragments(
        &mut self,
    ) -> impl Iterator<Item = (&(T, f64), MutableChromatogram<'_, V>)> {
        self.fragments.iter_mut_mzs()
    }

    pub fn rt_range_milis(&self) -> TupleRange<u32> {
        let min = self.ref_rt_ms.first().unwrap();
        let max = self.ref_rt_ms.last().unwrap();
        (*min, *max).try_into().expect("rt range should be sorted")
    }
}
