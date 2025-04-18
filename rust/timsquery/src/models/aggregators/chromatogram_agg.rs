use serde::Serialize;

use crate::errors::DataProcessingError;
use crate::models::base::{
    MutableChromatogram,
    MzMajorIntensityArray,
};
use crate::{
    ElutionGroup,
    KeyLike,
};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize)]
pub struct EGCAggregator<T: KeyLike> {
    pub eg: Arc<ElutionGroup<T>>,
    pub precursors: MzMajorIntensityArray<i8>,
    pub fragments: MzMajorIntensityArray<T>,
}

impl<'a, T: KeyLike> EGCAggregator<T> {
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
    ) -> impl Iterator<Item = (&(i8, f64), MutableChromatogram)> {
        self.precursors.iter_mut_mzs()
    }

    pub fn iter_mut_fragments(&mut self) -> impl Iterator<Item = (&(T, f64), MutableChromatogram)> {
        self.fragments.iter_mut_mzs()
    }

    pub fn unpack(
        self,
    ) -> (
        Arc<ElutionGroup<T>>,
        MzMajorIntensityArray<i8>,
        MzMajorIntensityArray<T>,
    ) {
        (self.eg, self.precursors, self.fragments)
    }
}
