use serde::Serialize;

use crate::{
    ElutionGroup,
    KeyLike,
};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize)]
pub struct EGSAggregator<T: KeyLike> {
    pub eg: Arc<ElutionGroup<T>>,
    precursors: Vec<f32>,
    fragments: Vec<f32>,
}

impl<'a, T: KeyLike> EGSAggregator<T> {
    pub fn new(eg: Arc<ElutionGroup<T>>) -> Self {
        let precursors = vec![0.0; eg.precursors.len()];
        let fragments = vec![0.0; eg.fragments.len()];
        Self {
            eg,
            precursors,
            fragments,
        }
    }

    pub fn iter_mut_precursors(&mut self) -> impl Iterator<Item = (&(i8, f64), &mut f32)> {
        self.eg.precursors.iter().zip(self.precursors.iter_mut())
    }

    pub fn iter_mut_fragments(&mut self) -> impl Iterator<Item = (&(T, f64), &mut f32)> {
        self.eg.fragments.iter().zip(self.fragments.iter_mut())
    }
}
