use serde::Serialize;

use crate::{
    ElutionGroup,
    KeyLike,
    ValueLike,
};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize)]
pub struct SpectralCollector<T: KeyLike, V: Default + ValueLike> {
    pub eg: Arc<ElutionGroup<T>>,
    precursors: Vec<V>,
    fragments: Vec<V>,
}

impl<T: KeyLike, V: ValueLike + Default> SpectralCollector<T, V> {
    pub fn new(eg: Arc<ElutionGroup<T>>) -> Self {
        let precursors = vec![Default::default(); eg.precursors.len()];
        let fragments = vec![Default::default(); eg.fragments.len()];
        Self {
            eg,
            precursors,
            fragments,
        }
    }

    pub fn iter_mut_precursors(&mut self) -> impl Iterator<Item = (&(i8, f64), &mut V)> {
        self.eg.precursors.iter().zip(self.precursors.iter_mut())
    }

    pub fn iter_mut_fragments(&mut self) -> impl Iterator<Item = (&(T, f64), &mut V)> {
        self.eg.fragments.iter().zip(self.fragments.iter_mut())
    }
}
