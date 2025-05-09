use crate::{
    ElutionGroup,
    KeyLike,
};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct MobilogramSlice<const N: usize> {
    inner: [f32; N],
    start: f32,
    end: f32,
}

impl<const N: usize> MobilogramSlice<N> {
    pub fn new(start: f32, end: f32) -> Self {
        let inner = [0.0; N];
        Self { inner, start, end }
    }

    pub fn add(&mut self, loc: f32, value: f32) {
        let idx = ((loc - self.start) / (self.end - self.start) * N as f32).round() as usize;
        if idx < N {
            self.inner[idx] += value;
        }
    }

    fn many_new(start: f32, end: f32, n: usize) -> Vec<Self> {
        let mut slices = Vec::with_capacity(n);
        for _i in 0..n {
            slices.push(Self::new(start, end));
        }
        slices
    }
}

#[derive(Debug, Clone)]
pub struct EGSMAggregator<T: KeyLike, const N: usize> {
    pub eg: Arc<ElutionGroup<T>>,
    precursors: Vec<MobilogramSlice<N>>,
    fragments: Vec<MobilogramSlice<N>>,
}

impl<T: KeyLike, const N: usize> EGSMAggregator<T, N> {
    pub fn new(eg: Arc<ElutionGroup<T>>, start: f32, end: f32) -> Self {
        let precursors = MobilogramSlice::many_new(start, end, eg.precursors.len());
        let fragments = MobilogramSlice::many_new(start, end, eg.fragments.len());
        Self {
            eg,
            precursors,
            fragments,
        }
    }

    pub fn iter_mut_precursors(
        &mut self,
    ) -> impl Iterator<Item = (&(i8, f64), &mut MobilogramSlice<N>)> {
        self.eg.precursors.iter().zip(self.precursors.iter_mut())
    }

    pub fn iter_mut_fragments(
        &mut self,
    ) -> impl Iterator<Item = (&(T, f64), &mut MobilogramSlice<N>)> {
        self.eg.fragments.iter().zip(self.fragments.iter_mut())
    }
}
