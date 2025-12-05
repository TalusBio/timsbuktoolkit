use serde::Serialize;

use crate::errors::DataProcessingError;
use crate::models::base::{
    ArrayElement,
    Chromatogram,
    MzMajorIntensityArray,
};
use crate::{
    KeyLike,
    TimsElutionGroup,
    ValueLike,
};
use std::sync::Arc;
use timscentroid::utils::TupleRange;

#[derive(Debug, Clone, Serialize)]
pub struct ChromatogramCollector<T: KeyLike, V: ArrayElement + ValueLike> {
    pub eg: TimsElutionGroup<T>,
    pub precursors: MzMajorIntensityArray<i8, V>,
    pub fragments: MzMajorIntensityArray<T, V>,
    pub ref_rt_ms: Arc<[u32]>, // TODO: replace with TupleRange?
}

impl<T: KeyLike, V: ValueLike + ArrayElement> ChromatogramCollector<T, V> {
    pub fn new(
        eg: TimsElutionGroup<T>,
        ref_rt_ms: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        let precursors = MzMajorIntensityArray::try_new_empty(
            eg.iter_precursors().collect(),
            ref_rt_ms.len(),
            0,
        )?;
        let fragments = MzMajorIntensityArray::try_new_empty(
            eg.iter_fragments_refs()
                .map(|(k, v)| (k.clone(), *v))
                .collect(),
            ref_rt_ms.len(),
            0,
        )?;
        Ok(Self {
            eg,
            precursors,
            fragments,
            ref_rt_ms,
        })
    }

    pub fn iter_mut_precursors(
        &mut self,
    ) -> impl Iterator<Item = (&(i8, f64), Chromatogram<'_, V>)> {
        self.precursors.iter_mut_mzs()
    }

    pub fn iter_mut_fragments(&mut self) -> impl Iterator<Item = (&(T, f64), Chromatogram<'_, V>)> {
        self.fragments.iter_mut_mzs()
    }

    pub fn rt_range_milis(&self) -> TupleRange<u32> {
        let min = self.ref_rt_ms.first().unwrap();
        let max = self.ref_rt_ms.last().unwrap();
        (*min, *max).try_into().expect("rt range should be sorted")
    }

    /// Filter ions using a predicate closure.
    ///
    /// The predicate receives a chromatogram slice and returns true to KEEP.
    /// Returns an iterator over the removed keys + m/z pairs.
    pub fn drain_nonmatching_fragments<F>(&mut self, predicate: F) -> impl Iterator<Item = (T, f64)>
    where
        F: Fn(&[V]) -> bool + Copy,
    {
        self.fragments.drain_non_matching(predicate).map(|(k, mz)| {
            self.eg.try_drop_fragment(&k).expect("key should exist");
            (k, mz)
        })
    }

    /// Filter ions using a predicate closure.
    ///
    /// The predicate receives a chromatogram slice and returns true to KEEP.
    /// Returns an iterator over the removed keys + m/z pairs.
    pub fn drain_nonmatching_precursors<F>(
        &mut self,
        predicate: F,
    ) -> impl Iterator<Item = (i8, f64)>
    where
        F: Fn(&[V]) -> bool + Copy,
    {
        self.precursors
            .drain_non_matching(predicate)
            .map(|(k, mz)| {
                self.eg.try_drop_precursor(k).expect("key should exist");
                (k, mz)
            })
    }

    /// Check if all fragments have been removed (no signal found)
    pub fn is_empty(&self) -> bool {
        self.fragments.num_ions() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tinyvec::tiny_vec;

    #[test]
    fn test_filter_ions_with_custom_predicate() {
        let eg = TimsElutionGroup::builder()
            .id(1)
            .mobility_ook0(0.8)
            .rt_seconds(100.0)
            .precursor_mzs(vec![400.0])
            .precursor_labels(tiny_vec!(0))
            .fragment_mzs(vec![600.0, 700.0, 800.0])
            .fragment_labels(tiny_vec!(0, 1, 2))
            .try_build()
            .expect("I passed valid vec lengths!");

        let rt_ms: Arc<[u32]> = vec![10, 20].into();
        let mut collector = ChromatogramCollector::<usize, f32>::new(eg, rt_ms).unwrap();

        // Set different intensities
        collector
            .precursors
            .arr
            .try_replace_row_with(0, &[5.0, 5.0])
            .unwrap(); // sum = 10.0
        collector
            .fragments
            .arr
            .try_replace_row_with(0, &[8.0, 2.0])
            .unwrap(); // sum = 10.0
        collector
            .fragments
            .arr
            .try_replace_row_with(1, &[3.0, 2.0])
            .unwrap(); // sum = 5.0
        collector
            .fragments
            .arr
            .try_replace_row_with(2, &[12.0, 3.0])
            .unwrap(); // sum = 15.0

        // Filter: keep only ions with sum >= 10.0
        let removed_keys: Vec<_> = collector
            .drain_nonmatching_fragments(|chrom| chrom.iter().sum::<f32>() >= 10.0)
            .collect();

        assert_eq!(removed_keys.len(), 1); // Only fragment 1 removed
        assert_eq!(collector.eg.precursor_count(), 1);
        assert_eq!(collector.eg.fragment_count(), 2);
        assert_eq!(collector.eg.iter_precursors().next().unwrap(), (0, 400.0));
        let mut frag_iter = collector.eg.iter_fragments();
        let _first_frag = frag_iter.next().unwrap();
        assert_eq!(frag_iter.next().unwrap(), (&2, 800.0));
    }
}
