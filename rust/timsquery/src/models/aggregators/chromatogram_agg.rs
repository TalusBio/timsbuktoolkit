use serde::Serialize;

use crate::errors::DataProcessingError;
use crate::models::base::{
    ArrayElement,
    Chromatogram,
    MzMajorIntensityArray,
};
use crate::traits::queriable_data::HasQueryData;
use crate::{
    KeyLike,
    TimsElutionGroup,
    ValueLike,
};
use timscentroid::rt_mapping::{
    CycleToRTMapping,
    MS1CycleIndex,
    RTIndex, // Trait needed for the index method.
};
use timscentroid::utils::TupleRange;

// TODO: rename to `ChromatogramAccumulator` — struct carries query scalars
// (id, mobility_ook0, etc.) alongside the accumulated chromatograms, but the
// "Collector" name dates from when it owned a full TimsElutionGroup. The query
// and accumulator roles are now structurally separated; a rename would match.
#[derive(Debug, Clone, Serialize)]
pub struct ChromatogramCollector<T: KeyLike, V: ArrayElement + ValueLike> {
    // Query scalars carried from the eg at reset time.
    pub id: u64,
    pub mobility_ook0: f32,
    pub rt_seconds: f32,
    pub precursor_mono_mz: f64,
    pub precursor_charge: u8,
    /// Cached from `TimsElutionGroup::get_precursor_mz_limits()` at reset
    /// (skips negative-isotope labels — do not derive from mono_mz + charge alone).
    pub precursor_mz_limits: (f64, f64),

    // mz_order inside each array IS the (key, mz) list we need for iteration —
    // labels/mzs are NOT duplicated on the collector.
    pub precursors: MzMajorIntensityArray<i8, V>,
    pub fragments: MzMajorIntensityArray<T, V>,
    pub rt_range_ms: TupleRange<u32>,

    /// MS1 peaks written into any precursor chromatogram cell during the
    /// most recent `add_query`. Informational only — downstream fast-path
    /// decisions key off `n_fragment_peaks_added`.
    pub n_precursor_peaks_added: u64,

    /// MS2 peaks written into any fragment chromatogram cell during the
    /// most recent `add_query`. Bumped inside the inner peak-insert loop.
    /// `== 0` lets downstream code skip `filter_zero_intensity_ions` +
    /// `find_apex` work for peptides that did not match any MS2 signal.
    pub n_fragment_peaks_added: u64,

    /// Quad isolation windows that intersected the query's precursor range
    /// during the most recent `add_query`. Distinguishes "fragments fall
    /// outside the run's quad schedule entirely" (counter == 0) from
    /// "fragments in window but no observed peaks" (counter > 0 &&
    /// n_fragment_peaks_added == 0). Note: exact value is impl-defined
    /// (eager counts windows once, lazy counts windows×fragments); only
    /// the zero / nonzero distinction is meaningful.
    pub n_quad_windows_matched: u32,
}

impl<T: KeyLike, V: ValueLike + ArrayElement> ChromatogramCollector<T, V> {
    pub fn new(
        eg: &TimsElutionGroup<T>,
        rt_range_ms: TupleRange<u32>,
        ref_rt_ms: &CycleToRTMapping<MS1CycleIndex>,
    ) -> Result<Self, DataProcessingError> {
        let start = ref_rt_ms.ms_to_closest_index(rt_range_ms.start());
        let end = ref_rt_ms.ms_to_closest_index(rt_range_ms.end());
        let num_cycles = end.index() - start.index() + 1;

        let precursor_order: Vec<_> = eg.iter_precursors().collect();
        let fragment_order: Vec<_> = eg
            .iter_fragments_refs()
            .map(|(k, v)| (k.clone(), *v))
            .collect();
        if precursor_order.is_empty() && fragment_order.is_empty() {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        if num_cycles == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }

        let precursors =
            MzMajorIntensityArray::try_new_empty(precursor_order, num_cycles, start.index())
                .expect("Already checked non-empty");
        let fragments =
            MzMajorIntensityArray::try_new_empty(fragment_order, num_cycles, start.index())?;
        Ok(Self {
            id: eg.id(),
            mobility_ook0: eg.mobility_ook0(),
            rt_seconds: eg.rt_seconds(),
            precursor_mono_mz: eg.mono_precursor_mz(),
            precursor_charge: eg.precursor_charge(),
            precursor_mz_limits: eg.get_precursor_mz_limits(),
            precursors,
            fragments,
            rt_range_ms,
            n_precursor_peaks_added: 0,
            n_fragment_peaks_added: 0,
            n_quad_windows_matched: 0,
        })
    }

    pub fn try_reset_with(
        &mut self,
        eg: &TimsElutionGroup<T>,
        rt_range_ms: TupleRange<u32>,
        ref_rt_ms: &CycleToRTMapping<MS1CycleIndex>,
    ) -> Result<(), DataProcessingError> {
        self.try_reset_with_overrides(eg, None, None, rt_range_ms, ref_rt_ms)
    }

    /// Like `try_reset_with` but lets callers override `rt_seconds` / `mobility_ook0`
    /// without rebuilding the source eg — replaces the `eg.clone().with_rt_seconds(..)`
    /// and `eg.clone().with_mobility(..)` clone-then-mutate pattern.
    pub fn try_reset_with_overrides(
        &mut self,
        eg: &TimsElutionGroup<T>,
        rt_override: Option<f32>,
        mobility_override: Option<f32>,
        rt_range_ms: TupleRange<u32>,
        ref_rt_ms: &CycleToRTMapping<MS1CycleIndex>,
    ) -> Result<(), DataProcessingError> {
        let start = ref_rt_ms.ms_to_closest_index(rt_range_ms.start());
        let end = ref_rt_ms.ms_to_closest_index(rt_range_ms.end());
        let num_cycles = end.index() - start.index() + 1;

        let precursor_order: Vec<_> = eg.iter_precursors().collect();
        let fragment_order: Vec<_> = eg
            .iter_fragments_refs()
            .map(|(k, v)| (k.clone(), *v))
            .collect();
        if precursor_order.is_empty() && fragment_order.is_empty() {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        if num_cycles == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }

        self.id = eg.id();
        self.mobility_ook0 = mobility_override.unwrap_or_else(|| eg.mobility_ook0());
        self.rt_seconds = rt_override.unwrap_or_else(|| eg.rt_seconds());
        self.precursor_mono_mz = eg.mono_precursor_mz();
        self.precursor_charge = eg.precursor_charge();
        self.precursor_mz_limits = eg.get_precursor_mz_limits();
        self.rt_range_ms = rt_range_ms;
        self.n_precursor_peaks_added = 0;
        self.n_fragment_peaks_added = 0;
        self.n_quad_windows_matched = 0;
        self.precursors
            .clear_with_order(precursor_order, num_cycles, start.index());
        self.fragments
            .clear_with_order(fragment_order, num_cycles, start.index());
        Ok(())
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
        self.rt_range_ms
    }

    /// Filter ions using a predicate closure.
    ///
    /// The predicate receives a chromatogram slice and returns true to KEEP.
    /// Returns an iterator over the removed keys + m/z pairs.
    pub fn drain_nonmatching_fragments<F>(&mut self, predicate: F) -> impl Iterator<Item = (T, f64)>
    where
        F: Fn(&[V]) -> bool + Copy,
    {
        self.fragments.drain_non_matching(predicate)
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
        self.precursors.drain_non_matching(predicate)
    }

    /// Check if all fragments have been removed (no signal found)
    pub fn is_empty(&self) -> bool {
        self.fragments.num_ions() == 0
    }

    pub fn num_cycles(&self) -> usize {
        assert_eq!(
            self.precursors.num_cycles(),
            self.fragments.num_cycles(),
            "Precursors and fragments must have the same number of cycles"
        );
        self.precursors.num_cycles()
    }

    pub fn cycle_offset(&self) -> usize {
        assert_eq!(
            self.precursors.cycle_offset, self.fragments.cycle_offset,
            "Precursors and fragments must have the same cycle offset"
        );
        self.precursors.cycle_offset
    }
}

impl<T: KeyLike, V: ArrayElement + ValueLike> HasQueryData<T> for ChromatogramCollector<T, V> {
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
        self.precursors.iter_mzs().map(|((k, mz), _)| (*k, *mz))
    }

    fn iter_fragments<'a>(&'a self) -> impl Iterator<Item = (&'a T, f64)> + 'a
    where
        T: 'a,
    {
        self.fragments.iter_mzs().map(|((k, mz), _)| (k, *mz))
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
            .precursor(400.0, 1u8)
            .precursor_labels(tiny_vec!(0))
            .fragment_mzs(vec![600.0, 700.0, 800.0])
            .fragment_labels(tiny_vec!(0, 1, 2))
            .try_build()
            .expect("I passed valid vec lengths!");

        let rt_ms = CycleToRTMapping::new(vec![10, 20]);
        let mut collector = ChromatogramCollector::<usize, f32>::new(
            &eg,
            TupleRange::try_new(9, 20).unwrap(),
            &rt_ms,
        )
        .unwrap();

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
        assert_eq!(collector.precursors.num_ions(), 1);
        assert_eq!(collector.fragments.num_ions(), 2);
        let mut prec_iter = collector.precursors.iter_mzs();
        assert_eq!(prec_iter.next().unwrap().0, &(0, 400.0));
        let mut frag_iter = collector.fragments.iter_mzs();
        let _first_frag = frag_iter.next().unwrap();
        assert_eq!(frag_iter.next().unwrap().0, &(2usize, 800.0));
    }
}
