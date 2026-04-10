//! Shared extraction builder -- used by both CLI (Scorer) and viewer.

use crate::models::ExpectedIntensities;
use crate::scoring::apex_finding::Extraction;
use crate::traits::MappableRTCycles;
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    KeyLike,
    OptionallyRestricted,
    QueriableData,
    Tolerance,
};

use super::pipeline::SkippingReason;

/// Build an Extraction from an elution group + index query.
///
/// RT range derived internally from elution_group.rt_seconds() + tolerance,
/// clamped to the index's cycle_mapping range.
///
/// top_n_fragments:
///   Some(n) -> filter_zero_intensity_ions + select_top_n_fragments(n)
///   None    -> no filtering, all ions kept
pub(crate) fn build_extraction<T, I>(
    elution_group: &timsquery::TimsElutionGroup<T>,
    mut expected_intensities: ExpectedIntensities<T>,
    index: &I,
    tolerance: &Tolerance,
    top_n_fragments: Option<usize>,
) -> Result<Extraction<T>, SkippingReason>
where
    T: KeyLike,
    I: QueriableData<ChromatogramCollector<T, f32>> + MappableRTCycles,
{
    let cycle_mapping = index.ms1_cycle_mapping();
    let max_range = cycle_mapping.range_milis();
    let max_range = TupleRange::try_new(max_range.0, max_range.1)
        .expect("Reference RTs should be sorted and valid");

    let rt_range = match tolerance.rt_range_as_milis(elution_group.rt_seconds()) {
        OptionallyRestricted::Unrestricted => max_range,
        OptionallyRestricted::Restricted(r) => r,
    };

    if !max_range.intersects(rt_range) {
        return Err(SkippingReason::RetentionTimeOutOfBounds);
    }

    let mut agg = ChromatogramCollector::new(
        elution_group.clone(),
        rt_range,
        cycle_mapping,
    )
    .map_err(|_| SkippingReason::RetentionTimeOutOfBounds)?;

    index.add_query(&mut agg, tolerance);

    if let Some(n) = top_n_fragments {
        super::pipeline::filter_zero_intensity_ions(&mut agg, &mut expected_intensities);
        super::pipeline::select_top_n_fragments(&mut agg, &mut expected_intensities, n);
    }

    Ok(Extraction {
        expected_intensities,
        chromatograms: agg,
    })
}
