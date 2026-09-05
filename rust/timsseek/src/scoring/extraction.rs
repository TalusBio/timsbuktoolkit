//! Shared extraction builder -- used by both CLI (Scorer) and viewer.

use crate::data_sources::reference_library::{
    ExpectedIntensity,
    RefQuery,
};
use crate::models::ExpectedIntensities;
use crate::scoring::apex_finding::Extraction;
use crate::traits::MappableRTCycles;
use timsquery::traits::QueryGeom;
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    IonAnnot,
    KeyLike,
    OptionallyRestricted,
    QueriableData,
    Tolerance,
};

use super::skip::SkipReason;

/// Build an extraction from owned query inputs.
///
/// RT range derived internally from query.rt_seconds() + tolerance,
/// clamped to the index's cycle_mapping range.
///
/// top_n_fragments:
///   Some(n) -> filter_zero_intensity_ions + select_top_n_fragments(n)
///   None    -> no filtering, all ions kept
pub fn build_extraction<T, I>(
    query: &timsquery::Target<T>,
    mut expected_intensities: ExpectedIntensities<T>,
    index: &I,
    tolerance: &Tolerance,
    top_n_fragments: Option<usize>,
) -> Result<Extraction<T>, SkipReason>
where
    T: KeyLike,
    I: QueriableData<ChromatogramCollector<T, f32>> + MappableRTCycles,
{
    let cycle_mapping = index.ms1_cycle_mapping();
    let max_range = cycle_mapping.range_milis();
    let max_range = TupleRange::try_new(max_range.0, max_range.1)
        .expect("Reference RTs should be sorted and valid");

    let rt_range = match tolerance.rt_range_as_milis(query.rt_seconds()) {
        OptionallyRestricted::Unrestricted => max_range,
        OptionallyRestricted::Restricted(r) => r,
    };

    if !max_range.intersects(rt_range) {
        return Err(SkipReason::RetentionTimeOutOfBounds);
    }

    let mut agg = ChromatogramCollector::new(query, rt_range, cycle_mapping)
        .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?;

    index.add_query(&mut agg, tolerance);

    classify_post_add_query(&agg)?;

    if let Some(n) = top_n_fragments {
        super::pipeline::filter_zero_intensity_ions(&mut agg, &mut expected_intensities);
        super::pipeline::select_top_n_fragments(&mut agg, &mut expected_intensities, n);
    }

    Ok(Extraction {
        expected_intensities,
        chromatograms: agg,
    })
}

/// Fast-path classification of a freshly-populated collector before any
/// per-cycle scoring runs. Uses the counters bumped inside `add_query`:
///
/// - `n_quad_windows_matched == 0` → fragments could not possibly match any
///   peak (scan-schedule / library mismatch), distinct from a mere absence.
/// - `n_peaks_added == 0` with nonzero quad matches → peptide is absent or
///   below LoD; no reason to run `find_apex` / `filter_zero_intensity_ions`.
fn classify_post_add_query<T: KeyLike>(
    agg: &ChromatogramCollector<T, f32>,
) -> Result<(), SkipReason> {
    if agg.n_fragment_peaks_added == 0 {
        if agg.n_quad_windows_matched == 0 {
            return Err(SkipReason::FragmentsOutsideScanRange);
        }
        return Err(SkipReason::NoObservedSignal);
    }
    Ok(())
}

/// Extract directly from a validated library query, reusing collector and
/// expected-intensity buffers. Reference values are copied only into the
/// extraction's mutable buffers, where zero-signal and top-N filtering occur.
///
/// On first call the `scratch` slot is `None` and a fresh `Extraction` is
/// allocated. Subsequent calls reset the existing one in place.
pub fn build_extraction_into<I>(
    scratch: &mut Option<Extraction<IonAnnot>>,
    query: &RefQuery<'_>,
    rt_override: Option<f32>,
    index: &I,
    tolerance: &Tolerance,
    top_n_fragments: Option<usize>,
) -> Result<(), SkipReason>
where
    I: QueriableData<ChromatogramCollector<IonAnnot, f32>> + MappableRTCycles,
{
    let cycle_mapping = index.ms1_cycle_mapping();
    let max_range = cycle_mapping.range_milis();
    let max_range = TupleRange::try_new(max_range.0, max_range.1)
        .expect("Reference RTs should be sorted and valid");

    let query_rt = rt_override.unwrap_or_else(|| query.rt_seconds());
    let rt_range = match tolerance.rt_range_as_milis(query_rt) {
        OptionallyRestricted::Unrestricted => max_range,
        OptionallyRestricted::Restricted(r) => r,
    };

    if !max_range.intersects(rt_range) {
        return Err(SkipReason::RetentionTimeOutOfBounds);
    }

    match scratch {
        Some(extr) => {
            extr.chromatograms
                .try_reset_with_overrides(query, rt_override, None, rt_range, cycle_mapping)
                .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?;
        }
        None => {
            let mut agg = ChromatogramCollector::new(query, rt_range, cycle_mapping)
                .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?;
            if let Some(rt) = rt_override {
                agg.rt_seconds = rt;
            }
            *scratch = Some(Extraction {
                expected_intensities: ExpectedIntensities::default(),
                chromatograms: agg,
            });
        }
    }

    let extr = scratch
        .as_mut()
        .expect("extraction set by build_extraction_into");
    refill_expected(extr, query);
    index.add_query(&mut extr.chromatograms, tolerance);

    classify_post_add_query(&extr.chromatograms)?;

    if let Some(n) = top_n_fragments {
        super::pipeline::filter_zero_intensity_ions(
            &mut extr.chromatograms,
            &mut extr.expected_intensities,
        );
        super::pipeline::select_top_n_fragments(
            &mut extr.chromatograms,
            &mut extr.expected_intensities,
            n,
        );
    }

    Ok(())
}

/// Copy validated reference values only into the buffers filtering will mutate.
fn refill_expected(extraction: &mut Extraction<IonAnnot>, query: &RefQuery<'_>) {
    let expected = &mut extraction.expected_intensities;
    expected.fragment_intensities.clear();
    expected
        .fragment_intensities
        .extend(query.iter_expected_fragments());
    expected.precursor_intensities.clear();
    expected
        .precursor_intensities
        .extend(query.expected_precursor_envelope());
    assert_eq!(
        expected.fragment_len(),
        extraction.chromatograms.fragments.num_ions(),
        "expected fragment intensities must match extraction geometry"
    );
    assert_eq!(
        expected.precursor_len(),
        extraction.chromatograms.precursors.num_ions(),
        "expected precursor intensities must match extraction geometry"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_sources::reference_library::ReferenceLibrary;
    use timscentroid::rt_mapping::{
        CycleToRTMapping,
        MS1CycleIndex,
    };
    use timsquery::models::capabilities::TargetCapabilities;
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };

    struct TestIndex(CycleToRTMapping<MS1CycleIndex>);

    impl MappableRTCycles for TestIndex {
        fn ms1_cycle_mapping(&self) -> &CycleToRTMapping<MS1CycleIndex> {
            &self.0
        }

        fn mobility_kind(&self) -> &timscentroid::MobilityKind {
            &timscentroid::MobilityKind::Ook0
        }
    }

    impl QueriableData<ChromatogramCollector<IonAnnot, f32>> for TestIndex {
        fn add_query(&self, agg: &mut ChromatogramCollector<IonAnnot, f32>, _: &Tolerance) {
            // First fragment and all precursors remain zero, exercising removal.
            for row in 1..agg.fragments.num_ions() {
                agg.fragments
                    .arr
                    .try_replace_row_with(row, &[1.0, 2.0])
                    .unwrap();
            }
            agg.n_quad_windows_matched = 1;
            agg.n_fragment_peaks_added = 1;
        }
    }

    #[test]
    fn direct_extraction_matches_materialized_queries_and_reuses_buffers() {
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        let fragments: Vec<_> = (1..=20)
            .map(|i| {
                (
                    IonAnnot::try_from(format!("y{i}").as_str()).unwrap(),
                    300.0 + i as f64,
                )
            })
            .collect();
        let mut intensities = Vec::new();
        for count in [20, 2] {
            builder.push_row(Row {
                precursor_mz: 500.0,
                charge: 2,
                frags: &fragments[..count],
                seq_strip: "PEPTIDEK",
                seq_mod: "PEPTIDEK",
                ..Default::default()
            });
            intensities.extend((0..count).map(|i| (i + 1) as f32));
        }
        let lib = ReferenceLibrary::try_from(timsquery::serde::TargetTable::Mzpaf {
            geom: builder.seal(crate::models::DecoyPolicy::Force).unwrap(),
            frag_intens: Some(intensities),
        })
        .unwrap();

        let index = TestIndex(CycleToRTMapping::new(vec![0, 20]));
        let tolerance = Tolerance::default();
        let mut slot = None;
        let mut allocation = None;
        for query in lib.iter() {
            // Compare both unfiltered extraction and zero-signal/top-N filtering.
            for top_n in [None, Some(8)] {
                let mut materialized = timsquery::Target::empty_like();
                materialized.reset_from(&query);
                let mut shifted_direct = timsquery::Target::empty_like();
                let mut shifted_owned = timsquery::Target::empty_like();
                crate::utils::elution_group_ops::apply_isotope_offset_fragments_into(
                    &mut shifted_direct,
                    &query,
                    1,
                );
                crate::utils::elution_group_ops::apply_isotope_offset_fragments_into(
                    &mut shifted_owned,
                    &materialized,
                    1,
                );
                assert_eq!(
                    serde_json::to_value(&shifted_direct).unwrap(),
                    serde_json::to_value(&shifted_owned).unwrap()
                );
                let expected = ExpectedIntensities::try_from_pairs(
                    query.iter_expected_fragments(),
                    query.expected_precursor_envelope(),
                )
                .unwrap();
                let baseline =
                    build_extraction(&materialized, expected, &index, &tolerance, top_n).unwrap();
                build_extraction_into(&mut slot, &query, None, &index, &tolerance, top_n).unwrap();
                let actual = slot.as_ref().unwrap();
                assert_eq!(
                    actual.expected_intensities.fragment_intensities,
                    baseline.expected_intensities.fragment_intensities
                );
                assert_eq!(
                    actual.expected_intensities.precursor_intensities,
                    baseline.expected_intensities.precursor_intensities
                );
                assert_eq!(
                    serde_json::to_value(&actual.chromatograms).unwrap(),
                    serde_json::to_value(&baseline.chromatograms).unwrap()
                );
                let ptr = actual.expected_intensities.fragment_intensities.as_ptr();
                if let Some(previous) = allocation {
                    assert_eq!(ptr, previous);
                }
                allocation = Some(ptr);
            }
            build_extraction_into(&mut slot, &query, Some(0.01), &index, &tolerance, None).unwrap();
            assert_eq!(slot.as_ref().unwrap().chromatograms.rt_seconds, 0.01);
            assert_eq!(
                query.rt_seconds(),
                0.0,
                "calibration must not change library RT"
            );
        }
    }
}
