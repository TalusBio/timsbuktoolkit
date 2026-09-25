//! Shared extraction builder -- used by both CLI (Scorer) and viewer.

use crate::data_sources::reference_library::ExpectedIntensity;
use crate::models::ExpectedIntensities;
use crate::scoring::apex_finding::Extraction;
use crate::traits::MappableRTCycles;
#[cfg(test)]
use timsquery::utils::TupleRange;
use timsquery::{
    ChromatogramCollector,
    ExtractionQuery,
    IonAnnot,
    KeyLike,
    PeakTolerance,
    QueriableData,
};

use super::skip::SkipReason;

/// Build an extraction from borrowed source geometry and resolved acquisition RT.
///
/// top_n_fragments:
///   Some(n) -> filter_zero_intensity_ions + select_top_n_fragments(n)
///   None    -> no filtering, all ions kept
pub fn build_extraction<T, I>(
    query: &ExtractionQuery<'_, impl timsquery::Target<Label = T>>,
    mut expected_intensities: ExpectedIntensities<T>,
    index: &I,
    tolerance: &PeakTolerance,
    top_n_fragments: Option<usize>,
) -> Result<Extraction<T>, SkipReason>
where
    T: KeyLike,
    I: QueriableData<ChromatogramCollector<T, f32>> + MappableRTCycles,
{
    let cycle_mapping = index.ms1_cycle_mapping();

    let mut agg = ChromatogramCollector::new(query, cycle_mapping)
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
pub fn build_extraction_into<I, Q>(
    scratch: &mut Option<Extraction<IonAnnot>>,
    query: &ExtractionQuery<'_, Q>,
    index: &I,
    tolerance: &PeakTolerance,
    top_n_fragments: Option<usize>,
) -> Result<(), SkipReason>
where
    Q: timsquery::Target<Label = IonAnnot> + ExpectedIntensity,
    I: QueriableData<ChromatogramCollector<IonAnnot, f32>> + MappableRTCycles,
{
    let cycle_mapping = index.ms1_cycle_mapping();

    match scratch {
        Some(extr) => {
            extr.chromatograms
                .try_reset_with(query, cycle_mapping)
                .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?;
        }
        None => {
            let agg = ChromatogramCollector::new(query, cycle_mapping)
                .map_err(|_| SkipReason::RetentionTimeOutOfBounds)?;
            *scratch = Some(Extraction {
                expected_intensities: ExpectedIntensities::default(),
                chromatograms: agg,
            });
        }
    }

    let extr = scratch
        .as_mut()
        .expect("extraction set by build_extraction_into");
    refill_expected(extr, query.source());
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
fn refill_expected(extraction: &mut Extraction<IonAnnot>, query: &impl ExpectedIntensity) {
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
    use timsquery::Tolerance;
    use timsquery::models::capabilities::TargetCapabilities;
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };

    #[test]
    fn extraction_preserves_source_coordinate() {
        let axis = timsquery::RtAxis::NormalizedIndex {
            scale: Some("anchors".into()),
        };
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        let fragments = [(IonAnnot::try_from("y1").unwrap(), 300.0)];
        builder.push_row(Row {
            precursor_mz: 500.0,
            charge: 2,
            frags: &fragments,
            rt: Some(timsquery::RtCoordinate {
                value: timsquery::LibraryRT(-20.0),
                axis: &axis,
            }),
            ..Default::default()
        });
        let lib = ReferenceLibrary::try_from(timsquery::serde::TargetTable::Mzpaf {
            geom: builder.seal(crate::models::DecoyPolicy::Never).unwrap(),
            frag_intens: Some(vec![1.0]),
        })
        .unwrap();
        let source = lib.iter().next().unwrap();
        let rt = timsquery::ResolvedRt::resolve(
            timsquery::RtSelection::Centered(timsquery::ObservedRTSeconds(120.0)),
            &timsquery::models::tolerance::RtTolerance::Minutes((0.5, 0.5)),
            TupleRange::try_new(
                timsquery::ObservedRTSeconds(0.0),
                timsquery::ObservedRTSeconds(300.0),
            )
            .unwrap(),
        )
        .unwrap();
        let query = ExtractionQuery::new(&source, rt);
        assert_eq!(
            query.rt().range_millis(),
            TupleRange::try_new(90_000, 150_000).unwrap()
        );
        assert_eq!(
            timsquery::Target::library_rt(query.source()),
            Some(timsquery::LibraryRT(-20.0))
        );
        assert!(std::ptr::eq(query.source(), &source));
    }

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
        fn add_query(
            &self,
            agg: &mut ChromatogramCollector<IonAnnot, f32>,
            _: &timsquery::PeakTolerance,
        ) {
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
                rt: Some(timsquery::RtCoordinate::seconds(0.0)),
                charge: 2,
                frags: &fragments[..count],
                analyte: timsquery::chemistry::analyte::Analyte::from_sequence("PEPTIDEK")
                    .as_input(),
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
                let mut materialized = timsquery::OwnedTarget::empty_like();
                materialized.reset_from(&query);
                let rt = timsquery::ResolvedRt::from_mapping(
                    timsquery::RtSelection::FullRun,
                    &tolerance.rt,
                    &index.0,
                )
                .unwrap();
                let direct_query = ExtractionQuery::new(&query, rt);
                let owned_query = ExtractionQuery::new(&materialized, rt);
                let mut shifted_direct = timsquery::SpectralCollector::new(&direct_query);
                let mut shifted_owned = timsquery::SpectralCollector::new(&owned_query);
                crate::utils::elution_group_ops::shift_fragment_isotopes(&mut shifted_direct, 1);
                crate::utils::elution_group_ops::shift_fragment_isotopes(&mut shifted_owned, 1);
                assert_eq!(
                    serde_json::to_value(&shifted_direct).unwrap(),
                    serde_json::to_value(&shifted_owned).unwrap()
                );
                let expected = ExpectedIntensities::try_from_pairs(
                    query.iter_expected_fragments(),
                    query.expected_precursor_envelope(),
                )
                .unwrap();
                let baseline = build_extraction(
                    &owned_query,
                    expected,
                    &index,
                    &tolerance.peak_tolerance(),
                    top_n,
                )
                .unwrap();
                build_extraction_into(
                    &mut slot,
                    &direct_query,
                    &index,
                    &tolerance.peak_tolerance(),
                    top_n,
                )
                .unwrap();
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
            build_extraction_into(
                &mut slot,
                &ExtractionQuery::new(
                    &query,
                    timsquery::ResolvedRt::from_mapping(
                        timsquery::RtSelection::Centered(timsquery::ObservedRTSeconds(0.01)),
                        &tolerance.rt,
                        &index.0,
                    )
                    .unwrap(),
                ),
                &index,
                &tolerance.peak_tolerance(),
                None,
            )
            .unwrap();
            assert_eq!(
                slot.as_ref().unwrap().chromatograms.rt.center().unwrap().0,
                0.01
            );
            assert_eq!(
                timsquery::Target::library_rt(&query).unwrap().0,
                0.0,
                "calibration must not change library RT"
            );
        }
    }
}
