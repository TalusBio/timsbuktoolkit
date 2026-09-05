use crate::cli::PossibleAggregator;
use crate::commands::JsonStreamSerializer;
use crate::error::CliError;
use rayon::iter::{
    IndexedParallelIterator,
    IntoParallelRefIterator,
    ParallelDrainRange,
    ParallelIterator,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::fmt::Display;
use std::io::Write;
use timscentroid::IndexedTimstofPeaks;
use timscentroid::rt_mapping::{
    CycleToRTMapping,
    MS1CycleIndex,
};
use timsquery::models::QueryRef;
use timsquery::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::ChromatogramOutput;
use timsquery::traits::{
    DecoyShift,
    QueryGeom,
};
use timsquery::{
    KeyLike,
    QueriableData,
};
use tracing::{
    error,
    warn,
};

/// Represents the output format for an aggregated spectrum.
#[derive(Debug, Serialize, Deserialize)]
pub struct SpectrumOutput {
    id: timsquery::models::OwnedSourceId,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mz: f64,
    precursor_charge: i8,
    fragment_mzs: Vec<f64>,
    fragment_intensities: Vec<f32>,
    precursor_intensities: Vec<f32>,
    precursor_labels: Vec<i8>,
}

impl SpectrumOutput {
    fn new<T: KeyLike + Display>(
        agg: &SpectralCollector<T, f32>,
        source_id: timsquery::models::SourceId<'_>,
    ) -> Self {
        let (fragment_mzs, fragment_intensities) = agg
            .iter_fragments()
            .map(|((_idx, mz), inten)| (mz, inten))
            .unzip();

        let (precursor_labels, precursor_intensities) = agg
            .iter_precursors()
            .map(|((idx, _mz), inten)| (idx, inten))
            .unzip();

        SpectrumOutput {
            id: source_id.to_owned_id(),
            mobility_ook0: agg.mobility_ook0,
            rt_seconds: agg.rt_seconds,
            precursor_mz: agg.precursor_mono_mz,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
            precursor_labels,
            precursor_charge: agg.precursor_charge as i8,
        }
    }
}

pub enum AggregatorContainer<T: KeyLike + Display> {
    Point(Vec<PointIntensityAggregator<T>>),
    Chromatogram(Vec<ChromatogramCollector<T, f32>>),
    Spectrum(Vec<SpectralCollector<T, f32>>),
}

impl<T: KeyLike + Display> AggregatorContainer<T> {
    pub fn new(
        queries: &[QueryRef<'_, T>],
        aggregator: PossibleAggregator,
        ref_rts: &CycleToRTMapping<MS1CycleIndex>,
        tolerance: &Tolerance,
    ) -> Result<Self, CliError>
    where
        T: DecoyShift,
    {
        Ok(match aggregator {
            PossibleAggregator::PointIntensity => AggregatorContainer::Point(
                queries.iter().map(PointIntensityAggregator::new).collect(),
            ),
            PossibleAggregator::Chromatogram => {
                let collectors = queries
                    .iter()
                    .map(|q| {
                        let rt_range = match tolerance.rt_range_as_milis(q.rt_seconds()) {
                            timsquery::OptionallyRestricted::Unrestricted => {
                                let range = ref_rts.range_milis();
                                timsquery::TupleRange::try_new(range.0, range.1)
                                    .expect("Reference RTs should be sorted and valid")
                            }
                            timsquery::OptionallyRestricted::Restricted(r) => r,
                        };
                        ChromatogramCollector::new(q, rt_range, ref_rts)
                            .map_err(|e| CliError::DataProcessing(format!("{:?}", e)))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                AggregatorContainer::Chromatogram(collectors)
            }
            PossibleAggregator::Spectrum => {
                AggregatorContainer::Spectrum(queries.iter().map(SpectralCollector::new).collect())
            }
        })
    }

    pub fn add_query(&mut self, index: &IndexedTimstofPeaks, tolerance: &Tolerance) {
        match self {
            AggregatorContainer::Point(aggregators) => {
                let n = aggregators.len();
                index.par_add_query_multi(aggregators, rayon::iter::repeat_n(tolerance, n));
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                let n = aggregators.len();
                index.par_add_query_multi(aggregators, rayon::iter::repeat_n(tolerance, n));
            }
            AggregatorContainer::Spectrum(aggregators) => {
                let n = aggregators.len();
                index.par_add_query_multi(aggregators, rayon::iter::repeat_n(tolerance, n));
            }
        }
    }

    /// `queries` must be the original queries in collector order. IDs are resolved
    /// here, before any empty chromatograms are filtered out.
    pub fn serialize_to_seq<W: Write>(
        &mut self,
        seq: &mut JsonStreamSerializer<W>,
        ref_rts: &CycleToRTMapping<MS1CycleIndex>,
        queries: &[QueryRef<'_, T>],
    ) -> Result<(), CliError>
    where
        T: DecoyShift,
    {
        let n = match self {
            Self::Point(a) => a.len(),
            Self::Chromatogram(a) => a.len(),
            Self::Spectrum(a) => a.len(),
        };
        if n != queries.len() {
            return Err(CliError::DataProcessing(
                "collector/query count mismatch".into(),
            ));
        }
        match self {
            AggregatorContainer::Point(aggregators) => {
                for agg in aggregators {
                    seq.serialize(agg).unwrap();
                }
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                let converted_results: Vec<ChromatogramOutput> = aggregators
                    .par_drain(..)
                    // Pair before filtering so skipped rows cannot shift output identities.
                    .zip(queries.par_iter())
                    .filter_map(|(mut agg, query)| {
                        let source_id = query.output_id();
                        match ChromatogramOutput::try_new(&mut agg, ref_rts, source_id) {
                            Ok(output) => Some(output),
                            Err(e) => {
                                match e {
                                    timsquery::errors::DataProcessingError::ExpectedNonEmptyData => {
                                        warn!(
                                            "Skipping chromatogram for elution group id {}: {:?}",
                                            source_id, e
                                        );
                                        None
                                    }
                                    _ => {
                                        error!(
                                            "Error generating chromatogram for elution group id {}: {:?}",
                                            source_id, e
                                        );
                                        panic!("Terminating due to unexpected error");
                                    }
                                }
                            }
                        }
                    })
                    .collect();

                for ser_agg in converted_results {
                    seq.serialize(&ser_agg).unwrap();
                }
            }
            AggregatorContainer::Spectrum(aggregators) => {
                for (agg, query) in aggregators.iter().zip(queries) {
                    let ser_agg = SpectrumOutput::new(agg, query.output_id());
                    seq.serialize(&ser_agg).unwrap();
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::SerializationFormat;
    use timsquery::TupleRange;
    use timsquery::serde::TargetTable;

    #[test]
    fn output_ids_survive_conversion_and_skipped_chromatograms() {
        // Exercise both JSON ID shapes through the actual reader and writer.
        for ids in [
            serde_json::json!([7, 8, 9]),
            serde_json::json!(["a", "b", "c"]),
        ] {
            let rows: Vec<_> = ids
                .as_array()
                .unwrap()
                .iter()
                .map(|id| {
                    serde_json::json!({
                        "id": id,
                        "mobility": 1.0,
                        "rt_seconds": 0.01,
                        "precursor_mz": 400.0,
                        "precursor_charge": 2,
                        "precursor_isotopes": [0],
                        "fragment_labels": ["y1"],
                        "fragments": [600.0]
                    })
                })
                .collect();
            let file = tempfile::NamedTempFile::new().unwrap();
            serde_json::to_writer(&file, &rows).unwrap();
            let TargetTable::Mzpaf { geom, .. } =
                timsquery::serde::read_targets(file.path()).unwrap()
            else {
                panic!("expected ion annotations");
            };
            let queries: Vec<_> = geom.flats().map(|f| geom.item_at(f)).collect();
            let cycles = CycleToRTMapping::new(vec![10, 20]);
            let mut spectra =
                AggregatorContainer::Spectrum(queries.iter().map(SpectralCollector::new).collect());
            let mut bytes = Vec::new();
            let mut ser = JsonStreamSerializer::new(&mut bytes, SerializationFormat::Json);
            assert!(
                spectra
                    .serialize_to_seq(&mut ser, &cycles, &queries[..2])
                    .is_err()
            );
            spectra
                .serialize_to_seq(&mut ser, &cycles, &queries)
                .unwrap();
            ser.finish().unwrap();
            let output: Vec<serde_json::Value> = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(
                output.iter().map(|r| r["id"].clone()).collect::<Vec<_>>(),
                *ids.as_array().unwrap()
            );

            let collectors = queries
                .iter()
                .enumerate()
                .map(|(i, q)| {
                    let mut agg =
                        ChromatogramCollector::new(q, TupleRange::try_new(9, 20).unwrap(), &cycles)
                            .unwrap();
                    // The middle query produces no signal and must be skipped.
                    if i != 1 {
                        agg.precursors
                            .arr
                            .try_replace_row_with(0, &[1.0, 2.0])
                            .unwrap();
                    }
                    agg
                })
                .collect();
            let mut chromatograms = AggregatorContainer::Chromatogram(collectors);
            let mut bytes = Vec::new();
            let mut ser = JsonStreamSerializer::new(&mut bytes, SerializationFormat::Json);
            chromatograms
                .serialize_to_seq(&mut ser, &cycles, &queries)
                .unwrap();
            ser.finish().unwrap();
            let output: Vec<serde_json::Value> = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(
                output.iter().map(|r| r["id"].clone()).collect::<Vec<_>>(),
                vec![ids[0].clone(), ids[2].clone()]
            );
        }
    }
}
