use crate::cli::PossibleAggregator;
use crate::error::CliError;
use rayon::iter::{
    ParallelDrainRange,
    ParallelIterator,
};
use serde::ser::SerializeSeq;
use serde::{
    Deserialize,
    Serialize,
};
use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::QueriableData;
use timsquery::models::aggregators::{
    ChromatogramCollector,
    PointIntensityAggregator,
    SpectralCollector,
};
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::serde::ChromatogramOutput;
use tracing::{
    error,
    warn,
};

/// Represents the output format for an aggregated spectrum.
#[derive(Debug, Serialize, Deserialize)]
pub struct SpectrumOutput {
    id: u64,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mz: f64,
    precursor_charge: i8,
    fragment_mzs: Vec<f64>,
    fragment_intensities: Vec<f32>,
    precursor_intensities: Vec<f32>,
    precursor_labels: Vec<i8>,
}

impl From<&SpectralCollector<u8, f32>> for SpectrumOutput {
    fn from(agg: &SpectralCollector<u8, f32>) -> Self {
        let (fragment_mzs, fragment_intensities) = agg
            .iter_fragments()
            .map(|((_idx, mz), inten)| (mz, inten))
            .unzip();

        let (precursor_labels, precursor_intensities) = agg
            .iter_precursors()
            .map(|((idx, _mz), inten)| (idx, inten))
            .unzip();

        SpectrumOutput {
            id: agg.eg.id(),
            mobility_ook0: agg.eg.mobility_ook0(),
            rt_seconds: agg.eg.rt_seconds(),
            precursor_mz: agg.eg.precursor_mz(),
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
            precursor_labels,
            precursor_charge: agg.eg.precursor_charge() as i8,
        }
    }
}

pub enum AggregatorContainer {
    Point(Vec<PointIntensityAggregator<u8>>),
    Chromatogram(Vec<ChromatogramCollector<u8, f32>>),
    Spectrum(Vec<SpectralCollector<u8, f32>>),
}

impl AggregatorContainer {
    pub fn new(
        queries: Vec<TimsElutionGroup<u8>>,
        aggregator: PossibleAggregator,
        ref_rts: Arc<[u32]>,
        tolerance: &Tolerance,
    ) -> Result<Self, CliError> {
        Ok(match aggregator {
            PossibleAggregator::PointIntensityAggregator => AggregatorContainer::Point(
                queries
                    .into_iter()
                    .map(|x| PointIntensityAggregator::new_with_elution_group(x.into()))
                    .collect(),
            ),
            PossibleAggregator::ChromatogramAggregator => {
                let collectors = queries
                    .into_iter()
                    .map(|x| {
                        let rt_range = match tolerance.rt_range_as_milis(x.rt_seconds()) {
                            timsquery::OptionallyRestricted::Unrestricted => {
                                timsquery::TupleRange::try_new(
                                    *ref_rts.first().unwrap(),
                                    *ref_rts.last().unwrap(),
                                )
                                .expect("Reference RTs should be sorted and valid")
                            }
                            timsquery::OptionallyRestricted::Restricted(r) => r,
                        };
                        ChromatogramCollector::new(x, rt_range, &ref_rts)
                            .map_err(|e| CliError::DataProcessing(format!("{:?}", e)))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                AggregatorContainer::Chromatogram(collectors)
            }
            PossibleAggregator::SpectrumAggregator => AggregatorContainer::Spectrum(
                queries.into_iter().map(SpectralCollector::new).collect(),
            ),
        })
    }

    pub fn add_query(&mut self, index: &IndexedTimstofPeaks, tolerance: &Tolerance) {
        match self {
            AggregatorContainer::Point(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            }
            AggregatorContainer::Spectrum(aggregators) => {
                index.par_add_query_multi(aggregators, tolerance);
            }
        }
    }

    pub fn serialize_to_seq<S>(&mut self, seq: &mut S, ref_rts: &[u32]) -> Result<(), CliError>
    where
        S: SerializeSeq,
    {
        match self {
            AggregatorContainer::Point(aggregators) => {
                for agg in aggregators {
                    seq.serialize_element(agg).unwrap();
                }
            }
            AggregatorContainer::Chromatogram(aggregators) => {
                let converted_results: Vec<ChromatogramOutput> = aggregators
                    .par_drain(..)
                    .filter_map(|agg| {
                        let agg_id = agg.eg.id();
                        match ChromatogramOutput::try_new(agg, ref_rts) {
                            Ok(output) => Some(output),
                            Err(e) => {
                                // if !matches!(e, CliError::EmptyChromatogram(_)) {
                                //     warn!(
                                //         "Skipping chromatogram for elution group id {}: {}",
                                //         agg_id, e
                                //     );
                                // }
                                // None
                                match e {
                                    timsquery::errors::DataProcessingError::ExpectedNonEmptyData => {
                                        warn!(
                                            "Skipping chromatogram for elution group id {}: {:?}",
                                            agg_id, e
                                        );
                                        None
                                    }
                                    _ => {
                                        error!(
                                            "Error generating chromatogram for elution group id {}: {:?}",
                                            agg_id, e
                                        );
                                        panic!("Terminating due to unexpected error");
                                    }
                                }
                            }
                        }
                    })
                    .collect();

                for ser_agg in converted_results {
                    seq.serialize_element(&ser_agg).unwrap();
                }
            }
            AggregatorContainer::Spectrum(aggregators) => {
                for agg in aggregators.iter() {
                    let ser_agg = SpectrumOutput::from(agg);
                    seq.serialize_element(&ser_agg).unwrap();
                }
            }
        }
        Ok(())
    }
}
