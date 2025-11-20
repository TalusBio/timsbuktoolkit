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
use timsquery::models::elution_group::ElutionGroup;
use timsquery::models::tolerance::Tolerance;
use tracing::warn;

/// Represents a user-friendly format for specifying elution groups in an input file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElutionGroupInput {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursors: Vec<f64>,
    pub fragments: Vec<f64>,
}

impl From<ElutionGroupInput> for ElutionGroup<usize> {
    fn from(val: ElutionGroupInput) -> Self {
        let precursors: Arc<[(i8, f64)]> = Arc::from(
            val.precursors
                .into_iter()
                .enumerate()
                .map(|(i, mz)| (i as i8, mz))
                .collect::<Vec<(i8, f64)>>(),
        );
        let fragments: Arc<[(usize, f64)]> = Arc::from(
            val.fragments
                .into_iter()
                .enumerate()
                .collect::<Vec<(usize, f64)>>(),
        );
        ElutionGroup {
            id: val.id,
            mobility: val.mobility,
            rt_seconds: val.rt_seconds,
            precursors,
            fragments,
        }
    }
}

/// Represents the output format for an aggregated spectrum.
#[derive(Debug, Serialize, Deserialize)]
pub struct SpectrumOutput {
    id: u64,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mzs: Vec<f64>,
    fragment_mzs: Vec<f64>,
    precursor_intensities: Vec<f32>,
    fragment_intensities: Vec<f32>,
}

/// Represents the output format for an aggregated chromatogram.
#[derive(Debug, Serialize, Deserialize)]
pub struct ChromatogramOutput {
    id: u64,
    mobility_ook0: f32,
    rt_seconds: f32,
    precursor_mzs: Vec<f64>,
    fragment_mzs: Vec<f64>,
    precursor_intensities: Vec<Vec<f32>>,
    fragment_intensities: Vec<Vec<f32>>,
    retention_time_results_seconds: Vec<f32>,
}

impl TryFrom<ChromatogramCollector<usize, f32>> for ChromatogramOutput {
    type Error = CliError;

    fn try_from(
        mut value: ChromatogramCollector<usize, f32>,
    ) -> Result<ChromatogramOutput, Self::Error> {
        let mut non_zero_min_idx = value.ref_rt_ms.len();
        let mut non_zero_max_idx = 0usize;

        value.iter_mut_precursors().for_each(|((_idx, _mz), cmg)| {
            let slc = cmg.as_slice();
            for (i, &inten) in slc.iter().enumerate() {
                if inten > 0.0 {
                    let abs_idx = i;
                    if abs_idx < non_zero_min_idx {
                        non_zero_min_idx = abs_idx;
                    }
                    if abs_idx > non_zero_max_idx {
                        non_zero_max_idx = abs_idx;
                    }
                }
            }
        });

        value.iter_mut_fragments().for_each(|((_idx, _mz), cmg)| {
            let slc = cmg.as_slice();
            for (i, &inten) in slc.iter().enumerate() {
                if inten > 0.0 {
                    let abs_idx = i;
                    if abs_idx < non_zero_min_idx {
                        non_zero_min_idx = abs_idx;
                    }
                    if abs_idx > non_zero_max_idx {
                        non_zero_max_idx = abs_idx;
                    }
                }
            }
        });

        if non_zero_min_idx > non_zero_max_idx {
            return Err(CliError::EmptyChromatogram(value.eg.id));
        }

        let (precursor_mzs, precursor_intensities): (Vec<f64>, Vec<Vec<f32>>) = value
            .iter_mut_precursors()
            .filter_map(|(&(idx, mz), cmg)| {
                let out_vec = match cmg.try_get_slice(non_zero_min_idx, non_zero_max_idx + 1) {
                    Some(slc) => slc.to_vec(),
                    None => {
                        return Some(Err(CliError::NonRecoverableError(format!(
                            "Failed to get slice for precursor mz {} in chromatogram id {}",
                            mz, idx,
                        ))));
                    }
                };
                if out_vec.iter().all(|&x| x == 0.0) {
                    return None;
                }
                Some(Ok((mz, out_vec)))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .unzip();

        let (fragment_mzs, fragment_intensities): (Vec<f64>, Vec<Vec<f32>>) = value
            .iter_mut_fragments()
            .filter_map(|(&(idx, mz), cmg)| {
                let out_vec = match cmg.try_get_slice(non_zero_min_idx, non_zero_max_idx + 1) {
                    Some(slc) => slc.to_vec(),
                    None => {
                        return Some(Err(CliError::NonRecoverableError(format!(
                            "Failed to get slice for fragment mz {} in chromatogram id {}",
                            mz, idx,
                        ))));
                    }
                };
                if out_vec.iter().all(|&x| x == 0.0) {
                    return None;
                }
                Some(Ok((mz, out_vec)))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .unzip();

        Ok(ChromatogramOutput {
            id: value.eg.id,
            mobility_ook0: value.eg.mobility,
            rt_seconds: value.eg.rt_seconds,
            precursor_mzs,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
            retention_time_results_seconds: value.ref_rt_ms[non_zero_min_idx..=non_zero_max_idx]
                .iter()
                .map(|&x| x as f32 / 1000.0)
                .collect(),
        })
    }
}

impl From<&SpectralCollector<usize, f32>> for SpectrumOutput {
    fn from(agg: &SpectralCollector<usize, f32>) -> Self {
        let (precursor_mzs, precursor_intensities): (Vec<f64>, Vec<f32>) = agg
            .iter_precursors()
            .map(|((_idx, mz), inten)| (mz, inten))
            .unzip();
        let (fragment_mzs, fragment_intensities) = agg
            .iter_fragments()
            .map(|((_idx, mz), inten)| (mz, inten))
            .unzip();

        SpectrumOutput {
            id: agg.eg.id,
            mobility_ook0: agg.eg.mobility,
            rt_seconds: agg.eg.rt_seconds,
            precursor_mzs,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
        }
    }
}

pub enum AggregatorContainer {
    Point(Vec<PointIntensityAggregator<usize>>),
    Chromatogram(Vec<ChromatogramCollector<usize, f32>>),
    Spectrum(Vec<SpectralCollector<usize, f32>>),
}

impl AggregatorContainer {
    pub fn new(
        queries: Vec<ElutionGroup<usize>>,
        aggregator: PossibleAggregator,
        ref_rts: Arc<[u32]>,
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
                        ChromatogramCollector::new(x, ref_rts.clone())
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

    pub fn serialize_to_seq<S>(&mut self, seq: &mut S) -> Result<(), CliError>
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
                        let agg_id = agg.eg.id;
                        match ChromatogramOutput::try_from(agg) {
                            Ok(output) => Some(output),
                            Err(e) => {
                                if !matches!(e, CliError::EmptyChromatogram(_)) {
                                    warn!(
                                        "Skipping chromatogram for elution group id {}: {}",
                                        agg_id, e
                                    );
                                }
                                None
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
