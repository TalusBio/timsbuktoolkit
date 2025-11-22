use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::QueriableData;
use timsquery::models::aggregators::ChromatogramCollector;
use timsquery::models::elution_group::ElutionGroup;
use timsquery::models::tolerance::Tolerance;

use crate::error::ViewerError;

/// Represents the output format for an aggregated chromatogram
/// (Copied from timsquery_cli for compatibility)
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ChromatogramOutput {
    pub id: u64,
    pub mobility_ook0: f32,
    pub rt_seconds: f32,
    pub precursor_mzs: Vec<f64>,
    pub fragment_mzs: Vec<f64>,
    pub precursor_intensities: Vec<Vec<f32>>,
    pub fragment_intensities: Vec<Vec<f32>>,
    pub retention_time_results_seconds: Vec<f32>,
}

impl TryFrom<ChromatogramCollector<usize, f32>> for ChromatogramOutput {
    type Error = ViewerError;

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
            return Err(ViewerError::General(format!(
                "Empty chromatogram for elution group id {}",
                value.eg.id
            )));
        }

        let (precursor_mzs, precursor_intensities): (Vec<f64>, Vec<Vec<f32>>) = value
            .iter_mut_precursors()
            .filter_map(|(&(idx, mz), cmg)| {
                let out_vec = match cmg.try_get_slice(non_zero_min_idx, non_zero_max_idx + 1) {
                    Some(slc) => slc.to_vec(),
                    None => {
                        return Some(Err(ViewerError::General(format!(
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
                        return Some(Err(ViewerError::General(format!(
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

/// Generates a chromatogram for a single elution group
pub fn generate_chromatogram(
    elution_group: &ElutionGroup<usize>,
    index: &IndexedTimstofPeaks,
    ms1_rts: Arc<[u32]>,
    tolerance: &Tolerance,
) -> Result<ChromatogramOutput, ViewerError> {
    // Create the chromatogram collector
    let mut collector = ChromatogramCollector::new(elution_group.clone(), ms1_rts)
        .map_err(|e| ViewerError::General(format!("Failed to create collector: {:?}", e)))?;

    // Query the index
    index.add_query(&mut collector, tolerance);

    // Convert to output format
    ChromatogramOutput::try_from(collector)
}
