use crate::{
    ChromatogramCollector,
    KeyLike,
};

/// Represents the output format for an aggregated chromatogram
/// It is pretty ugly performance-wise but I am keeping it as is because
/// we are meant to have only one for the whole runtime of the viewer.
/// and not that many in the cli (RN max == batch size) so it should be fine.
/// I can optimize it as needed but probably the biggest gain will be
/// in the logic so we keep them as buffers, so only 1 exists at a time.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ChromatogramOutput {
    pub id: u64,
    pub mobility_ook0: f32,
    pub rt_seconds: f32,
    pub precursor_mzs: Vec<f64>,
    pub fragment_mzs: Vec<f64>,
    pub precursor_intensities: Vec<Vec<f32>>,
    pub fragment_intensities: Vec<Vec<f32>>,
    pub fragment_labels: Vec<String>,
    pub retention_time_results_seconds: Vec<f32>,
    /// Library (predicted) fragment intensities aligned with `fragment_labels`.
    /// Used for mirror plot visualization to compare observed vs predicted spectra.
    /// Populated from library sidecar when available.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub library_fragment_intensities: Option<Vec<f32>>,
}

impl ChromatogramOutput {
    pub fn try_new<T: KeyLike + std::fmt::Display>(
        mut collector: ChromatogramCollector<T, f32>,
        reference_cycles: &[u32],
    ) -> Result<Self, crate::errors::DataProcessingError> {
        let mut local_non_zero_min_idx = collector.num_cycles();
        let mut local_non_zero_max_idx = 0usize;

        collector
            .iter_mut_precursors()
            .for_each(|((_idx, _mz), cmg)| {
                let slc = cmg.as_slice();
                for (i, &inten) in slc.iter().enumerate() {
                    if inten > 0.0 {
                        let abs_idx = i;
                        if abs_idx < local_non_zero_min_idx {
                            local_non_zero_min_idx = abs_idx;
                        }
                        if abs_idx > local_non_zero_max_idx {
                            local_non_zero_max_idx = abs_idx;
                        }
                    }
                }
            });

        collector
            .iter_mut_fragments()
            .for_each(|((_idx, _mz), cmg)| {
                let slc = cmg.as_slice();
                for (i, &inten) in slc.iter().enumerate() {
                    if inten > 0.0 {
                        let abs_idx = i;
                        if abs_idx < local_non_zero_min_idx {
                            local_non_zero_min_idx = abs_idx;
                        }
                        if abs_idx > local_non_zero_max_idx {
                            local_non_zero_max_idx = abs_idx;
                        }
                    }
                }
            });

        if local_non_zero_min_idx > local_non_zero_max_idx {
            return Err(crate::errors::DataProcessingError::ExpectedNonEmptyData);
        }

        let cycle_offset = collector.cycle_offset();
        let (precursor_mzs, precursor_intensities): (Vec<f64>, Vec<Vec<f32>>) = collector
            .iter_mut_precursors()
            .filter_map(|(&(idx, mz), cmg)| {
                let out_vec = cmg
                    .try_get_slice(
                        local_non_zero_min_idx + cycle_offset,
                        local_non_zero_max_idx + 1 + cycle_offset,
                    )
                    .unwrap();

                if out_vec.iter().all(|&x| x == 0.0) {
                    return None;
                }
                Some((mz, out_vec.to_vec()))
            })
            .unzip();

        let non_zero_min_idx = local_non_zero_min_idx + collector.cycle_offset();
        let non_zero_max_idx = local_non_zero_max_idx + collector.cycle_offset();

        let ((fragment_mzs, fragment_intensities), fragment_labels): (
            (Vec<f64>, Vec<Vec<f32>>),
            Vec<String>,
        ) = collector
            .iter_mut_fragments()
            .filter_map(|(&(ref idx, mz), cmg)| {
                let out_slc = cmg
                    .try_get_slice(non_zero_min_idx, non_zero_max_idx + 1)
                    .expect("Failed to get slice from chromatogram");
                if out_slc.iter().all(|&x| x == 0.0) {
                    return None;
                }
                let out_vec = out_slc.to_vec();
                Some(Ok(((mz, out_vec), format!("{}", idx))))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .unzip();

        Ok(ChromatogramOutput {
            id: collector.eg.id(),
            mobility_ook0: collector.eg.mobility_ook0(),
            rt_seconds: collector.eg.rt_seconds(),
            precursor_mzs,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
            fragment_labels,
            retention_time_results_seconds: reference_cycles[non_zero_min_idx..=non_zero_max_idx]
                .iter()
                .map(|&x| x as f32 / 1000.0)
                .collect(),
            library_fragment_intensities: None,
        })
    }
}
