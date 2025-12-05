use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
use timsquery::models::aggregators::ChromatogramCollector;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::{
    KeyLike,
    QueriableData,
};

use crate::error::ViewerError;
use crate::plot_renderer::MS2Spectrum;
use tracing::instrument;

/// Smoothing method configuration
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize, PartialEq)]
pub enum SmoothingMethod {
    None,
    SavitzkyGolay { window: usize, polynomial: usize },
    Gaussian { sigma: f32 },
}

impl SmoothingMethod {
    /// Default sigma value for Gaussian smoothing.
    pub const DEFAULT_GAUSSIAN_SIGMA: f32 = 2.0;
    /// Default polynomial order for Savitzky-Golay smoothing.
    pub const DEFAULT_SG_POLYNOMIAL: usize = 2;
    /// Default window size for Savitzky-Golay smoothing.
    pub const DEFAULT_SG_WINDOW: usize = 5;

    pub fn default_savitzky_golay() -> Self {
        Self::SavitzkyGolay {
            window: Self::DEFAULT_SG_WINDOW,
            polynomial: Self::DEFAULT_SG_POLYNOMIAL,
        }
    }

    pub fn default_gaussian() -> Self {
        Self::Gaussian {
            sigma: Self::DEFAULT_GAUSSIAN_SIGMA,
        }
    }
}

impl Default for SmoothingMethod {
    fn default() -> Self {
        Self::None
    }
}

pub fn apply_smoothing(intensities: &[f32], method: &SmoothingMethod) -> Option<Vec<f32>> {
    match method {
        SmoothingMethod::None => None,
        SmoothingMethod::SavitzkyGolay { window, polynomial } => {
            Some(savitzky_golay_smooth(intensities, *window, *polynomial))
        }
        SmoothingMethod::Gaussian { sigma } => Some(gaussian_smooth(intensities, *sigma)),
    }
}

/// Savitzky-Golay smoothing filter
/// Uses a simplified weighted moving average approximation rather than full least squares
fn savitzky_golay_smooth(data: &[f32], window: usize, polynomial: usize) -> Vec<f32> {
    if data.len() < window || window < 3 || window.is_multiple_of(2) {
        return data.to_vec();
    }

    if polynomial >= window {
        return data.to_vec();
    }

    let half_window = window / 2;
    let mut smoothed = Vec::with_capacity(data.len());

    let weights = compute_savitzky_golay_weights(window, polynomial);

    for point_idx in 0..data.len() {
        let mut weighted_sum = 0.0;
        let mut weight_total = 0.0;

        for window_offset in 0..window {
            let data_idx = (point_idx + window_offset).saturating_sub(half_window);
            if data_idx < data.len() {
                weighted_sum += data[data_idx] * weights[window_offset];
                weight_total += weights[window_offset];
            }
        }

        smoothed.push(weighted_sum / weight_total);
    }

    smoothed
}

/// Compute Savitzky-Golay filter weights.
///
/// Simplified version using polynomial-based weight decay from center,
/// rather than proper least-squares fitting. Suitable for visualization
/// smoothing but not for precise quantitative analysis.
fn compute_savitzky_golay_weights(window_size: usize, polynomial_order: usize) -> Vec<f32> {
    let half_window = window_size / 2;
    let mut weights = vec![0.0; window_size];

    for position in 0..window_size {
        let distance_from_center = ((position as isize) - (half_window as isize)).abs() as f32;
        let normalized_distance = distance_from_center / (half_window as f32);

        weights[position] = match polynomial_order {
            0 | 1 => 1.0 - normalized_distance,
            2 => (1.0 - normalized_distance * normalized_distance).max(0.0),
            3 => (1.0 - normalized_distance.powi(3)).max(0.0),
            _ => (1.0 - normalized_distance.powi(polynomial_order as i32)).max(0.0),
        };
    }

    weights
}

/// Gaussian smoothing filter
fn gaussian_smooth(data: &[f32], sigma: f32) -> Vec<f32> {
    if sigma <= 0.0 || data.len() < 3 {
        return data.to_vec();
    }

    let kernel_radius = (3.0 * sigma).ceil() as usize;
    let kernel_size = 2 * kernel_radius + 1;
    let mut kernel = vec![0.0; kernel_size];

    // Compute Gaussian kernel weights
    let mut normalization_sum = 0.0;
    for kernel_pos in 0..kernel_size {
        let distance_from_center = (kernel_pos as f32) - (kernel_radius as f32);
        let gaussian_value = (-0.5 * (distance_from_center / sigma).powi(2)).exp();
        kernel[kernel_pos] = gaussian_value;
        normalization_sum += gaussian_value;
    }

    // Normalize kernel so weights sum to 1.0
    for weight in kernel.iter_mut() {
        *weight /= normalization_sum;
    }

    let mut smoothed = Vec::with_capacity(data.len());

    // Apply Gaussian kernel convolution
    for point_idx in 0..data.len() {
        let mut weighted_sum = 0.0;
        let mut weight_total = 0.0;

        for kernel_offset in 0..kernel_size {
            let data_idx = (point_idx + kernel_offset).saturating_sub(kernel_radius);
            if data_idx < data.len() {
                weighted_sum += data[data_idx] * kernel[kernel_offset];
                weight_total += kernel[kernel_offset];
            }
        }

        smoothed.push(weighted_sum / weight_total);
    }

    smoothed
}

/// Represents the output format for an aggregated chromatogram
/// It is pretty ugly performance-wise but I am keeping it as is because
/// we are meant to have only one for the whole runtime of the viewer.
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
}

impl<T: KeyLike + std::fmt::Display> TryFrom<ChromatogramCollector<T, f32>> for ChromatogramOutput {
    type Error = ViewerError;

    fn try_from(
        mut value: ChromatogramCollector<T, f32>,
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
                value.eg.id()
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

        let ((fragment_mzs, fragment_intensities), fragment_labels): (
            (Vec<f64>, Vec<Vec<f32>>),
            Vec<String>,
        ) = value
            .iter_mut_fragments()
            .filter_map(|(&(ref idx, mz), cmg)| {
                let out_vec = match cmg.try_get_slice(non_zero_min_idx, non_zero_max_idx + 1) {
                    Some(slc) => slc.to_vec(),
                    None => {
                        return Some(Err(ViewerError::General(format!(
                            "Failed to get slice for fragment mz {} in chromatogram id {:#?}",
                            mz,
                            idx.clone(),
                        ))));
                    }
                };
                if out_vec.iter().all(|&x| x == 0.0) {
                    return None;
                }
                Some(Ok(((mz, out_vec), format!("{}", idx))))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .unzip();

        Ok(ChromatogramOutput {
            id: value.eg.id(),
            mobility_ook0: value.eg.mobility_ook0(),
            rt_seconds: value.eg.rt_seconds(),
            precursor_mzs,
            fragment_mzs,
            precursor_intensities,
            fragment_intensities,
            fragment_labels,
            retention_time_results_seconds: value.ref_rt_ms[non_zero_min_idx..=non_zero_max_idx]
                .iter()
                .map(|&x| x as f32 / 1000.0)
                .collect(),
        })
    }
}

impl ChromatogramOutput {
    /// Apply smoothing to all intensity traces in the chromatogram
    pub fn apply_smoothing(&mut self, method: &SmoothingMethod) {
        if matches!(method, SmoothingMethod::None) {
            return;
        }

        // Smooth precursor intensities
        for intensities in self.precursor_intensities.iter_mut() {
            if let Some(smoothed) = apply_smoothing(intensities, method) {
                *intensities = smoothed;
            }
        }

        // Smooth fragment intensities
        for intensities in self.fragment_intensities.iter_mut() {
            if let Some(smoothed) = apply_smoothing(intensities, method) {
                *intensities = smoothed;
            }
        }
    }
}

#[instrument(skip_all)]
pub fn generate_chromatogram<T: KeyLike + std::fmt::Display>(
    elution_group: &TimsElutionGroup<T>,
    index: &IndexedTimstofPeaks,
    ms1_rts: Arc<[u32]>,
    tolerance: &Tolerance,
    smoothing: &SmoothingMethod,
) -> Result<ChromatogramOutput, ViewerError> {
    let mut collector = ChromatogramCollector::new(elution_group.clone(), ms1_rts)
        .map_err(|e| ViewerError::General(format!("Failed to create collector: {:?}", e)))?;

    index.add_query(&mut collector, tolerance);

    let mut output = ChromatogramOutput::try_from(collector)?;

    output.apply_smoothing(smoothing);

    Ok(output)
}

#[instrument(skip(chromatogram))]
pub fn extract_ms2_spectrum_from_chromatogram(
    chromatogram: &ChromatogramOutput,
    rt_seconds: f64,
) -> Result<MS2Spectrum, ViewerError> {
    let rt_index = chromatogram
        .retention_time_results_seconds
        .partition_point(|&rt| (rt as f64) < rt_seconds)
        .min(
            chromatogram
                .retention_time_results_seconds
                .len()
                .saturating_sub(1),
        );

    let mut mz_values = Vec::new();
    let mut intensities = Vec::new();
    let mut labels = Vec::new();

    for (frag_idx, (&mz, intensity_vec)) in chromatogram
        .fragment_mzs
        .iter()
        .zip(chromatogram.fragment_intensities.iter())
        .enumerate()
    {
        if rt_index < intensity_vec.len() {
            mz_values.push(mz);
            intensities.push(intensity_vec[rt_index]);
            labels.push(
                chromatogram
                    .fragment_labels
                    .get(frag_idx)
                    .cloned()
                    .unwrap_or_else(|| format!("Fragment {}", frag_idx + 1)),
            );
        }
    }

    let actual_rt = chromatogram
        .retention_time_results_seconds
        .get(rt_index)
        .copied()
        .unwrap_or(rt_seconds as f32) as f64;

    Ok(MS2Spectrum {
        mz_values,
        intensities,
        rt_seconds: actual_rt,
        fragment_labels: labels,
    })
}
