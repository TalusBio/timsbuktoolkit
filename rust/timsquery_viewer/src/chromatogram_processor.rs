use std::sync::Arc;
use timscentroid::IndexedTimstofPeaks;
pub use timsquery::models::aggregators::ChromatogramCollector;
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::models::tolerance::Tolerance;
use timsquery::{
    KeyLike,
    QueriableData,
};

use crate::error::ViewerError;
use crate::plot_renderer::MS2Spectrum;
pub use timsquery::serde::ChromatogramOutput;
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

/// Apply smoothing to all intensity traces in the chromatogram
pub fn apply_smoothing_chromatogram(
    chromatogram: &mut ChromatogramOutput,
    method: &SmoothingMethod,
) {
    if matches!(method, SmoothingMethod::None) {
        return;
    }

    // Smooth precursor intensities
    for intensities in chromatogram.precursor_intensities.iter_mut() {
        if let Some(smoothed) = apply_smoothing(intensities, method) {
            *intensities = smoothed;
        }
    }

    // Smooth fragment intensities
    for intensities in chromatogram.fragment_intensities.iter_mut() {
        if let Some(smoothed) = apply_smoothing(intensities, method) {
            *intensities = smoothed;
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
    let rt_range_ms = match tolerance.rt_range_as_milis(elution_group.rt_seconds()) {
        timsquery::OptionallyRestricted::Unrestricted => {
            timsquery::TupleRange::try_new(*ms1_rts.first().unwrap(), *ms1_rts.last().unwrap())
                .expect("Reference RTs should be sorted and valid")
        }
        timsquery::OptionallyRestricted::Restricted(r) => r,
    };
    let mut collector = ChromatogramCollector::new(elution_group.clone(), rt_range_ms, &ms1_rts)
        .map_err(|e| ViewerError::General(format!("Failed to create collector: {:?}", e)))?;

    index.add_query(&mut collector, tolerance);

    let mut output = match ChromatogramOutput::try_new(collector, &ms1_rts) {
        Ok(cmg) => cmg,
        Err(e) => {
            return Err(ViewerError::General(format!(
                "Failed to generate chromatogram output: {:?}",
                e
            )));
        }
    };

    apply_smoothing_chromatogram(&mut output, smoothing);

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
