use pyo3::prelude::*;
use timsquery::{
    MzMobilityStatsCollector,
    SpectralCollector,
};

/// Result of a spectral query — total summed intensity per ion.
///
/// Each precursor/fragment gets a single f32 intensity value (summed
/// across all matching peaks within the tolerance window).
///
/// NOTE: Uses `usize` fragment keys and `f32` intensities.
#[pyclass(frozen)]
pub struct PySpectralResult {
    collector: SpectralCollector<usize, f32>,
}

impl PySpectralResult {
    pub fn new(collector: SpectralCollector<usize, f32>) -> Self {
        Self { collector }
    }
}

#[pymethods]
impl PySpectralResult {
    /// Total intensity per precursor isotope.
    #[getter]
    fn precursor_intensities(&self) -> Vec<f32> {
        self.collector
            .iter_precursors()
            .map(|(_, val)| *val)
            .collect()
    }

    /// Total intensity per fragment ion.
    #[getter]
    fn fragment_intensities(&self) -> Vec<f32> {
        self.collector
            .iter_fragments()
            .map(|(_, val)| *val)
            .collect()
    }

    /// List of (isotope_label, mz) tuples for each precursor.
    #[getter]
    fn precursor_labels(&self) -> Vec<(i8, f64)> {
        self.collector
            .iter_precursors()
            .map(|((label, mz), _)| (label, mz))
            .collect()
    }

    /// List of (fragment_label, mz) tuples for each fragment.
    #[getter]
    fn fragment_labels(&self) -> Vec<(usize, f64)> {
        self.collector
            .iter_fragments()
            .map(|((label, mz), _)| (*label, *mz))
            .collect()
    }

    /// The elution group id.
    #[getter]
    fn id(&self) -> u64 {
        self.collector.eg.id()
    }

    fn __repr__(&self) -> String {
        let prec_sum: f32 = self.collector.iter_precursors().map(|(_, v)| v).sum();
        let frag_sum: f32 = self.collector.iter_fragments().map(|(_, v)| v).sum();
        format!(
            "SpectralResult(id={}, precursors={}, fragments={}, prec_total={:.1}, frag_total={:.1})",
            self.collector.eg.id(),
            self.collector.eg.precursor_count(),
            self.collector.eg.fragment_count(),
            prec_sum,
            frag_sum,
        )
    }
}

/// Stats for a single ion from MzMobilityStatsCollector.
///
/// Exposed as a tuple: (weight, mean_mz, mean_mobility).
/// If no peaks were found, mean_mz and mean_mobility are NaN.
fn stats_to_tuple(stats: &MzMobilityStatsCollector) -> (f64, f64, f64) {
    (
        stats.weight(),
        stats.mean_mz().unwrap_or(f64::NAN),
        stats.mean_mobility().unwrap_or(f64::NAN),
    )
}

/// Result of an m/z + mobility stats query.
///
/// Each precursor/fragment gets intensity-weighted running statistics:
///   - weight: total accumulated intensity
///   - mean_mz: intensity-weighted mean m/z
///   - mean_mobility: intensity-weighted mean ion mobility (1/K0)
///
/// Stats are returned as (weight, mean_mz, mean_mobility) tuples.
/// NaN values indicate no peaks were found for that ion.
///
/// NOTE: Uses `usize` fragment keys.
#[pyclass(frozen)]
pub struct PyMzMobilityResult {
    collector: SpectralCollector<usize, MzMobilityStatsCollector>,
}

impl PyMzMobilityResult {
    pub fn new(collector: SpectralCollector<usize, MzMobilityStatsCollector>) -> Self {
        Self { collector }
    }
}

#[pymethods]
impl PyMzMobilityResult {
    /// Stats per precursor isotope: list of (weight, mean_mz, mean_mobility).
    #[getter]
    fn precursor_stats(&self) -> Vec<(f64, f64, f64)> {
        self.collector
            .iter_precursors()
            .map(|(_, stats)| stats_to_tuple(stats))
            .collect()
    }

    /// Stats per fragment ion: list of (weight, mean_mz, mean_mobility).
    #[getter]
    fn fragment_stats(&self) -> Vec<(f64, f64, f64)> {
        self.collector
            .iter_fragments()
            .map(|(_, stats)| stats_to_tuple(stats))
            .collect()
    }

    /// List of (isotope_label, mz) tuples for each precursor.
    #[getter]
    fn precursor_labels(&self) -> Vec<(i8, f64)> {
        self.collector
            .iter_precursors()
            .map(|((label, mz), _)| (label, mz))
            .collect()
    }

    /// List of (fragment_label, mz) tuples for each fragment.
    #[getter]
    fn fragment_labels(&self) -> Vec<(usize, f64)> {
        self.collector
            .iter_fragments()
            .map(|((label, mz), _)| (*label, *mz))
            .collect()
    }

    /// The elution group id.
    #[getter]
    fn id(&self) -> u64 {
        self.collector.eg.id()
    }

    fn __repr__(&self) -> String {
        format!(
            "MzMobilityResult(id={}, precursors={}, fragments={})",
            self.collector.eg.id(),
            self.collector.eg.precursor_count(),
            self.collector.eg.fragment_count(),
        )
    }
}
