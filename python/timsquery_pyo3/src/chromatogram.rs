use numpy::PyArray2;
use pyo3::prelude::*;
use timsquery::ChromatogramCollector;

use crate::numpy_utils::array2d_to_numpy;

/// Result of a chromatogram query.
///
/// Contains 2D numpy arrays for precursor and fragment intensity traces
/// across retention time cycles. Shaped (n_ions, n_cycles).
#[pyclass]
pub struct PyChromatogramResult {
    pub(crate) collector: ChromatogramCollector<usize, f32>,
}

impl PyChromatogramResult {
    pub fn new(collector: ChromatogramCollector<usize, f32>) -> Self {
        Self { collector }
    }
}

#[pymethods]
impl PyChromatogramResult {
    #[getter]
    fn precursor_intensities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        array2d_to_numpy(py, &self.collector.precursors.arr)
    }

    #[getter]
    fn fragment_intensities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        array2d_to_numpy(py, &self.collector.fragments.arr)
    }

    #[getter]
    fn precursor_labels(&self) -> Vec<(i8, f64)> {
        self.collector
            .precursors
            .mz_order
            .iter()
            .map(|(k, mz)| (*k, *mz))
            .collect()
    }

    #[getter]
    fn fragment_labels(&self) -> Vec<(usize, f64)> {
        self.collector
            .fragments
            .mz_order
            .iter()
            .map(|(k, mz)| (*k, *mz))
            .collect()
    }

    #[getter]
    fn rt_range_ms(&self) -> (u32, u32) {
        let r = self.collector.rt_range_milis();
        (r.start(), r.end())
    }

    #[getter]
    fn num_cycles(&self) -> usize {
        self.collector.num_cycles()
    }

    #[getter]
    fn id(&self) -> u64 {
        self.collector.eg.id()
    }

    fn __repr__(&self) -> String {
        format!(
            "ChromatogramResult(id={}, precursors={}x{}, fragments={}x{}, rt_ms=({}, {}))",
            self.collector.eg.id(),
            self.collector.precursors.arr.nrows(),
            self.collector.precursors.arr.ncols(),
            self.collector.fragments.arr.nrows(),
            self.collector.fragments.arr.ncols(),
            self.collector.rt_range_milis().start(),
            self.collector.rt_range_milis().end(),
        )
    }
}
