use std::collections::VecDeque;
use std::sync::Arc;

use pyo3::prelude::*;
use timsquery::serde::IndexedPeaksHandle;
use timsquery::traits::queriable_data::QueriableData;
use timsquery::{
    ChromatogramCollector,
    Tolerance,
};

use crate::elution_group::PyElutionGroup;
use crate::index::rt_range_ms_for_chromatogram;
use crate::numpy_utils::array2d_to_numpy;
use crate::tolerance::PyTolerance;

/// Source of tolerances: either one shared or one per query from a Python iterator.
pub enum ToleranceSource {
    Single(Tolerance),
    PerQuery(PyObject),
}

/// Lightweight result yielded by the streaming iterator.
///
/// Owns materialized numpy arrays and metadata. The iterator's internal
/// collector pool is never exposed — it reuses Rust-side buffers across chunks.
#[pyclass(frozen)]
pub struct PyChromatogramArrays {
    #[pyo3(get)]
    id: u64,
    precursor_intensities: PyObject,
    fragment_intensities: PyObject,
    #[pyo3(get)]
    precursor_labels: Vec<(i8, f64)>,
    #[pyo3(get)]
    fragment_labels: Vec<(usize, f64)>,
    #[pyo3(get)]
    rt_range_ms: (u32, u32),
    #[pyo3(get)]
    num_cycles: usize,
}

#[pymethods]
impl PyChromatogramArrays {
    #[getter]
    fn precursor_intensities<'py>(&self, py: Python<'py>) -> Bound<'py, PyAny> {
        self.precursor_intensities.clone_ref(py).into_bound(py)
    }

    #[getter]
    fn fragment_intensities<'py>(&self, py: Python<'py>) -> Bound<'py, PyAny> {
        self.fragment_intensities.clone_ref(py).into_bound(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "ChromatogramArrays(id={}, precursors={}, fragments={}, cycles={})",
            self.id,
            self.precursor_labels.len(),
            self.fragment_labels.len(),
            self.num_cycles,
        )
    }
}

fn extract_arrays(
    py: Python<'_>,
    collector: &ChromatogramCollector<usize, f32>,
) -> PyResult<PyChromatogramArrays> {
    let prec_np = array2d_to_numpy(py, &collector.precursors.arr)?;
    let frag_np = array2d_to_numpy(py, &collector.fragments.arr)?;
    let rt = collector.rt_range_milis();

    Ok(PyChromatogramArrays {
        id: collector.id,
        precursor_intensities: prec_np.into_any().unbind(),
        fragment_intensities: frag_np.into_any().unbind(),
        precursor_labels: collector
            .precursors
            .mz_order
            .iter()
            .map(|(k, mz)| (*k, *mz))
            .collect(),
        fragment_labels: collector
            .fragments
            .mz_order
            .iter()
            .map(|(k, mz)| (*k, *mz))
            .collect(),
        rt_range_ms: (rt.start(), rt.end()),
        num_cycles: collector.num_cycles(),
    })
}

/// Streaming chromatogram iterator with internal collector reuse.
#[pyclass]
pub struct PyChromatogramIterator {
    handle: Arc<IndexedPeaksHandle>,
    tol_source: ToleranceSource,
    eg_source: PyObject,
    pool: Vec<ChromatogramCollector<usize, f32>>,
    chunk_tolerances: Vec<Tolerance>,
    buffer: VecDeque<PyChromatogramArrays>,
    chunk_size: usize,
    exhausted: bool,
}

impl PyChromatogramIterator {
    pub fn new(
        handle: Arc<IndexedPeaksHandle>,
        tol_source: ToleranceSource,
        eg_source: PyObject,
        chunk_size: usize,
    ) -> Self {
        Self {
            handle,
            tol_source,
            eg_source,
            pool: Vec::with_capacity(chunk_size),
            chunk_tolerances: Vec::with_capacity(chunk_size),
            buffer: VecDeque::with_capacity(chunk_size),
            chunk_size,
            exhausted: false,
        }
    }

    fn fill_buffer(&mut self, py: Python<'_>) -> PyResult<()> {
        let ref_rt = self.handle.ms1_cycle_mapping();
        let mut n_this_chunk = 0;
        self.chunk_tolerances.clear();

        for i in 0..self.chunk_size {
            let next_result = self.eg_source.call_method0(py, "__next__");
            match next_result {
                Ok(obj) => {
                    let eg_ref: PyRef<'_, PyElutionGroup> = obj.extract(py)?;
                    let eg = eg_ref.inner.clone();

                    let tol = match &self.tol_source {
                        ToleranceSource::Single(t) => t.clone(),
                        ToleranceSource::PerQuery(iter_obj) => {
                            let tol_obj = iter_obj.call_method0(py, "__next__").map_err(|e| {
                                if e.is_instance_of::<pyo3::exceptions::PyStopIteration>(py) {
                                    PyErr::new::<pyo3::exceptions::PyValueError, _>(
                                        "tolerance iterator exhausted before elution_groups iterator",
                                    )
                                } else {
                                    e
                                }
                            })?;
                            let tol_ref: PyRef<'_, PyTolerance> = tol_obj.extract(py)?;
                            tol_ref.inner.clone()
                        }
                    };

                    let rt_range_ms =
                        rt_range_ms_for_chromatogram(&tol, eg.rt_seconds(), &self.handle)?;

                    if i < self.pool.len() {
                        self.pool[i]
                            .try_reset_with(&eg, rt_range_ms, ref_rt)
                            .map_err(|e| {
                                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{e:?}"))
                            })?;
                    } else {
                        let collector =
                            ChromatogramCollector::<usize, f32>::new(&eg, rt_range_ms, ref_rt)
                                .map_err(|e| {
                                    PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                                        "{e:?}"
                                    ))
                                })?;
                        self.pool.push(collector);
                    }
                    self.chunk_tolerances.push(tol);
                    n_this_chunk += 1;
                }
                Err(err) if err.is_instance_of::<pyo3::exceptions::PyStopIteration>(py) => {
                    self.exhausted = true;
                    break;
                }
                Err(err) => return Err(err),
            }
        }

        if n_this_chunk == 0 {
            return Ok(());
        }

        // Release GIL for the parallel query work
        let pool_slice = &mut self.pool[..n_this_chunk];
        let tol_slice = &self.chunk_tolerances[..];
        let handle = &self.handle;
        py.allow_threads(|| {
            handle.par_add_query_multi(pool_slice, tol_slice);
        });

        for collector in &self.pool[..n_this_chunk] {
            self.buffer.push_back(extract_arrays(py, collector)?);
        }

        Ok(())
    }
}

#[pymethods]
impl PyChromatogramIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PyChromatogramArrays>> {
        if let Some(arrays) = self.buffer.pop_front() {
            return Ok(Some(arrays));
        }

        if self.exhausted {
            return Ok(None);
        }

        self.fill_buffer(py)?;
        Ok(self.buffer.pop_front())
    }
}
