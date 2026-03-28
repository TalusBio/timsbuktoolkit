use std::sync::Arc;

use pyo3::prelude::*;
use timscentroid::rt_mapping::RTIndex;
use timsquery::serde::IndexedPeaksHandle;
use timsquery::traits::queriable_data::QueriableData;
use timsquery::{
    ChromatogramCollector, MzMobilityStatsCollector, OptionallyRestricted, SpectralCollector,
    Tolerance,
};

use crate::chromatogram::PyChromatogramResult;
use crate::elution_group::PyElutionGroup;
use crate::iterator::PyChromatogramIterator;
use crate::spectrum::{PyMzMobilityResult, PySpectralResult};
use crate::tolerance::PyTolerance;

/// Compute the RT range in milliseconds for a chromatogram query.
///
/// If the tolerance is restricted, returns the tolerance-derived range.
/// If unrestricted, falls back to the full acquisition RT range from the
/// cycle mapping (giving a chromatogram spanning the entire run).
pub(crate) fn rt_range_ms_for_chromatogram(
    tol: &Tolerance,
    rt_seconds: f32,
    handle: &IndexedPeaksHandle,
) -> PyResult<timsquery::TupleRange<u32>> {
    match tol.rt_range_as_milis(rt_seconds) {
        OptionallyRestricted::Restricted(range) => Ok(range),
        OptionallyRestricted::Unrestricted => {
            let (start, end) = handle.ms1_cycle_mapping().range_milis();
            timsquery::TupleRange::try_new(start, end).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Empty RT range in index: {:?}",
                    e
                ))
            })
        }
    }
}

/// Resolved tolerances: either one shared or one per query.
pub(crate) enum ResolvedTolerances {
    Single(Tolerance),
    PerQuery(Vec<Tolerance>),
}

impl ResolvedTolerances {
    /// Extract from a Python object: either a PyTolerance or a list of PyTolerance.
    pub fn from_py(py: Python<'_>, obj: &PyObject, expected_len: Option<usize>) -> PyResult<Self> {
        // Try single PyTolerance first
        if let Ok(tol) = obj.extract::<PyRef<'_, PyTolerance>>(py) {
            return Ok(Self::Single(tol.inner.clone()));
        }

        // Try list of PyTolerance
        let list: Vec<PyRef<'_, PyTolerance>> = obj.extract(py).map_err(|_| {
            PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "tolerance must be a PyTolerance or a list of PyTolerance",
            )
        })?;

        if let Some(expected) = expected_len {
            if list.len() != expected {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "tolerance list length ({}) must match elution_groups length ({})",
                    list.len(),
                    expected,
                )));
            }
        }

        Ok(Self::PerQuery(
            list.iter().map(|t| t.inner.clone()).collect(),
        ))
    }

    /// Get the tolerance for query at index i.
    pub fn get(&self, i: usize) -> &Tolerance {
        match self {
            Self::Single(tol) => tol,
            Self::PerQuery(tols) => &tols[i],
        }
    }
}

/// A loaded timsTOF index, ready for querying.
///
/// Wraps either an eager (fully in-memory) or lazy (on-demand parquet)
/// indexed peaks handle.
#[pyclass]
pub struct PyTimsIndex {
    pub(crate) handle: Arc<IndexedPeaksHandle>,
}

#[pymethods]
impl PyTimsIndex {
    /// Load a timsTOF index from a .d directory or .d.idx cache.
    ///
    /// Args:
    ///     path: Path to the .d file or .d.idx cached index.
    ///     prefer_lazy: If True, prefer lazy loading for cached indexes (default: False).
    #[new]
    #[pyo3(signature = (path, prefer_lazy=false))]
    fn new(path: &str, prefer_lazy: bool) -> PyResult<Self> {
        let config = timsquery::serde::IndexLoadConfig {
            prefer_lazy,
            ..Default::default()
        };
        let handle = timsquery::serde::load_index_auto(path, Some(config))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{:?}", e)))?;
        Ok(Self {
            handle: Arc::new(handle),
        })
    }

    /// Query a single elution group and return a ChromatogramResult.
    ///
    /// Args:
    ///     elution_group: The query target.
    ///     tolerance: Search tolerances.
    ///
    /// Returns:
    ///     ChromatogramResult with precursor and fragment intensity arrays.
    fn query_chromatogram(
        &self,
        elution_group: &PyElutionGroup,
        tolerance: &PyTolerance,
    ) -> PyResult<PyChromatogramResult> {
        let eg = elution_group.inner.clone();
        let tol = &tolerance.inner;
        let rt_range_ms = rt_range_ms_for_chromatogram(tol, eg.rt_seconds(), &self.handle)?;
        let ref_rt = self.handle.ms1_cycle_mapping();

        let mut collector = ChromatogramCollector::<usize, f32>::new(eg, rt_range_ms, ref_rt)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{:?}", e)))?;

        self.handle.add_query(&mut collector, tol);

        Ok(PyChromatogramResult::new(collector))
    }

    /// Re-query into an existing ChromatogramResult, reusing its allocation.
    ///
    /// The internal arrays are reset and refilled. The Vec capacity is preserved,
    /// so repeated calls with similarly-sized elution groups avoid reallocation.
    ///
    /// Args:
    ///     result: A previously returned ChromatogramResult (mutated in place).
    ///     elution_group: The new query target.
    ///     tolerance: Search tolerances.
    fn query_chromatogram_into(
        &self,
        result: &mut PyChromatogramResult,
        elution_group: &PyElutionGroup,
        tolerance: &PyTolerance,
    ) -> PyResult<()> {
        let eg = elution_group.inner.clone();
        let tol = &tolerance.inner;
        let rt_range_ms = rt_range_ms_for_chromatogram(tol, eg.rt_seconds(), &self.handle)?;
        let ref_rt = self.handle.ms1_cycle_mapping();

        result
            .collector
            .try_reset_with(eg, rt_range_ms, ref_rt)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{:?}", e)))?;

        self.handle.add_query(&mut result.collector, tol);

        Ok(())
    }

    /// Query a single elution group for total summed intensity per ion.
    ///
    /// Returns one f32 intensity value per precursor/fragment (no RT dimension).
    ///
    /// Args:
    ///     elution_group: The query target.
    ///     tolerance: Search tolerances.
    fn query_spectrum(
        &self,
        elution_group: &PyElutionGroup,
        tolerance: &PyTolerance,
    ) -> PyResult<PySpectralResult> {
        let eg = elution_group.inner.clone();
        let tol = &tolerance.inner;
        let mut collector = SpectralCollector::<usize, f32>::new(eg);
        self.handle.add_query(&mut collector, tol);
        Ok(PySpectralResult::new(collector))
    }

    /// Query a single elution group for intensity-weighted mean m/z and mobility.
    ///
    /// Returns (weight, mean_mz, mean_mobility) per precursor/fragment.
    /// NaN values indicate no peaks were found for that ion.
    ///
    /// Args:
    ///     elution_group: The query target.
    ///     tolerance: Search tolerances.
    fn query_mz_mobility(
        &self,
        elution_group: &PyElutionGroup,
        tolerance: &PyTolerance,
    ) -> PyResult<PyMzMobilityResult> {
        let eg = elution_group.inner.clone();
        let tol = &tolerance.inner;
        let mut collector = SpectralCollector::<usize, MzMobilityStatsCollector>::new(eg);
        self.handle.add_query(&mut collector, tol);
        Ok(PyMzMobilityResult::new(collector))
    }

    /// Query multiple elution groups in parallel (via rayon).
    ///
    /// Args:
    ///     elution_groups: List of query targets.
    ///     tolerance: A single PyTolerance (shared) or a list of PyTolerance
    ///                (one per elution group, must match length).
    ///
    /// Returns:
    ///     List of ChromatogramResult, one per elution group.
    fn query_chromatograms_batch(
        &self,
        py: Python<'_>,
        elution_groups: Vec<PyRef<'_, PyElutionGroup>>,
        tolerance: PyObject,
    ) -> PyResult<Vec<PyChromatogramResult>> {
        let n = elution_groups.len();
        let tolerances = ResolvedTolerances::from_py(py, &tolerance, Some(n))?;
        let ref_rt = self.handle.ms1_cycle_mapping();

        let mut collectors: Vec<ChromatogramCollector<usize, f32>> = elution_groups
            .iter()
            .enumerate()
            .map(|(i, eg)| {
                let inner = eg.inner.clone();
                let tol = tolerances.get(i);
                let rt_range_ms =
                    rt_range_ms_for_chromatogram(tol, inner.rt_seconds(), &self.handle)?;
                ChromatogramCollector::<usize, f32>::new(inner, rt_range_ms, ref_rt).map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{:?}", e))
                })
            })
            .collect::<PyResult<Vec<_>>>()?;

        // Release GIL for the parallel query work
        let handle = &*self.handle;
        py.allow_threads(|| match &tolerances {
            ResolvedTolerances::Single(tol) => {
                handle.par_add_query_multi(
                    &mut collectors[..],
                    rayon::iter::repeat_n(tol, n),
                );
            }
            ResolvedTolerances::PerQuery(tols) => {
                handle.par_add_query_multi(&mut collectors[..], &tols[..]);
            }
        });

        Ok(collectors
            .into_iter()
            .map(PyChromatogramResult::new)
            .collect())
    }

    /// Streaming query over an iterator of elution groups.
    ///
    /// Internally reuses collector allocations and processes queries in parallel
    /// chunks via rayon. Yields lightweight ChromatogramArrays with owned numpy
    /// arrays.
    ///
    /// Args:
    ///     elution_groups: Any Python iterable of ElutionGroup objects.
    ///     tolerance: A single PyTolerance (shared) or a Python iterable of
    ///                PyTolerance (one per elution group, consumed in lockstep).
    ///     chunk_size: Number of queries per parallel batch (default: 256).
    ///
    /// Returns:
    ///     An iterator yielding ChromatogramArrays.
    #[pyo3(signature = (elution_groups, tolerance, chunk_size=None))]
    fn query_chromatograms_iter(
        &self,
        py: Python<'_>,
        elution_groups: PyObject,
        tolerance: PyObject,
        chunk_size: Option<usize>,
    ) -> PyResult<PyChromatogramIterator> {
        let eg_iter = elution_groups.call_method0(py, "__iter__")?;

        // Check if tolerance is a single PyTolerance or an iterable
        let tol_source = match tolerance.extract::<PyRef<'_, PyTolerance>>(py) {
            Ok(tol_ref) => {
                crate::iterator::ToleranceSource::Single(tol_ref.inner.clone())
            }
            Err(_) => {
                let tol_iter = tolerance.call_method0(py, "__iter__").map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "tolerance must be a PyTolerance or an iterable of PyTolerance",
                    )
                })?;
                crate::iterator::ToleranceSource::PerQuery(tol_iter)
            }
        };

        Ok(PyChromatogramIterator::new(
            Arc::clone(&self.handle),
            tol_source,
            eg_iter,
            chunk_size.unwrap_or(256),
        ))
    }

    /// Convert a retention time in seconds to the nearest cycle index.
    fn rt_seconds_to_cycle_index(&self, rt_seconds: f32) -> usize {
        let rt_ms = (rt_seconds * 1000.0) as u32;
        self.handle
            .ms1_cycle_mapping()
            .ms_to_closest_index(rt_ms)
            .index()
    }

    /// Convert a cycle index to retention time in milliseconds.
    ///
    /// Raises IndexError if the index is out of bounds.
    fn cycle_index_to_rt_ms(&self, index: u32) -> PyResult<u32> {
        use timscentroid::rt_mapping::MS1CycleIndex;
        let idx = MS1CycleIndex::new(index);
        self.handle
            .ms1_cycle_mapping()
            .rt_milis_for_index(&idx)
            .map_err(|_| {
                PyErr::new::<pyo3::exceptions::PyIndexError, _>(format!(
                    "Cycle index {} out of bounds (num_cycles={})",
                    index,
                    self.handle.ms1_cycle_mapping().len(),
                ))
            })
    }

    /// All cycle retention times in milliseconds, as a list.
    ///
    /// Index i corresponds to cycle i. Useful for building an RT axis
    /// to align with chromatogram arrays.
    #[getter]
    fn rt_values_ms(&self) -> Vec<u32> {
        let mapping = self.handle.ms1_cycle_mapping();
        (0..mapping.len() as u32)
            .map(|i| {
                use timscentroid::rt_mapping::MS1CycleIndex;
                mapping
                    .rt_milis_for_index(&MS1CycleIndex::new(i))
                    .unwrap()
            })
            .collect()
    }

    /// Total number of MS1 cycles in the acquisition.
    #[getter]
    fn num_cycles(&self) -> usize {
        self.handle.ms1_cycle_mapping().len()
    }

    /// Full acquisition RT range as (start_ms, end_ms).
    #[getter]
    fn rt_range_ms(&self) -> (u32, u32) {
        self.handle.ms1_cycle_mapping().range_milis()
    }

    /// Whether this index is lazily loaded.
    #[getter]
    fn is_lazy(&self) -> bool {
        self.handle.is_lazy()
    }

    fn __repr__(&self) -> String {
        let mode = if self.handle.is_lazy() {
            "lazy"
        } else {
            "eager"
        };
        format!("TimsIndex(mode={})", mode)
    }
}
