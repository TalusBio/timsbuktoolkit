use pyo3::prelude::*;
use timsquery::TimsElutionGroup;
use timsquery::tinyvec::tiny_vec;

/// An elution group defines a query target: one precursor and its fragments.
///
/// NOTE: Fragment labels are `usize` only in this binding. This is a deliberate
/// simplification — the Rust side is generic over `T: KeyLike` but we monomorphize
/// to `usize` here for a clean Python interface. Other key types (e.g. `IonAnnot`)
/// may be added in future versions.
#[pyclass]
#[derive(Debug, Clone)]
pub struct PyElutionGroup {
    pub(crate) inner: TimsElutionGroup<usize>,
}

#[pymethods]
impl PyElutionGroup {
    /// Create a new ElutionGroup.
    ///
    /// Args:
    ///     id: Unique identifier.
    ///     precursor_mz: Monoisotopic precursor m/z.
    ///     precursor_charge: Charge state.
    ///     rt_seconds: Expected retention time in seconds.
    ///     mobility: Expected ion mobility (1/K0).
    ///     fragment_mzs: List of fragment m/z values.
    ///     fragment_labels: List of integer labels (one per fragment, same length as fragment_mzs).
    ///     precursor_labels: Isotope offset labels (e.g. [0, 1, -1] for M0, M+1, M-1).
    #[new]
    #[pyo3(signature = (id, precursor_mz, precursor_charge, rt_seconds, mobility, fragment_mzs, fragment_labels, precursor_labels=None))]
    fn new(
        id: u64,
        precursor_mz: f64,
        precursor_charge: u8,
        rt_seconds: f32,
        mobility: f32,
        fragment_mzs: Vec<f64>,
        fragment_labels: Vec<usize>,
        precursor_labels: Option<Vec<i8>>,
    ) -> PyResult<Self> {
        let precursor_labels_tv = match precursor_labels {
            Some(labels) => labels.into_iter().collect(),
            None => tiny_vec![0i8],
        };
        let fragment_labels_tv = fragment_labels.into_iter().collect();

        let eg = TimsElutionGroup::builder()
            .id(id)
            .precursor(precursor_mz, precursor_charge)
            .mobility_ook0(mobility)
            .rt_seconds(rt_seconds)
            .fragment_mzs(fragment_mzs)
            .fragment_labels(fragment_labels_tv)
            .precursor_labels(precursor_labels_tv)
            .try_build()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{:?}", e)))?;

        Ok(Self { inner: eg })
    }

    #[getter]
    fn id(&self) -> u64 {
        self.inner.id()
    }

    #[getter]
    fn precursor_mz(&self) -> f64 {
        self.inner.precursor_mz()
    }

    #[getter]
    fn precursor_charge(&self) -> u8 {
        self.inner.precursor_charge()
    }

    #[getter]
    fn rt_seconds(&self) -> f32 {
        self.inner.rt_seconds()
    }

    #[getter]
    fn mobility(&self) -> f32 {
        self.inner.mobility_ook0()
    }

    #[getter]
    fn num_fragments(&self) -> usize {
        self.inner.fragment_count()
    }

    #[getter]
    fn num_precursors(&self) -> usize {
        self.inner.precursor_count()
    }

    fn __repr__(&self) -> String {
        format!(
            "ElutionGroup(id={}, mz={:.4}, charge={}, rt={:.1}s, mob={:.3}, frags={}, precs={})",
            self.inner.id(),
            self.inner.precursor_mz(),
            self.inner.precursor_charge(),
            self.inner.rt_seconds(),
            self.inner.mobility_ook0(),
            self.inner.fragment_count(),
            self.inner.precursor_count(),
        )
    }
}
