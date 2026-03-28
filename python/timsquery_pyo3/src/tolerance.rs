use pyo3::prelude::*;

#[pyclass(frozen)]
#[derive(Debug, Clone)]
pub struct PyMzTolerance {
    pub(crate) inner: timsquery::models::tolerance::MzTolerance,
}

#[pymethods]
impl PyMzTolerance {
    /// Create a parts-per-million tolerance: ±(low_ppm, high_ppm).
    #[staticmethod]
    fn ppm(low: f64, high: f64) -> Self {
        Self {
            inner: timsquery::models::tolerance::MzTolerance::Ppm((low, high)),
        }
    }

    /// Create an absolute (dalton) tolerance: ±(low_da, high_da).
    #[staticmethod]
    fn absolute(low: f64, high: f64) -> Self {
        Self {
            inner: timsquery::models::tolerance::MzTolerance::Absolute((low, high)),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(frozen)]
#[derive(Debug, Clone)]
pub struct PyRtTolerance {
    pub(crate) inner: timsquery::models::tolerance::RtTolerance,
}

#[pymethods]
impl PyRtTolerance {
    /// Create a fixed-minutes tolerance: ±(low_min, high_min).
    #[staticmethod]
    fn minutes(low: f32, high: f32) -> Self {
        Self {
            inner: timsquery::models::tolerance::RtTolerance::Minutes((low, high)),
        }
    }

    /// Create a percentage tolerance: ±(low_pct, high_pct).
    #[staticmethod]
    fn pct(low: f32, high: f32) -> Self {
        Self {
            inner: timsquery::models::tolerance::RtTolerance::Pct((low, high)),
        }
    }

    /// No restriction on retention time.
    #[staticmethod]
    fn unrestricted() -> Self {
        Self {
            inner: timsquery::models::tolerance::RtTolerance::Unrestricted,
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(frozen)]
#[derive(Debug, Clone)]
pub struct PyMobilityTolerance {
    pub(crate) inner: timsquery::models::tolerance::MobilityTolerance,
}

#[pymethods]
impl PyMobilityTolerance {
    /// Create an absolute (1/K0) tolerance: ±(low, high).
    #[staticmethod]
    fn absolute(low: f32, high: f32) -> Self {
        Self {
            inner: timsquery::models::tolerance::MobilityTolerance::Absolute((low, high)),
        }
    }

    /// Create a percentage tolerance: ±(low_pct, high_pct).
    #[staticmethod]
    fn pct(low: f32, high: f32) -> Self {
        Self {
            inner: timsquery::models::tolerance::MobilityTolerance::Pct((low, high)),
        }
    }

    /// No restriction on ion mobility.
    #[staticmethod]
    fn unrestricted() -> Self {
        Self {
            inner: timsquery::models::tolerance::MobilityTolerance::Unrestricted,
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(frozen)]
#[derive(Debug, Clone)]
pub struct PyQuadTolerance {
    pub(crate) inner: timsquery::models::tolerance::QuadTolerance,
}

#[pymethods]
impl PyQuadTolerance {
    /// Create an absolute (dalton) quadrupole tolerance: ±(low_da, high_da).
    #[staticmethod]
    fn absolute(low: f32, high: f32) -> Self {
        Self {
            inner: timsquery::models::tolerance::QuadTolerance::Absolute((low, high)),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(frozen)]
#[derive(Debug, Clone)]
pub struct PyTolerance {
    pub(crate) inner: timsquery::Tolerance,
}

#[pymethods]
impl PyTolerance {
    /// Construct a Tolerance from per-dimension tolerance objects.
    #[new]
    fn new(
        mz: &PyMzTolerance,
        rt: &PyRtTolerance,
        mobility: &PyMobilityTolerance,
        quad: &PyQuadTolerance,
    ) -> Self {
        Self {
            inner: timsquery::Tolerance {
                ms: mz.inner.clone(),
                rt: rt.inner.clone(),
                mobility: mobility.inner.clone(),
                quad: quad.inner.clone(),
            },
        }
    }

    /// Default tolerance: 20 ppm m/z, ±5 min RT, 3% mobility, 0.1 Da quad.
    #[staticmethod]
    fn default() -> Self {
        Self {
            inner: timsquery::Tolerance::default(),
        }
    }

    /// Return a new Tolerance with the m/z tolerance replaced.
    fn with_mz(&self, mz: &PyMzTolerance) -> Self {
        Self {
            inner: timsquery::Tolerance {
                ms: mz.inner.clone(),
                ..self.inner.clone()
            },
        }
    }

    /// Return a new Tolerance with the RT tolerance replaced.
    fn with_rt(&self, rt: &PyRtTolerance) -> Self {
        Self {
            inner: timsquery::Tolerance {
                rt: rt.inner.clone(),
                ..self.inner.clone()
            },
        }
    }

    /// Return a new Tolerance with the mobility tolerance replaced.
    fn with_mobility(&self, mobility: &PyMobilityTolerance) -> Self {
        Self {
            inner: timsquery::Tolerance {
                mobility: mobility.inner.clone(),
                ..self.inner.clone()
            },
        }
    }

    /// Return a new Tolerance with the quad tolerance replaced.
    fn with_quad(&self, quad: &PyQuadTolerance) -> Self {
        Self {
            inner: timsquery::Tolerance {
                quad: quad.inner.clone(),
                ..self.inner.clone()
            },
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
