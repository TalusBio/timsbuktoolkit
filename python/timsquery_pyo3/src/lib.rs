//! PyO3 bindings for timsquery.
//!
//! This crate exposes a Python-friendly interface to the timsquery library
//! for querying timsTOF mass spectrometry data.
//!
//! # Key simplifications
//!
//! - **Fragment key type is `usize` only.** The Rust library is generic over
//!   `T: KeyLike`, but this binding monomorphizes to `usize` for simplicity.
//!   Other key types (e.g. `IonAnnot`) may be added in future versions.
//!
//! - **Intensity type is `f32` only.** Chromatogram and spectral (f32) results
//!   use `f32` intensities. The MzMobility variant uses `MzMobilityStatsCollector`.

mod chromatogram;
mod index;
pub(crate) mod iterator;
mod numpy_utils;
mod spectrum;
mod target;
mod tolerance;

use pyo3::prelude::*;
use timsquery::models::OwnedSourceId;

/// A source id keeps the shape the library gave it, so Python sees an `int`
/// for a numeric id and a `str` for a text one (DIA-NN's
/// `transition_group_id`) rather than one coerced into the other.
pub(crate) fn source_id_to_py<'py>(
    py: Python<'py>,
    id: &OwnedSourceId,
) -> PyResult<Bound<'py, PyAny>> {
    use pyo3::IntoPyObject;
    match id {
        OwnedSourceId::Numeric(n) => Ok(n.into_pyobject(py)?.into_any()),
        OwnedSourceId::Text(s) => Ok(s.into_pyobject(py)?.into_any()),
    }
}

#[pymodule]
fn timsquery_pyo3(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Tolerance types
    m.add_class::<tolerance::PyMzTolerance>()?;
    m.add_class::<tolerance::PyRtTolerance>()?;
    m.add_class::<tolerance::PyMobilityTolerance>()?;
    m.add_class::<tolerance::PyQuadTolerance>()?;
    m.add_class::<tolerance::PyTolerance>()?;

    // Query definition
    m.add_class::<target::PyTarget>()?;

    // Index
    m.add_class::<index::PyTimsIndex>()?;

    // Chromatogram results
    m.add_class::<chromatogram::PyChromatogramResult>()?;

    // Spectral results
    m.add_class::<spectrum::PySpectralResult>()?;
    m.add_class::<spectrum::PyMzMobilityResult>()?;

    // Streaming iterator types
    m.add_class::<iterator::PyChromatogramArrays>()?;
    m.add_class::<iterator::PyChromatogramIterator>()?;

    Ok(())
}
