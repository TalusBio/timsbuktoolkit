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
mod elution_group;
mod index;
pub(crate) mod iterator;
mod numpy_utils;
mod spectrum;
mod tolerance;

use pyo3::prelude::*;

#[pymodule]
fn timsquery_pyo3(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Tolerance types
    m.add_class::<tolerance::PyMzTolerance>()?;
    m.add_class::<tolerance::PyRtTolerance>()?;
    m.add_class::<tolerance::PyMobilityTolerance>()?;
    m.add_class::<tolerance::PyQuadTolerance>()?;
    m.add_class::<tolerance::PyTolerance>()?;

    // Query definition
    m.add_class::<elution_group::PyElutionGroup>()?;

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
