use numpy::{
    PyArray1,
    PyArray2,
    PyArrayMethods,
};
use pyo3::prelude::*;
use timsquery::Array2D;

/// Copy an Array2D<f32> into a new numpy 2D array (safe, no unsafe blocks).
pub fn array2d_to_numpy<'py>(
    py: Python<'py>,
    arr: &Array2D<f32>,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    let flat = PyArray1::from_slice(py, arr.as_flat_slice());
    flat.reshape((arr.nrows(), arr.ncols()))
}
