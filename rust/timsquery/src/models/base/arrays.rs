use std::fmt::Debug;

use serde::Serialize;
use serde::ser::SerializeSeq;

use crate::errors::DataProcessingError;

use std::ops::Range;

pub trait ArrayElement:
    Clone
    + Copy
    + Default
    + std::ops::Mul
    + std::fmt::Display
    + std::fmt::Debug
    + std::ops::Add<Self, Output = Self>
    + std::ops::AddAssign<Self>
{
}

/// Blanket trait implementation on elements that
/// can be used as a value in the array.
impl<
    T: Clone
        + Copy
        + Default
        + std::ops::Mul
        + std::fmt::Display
        + std::fmt::Debug
        + std::ops::Add<T, Output = T>
        + std::ops::AddAssign<T>,
> ArrayElement for T
{
}

/// Implements a way to represent an array of
/// dimensions x-y that will be later used to
/// implement an mz-major and a rt-major array
/// representation.
///
/// Simple 2D array
///
/// `values` is a flattened array of values
/// `n_col` is the number of values in each row
/// `n_row` is the number of rows
///
/// Note on memory layout:
///
/// Values that belong to the same row are adjacent
/// in memory.
#[derive(Debug, Clone, PartialEq)]
pub struct Array2D<T: ArrayElement> {
    pub(super) values: Vec<T>,
    pub(super) n_col: usize,
    pub(super) n_row: usize,
}

struct Array2DView<'a, T: ArrayElement> {
    data: &'a Array2D<T>,
    col_range: Range<usize>,
    row_range: Range<usize>,
}

impl<T: Serialize + ArrayElement> Serialize for Array2DView<'_, T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(self.row_range.len()))?;
        for (ir, row) in self.data.iter_rows().enumerate() {
            if !self.row_range.contains(&ir) {
                continue;
            }
            let row_slice = &row[self.col_range.clone()];
            seq.serialize_element(row_slice)?;
        }
        seq.end()
    }
}

impl<T: Serialize + ArrayElement> Serialize for Array2D<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(self.n_row))?;
        for row in self.iter_rows() {
            seq.serialize_element(row)?;
        }
        seq.end()
    }
}

impl<T: ArrayElement> Array2D<T> {
    pub fn new<S: AsRef<[T]>, C: AsRef<[S]>>(values: C) -> Result<Array2D<T>, DataProcessingError> {
        let nrows = values.as_ref().len();
        if nrows == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        let ncols = values.as_ref()[0].as_ref().len();
        if ncols == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }

        let expected_size = nrows * ncols;
        let values: Vec<T> = values
            .as_ref()
            .iter()
            .flat_map(|x| x.as_ref())
            .cloned()
            .collect();

        if values.len() != expected_size {
            return Err(DataProcessingError::ExpectedVectorSameLength);
        }

        Ok(Array2D {
            values,
            n_col: ncols,
            n_row: nrows,
        })
    }

    pub fn new_transposed<S: AsRef<[T]>, C: AsRef<[S]>>(
        values: C,
    ) -> Result<Array2D<T>, DataProcessingError> {
        let ncols = values.as_ref().len();
        if ncols == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        let nrows = values.as_ref()[0].as_ref().len();
        if nrows == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }

        let expected_size = nrows * ncols;
        let mut out_values = vec![None; expected_size];

        for (ci, col) in values.as_ref().iter().enumerate() {
            if col.as_ref().len() != nrows {
                return Err(DataProcessingError::ExpectedVectorSameLength);
            }
            for (ri, val) in col.as_ref().iter().enumerate() {
                let idx = ri * ncols + ci; // Changed indexing for row-major order
                out_values[idx] = Some(*val);
            }
        }

        let unwrapped_values: Vec<T> = out_values.into_iter().flatten().collect();

        Ok(Array2D {
            values: unwrapped_values,
            n_col: ncols,
            n_row: nrows,
        })
    }

    pub fn from_flat_vector(
        values: Vec<T>,
        nrows: usize,
        ncols: usize,
    ) -> Result<Array2D<T>, DataProcessingError> {
        if values.len() != nrows * ncols {
            return Err(DataProcessingError::ExpectedVectorSameLength);
        }
        Ok(Array2D {
            values,
            n_col: ncols,
            n_row: nrows,
        })
    }

    /// Apply a function to each row of the array
    ///
    /// Example:
    /// ```
    /// use timsquery::Array2D;
    /// let array = Array2D::new(vec![vec![1, 2, 3], vec![4, 5, 6]]).unwrap();
    /// let result: Vec<u32> = array.row_apply(|x| x.iter().sum()).collect();
    /// assert_eq!(result, vec![6, 15]);
    ///
    /// let array = Array2D::new_transposed(vec![vec![1, 2, 3], vec![4, 5, 6]]).unwrap();
    /// let result: Vec<u32> = array.row_apply(|x| x.iter().sum()).collect();
    /// assert_eq!(result, vec![5, 7, 9]);
    /// ```
    pub fn row_apply<'a: 'b, 'b, W, F: FnMut(&[T]) -> W + 'b>(
        &'a self,
        f: F,
    ) -> impl Iterator<Item = W> + 'b {
        self.iter_rows().map(f)
    }

    pub fn iter_mut_rows(&mut self) -> impl Iterator<Item = &mut [T]> {
        self.values.chunks_mut(self.n_col)
    }

    pub fn iter_rows(&self) -> impl Iterator<Item = &[T]> {
        self.values.chunks(self.n_col)
    }

    /// Apply a function to each pair of rows of the array
    ///
    /// For example if I want to calculate the MAE between each row of the array
    /// I can do:
    /// ```
    /// use timsquery::Array2D;
    /// let array: Array2D<f64> = Array2D::new(vec![vec![1., 2., 3.], vec![4., 5., 6.]]).unwrap();
    /// let result: Vec<f64> = array.outer_row_apply(|x, y| (x.iter().zip(y.iter()).map(|(a, b)| (a - b).abs() as f64).sum::<f64>()) / x.len() as f64);
    /// assert_eq!(result, vec![3.0]);
    ///
    /// let array: Array2D<f64> = Array2D::new_transposed(vec![vec![1.,2., 3.], vec![4., 5., 6.]]).unwrap();
    /// let result: Vec<f64> = array.outer_row_apply(|x, y| (x.iter().zip(y.iter()).map(|(a, b)| (a - b).abs() as f64).sum::<f64>()) / x.len() as f64);
    /// assert_eq!(result, vec![1.0, 2.0, 1.0]);
    /// ```
    pub fn outer_row_apply<W, F: FnMut(&[T], &[T]) -> W>(&self, mut f: F) -> Vec<W> {
        let mut result = Vec::new();
        for i in 0..self.n_row {
            let row = self.get_row(i).expect("Failed to get row, malformed array");
            for j in 0..self.n_row {
                if j >= i {
                    continue;
                }
                let other_row = self.get_row(j).expect("Failed to get row, malformed array");
                result.push(f(row, other_row));
            }
        }
        result
    }

    pub fn insert(&mut self, row_idx: usize, col_idx: usize, value: T) {
        let idx = row_idx * self.n_col + col_idx;
        match self.values.get_mut(idx) {
            Some(v) => *v = value,
            None => panic!(
                "Index out of bounds ({}/{}): row {}/{}, col {}/{}",
                idx,
                self.values.len(),
                row_idx,
                self.n_row,
                col_idx,
                self.n_col,
            ),
        };
    }

    fn get_row_limits(&self, index: usize) -> Option<Range<usize>> {
        let start = index * self.n_col;
        let end = start + self.n_col;
        if end > self.values.len() || start >= self.values.len() {
            return None;
        }
        Some(start..end)
    }

    pub fn get_row(&self, index: usize) -> Option<&[T]> {
        Some(&self.values[self.get_row_limits(index)?])
    }

    pub fn get_row_mut(&mut self, index: usize) -> Result<&mut [T], DataProcessingError> {
        let range = self
            .get_row_limits(index)
            .ok_or(DataProcessingError::IndexOutOfBoundsError(index))?;
        Ok(&mut self.values[range])
    }

    pub fn iter_column(&self, index: usize) -> impl '_ + Iterator<Item = T> {
        let out = self.iter_rows().map(move |row| row[index]);
        out
    }

    pub fn try_swap_rows(&mut self, row1: usize, row2: usize) -> Result<(), DataProcessingError> {
        let range_1 = self
            .get_row_limits(row1)
            .ok_or(DataProcessingError::IndexOutOfBoundsError(row1))?;
        let range_2 = self
            .get_row_limits(row2)
            .ok_or(DataProcessingError::IndexOutOfBoundsError(row2))?;
        for (i, j) in range_1.zip(range_2) {
            self.values.swap(i, j);
        }
        Ok(())
    }

    pub fn try_replace_row_with(
        &mut self,
        row_idx: usize,
        row: &[T],
    ) -> Result<(), DataProcessingError> {
        let range = self
            .get_row_limits(row_idx)
            .ok_or(DataProcessingError::IndexOutOfBoundsError(row_idx))?;
        self.values[range].copy_from_slice(row);
        Ok(())
    }

    pub fn nrows(&self) -> usize {
        self.n_row
    }

    pub fn ncols(&self) -> usize {
        self.n_col
    }

    pub fn transpose_clone(&self) -> Array2D<T> {
        // Swap major_dim and minor_dim
        let col_dim = self.n_col;
        let row_dim = self.n_row;

        let mut result = vec![self.values[0]; row_dim * col_dim];

        for (i, crow) in self.iter_rows().enumerate() {
            for (j, v) in crow.iter().enumerate() {
                let idx = (j * row_dim) + i;
                result[idx] = *v;
            }
        }

        Array2D {
            values: result,
            n_col: row_dim,
            n_row: col_dim,
        }
    }

    pub fn reset_with_value(&mut self, ncols: usize, nrows: usize, value: T) {
        self.values.clear();
        self.values.resize(ncols * nrows, value);
        self.n_col = ncols;
        self.n_row = nrows;
    }
}

impl<T: ArrayElement + std::ops::Mul<T, Output = T>> Array2D<T> {
    /// RowConvolve
    /// Apply a a convolution on each row of the array separately.
    /// It is equivalent to applying the passed kernel on each row of the array.
    /// and padding with zeros.
    pub fn row_convolve(&self, kernel: &[T], default_value: T) -> Array2D<T> {
        let mut result = vec![default_value; self.n_row * self.n_col];
        let offset_size = (kernel.len() - 1) / 2;

        for i in 0..self.n_row {
            let row_offset = i * (self.n_col);
            let row = self.get_row(i).expect("Failed to get row, malformed array");
            row.windows(kernel.len())
                .enumerate()
                .for_each(|(ii, window)| {
                    window
                        .iter()
                        .zip(kernel.iter())
                        .map(|(&a, &b)| a * b)
                        .for_each(|prod| {
                            result[ii + offset_size + row_offset] += prod;
                        });
                })
        }

        Array2D::from_flat_vector(result, self.nrows(), self.ncols()).unwrap()
    }

    pub fn convolve_fold(
        &self,
        kernel: &[T],
        default_value: T,
        fold_func: impl Fn(T, T) -> T,
    ) -> Vec<T> {
        let mut result = vec![default_value; self.n_col];
        let offset_size = (kernel.len() - 1) / 2;

        for i in 0..self.n_row {
            let row = self.get_row(i).expect("Failed to get row, malformed array");
            row.windows(kernel.len())
                .enumerate()
                .for_each(|(ii, window)| {
                    window
                        .iter()
                        .zip(kernel.iter())
                        .map(|(&a, &b)| a * b)
                        .for_each(|prod| {
                            result[ii + offset_size] = fold_func(result[ii + offset_size], prod);
                        });
                })
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_array2d_new() {
        // Test creating a 2x3 array
        let values = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let array = Array2D::new(&values).unwrap();

        // Check dimensions
        assert_eq!(array.n_col, 3); // columns
        assert_eq!(array.n_row, 2); // rows

        // Check memory layout - values in same row should be adjacent
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn test_array2d_new_transposed() {
        // Test creating a 3x2 array from columns
        let columns = vec![
            vec![1, 4], // first column
            vec![2, 5], // second column
            vec![3, 6], // third column
        ];
        let array = Array2D::new_transposed(&columns).unwrap();

        // Check dimensions
        assert_eq!(array.n_col, 3); // columns
        assert_eq!(array.n_row, 2); // rows

        // Check memory layout - values should be arranged row-major
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn test_array2d_error_handling() {
        // Test with inconsistent row lengths
        let invalid_values = vec![
            vec![1, 2, 3],
            vec![4, 5], // Missing one value
        ];
        assert!(Array2D::new(&invalid_values).is_err());

        // Test with empty array
        let empty_values: Vec<Vec<i32>> = vec![];
        assert!(Array2D::new(&empty_values).is_err());

        // Test transposed with inconsistent column lengths
        let invalid_columns = vec![
            vec![1, 4],
            vec![2], // Missing one value
            vec![3, 6],
        ];
        assert!(Array2D::new_transposed(&invalid_columns).is_err());
    }

    #[test]
    fn test_array2d_large() {
        // Test with a larger array to verify memory efficiency
        let size = 1000;
        let values: Vec<Vec<f64>> = (0..size)
            .map(|i| (0..size).map(|j| (i * size + j) as f64).collect())
            .collect();

        let array = Array2D::new(&values).unwrap();

        // Verify dimensions
        assert_eq!(array.n_col, size);
        assert_eq!(array.n_row, size);

        // Verify first and last values
        assert_eq!(array.values[0], 0.0);
        assert_eq!(array.values[size * size - 1], (size * size - 1) as f64);
    }

    #[test]
    fn test_array2d_transpose() {
        let values = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let array = Array2D::new(&values).unwrap();

        // Check dimensions
        assert_eq!(array.n_col, 3); // columns
        assert_eq!(array.n_row, 2); // rows

        // Check memory layout - values in same row should be adjacent
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);

        let transposed = array.transpose_clone();

        // Check dimensions
        assert_eq!(transposed.n_col, 2); // columns
        assert_eq!(transposed.n_row, 3); // rows

        // Check memory layout - values should be arranged row-major
        assert_eq!(transposed.values, vec![1, 4, 2, 5, 3, 6]);
    }

    #[test]
    fn test_convolve() {
        let values = vec![vec![1, 2, 3, 4, 5, 6], vec![4, 5, 6, 7, 8, 9]];
        let array = Array2D::new(&values).unwrap();
        let array2 = Array2D {
            values: vec![1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 9],
            n_col: 6,
            n_row: 2,
        };
        assert_eq!(array, array2);

        let convolved = array.row_convolve(&[1, 1, 1], 0);

        // 6 = 1 + 2 + 3
        // 9 = 2 + 3 + 4
        let expect = Array2D {
            values: vec![0, 6, 9, 12, 15, 0, 0, 15, 18, 21, 24, 0],
            n_col: 6,
            n_row: 2,
        };
        assert_eq!(convolved, expect);
    }

    #[test]
    fn test_convolve_reduce() {
        let values = vec![vec![1, 2, 3, 4, 5, 6], vec![4, 5, 6, 7, 8, 9]];
        let array = Array2D::new(&values).unwrap();
        let array2 = Array2D {
            values: vec![1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 9],
            n_col: 6,
            n_row: 2,
        };
        assert_eq!(array, array2);

        let convolved = array.convolve_fold(&[1, 1, 1], 0, |a, b| a + b);

        // 21 = 1 + 2 + 3 + 4 + 5 + 6
        let expect = vec![0, 21, 27, 33, 39, 0];
        assert_eq!(convolved, expect);
    }

    #[test]
    fn test_reset_with_default() {
        let mut array = Array2D::new(vec![vec![1, 2, 3], vec![4, 5, 6]]).unwrap();
        assert_eq!(array.ncols(), 3);
        assert_eq!(array.nrows(), 2);

        array.reset_with_value(2, 3, 0i32);
        assert_eq!(array.values, vec![0, 0, 0, 0, 0, 0]);
        assert_eq!(array.ncols(), 2);
        assert_eq!(array.nrows(), 3);
        assert_eq!(array.get_row(0), Some(vec![0i32, 0i32].as_ref()));
        assert_eq!(array.get_row(1), Some(vec![0i32, 0i32].as_ref()));
        assert_eq!(array.get_row(2), Some(vec![0i32, 0i32].as_ref()));
    }

    #[test]
    fn test_insertion() {
        let mut array = Array2D::new(vec![vec![1, 2, 3], vec![4, 5, 6]]).unwrap();
        array.insert(0, 0, 7);
        assert_eq!(array.values, vec![7, 2, 3, 4, 5, 6]);
        array.insert(1, 2, 8);
        assert_eq!(array.values, vec![7, 2, 3, 4, 5, 8,]);
    }

    #[test]
    fn test_swap_rows() {
        let mut array = Array2D::new(vec![vec![1, 2, 3], vec![4, 5, 6]]).unwrap();
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);
        array.try_swap_rows(0, 1).unwrap();
        assert_eq!(array.values, vec![4, 5, 6, 1, 2, 3]);
        if array.try_swap_rows(1, 2).is_ok() {
            panic!("Should not have succeeded")
        };
    }
}
