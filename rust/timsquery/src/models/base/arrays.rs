use std::fmt::Debug;

use crate::errors::{
    DataProcessingError,
    Result,
};

/// Implements a way to represent an array of
/// dimensions x-y that will be later used to
/// implement an mz-major and a rt-major array
/// representation.
///
/// Simple 2D array
///
/// `values` is a flattened array of values
/// `major_dim` is the number of values in each row
/// `minor_dim` is the number of rows
///
/// Note on memory layout:
///
/// Values that belong to the same row are adjacent
/// in memory.
#[derive(Debug, Clone, PartialEq)]
pub struct Array2D<T: Clone + Copy + std::ops::Mul> {
    pub(super) values: Vec<T>,
    pub(super) major_dim: usize,
    pub(super) minor_dim: usize,
}

impl<
    T: Clone
        + Copy
        + std::ops::Mul<T, Output = T>
        + std::ops::Add<T, Output = T>
        + std::ops::AddAssign<T>,
> Array2D<T>
{
    pub fn new<S: AsRef<[T]>, C: AsRef<[S]>>(values: C) -> Result<Array2D<T>> {
        let nrows = values.as_ref().len();
        if nrows == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData 
            .into());
        }
        let ncols = values.as_ref()[0].as_ref().len();
        if ncols == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData  .into());
        }

        let expected_size = nrows * ncols;
        let values: Vec<T> = values
            .as_ref()
            .iter()
            .flat_map(|x| x.as_ref())
            .cloned()
            .collect();

        if values.len() != expected_size {
            return Err(DataProcessingError::ExpectedVectorSameLength 
            .into());
        }

        Ok(Array2D {
            values,
            major_dim: ncols,
            minor_dim: nrows,
        })
    }

    pub fn new_transposed<S: AsRef<[T]>, C: AsRef<[S]>>(values: C) -> Result<Array2D<T>> {
        let ncols = values.as_ref().len();
        if ncols == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData 
            .into());
        }
        let nrows = values.as_ref()[0].as_ref().len();
        if nrows == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData 
            .into());
        }

        let expected_size = nrows * ncols;
        let mut out_values = vec![None; expected_size];

        for (ci, col) in values.as_ref().iter().enumerate() {
            if col.as_ref().len() != nrows {
                return Err(DataProcessingError::ExpectedVectorSameLength .into());
            }
            for (ri, val) in col.as_ref().iter().enumerate() {
                let idx = ri * ncols + ci; // Changed indexing for row-major order
                out_values[idx] = Some(*val);
            }
        }

        let unwrapped_values: Vec<T> = out_values.into_iter().flatten().collect();

        Ok(Array2D {
            values: unwrapped_values,
            major_dim: ncols,
            minor_dim: nrows,
        })
    }

    pub fn from_flat_vector(values: Vec<T>, nrows: usize, ncols: usize) -> Result<Array2D<T>> {
        if values.len() != nrows * ncols {
            return Err(DataProcessingError::ExpectedVectorSameLength.into());
        }
        Ok(Array2D {
            values,
            major_dim: ncols,
            minor_dim: nrows,
        })
    }

    /// Apply a function to each row of the array
    ///
    /// Example:
    /// ```
    /// use timsseek::models::Array2D;
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
        self.values.chunks(self.major_dim).map(f)
    }

    pub fn iter_mut_rows(&mut self) -> impl Iterator<Item = &mut [T]> {
        self.values.chunks_mut(self.major_dim)
    }

    /// RowConvolve
    /// Apply a a convolution on each row of the array separately.
    /// It is equivalent to applying the passed kernel on each row of the array.
    /// and padding with zeros.
    pub fn row_convolve(&self, kernel: &[T], default_value: T) -> Array2D<T> {
        let mut result = vec![default_value; self.minor_dim * self.major_dim];
        let offset_size = (kernel.len() - 1) / 2;

        for i in 0..self.minor_dim {
            let row_offset = i * (self.major_dim);
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
        let mut result = vec![default_value; self.major_dim];
        let offset_size = (kernel.len() - 1) / 2;

        for i in 0..self.minor_dim {
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

    /// Apply a function to each pair of rows of the array
    ///
    /// For example if I want to calculate the MAE between each row of the array
    /// I can do:
    /// ```
    /// use timsseek::models::Array2D;
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
        for i in 0..self.minor_dim {
            let row = self.get_row(i).expect("Failed to get row, malformed array");
            for j in 0..self.minor_dim {
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
        let idx = row_idx * self.major_dim + col_idx;
        self.values[idx] = value;
    }

    pub fn get_row(&self, index: usize) -> Option<&[T]> {
        let start = index * self.major_dim;
        let end = start + self.major_dim;
        if end > self.values.len() || start >= self.values.len() {
            return None;
        }
        Some(&self.values[start..end])
    }

    pub fn nrows(&self) -> usize {
        self.minor_dim
    }

    pub fn ncols(&self) -> usize {
        self.major_dim
    }

    pub fn transpose(self) -> Array2D<T> {
        // Swap major_dim and minor_dim
        let col_dim = self.major_dim;
        let row_dim = self.minor_dim;
        let vals = self.values;

        let mut result = vec![vals[0]; row_dim * col_dim];

        for i in 0..row_dim {
            for j in 0..col_dim {
                result[j * row_dim + i] = vals[i * col_dim + j];
            }
        }

        Array2D {
            values: result,
            major_dim: row_dim,
            minor_dim: col_dim,
        }
    }

    pub fn reset_with_default(&mut self, ncols: usize, nrows: usize, value: T) {
        self.values.clear();
        self.values.resize(ncols * nrows, value);
        self.major_dim = ncols;
        self.minor_dim = nrows;
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_array2d_new() -> Result<()> {
        // Test creating a 2x3 array
        let values = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let array = Array2D::new(&values)?;

        // Check dimensions
        assert_eq!(array.major_dim, 3); // columns
        assert_eq!(array.minor_dim, 2); // rows

        // Check memory layout - values in same row should be adjacent
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);

        Ok(())
    }

    #[test]
    fn test_array2d_new_transposed() -> Result<()> {
        // Test creating a 3x2 array from columns
        let columns = vec![
            vec![1, 4], // first column
            vec![2, 5], // second column
            vec![3, 6], // third column
        ];
        let array = Array2D::new_transposed(&columns)?;

        // Check dimensions
        assert_eq!(array.major_dim, 3); // columns
        assert_eq!(array.minor_dim, 2); // rows

        // Check memory layout - values should be arranged row-major
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);

        Ok(())
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
    fn test_array2d_large() -> Result<()> {
        // Test with a larger array to verify memory efficiency
        let size = 1000;
        let values: Vec<Vec<f64>> = (0..size)
            .map(|i| (0..size).map(|j| (i * size + j) as f64).collect())
            .collect();

        let array = Array2D::new(&values)?;

        // Verify dimensions
        assert_eq!(array.major_dim, size);
        assert_eq!(array.minor_dim, size);

        // Verify first and last values
        assert_eq!(array.values[0], 0.0);
        assert_eq!(array.values[size * size - 1], (size * size - 1) as f64);

        Ok(())
    }

    #[test]
    fn test_array2d_transpose() -> Result<()> {
        let values = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let array = Array2D::new(&values)?;

        // Check dimensions
        assert_eq!(array.major_dim, 3); // columns
        assert_eq!(array.minor_dim, 2); // rows

        // Check memory layout - values in same row should be adjacent
        assert_eq!(array.values, vec![1, 2, 3, 4, 5, 6]);

        let transposed = array.transpose();

        // Check dimensions
        assert_eq!(transposed.major_dim, 2); // columns
        assert_eq!(transposed.minor_dim, 3); // rows

        // Check memory layout - values should be arranged row-major
        assert_eq!(transposed.values, vec![1, 4, 2, 5, 3, 6]);

        Ok(())
    }

    #[test]
    fn test_convolve() {
        let values = vec![vec![1, 2, 3, 4, 5, 6], vec![4, 5, 6, 7, 8, 9]];
        let array = Array2D::new(&values).unwrap();
        let array2 = Array2D {
            values: vec![1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 9],
            major_dim: 6,
            minor_dim: 2,
        };
        assert_eq!(array, array2);

        let convolved = array.row_convolve(&[1, 1, 1], 0);

        // 6 = 1 + 2 + 3
        // 9 = 2 + 3 + 4
        let expect = Array2D {
            values: vec![0, 6, 9, 12, 15, 0, 0, 15, 18, 21, 24, 0],
            major_dim: 6,
            minor_dim: 2,
        };
        assert_eq!(convolved, expect);
    }

    #[test]
    fn test_convolve_reduce() {
        let values = vec![vec![1, 2, 3, 4, 5, 6], vec![4, 5, 6, 7, 8, 9]];
        let array = Array2D::new(&values).unwrap();
        let array2 = Array2D {
            values: vec![1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 9],
            major_dim: 6,
            minor_dim: 2,
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

        array.reset_with_default(2, 3, 0i32);
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
}
