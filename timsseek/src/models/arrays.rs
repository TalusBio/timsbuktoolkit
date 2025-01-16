use serde::Serialize;
use std::fmt::Debug;
use std::hash::Hash;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::arrays::PartitionedCMGArrayStats;

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
#[derive(Debug, Clone)]
pub struct Array2D<T: Clone + Copy> {
    values: Vec<T>,
    major_dim: usize,
    minor_dim: usize,
}

impl<T: Clone + Copy> Array2D<T> {
    pub fn new<S: AsRef<[T]>, C: AsRef<[S]>>(values: C) -> Result<Array2D<T>> {
        let nrows = values.as_ref().len();
        if nrows == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            }
            .into());
        }
        let ncols = values.as_ref()[0].as_ref().len();
        if ncols == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero columns".to_string()),
            }
            .into());
        }

        let expected_size = nrows * ncols;
        let values: Vec<T> = values
            .as_ref()
            .iter()
            .flat_map(|x| x.as_ref())
            .cloned()
            .collect();

        if values.len() != expected_size {
            return Err(DataProcessingError::ExpectedSlicesSameLength {
                expected: expected_size,
                other: values.len(),
                context: "Expected array with nrows * ncols values, got different size".to_string(),
            }
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
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero columns".to_string()),
            }
            .into());
        }
        let nrows = values.as_ref()[0].as_ref().len();
        if nrows == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            }
            .into());
        }

        let expected_size = nrows * ncols;
        let mut out_values = vec![None; expected_size];

        for (ci, col) in values.as_ref().iter().enumerate() {
            if col.as_ref().len() != nrows {
                return Err(DataProcessingError::ExpectedSlicesSameLength {
                    expected: nrows,
                    other: col.as_ref().len(),
                    context: "All columns must have the same length".to_string(),
                }
                .into());
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
            return Err(DataProcessingError::ExpectedSlicesSameLength {
                expected: nrows * ncols,
                other: values.len(),
                context: "Expected array with nrows * ncols values, got different size".to_string(),
            }
            .into());
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
    /// let result: Vec<u32> = array.row_apply(|x| x.iter().sum());
    /// assert_eq!(result, vec![6, 15]);
    ///
    /// let array = Array2D::new_transposed(vec![vec![1, 2, 3], vec![4, 5, 6]]).unwrap();
    /// let result: Vec<u32> = array.row_apply(|x| x.iter().sum());
    /// assert_eq!(result, vec![5, 7, 9]);
    /// ```
    pub fn row_apply<W, F: FnMut(&[T]) -> W>(&self, f: F) -> Vec<W> {
        self.values
            .chunks(self.major_dim)
            .map(f)
            .collect::<Vec<W>>()
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
}

/// Array representation of a series of chromatograms
/// In this representation all elements with the same retention time
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct RTMajorIntensityArray {
    pub arr: Array2D<f32>,
}

/// Array representation of a series of chromatograms
/// In this representation all elements with the same m/z
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct MzMajorIntensityArray {
    pub arr: Array2D<f32>,
}

impl MzMajorIntensityArray {
    pub fn new<FH: Clone + Eq + Serialize + Hash + Send + Sync + Debug>(
        array: &PartitionedCMGArrayStats<FH>,
        order: Option<&[FH]>,
    ) -> std::result::Result<Self, DataProcessingError> {
        // TODO: Do I need to check if the order has no duplicates?
        let major_dim = match order {
            Some(order) => order.len(),
            None => array.intensities.len(),
        };
        let minor_dim = array.retention_time_miliseconds.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            });
        }
        let mut vals: Vec<f32> = vec![0.0; major_dim * minor_dim];
        match order {
            Some(order) => {
                for (j, fh) in order.iter().enumerate() {
                    let tmp = array.intensities.get(fh);
                    if tmp.is_none() {
                        continue;
                    }
                    let inten_vec = tmp.unwrap();
                    for (i, inten) in inten_vec.iter().enumerate() {
                        let idx = j * minor_dim + i;
                        vals[idx] = *inten as f32;
                    }
                }
            }
            None => {
                for (j, inten_vec) in array.intensities.values().enumerate() {
                    for (i, inten) in inten_vec.iter().enumerate() {
                        let idx = j * minor_dim + i;
                        vals[idx] = *inten as f32;
                    }
                }
            }
        }

        let out_arr = Array2D::from_flat_vector(vals, major_dim, minor_dim)
            .expect("CMG array should already be size checked");

        // This assertion just makes sure I did not mix up rows and columns
        // in this implementation.
        assert!(array.intensities.len() <= out_arr.nrows());

        Ok(Self { arr: out_arr })
    }
}

// TODO: Actually make the PartitionedCMGArrayStats store its data as the
//       mz-major array. and implement the serialization so it preserves the
//       output (I like the current serialization but its not great as a mem
//       layout)
impl RTMajorIntensityArray {
    pub fn new<FH: Clone + Eq + Serialize + Hash + Send + Sync>(
        array: &PartitionedCMGArrayStats<FH>,
        order: Option<&[FH]>,
    ) -> core::result::Result<Self, DataProcessingError> {
        let major_dim = match order {
            Some(order) => order.len(),
            None => array.intensities.len(),
        };
        let minor_dim = array.retention_time_miliseconds.len();
        if major_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero columns".to_string()),
            });
        }
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            });
        }

        let mut vals: Vec<f32> = vec![0.0; minor_dim * major_dim];
        match order {
            Some(order) => {
                for (j, fh) in order.iter().enumerate() {
                    let tmp = array.intensities.get(fh);
                    if tmp.is_none() {
                        // This in theory can happen if there are no peaks
                        // that match that specific transition.
                        continue;
                    }
                    for (i, inten) in tmp.unwrap().iter().enumerate() {
                        let idx = i * major_dim + j;
                        vals[idx] = *inten as f32;
                    }
                }
            }
            None => {
                for (j, inten_vec) in array.intensities.values().enumerate() {
                    for (i, inten) in inten_vec.iter().enumerate() {
                        let idx = i * major_dim + j;
                        vals[idx] = *inten as f32;
                    }
                }
            }
        }

        let out_arr = Array2D::from_flat_vector(vals, minor_dim, major_dim)
            .expect("CMG array should already be size checked");

        // This assertion just makes sure I did not mix up rows and columns
        // in this implementation.
        assert_eq!(out_arr.values.len(), minor_dim * major_dim);

        Ok(Self { arr: out_arr })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

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

    fn sample_cmg_array() -> PartitionedCMGArrayStats<u8> {
        let retention_time_miliseconds = vec![0, 1, 2];
        let mut intensities: HashMap<u8, Vec<u64>> = HashMap::new();
        intensities.insert(0, vec![1, 2, 3]);
        intensities.insert(1, vec![4, 5, 6]);

        let summed_intensity_vec: Vec<u64> = vec![5, 7, 9];

        let mut ims = HashMap::new();
        let mut mz = HashMap::new();
        ims.insert(0, vec![1., 2., 3.]);
        ims.insert(1, vec![4., 5., 6.]);
        mz.insert(0, vec![1., 2., 3.]);
        mz.insert(1, vec![4., 5., 6.]);

        let weighted_ims_mean = vec![2.0, 4.0, 5.0];

        PartitionedCMGArrayStats {
            summed_intensity: summed_intensity_vec,
            intensities,
            weighted_ims_mean,
            retention_time_miliseconds,
            ims_means: ims,
            mz_means: mz,
        }
    }

    #[test]
    fn test_array2d_from_cmg_int_mz_major() {
        let array = sample_cmg_array();
        let arr = MzMajorIntensityArray::new(&array, Some(&[0, 1])).unwrap();
        assert_eq!(arr.arr.values, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let first_row = arr.arr.get_row(0).unwrap();
        assert_eq!(first_row, &[1.0, 2.0, 3.0]);
        assert_eq!(arr.arr.major_dim, 3);
        assert_eq!(arr.arr.minor_dim, 2);
    }

    #[test]
    fn test_array2d_from_cmg_int_rt_major() {
        let array = sample_cmg_array();
        let arr = RTMajorIntensityArray::new(&array, Some(&[0, 1])).unwrap();
        assert_eq!(arr.arr.values, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        let first_row = arr.arr.get_row(0).unwrap();
        assert_eq!(first_row, &[1.0, 4.0]);
        assert_eq!(arr.arr.major_dim, 2);
        assert_eq!(arr.arr.minor_dim, 3);
    }
}
