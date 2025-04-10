use super::Array2D;
use crate::errors::DataProcessingError;
use std::sync::Arc;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::arrays::PartitionedCMGArrayStats;
use timsquery::traits::key_like::KeyLike;

/// Array representation of a series of chromatograms
/// In this representation all elements with the same retention time
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct RTMajorIntensityArray<K: Clone> {
    pub arr: Array2D<f32>,
    pub order: Arc<[K]>,
    pub rts_ms: Arc<[u32]>,
}

/// Array representation of a series of chromatograms
/// In this representation all elements with the same m/z
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct MzMajorIntensityArray<K: Clone> {
    pub arr: Array2D<f32>,
    pub order: Arc<[K]>,
    pub rts_ms: Arc<[u32]>,
}

impl<K: KeyLike> MzMajorIntensityArray<K> {
    pub fn new(
        array: &PartitionedCMGArrayStats<K>,
        order: Arc<[K]>,
    ) -> std::result::Result<Self, DataProcessingError> {
        // TODO: Do I need to check if the order has no duplicates?
        let major_dim = order.len();
        let minor_dim = array.retention_time_miliseconds.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            });
        }
        let vals: Vec<f32> = vec![0.0; major_dim * minor_dim];
        let out_arr = Array2D::from_flat_vector(vals, major_dim, minor_dim)
            .expect("CMG array should already be size checked");
        let mut out = Self {
            arr: out_arr,
            order: order.clone(),
            rts_ms: array.retention_time_miliseconds.clone(),
        };
        out.fill_with(array, order);

        Ok(out)
    }

    pub fn new_empty(order: Arc<[K]>, rts_ms: Arc<[u32]>) -> Result<Self, DataProcessingError> {
        let major_dim = order.len();
        let minor_dim = rts_ms.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            });
        }
        let vals: Vec<f32> = vec![0.0; major_dim * minor_dim];
        let out_arr = Array2D::from_flat_vector(vals, major_dim, minor_dim)
            .expect("CMG array should already be size checked");
        Ok(Self {
            arr: out_arr,
            order: order.clone(),
            rts_ms: rts_ms.clone(),
        })
    }

    // The main purpose of this function is to preserve the allocation
    // of the array but replace the data contained in it.
    pub fn reset_with(
        &mut self,
        array: &PartitionedCMGArrayStats<K>,
        order: Arc<[K]>,
    ) -> Result<(), DataProcessingError> {
        // TODO consider if its worth abstracting these 5 lines ...
        let major_dim = order.len();
        let minor_dim = array.retention_time_miliseconds.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData {
                context: Some("Cannot create array with zero rows".to_string()),
            });
        }
        self.arr.reset_with_default(minor_dim, major_dim, 0.0);
        self.fill_with(array, order);
        Ok(())
    }

    fn fill_with(&mut self, array: &PartitionedCMGArrayStats<K>, order: Arc<[K]>) {
        for (j, fh) in order.iter().enumerate() {
            let tmp = array.intensities.get(fh);
            if tmp.is_none() {
                continue;
            }
            let inten_vec = tmp.unwrap();
            for (i, inten) in inten_vec.iter().enumerate() {
                self.arr.insert(j, i, *inten as f32);
            }
        }
        self.order = order;
        self.rts_ms = array.retention_time_miliseconds.clone();
    }
}

// TODO: Actually make the PartitionedCMGArrayStats store its data as the
//       mz-major array. and implement the serialization so it preserves the
//       output (I like the current serialization but its not great as a mem
//       layout)
impl<FH: KeyLike> RTMajorIntensityArray<FH> {
    pub fn new(
        array: &PartitionedCMGArrayStats<FH>,
        order: Arc<[FH]>,
    ) -> core::result::Result<Self, DataProcessingError> {
        let major_dim = order.len();
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

        let vals: Vec<f32> = vec![0.0; minor_dim * major_dim];
        let out_arr = Array2D::from_flat_vector(vals, minor_dim, major_dim)
            .expect("CMG array should already be size checked");
        let mut out = Self {
            arr: out_arr,
            order: order.clone(), // TODO fix this un-necessary clone later ...
            rts_ms: array.retention_time_miliseconds.clone(),
        };
        out.fill_with(array, order);

        Ok(out)
    }

    pub fn new_empty(order: Arc<[FH]>, rts_ms: Arc<[u32]>) -> Result<Self, DataProcessingError> {
        let minor_dim = rts_ms.len();
        let major_dim = order.len();
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
        let vals: Vec<f32> = vec![0.0; minor_dim * major_dim];
        let out_arr = Array2D::from_flat_vector(vals, minor_dim, major_dim)
            .expect("CMG array should already be size checked");
        Ok(Self {
            arr: out_arr,
            order: order.clone(),
            rts_ms: rts_ms.clone(),
        })
    }

    pub fn reset_with(
        &mut self,
        array: &PartitionedCMGArrayStats<FH>,
        order: Arc<[FH]>,
    ) -> Result<(), DataProcessingError> {
        // TODO consider if its worth abstracting these 5 lines ...
        let minor_dim = order.len();
        let major_dim = array.retention_time_miliseconds.len();
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
        self.arr.reset_with_default(minor_dim, major_dim, 0.0);
        self.fill_with(array, order);
        Ok(())
    }

    fn fill_with(&mut self, array: &PartitionedCMGArrayStats<FH>, order: Arc<[FH]>) {
        for (j, fh) in order.iter().enumerate() {
            let tmp = array.intensities.get(fh);
            if tmp.is_none() {
                // This in theory can happen if there are no peaks
                // that match that specific transition.
                continue;
            }
            for (i, inten) in tmp.unwrap().iter().enumerate() {
                self.arr.insert(i, j, *inten as f32);
            }
        }
        self.order = order;
        self.rts_ms = array.retention_time_miliseconds.clone();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn sample_cmg_array() -> PartitionedCMGArrayStats<u8> {
        let retention_time_miliseconds: Arc<[u32]> = vec![0, 1, 2].into();
        let mut intensities: HashMap<u8, Vec<u64>> = HashMap::new();
        intensities.insert(0, vec![1, 2, 3]);
        intensities.insert(1, vec![4, 5, 6]);

        let mut ims = HashMap::new();
        let mut mz = HashMap::new();
        ims.insert(0, vec![1., 2., 3.]);
        ims.insert(1, vec![4., 5., 6.]);
        mz.insert(0, vec![1., 2., 3.]);
        mz.insert(1, vec![4., 5., 6.]);

        let weighted_ims_mean = vec![2.0, 4.0, 5.0];

        PartitionedCMGArrayStats {
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
        let arr = MzMajorIntensityArray::new(&array, [0, 1].into()).unwrap();
        assert_eq!(arr.arr.values, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let first_row = arr.arr.get_row(0).unwrap();
        assert_eq!(first_row, &[1.0, 2.0, 3.0]);
        assert_eq!(arr.arr.major_dim, 3);
        assert_eq!(arr.arr.minor_dim, 2);
    }

    #[test]
    fn test_array2d_from_cmg_int_rt_major() {
        let array = sample_cmg_array();
        let arr = RTMajorIntensityArray::new(&array, [0, 1].into()).unwrap();
        assert_eq!(arr.arr.values, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        let first_row = arr.arr.get_row(0).unwrap();
        assert_eq!(first_row, &[1.0, 4.0]);
        assert_eq!(arr.arr.major_dim, 2);
        assert_eq!(arr.arr.minor_dim, 3);
    }
}
