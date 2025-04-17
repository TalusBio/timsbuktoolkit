use super::Array2D;
use crate::errors::DataProcessingError;
use std::sync::Arc;
use crate::KeyLike;

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
    pub fn new_empty(order: Arc<[K]>, rts_ms: Arc<[u32]>) -> Result<Self, DataProcessingError> {
        let major_dim = order.len();
        let minor_dim = rts_ms.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData );
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

    pub fn iter_mut_mzs(&mut self) -> impl Iterator<Item = (&K, &mut[f32])> {
        assert_eq!(self.arr.nrows(), self.order.len());
        self.order.iter().zip(self.arr.iter_mut_rows())
    }

    // // The main purpose of this function is to preserve the allocation
    // // of the array but replace the data contained in it.
    // pub fn reset_with(
    //     &mut self,
    //     array: &PartitionedCMGArrayStats<K>,
    //     order: Arc<[K]>,
    // ) -> Result<(), DataProcessingError> {
    //     // TODO consider if its worth abstracting these 5 lines ...
    //     let major_dim = order.len();
    //     let minor_dim = array.retention_time_miliseconds.len();
    //     if minor_dim == 0 {
    //         return Err(DataProcessingError::ExpectedNonEmptyData {
    //             context: Some("Cannot create array with zero rows".to_string()),
    //         });
    //     }
    //     self.arr.reset_with_default(minor_dim, major_dim, 0.0);
    //     self.fill_with(array, order);
    //     Ok(())
    // }

    // fn fill_with(&mut self, array: &PartitionedCMGArrayStats<K>, order: Arc<[K]>) {
    //     for (j, fh) in order.iter().enumerate() {
    //         let tmp = array.intensities.get(fh);
    //         if tmp.is_none() {
    //             continue;
    //         }
    //         let inten_vec = tmp.unwrap();
    //         for (i, inten) in inten_vec.iter().enumerate() {
    //             self.arr.insert(j, i, *inten as f32);
    //         }
    //     }
    //     self.order = order;
    //     self.rts_ms = array.retention_time_miliseconds.clone();
    // }
}

impl<FH: KeyLike> RTMajorIntensityArray<FH> {
    pub fn new_empty(order: Arc<[FH]>, rts_ms: Arc<[u32]>) -> Result<Self, DataProcessingError> {
        let minor_dim = rts_ms.len();
        let major_dim = order.len();
        if major_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData );
        }
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData );
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

    // pub fn reset_with(
    //     &mut self,
    //     array: &PartitionedCMGArrayStats<FH>,
    //     order: Arc<[FH]>,
    // ) -> Result<(), DataProcessingError> {
    //     // TODO consider if its worth abstracting these 5 lines ...
    //     let minor_dim = order.len();
    //     let major_dim = array.retention_time_miliseconds.len();
    //     if major_dim == 0 {
    //         return Err(DataProcessingError::ExpectedNonEmptyData {
    //             context: Some("Cannot create array with zero columns".to_string()),
    //         });
    //     }
    //     if minor_dim == 0 {
    //         return Err(DataProcessingError::ExpectedNonEmptyData {
    //             context: Some("Cannot create array with zero rows".to_string()),
    //         });
    //     }
    //     self.arr.reset_with_default(minor_dim, major_dim, 0.0);
    //     self.fill_with(array, order);
    //     Ok(())
    // }

}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_array2d_from_cmg_int_mz_major() {
        todo!();
        // let array = sample_cmg_array();
        // let arr = MzMajorIntensityArray::new(&array, [0, 1].into()).unwrap();
        // assert_eq!(arr.arr.values, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        // let first_row = arr.arr.get_row(0).unwrap();
        // assert_eq!(first_row, &[1.0, 2.0, 3.0]);
        // assert_eq!(arr.arr.major_dim, 3);
        // assert_eq!(arr.arr.minor_dim, 2);
    }

    #[test]
    fn test_array2d_from_cmg_int_rt_major() {
        todo!();
        // let array = sample_cmg_array();
        // let arr = RTMajorIntensityArray::new(&array, [0, 1].into()).unwrap();
        // assert_eq!(arr.arr.values, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        // let first_row = arr.arr.get_row(0).unwrap();
        // assert_eq!(first_row, &[1.0, 4.0]);
        // assert_eq!(arr.arr.major_dim, 2);
        // assert_eq!(arr.arr.minor_dim, 3);
    }
}
