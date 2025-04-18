use super::Array2D;
use crate::KeyLike;
use crate::errors::DataProcessingError;
use crate::models::aggregators::EGCAggregator;
use serde::ser::{
    Serialize,
    SerializeStruct,
    Serializer,
};
use std::sync::Arc;

/// Array representation of a series of chromatograms
/// In this representation all elements with the same retention time
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct RTMajorIntensityArray<K: Clone> {
    pub arr: Array2D<f32>,
    pub order: Arc<[(K, f64)]>,
    pub rts_ms: Arc<[u32]>,
}

/// Array representation of a series of chromatograms
/// In this representation all elements with the same m/z
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct MzMajorIntensityArray<K: Clone + Eq> {
    pub arr: Array2D<f32>,
    pub order_mz: Arc<[(K, f64)]>,
    pub rts_ms: Arc<[u32]>,
}

impl<K: KeyLike> Serialize for MzMajorIntensityArray<K> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("MzMajorIntensityArray", 3)?;
        println!("FINISH IMPLEMENTING THISSS!!!");
        // The values should go in different fields ...
        state.serialize_field("arr", &self.arr.values)?;
        state.serialize_field("order_mz", &self.order_mz)?;
        state.serialize_field("rts_ms", &self.rts_ms)?;
        state.end()
    }
}

#[derive(Debug)]
pub struct MutableChromatogram<'a> {
    slc: &'a mut [f32],
    rts: &'a [u32],
}

impl<'a> MutableChromatogram<'a> {
    pub fn add_at_close_rt(&mut self, rt_ms: u32, intensity: f32) {
        let loc = self.rts.partition_point(|&rt| rt < rt_ms);
        self.slc[loc] += intensity;
    }
}

impl<K: KeyLike> MzMajorIntensityArray<K> {
    pub fn try_new_empty(
        order_mz: Arc<[(K, f64)]>,
        rts_ms: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        let major_dim = order_mz.len();
        let minor_dim = rts_ms.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        let vals: Vec<f32> = vec![0.0; major_dim * minor_dim];
        let out_arr = Array2D::from_flat_vector(vals, major_dim, minor_dim)
            .expect("CMG array should already be size checked");
        Ok(Self {
            arr: out_arr,
            order_mz,
            rts_ms: rts_ms.clone(),
        })
    }

    // Not sure if "row" makes that much sense ... but I feel like get_mz is even less clear ...
    pub fn get_row(&self, key: &K) -> Option<&[f32]> {
        let idx = self.order_mz.iter().position(|(k, _)| k == key)?;
        self.arr.get_row(idx)
    }

    pub fn iter_mut_mzs(&mut self) -> impl Iterator<Item = (&(K, f64), MutableChromatogram<'_>)> {
        assert_eq!(self.arr.nrows(), self.order_mz.len());
        self.order_mz
            .iter()
            .zip(self.arr.iter_mut_rows().map(|slc| MutableChromatogram {
                slc,
                rts: &self.rts_ms,
            }))
    }

    pub fn transpose_clone(&self) -> RTMajorIntensityArray<K> {
        RTMajorIntensityArray {
            arr: self.arr.transpose_clone(),
            order: self.order_mz.clone(),
            rts_ms: self.rts_ms.clone(),
        }
    }

    pub fn reorder_with(&mut self, order: Arc<[(K, f64)]>) {
        // This is essentially bubble sort ...
        assert_eq!(self.arr.nrows(), self.order_mz.len());
        // Should I check there are no dupes?
        //
        let mut local_order: Vec<_> = self.order_mz.iter().cloned().collect();
        for (i, (k, _mz)) in order.iter().enumerate() {
            let j = local_order.iter().position(|(k2, _)| k2 == k).unwrap();
            if i != j {
                local_order.swap(i, j);
                self.arr.try_swap_rows(i, j).expect("Array in bounds");
            }
        }

        for (l, r) in local_order.iter().zip(order.iter()) {
            assert_eq!(l, r);
        }
        self.order_mz = order;
    }

    // The main purpose of this function is to preserve the allocation
    // of the array but replace the data contained in it.
    pub fn try_reset_with(
        &mut self,
        array: &MzMajorIntensityArray<K>,
    ) -> Result<(), DataProcessingError> {
        // TODO consider if its worth abstracting these 5 lines ...
        let major_dim = array.order_mz.len();
        let minor_dim = array.rts_ms.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        self.arr.reset_with_default(minor_dim, major_dim, 0.0);
        self.fill_with(array);
        Ok(())
    }

    fn fill_with(&mut self, array: &MzMajorIntensityArray<K>) {
        let order = array.order_mz.clone();

        for (j, (fh, fh_mz)) in order.iter().enumerate() {
            let tmp = array.get_row(fh);
            if tmp.is_none() {
                panic!("No row for..."); // , fh, fh_mz);
                continue;
            }
            let inten_slc = tmp.unwrap();
            self.arr
                .try_replace_row_with(j, inten_slc)
                .expect("Sufficient Capacity");
        }
        self.order_mz = order;
        // I am pretty sure this is not needed ...
        // Should I be checking for it??
        self.rts_ms = array.rts_ms.clone();
    }
}

impl<FH: KeyLike> RTMajorIntensityArray<FH> {
    pub fn try_new_empty(
        order: Arc<[(FH, f64)]>,
        rts_ms: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        let minor_dim = rts_ms.len();
        let major_dim = order.len();
        if major_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
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

    pub fn try_reset_with(
        &mut self,
        array: &MzMajorIntensityArray<FH>,
    ) -> Result<(), DataProcessingError> {
        // TODO consider if its worth abstracting these 5 lines ...
        let minor_dim = array.order_mz.len();
        let major_dim = array.rts_ms.len();
        if major_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        self.arr.reset_with_default(minor_dim, major_dim, 0.0);
        self.fill_with(array);
        Ok(())
    }

    fn fill_with(&mut self, array: &MzMajorIntensityArray<FH>) {
        let order = array.order_mz.clone();
        for (j, (fh, _mz)) in order.iter().enumerate() {
            let tmp = array.get_row(fh);
            if tmp.is_none() {
                continue;
            }
            let inten_slc = tmp.unwrap();
            for (i, _) in inten_slc.iter().enumerate() {
                self.arr.insert(i, j, inten_slc[i]);
            }
        }
        self.order = order;
        self.rts_ms = array.rts_ms.clone();
    }
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
