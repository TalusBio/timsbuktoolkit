use super::{
    Array2D,
    ArrayElement,
};
use crate::KeyLike;
use crate::errors::DataProcessingError;
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
pub struct RTMajorIntensityArray<K: Clone, V: ArrayElement> {
    pub arr: Array2D<V>,
    pub mz_order: Arc<[(K, f64)]>,
    pub rts_ms: Arc<[u32]>,
}

/// Array representation of a series of chromatograms
/// In this representation all elements with the same m/z
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct MzMajorIntensityArray<K: Clone + Eq, V: ArrayElement> {
    pub arr: Array2D<V>,
    pub mz_order: Arc<[(K, f64)]>,
    pub rts_ms: Arc<[u32]>,
}

impl<K: KeyLike, V: ArrayElement + Serialize> Serialize for MzMajorIntensityArray<K, V> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("MzMajorIntensityArray", 3)?;
        // The values should go in different fields ...
        state.serialize_field("arr", &self.arr)?;
        state.serialize_field("mz_order", &self.mz_order)?;
        state.serialize_field("rts_ms", &self.rts_ms)?;
        state.end()
    }
}

#[derive(Debug)]
pub struct MutableChromatogram<'a, V: ArrayElement> {
    slc: &'a mut [V],
    rts: &'a [u32],
}

impl<V: ArrayElement> MutableChromatogram<'_, V> {
    /// Add the passed value at the first position higher
    /// than the internal reference retention time.
    ///
    /// For instance if the reference positions are [10, 20, 30]
    /// any value with rt <10 will be assigned to the first element
    /// any between 10-20 to the second, and any higher to the third.
    pub fn add_at_close_rt<T>(&mut self, rt_ms: u32, intensity: T)
    where
        V: std::ops::AddAssign<T>,
    {
        // For now clamping at the length seems like the right move ...
        // but in the future I might make a more explicit behabior returning an out of bounds
        // error.
        let loc = self
            .rts
            .partition_point(|&rt| rt < rt_ms)
            .min(self.slc.len() - 1);
        self.slc[loc] += intensity;
    }
}

impl<K: KeyLike, V: ArrayElement> MzMajorIntensityArray<K, V> {
    /// Attempt to create a new empty one.
    ///
    /// Errors will be returned if the data is empty.
    pub fn try_new_empty(
        mz_order: Arc<[(K, f64)]>,
        rts_ms: Arc<[u32]>,
    ) -> Result<Self, DataProcessingError> {
        let major_dim = mz_order.len();
        let minor_dim = rts_ms.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        let vals: Vec<_> = vec![V::default(); major_dim * minor_dim];
        let out_arr = Array2D::from_flat_vector(vals, major_dim, minor_dim)
            .expect("CMG array should already be size checked");
        Ok(Self {
            arr: out_arr,
            mz_order,
            rts_ms: rts_ms.clone(),
        })
    }

    /// Get a single 'row', all the chromatogram for a single mz
    pub fn get_row(&self, key: &K) -> Option<&[V]> {
        // Not sure if "row" makes that much sense ... but I feel like get_mz is even less clear ...
        let idx = self.mz_order.iter().position(|(k, _)| k == key)?;
        self.get_row_idx(idx)
    }

    /// Iterate over the values of a single column (all values for a RT)
    pub fn iter_column_idx(&self, index: usize) -> impl '_ + Iterator<Item = (V, &K)> {
        let vals = self.arr.iter_column(index);
        let out = vals.zip(self.mz_order.iter().map(|(k, _mz)| k));
        out
    }

    /// Get a single 'row', all the chromatogram for a single mz
    pub fn get_row_idx(&self, idx: usize) -> Option<&[V]> {
        // Not sure if "row" makes that much sense ... but I feel like get_mz is even less clear ...
        self.arr.get_row(idx)
    }

    /// Get a single 'row', all the chromatogram for a single mz mutably
    fn get_row_mut_idx(&mut self, idx: usize) -> Result<&mut [V], DataProcessingError> {
        // Not sure if "row" makes that much sense ... but I feel like get_mz is even less clear ...
        self.arr.get_row_mut(idx)
    }

    pub fn iter_mut_mzs(
        &mut self,
    ) -> impl Iterator<Item = (&(K, f64), MutableChromatogram<'_, V>)> {
        assert_eq!(self.arr.nrows(), self.mz_order.len());
        self.mz_order
            .iter()
            .zip(self.arr.iter_mut_rows().map(|slc| MutableChromatogram {
                slc,
                rts: &self.rts_ms,
            }))
    }

    /// Create a clone of the array but transposed,
    /// instead of every chromatogram beging adjacent in memory,
    /// every retention time is adjacent in memory.
    pub fn transpose_clone(&self) -> RTMajorIntensityArray<K, V> {
        RTMajorIntensityArray {
            arr: self.arr.transpose_clone(),
            mz_order: self.mz_order.clone(),
            rts_ms: self.rts_ms.clone(),
        }
    }

    /// In-place reorder the array by swapping the chromatograms
    /// to match the passed order.
    pub fn reorder_with(&mut self, order: Arc<[(K, f64)]>) {
        // TODO: make this an error and consume the `self`
        //
        // This is essentially bubble sort ...
        assert_eq!(self.arr.nrows(), self.mz_order.len());
        // Q: Should I check there are no dupes? 2025-Apr-20
        let mut local_order: Vec<_> = self.mz_order.iter().cloned().collect();
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
        self.mz_order = order;
    }

    /// The main purpose of this function is to preserve the allocation
    /// of the array but replace the data contained in it.
    pub fn try_reset_with(
        &mut self,
        array: &MzMajorIntensityArray<K, V>,
    ) -> Result<(), DataProcessingError> {
        // TODO consider if its worth abstracting these 5 lines ...
        let major_dim = array.mz_order.len();
        let minor_dim = array.rts_ms.len();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        self.mz_order = array.mz_order.clone();
        self.rts_ms = array.rts_ms.clone();
        self.arr
            .reset_with_value(minor_dim, major_dim, V::default());
        self.fill_with(array);
        Ok(())
    }

    /// Replace all the values of this array using the values
    /// This is meant to be used after extending the allocation.
    fn fill_with(&mut self, array: &MzMajorIntensityArray<K, V>) {
        let order = array.mz_order.clone();

        for (j, (fh, _fh_mz)) in order.iter().enumerate() {
            let tmp = array.get_row(fh);
            if tmp.is_none() {
                panic!("No row for... {:?}", fh); // , fh, fh_mz);
                // TODO: make sure I dont need this continue ...
                // I think It made sense before but not in the new
                // implementation.
                // continue;
            }
            let inten_slc = tmp.unwrap();
            self.arr
                .try_replace_row_with(j, inten_slc)
                .expect("Sufficient Capacity");
        }
        self.mz_order = order;
        // I am pretty sure this is not needed ...
        // Should I be checking for it??
        self.rts_ms = array.rts_ms.clone();
    }

    pub fn num_ions(&self) -> usize {
        self.mz_order.len()
    }
}

impl<FH: KeyLike, V: ArrayElement> RTMajorIntensityArray<FH, V> {
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
        let vals: Vec<_> = vec![V::default(); minor_dim * major_dim];
        let out_arr = Array2D::from_flat_vector(vals, minor_dim, major_dim)
            .expect("CMG array should already be size checked");
        Ok(Self {
            arr: out_arr,
            mz_order: order.clone(),
            rts_ms: rts_ms.clone(),
        })
    }

    pub fn try_reset_with(
        &mut self,
        array: &MzMajorIntensityArray<FH, V>,
    ) -> Result<(), DataProcessingError> {
        // Q: wtf am I using this for??
        // TODO consider if its worth abstracting these 5 lines ...
        let n_mz = array.mz_order.len();
        let n_rt = array.rts_ms.len();
        if n_mz == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        if n_rt == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        self.mz_order = array.mz_order.clone();
        self.rts_ms = array.rts_ms.clone();
        self.arr.reset_with_value(n_mz, n_rt, V::default());
        self.fill_with(array);
        Ok(())
    }

    fn fill_with(&mut self, array: &MzMajorIntensityArray<FH, V>) {
        let order = array.mz_order.clone();
        for (j, (fh, _mz)) in order.iter().enumerate() {
            let tmp = array.get_row(fh);
            if tmp.is_none() {
                println!("Could not get key");
                continue;
            }
            let inten_slc = tmp.unwrap();
            for (i, &v) in inten_slc.iter().enumerate() {
                // Note that in essence here we are doing a transposition
                // so row indices become col indices
                self.arr.insert(i, j, v);
            }
        }
        self.mz_order = order;
        self.rts_ms = array.rts_ms.clone();
    }
}

#[cfg(test)]
mod tests {

    use crate::models::Array2D;
    use std::sync::Arc;

    use super::{
        MzMajorIntensityArray,
        RTMajorIntensityArray,
    };

    #[test]
    fn test_array2d_int_mz_major() {
        // let array = sample_cmg_array();
        let rt_ms: Arc<[u32]> = [10u32, 20, 30, 40].into();
        let order: Arc<[(String, f64)]> =
            [("FOO".into(), 123.5f64), ("BAR".into(), 456.3f64)].into();
        let mut arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order, rt_ms).unwrap();
        assert_eq!(arr.arr.values, vec![0.0; 8]);
        // I am expecting that adding under bounds gets added to the lowest.
        let to_insert = [(0, 4.0), (21, 2.0)];
        for (x, (_y, mut yc)) in to_insert.iter().zip(arr.iter_mut_mzs()) {
            yc.add_at_close_rt(x.0, x.1);
            yc.add_at_close_rt(x.0, x.1);
            // I am expecting that out of bounds gets added to the last.
            yc.add_at_close_rt(x.0 + 50, x.1);
        }
        assert_eq!(arr.arr.values, vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0]);
    }

    #[test]
    fn test_array2d_int_rt_major() {
        let rts_ms: Arc<[u32]> = [10u32, 20, 30, 40].into();
        let mz_order: Arc<[(String, f64)]> =
            [("FOO".into(), 123.5f64), ("BAR".into(), 456.3f64)].into();
        let vals = vec![
            8.0, 0.0, 0.0, 4.0, // Foo
            0.0, 0.0, 4.0, 2.0, // bar
        ];
        let arr = MzMajorIntensityArray {
            arr: Array2D::from_flat_vector(vals.clone(), mz_order.len(), rts_ms.len()).unwrap(),
            mz_order,
            rts_ms,
        };
        assert_eq!(arr.arr.values, vals);
        let t_arr: RTMajorIntensityArray<String, f32> = arr.transpose_clone();
        let expect = vec![
            8.0, 0.0, // rt < 10
            0.0, 0.0, // rt 10-20
            0.0, 4.0, // rt 20-30
            4.0, 2.0, // rt >30
        ];
        assert_eq!(t_arr.arr.values, expect);
    }

    #[test]
    fn test_array2d_reset_resize_larger() {
        let rt_ms: Arc<[u32]> = [10u32, 20, 30, 40].into();
        let order: Arc<[(String, f64)]> =
            [("FOO".into(), 123.5f64), ("BAR".into(), 456.3f64)].into();
        let mut arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order, rt_ms).unwrap();

        let rt_ms_tgt: Arc<[u32]> = [8678309u32].into();
        let order_tgt: Arc<[(String, f64)]> = [("TOMATO".into(), 123.5f64)].into();
        let mut target_arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order_tgt, rt_ms_tgt).unwrap();

        assert_eq!(arr.arr.values, vec![0.0; 8]);
        assert_eq!(target_arr.arr.values, vec![0.0]);

        // I am expecting that adding under bounds gets added to the lowest.
        let to_insert = [(0, 4.0), (21, 2.0)];
        for (x, (_y, mut yc)) in to_insert.iter().zip(arr.iter_mut_mzs()) {
            yc.add_at_close_rt(x.0, x.1);
            yc.add_at_close_rt(x.0, x.1);
            // I am expecting that out of bounds gets added to the last.
            yc.add_at_close_rt(x.0 + 50, x.1);
        }
        assert_eq!(arr.arr.values, vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0]);

        target_arr.try_reset_with(&arr).unwrap();
        assert_eq!(
            target_arr.arr.values,
            vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0]
        );
    }

    #[test]
    fn test_array2d_reset_resize_smaller() {
        let rt_ms: Arc<[u32]> = [10u32, 20, 30, 40].into();
        let order: Arc<[(String, f64)]> =
            [("FOO".into(), 123.5f64), ("BAR".into(), 456.3f64)].into();
        let mut arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order, rt_ms).unwrap();

        let rt_ms_tgt: Arc<[u32]> = [8678309u32].into();
        let order_tgt: Arc<[(String, f64)]> = [("TOMATO".into(), 123.5f64)].into();
        let target_arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order_tgt, rt_ms_tgt).unwrap();

        assert_eq!(arr.arr.values, vec![0.0; 8]);
        assert_eq!(target_arr.arr.values, vec![0.0]);

        arr.try_reset_with(&target_arr).unwrap();
        assert_eq!(arr.arr.values, vec![0.0]);
    }

    #[test]
    fn test_array2d_reset_resize_larger_to_rt() {
        let rt_ms: Arc<[u32]> = [10u32, 20, 30, 40].into();
        let order: Arc<[(String, f64)]> =
            [("FOO".into(), 123.5f64), ("BAR".into(), 456.3f64)].into();
        let mut arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order, rt_ms).unwrap();

        let rt_ms_tgt: Arc<[u32]> = [8678309u32].into();
        let order_tgt: Arc<[(String, f64)]> = [("TOMATO".into(), 123.5f64)].into();
        let mut target_arr: RTMajorIntensityArray<String, f32> =
            RTMajorIntensityArray::try_new_empty(order_tgt, rt_ms_tgt).unwrap();

        assert_eq!(arr.arr.values, vec![0.0; 8]);
        assert_eq!(target_arr.arr.values, vec![0.0]);

        // I am expecting that adding under bounds gets added to the lowest.
        let to_insert = [(0, 4.0), (21, 2.0)];
        for (x, (_y, mut yc)) in to_insert.iter().zip(arr.iter_mut_mzs()) {
            yc.add_at_close_rt(x.0, x.1);
            yc.add_at_close_rt(x.0, x.1);
            // I am expecting that out of bounds gets added to the last.
            yc.add_at_close_rt(x.0 + 50, x.1);
        }

        target_arr.try_reset_with(&arr).unwrap();
        // Note the change in position
        assert_eq!(arr.arr.values, vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0]);
        assert_eq!(
            target_arr.arr.values,
            vec![8.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 2.0]
        );
    }
}
