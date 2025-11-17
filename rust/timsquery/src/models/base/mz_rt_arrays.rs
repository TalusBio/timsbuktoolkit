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
    pub cycle_offset: usize,
}

/// Array representation of a series of chromatograms
/// In this representation all elements with the same m/z
/// are adjacent in memory.
#[derive(Debug, Clone)]
pub struct MzMajorIntensityArray<K: Clone + Eq, V: ArrayElement> {
    pub arr: Array2D<V>,
    pub mz_order: Arc<[(K, f64)]>,
    pub cycle_offset: usize,
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
        state.end()
    }
}

/// A mutable view into a single chromatogram within a MzMajorIntensityArray
///
/// The meaning of the cycle_offset is that the indices passed to add_at_index
/// are adjusted by subtracting the cycle_offset before accessing the underlying slice.
#[derive(Debug)]
pub struct Chromatogram<'a, V: ArrayElement> {
    /// The slice of intensities for this chromatogram
    slc: &'a mut [V],

    /// The offset to apply to indices when adding intensities.
    /// In essence it means "how much do I need to add to the index in this slice to match
    /// the global cyle (retention time) index".
    cycle_offset: usize,
}

impl<V: ArrayElement> Chromatogram<'_, V> {
    fn new<'a>(slc: &'a mut [V], cycle_offset: usize) -> Chromatogram<'a, V> {
        Chromatogram { slc, cycle_offset }
    }

    pub fn as_slice(&self) -> &[V] {
        self.slc
    }

    pub fn try_get_slice(&self, start: usize, end: usize) -> Option<&[V]> {
        // add index offsetting
        let start = start.saturating_sub(self.cycle_offset);
        let end = end.saturating_sub(self.cycle_offset);
        if start >= self.slc.len() || end > self.slc.len() || start >= end {
            return None;
        }
        Some(&self.slc[start..end])
    }

    pub fn len(&self) -> usize {
        self.slc.len()
    }

    /// Add intensity at the given index, adjusting for cycle offset.
    /// If the index is out of bounds, add to the last element.
    pub fn add_at_index<T>(&mut self, idx: u32, intensity: T)
    where
        V: std::ops::AddAssign<T>,
    {
        let idx = idx.saturating_sub(self.cycle_offset as u32);
        if idx as usize >= self.slc.len() {
            self.slc[self.slc.len() - 1] += intensity;
            return;
        }
        self.slc[idx as usize] += intensity;
    }
}

impl<K: KeyLike, V: ArrayElement> MzMajorIntensityArray<K, V> {
    /// Attempt to create a new empty one.
    ///
    /// Errors will be returned if the data is empty.
    pub fn try_new_empty(
        mz_order: Arc<[(K, f64)]>,
        rt_cycles: usize,
        cycle_offset: usize,
    ) -> Result<Self, DataProcessingError> {
        let major_dim = mz_order.len();
        let minor_dim = rt_cycles;
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        let vals: Vec<_> = vec![V::default(); major_dim * minor_dim];
        let out_arr = Array2D::from_flat_vector(vals, major_dim, minor_dim)
            .expect("CMG array should already be size checked");
        Ok(Self {
            arr: out_arr,
            mz_order,
            cycle_offset,
        })
    }

    /// Get a single 'row', all the chromatogram for a single mz
    pub fn get_row(&self, key: &K) -> Option<&[V]> {
        // Not sure if "row" makes that much sense ... but I feel like get_mz is even less clear ...
        let idx = self.mz_order.iter().position(|(k, _)| k == key)?;
        self.get_row_idx(idx)
    }

    /// Iterate over the values of a single column (all values for a RT)
    pub fn iter_column_idx(&self, index: usize) -> impl '_ + Iterator<Item = (&K, V)> {
        let vals = self.arr.iter_column(index);
        let out = self.mz_order.iter().map(|(k, _mz)| k).zip(vals);
        out
    }

    /// Get a single 'row', all the chromatogram for a single mz
    pub fn get_row_idx(&self, idx: usize) -> Option<&[V]> {
        // Not sure if "row" makes that much sense ... but I feel like get_mz is even less clear ...
        self.arr.get_row(idx)
    }

    pub fn iter_mut_mzs(&mut self) -> impl Iterator<Item = (&(K, f64), Chromatogram<'_, V>)> {
        assert_eq!(self.arr.nrows(), self.mz_order.len());
        self.mz_order.iter().zip(
            self.arr
                .iter_mut_rows()
                .map(|slc| Chromatogram::new(slc, self.cycle_offset)),
        )
    }

    pub fn iter_mzs(&self) -> impl Iterator<Item = (&(K, f64), &[V])> {
        assert_eq!(self.arr.nrows(), self.mz_order.len());
        self.mz_order.iter().zip(self.arr.iter_rows())
    }

    /// Create a clone of the array but transposed,
    /// instead of every chromatogram beging adjacent in memory,
    /// every retention time is adjacent in memory.
    pub fn transpose_clone(&self) -> RTMajorIntensityArray<K, V> {
        RTMajorIntensityArray {
            arr: self.arr.transpose_clone(),
            mz_order: self.mz_order.clone(),
            cycle_offset: self.cycle_offset,
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
        let minor_dim = array.arr.ncols();
        if minor_dim == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        self.mz_order = array.mz_order.clone();
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
        self.cycle_offset = array.cycle_offset;
    }

    pub fn num_ions(&self) -> usize {
        self.mz_order.len()
    }
}

impl<FH: KeyLike, V: ArrayElement> RTMajorIntensityArray<FH, V> {
    pub fn try_new_empty(
        order: Arc<[(FH, f64)]>,
        num_cycles: usize,
        cycle_offset: usize,
    ) -> Result<Self, DataProcessingError> {
        let minor_dim = num_cycles;
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
            cycle_offset,
        })
    }

    pub fn try_reset_with(
        &mut self,
        array: &MzMajorIntensityArray<FH, V>,
    ) -> Result<(), DataProcessingError> {
        // Q: wtf am I using this for??
        // TODO consider if its worth abstracting these 5 lines ...
        let n_mz = array.mz_order.len();
        let n_rt = array.arr.ncols();
        if n_mz == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        if n_rt == 0 {
            return Err(DataProcessingError::ExpectedNonEmptyData);
        }
        self.mz_order = array.mz_order.clone();
        self.cycle_offset = array.cycle_offset;
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
        self.cycle_offset = array.cycle_offset;
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
            MzMajorIntensityArray::try_new_empty(order, rt_ms.len(), 0).unwrap();
        assert_eq!(arr.arr.values, vec![0.0; 8]);
        // I am expecting that adding under bounds gets added to the lowest.
        let to_insert = [(0, 4.0), (21, 2.0)];
        for (x, (_y, mut yc)) in to_insert.iter().zip(arr.iter_mut_mzs()) {
            // On each iteration we get [0 ,0, 0, 0]
            yc.add_at_index(x.0, x.1);
            yc.add_at_index(x.0, x.1);
            // here should be [8.0 ,0, 0, 0  ] on iteration 1
            // here should be [0   ,0, 0, 4.0] on iteration 2
            // I am expecting that out of bounds gets added to the last.
            yc.add_at_index(x.0 + 50, x.1);
            // here should be [8.0 ,0, 0, 4.0] on iteration 1
            // here should be [0   ,0, 0, 6.0] on iteration 2
        }
        assert_eq!(arr.arr.values, vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 6.0]);
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
            cycle_offset: 0,
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
            MzMajorIntensityArray::try_new_empty(order, rt_ms.len(), 0).unwrap();

        let rt_ms_tgt: Arc<[u32]> = [8678309u32].into();
        let order_tgt: Arc<[(String, f64)]> = [("TOMATO".into(), 123.5f64)].into();
        let mut target_arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order_tgt, rt_ms_tgt.len(), 0).unwrap();

        assert_eq!(arr.arr.values, vec![0.0; 8]);
        assert_eq!(target_arr.arr.values, vec![0.0]);

        // I am expecting that adding under bounds gets added to the lowest.
        // let to_insert = [(0, 4.0), (21, 2.0)]; // Old versions added based on RT, new ones based
        // on index
        let to_insert = [(0, 4.0), (1, 2.0)];
        for (x, (_y, mut yc)) in to_insert.iter().zip(arr.iter_mut_mzs()) {
            // On each iteration we get [0 ,0 ,0, 0]
            yc.add_at_index(x.0, x.1);
            yc.add_at_index(x.0, x.1);
            // here should be [8.0 ,0 ,0, 0] on iteration 1
            // here should be [0 ,4.0 ,0, 0] on iteration 2
            // I am expecting that out of bounds gets added to the last.
            yc.add_at_index(x.0 + 50, x.1);
            // here should be [8.0 ,0 ,0, 4.0] on iteration 1
            // here should be [0 ,4.0, 0.0, 2.0] on iteration 2
        }
        assert_eq!(arr.arr.values, vec![8.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 2.0]);

        // Since we are resetting two mz major arrays the resulting data should be
        // the same.
        target_arr.try_reset_with(&arr).unwrap();
        assert_eq!(
            target_arr.arr.values,
            vec![8.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 2.0]
        );
    }

    #[test]
    fn test_array2d_reset_resize_smaller() {
        let rt_ms: Arc<[u32]> = [10u32, 20, 30, 40].into();
        let order: Arc<[(String, f64)]> =
            [("FOO".into(), 123.5f64), ("BAR".into(), 456.3f64)].into();
        let mut arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order, rt_ms.len(), 0).unwrap();

        let rt_ms_tgt: Arc<[u32]> = [8678309u32].into();
        let order_tgt: Arc<[(String, f64)]> = [("TOMATO".into(), 123.5f64)].into();
        let target_arr: MzMajorIntensityArray<String, f32> =
            MzMajorIntensityArray::try_new_empty(order_tgt, rt_ms_tgt.len(), 0).unwrap();

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
            MzMajorIntensityArray::try_new_empty(order, rt_ms.len(), 0).unwrap();

        let rt_ms_tgt: Arc<[u32]> = [8678309u32].into();
        let order_tgt: Arc<[(String, f64)]> = [("TOMATO".into(), 123.5f64)].into();
        let mut target_arr: RTMajorIntensityArray<String, f32> =
            RTMajorIntensityArray::try_new_empty(order_tgt, rt_ms_tgt.len(), 0).unwrap();

        assert_eq!(arr.arr.values, vec![0.0; 8]);
        assert_eq!(target_arr.arr.values, vec![0.0]);

        arr.arr.values = vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0];
        target_arr.try_reset_with(&arr).unwrap();

        // Note the change in position
        assert_eq!(arr.arr.values, vec![8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 2.0]);
        assert_eq!(
            target_arr.arr.values,
            vec![8.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 2.0]
        );
    }
}
