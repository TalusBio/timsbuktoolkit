use crate::traits::KeyLike;
use serde::{
    Deserialize,
    Serialize,
};
use std::sync::Arc;

/// A struct that represents an elution group.
///
/// The elution group is a single precursor ion that is framented.
/// The fragments m/z values are stored in a hashmap where the key is
/// the generic type `T` and the value is the fragment m/z.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ElutionGroup<T: KeyLike> {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    // This i8 might not be the best for performance due to alignment ...
    pub precursors: Arc<[(i8, f64)]>,
    pub fragments: Arc<[(T, f64)]>,
}

impl<T: KeyLike> ElutionGroup<T> {
    pub fn get_precursor_mz_limits(&self) -> (f64, f64) {
        let mut min_precursor_mz = f64::MAX;
        let mut max_precursor_mz = f64::MIN;
        for (isotope_index, precursor_mz) in self.precursors.iter() {
            if *isotope_index < 0 {
                continue;
            }
            if *precursor_mz < min_precursor_mz {
                min_precursor_mz = *precursor_mz;
            }
            if *precursor_mz > max_precursor_mz {
                max_precursor_mz = *precursor_mz;
            }
        }
        (min_precursor_mz, max_precursor_mz)
    }

    pub fn get_monoisotopic_precursor_mz(&self) -> Option<f64> {
        self.precursors
            .iter()
            .find(|(isotope_index, _)| *isotope_index == 0)
            .map(|(_, precursor_mz)| *precursor_mz)
    }

    pub fn iter_precursors(&self) -> impl Iterator<Item = &(i8, f64)> {
        self.precursors.iter()
    }

    pub fn iter_fragments(&self) -> impl Iterator<Item = &(T, f64)> {
        self.fragments.iter()
    }
}
