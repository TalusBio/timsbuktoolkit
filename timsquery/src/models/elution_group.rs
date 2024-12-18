use serde::{
    Deserialize,
    Serialize,
};
use std::collections::HashMap;
use std::hash::Hash;

/// A struct that represents an elution group.
///
/// The elution group is a single precursor ion that is framented.
/// The fragments m/z values are stored in a hashmap where the key is
/// the generic type `T` and the value is the fragment m/z.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ElutionGroup<T: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursor_mzs: Vec<f64>,
    pub fragment_mzs: HashMap<T, f64>,
    pub expected_fragment_intensity: Option<HashMap<T, f32>>,
    pub expected_precursor_intensity: Option<Vec<f32>>,
}

impl<T: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> ElutionGroup<T> {
    pub fn get_precursor_mz_limits(&self) -> (f64, f64) {
        let mut min_precursor_mz = f64::MAX;
        let mut max_precursor_mz = f64::MIN;
        for precursor_mz in self.precursor_mzs.iter() {
            if *precursor_mz < min_precursor_mz {
                min_precursor_mz = *precursor_mz;
            }
            if *precursor_mz > max_precursor_mz {
                max_precursor_mz = *precursor_mz;
            }
        }
        (min_precursor_mz, max_precursor_mz)
    }
}
