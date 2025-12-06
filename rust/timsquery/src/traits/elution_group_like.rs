use super::KeyLike;

/// A trait representing the common behavior of elution groups.
///
/// This trait is meant to abstract the behavior of elution groups,
/// allowing for different implementations to be used interchangeably.
/// (The main idea for that is to allow different underlying data structures
/// to represent elution groups while providing a consistent interface.
/// Which should be nice to reduce strain on serialization/deserialization,
/// and have flexibility on performance optimizations.)
pub trait TimsElutionGroupLike<T: KeyLike> {
    fn id(&self) -> u64;
    fn precursor_count(&self) -> usize;
    fn fragment_count(&self) -> usize;
    fn rt_seconds(&self) -> f32;
    fn mobility_ook0(&self) -> f32;

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> + '_;
    fn iter_fragments<'a>(&'a self) -> impl Iterator<Item = (&'a T, f64)> + 'a
    where
        T: 'a;
}
