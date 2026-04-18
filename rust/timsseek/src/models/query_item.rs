use super::DigestSlice;
use micromzpaf::IonAnnot;
use serde::{
    Deserialize,
    Serialize,
};
use timsquery::tinyvec::{
    TinyVec,
    tiny_vec,
};
use timsquery::{
    KeyLike,
    TimsElutionGroup,
};

/// Inline capacity for fragment and precursor intensity pairs. Mirrors the
/// `fragment_labels`/`precursor_labels` capacity used by `TimsElutionGroup`
/// so typical peptides stay fully stack-resident and no heap allocation
/// happens on clone/mutate.
pub const INLINE_FRAG_CAPACITY: usize = 13;
pub const INLINE_PREC_CAPACITY: usize = 13;

pub type FragmentIntensityVec<T> = TinyVec<[(T, f32); INLINE_FRAG_CAPACITY]>;
pub type PrecursorIntensityVec = TinyVec<[(i8, f32); INLINE_PREC_CAPACITY]>;

/// Linear lookup for a `(key, value)` slice. Used throughout scoring for
/// `ExpectedIntensities`-shaped arrays whose length is bounded by
/// [`INLINE_FRAG_CAPACITY`] / [`INLINE_PREC_CAPACITY`].
#[inline]
pub fn linear_get<K: PartialEq, V: Copy>(entries: &[(K, V)], key: &K) -> Option<V> {
    entries.iter().find(|(k, _)| k == key).map(|(_, v)| *v)
}

/// Expected (theoretical / predicted) precursor and fragment intensities for
/// a peptide query.
///
/// Stored as `TinyVec<[(K, f32); 13]>` per field: each entry keeps its key
/// alongside the intensity for diagnostics and by-key lookup. For typical
/// peptides (<=13 fragments / <=13 precursors) the backing storage is inline,
/// so `clone()` is a stack memcpy and mutation is in-place swap-remove —
/// zero heap allocations on the hot path.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedIntensities<T: KeyLike + Default> {
    pub fragment_intensities: FragmentIntensityVec<T>,
    pub precursor_intensities: PrecursorIntensityVec,
}

impl<T: KeyLike + Default> Default for ExpectedIntensities<T> {
    fn default() -> Self {
        Self {
            fragment_intensities: TinyVec::new(),
            precursor_intensities: TinyVec::new(),
        }
    }
}

impl<T: KeyLike + Default> ExpectedIntensities<T> {
    pub fn fragment_len(&self) -> usize {
        self.fragment_intensities.len()
    }

    pub fn precursor_len(&self) -> usize {
        self.precursor_intensities.len()
    }

    pub fn get_fragment(&self, key: &T) -> Option<f32> {
        linear_get(&self.fragment_intensities, key)
    }

    pub fn get_precursor(&self, key: i8) -> Option<f32> {
        linear_get(&self.precursor_intensities, &key)
    }

    /// Swap-remove by key. Returns the removed value if present.
    pub fn remove_fragment(&mut self, key: &T) -> Option<f32> {
        let idx = self
            .fragment_intensities
            .iter()
            .position(|(k, _)| k == key)?;
        Some(self.fragment_intensities.swap_remove(idx).1)
    }

    pub fn remove_precursor(&mut self, key: i8) -> Option<f32> {
        let idx = self
            .precursor_intensities
            .iter()
            .position(|(k, _)| *k == key)?;
        Some(self.precursor_intensities.swap_remove(idx).1)
    }
}

#[derive(Debug, Clone)]
pub struct QueryItemToScore {
    // Kinda hate this
    pub digest: DigestSlice,
    pub query: TimsElutionGroup<IonAnnot>,
    pub expected_intensity: ExpectedIntensities<IonAnnot>,
}

impl QueryItemToScore {
    pub fn sample() -> Self {
        let eg = TimsElutionGroup::builder()
            .id(42)
            .mobility_ook0(0.75)
            .rt_seconds(123.4)
            .fragment_labels(
                [
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                    IonAnnot::try_from("y3").unwrap(),
                    IonAnnot::try_from("y4").unwrap(),
                ]
                .as_slice()
                .into(),
            )
            .fragment_mzs(vec![450.0, 650.5, 751.0, 751.5])
            .precursor_labels(tiny_vec!(-1, 0, 1, 2))
            .precursor(450.5, 2)
            .try_build()
            .unwrap();

        let ei = ExpectedIntensities {
            fragment_intensities: tiny_vec!(
                [(IonAnnot, f32); INLINE_FRAG_CAPACITY] =>
                (IonAnnot::try_from("y1").unwrap(), 1.0),
                (IonAnnot::try_from("y2").unwrap(), 1.0),
                (IonAnnot::try_from("y3").unwrap(), 1.0),
                (IonAnnot::try_from("y4").unwrap(), 1.0),
            ),
            precursor_intensities: tiny_vec!(
                [(i8, f32); INLINE_PREC_CAPACITY] =>
                (-1, 0.5),
                (0, 1.0),
                (1, 0.8),
                (2, 0.3),
            ),
        };
        let pepseq = "PEPTIDEPINKPEPTIDE".into();
        let digest = DigestSlice::from_string(pepseq, false, 1);
        let query = eg;
        let expected_intensity = ei;
        QueryItemToScore {
            digest,
            query,
            expected_intensity,
        }
    }
}
