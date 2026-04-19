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

/// Rejection returned by [`ExpectedIntensities::try_from_pairs`] when the
/// input contains a repeated fragment or precursor key. Keys are expected
/// to be unique — `linear_get`/`remove_*` return only the first match and
/// silently masking a duplicate corrupts downstream scoring.
#[derive(Debug, Clone)]
pub struct DuplicateKeyError {
    pub which: &'static str,
    pub key: String,
}

impl std::fmt::Display for DuplicateKeyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "duplicate {} key: {}", self.which, self.key)
    }
}

impl std::error::Error for DuplicateKeyError {}

/// Expected (theoretical / predicted) precursor and fragment intensities for
/// a peptide query.
///
/// Stored as `TinyVec<[(K, f32); 13]>` per field: each entry keeps its key
/// alongside the intensity for diagnostics and by-key lookup. For typical
/// peptides (<=13 fragments / <=13 precursors) the backing storage is inline,
/// so `clone()` is a stack memcpy and mutation is in-place swap-remove —
/// zero heap allocations on the hot path.
///
/// # Invariants
///
/// Keys within each field must be unique. [`linear_get`] and
/// [`ExpectedIntensities::remove_fragment`]/[`remove_precursor`] return only
/// the first match; a duplicate silently hides data. All construction paths
/// go through [`ExpectedIntensities::try_from_pairs`] which enforces this.
/// Direct struct-literal construction bypasses the check — avoid it outside
/// tests where uniqueness is obvious.
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

impl<T: KeyLike + Default + std::fmt::Debug> ExpectedIntensities<T> {
    /// Construct from fragment and precursor pair iterators, erroring on any
    /// repeated key in either input. Preferred entry point for all library
    /// load paths (speclib ndjson/msgpack, DIA-NN/Spectronaut/Skyline TSV).
    pub fn try_from_pairs<FI, PI>(frags: FI, precs: PI) -> Result<Self, DuplicateKeyError>
    where
        FI: IntoIterator<Item = (T, f32)>,
        PI: IntoIterator<Item = (i8, f32)>,
    {
        let mut out = Self::default();
        for (k, v) in frags {
            if out.get_fragment(&k).is_some() {
                return Err(DuplicateKeyError {
                    which: "fragment",
                    key: format!("{:?}", k),
                });
            }
            out.fragment_intensities.push((k, v));
        }
        for (k, v) in precs {
            if out.get_precursor(k).is_some() {
                return Err(DuplicateKeyError {
                    which: "precursor",
                    key: k.to_string(),
                });
            }
            out.precursor_intensities.push((k, v));
        }
        Ok(out)
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

    /// Keep only the `n` fragments with largest predicted intensity,
    /// partial-sorting in place and truncating. Zero alloc. O(len).
    /// No-op if `len <= n`. Surviving order is unspecified.
    pub fn retain_top_n_fragments(&mut self, n: usize) {
        if self.fragment_intensities.len() <= n {
            return;
        }
        // Partition by intensity descending: the `n`th element splits top-n / rest.
        self.fragment_intensities
            .select_nth_unstable_by(n, |a, b| b.1.total_cmp(&a.1));
        self.fragment_intensities.truncate(n);
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

#[cfg(test)]
mod tests {
    use super::*;

    fn yi(s: &str) -> IonAnnot {
        IonAnnot::try_from(s).unwrap()
    }

    #[test]
    fn try_from_pairs_empty_is_ok() {
        let out: ExpectedIntensities<IonAnnot> =
            ExpectedIntensities::try_from_pairs(std::iter::empty(), std::iter::empty()).unwrap();
        assert_eq!(out.fragment_len(), 0);
        assert_eq!(out.precursor_len(), 0);
    }

    #[test]
    fn try_from_pairs_unique_ok() {
        let out: ExpectedIntensities<IonAnnot> = ExpectedIntensities::try_from_pairs(
            [(yi("y1"), 0.5), (yi("y2"), 0.7)],
            [(0i8, 1.0), (1i8, 0.3)],
        )
        .unwrap();
        assert_eq!(out.get_fragment(&yi("y1")), Some(0.5));
        assert_eq!(out.get_precursor(1), Some(0.3));
    }

    #[test]
    fn try_from_pairs_duplicate_fragment_errs() {
        let err = ExpectedIntensities::<IonAnnot>::try_from_pairs(
            [(yi("y1"), 0.5), (yi("y1"), 0.9)],
            std::iter::empty::<(i8, f32)>(),
        )
        .expect_err("duplicate fragment must error");
        assert_eq!(err.which, "fragment");
    }

    #[test]
    fn try_from_pairs_duplicate_precursor_errs() {
        let err = ExpectedIntensities::<IonAnnot>::try_from_pairs(
            std::iter::empty::<(IonAnnot, f32)>(),
            [(0i8, 1.0), (0i8, 2.0)],
        )
        .expect_err("duplicate precursor must error");
        assert_eq!(err.which, "precursor");
    }

    #[test]
    fn retain_top_n_noop_when_under_cap() {
        let mut ei = ExpectedIntensities::<IonAnnot>::try_from_pairs(
            [(yi("y1"), 0.5), (yi("y2"), 0.7)],
            std::iter::empty::<(i8, f32)>(),
        )
        .unwrap();
        ei.retain_top_n_fragments(8);
        assert_eq!(ei.fragment_len(), 2);
    }

    #[test]
    fn retain_top_n_keeps_highest() {
        let mut ei = ExpectedIntensities::<IonAnnot>::try_from_pairs(
            [
                (yi("y1"), 0.1),
                (yi("y2"), 0.9),
                (yi("y3"), 0.5),
                (yi("y4"), 0.2),
            ],
            std::iter::empty::<(i8, f32)>(),
        )
        .unwrap();
        ei.retain_top_n_fragments(2);
        assert_eq!(ei.fragment_len(), 2);
        // Both survivors must have intensity >= 0.5
        for (_, v) in ei.fragment_intensities.iter() {
            assert!(*v >= 0.5, "unexpected survivor {}", v);
        }
        // y2 (0.9) and y3 (0.5) must both be present
        assert!(ei.get_fragment(&yi("y2")).is_some());
        assert!(ei.get_fragment(&yi("y3")).is_some());
    }

    #[test]
    fn try_from_pairs_checks_fragments_before_precursors() {
        // Fragment dup should surface even when precursor block also has a dup.
        let err = ExpectedIntensities::<IonAnnot>::try_from_pairs(
            [(yi("y1"), 0.5), (yi("y1"), 0.5)],
            [(0i8, 1.0), (0i8, 2.0)],
        )
        .expect_err("should fail on fragment dup first");
        assert_eq!(err.which, "fragment");
    }
}
