use serde::{
    Deserialize,
    Serialize,
};
use timsquery::KeyLike;

/// Fragment and precursor intensity pairs. Reused in place on the scoring hot
/// path, so the buffer is allocated once per worker rather than per item.
pub type FragmentIntensityVec<T> = Vec<(T, f32)>;
pub type PrecursorIntensityVec = Vec<(i8, f32)>;

/// Linear lookup for a `(key, value)` slice. Used throughout scoring for
/// `ExpectedIntensities`-shaped arrays, which hold a dozen entries in a
/// typical library, so a scan beats any index.
#[inline]
pub fn linear_get<K: PartialEq, V: Copy>(entries: &[(K, V)], key: &K) -> Option<V> {
    entries.iter().find(|(k, _)| k == key).map(|(_, v)| *v)
}

/// Rejection returned by [`ExpectedIntensities::try_from_pairs`] when the
/// input contains a repeated fragment or precursor key. Keys are expected
/// to be unique -- `linear_get`/`remove_*` return only the first match and
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
/// A `Vec` per field, each entry keeping its key alongside the intensity for
/// diagnostics and by-key lookup.
///
/// On the scoring path the buffers are sized once per worker via
/// [`with_capacity`](Self::with_capacity) and reused per item through
/// [`refill_from_pairs`](Self::refill_from_pairs), so scoring does not
/// allocate.
///
/// # Invariants
///
/// Keys within each field must be unique. [`linear_get`] and
/// [`ExpectedIntensities::remove_fragment`]/`remove_precursor` return only
/// the first match; a duplicate silently hides data. All construction paths
/// go through [`ExpectedIntensities::try_from_pairs`] which enforces this.
/// Direct struct-literal construction bypasses the check -- avoid it outside
/// tests where uniqueness is obvious.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedIntensities<T: KeyLike + Default> {
    pub fragment_intensities: FragmentIntensityVec<T>,
    pub precursor_intensities: PrecursorIntensityVec,
}

/// Precursor isotope slots. `isotope_dist_or_averagine` returns `[f32; 3]` and
/// every `IsotopeStrategy::FromComposition` site sets `n_isotopes: 3`, so the
/// envelope is this long by construction.
pub const PRECURSOR_ENVELOPE_LEN: usize = 3;

impl<T: KeyLike + Default> ExpectedIntensities<T> {
    /// Allocate both buffers up front. Callers on the scoring path pass the
    /// library's maximum fragment count so that
    /// [`refill_from_pairs`](Self::refill_from_pairs) never grows a buffer
    /// mid-run.
    pub fn with_capacity(fragments: usize) -> Self {
        Self {
            fragment_intensities: Vec::with_capacity(fragments),
            precursor_intensities: Vec::with_capacity(PRECURSOR_ENVELOPE_LEN),
        }
    }
}

impl<T: KeyLike + Default> Default for ExpectedIntensities<T> {
    fn default() -> Self {
        Self {
            fragment_intensities: Vec::new(),
            precursor_intensities: Vec::new(),
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
        out.refill_from_pairs(frags, precs)?;
        Ok(out)
    }

    /// [`try_from_pairs`](Self::try_from_pairs) in place, keeping the buffers
    /// already allocated. Same uniqueness check.
    ///
    /// Use this on the scoring hot path: assigning a fresh value drops the
    /// buffer the next item would have reused.
    ///
    /// On a duplicate key both fields are left empty, so a rejected refill
    /// never leaves a partially filled set behind.
    pub fn refill_from_pairs<FI, PI>(
        &mut self,
        frags: FI,
        precs: PI,
    ) -> Result<(), DuplicateKeyError>
    where
        FI: IntoIterator<Item = (T, f32)>,
        PI: IntoIterator<Item = (i8, f32)>,
    {
        self.fragment_intensities.clear();
        self.precursor_intensities.clear();
        for (k, v) in frags {
            if self.get_fragment(&k).is_some() {
                let key = format!("{:?}", k);
                self.fragment_intensities.clear();
                return Err(DuplicateKeyError {
                    which: "fragment",
                    key,
                });
            }
            self.fragment_intensities.push((k, v));
        }
        for (k, v) in precs {
            if self.get_precursor(k).is_some() {
                let key = k.to_string();
                self.fragment_intensities.clear();
                self.precursor_intensities.clear();
                return Err(DuplicateKeyError {
                    which: "precursor",
                    key,
                });
            }
            self.precursor_intensities.push((k, v));
        }
        Ok(())
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

    /// In-place copy from `other`, reusing each buffer when it is already big
    /// enough. Prefer this over `*self = other.clone()` on the scoring hot
    /// path: the default `Clone::clone_from` drops the destination buffer
    /// before reallocating.
    pub fn clone_from_ref(&mut self, other: &Self)
    where
        T: Clone,
    {
        self.fragment_intensities.clear();
        self.fragment_intensities
            .extend(other.fragment_intensities.iter().cloned());
        self.precursor_intensities.clear();
        self.precursor_intensities
            .extend(other.precursor_intensities.iter().cloned());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use micromzpaf::IonAnnot;

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
    fn with_capacity_survives_a_refill_without_growing() {
        let mut ei = ExpectedIntensities::<IonAnnot>::with_capacity(32);
        let cap = ei.fragment_intensities.capacity();
        let frags: Vec<_> = (1..=32).map(|i| (yi(&format!("y{i}")), i as f32)).collect();
        ei.refill_from_pairs(frags, [(0i8, 0.5), (1, 0.3), (2, 0.1)])
            .unwrap();
        assert_eq!(
            ei.fragment_intensities.capacity(),
            cap,
            "a library-sized refill must not reallocate"
        );
        assert_eq!(ei.precursor_intensities.capacity(), PRECURSOR_ENVELOPE_LEN);
    }

    #[test]
    fn refill_replaces_previous_contents() {
        let mut ei = ExpectedIntensities::try_from_pairs(
            [(yi("y1"), 1.0), (yi("y2"), 2.0), (yi("y3"), 3.0)],
            [(0i8, 0.5)],
        )
        .unwrap();
        ei.refill_from_pairs([(yi("b1"), 9.0)], [(2i8, 0.1)])
            .unwrap();
        assert_eq!(ei.fragment_len(), 1);
        assert_eq!(ei.precursor_len(), 1);
        assert_eq!(ei.get_fragment(&yi("b1")), Some(9.0));
        assert_eq!(ei.get_fragment(&yi("y1")), None, "stale key must be gone");
        assert_eq!(ei.get_precursor(2), Some(0.1));
        assert_eq!(ei.get_precursor(0), None, "stale isotope must be gone");
    }

    #[test]
    fn refill_leaves_nothing_behind_on_duplicate() {
        let mut ei = ExpectedIntensities::try_from_pairs([(yi("y1"), 1.0)], [(0i8, 0.5)]).unwrap();
        let err = ei
            .refill_from_pairs([(yi("b1"), 1.0), (yi("b1"), 2.0)], [(0i8, 0.5)])
            .unwrap_err();
        assert_eq!(err.which, "fragment");
        assert_eq!(ei.fragment_len(), 0, "a rejected refill must not half-fill");
        assert_eq!(ei.precursor_len(), 0);

        let err = ei
            .refill_from_pairs([(yi("b1"), 1.0)], [(0i8, 0.5), (0i8, 0.6)])
            .unwrap_err();
        assert_eq!(err.which, "precursor");
        assert_eq!(ei.fragment_len(), 0);
        assert_eq!(ei.precursor_len(), 0);
    }

    /// The point of the method: a spilled buffer is reused rather than freed
    /// and reallocated. Assigning a fresh value would drop this capacity.
    #[test]
    fn refill_reuses_the_buffer() {
        let many: Vec<_> = (1..=18).map(|i| (yi(&format!("y{i}")), i as f32)).collect();
        let mut ei = ExpectedIntensities::<IonAnnot>::default();
        ei.refill_from_pairs(many.iter().copied(), []).unwrap();
        let cap = ei.fragment_intensities.capacity();
        assert!(cap >= many.len());

        ei.refill_from_pairs([(yi("y1"), 1.0)], []).unwrap();
        assert_eq!(
            ei.fragment_intensities.capacity(),
            cap,
            "refill must keep the buffer it already had"
        );
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
