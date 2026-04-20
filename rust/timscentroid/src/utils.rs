use half::f16;
use std::ops::RangeInclusive;
use thiserror::Error;

/// Non-negative non-NaN f16 mobility, stored as raw bits.
///
/// Invariant: `bits` encodes a non-negative, non-NaN f16 value.
/// Enforced at construction via `try_new` / `from_f16` — every `MobInt`
/// in flight has already passed the check.
///
/// For values in the valid range, `u16` bit comparison matches `f16`
/// comparison, so `Ord` can be derived from the u16 representation
/// (no f32 widening required). This is what makes the SoA filter
/// hot path SIMD-trivial.
///
/// `#[repr(transparent)]` ensures `Vec<MobInt>` is byte-identical to
/// `Vec<u16>` — zero-cost reinterpretation in either direction when
/// needed (though construction must still go through `try_new`).
#[derive(
    Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
#[repr(transparent)]
pub struct MobInt(u16);

#[derive(Error, Debug)]
pub enum MobIntError {
    #[error("mobility is NaN (bits = {bits:#06x})")]
    Nan { bits: u16 },
    #[error("mobility is negative (bits = {bits:#06x})")]
    Negative { bits: u16 },
}

impl MobInt {
    /// Construct from raw f16 bits. Validates non-negative + non-NaN
    /// via direct bit tests (no f16 widen needed).
    #[inline]
    pub fn try_new(bits: u16) -> Result<Self, MobIntError> {
        // IEEE 754 binary16:
        //   sign: bit 15
        //   exp:  bits 10-14 (0x7C00 mask)
        //   mantissa: bits 0-9 (0x03FF mask)
        // Negative: sign bit set. Mobility must be non-negative.
        if bits & 0x8000 != 0 {
            return Err(MobIntError::Negative { bits });
        }
        // NaN: exp all-ones AND mantissa != 0.
        if (bits & 0x7C00) == 0x7C00 && (bits & 0x03FF) != 0 {
            return Err(MobIntError::Nan { bits });
        }
        Ok(Self(bits))
    }

    #[inline]
    pub fn from_f16(m: f16) -> Result<Self, MobIntError> {
        Self::try_new(m.to_bits())
    }

    /// Free reinterpret — same bits, same register.
    #[inline(always)]
    pub fn to_f16(self) -> f16 {
        f16::from_bits(self.0)
    }

    #[inline(always)]
    pub fn bits(self) -> u16 {
        self.0
    }
}

impl From<MobInt> for f16 {
    #[inline(always)]
    fn from(m: MobInt) -> Self {
        m.to_f16()
    }
}

/// Finds the index range of elements in a sorted slice whose keys fall within the specified range.
///
/// This function performs a binary search on a slice that is sorted by the result of applying
/// the key function to each element. It returns a `Range<usize>` that can be used to slice
/// the original array to get all elements whose keys fall within the inclusive range.
///
/// # Arguments
///
/// * `slice` - A slice sorted by the result of applying `key_fn` to each element
/// * `key_range` - An inclusive range of key values to search for
/// * `key_fn` - A function that extracts the key from each element for comparison
///
/// # Returns
///
/// A `Range<usize>` representing the indices of elements whose keys fall within `key_range`.
/// The range can be empty if no elements match. The returned range can be used directly
/// with slice indexing: `&slice[result_range]`.
///
/// # Time Complexity
///
/// O(log n) where n is the length of the slice.
///
/// # Panics
///
/// This function will panic if the slice is not sorted by the key function.
///
/// # Examples
///
/// ```
/// use timscentroid::utils::binary_search_range_by_key;
/// use std::ops::RangeInclusive;
///
/// #[derive(Debug, PartialEq)]
/// struct Person {
///     name: &'static str,
///     age: u32,
/// }
///
/// // Note: Sorted by age
/// let people = vec![
///     Person { name: "Alice", age: 20 },
///     Person { name: "Bob", age: 25 },
///     Person { name: "Charlie", age: 30 },
///     Person { name: "David", age: 35 },
///     Person { name: "Eve", age: 40 },
/// ];
///
/// // Find people aged 25-35 (inclusive)
/// let age_range = 25..=35;
/// let result_range = binary_search_range_by_key(&people, age_range, |p| p.age);
/// let matching_people = &people[result_range];
///
/// assert_eq!(matching_people.len(), 3);
/// assert_eq!(matching_people[0].name, "Bob");
/// assert_eq!(matching_people[1].name, "Charlie");
/// assert_eq!(matching_people[2].name, "David");
///
/// // Empty range when no matches
/// let empty_range = binary_search_range_by_key(&people, 50..=60, |p| p.age);
/// assert!(people[empty_range].is_empty());
///
/// // Single element match
/// let single_range = binary_search_range_by_key(&people, 20..=20, |p| p.age);
/// assert_eq!(people[single_range.clone()].len(), 1);
/// assert_eq!(people[single_range.clone()][0].name, "Alice");
/// ```
pub fn binary_search_range_by_key<T, K, F>(
    slice: &[T],
    key_range: RangeInclusive<K>,
    key_fn: F,
) -> std::ops::Range<usize>
where
    F: Fn(&T) -> K,
    K: Ord,
{
    let start_idx = slice.partition_point(|x| key_fn(x) < *key_range.start());
    let end_idx = start_idx + slice[start_idx..].partition_point(|x| key_fn(x) <= *key_range.end());

    start_idx..end_idx
}

/// TupleRange represents a range defined by a tuple of two elements (T, T).
///
/// It represents a range as closed-closed [a, b], meaning both endpoints are inclusive.
/// Importantly, it ensures that the first element is always less than or equal to the second,
/// to enforce early exit in cases where a range would be propagated poinlessly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct TupleRange<T: Copy + PartialOrd>(T, T);

#[derive(Error, Debug)]
pub enum TupleRangeError<T: Copy + PartialOrd + std::fmt::Debug> {
    #[error(
        "Expected the first element to be less than or equal to the second, got ({0:?}, {1:?})"
    )]
    ExpectedOrderedRange(T, T),
}

impl<T: Copy + PartialOrd + std::fmt::Debug> TupleRange<T> {
    /// Creates a new `TupleRange` ensuring that the first element
    /// is less than or equal to the second.
    pub fn try_new(left: T, right: T) -> Result<Self, TupleRangeError<T>> {
        // Swap if left > right
        if left > right {
            Err(TupleRangeError::ExpectedOrderedRange(left, right))
        } else {
            Ok(Self(left, right))
        }
    }

    pub fn as_tuple(&self) -> (T, T) {
        (self.0, self.1)
    }

    pub fn as_inclusive_range(&self) -> std::ops::RangeInclusive<T> {
        self.0..=self.1
    }

    pub fn contains(&self, x: T) -> bool {
        self.0 <= x && x <= self.1
    }

    pub fn start(&self) -> T {
        self.0
    }

    pub fn end(&self) -> T {
        self.1
    }

    pub fn map_elems<W: Copy + PartialOrd + std::fmt::Debug>(
        self,
        f: impl Fn(T) -> W,
    ) -> Result<TupleRange<W>, TupleRangeError<W>> {
        TupleRange::try_new(f(self.0), f(self.1))
    }

    pub fn intersects(&self, other: Self) -> bool {
        !(self.end() < other.start() || other.end() < self.start())
    }

    /// True iff `other` is fully contained within `self` (closed interval).
    pub fn encloses(&self, other: Self) -> bool {
        self.start() <= other.start() && other.end() <= self.end()
    }

    pub fn try_intercept(&self, other: Self) -> Option<Self> {
        let left = match self.start().partial_cmp(&other.start()) {
            Some(std::cmp::Ordering::Equal) => self.start(),
            Some(std::cmp::Ordering::Less) => other.start(),
            Some(std::cmp::Ordering::Greater) => self.start(),
            None => self.start(),
        };
        let right = match self.end().partial_cmp(&other.end()) {
            Some(std::cmp::Ordering::Equal) => self.end(),
            Some(std::cmp::Ordering::Less) => self.end(),
            Some(std::cmp::Ordering::Greater) => other.end(),
            None => self.end(),
        };
        if left > right {
            None
        } else {
            Some(Self(left, right))
        }
    }
}

impl<T> AsRef<TupleRange<T>> for TupleRange<T>
where
    T: Copy + PartialOrd,
{
    fn as_ref(&self) -> &TupleRange<T> {
        self
    }
}

impl<T> TryInto<TupleRange<T>> for (T, T)
where
    T: Copy + PartialOrd + std::fmt::Debug,
{
    type Error = TupleRangeError<T>;

    fn try_into(self) -> Result<TupleRange<T>, Self::Error> {
        TupleRange::try_new(self.0, self.1)
    }
}

/// An enum representing a value that can either be restricted to a specific value
/// or be unrestricted (i.e., no restriction).
///
/// In essence is the same as `Option<T>` but with different semantics,
/// since `Option<T>` could mean either no restriction OR no value allowed.
///
/// Mainly meant to be used in cases where I can filter based on a range
/// This allows for cleaner code when dealing with optional restrictions.
#[derive(Debug, Clone, Copy)]
pub enum OptionallyRestricted<T: Copy> {
    Restricted(T),
    Unrestricted,
}

impl<T: Copy> OptionallyRestricted<T> {
    pub fn to_option(&self) -> Option<T> {
        match self {
            OptionallyRestricted::Restricted(x) => Some(*x),
            OptionallyRestricted::Unrestricted => None,
        }
    }

    pub fn map<T2: Copy>(self, f: impl FnOnce(T) -> T2) -> OptionallyRestricted<T2> {
        match self {
            OptionallyRestricted::Restricted(x) => OptionallyRestricted::Restricted(f(x)),
            OptionallyRestricted::Unrestricted => OptionallyRestricted::Unrestricted,
        }
    }

    pub fn is_unrestricted_or(&self, f: impl FnOnce(&T) -> bool) -> bool
    where
        T: PartialEq,
    {
        match self {
            OptionallyRestricted::Restricted(x) => f(x),
            OptionallyRestricted::Unrestricted => true,
        }
    }

    pub fn unwrap_or(self, default: T) -> T {
        match self {
            OptionallyRestricted::Restricted(x) => x,
            OptionallyRestricted::Unrestricted => default,
        }
    }

    pub fn unwrap(self) -> T {
        match self {
            OptionallyRestricted::Restricted(x) => x,
            OptionallyRestricted::Unrestricted => panic!("Called unwrap on an Unrestricted value"),
        }
    }

    pub fn expect(self, msg: &str) -> T {
        match self {
            OptionallyRestricted::Restricted(x) => x,
            OptionallyRestricted::Unrestricted => panic!("{}", msg),
        }
    }
}

impl<T> AsRef<OptionallyRestricted<T>> for OptionallyRestricted<T>
where
    T: Copy,
{
    fn as_ref(&self) -> &OptionallyRestricted<T> {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mi(v: f32) -> Result<MobInt, MobIntError> {
        MobInt::from_f16(f16::from_f32(v))
    }

    #[test]
    fn mobint_accepts_non_negative_non_nan() {
        for v in [0.0, f32::MIN_POSITIVE, 0.85, 1.5, 1000.0, f32::INFINITY] {
            mi(v).unwrap_or_else(|e| panic!("{v}: {e}"));
        }
    }

    #[test]
    fn mobint_rejects_negative() {
        for v in [-0.0, -1.0e-5, -0.85, -1.5, f32::NEG_INFINITY] {
            assert!(matches!(mi(v), Err(MobIntError::Negative { .. })), "{v}");
        }
    }

    #[test]
    fn mobint_rejects_nan() {
        assert!(matches!(mi(f32::NAN), Err(MobIntError::Nan { .. })));
    }

    #[test]
    fn mobint_f16_roundtrip() {
        let m = f16::from_f32(0.85);
        let mob = MobInt::from_f16(m).unwrap();
        assert_eq!(mob.to_f16(), m);
        assert_eq!(mob.bits(), m.to_bits());
    }

    #[test]
    fn mobint_ord_matches_f16() {
        // For validated (non-negative, non-NaN) f16, u16 bit order == f16 order.
        let values = [0.0f32, 0.25, 0.5, 0.9, 1.0, 1.25, 10.0, 1000.0];
        for (i, a) in values.iter().enumerate() {
            for b in &values[i..] {
                let fa = f16::from_f32(*a);
                let fb = f16::from_f32(*b);
                let ma = MobInt::from_f16(fa).unwrap();
                let mb = MobInt::from_f16(fb).unwrap();
                assert_eq!(ma.cmp(&mb), fa.partial_cmp(&fb).unwrap(), "{a} vs {b}");
            }
        }
    }

    #[test]
    fn test_slice_search() {
        let input = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let range = 3..=7;
        let result = binary_search_range_by_key(&input, range, |&x| x);
        assert_eq!(result, 2..7);

        assert_eq!(&input[result], &[3, 4, 5, 6, 7]);
    }

    #[test]
    fn test_slice_search_repeats() {
        let input = vec![
            1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 7, 7, 7, 7, 8, 9, 10,
        ];
        let range = 3..=7;
        let result = binary_search_range_by_key(&input, range, |&x| x);
        assert_eq!(result, 2..18);

        assert_eq!(
            &input[result],
            &[3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 7, 7, 7, 7]
        );
    }
}
