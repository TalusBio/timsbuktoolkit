//! A minimal fixed-size bitset. Hand-rolled rather than pulling in
//! `fixedbitset`: three operations, and the crate's dependency list is a
//! deliberate constraint.

pub struct BitSet {
    words: Vec<u64>,
    bits: usize,
}

impl BitSet {
    pub fn new(bits: usize) -> Self {
        Self {
            words: vec![0; bits.div_ceil(64)],
            bits,
        }
    }

    pub fn clear(&mut self) {
        self.words.fill(0);
    }

    pub fn set(&mut self, i: usize) {
        if i < self.bits {
            self.words[i / 64] |= 1u64 << (i % 64);
        }
    }

    pub fn get(&self, i: usize) -> bool {
        i < self.bits && (self.words[i / 64] >> (i % 64)) & 1 == 1
    }

    pub fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    #[cfg(test)]
    pub fn word_capacity(&self) -> usize {
        self.words.capacity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_and_get_round_trip_across_word_boundaries() {
        let mut b = BitSet::new(130);
        for i in [0usize, 63, 64, 65, 129] {
            b.set(i);
        }
        for i in [0usize, 63, 64, 65, 129] {
            assert!(b.get(i), "bit {i} should be set");
        }
        assert!(!b.get(1));
        assert!(!b.get(128));
        assert_eq!(b.count_ones(), 5);
    }

    #[test]
    fn clear_resets_without_shrinking() {
        let mut b = BitSet::new(130);
        let cap = b.word_capacity();
        b.set(100);
        b.clear();
        assert!(!b.get(100));
        assert_eq!(b.count_ones(), 0);
        assert_eq!(b.word_capacity(), cap);
    }

    #[test]
    fn out_of_range_get_is_false_not_a_panic() {
        let b = BitSet::new(10);
        assert!(!b.get(1000));
    }
}
