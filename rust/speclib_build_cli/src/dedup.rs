use bloomfilter::Bloom;
use std::collections::HashMap;
use timsseek::DigestSlice;

pub struct PeptideDedup;

impl PeptideDedup {
    /// Deduplicate DigestSlices using bloom filter fast-pass + bucket HashMap.
    /// `estimated_count` used to preallocate bloom + HashMap capacity.
    pub fn dedup(slices: Vec<DigestSlice>, estimated_count: usize) -> Vec<DigestSlice> {
        let est = estimated_count.max(64);
        let mut bloom = Bloom::new_for_fp_rate(est, 0.01);
        let mut buckets: HashMap<(u8, u16), Vec<DigestSlice>> = HashMap::with_capacity(est / 20);

        for slice in slices {
            let seq = slice.as_str().as_bytes();
            let key = (seq[0], seq.len() as u16);

            if !bloom.check(seq) {
                // Definitely new
                bloom.set(seq);
                buckets.entry(key).or_default().push(slice);
            } else {
                // Maybe seen — check bucket
                let bucket = buckets.entry(key).or_default();
                if !bucket.iter().any(|existing| existing.as_str().as_bytes() == seq) {
                    bucket.push(slice);
                }
            }
        }

        buckets.into_values().flat_map(|v| v.into_iter()).collect()
    }

    /// Estimate unique peptide count from protein lengths for preallocation.
    pub fn estimate_from_proteins(total_aa: usize, missed_cleavages: usize) -> usize {
        (total_aa / 10) * (1 + missed_cleavages)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use timsseek::models::DecoyMarking;

    fn make_slice(seq: &str) -> DigestSlice {
        let s: Arc<str> = seq.into();
        let len = s.len() as u16;
        DigestSlice::new(s, 0..len, DecoyMarking::Target, 0)
    }

    #[test]
    fn test_dedup_removes_duplicates() {
        // 5 slices, 2 dupes → 3 unique
        let slices = vec![
            make_slice("PEPTIDE"),
            make_slice("LAGER"),
            make_slice("PEPTIDE"), // dupe
            make_slice("TOMATO"),
            make_slice("LAGER"), // dupe
        ];
        let result = PeptideDedup::dedup(slices, 10);
        assert_eq!(result.len(), 3);
        // All three unique sequences are present
        let seqs: Vec<&str> = result.iter().map(|s| s.as_str()).collect();
        assert!(seqs.contains(&"PEPTIDE"));
        assert!(seqs.contains(&"LAGER"));
        assert!(seqs.contains(&"TOMATO"));
    }

    #[test]
    fn test_dedup_same_len_different_seq() {
        // Same length, different content → both kept
        let slices = vec![
            make_slice("ABCDE"),
            make_slice("VWXYZ"),
        ];
        let result = PeptideDedup::dedup(slices, 10);
        assert_eq!(result.len(), 2);
        let seqs: Vec<&str> = result.iter().map(|s| s.as_str()).collect();
        assert!(seqs.contains(&"ABCDE"));
        assert!(seqs.contains(&"VWXYZ"));
    }

    #[test]
    fn test_dedup_shared_protein_backbone() {
        // Two slices from same Arc, same range → deduped
        let backbone: Arc<str> = "PEPTIDEPINKTOMATO".into();
        let len = backbone.len() as u16;
        let slice1 = DigestSlice::new(backbone.clone(), 0..len, DecoyMarking::Target, 0);
        let slice2 = DigestSlice::new(backbone.clone(), 0..len, DecoyMarking::Target, 0);
        let slices = vec![slice1, slice2];
        let result = PeptideDedup::dedup(slices, 10);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].as_str(), "PEPTIDEPINKTOMATO");
    }

    #[test]
    fn test_dedup_empty() {
        // Empty input → empty output
        let result = PeptideDedup::dedup(vec![], 0);
        assert!(result.is_empty());
    }
}
