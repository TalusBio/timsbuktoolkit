use super::decoy::DecoyMarking;
use serde::Serialize;
use std::collections::HashSet;
use std::fmt::Display as FmtDisplay;
use std::ops::Range;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinSlice {
    ref_seq: Arc<str>,
    range: Range<u16>,
    pub decoy: DecoyMarking,
    pub decoy_group: u32,
}

impl Serialize for ProteinSlice {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let local_str = Into::<String>::into(self.clone());
        serializer.serialize_str(local_str.as_str())
    }
}

impl ProteinSlice {
    pub fn new(
        ref_seq: Arc<str>,
        range: Range<u16>,
        decoy: DecoyMarking,
        decoy_group: u32,
    ) -> Self {
        Self {
            ref_seq,
            range,
            decoy,
            decoy_group,
        }
    }

    pub fn from_string(seq: String, decoy: bool, decoy_group: u32) -> ProteinSlice {
        let len = seq.len();
        ProteinSlice {
            ref_seq: seq.into(),
            range: 0..(len as u16),
            decoy: if decoy {
                DecoyMarking::ReversedDecoy
            } else {
                DecoyMarking::Target
            },
            decoy_group,
        }
    }

    pub fn len(&self) -> usize {
        self.range.len()
    }

    pub fn is_empty(&self) -> bool {
        self.range.is_empty()
    }

    /// Returns the peptide sequence as a string slice without allocation.
    pub fn as_str(&self) -> &str {
        &self.ref_seq[self.range.start as usize..self.range.end as usize]
    }

    pub fn is_decoy(&self) -> bool {
        self.decoy.is_decoy()
    }
}

impl FmtDisplay for ProteinSlice {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let local_str: String = self.clone().into();
        write!(f, "{}", local_str)
    }
}

impl From<ProteinSlice> for String {
    fn from(x: ProteinSlice) -> Self {
        x.ref_seq.as_ref()[(x.range.start as usize)..(x.range.end as usize)].to_string()
    }
}

pub fn deduplicate_digests(mut digest_slices: Vec<ProteinSlice>) -> Vec<ProteinSlice> {
    let mut seen = HashSet::new();
    digest_slices.retain(|x| {
        let local_str: String = x.clone().into();
        let is_first = !seen.contains(&local_str);
        seen.insert(local_str);
        is_first
    });
    digest_slices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deduplicate_digests() {
        let seq: Arc<str> = "PEPTIDEPINKTOMATOTOMATO".into();
        let seq2: Arc<str> = "PEPTIDEPINKTOMATO".into();
        let seq2_rep: Arc<str> = "PEPTIDEPINKTOMATO".into();
        let digests: Vec<ProteinSlice> = vec![
            ProteinSlice {
                ref_seq: seq.clone(),
                range: (0_u16..seq.as_ref().len() as u16),
                decoy: DecoyMarking::Target,
                decoy_group: 1,
            },
            ProteinSlice {
                ref_seq: seq.clone(),
                range: (0_u16..seq2.as_ref().len() as u16), // Note the short length
                decoy: DecoyMarking::Target,
                decoy_group: 1,
            },
            ProteinSlice {
                ref_seq: seq2.clone(),
                range: (0_u16..seq2.as_ref().len() as u16),
                decoy: DecoyMarking::Target,
                decoy_group: 1,
            },
            ProteinSlice {
                ref_seq: seq2_rep.clone(),
                range: (0_u16..seq2_rep.as_ref().len() as u16),
                decoy: DecoyMarking::Target,
                decoy_group: 1,
            },
        ];
        let deduped = deduplicate_digests(digests);
        assert_eq!(deduped.len(), 2);
        assert_eq!(deduped[0].len(), seq.as_ref().len());
        assert_eq!(deduped[1].len(), seq2.as_ref().len());
    }

    #[test]
    fn test_as_str() {
        let seq: Arc<str> = "PEPTIDEPINKTOMATO".into();
        let slice = ProteinSlice::new(seq.clone(), 0..7, DecoyMarking::Target, 0);
        assert_eq!(slice.as_str(), "PEPTIDE");
        assert_eq!(slice.as_str().as_bytes(), b"PEPTIDE");
    }

    #[test]
    fn test_from_string() {
        let seq = "PEPTIDEPINKTOMATO".to_string();
        let expect_len = seq.len();
        let digests = ProteinSlice::from_string(seq.clone(), false, 1);
        let digests_decoy = ProteinSlice::from_string(seq.clone(), true, 2);
        assert_eq!(digests.len(), expect_len);
        assert_eq!(digests.decoy, DecoyMarking::Target);
        assert_eq!(digests_decoy.len(), expect_len);
        assert_eq!(digests_decoy.decoy, DecoyMarking::ReversedDecoy);

        // Check that converting back to string does not reverse the sequence
        assert_eq!(Into::<String>::into(digests.clone()), seq);
        assert_eq!(Into::<String>::into(digests_decoy.clone()), seq);
    }
}
