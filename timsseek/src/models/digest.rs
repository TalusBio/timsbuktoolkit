use super::decoy::{
    DecoyMarking,
    as_decoy_string,
};
use serde::Serialize;
use std::collections::HashSet;
use std::ops::Range;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DigestSlice {
    ref_seq: Arc<str>,
    range: Range<usize>,
    pub decoy: DecoyMarking,
}

impl Serialize for DigestSlice {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let local_str = Into::<String>::into(self.clone());
        serializer.serialize_str(local_str.as_str())
    }
}

impl DigestSlice {
    pub fn new(ref_seq: Arc<str>, range: Range<usize>, decoy: DecoyMarking) -> Self {
        Self {
            ref_seq,
            range,
            decoy,
        }
    }

    pub fn from_string(seq: String, decoy: bool) -> DigestSlice {
        let len = seq.len();
        DigestSlice {
            ref_seq: seq.into(),
            range: 0..len,
            decoy: if decoy {
                DecoyMarking::ReversedDecoy
            } else {
                DecoyMarking::Target
            },
        }
    }

    pub fn as_decoy(&self) -> DigestSlice {
        DigestSlice {
            ref_seq: self.ref_seq.clone(),
            range: self.range.clone(),
            decoy: DecoyMarking::Decoy,
        }
    }

    pub fn as_decoy_string(&self) -> String {
        as_decoy_string(&self.ref_seq.as_ref()[self.range.clone()])
    }

    pub fn len(&self) -> usize {
        self.range.len()
    }

    pub fn is_empty(&self) -> bool {
        self.range.is_empty()
    }
}

impl From<DigestSlice> for String {
    fn from(x: DigestSlice) -> Self {
        let tmp = &x.ref_seq.as_ref()[x.range.clone()];

        match x.decoy {
            DecoyMarking::Target => tmp.to_string(),
            DecoyMarking::ReversedDecoy => tmp.to_string(),
            DecoyMarking::Decoy => as_decoy_string(tmp),
        }
    }
}

pub fn deduplicate_digests(mut digest_slices: Vec<DigestSlice>) -> Vec<DigestSlice> {
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
        let digests: Vec<DigestSlice> = vec![
            DigestSlice {
                ref_seq: seq.clone(),
                range: 0..seq.as_ref().len(),
                decoy: DecoyMarking::Target,
            },
            DigestSlice {
                ref_seq: seq.clone(),
                range: 0..seq2.as_ref().len(), // Note the short length
                decoy: DecoyMarking::Target,
            },
            DigestSlice {
                ref_seq: seq2.clone(),
                range: 0..seq2.as_ref().len(),
                decoy: DecoyMarking::Target,
            },
            DigestSlice {
                ref_seq: seq2_rep.clone(),
                range: 0..seq2_rep.as_ref().len(),
                decoy: DecoyMarking::Target,
            },
        ];
        let deduped = deduplicate_digests(digests);
        assert_eq!(deduped.len(), 2);
        assert_eq!(deduped[0].len(), seq.as_ref().len());
        assert_eq!(deduped[1].len(), seq2.as_ref().len());
    }

    #[test]
    fn test_from_string() {
        let seq = "PEPTIDEPINKTOMATO".to_string();
        let expect_len = seq.len();
        let digests = DigestSlice::from_string(seq.clone(), false);
        let digests_decoy = DigestSlice::from_string(seq.clone(), true);
        assert_eq!(digests.len(), expect_len);
        assert_eq!(digests.decoy, DecoyMarking::Target);
        assert_eq!(digests_decoy.len(), expect_len);
        assert_eq!(digests_decoy.decoy, DecoyMarking::ReversedDecoy);

        // Check that converting back to string does not reverse the sequence
        assert_eq!(Into::<String>::into(digests.clone()), seq);
        assert_eq!(Into::<String>::into(digests_decoy.clone()), seq);
    }
}
