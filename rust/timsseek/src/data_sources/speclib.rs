use crate::errors::{
    LibraryReadingError,
    TimsSeekError,
};
use crate::models::{
    DecoyMarking,
    DigestSlice,
};
use crate::{
    ExpectedIntensities,
    IonAnnot,
    QueryItemToScore,
};
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};
use std::io::{
    BufRead,
    BufReader,
};
use std::ops::Deref;
use std::path::{
    self,
    Path,
    PathBuf,
};
use std::sync::Arc;
use timsquery::models::elution_group::ElutionGroup;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SerSpeclibElement {
    precursor: PrecursorEntry,
    elution_group: ReferenceEG,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PrecursorEntry {
    sequence: String,
    charge: u8,
    decoy: bool,
}

impl From<PrecursorEntry> for DigestSlice {
    fn from(x: PrecursorEntry) -> Self {
        let decoy = if x.decoy {
            DecoyMarking::ReversedDecoy
        } else {
            DecoyMarking::Target
        };
        let seq: Arc<str> = x.sequence.clone().into();
        if seq.as_ref().len() >= u16::MAX as usize {
            panic!("Sequence too long (gt: {}): {}", u16::MAX, seq);
        }
        let range = 0u16..seq.as_ref().len() as u16;
        DigestSlice::new(seq, range, decoy)
    }
}

impl From<SerSpeclibElement> for QueryItemToScore {
    fn from(x: SerSpeclibElement) -> Self {
        let charge = x.precursor.charge;
        let precursor = x.precursor.into();
        let elution_group = x.elution_group;
        QueryItemToScore {
            expected_intensity: elution_group.expected_intensities,
            charge,
            query: elution_group.elution_group,
            digest: precursor,
        }
    }
}

impl From<QueryItemToScore> for SerSpeclibElement {
    fn from(x: QueryItemToScore) -> Self {
        let precursor = PrecursorEntry {
            sequence: x.digest.clone().into(),
            charge: x.charge,
            decoy: x.digest.is_decoy(),
        };
        let elution_group = ReferenceEG {
            elution_group: x.query.clone(),
            expected_intensities: x.expected_intensity,
        };
        SerSpeclibElement {
            precursor,
            elution_group,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceEG {
    #[serde(flatten)]
    pub elution_group: Arc<ElutionGroup<IonAnnot>>,
    #[serde(flatten)]
    pub expected_intensities: ExpectedIntensities,
}

#[derive(Debug, Clone)]
pub struct Speclib {
    elems: Vec<QueryItemToScore>,
    idx: usize,
}

impl Serialize for Speclib {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let speclib_ser: Vec<SerSpeclibElement> =
            self.elems.iter().map(|x| x.clone().into()).collect();
        speclib_ser.serialize(serializer)
    }
}

impl Speclib {
    pub fn from_json(json: &str) -> Result<Self, LibraryReadingError> {
        let speclib_ser: Vec<SerSpeclibElement> = match serde_json::from_str(json) {
            Ok(x) => x,
            Err(e) => {
                return Err(LibraryReadingError::SpeclibParsingError {
                    source: e,
                    context: "Error parsing JSON",
                });
            }
        };

        let speclib = speclib_ser.into_iter().map(|x| x.into()).collect();

        Ok(Self {
            elems: speclib,
            idx: 0,
        })
    }

    pub fn from_ndjson(json: &str) -> Result<Self, LibraryReadingError> {
        // Split on newlines and parse each ...
        // In the future I want to make this lazy but this will do for now
        let lines: Vec<&str> = json.split('\n').collect();
        Self::from_ndjson_elems(&lines)
    }

    fn from_ndjson_elems(lines: &[&str]) -> Result<Self, LibraryReadingError> {
        let tmp = lines
            .into_par_iter()
            .filter(|line| !line.is_empty())
            .map(|line| {
                let elem: SerSpeclibElement = match serde_json::from_str(line) {
                    Ok(x) => x,
                    Err(e) => {
                        return Err(LibraryReadingError::SpeclibParsingError {
                            source: e,
                            context: "Error parsing line in NDJSON",
                        });
                    }
                };

                Ok(elem.into())
            })
            .collect::<Result<Vec<_>, LibraryReadingError>>()?;

        Ok(Self { elems: tmp, idx: 0 })
    }

    pub fn from_ndjson_file(path: &path::Path) -> Result<Self, LibraryReadingError> {
        let file =
            std::fs::File::open(path).map_err(|e| LibraryReadingError::FileReadingError {
                source: e,
                context: "Error opening NDJSON file",
                path: PathBuf::from(path),
            })?;

        let mut reader = BufReader::new(file);
        let mut current_chunk = String::new();
        let mut curr_nlines = 0;
        let mut results = Self {
            elems: Vec::with_capacity(100_001),
            idx: 0,
        };

        while let Ok(len) = reader.read_line(&mut current_chunk) {
            curr_nlines += 1;
            if len == 0 {
                break;
            }

            if curr_nlines % 100_000 == 0 {
                let chunk_result = Self::from_ndjson(&current_chunk)?;
                results = results.fold(chunk_result);
                current_chunk.clear();
                curr_nlines = 0;
            }
        }

        // Process remaining lines
        if !current_chunk.is_empty() {
            let chunk_result = Self::from_ndjson(&current_chunk)?;
            results = results.fold(chunk_result);
        }

        Ok(results)
    }

    pub fn sample() -> Self {
        Self {
            elems: vec![QueryItemToScore::sample()],
            idx: 0,
        }
    }

    fn fold(self, other: Self) -> Self {
        assert_eq!(self.idx, 0);
        assert_eq!(other.idx, 0);
        Self {
            elems: self.elems.into_iter().chain(other.elems).collect(),
            idx: 0,
        }
    }

    pub fn as_slice(&self) -> &[QueryItemToScore] {
        &self.elems
    }
}

impl Iterator for Speclib {
    type Item = QueryItemToScore;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.elems.len() {
            return None;
        }
        // TODO: reimplement as a real stream ...
        // this is probably not the bottleneck but would be nice
        // to check how expensive the clone actually is.
        let elem = self.elems[self.idx].clone();
        self.idx += 1;
        Some(elem)
    }
}

impl ExactSizeIterator for Speclib {
    fn len(&self) -> usize {
        self.elems.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_speclib() {
        let json = r#"[
            {
                "precursor": {
                    "sequence": "PEPTIDEPINK",
                    "charge": 2,
                    "decoy": false
                },
                "elution_group": {
                    "id": 0,
                    "precursor_mzs": [
                        1810.917339999999,
                        1810.917339999999
                    ],
                    "fragment_mzs": {
                        "a1": 123.0,
                        "b1": 123.0,
                        "c1^2": 123.0
                    },
                    "precursor_charge": 2,
                    "mobility": 0.8,
                    "rt_seconds": 0.0,
                    "precursor_intensities": [
                        1.0,
                        1.0
                    ],
                    "fragment_intensities": {
                        "a1": 1.0,
                        "b1": 1.0,
                        "c1^2": 1.0
                    }
                }
            }
        ]"#;
        let speclib = Speclib::from_json(json).unwrap();
        assert_eq!(speclib.len(), 1);
        assert_eq!(speclib.elems.len(), 1);
        println!("{:?}", speclib);

        assert_eq!(speclib.elems[0].digest.decoy, DecoyMarking::Target);
        assert_eq!(speclib.elems[0].digest.len(), "PEPTIDEPINK".len());
        assert_eq!(speclib.elems[0].query.fragments.len(), 3);
    }
}
