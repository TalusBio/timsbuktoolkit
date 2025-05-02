use crate::errors::TimsSeekError;
use crate::fragment_mass::IonAnnot;
use crate::models::{
    DecoyMarking,
    DigestSlice,
};
use crate::scoring::QueryItemToScore;
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};
use std::collections::HashMap;
use std::io::{
    BufRead,
    BufReader,
};
use std::path;
use std::sync::Arc;
use timsquery::models::elution_group::ElutionGroup;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SerSpeclibElement {
    precursor: PrecursorEntry,
    elution_group: ReferenceEG,
}

#[derive(Debug, Clone)]
struct SpeclibElement {
    precursor: DigestSlice,
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
        let range = 0 as u16..seq.as_ref().len() as u16;
        DigestSlice::new(seq, range, decoy)
    }
}

impl From<SerSpeclibElement> for SpeclibElement {
    fn from(x: SerSpeclibElement) -> Self {
        let precursor = x.precursor.into();
        let elution_group = x.elution_group;
        SpeclibElement {
            precursor,
            elution_group,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedIntensities {
    pub fragment_intensities: HashMap<IonAnnot, f32>,
    pub precursor_intensities: Vec<f32>,
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
    elems: Vec<SpeclibElement>,
}

impl Speclib {
    pub fn from_json(json: &str) -> Self {
        let speclib_ser: Vec<SerSpeclibElement> = serde_json::from_str(json).unwrap();

        let speclib = speclib_ser.into_iter().map(|x| {
            let precursor = x.precursor.into();
            let elution_group = x.elution_group;
            SpeclibElement {
                precursor,
                elution_group,
            }
        }).collect();

        Self {
            elems: speclib
        }
    }

    pub fn from_ndjson(json: &str) -> Self {
        // Split on newlines and parse each ...
        let lines: Vec<&str> = json.split('\n').collect();
        Self::from_ndjson_elems(&lines)
    }

    fn from_ndjson_elems(lines: &[&str]) -> Self {
        let tmp= lines
            .into_par_iter()
            .filter(|line| !line.is_empty())
            .map(|line| {
                let elem: SerSpeclibElement = match serde_json::from_str(line) {
                    Ok(x) => x,
                    Err(e) => {
                        panic!("Error parsing line: {:?}; error: {:?}", line, e);
                        // TODO: Make this a proper error
                        // return Err(TimsSeekError::TimsRust(TimsRustError::Serde(e)));
                    }
                };

                elem.into()
            }).collect::<Vec<SpeclibElement>>();

        Self { elems: tmp }
    }

    pub fn from_ndjson_file(path: &path::Path) -> Result<Self, TimsSeekError> {
        let file = std::fs::File::open(path).map_err(|e| TimsSeekError::Io {
            source: e,
            path: Some(path.to_path_buf()),
        })?;

        let mut reader = BufReader::new(file);
        let mut current_chunk = String::new();
        let mut curr_nlines = 0;
        let mut results = Self {
            elems: Vec::with_capacity(100_001),
        };

        while let Ok(len) = reader.read_line(&mut current_chunk) {
            curr_nlines += 1;
            if len == 0 {
                break;
            }

            if curr_nlines % 100_000 == 0 {
                let chunk_result = Self::from_ndjson(&current_chunk);
                results = results.fold(chunk_result);
                current_chunk.clear();
                curr_nlines = 0;
            }
        }

        // Process remaining lines
        if !current_chunk.is_empty() {
            let chunk_result = Self::from_ndjson(&current_chunk);
            results = results.fold(chunk_result);
        }

        Ok(results)
    }

    fn fold(self, other: Self) -> Self {
        Self {
            elems: self.elems.into_iter().chain(other.elems).collect(),
        }
    }

    pub fn len(&self) -> usize {
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
        let speclib = Speclib::from_json(json);
        assert_eq!(speclib.len(), 1);
        assert_eq!(speclib.elems.len(), 1);
        println!("{:?}", speclib);

        assert_eq!(speclib.elems[0].precursor.decoy, DecoyMarking::Target);
        assert_eq!(speclib.elems[0].precursor.len(), "PEPTIDEPINK".len());
        assert_eq!(speclib.elems[0].elution_group.elution_group.fragments.len(), 3);
    }
}
