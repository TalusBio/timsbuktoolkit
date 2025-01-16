use crate::errors::TimsSeekError;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::models::{
    DecoyMarking,
    DigestSlice,
    NamedQueryChunk,
};
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
pub struct ExpectedIntensities {
    pub fragment_intensities: HashMap<SafePosition, f32>,
    pub precursor_intensities: Vec<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceEG {
    #[serde(flatten)]
    pub elution_group: ElutionGroup<SafePosition>,
    #[serde(flatten)]
    pub expected_intensities: ExpectedIntensities,
}

#[derive(Debug, Clone)]
pub struct Speclib {
    digests: Vec<DigestSlice>,
    charges: Vec<u8>,
    queries: Vec<ElutionGroup<SafePosition>>,
    expected_intensities: Vec<ExpectedIntensities>,
}

pub struct SpeclibIterator {
    speclib: Speclib,
    chunk_size: usize,
    max_iterations: usize,
    iteration_index: usize,
}

impl SpeclibIterator {
    pub fn new(speclib: Speclib, chunk_size: usize) -> Self {
        let max_iters = speclib.digests.len() / chunk_size;
        Self {
            speclib,
            chunk_size,
            max_iterations: max_iters,
            iteration_index: 0,
        }
    }
}

impl Iterator for SpeclibIterator {
    type Item = NamedQueryChunk;

    fn next(&mut self) -> Option<Self::Item> {
        // No need to make decoys when we have a speclib!!
        let out = self
            .speclib
            .get_chunk(self.iteration_index, self.chunk_size);
        self.iteration_index += 1;
        out
    }
}

impl ExactSizeIterator for SpeclibIterator {
    fn len(&self) -> usize {
        self.max_iterations
    }
}

impl Speclib {
    pub fn from_json(json: &str) -> Self {
        let speclib: Vec<SpeclibElement> = serde_json::from_str(json).unwrap();

        let (exp_int, (queries, (charges, digests))): (
            Vec<ExpectedIntensities>,
            (Vec<ElutionGroup<SafePosition>>, (Vec<u8>, Vec<DigestSlice>)),
        ) = speclib
            .into_par_iter()
            .map(|x| {
                let charge = x.precursor.charge;
                let elution_group = x.elution_group.elution_group;
                let expected_intensities = x.elution_group.expected_intensities;
                let digest = x.precursor.into();
                (expected_intensities, (elution_group, (charge, digest)))
            })
            .unzip();

        Self {
            digests,
            charges,
            queries,
            expected_intensities: exp_int,
        }
    }

    pub fn from_ndjson(json: &str) -> Self {
        // Split on newlines and parse each ...
        let lines: Vec<&str> = json.split('\n').collect();
        Self::from_ndjson_elems(&lines)
    }

    fn from_ndjson_elems(lines: &[&str]) -> Self {
        let ((charges, digests), (queries, expected_intensities)): (
            (Vec<_>, Vec<_>),
            (Vec<_>, Vec<_>),
        ) = lines
            .into_par_iter()
            .filter(|line| !line.is_empty())
            .map(|line| {
                let elem: SpeclibElement = match serde_json::from_str(line) {
                    Ok(x) => x,
                    Err(e) => {
                        panic!("Error parsing line: {:?}; error: {:?}", line, e);
                        // TODO: Make this a proper error
                        // return Err(TimsSeekError::TimsRust(TimsRustError::Serde(e)));
                    }
                };

                (
                    (elem.precursor.charge, elem.precursor.into()),
                    (
                        elem.elution_group.elution_group,
                        elem.elution_group.expected_intensities,
                    ),
                )
            })
            .unzip();

        if digests.is_empty() {
            panic!("No digests found in speclib file");
        }

        Self {
            digests,
            charges,
            queries,
            expected_intensities,
        }
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
            digests: Vec::with_capacity(100_001),
            charges: Vec::with_capacity(100_001),
            queries: Vec::with_capacity(100_001),
            expected_intensities: Vec::with_capacity(100_001),
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

    fn get_chunk(&self, chunk_index: usize, chunk_size: usize) -> Option<NamedQueryChunk> {
        let start = chunk_index * chunk_size;
        if start >= self.digests.len() {
            return None;
        }
        let end = start + chunk_size;
        let end = if end > self.digests.len() {
            self.digests.len()
        } else {
            end
        };
        let digests = &self.digests[start..end];
        let charges = &self.charges[start..end];
        let queries = &self.queries[start..end];
        let expected_intensities = &self.expected_intensities[start..end];
        Some(NamedQueryChunk::new(
            digests.to_vec(),
            charges.to_vec(),
            queries.to_vec(),
            expected_intensities.to_vec(),
        ))
    }

    fn fold(self, other: Self) -> Self {
        let mut digests = self.digests;
        digests.extend(other.digests);
        let mut charges = self.charges;
        charges.extend(other.charges);
        let mut queries = self.queries;
        queries.extend(other.queries);
        let mut expected_intensities = self.expected_intensities;
        expected_intensities.extend(other.expected_intensities);
        Self {
            digests,
            charges,
            queries,
            expected_intensities,
        }
    }

    pub fn as_iterator(self, chunk_size: usize) -> SpeclibIterator {
        // TODO make an iterable version of this where the file can be "read" lazily
        SpeclibIterator::new(self, chunk_size)
    }

    pub fn len(&self) -> usize {
        self.digests.len()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpeclibElement {
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
        let range = 0..seq.as_ref().len();
        DigestSlice::new(seq, range, decoy)
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
        assert_eq!(speclib.digests.len(), 1);
        assert_eq!(speclib.charges.len(), 1);
        assert_eq!(speclib.queries.len(), 1);
        println!("{:?}", speclib);

        assert_eq!(speclib.digests[0].decoy, DecoyMarking::Target);
        assert_eq!(speclib.digests[0].len(), "PEPTIDEPINK".len());
        assert_eq!(speclib.queries[0].fragment_mzs.len(), 3);
    }
}
