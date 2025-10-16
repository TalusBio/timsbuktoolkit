use crate::errors::LibraryReadingError;
use crate::models::{
    DecoyMarking,
    DigestSlice,
};
use crate::{
    ExpectedIntensities,
    IonAnnot,
    QueryItemToScore,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::io::{
    BufRead,
    BufReader,
    Read,
};
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;
use timsquery::models::elution_group::ElutionGroup;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerSpeclibElement {
    precursor: PrecursorEntry,
    elution_group: ReferenceEG,
}

impl SerSpeclibElement {
    pub fn sample() -> Self {
        QueryItemToScore::sample().into()
    }
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

#[derive(Debug, Clone, Copy)]
pub enum SpeclibFormat {
    NdJson,
    NdJsonZstd,
    MessagePack,
    MessagePackZstd,
}

impl SpeclibFormat {
    pub fn detect_from_path(path: &Path) -> Result<Self, LibraryReadingError> {
        let path_str = path.to_string_lossy().to_lowercase();

        if path_str.ends_with(".msgpack.zst") {
            Ok(SpeclibFormat::MessagePackZstd)
        } else if path_str.ends_with(".msgpack") {
            Ok(SpeclibFormat::MessagePack)
        } else if path_str.ends_with(".ndjson.zst") {
            Ok(SpeclibFormat::NdJsonZstd)
        } else if path_str.ends_with(".ndjson") {
            Ok(SpeclibFormat::NdJson)
        } else {
            // Try to detect by reading first few bytes
            Self::detect_from_content(path)
        }
    }

    fn detect_from_content(path: &Path) -> Result<Self, LibraryReadingError> {
        let file =
            std::fs::File::open(path).map_err(|e| LibraryReadingError::FileReadingError {
                source: e,
                context: "Error opening file for format detection",
                path: PathBuf::from(path),
            })?;

        let mut reader = BufReader::new(file);
        let mut buffer = [0u8; 8];

        match reader.read(&mut buffer) {
            Ok(bytes_read) if bytes_read >= 4 => {
                if buffer[0..4] == [0x28, 0xB5, 0x2F, 0xFD] {
                    Ok(SpeclibFormat::MessagePackZstd)
                } else if buffer[0] == b'{' {
                    Ok(SpeclibFormat::NdJson)
                } else {
                    Ok(SpeclibFormat::MessagePack)
                }
            }
            _ => Ok(SpeclibFormat::NdJson),
        }
    }
}

pub struct SpeclibReader<'a> {
    inner: Box<dyn Iterator<Item = Result<QueryItemToScore, LibraryReadingError>> + Send + 'a>,
}

impl<'a> SpeclibReader<'a> {
    pub fn new<R: Read + Send + 'a>(
        reader: R,
        format: SpeclibFormat,
    ) -> Result<Self, LibraryReadingError> {
        let inner: Box<dyn Iterator<Item = Result<QueryItemToScore, LibraryReadingError>> + Send> =
            match format {
                SpeclibFormat::NdJson => Box::new(NdJsonReader::new(BufReader::new(reader))),
                SpeclibFormat::NdJsonZstd => {
                    let decoder = zstd::Decoder::new(reader).map_err(|e| {
                        LibraryReadingError::SpeclibParsingError {
                            source: serde_json::Error::io(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                e,
                            )),
                            context: "Error creating ZSTD decoder",
                        }
                    })?;
                    Box::new(NdJsonReader::new(BufReader::new(decoder)))
                }
                SpeclibFormat::MessagePack => Box::new(MessagePackReader::new(reader)),
                SpeclibFormat::MessagePackZstd => {
                    let decoder = zstd::Decoder::new(reader).map_err(|e| {
                        LibraryReadingError::SpeclibParsingError {
                            source: serde_json::Error::io(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                e,
                            )),
                            context: "Error creating ZSTD decoder",
                        }
                    })?;
                    Box::new(MessagePackReader::new(decoder))
                }
            };

        Ok(SpeclibReader { inner })
    }
}

impl Iterator for SpeclibReader<'_> {
    type Item = Result<QueryItemToScore, LibraryReadingError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

struct NdJsonReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> NdJsonReader<R> {
    fn new(reader: R) -> Self {
        Self { reader }
    }
}

impl<R: BufRead> Iterator for NdJsonReader<R> {
    type Item = Result<QueryItemToScore, LibraryReadingError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => None, // EOF
            Ok(_) => {
                if line.trim().is_empty() {
                    return self.next(); // Skip empty lines
                }

                let elem: SerSpeclibElement = match serde_json::from_str(&line) {
                    Ok(x) => x,
                    Err(e) => {
                        return Some(Err(LibraryReadingError::SpeclibParsingError {
                            source: e,
                            context: "Error parsing NDJSON line",
                        }));
                    }
                };

                Some(Ok(elem.into()))
            }
            Err(e) => Some(Err(LibraryReadingError::FileReadingError {
                source: e,
                context: "Error reading line",
                path: PathBuf::new(),
            })),
        }
    }
}

struct MessagePackReader<R: Read> {
    deserializer: rmp_serde::Deserializer<rmp_serde::decode::ReadReader<R>>,
}

impl<R: Read> MessagePackReader<R> {
    fn new(reader: R) -> Self {
        Self {
            deserializer: rmp_serde::Deserializer::new(reader),
        }
    }
}

impl<R: Read> Iterator for MessagePackReader<R> {
    type Item = Result<QueryItemToScore, LibraryReadingError>;

    fn next(&mut self) -> Option<Self::Item> {
        use serde::Deserialize;

        match SerSpeclibElement::deserialize(&mut self.deserializer) {
            Ok(elem) => Some(Ok(elem.into())),
            Err(rmp_serde::decode::Error::InvalidMarkerRead(ref io_err))
                if io_err.kind() == std::io::ErrorKind::UnexpectedEof =>
            {
                None
            } // EOF
            Err(rmp_serde::decode::Error::InvalidDataRead(ref io_err))
                if io_err.kind() == std::io::ErrorKind::UnexpectedEof =>
            {
                None
            } // EOF
            Err(e) => Some(Err(LibraryReadingError::SpeclibParsingError {
                source: serde_json::Error::io(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    e,
                )),
                context: "Error reading MessagePack",
            })),
        }
    }
}

impl Speclib {
    pub fn from_file(path: &Path) -> Result<Self, LibraryReadingError> {
        let format = SpeclibFormat::detect_from_path(path)?;
        Self::from_file_with_format(path, format)
    }

    pub fn from_file_with_format(
        path: &Path,
        format: SpeclibFormat,
    ) -> Result<Self, LibraryReadingError> {
        let file =
            std::fs::File::open(path).map_err(|e| LibraryReadingError::FileReadingError {
                source: e,
                context: "Error opening speclib file",
                path: PathBuf::from(path),
            })?;

        let reader = SpeclibReader::new(file, format)?;
        let elements: Result<Vec<_>, _> = reader.collect();

        Ok(Self {
            elems: elements?,
            idx: 0,
        })
    }

    // pub(crate) fn from_ndjson(json: &str) -> Result<Self, LibraryReadingError> {
    //     // Split on newlines and parse each ...
    //     // In the future I want to make this lazy but this will do for now
    //     let lines: Vec<&str> = json.split('\n').collect();
    //     Self::from_ndjson_elems(&lines)
    // }

    // fn from_ndjson_elems(lines: &[&str]) -> Result<Self, LibraryReadingError> {
    //     let tmp = lines
    //         .into_par_iter()
    //         .filter(|line| !line.is_empty())
    //         .map(|line| {
    //             let elem: SerSpeclibElement = match serde_json::from_str(line) {
    //                 Ok(x) => x,
    //                 Err(e) => {
    //                     return Err(LibraryReadingError::SpeclibParsingError {
    //                         source: e,
    //                         context: "Error parsing line in NDJSON",
    //                     });
    //                 }
    //             };

    //             Ok(elem.into())
    //         })
    //         .collect::<Result<Vec<_>, LibraryReadingError>>()?;

    //     Ok(Self { elems: tmp, idx: 0 })
    // }

    // pub(crate) fn from_ndjson_file(path: &Path) -> Result<Self, LibraryReadingError> {
    //     let file =
    //         std::fs::File::open(path).map_err(|e| LibraryReadingError::FileReadingError {
    //             source: e,
    //             context: "Error opening NDJSON file",
    //             path: PathBuf::from(path),
    //         })?;

    //     let mut reader = BufReader::new(file);
    //     let mut current_chunk = String::new();
    //     let mut curr_nlines = 0;
    //     let mut results = Self {
    //         elems: Vec::with_capacity(100_001),
    //         idx: 0,
    //     };

    //     while let Ok(len) = reader.read_line(&mut current_chunk) {
    //         curr_nlines += 1;
    //         if len == 0 {
    //             break;
    //         }

    //         if curr_nlines % 100_000 == 0 {
    //             let chunk_result = Self::from_ndjson(&current_chunk)?;
    //             results = results.fold(chunk_result);
    //             current_chunk.clear();
    //             curr_nlines = 0;
    //         }
    //     }

    //     // Process remaining lines
    //     if !current_chunk.is_empty() {
    //         let chunk_result = Self::from_ndjson(&current_chunk)?;
    //         results = results.fold(chunk_result);
    //     }

    //     Ok(results)
    // }

    pub fn sample() -> Self {
        Self {
            elems: vec![QueryItemToScore::sample()],
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

    impl Speclib {
        // Technically its used in testing ...
        fn from_json(json: &str) -> Result<Self, LibraryReadingError> {
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
    }

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
                    "precursors": [
                        [0, 1810.917339999999],
                        [1, 1810.917339999999]
                    ],
                    "fragments": [
                        ["a1", 123.0],
                        ["b1", 123.0],
                        ["c1^2", 123.0]
                    ],
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
