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
use timsquery::models::elution_group::TimsElutionGroup;

/// This is meant to the be the serializable version of the speclib element
/// so ... in general should be backwards compatible and implement `Into<QueryItemToScore>`
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerSpeclibElement {
    precursor: PrecursorEntry,
    elution_group: ReferenceEG,
}

impl SerSpeclibElement {
    pub fn sample() -> Self {
        // let eg = TimsElutionGroup::builder()
        //     .id(42)
        //     .mobility_ook0(0.75)
        //     .rt_seconds(123.4)
        //     .fragment_labels([
        //         IonAnnot::try_from("y1").unwrap(),
        //         IonAnnot::try_from("y2").unwrap(),
        //         IonAnnot::try_from("y3").unwrap(),
        //         IonAnnot::try_from("y4").unwrap(),
        //         ].as_slice().into())
        //     .fragment_mzs(
        //         vec![450.0,
        //             650.5,
        //             751.0,
        //             751.5,]
        //     )
        //     .precursor_labels(tiny_vec!(-1,0,1,2))
        //     .precursor_mzs(vec![450.0, 450.5, 451.0, 451.5]).try_build().unwrap();

        // let ei = ExpectedIntensities {
        //     fragment_intensities: HashMap::from_iter(
        //         [
        //             (IonAnnot::try_from("y1").unwrap(), 1.0),
        //             (IonAnnot::try_from("y2").unwrap(), 1.0),
        //             (IonAnnot::try_from("y3").unwrap(), 1.0),
        //             (IonAnnot::try_from("y4").unwrap(), 1.0),
        //         ]
        //         .iter()
        //         .cloned(),
        //     ),
        //     precursor_intensities: HashMap::from_iter(
        //         [(-1, 0.5), (0, 1.0), (1, 0.8), (2, 0.3)].iter().cloned(),
        //     ),
        // };
        // let pepseq = "PEPTIDEPINKPEPTIDE".into();
        // let digest = DigestSlice::from_string(pepseq, false, 1);
        // let charge = 2;
        // let query = eg;
        // let expected_intensity = ei;
        // QueryItemToScore {
        //     digest,
        //     charge,
        //     query,
        //     expected_intensity,
        // }
        SerSpeclibElement {
            precursor: PrecursorEntry {
                sequence: "PEPTIDESEK".into(),
                charge: 2,
                decoy: false,
                decoy_group: 32,
            },
            elution_group: ReferenceEG {
                id: 32,
                precursor_mz: 512.2,
                precursor_labels: vec![0, 2],
                fragment_mzs: vec![312.2, 675.7],
                fragment_labels: vec![
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                ],
                precursor_intensities: vec![1.0, 0.5],
                fragment_intensities: vec![0.8, 0.3],
                mobility_ook0: 0.75,
                rt_seconds: 120.0,
            },
        }
    }

    pub fn sample_json() -> &'static str {
        r#"{
            "precursor": {
                "sequence": "PEPTIDEPINK",
                "charge": 2,
                "decoy": false,
                "decoy_group": 0
            },
            "elution_group": {
                "id": 0,
                "precursor_mz": 876.5432,
                "precursor_labels": [ 0, 1 ],
                "fragment_mzs": [ 123.0, 123.0, 123.0 ],
                "fragment_labels": ["a1", "b1", "c1^2"],
                "precursor_intensities": [1.0, 1.0],
                "fragment_intensities": [1.0, 1.0, 1.0],
                "precursor_charge": 2,
                "mobility_ook0": 0.8,
                "rt_seconds": 0.0
            }
        }"#
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PrecursorEntry {
    sequence: String,
    charge: u8,
    decoy: bool,
    decoy_group: u32,
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
        DigestSlice::new(seq, range, decoy, x.decoy_group)
    }
}

impl From<SerSpeclibElement> for QueryItemToScore {
    fn from(x: SerSpeclibElement) -> Self {
        let charge = x.precursor.charge;
        let precursor = x.precursor.into();
        let ref_eg = x.elution_group;
        QueryItemToScore {
            expected_intensity: ExpectedIntensities {
                fragment_intensities: ref_eg
                    .fragment_labels
                    .iter()
                    .cloned()
                    .zip(ref_eg.fragment_intensities.iter().cloned())
                    .collect(),
                precursor_intensities: ref_eg
                    .precursor_labels
                    .iter()
                    .cloned()
                    .zip(ref_eg.precursor_intensities.iter().cloned())
                    .collect(),
            },
            query: TimsElutionGroup::builder()
                .id(ref_eg.id as u64)
                .mobility_ook0(ref_eg.mobility_ook0)
                .rt_seconds(ref_eg.rt_seconds)
                .precursor_labels(ref_eg.precursor_labels.as_slice().into())
                .precursor(ref_eg.precursor_mz, charge)
                .fragment_labels(ref_eg.fragment_labels.as_slice().into())
                .fragment_mzs(ref_eg.fragment_mzs)
                .try_build()
                .unwrap(),
            digest: precursor,
        }
    }
}

// impl From<QueryItemToScore> for SerSpeclibElement {
//     fn from(x: QueryItemToScore) -> Self {
//         let precursor = PrecursorEntry {
//             sequence: x.digest.clone().into(),
//             charge: x.charge,
//             decoy: x.digest.is_decoy(),
//             decoy_group: x.digest.decoy_group,
//         };
//         let mut precursor_mzs = Vec::with_capacity(x.query.precursors.len());
//         let mut precursor_labels = Vec::with_capacity(x.query.precursors.len());
//         for (label, mz) in x.query.precursors.iter() {
//             precursor_labels.push(*label);
//             precursor_mzs.push(*mz);
//         }
//         let mut fragment_mzs = Vec::with_capacity(x.query.fragments.len());
//         let mut fragment_labels = Vec::with_capacity(x.query.fragments.len());
//         for (label, mz) in x.query.fragments.iter() {
//             fragment_labels.push(*label);
//             fragment_mzs.push(*mz);
//         }
//         let precursor_intensities: Vec<f32> = precursor_labels
//             .iter()
//             .map(|l| {
//                 *x.expected_intensity
//                     .precursor_intensities
//                     .get(l)
//                     .expect("Correctly built reference eg should have all keys")
//             })
//             .collect();
//         let fragment_intensities: Vec<f32> = fragment_labels
//             .iter()
//             .map(|l| {
//                 *x.expected_intensity
//                     .fragment_intensities
//                     .get(l)
//                     .expect("Correctly built reference eg should have all keys")
//             })
//             .collect();
//         let elution_group = ReferenceEG {
//             id: x.query.id() as u32,
//             precursor_mzs,
//             precursor_labels,
//             fragment_mzs,
//             fragment_labels,
//             precursor_intensities,
//             fragment_intensities,
//             precursor_charge: x.charge,
//             mobility_ook0: x.query.mobility_ook0(),
//             rt_seconds: x.query.rt_seconds(),
//         };
//         SerSpeclibElement {
//             precursor,
//             elution_group,
//         }
//     }
// }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceEG {
    id: u32,
    precursor_mz: f64,
    precursor_labels: Vec<i8>,
    #[serde(alias = "fragment_mz")]
    fragment_mzs: Vec<f64>,
    fragment_labels: Vec<IonAnnot>,
    precursor_intensities: Vec<f32>,
    fragment_intensities: Vec<f32>,
    #[serde(alias = "mobility")]
    mobility_ook0: f32,
    rt_seconds: f32,
}

#[derive(Debug, Clone)]
pub struct Speclib {
    pub elems: Vec<QueryItemToScore>,
}

struct SpeclibIterator<'a> {
    speclib: &'a Speclib,
    idx: usize,
}

impl Iterator for SpeclibIterator<'_> {
    type Item = QueryItemToScore;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.speclib.elems.len() {
            return None;
        }
        let elem = self.speclib.elems[self.idx].clone();
        self.idx += 1;
        Some(elem)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.speclib.elems.len() - self.idx;
        (remaining, Some(remaining))
    }
}

// impl Serialize for Speclib {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: serde::Serializer,
//     {
//         let speclib_ser: Vec<SerSpeclibElement> =
//             self.elems.iter().map(|x| x.clone().into()).collect();
//         speclib_ser.serialize(serializer)
//     }
// }

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
    pub fn iter<'a>(&'a self) -> impl Iterator<Item = QueryItemToScore> + 'a {
        SpeclibIterator {
            speclib: self,
            idx: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.elems.len()
    }

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
        let elems = elements?;

        Ok(Self { elems })
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
        }
    }

    pub fn as_slice(&self) -> &[QueryItemToScore] {
        &self.elems
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

            Ok(Self { elems: speclib })
        }
    }

    #[test]
    fn test_speclib() {
        let json = r#"[
            {
                "precursor": {
                    "sequence": "PEPTIDEPINK",
                    "charge": 2,
                    "decoy": false,
                    "decoy_group": 0
                },
                "elution_group": {
                    "id": 0,
                    "precursor_mz": 876.5432,
                    "precursor_labels": [ 0, 1 ],
                    "fragment_mzs": [ 123.0, 123.0, 123.0 ],
                    "fragment_labels": ["a1", "b1", "c1^2"],
                    "precursor_intensities": [1.0, 1.0],
                    "fragment_intensities": [1.0, 1.0, 1.0],
                    "precursor_charge": 2,
                    "mobility_ook0": 0.8,
                    "rt_seconds": 0.0
                }
            }
        ]"#;
        let speclib = Speclib::from_json(json).unwrap();
        assert_eq!(speclib.elems.len(), 1);
        println!("{:?}", speclib);

        assert_eq!(speclib.elems[0].digest.decoy, DecoyMarking::Target);
        assert_eq!(speclib.elems[0].digest.len(), "PEPTIDEPINK".len());
        assert_eq!(speclib.elems[0].query.fragment_count(), 3);
    }

    #[test]
    fn test_sample_elution_group_deserialization() {
        let json = SerSpeclibElement::sample_json();
        let elem: SerSpeclibElement = serde_json::from_str(json).unwrap();
        let query_item: QueryItemToScore = elem.into();

        assert_eq!(query_item.digest.len(), "PEPTIDEPINK".len());
        assert_eq!(query_item.query.fragment_count(), 3);
    }
}
