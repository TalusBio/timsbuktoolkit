use crate::errors::LibraryReadingError;
use crate::fragment_mass::elution_group_converter::isotope_dist_from_seq;
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
use std::collections::HashMap;
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
use timsquery::serde::{
    DiannPrecursorExtras,
    ElutionGroupCollection,
    FileReadingExtras,
    read_library_file as read_timsquery_library,
};

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

/// Convert a DIA-NN library entry to a QueryItemToScore
///
/// This handles conversion from TimsElutionGroup + DiannPrecursorExtras
/// to the format expected by the timsseek scoring pipeline.
fn convert_diann_to_query_item(
    eg: TimsElutionGroup<IonAnnot>,
    extra: DiannPrecursorExtras,
    decoy_group_id: u32,
) -> Result<QueryItemToScore, LibraryReadingError> {
    let digest = DigestSlice::from_string(
        extra.stripped_peptide.clone(),
        extra.is_decoy,
        decoy_group_id,
    );

    let fragment_intensities: HashMap<IonAnnot, f32> =
        extra.relative_intensities.into_iter().collect();

    let (rebuilt_eg, precursor_intensities) = match isotope_dist_from_seq(&extra.stripped_peptide) {
        Ok(isotope_ratios) => {
            let fragment_labels: Vec<IonAnnot> =
                eg.iter_fragments().map(|(label, _mz)| *label).collect();
            let fragment_mzs: Vec<f64> = eg.iter_fragments().map(|(_label, mz)| mz).collect();

            let rebuilt = TimsElutionGroup::builder()
                .id(eg.id())
                .mobility_ook0(eg.mobility_ook0())
                .rt_seconds(eg.rt_seconds())
                .precursor(eg.precursor_mz(), eg.precursor_charge())
                .precursor_labels([0i8, 1i8, 2i8].as_slice().into())
                .fragment_labels(fragment_labels.as_slice().into())
                .fragment_mzs(fragment_mzs)
                .try_build()
                .map_err(|e| LibraryReadingError::UnsupportedFormat {
                    message: format!("Failed to rebuild elution group: {:?}", e),
                })?;

            let intensities: HashMap<i8, f32> = vec![
                (0i8, isotope_ratios[0]),
                (1i8, isotope_ratios[1]),
                (2i8, isotope_ratios[2]),
            ]
            .into_iter()
            .collect();

            (rebuilt, intensities)
        }
        Err(e) => {
            tracing::debug!(
                "Failed to calculate isotope distribution for {}: {}. Using monoisotope only.",
                extra.stripped_peptide,
                e
            );
            (eg, HashMap::from_iter([(0i8, 1.0)]))
        }
    };

    let expected_intensity = ExpectedIntensities {
        fragment_intensities,
        precursor_intensities,
    };

    Ok(QueryItemToScore {
        digest,
        query: rebuilt_eg,
        expected_intensity,
    })
}

/// Generate a mass-shifted decoy from a target QueryItemToScore
///
/// Creates a decoy by shifting the precursor m/z by a fixed mass difference.
/// All other properties (RT, mobility, fragments, intensities) remain the same.
fn create_mass_shifted_decoy(
    target: &QueryItemToScore,
    decoy_group_id: u32,
    mass_shift_da: f64,
) -> Result<QueryItemToScore, LibraryReadingError> {
    let target_sequence: String = target.digest.clone().into();
    let decoy_digest = DigestSlice::from_string(target_sequence, true, decoy_group_id);

    let charge = target.query.precursor_charge();
    let mz_shift = mass_shift_da / charge as f64;
    let shifted_precursor_mz = target.query.precursor_mz() + mz_shift;

    let fragment_labels: Vec<IonAnnot> = target
        .query
        .iter_fragments()
        .map(|(label, _)| *label)
        .collect();
    let fragment_mzs: Vec<f64> = target
        .query
        .iter_fragments()
        .map(|(lab, mz)| {
            if let Some(ord) = lab.try_get_ordinal() {
                // Only shift masses if the ordinal is more than position
                // 2
                if ord <= 2 {
                    return mz;
                }
            }
            mz + (mass_shift_da / lab.get_charge() as f64)
        })
        .collect();

    let precursor_labels: Vec<i8> = target
        .query
        .iter_precursors()
        .map(|(label, _)| label)
        .collect();

    let shifted_eg = TimsElutionGroup::builder()
        .id(target.query.id())
        .mobility_ook0(target.query.mobility_ook0())
        .rt_seconds(target.query.rt_seconds())
        .precursor(shifted_precursor_mz, charge)
        .precursor_labels(precursor_labels.as_slice().into())
        .fragment_labels(fragment_labels.as_slice().into())
        .fragment_mzs(fragment_mzs)
        .try_build()
        .map_err(|e| LibraryReadingError::UnsupportedFormat {
            message: format!("Failed to build mass-shifted decoy elution group: {:?}", e),
        })?;

    Ok(QueryItemToScore {
        digest: decoy_digest,
        query: shifted_eg,
        expected_intensity: target.expected_intensity.clone(),
    })
}

/// Convert an ElutionGroupCollection from timsquery to a Speclib (without applying strategy)
fn convert_elution_group_collection(
    collection: ElutionGroupCollection,
) -> Result<Speclib, LibraryReadingError> {
    match collection {
        ElutionGroupCollection::MzpafLabels(egs, Some(FileReadingExtras::Diann(extras))) => {
            if egs.len() != extras.len() {
                return Err(LibraryReadingError::UnsupportedFormat {
                    message: format!(
                        "Mismatch between elution groups ({}) and extras ({})",
                        egs.len(),
                        extras.len()
                    ),
                });
            }

            tracing::info!(
                "Converting {} DIA-NN library entries to Speclib format",
                egs.len()
            );

            // Convert everything as-is, preserving existing decoy groups
            let elems: Vec<_> = egs
                .into_iter()
                .zip(extras)
                .enumerate()
                .map(|(idx, (eg, extra))| {
                    let decoy_group = idx as u32;
                    convert_diann_to_query_item(eg, extra, decoy_group)
                })
                .collect::<Result<Vec<_>, _>>()?;

            Ok(Speclib { elems })
        }
        ElutionGroupCollection::MzpafLabels(_, None) => {
            Err(LibraryReadingError::UnsupportedFormat {
                message: "MzpafLabels variant without DiannPrecursorExtras is not supported"
                    .to_string(),
            })
        }
        _ => {
            tracing::warn!("Unsupported ElutionGroupCollection variant for timsseek");
            Err(LibraryReadingError::UnsupportedFormat {
                message: "Only DIA-NN TSV/Parquet libraries are currently supported via timsquery"
                    .to_string(),
            })
        }
    }
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

    pub fn from_file(
        path: &Path,
        decoy_strategy: crate::models::DecoyStrategy,
    ) -> Result<Self, LibraryReadingError> {
        // First try to detect if it's one of our native formats (msgpack/ndjson)
        match SpeclibFormat::detect_from_path(path) {
            Ok(format) => {
                // Try to load with native format
                match Self::from_file_with_format(path, format, decoy_strategy) {
                    Ok(speclib) => return Ok(speclib),
                    Err(e) => {
                        tracing::debug!(
                            "Failed to load as native format, trying timsquery fallback: {:?}",
                            e
                        );
                    }
                }
            }
            Err(_) => {
                tracing::debug!("Could not detect native format, trying timsquery fallback");
            }
        }

        tracing::info!(
            "Attempting to load library using timsquery format detection: {}",
            path.display()
        );
        let collection = read_timsquery_library(path)?;
        let speclib = convert_elution_group_collection(collection)?;

        // Apply decoy strategy to the loaded speclib
        Self::apply_decoy_strategy(speclib, decoy_strategy)
    }

    /// Apply decoy strategy to an already-loaded speclib (works for any format)
    fn apply_decoy_strategy(
        speclib: Speclib,
        strategy: crate::models::DecoyStrategy,
    ) -> Result<Speclib, LibraryReadingError> {
        use crate::models::{
            DecoyMarking,
            DecoyStrategy,
        };

        // Separate targets from decoys
        let (targets, decoys): (Vec<_>, Vec<_>) = speclib
            .elems
            .into_iter()
            .partition(|item| item.digest.decoy == DecoyMarking::Target);

        let has_decoys = !decoys.is_empty();

        match strategy {
            DecoyStrategy::Never => {
                // Return as-is (targets + existing decoys)
                tracing::info!(
                    "Decoy strategy: Never - using library as-is ({} targets, {} decoys)",
                    targets.len(),
                    decoys.len()
                );
                Ok(Speclib {
                    elems: targets.into_iter().chain(decoys).collect(),
                })
            }

            DecoyStrategy::Force => {
                // Drop existing decoys and generate new mass-shift decoys
                if has_decoys {
                    tracing::warn!(
                        "Decoy strategy: Force - dropping {} existing decoys and regenerating mass-shift decoys",
                        decoys.len()
                    );
                } else {
                    tracing::info!("Decoy strategy: Force - generating mass-shift decoys");
                }

                // const CARBON_MASS: f64 = 12.0;
                // This seems to leak data.
                const CH2_MASS: f64 = 14.0;
                let mut all_entries = Vec::with_capacity(targets.len() * 3);

                for (idx, target) in targets.into_iter().enumerate() {
                    let decoy_group_id = idx as u32;

                    // Update target's decoy_group
                    let mut target = target;
                    target.digest.decoy_group = decoy_group_id;

                    all_entries.push(target.clone());

                    let plus_decoy = create_mass_shifted_decoy(&target, decoy_group_id, CH2_MASS)?;
                    all_entries.push(plus_decoy);

                    let minus_decoy =
                        create_mass_shifted_decoy(&target, decoy_group_id, -CH2_MASS)?;
                    all_entries.push(minus_decoy);
                }

                tracing::info!(
                    "Generated {} total entries ({} targets + {} decoys)",
                    all_entries.len(),
                    all_entries.len() / 3,
                    all_entries.len() * 2 / 3
                );

                Ok(Speclib { elems: all_entries })
            }

            DecoyStrategy::IfMissing => {
                if has_decoys {
                    tracing::info!(
                        "Library contains {} targets and {} decoys. Using existing decoys.",
                        targets.len(),
                        decoys.len()
                    );
                    Ok(Speclib {
                        elems: targets.into_iter().chain(decoys).collect(),
                    })
                } else {
                    tracing::warn!(
                        "Library contains NO decoys. Will generate synthetic mass-shift decoys:\n\
                         - Creating 2 decoys per target (±12.0 Da / C12 mass)\n\
                         - Total library size will be 3x original (1 target + 2 decoys)\n\
                         - Each target-decoy triplet will share a decoy_group for FDR estimation"
                    );

                    const CARBON_MASS: f64 = 12.0;
                    let mut all_entries = Vec::with_capacity(targets.len() * 3);

                    for (idx, target) in targets.into_iter().enumerate() {
                        let decoy_group_id = idx as u32;

                        // Update target's decoy_group
                        let mut target = target;
                        target.digest.decoy_group = decoy_group_id;

                        all_entries.push(target.clone());

                        let plus_decoy =
                            create_mass_shifted_decoy(&target, decoy_group_id, CARBON_MASS)?;
                        all_entries.push(plus_decoy);

                        let minus_decoy =
                            create_mass_shifted_decoy(&target, decoy_group_id, -CARBON_MASS)?;
                        all_entries.push(minus_decoy);
                    }

                    tracing::info!(
                        "Generated {} total entries ({} targets + {} decoys)",
                        all_entries.len(),
                        all_entries.len() / 3,
                        all_entries.len() * 2 / 3
                    );

                    Ok(Speclib { elems: all_entries })
                }
            }
        }
    }

    pub fn from_file_with_format(
        path: &Path,
        format: SpeclibFormat,
        decoy_strategy: crate::models::DecoyStrategy,
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

        let speclib = Self { elems };

        // Apply decoy strategy to the loaded speclib
        Self::apply_decoy_strategy(speclib, decoy_strategy)
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

    #[test]
    fn test_load_diann_tsv_library() {
        // Use the test file from timsquery tests
        // Note: sample_lib.tsv is in Skyline format and won't load as DIA-NN
        // So we test with sample_lib.txt which is in DIA-NN TSV format
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        assert!(
            test_file.exists(),
            "Test file should exist at {:?}",
            test_file
        );

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        // The sample file has 2 unique precursors with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys
        assert_eq!(
            speclib.len(),
            6,
            "Expected 6 entries (2 targets + 4 decoys)"
        );

        // Verify first target entry structure (targets come first)
        let first_target = speclib
            .elems
            .iter()
            .find(|e| !e.digest.is_decoy())
            .expect("Should have at least one target");
        assert_eq!(
            first_target.query.fragment_count(),
            5,
            "First entry should have 5 fragments"
        );

        // Verify isotope envelope was added (should have 3 isotopes: 0, 1, 2)
        let precursor_count: usize = first_target.query.iter_precursors().count();
        assert_eq!(
            precursor_count, 3,
            "Should have 3 isotopes in precursor envelope"
        );

        // Verify precursor intensities match isotope count
        assert_eq!(
            first_target.expected_intensity.precursor_intensities.len(),
            3,
            "Should have intensities for all 3 isotopes"
        );

        // Verify all intensities are present
        assert!(
            first_target
                .expected_intensity
                .precursor_intensities
                .contains_key(&0),
            "Should have intensity for isotope 0"
        );
        assert!(
            first_target
                .expected_intensity
                .precursor_intensities
                .contains_key(&1),
            "Should have intensity for isotope 1"
        );
        assert!(
            first_target
                .expected_intensity
                .precursor_intensities
                .contains_key(&2),
            "Should have intensity for isotope 2"
        );
    }

    #[test]
    fn test_load_diann_txt_library() {
        // Use the test file from timsquery tests
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        assert!(
            test_file.exists(),
            "Test file should exist at {:?}",
            test_file
        );

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load TXT library");

        // The sample file has 2 unique precursors with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys
        assert_eq!(
            speclib.len(),
            6,
            "Expected 6 entries (2 targets + 4 decoys)"
        );

        // Verify we have both targets and decoys
        let n_targets = speclib
            .elems
            .iter()
            .filter(|e| !e.digest.is_decoy())
            .count();
        let n_decoys = speclib.elems.iter().filter(|e| e.digest.is_decoy()).count();

        assert_eq!(n_targets, 2, "Should have 2 targets");
        assert_eq!(n_decoys, 4, "Should have 4 decoys");
    }

    #[test]
    fn test_load_diann_parquet_library() {
        // Use the test file from timsquery tests
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_pq_speclib.parquet");

        assert!(
            test_file.exists(),
            "Test file should exist at {:?}",
            test_file
        );

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load Parquet library");

        // The sample parquet file has 3 unique precursors with no decoys
        // Should generate 3x entries: 3 targets + 6 decoys
        assert_eq!(
            speclib.len(),
            9,
            "Expected 9 entries (3 targets + 6 decoys)"
        );

        // Verify we have both targets and decoys
        let n_targets = speclib
            .elems
            .iter()
            .filter(|e| !e.digest.is_decoy())
            .count();
        let n_decoys = speclib.elems.iter().filter(|e| e.digest.is_decoy()).count();

        assert_eq!(n_targets, 3, "Should have 3 targets");
        assert_eq!(n_decoys, 6, "Should have 6 decoys");

        // Verify isotope envelope for targets
        for entry in speclib.elems.iter().filter(|e| !e.digest.is_decoy()) {
            let precursor_count: usize = entry.query.iter_precursors().count();
            assert_eq!(
                precursor_count, 3,
                "Each entry should have 3 isotopes in precursor envelope"
            );
        }
    }

    #[test]
    fn test_isotope_envelope_calculation() {
        // Use the DIA-NN TSV test file
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        // Check that isotope intensities are normalized (M0 should be 1.0)
        for entry in &speclib.elems {
            let m0_intensity = entry
                .expected_intensity
                .precursor_intensities
                .get(&0)
                .expect("Should have M0 intensity");

            // M0 should be the highest intensity (normalized to 1.0)
            assert!(
                (*m0_intensity - 1.0).abs() < 0.01,
                "M0 should be normalized to ~1.0, got {}",
                m0_intensity
            );

            // M+1 and M+2 should be less than M0
            let m1_intensity = entry.expected_intensity.precursor_intensities.get(&1);
            let m2_intensity = entry.expected_intensity.precursor_intensities.get(&2);

            if let Some(m1) = m1_intensity {
                assert!(*m1 <= 1.0, "M+1 intensity should be <= 1.0, got {}", m1);
            }
            if let Some(m2) = m2_intensity {
                assert!(*m2 <= 1.0, "M+2 intensity should be <= 1.0, got {}", m2);
            }
        }
    }

    #[test]
    fn test_decoy_generation_for_library_without_decoys() {
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        // Test file has 2 targets with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys (2 per target)
        assert_eq!(
            speclib.len(),
            6,
            "Should have 6 entries (2 targets + 4 decoys)"
        );

        // Count targets vs decoys
        let n_targets = speclib
            .elems
            .iter()
            .filter(|e| !e.digest.is_decoy())
            .count();
        let n_decoys = speclib.elems.iter().filter(|e| e.digest.is_decoy()).count();

        assert_eq!(n_targets, 2, "Should have 2 target entries");
        assert_eq!(n_decoys, 4, "Should have 4 decoy entries (2 per target)");

        // Each decoy_group should have exactly 3 entries (1 target + 2 decoys)
        let mut groups: std::collections::HashMap<u32, Vec<&crate::QueryItemToScore>> =
            std::collections::HashMap::new();
        for entry in &speclib.elems {
            groups
                .entry(entry.digest.decoy_group)
                .or_insert_with(Vec::new)
                .push(entry);
        }

        assert_eq!(groups.len(), 2, "Should have 2 unique decoy groups");

        for (group_id, entries) in groups {
            assert_eq!(
                entries.len(),
                3,
                "Decoy group {} should have exactly 3 entries",
                group_id
            );

            // Verify 1 target + 2 decoys per group
            let targets_in_group = entries.iter().filter(|e| !e.digest.is_decoy()).count();
            let decoys_in_group = entries.iter().filter(|e| e.digest.is_decoy()).count();

            assert_eq!(targets_in_group, 1, "Should have 1 target in group");
            assert_eq!(decoys_in_group, 2, "Should have 2 decoys in group");
        }
    }

    #[test]
    fn test_mass_shift_decoys() {
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        // Group by decoy_group_id
        let mut groups: std::collections::HashMap<u32, Vec<&crate::QueryItemToScore>> =
            std::collections::HashMap::new();
        for entry in &speclib.elems {
            groups
                .entry(entry.digest.decoy_group)
                .or_insert_with(Vec::new)
                .push(entry);
        }

        const CARBON_MASS: f64 = 12.0;

        for (_group_id, entries) in groups {
            // Find target and decoys
            let target = entries
                .iter()
                .find(|e| !e.digest.is_decoy())
                .expect("Should have a target");
            let decoys: Vec<_> = entries.iter().filter(|e| e.digest.is_decoy()).collect();

            assert_eq!(decoys.len(), 2, "Should have 2 decoys");

            let target_mz = target.query.precursor_mz();
            let charge = target.query.precursor_charge();
            let mz_shift = CARBON_MASS / charge as f64;

            // Check that decoys have expected m/z shifts
            let decoy_mzs: Vec<f64> = decoys.iter().map(|d| d.query.precursor_mz()).collect();

            let expected_plus = target_mz + mz_shift;
            let expected_minus = target_mz - mz_shift;

            // One decoy should be +C12, one should be -C12
            let has_plus = decoy_mzs
                .iter()
                .any(|mz| (mz - expected_plus).abs() < 0.001);
            let has_minus = decoy_mzs
                .iter()
                .any(|mz| (mz - expected_minus).abs() < 0.001);

            assert!(
                has_plus,
                "Should have a +C12 decoy at m/z {}, found {:?}",
                expected_plus, decoy_mzs
            );
            assert!(
                has_minus,
                "Should have a -C12 decoy at m/z {}, found {:?}",
                expected_minus, decoy_mzs
            );

            // Verify other properties are preserved
            for decoy in &decoys {
                assert_eq!(
                    decoy.query.rt_seconds(),
                    target.query.rt_seconds(),
                    "RT should be preserved"
                );
                assert_eq!(
                    decoy.query.mobility_ook0(),
                    target.query.mobility_ook0(),
                    "Mobility should be preserved"
                );
                assert_eq!(
                    decoy.query.fragment_count(),
                    target.query.fragment_count(),
                    "Fragment count should be preserved"
                );
            }
        }
    }

    #[test]
    fn test_fragment_intensities_preserved() {
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        for entry in &speclib.elems {
            // Number of fragment intensities should match number of fragments
            assert_eq!(
                entry.expected_intensity.fragment_intensities.len(),
                entry.query.fragment_count(),
                "Fragment intensity count should match fragment count"
            );

            // All fragment intensities should be positive
            for (_label, intensity) in &entry.expected_intensity.fragment_intensities {
                assert!(*intensity > 0.0, "Fragment intensities should be positive");
            }
        }
    }
}
