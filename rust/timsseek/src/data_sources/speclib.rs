use crate::data_sources::reference_library::{
    RefQuery,
    ReferenceLibrary,
};
use crate::errors::LibraryReadingError;
use crate::fragment_mass::elution_group_converter::isotope_dist_from_seq;
use crate::fragment_mass::{
    IsotopeSource,
    isotope_dist_or_averagine,
};
use crate::models::sequence::{
    Peptide,
    SeqFormat,
    SpeclibMeta,
    normalize_to_proforma,
    parse_sequence,
};
use crate::models::{
    DecoyMarking,
    map_decoy_strategy,
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
use timsquery::models::QueryCollection;
use timsquery::models::capabilities::{
    LibCapabilities,
    SeqFeatureState,
};
use timsquery::models::elution_group::TimsElutionGroup;
use timsquery::serde::{
    DiannPrecursorExtras,
    ElutionGroupCollection,
    FileReadingExtras,
    SkylinePrecursorExtras,
    read_library_file as read_timsquery_library,
};
use timsquery::utils::constants::PROTON_MASS;

/// This is meant to the be the serializable version of the speclib element
/// so ... in general should be backwards compatible and implement `Into<QueryItemToScore>`
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerSpeclibElement {
    precursor: PrecursorEntry,
    elution_group: ReferenceEG,
}

impl SerSpeclibElement {
    pub fn new(precursor: PrecursorEntry, elution_group: ReferenceEG) -> Self {
        Self {
            precursor,
            elution_group,
        }
    }

    pub fn sample() -> Self {
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
pub struct PrecursorEntry {
    sequence: String,
    charge: u8,
    decoy: bool,
    decoy_group: u32,
}

impl PrecursorEntry {
    pub fn new(sequence: String, charge: u8, decoy: bool, decoy_group: u32) -> Self {
        Self {
            sequence,
            charge,
            decoy,
            decoy_group,
        }
    }
}

impl TryFrom<SerSpeclibElement> for QueryItemToScore {
    type Error = LibraryReadingError;

    fn try_from(x: SerSpeclibElement) -> Result<Self, Self::Error> {
        let charge = x.precursor.charge;
        let raw: Arc<str> = x.precursor.sequence.clone().into();
        let normalized = normalize_to_proforma(&raw);
        let parsed = parse_sequence(&normalized);
        let precursor = Peptide {
            raw,
            parsed,
            decoy: if x.precursor.decoy {
                DecoyMarking::ReversedDecoy
            } else {
                DecoyMarking::Target
            },
            decoy_group: x.precursor.decoy_group,
        };
        let ref_eg = x.elution_group;
        let expected_intensity = ExpectedIntensities::try_from_pairs(
            ref_eg
                .fragment_labels
                .iter()
                .cloned()
                .zip(ref_eg.fragment_intensities.iter().cloned()),
            ref_eg
                .precursor_labels
                .iter()
                .cloned()
                .zip(ref_eg.precursor_intensities.iter().cloned()),
        )
        .map_err(|e| LibraryReadingError::UnsupportedFormat {
            message: format!("speclib entry has {}", e),
        })?;
        Ok(QueryItemToScore {
            expected_intensity,
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
        })
    }
}

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

impl ReferenceEG {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: u32,
        precursor_mz: f64,
        precursor_labels: Vec<i8>,
        fragment_mzs: Vec<f64>,
        fragment_labels: Vec<IonAnnot>,
        precursor_intensities: Vec<f32>,
        fragment_intensities: Vec<f32>,
        mobility_ook0: f32,
        rt_seconds: f32,
    ) -> Self {
        Self {
            id,
            precursor_mz,
            precursor_labels,
            fragment_mzs,
            fragment_labels,
            precursor_intensities,
            fragment_intensities,
            mobility_ook0,
            rt_seconds,
        }
    }
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
    // Prefer modified form for raw; fall back to stripped if empty.
    let raw_src = if extra.modified_peptide.is_empty() {
        extra.stripped_peptide.clone()
    } else {
        extra.modified_peptide.clone()
    };
    let raw: Arc<str> = raw_src.into();
    let normalized = normalize_to_proforma(&raw);
    let parsed = parse_sequence(&normalized);
    let digest = Peptide {
        raw,
        parsed,
        decoy: if extra.is_decoy {
            DecoyMarking::ReversedDecoy
        } else {
            DecoyMarking::Target
        },
        decoy_group: decoy_group_id,
    };

    let (rebuilt_eg, precursor_pairs) = match isotope_dist_from_seq(&extra.stripped_peptide) {
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

            let pairs: Vec<(i8, f32)> = vec![
                (0i8, isotope_ratios[0]),
                (1i8, isotope_ratios[1]),
                (2i8, isotope_ratios[2]),
            ];
            (rebuilt, pairs)
        }
        Err(e) => {
            tracing::debug!(
                "Failed to calculate isotope distribution for {}: {}. Using monoisotope only.",
                extra.stripped_peptide,
                e
            );
            (eg, vec![(0i8, 1.0)])
        }
    };

    let expected_intensity =
        ExpectedIntensities::try_from_pairs(extra.relative_intensities, precursor_pairs).map_err(
            |e| LibraryReadingError::UnsupportedFormat {
                message: format!("DIA-NN entry {}: {}", extra.stripped_peptide, e),
            },
        )?;

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
    let decoy_digest = Peptide {
        raw: target.digest.raw.clone(),
        parsed: target.digest.parsed.clone(),
        decoy: DecoyMarking::MassShiftedDecoy,
        decoy_group: decoy_group_id,
    };

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

/// `SkylinePrecursorExtras` is structurally identical to `DiannPrecursorExtras`
/// (peptide + protein + decoy flag + fragment intensity pairs), so we adapt it
/// into the DIA-NN shape and reuse `convert_diann_to_query_item`.
fn skyline_to_diann_extras(sky: SkylinePrecursorExtras) -> DiannPrecursorExtras {
    DiannPrecursorExtras {
        modified_peptide: sky.modified_peptide,
        stripped_peptide: sky.stripped_peptide,
        protein_id: sky.protein_id,
        is_decoy: sky.is_decoy,
        relative_intensities: sky.relative_intensities,
    }
}

/// Map library rows to `elems`, in parallel when the `rayon` feature is on.
///
/// `f` maps `(row_index, eg, extra) -> Result<QueryItemToScore>`. Order and the
/// row index are preserved in both builds (indexed par-iter + order-preserving
/// `Result` collect), so `decoy_group = idx` and every downstream flat-index
/// assumption match the serial result exactly. The per-row work (`parse_sequence`
/// + `isotope_dist_from_seq`) is pure, so parallelizing is safe.
fn convert_rows<E, F>(
    egs: Vec<TimsElutionGroup<IonAnnot>>,
    extras: Vec<E>,
    f: F,
) -> Result<Vec<QueryItemToScore>, LibraryReadingError>
where
    E: Send,
    F: Fn(usize, TimsElutionGroup<IonAnnot>, E) -> Result<QueryItemToScore, LibraryReadingError>
        + Sync
        + Send,
{
    #[cfg(feature = "rayon")]
    {
        use rayon::prelude::*;
        egs.into_par_iter()
            .zip(extras)
            .enumerate()
            .map(|(idx, (eg, extra))| f(idx, eg, extra))
            .collect()
    }
    #[cfg(not(feature = "rayon"))]
    {
        egs.into_iter()
            .zip(extras)
            .enumerate()
            .map(|(idx, (eg, extra))| f(idx, eg, extra))
            .collect()
    }
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
            let elems = convert_rows(egs, extras, |idx, eg, extra| {
                convert_diann_to_query_item(eg, extra, idx as u32)
            })?;

            Ok(Speclib::Materialized {
                elems,
                meta: SpeclibMeta::default(),
            })
        }
        ElutionGroupCollection::MzpafLabels(egs, Some(FileReadingExtras::Skyline(extras))) => {
            if egs.len() != extras.len() {
                return Err(LibraryReadingError::UnsupportedFormat {
                    message: format!(
                        "Mismatch between Skyline elution groups ({}) and extras ({})",
                        egs.len(),
                        extras.len()
                    ),
                });
            }

            tracing::info!(
                "Converting {} Skyline transition-list entries to Speclib format",
                egs.len()
            );

            let elems = convert_rows(egs, extras, |idx, eg, sky_extra| {
                convert_diann_to_query_item(eg, skyline_to_diann_extras(sky_extra), idx as u32)
            })?;

            Ok(Speclib::Materialized {
                elems,
                meta: SpeclibMeta::default(),
            })
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
                message: "Only DIA-NN TSV/Parquet and Skyline CSV libraries are currently supported via timsquery"
                    .to_string(),
            })
        }
    }
}

/// Summary of a [`build_reference_library`] call, for load-time logging.
#[derive(Debug, Clone, Copy)]
pub struct LoadReport {
    pub n_targets: usize,
    pub n_averagine_fallback: usize,
    pub sequence_features: SeqFeatureState,
}

/// Whether any row in the collection's extras is already flagged as a decoy.
///
/// Used to decide (`map_decoy_strategy`) whether `IfMissing` should generate
/// lazy mass-shift decoys or defer to the file's own decoy rows.
fn collection_has_file_decoys(collection: &ElutionGroupCollection) -> bool {
    match collection {
        ElutionGroupCollection::MzpafLabels(_, Some(FileReadingExtras::Diann(extras))) => {
            extras.iter().any(|e| e.is_decoy)
        }
        ElutionGroupCollection::MzpafLabels(_, Some(FileReadingExtras::Skyline(extras))) => {
            extras.iter().any(|e| e.is_decoy)
        }
        ElutionGroupCollection::MzpafLabels(_, Some(FileReadingExtras::Spectronaut(extras))) => {
            extras.iter().any(|e| e.is_decoy)
        }
        _ => false,
    }
}

/// Build the lazy, columnar `ReferenceLibrary` arena directly from the
/// timsquery bridge output, without ever materializing a `QueryItemToScore`
/// per row. This is the memory-optimized path for the DEFAULT `.speclib`
/// load (see `speclib_data_flow.md`): it is what avoids the 9 GB peak RSS of
/// the fully-materialized target+2-decoy expansion.
///
/// Only the DIA-NN `MzpafLabels` + `FileReadingExtras::Diann` shape is
/// supported today (mirrors `convert_elution_group_collection`'s DIA-NN
/// branch); everything else stays on the materialized path in `from_file`.
pub fn build_reference_library(
    collection: ElutionGroupCollection,
    caps: LibCapabilities,
) -> Result<(ReferenceLibrary, LoadReport), LibraryReadingError> {
    let (egs, extras) = match collection {
        ElutionGroupCollection::MzpafLabels(egs, Some(FileReadingExtras::Diann(extras))) => {
            (egs, extras)
        }
        _ => {
            return Err(LibraryReadingError::UnsupportedFormat {
                message:
                    "build_reference_library only supports DIA-NN MzpafLabels collections (with extras)"
                        .to_string(),
            });
        }
    };

    if egs.len() != extras.len() {
        return Err(LibraryReadingError::UnsupportedFormat {
            message: format!(
                "Mismatch between elution groups ({}) and extras ({})",
                egs.len(),
                extras.len()
            ),
        });
    }

    let n_targets = egs.len();
    tracing::info!(
        "Building lazy ReferenceLibrary arena from {} DIA-NN library entries",
        n_targets
    );

    let mut geom = QueryCollection::with_capabilities(caps);
    let mut frag_intens: Vec<f32> = Vec::with_capacity(n_targets * 8);
    let mut n_averagine_fallback = 0usize;
    let mut all_parsable = true;

    for (eg, extra) in egs.into_iter().zip(extras.into_iter()) {
        let frags: Vec<(IonAnnot, f64)> = eg.iter_fragments().map(|(lab, mz)| (*lab, mz)).collect();

        // Ingest assertion: `frag_intens` is co-built with `frag_labels` from
        // the same source pair, in the same order. If DIA-NN ever reorders
        // one relative to the other, `RefQuery::iter_expected_fragments`
        // would silently zip mismatched (label, intensity) pairs downstream.
        debug_assert!(
            frags
                .iter()
                .map(|(l, _)| *l)
                .eq(extra.relative_intensities.iter().map(|(l, _)| *l)),
            "frag label order must match intensity key order at ingest"
        );

        let stripped = extra.stripped_peptide.as_str();
        let modified = if extra.modified_peptide.is_empty() {
            extra.stripped_peptide.as_str()
        } else {
            extra.modified_peptide.as_str()
        };

        if all_parsable {
            let normalized = normalize_to_proforma(modified);
            if parse_sequence(&normalized).is_none() {
                all_parsable = false;
            }
        }

        let charge = eg.precursor_charge();
        let neutral_mass = eg.precursor_mz() * charge as f64 - charge as f64 * PROTON_MASS;
        let (isotope_src, _envelope) = isotope_dist_or_averagine(stripped, neutral_mass);
        if isotope_src == IsotopeSource::Averagine {
            n_averagine_fallback += 1;
        }

        geom.push_target(
            eg.precursor_mz(),
            charge,
            eg.rt_seconds(),
            eg.mobility_ook0(),
            &frags,
            stripped,
            modified,
            &[],
        );

        frag_intens.extend(extra.relative_intensities.iter().map(|(_, i)| *i));
    }

    geom.seal();

    let sequence_features = if all_parsable {
        SeqFeatureState::Available
    } else {
        SeqFeatureState::Unavailable
    };
    geom.caps.sequence_features = sequence_features;

    if n_averagine_fallback > 0 {
        tracing::warn!(
            "{}/{} library entries fell back to averagine isotope envelopes",
            n_averagine_fallback,
            n_targets
        );
    }
    tracing::info!(
        "Lazy ReferenceLibrary arena built: {} targets, sequence_features={:?}",
        n_targets,
        sequence_features
    );

    let report = LoadReport {
        n_targets,
        n_averagine_fallback,
        sequence_features,
    };

    Ok((ReferenceLibrary { geom, frag_intens }, report))
}

#[derive(Debug, Clone)]
pub enum Speclib {
    /// Columnar arena path (see `speclib_data_flow.md`): default route for a
    /// `.speclib` load whose decoy strategy resolves to `LazyMassShift`.
    /// Scoring against this arm is wired up in Task 9/10 — for now several
    /// `Speclib` methods panic on this arm; see their docs.
    Lazy(ReferenceLibrary),
    /// Legacy fully-materialized path: one `QueryItemToScore` per row
    /// (target + every decoy variant). Kept for `Passthrough`/`None` decoy
    /// strategies and the native `SerSpeclibElement` formats until Task 11
    /// deletes it.
    Materialized {
        elems: Vec<QueryItemToScore>,
        meta: SpeclibMeta,
    },
}

struct SpeclibIterator<'a> {
    speclib: &'a Speclib,
    idx: usize,
}

impl SpeclibIterator<'_> {
    /// Both call sites need the materialized slice; this is the one place
    /// that panics on the `Lazy` arm today.
    ///
    /// `Speclib::iter()` yields owned `QueryItemToScore`, which the `Lazy`
    /// arm has none of (it hands out `RefQuery` flyweights instead via
    /// `Speclib::item_at`/`ReferenceLibrary::iter`). No current caller
    /// exercises `Speclib::iter()` on a `Lazy` lib; unifying the two
    /// iteration shapes is Task 9/10's job.
    fn elems(&self) -> &[QueryItemToScore] {
        match self.speclib {
            Speclib::Materialized { elems, .. } => elems,
            Speclib::Lazy(_) => unimplemented!(
                "Speclib::iter() over the Lazy arm is not supported; use ReferenceLibrary::iter (Task 9/10)"
            ),
        }
    }
}

impl Iterator for SpeclibIterator<'_> {
    type Item = QueryItemToScore;

    fn next(&mut self) -> Option<Self::Item> {
        let elems = self.elems();
        if self.idx >= elems.len() {
            return None;
        }
        let elem = elems[self.idx].clone();
        self.idx += 1;
        Some(elem)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.elems().len() - self.idx;
        (remaining, Some(remaining))
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
    /// Detect a native timsseek format by EXTENSION ONLY. Returns `None` for
    /// anything else (including `.speclib`), which routes to the timsquery
    /// bridge.
    ///
    /// Extension-only is deliberate: msgpack has no reliable magic byte, so a
    /// content sniff would misclaim raw binaries like `.speclib` as msgpack.
    pub fn detect_from_extension(path: &Path) -> Option<Self> {
        let path_str = path.to_string_lossy().to_lowercase();

        // Accept both `.zst` and `.zstd` — DIA-NN/user pipelines use either.
        if path_str.ends_with(".msgpack.zst") || path_str.ends_with(".msgpack.zstd") {
            Some(SpeclibFormat::MessagePackZstd)
        } else if path_str.ends_with(".msgpack") {
            Some(SpeclibFormat::MessagePack)
        } else if path_str.ends_with(".ndjson.zst") || path_str.ends_with(".ndjson.zstd") {
            Some(SpeclibFormat::NdJsonZstd)
        } else if path_str.ends_with(".ndjson") {
            Some(SpeclibFormat::NdJson)
        } else {
            None
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

                Some(elem.try_into())
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
            Ok(elem) => Some(elem.try_into()),
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

/// Inspect every `Peptide.parsed`. If any is `None`, zero them all, flip
/// `parsable_sequences` off, log counts. Returns the finalized meta.
///
/// MUST run AFTER all decoy generation (including `create_mass_shifted_decoy`
/// and DIA-NN target-decoy pairing). Otherwise decoys generated post-gate
/// inherit `parsed: Some(...)` from their target while the gate said off.
/// Whether every entry parsed. Cheap to compute on TARGETS before decoy
/// expansion (decoys clone their target's `parsed`, so the verdict is identical)
/// — doing it there avoids re-walking the 3x-larger, peak-RSS `elems`.
fn all_sequences_parsed(items: &[QueryItemToScore]) -> bool {
    !items.is_empty() && items.iter().all(|q| q.digest.parsed.is_some())
}

/// Apply the all-or-nothing parse gate given a precomputed `parsable` verdict
/// (see [`all_sequences_parsed`]). When parsable, this is O(1) — no full walk.
/// When not, it must still zero every entry's `parsed` (the degraded path).
fn apply_parse_gate(
    items: &mut [QueryItemToScore],
    expected_format: SeqFormat,
    parsable: bool,
) -> SpeclibMeta {
    let total = items.len();
    if !parsable {
        tracing::warn!(
            "speclib load: n_total={}; not all sequences parsable, sequence features disabled",
            total,
        );
        for q in items.iter_mut() {
            q.digest.parsed = None;
        }
    } else {
        tracing::info!(
            "speclib load: n_total={}, all parsed; sequence features enabled",
            total
        );
    }
    SpeclibMeta {
        parsable_sequences: parsable,
        sequence_format: expected_format,
    }
}

/// Target number of entries to sample for diagnostic stats. A full walk of a
/// 26M-entry library at peak RSS costs tens of seconds (paging) for numbers that
/// are only logged; a strided sample of this many is statistically equivalent.
const STATS_SAMPLE_TARGET: usize = 100_000;

impl Speclib {
    pub fn iter<'a>(&'a self) -> impl Iterator<Item = QueryItemToScore> + 'a {
        SpeclibIterator {
            speclib: self,
            idx: 0,
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Speclib::Materialized { elems, .. } => elems.len(),
            Speclib::Lazy(lib) => lib.len(),
        }
    }

    /// Whether every sequence in the library parsed (gates sequence-derived
    /// scoring features). The single call site that used to read
    /// `speclib.meta.parsable_sequences` directly now goes through this
    /// accessor so it works for both arms.
    pub fn parsable_sequences(&self) -> bool {
        match self {
            Speclib::Materialized { meta, .. } => meta.parsable_sequences,
            Speclib::Lazy(lib) => lib.geom.caps.sequence_features == SeqFeatureState::Available,
        }
    }

    /// Flat-index accessor into the `Lazy` arena, in `target, +, -` variant
    /// order (see `ReferenceLibrary::item_at`). The `Materialized` arm has no
    /// `RefQuery` flyweight to hand out yet — unifying both arms behind one
    /// scoring-facing accessor is Task 9/10's job.
    pub fn item_at(&self, flat_idx: usize) -> RefQuery<'_> {
        match self {
            Speclib::Lazy(lib) => lib.item_at(flat_idx),
            Speclib::Materialized { .. } => unimplemented!(
                "Speclib::item_at on the Materialized arm; scoring unification lands in Task 9/10"
            ),
        }
    }

    /// Mean number of fragments per entry (0.0 for an empty library).
    ///
    /// Estimated from a strided sample: a full walk of a multi-million-entry
    /// library touches gigabytes at peak RSS (seconds under paging) for a number
    /// that is only reported. The sample mean is statistically identical.
    pub fn avg_fragments(&self) -> f64 {
        match self {
            Speclib::Materialized { elems, .. } => {
                let n = elems.len();
                if n == 0 {
                    return 0.0;
                }
                let stride = (n / STATS_SAMPLE_TARGET).max(1);
                let (mut sum, mut cnt) = (0usize, 0usize);
                let mut i = 0;
                while i < n {
                    sum += elems[i].expected_intensity.fragment_len();
                    cnt += 1;
                    i += stride;
                }
                sum as f64 / cnt as f64
            }
            // Every variant of a target shares the same fragment set, so the
            // per-target arena length is exact, not a sample — cheaper than
            // the materialized path's stride, not an approximation of it.
            Speclib::Lazy(lib) => {
                let n_targets = lib.geom.n_targets();
                if n_targets == 0 {
                    return 0.0;
                }
                lib.geom.frag_labels.len() as f64 / n_targets as f64
            }
        }
    }

    pub fn from_file(
        path: &Path,
        decoy_strategy: crate::models::DecoyStrategy,
    ) -> Result<Self, LibraryReadingError> {
        // Native timsseek formats are matched by EXTENSION ONLY: a native
        // extension commits to the native reader and surfaces its error. A
        // `.speclib` matches no native extension and falls through to the
        // bridge -> timsquery registry -> binary reader.
        if let Some(format) = SpeclibFormat::detect_from_extension(path) {
            tracing::info!(
                "Loading native speclib format ({:?}) from {}",
                format,
                path.display()
            );
            return Self::from_file_with_format(path, format, decoy_strategy);
        }

        // Terminal source: bridge to the timsquery reader registry (DIA-NN
        // TSV/parquet/.speclib, Spectronaut, Skyline, JSON), then convert.
        tracing::info!(
            "Loading library via timsquery format detection: {}",
            path.display()
        );
        let collection = read_timsquery_library(path)?;

        // Route decision: the memory-optimized lazy arena is isolated to
        // exactly the workload with the 9 GB materialized-decoy problem
        // (`LazyMassShift`) AND the DIA-NN collection shape `build_reference_library`
        // actually knows how to ingest today. `Passthrough`/`None`, and any
        // non-DIA-NN shape (Skyline, Spectronaut, native formats), keep the
        // existing materialized path unchanged — see `map_decoy_strategy`.
        let has_file_decoys = collection_has_file_decoys(&collection);
        let is_diann_shape = matches!(
            collection,
            ElutionGroupCollection::MzpafLabels(_, Some(FileReadingExtras::Diann(_)))
        );
        use timsquery::models::capabilities::DecoyStrategy as TqDecoyStrategy;
        let mapped = map_decoy_strategy(decoy_strategy, has_file_decoys);
        match mapped {
            TqDecoyStrategy::LazyMassShift { offset, n_decoys } if is_diann_shape => {
                let mut caps = LibCapabilities::default_diann();
                caps.decoys = TqDecoyStrategy::LazyMassShift { offset, n_decoys };
                let (lib, report) = build_reference_library(collection, caps)?;
                tracing::info!(
                    "Lazy speclib arena ready: {} targets, sequence_features={:?}",
                    report.n_targets,
                    report.sequence_features
                );
                if report.n_averagine_fallback > 0 {
                    tracing::warn!(
                        "{}/{} library entries used averagine isotope fallback",
                        report.n_averagine_fallback,
                        report.n_targets
                    );
                }
                Ok(Speclib::Lazy(lib))
            }
            // `Passthrough`/`None`, or a `LazyMassShift` resolution on a
            // non-DIA-NN shape that `build_reference_library` can't ingest
            // yet (e.g. Skyline) — same materialized path as before.
            _ => {
                let speclib = convert_elution_group_collection(collection)?;

                // Decide the parse gate on targets, BEFORE the 3x decoy
                // expansion and peak RSS — decoys inherit their target's
                // `parsed`, so the verdict is identical but the walk is 3x
                // smaller and not under memory pressure.
                let Speclib::Materialized { elems, .. } = &speclib else {
                    unreachable!("convert_elution_group_collection always returns Materialized")
                };
                let parsable = all_sequences_parsed(elems);
                let speclib = Self::apply_decoy_strategy(speclib, decoy_strategy)?;
                let Speclib::Materialized { mut elems, .. } = speclib else {
                    unreachable!("apply_decoy_strategy always returns Materialized")
                };
                let meta = apply_parse_gate(&mut elems, SeqFormat::Modified, parsable);
                let speclib = Speclib::Materialized { elems, meta };
                speclib.log_entry_stats();
                Ok(speclib)
            }
        }
    }

    /// Apply decoy strategy to an already-loaded, materialized speclib.
    ///
    /// Only ever called from the `Passthrough`/`None` branch of `from_file`
    /// (and from `from_file_with_format`'s native-format path, which is
    /// entirely untouched by the lazy-arena cutover) — both of which produce
    /// `Speclib::Materialized`. The `Lazy` arm never reaches here.
    fn apply_decoy_strategy(
        speclib: Speclib,
        strategy: crate::models::DecoyStrategy,
    ) -> Result<Speclib, LibraryReadingError> {
        use crate::models::{
            DecoyMarking,
            DecoyStrategy,
        };

        let elems = match speclib {
            Speclib::Materialized { elems, .. } => elems,
            Speclib::Lazy(_) => {
                unreachable!("apply_decoy_strategy is only ever called on a Materialized speclib")
            }
        };

        // Separate targets from decoys
        let (targets, decoys): (Vec<_>, Vec<_>) = elems
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
                Ok(Speclib::Materialized {
                    elems: targets.into_iter().chain(decoys).collect(),
                    meta: SpeclibMeta::default(),
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

                Ok(Speclib::Materialized {
                    elems: all_entries,
                    meta: SpeclibMeta::default(),
                })
            }

            DecoyStrategy::IfMissing => {
                if has_decoys {
                    tracing::info!(
                        "Library contains {} targets and {} decoys. Using existing decoys.",
                        targets.len(),
                        decoys.len()
                    );
                    Ok(Speclib::Materialized {
                        elems: targets.into_iter().chain(decoys).collect(),
                        meta: SpeclibMeta::default(),
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

                    Ok(Speclib::Materialized {
                        elems: all_entries,
                        meta: SpeclibMeta::default(),
                    })
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

        // Gate verdict on targets before the decoy expansion (see from_file).
        let parsable = all_sequences_parsed(&elems);
        let speclib = Self::apply_decoy_strategy(
            Speclib::Materialized {
                elems,
                meta: SpeclibMeta::default(),
            },
            decoy_strategy,
        )?;
        let Speclib::Materialized { mut elems, .. } = speclib else {
            unreachable!("apply_decoy_strategy always returns Materialized")
        };
        let meta = apply_parse_gate(&mut elems, SeqFormat::Modified, parsable);
        let speclib = Speclib::Materialized { elems, meta };
        speclib.log_entry_stats();
        Ok(speclib)
    }

    pub fn sample() -> Self {
        Speclib::Materialized {
            elems: vec![QueryItemToScore::sample()],
            meta: SpeclibMeta::default(),
        }
    }

    /// Materialized-arm-only accessor. The `Lazy` arm has no materialized
    /// slice by design (that is the memory optimization) — callers on that
    /// path must move to `item_at`/`iter` (Task 9/10 wires up the CLI/scoring
    /// call sites that still call `as_slice()` on a lazily-loaded lib).
    pub fn as_slice(&self) -> &[QueryItemToScore] {
        match self {
            Speclib::Materialized { elems, .. } => elems,
            Speclib::Lazy(_) => unimplemented!(
                "lazy speclib has no materialized slice; use item_at/iter (Task 9/10)"
            ),
        }
    }

    /// Log a summary of per-entry fragment and precursor counts.
    ///
    /// Entries whose fragment count exceeds [`INLINE_FRAG_CAPACITY`] spill out
    /// of the inline storage for `ExpectedIntensities` and incur a heap
    /// allocation per clone, so the spillover count is worth surfacing.
    pub fn log_entry_stats(&self) {
        use crate::models::query_item::{
            INLINE_FRAG_CAPACITY,
            INLINE_PREC_CAPACITY,
        };

        let elems = match self {
            Speclib::Materialized { elems, .. } => elems,
            Speclib::Lazy(lib) => {
                // The arena has no per-entry `ExpectedIntensities` to sample
                // (that inline-capacity/spillover concept is specific to the
                // materialized layout) — report the shape we do have.
                tracing::info!(
                    "Speclib stats: lazy arena, {} targets ({} flat scoring entries, {} total fragment slots)",
                    lib.geom.n_targets(),
                    lib.len(),
                    lib.geom.frag_labels.len(),
                );
                return;
            }
        };

        let n = elems.len();
        if n == 0 {
            tracing::info!("Speclib stats: empty library");
            return;
        }

        let mut frag_max = 0usize;
        let mut prec_max = 0usize;
        let mut frag_sum = 0usize;
        let mut frag_spill = 0usize;
        let mut prec_spill = 0usize;

        // Strided sample — a full walk is a multi-GB, seconds-long pass at peak
        // RSS for numbers that are only logged (see `avg_fragments`).
        let stride = (n / STATS_SAMPLE_TARGET).max(1);
        let mut sampled = 0usize;
        let mut i = 0;
        while i < n {
            let item = &elems[i];
            let nf = item.expected_intensity.fragment_len();
            let np = item.expected_intensity.precursor_len();
            frag_sum += nf;
            frag_max = frag_max.max(nf);
            prec_max = prec_max.max(np);
            if nf > INLINE_FRAG_CAPACITY {
                frag_spill += 1;
            }
            if np > INLINE_PREC_CAPACITY {
                prec_spill += 1;
            }
            sampled += 1;
            i += stride;
        }

        // Scale sampled spill counts back to the full library for reporting.
        let frag_spill = frag_spill * stride;
        let prec_spill = prec_spill * stride;
        let frag_mean = frag_sum as f64 / sampled as f64;

        tracing::info!(
            "Speclib stats: {} entries; fragments per entry: mean={:.1} max={} (inline cap={}, spillover entries={}); precursors per entry: max={} (inline cap={}, spillover entries={})",
            n,
            frag_mean,
            frag_max,
            INLINE_FRAG_CAPACITY,
            frag_spill,
            prec_max,
            INLINE_PREC_CAPACITY,
            prec_spill,
        );

        if frag_spill > 0 {
            let pct = (frag_spill as f64 / n as f64) * 100.0;
            tracing::warn!(
                "{}/{} speclib entries ({:.2}%) exceed INLINE_FRAG_CAPACITY={} — those entries will heap-allocate on every clone. Consider raising the inline cap in `models::query_item` if this is a large fraction.",
                frag_spill,
                n,
                pct,
                INLINE_FRAG_CAPACITY,
            );
        }
    }
}

pub struct SpeclibWriter<W: std::io::Write> {
    inner: SpeclibWriterInner<W>,
}

enum SpeclibWriterInner<W: std::io::Write> {
    MsgpackZstd(zstd::Encoder<'static, W>),
}

impl<W: std::io::Write> SpeclibWriter<W> {
    pub fn new_msgpack_zstd(writer: W) -> Result<Self, std::io::Error> {
        let encoder = zstd::Encoder::new(writer, 3)?;
        Ok(Self {
            inner: SpeclibWriterInner::MsgpackZstd(encoder),
        })
    }

    pub fn append(&mut self, elem: &SerSpeclibElement) -> Result<(), LibraryReadingError> {
        match &mut self.inner {
            SpeclibWriterInner::MsgpackZstd(encoder) => {
                rmp_serde::encode::write(encoder, elem).map_err(|e| {
                    LibraryReadingError::SpeclibParsingError {
                        source: serde_json::Error::io(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            e,
                        )),
                        context: "Error writing MessagePack",
                    }
                })?;
            }
        }
        Ok(())
    }

    pub fn finish(self) -> Result<W, std::io::Error> {
        match self.inner {
            SpeclibWriterInner::MsgpackZstd(encoder) => encoder.finish(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_sources::reference_library::ExpectedIntensity;

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

            let speclib: Vec<QueryItemToScore> = speclib_ser
                .into_iter()
                .map(QueryItemToScore::try_from)
                .collect::<Result<_, _>>()?;

            Ok(Self::Materialized {
                elems: speclib,
                meta: SpeclibMeta::default(),
            })
        }
    }

    #[test]
    fn test_detect_native_format_by_extension() {
        use std::path::Path;
        // Both .zst and .zstd must map to the native zstd readers.
        for ext in ["lib.msgpack.zst", "lib.msgpack.zstd"] {
            assert!(matches!(
                SpeclibFormat::detect_from_extension(Path::new(ext)),
                Some(SpeclibFormat::MessagePackZstd)
            ));
        }
        for ext in ["lib.ndjson.zst", "lib.ndjson.zstd"] {
            assert!(matches!(
                SpeclibFormat::detect_from_extension(Path::new(ext)),
                Some(SpeclibFormat::NdJsonZstd)
            ));
        }
        assert!(matches!(
            SpeclibFormat::detect_from_extension(Path::new("lib.msgpack")),
            Some(SpeclibFormat::MessagePack)
        ));
        assert!(matches!(
            SpeclibFormat::detect_from_extension(Path::new("lib.ndjson")),
            Some(SpeclibFormat::NdJson)
        ));
        // A .speclib must NOT be claimed as native -> routes to the bridge.
        assert!(SpeclibFormat::detect_from_extension(Path::new("lib.speclib")).is_none());
        assert!(SpeclibFormat::detect_from_extension(Path::new("lib.tsv")).is_none());
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
        assert_eq!(speclib.as_slice().len(), 1);
        println!("{:?}", speclib);

        assert_eq!(speclib.as_slice()[0].digest.decoy, DecoyMarking::Target);
        assert_eq!(speclib.as_slice()[0].digest.len(), "PEPTIDEPINK".len());
        assert_eq!(speclib.as_slice()[0].query.fragment_count(), 3);
    }

    #[test]
    fn test_sample_elution_group_deserialization() {
        let json = SerSpeclibElement::sample_json();
        let elem: SerSpeclibElement = serde_json::from_str(json).unwrap();
        let query_item: QueryItemToScore = elem.try_into().unwrap();

        assert_eq!(query_item.digest.len(), "PEPTIDEPINK".len());
        assert_eq!(query_item.query.fragment_count(), 3);
    }

    /// Every fixture below has no shipped decoys, so under the default CLI
    /// `DecoyStrategy::IfMissing` they all now resolve to `LazyMassShift`
    /// (see `map_decoy_strategy`) and load as `Speclib::Lazy`, not the old
    /// materialized target+2-decoy `Vec<QueryItemToScore>`. That routing is
    /// the whole point of the cutover, so these tests assert against the
    /// `ReferenceLibrary` arena (`RefQuery`/`QueryGeom`) instead of
    /// `Speclib::as_slice()`.
    fn expect_lazy(speclib: &Speclib) -> &ReferenceLibrary {
        match speclib {
            Speclib::Lazy(lib) => lib,
            Speclib::Materialized { .. } => {
                panic!("expected a Lazy speclib: fixture has no shipped decoys under IfMissing")
            }
        }
    }

    #[test]
    fn test_load_diann_tsv_library() {
        use timsquery::traits::QueryGeom;

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
        // Should generate 3x flat entries: 2 targets + 4 mass-shift decoys
        assert_eq!(
            speclib.len(),
            6,
            "Expected 6 entries (2 targets + 4 decoys)"
        );

        let lib = expect_lazy(&speclib);

        // Verify first target entry structure (variant 0 == target)
        let first_target = lib
            .iter()
            .find(|q| q.geom().variant() == 0)
            .expect("Should have at least one target");
        assert_eq!(
            first_target.fragment_count(),
            5,
            "First entry should have 5 fragments"
        );

        // Verify isotope envelope was added (should have 3 isotopes: 0, 1, 2)
        let envelope = first_target.expected_precursor_envelope();
        assert_eq!(
            envelope.len(),
            3,
            "Should have 3 isotopes in precursor envelope"
        );
        assert_eq!(envelope[0].0, 0i8, "isotope 0 present");
        assert_eq!(envelope[1].0, 1i8, "isotope 1 present");
        assert_eq!(envelope[2].0, 2i8, "isotope 2 present");
    }

    #[test]
    fn test_diann_tsv_parsable_gate() {
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        assert!(
            speclib.parsable_sequences(),
            "DIA-NN sample fixture should all parse with modified_peptide"
        );
    }

    #[test]
    fn test_load_skyline_csv_library() {
        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("skyline_io_files")
            .join("sample_transition_list.csv");

        assert!(
            test_file.exists(),
            "Test file should exist at {:?}",
            test_file
        );

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load Skyline CSV library");

        // Skyline is not the DIA-NN collection shape `build_reference_library`
        // ingests (see the `is_diann_shape` gate in `Speclib::from_file`), so
        // it stays on the materialized path even though the resolved decoy
        // strategy is `LazyMassShift` (no shipped decoys, default IfMissing).
        // Fixture has 14 PRTC targets, no decoys -> 14 targets + 28 mass-shift decoys
        assert_eq!(
            speclib.len(),
            42,
            "Expected 42 entries (14 targets + 28 decoys)"
        );

        let n_targets = speclib
            .as_slice()
            .iter()
            .filter(|e| !e.digest.is_decoy())
            .count();
        let n_decoys = speclib
            .as_slice()
            .iter()
            .filter(|e| e.digest.is_decoy())
            .count();
        assert_eq!(n_targets, 14, "Should have 14 targets");
        assert_eq!(n_decoys, 28, "Should have 28 decoys");

        // Isotope envelope should have been attached (3 isotopes) for every target
        for entry in speclib.as_slice().iter().filter(|e| !e.digest.is_decoy()) {
            assert_eq!(
                entry.query.iter_precursors().count(),
                3,
                "Each target should have 3 isotopes in precursor envelope"
            );
            assert!(
                entry.query.fragment_count() > 0,
                "Each target should have at least one fragment"
            );
        }
    }

    #[test]
    fn test_load_diann_txt_library() {
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

        let lib = expect_lazy(&speclib);
        let n_targets = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();

        assert_eq!(n_targets, 2, "Should have 2 targets");
        assert_eq!(n_decoys, 4, "Should have 4 decoys");
    }

    #[test]
    fn test_load_diann_parquet_library() {
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

        let lib = expect_lazy(&speclib);
        let n_targets = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();

        assert_eq!(n_targets, 3, "Should have 3 targets");
        assert_eq!(n_decoys, 6, "Should have 6 decoys");

        // Verify isotope envelope for targets
        for q in lib.iter().filter(|q| q.geom().variant() == 0) {
            assert_eq!(
                q.expected_precursor_envelope().len(),
                3,
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

        let lib = expect_lazy(&speclib);

        // Check that isotope intensities are normalized (M0 should be 1.0),
        // for every flat entry (targets AND decoy variants — the envelope is
        // computed per-target and shared across variants, see
        // `decoy_variant_reuses_target_intensities` in reference_library.rs).
        for q in lib.iter() {
            let envelope = q.expected_precursor_envelope();
            let m0 = envelope[0].1;
            assert!(
                (m0 - 1.0).abs() < 0.01,
                "M0 should be normalized to ~1.0, got {}",
                m0
            );
            for &(_iso, intensity) in envelope.iter().skip(1) {
                assert!(
                    intensity <= 1.0,
                    "M+n intensity should be <= 1.0, got {}",
                    intensity
                );
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

        let lib = expect_lazy(&speclib);
        let n_targets = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();

        assert_eq!(n_targets, 2, "Should have 2 target entries");
        assert_eq!(n_decoys, 4, "Should have 4 decoy entries (2 per target)");

        // Each target index (== decoy_group in the old materialized scheme)
        // should have exactly 3 flat variants (1 target + 2 decoys) — this is
        // guaranteed structurally by `ReferenceLibrary::item_at`'s `t,+,-`
        // packing, so assert it directly rather than re-deriving groups.
        let n_target_indices = lib.geom.n_targets();
        assert_eq!(n_target_indices, 2, "Should have 2 unique targets");
        for tgt in 0..n_target_indices as u32 {
            let variants: Vec<u8> = (0..3)
                .map(|v| RefQuery::new(lib, tgt, v).geom().variant())
                .collect();
            assert_eq!(
                variants,
                vec![0, 1, 2],
                "target {tgt} should have exactly 1 target + 2 decoy variants"
            );
        }
    }

    #[test]
    fn test_mass_shift_decoys() {
        use timsquery::traits::QueryGeom;

        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        let lib = expect_lazy(&speclib);

        // Unified CH2 offset (see `map_decoy_strategy`), replacing the old
        // 12.0 (materialized `IfMissing`) / 14.0 (materialized `Force`) split.
        const CH2_MASS: f64 = 14.0;

        for tgt in 0..lib.geom.n_targets() as u32 {
            let target = RefQuery::new(lib, tgt, 0);
            let plus = RefQuery::new(lib, tgt, 1);
            let minus = RefQuery::new(lib, tgt, 2);

            let target_mz = target.mono_precursor_mz();
            let charge = target.precursor_charge();
            let mz_shift = CH2_MASS / charge as f64;

            assert!(
                (plus.mono_precursor_mz() - (target_mz + mz_shift)).abs() < 0.001,
                "+ decoy should be shifted by +CH2/charge"
            );
            assert!(
                (minus.mono_precursor_mz() - (target_mz - mz_shift)).abs() < 0.001,
                "- decoy should be shifted by -CH2/charge"
            );

            // Verify other properties are preserved
            for decoy in [&plus, &minus] {
                assert_eq!(
                    decoy.rt_seconds(),
                    target.rt_seconds(),
                    "RT should be preserved"
                );
                assert_eq!(
                    decoy.mobility_ook0(),
                    target.mobility_ook0(),
                    "Mobility should be preserved"
                );
                assert_eq!(
                    decoy.fragment_count(),
                    target.fragment_count(),
                    "Fragment count should be preserved"
                );
            }
        }
    }

    #[test]
    fn test_fragment_intensities_preserved() {
        use timsquery::traits::QueryGeom;

        let test_file = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("timsquery")
            .join("tests")
            .join("diann_io_files")
            .join("sample_lib.txt");

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyStrategy::default())
            .expect("Failed to load DIA-NN TSV library");

        let lib = expect_lazy(&speclib);
        for q in lib.iter() {
            let fragments: Vec<_> = q.iter_expected_fragments().collect();
            assert_eq!(
                fragments.len(),
                q.fragment_count(),
                "Fragment intensity count should match fragment count"
            );
            for (_label, intensity) in fragments {
                assert!(intensity > 0.0, "Fragment intensities should be positive");
            }
        }
    }

    #[test]
    fn test_ser_speclib_element_roundtrip() {
        let elem = SerSpeclibElement::new(
            PrecursorEntry::new("PEPTIDEK".to_string(), 2, false, 0),
            ReferenceEG::new(
                0,
                500.0,
                vec![0, 1, 2],
                vec![300.0, 400.0],
                vec![
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                ],
                vec![1.0, 0.5, 0.2],
                vec![0.8, 0.3],
                0.75,
                120.0,
            ),
        );
        let bytes = rmp_serde::to_vec(&elem).unwrap();
        let decoded: SerSpeclibElement = rmp_serde::from_slice(&bytes).unwrap();
        let qi: QueryItemToScore = decoded.try_into().unwrap();
        assert_eq!(qi.digest.len(), "PEPTIDEK".len());
        assert_eq!(qi.query.fragment_count(), 2);
    }

    #[test]
    fn test_speclib_writer_roundtrip() {
        let elem = SerSpeclibElement::sample();
        let mut buf = Vec::new();
        {
            let mut writer = SpeclibWriter::new_msgpack_zstd(&mut buf).unwrap();
            writer.append(&elem).unwrap();
            writer.append(&elem).unwrap();
            writer.finish().unwrap();
        }
        let reader =
            SpeclibReader::new(std::io::Cursor::new(&buf), SpeclibFormat::MessagePackZstd).unwrap();
        let items: Vec<_> = reader.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(items.len(), 2);
    }

    #[cfg(test)]
    fn fixture_query_item(
        raw_seq: &str,
        parsed: Option<crate::models::sequence::ParsedSequence>,
    ) -> QueryItemToScore {
        let elem = SerSpeclibElement::new(
            PrecursorEntry::new(raw_seq.to_string(), 2, false, 0),
            ReferenceEG::new(
                0,
                500.0,
                vec![0, 1, 2],
                vec![300.0, 400.0],
                vec![
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                ],
                vec![1.0, 0.5, 0.2],
                vec![0.8, 0.3],
                0.75,
                120.0,
            ),
        );
        let mut q: QueryItemToScore = elem.try_into().expect("fixture convert");
        q.digest.parsed = parsed;
        q
    }

    #[test]
    fn test_parse_gate_off_on_poisoned_row() {
        use crate::models::sequence::{
            ParsedSequence,
            SeqFormat,
        };
        use smallvec::SmallVec;

        let mut items = vec![
            fixture_query_item(
                "PEPTIDEK",
                Some(ParsedSequence {
                    residues: SmallVec::new(),
                    mods: SmallVec::new(),
                }),
            ),
            fixture_query_item("GARBAGE!!!", None),
        ];
        let parsable = all_sequences_parsed(&items);
        let meta = apply_parse_gate(&mut items, SeqFormat::Modified, parsable);
        assert!(!meta.parsable_sequences);
        assert!(items.iter().all(|q| q.digest.parsed.is_none()));
    }

    #[test]
    fn test_parse_gate_on_when_all_parsed() {
        use crate::models::sequence::{
            ParsedSequence,
            SeqFormat,
        };
        use smallvec::SmallVec;

        let mut items = vec![
            fixture_query_item(
                "PEPTIDEK",
                Some(ParsedSequence {
                    residues: SmallVec::new(),
                    mods: SmallVec::new(),
                }),
            ),
            fixture_query_item(
                "ABCDEK",
                Some(ParsedSequence {
                    residues: SmallVec::new(),
                    mods: SmallVec::new(),
                }),
            ),
        ];
        let parsable = all_sequences_parsed(&items);
        let meta = apply_parse_gate(&mut items, SeqFormat::Modified, parsable);
        assert!(meta.parsable_sequences);
        assert!(items.iter().all(|q| q.digest.parsed.is_some()));
    }

    /// One target, one DIA-NN extras row -> the arena should have exactly
    /// one target (3 flat entries: target + 2 lazy decoys), fragment m/z
    /// should round-trip through the CSR arena, and `frag_intens` should be
    /// parallel to `frag_labels` in the same order.
    #[test]
    fn build_reference_library_from_diann_collection() {
        use timsquery::traits::QueryGeom;

        let eg = TimsElutionGroup::<IonAnnot>::builder()
            .id(0)
            .mobility_ook0(1.0)
            .rt_seconds(100.0)
            .precursor(500.25, 2)
            .precursor_labels([0i8].as_slice().into())
            .fragment_labels(
                [
                    IonAnnot::try_from("y3").unwrap(),
                    IonAnnot::try_from("y5").unwrap(),
                ]
                .as_slice()
                .into(),
            )
            .fragment_mzs(vec![300.1, 500.2])
            .try_build()
            .expect("build fixture elution group");

        let extra = DiannPrecursorExtras {
            modified_peptide: "PEPTIDEK".to_string(),
            stripped_peptide: "PEPTIDEK".to_string(),
            protein_id: "PROT1".to_string(),
            is_decoy: false,
            relative_intensities: vec![
                (IonAnnot::try_from("y3").unwrap(), 1.0),
                (IonAnnot::try_from("y5").unwrap(), 0.5),
            ],
        };

        let collection = ElutionGroupCollection::MzpafLabels(
            vec![eg],
            Some(FileReadingExtras::Diann(vec![extra])),
        );

        let (lib, report) =
            build_reference_library(collection, LibCapabilities::default_diann()).expect("build");

        assert_eq!(lib.geom.n_targets(), 1);
        assert_eq!(lib.len(), 3, "1 target x (1 target + 2 decoys) variants");
        assert_eq!(lib.frag_intens.len(), lib.geom.frag_labels.len());
        assert_eq!(report.n_targets, 1);
        assert_eq!(report.n_averagine_fallback, 0);
        assert_eq!(report.sequence_features, SeqFeatureState::Available);

        // Fragment m/z round-trips through the arena for the target variant
        // (variant 0 is unshifted).
        let target = lib.item_at(0);
        let frags: Vec<_> = target
            .iter_fragments_refs()
            .map(|(l, mz)| (*l, mz))
            .collect();
        assert_eq!(frags.len(), 2);
        assert_eq!(frags[0].0, IonAnnot::try_from("y3").unwrap());
        assert!((frags[0].1 - 300.1).abs() < 1e-9);
        assert_eq!(frags[1].0, IonAnnot::try_from("y5").unwrap());
        assert!((frags[1].1 - 500.2).abs() < 1e-9);

        // frag_intens is parallel to frag_labels, same order as extras.
        assert!((lib.frag_intens[0] - 1.0).abs() < 1e-6);
        assert!((lib.frag_intens[1] - 0.5).abs() < 1e-6);
    }
}
