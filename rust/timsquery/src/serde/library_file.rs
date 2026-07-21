pub use super::diann_io::DiannPrecursorExtras;
use super::diann_io::{
    read_library_file as read_diann_tsv,
    read_parquet_library_file as read_diann_parquet,
    sniff_diann_library_file,
    sniff_diann_parquet_library_file,
};
use super::diann_speclib_io::{
    read_diann_speclib_library_file,
    sniff_diann_speclib_library_file,
};
use super::elution_group_inputs::{
    ElutionGroupInput,
    ElutionGroupInputError,
};
pub use super::skyline_io::SkylinePrecursorExtras;
use super::skyline_io::{
    read_library_file as read_skyline_csv,
    sniff_skyline_library_file,
};
pub use super::spectronaut_io::SpectronautPrecursorExtras;
use super::spectronaut_io::{
    read_library_file as read_spectronaut_tsv,
    sniff_spectronaut_library_file,
};
use crate::TimsElutionGroup;
use crate::ion::IonAnnot;
use crate::models::{
    LibCapabilities,
    QueryCollection,
};
use std::path::Path;
use std::sync::Arc;
use tracing::{
    debug,
    info,
    warn,
};

#[derive(Debug)]
pub enum LibraryReadingError {
    IoError(std::io::Error),
    SerdeJsonError(serde_json::Error),
    ElutionGroupInputError(ElutionGroupInputError),
    UnableToParseElutionGroups,
    /// A `.speclib` whose version is newer (more negative) than this reader
    /// supports.
    UnsupportedSpeclibVersion(i32),
    /// A `.speclib` structural desync / bad byte count / unexpected value, with
    /// context.
    SpeclibParse(String),
}

impl From<serde_json::Error> for LibraryReadingError {
    fn from(err: serde_json::Error) -> Self {
        LibraryReadingError::SerdeJsonError(err)
    }
}
impl From<ElutionGroupInputError> for LibraryReadingError {
    fn from(err: ElutionGroupInputError) -> Self {
        LibraryReadingError::ElutionGroupInputError(err)
    }
}

#[derive(Debug)]
pub enum FileReadingExtras {
    Diann(Vec<DiannPrecursorExtras>),
    Spectronaut(Vec<SpectronautPrecursorExtras>),
    Skyline(Vec<SkylinePrecursorExtras>),
}

#[derive(Debug)]
pub enum ElutionGroupCollection {
    StringLabels(Vec<TimsElutionGroup<String>>, Option<FileReadingExtras>),
    MzpafLabels(Vec<TimsElutionGroup<IonAnnot>>, Option<FileReadingExtras>),
    TinyIntLabels(Vec<TimsElutionGroup<u8>>, Option<FileReadingExtras>),
    IntLabels(Vec<TimsElutionGroup<u32>>, Option<FileReadingExtras>),
}

impl ElutionGroupCollection {
    pub fn len(&self) -> usize {
        match self {
            ElutionGroupCollection::StringLabels(egs, _) => egs.len(),
            ElutionGroupCollection::MzpafLabels(egs, _) => egs.len(),
            ElutionGroupCollection::TinyIntLabels(egs, _) => egs.len(),
            ElutionGroupCollection::IntLabels(egs, _) => egs.len(),
        }
    }

    fn try_read_json(path: &Path) -> Result<Self, LibraryReadingError> {
        let file_content = std::fs::read_to_string(path).map_err(LibraryReadingError::IoError)?;
        info!("Read file content from {}", path.display());
        // First try direct deserialization
        if let Ok(egs) = Self::try_deser_direct(&file_content) {
            info!("Successfully deserialized elution groups directly");
            return Ok(egs);
        }
        // Next try deserialization via inputed format
        match Self::try_deser_inputed(&file_content) {
            Ok(egs) => {
                info!("Successfully deserialized elution groups via inputed format");
                Ok(egs)
            }
            Err(_) => Err(LibraryReadingError::UnableToParseElutionGroups),
        }
    }

    fn try_deser_inputed(content: &str) -> Result<Self, LibraryReadingError> {
        // We can try from smallest to largest overhead
        // Here we try to do the deser into ElutionGroupInput variants first
        debug!("Attempting deserialization of elution group inputs");
        debug!("Attempting to deserialize elution group inputs with tiny int labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<u8>>>(content) {
            // Here we can handle filling the inputs if they are needed...
            let eg_inputs = if eg_inputs.first().is_some_and(|x| x.needs_fragment_labels()) {
                debug!("Filling missing fragment labels with tiny int labels");
                eg_inputs
                    .into_iter()
                    .map(|x| x.try_fill_labels_u8())
                    .collect::<Result<_, _>>()
            } else {
                Ok(eg_inputs)
            };

            let out: Result<Vec<TimsElutionGroup<u8>>, ElutionGroupInputError> = eg_inputs?
                .into_iter()
                .map(<ElutionGroupInput<u8> as TryInto<TimsElutionGroup<u8>>>::try_into)
                .collect();
            return Ok(ElutionGroupCollection::TinyIntLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with int labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<u32>>>(content) {
            let out: Result<Vec<TimsElutionGroup<u32>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::IntLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with mzpaf labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<IonAnnot>>>(content) {
            let out: Result<Vec<TimsElutionGroup<IonAnnot>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::MzpafLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with string labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<String>>>(content) {
            let out: Result<Vec<TimsElutionGroup<String>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::StringLabels(out?, None));
        }
        Err(LibraryReadingError::UnableToParseElutionGroups)
    }

    fn try_deser_direct(content: &str) -> Result<Self, LibraryReadingError> {
        // We can try from smallest to largest overhead
        // Here we try to do the direct deser into ElutionGroupCollection variants
        // u8 -> u32 -> IonAnnot -> String
        debug!("Attempting direct deserialization of elution groups");
        debug!("Attempting to deserialize elution groups with tiny int labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<u8>>>(content) {
            return Ok(ElutionGroupCollection::TinyIntLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with int labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<u32>>>(content) {
            return Ok(ElutionGroupCollection::IntLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with mzpaf labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<IonAnnot>>>(content) {
            return Ok(ElutionGroupCollection::MzpafLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with string labels");
        if let Ok(egs) = serde_json::from_str::<Vec<TimsElutionGroup<String>>>(content) {
            return Ok(ElutionGroupCollection::StringLabels(egs, None));
        }
        Err(LibraryReadingError::UnableToParseElutionGroups)
    }
}

/// The label-typed columnar library store returned by [`read_library_file`].
///
/// One funnel: every format lands in exactly one variant. DIA-NN family formats
/// (`.speclib`/TSV/parquet) carry ion-chemistry (`IonAnnot`) labels and land in
/// [`LibraryArena::Mzpaf`]; string-labelled JSON lands in [`LibraryArena::Str`].
///
/// `frag_intens` is the reference-intensity sidecar, parallel to
/// `geom.frag_labels`/`geom.frag_mzs` (same length). The columnar store itself
/// stays intensity-free (extraction does not need intensities); the DIA-NN
/// family readers (`.speclib`/TSV/parquet) populate the sidecar (`Some`) so the
/// timsseek bridge can zip in reference intensities, while the intensity-free
/// mzpaf path and string-labelled JSON leave it `None`. Extraction (cli)
/// ignores it.
pub enum LibraryArena {
    Mzpaf {
        geom: QueryCollection<IonAnnot>,
        frag_intens: Option<Vec<f32>>,
    },
    Str {
        geom: QueryCollection<Arc<str>>,
    },
}

/// The reader-extras fields the arena consumes, flattened out of the
/// otherwise-identical DIA-NN/Skyline/Spectronaut `*PrecursorExtras` structs.
struct PrecursorExtrasRow {
    modified: String,
    stripped: String,
    is_decoy: bool,
    relative_intensities: Vec<(IonAnnot, f32)>,
}

impl From<DiannPrecursorExtras> for PrecursorExtrasRow {
    fn from(e: DiannPrecursorExtras) -> Self {
        Self {
            modified: e.modified_peptide,
            stripped: e.stripped_peptide,
            is_decoy: e.is_decoy,
            relative_intensities: e.relative_intensities,
        }
    }
}

impl From<SkylinePrecursorExtras> for PrecursorExtrasRow {
    fn from(e: SkylinePrecursorExtras) -> Self {
        Self {
            modified: e.modified_peptide,
            stripped: e.stripped_peptide,
            is_decoy: e.is_decoy,
            relative_intensities: e.relative_intensities,
        }
    }
}

impl From<SpectronautPrecursorExtras> for PrecursorExtrasRow {
    fn from(e: SpectronautPrecursorExtras) -> Self {
        Self {
            modified: e.modified_peptide,
            stripped: e.stripped_peptide,
            is_decoy: e.is_decoy,
            relative_intensities: e.relative_intensities,
        }
    }
}

impl LibraryArena {
    /// Build an `Mzpaf` arena WITH the reference-intensity sidecar from
    /// `IonAnnot`-labelled groups plus their reader extras.
    ///
    /// Reference intensities are matched BY LABEL: the extras carry
    /// `(IonAnnot, intensity)` pairs, and for every fragment pushed into the
    /// arena (in `iter_fragments` order, which is the order that lands in
    /// `geom.frag_labels`) the intensity is looked up by that fragment's label.
    /// A fragment with no matching reference intensity fails the whole load
    /// loudly (never a silent zero). The resulting `frag_intens` is parallel to
    /// `geom.frag_labels`. Sequences and the per-row decoy flag are threaded in
    /// too, so `seal()` sees shipped decoys and the timsseek parse gate sees the
    /// modified sequence.
    fn mzpaf_with_intensities(
        egs: Vec<TimsElutionGroup<IonAnnot>>,
        extras: FileReadingExtras,
    ) -> Result<Self, LibraryReadingError> {
        let rows: Vec<PrecursorExtrasRow> = match extras {
            FileReadingExtras::Diann(v) => v.into_iter().map(PrecursorExtrasRow::from).collect(),
            FileReadingExtras::Skyline(v) => v.into_iter().map(PrecursorExtrasRow::from).collect(),
            FileReadingExtras::Spectronaut(v) => {
                v.into_iter().map(PrecursorExtrasRow::from).collect()
            }
        };

        if egs.len() != rows.len() {
            return Err(LibraryReadingError::SpeclibParse(format!(
                "elution groups ({}) and reader extras ({}) length mismatch",
                egs.len(),
                rows.len()
            )));
        }

        let mut geom =
            QueryCollection::with_capabilities(LibCapabilities::default_diann_no_decoys());
        let mut frag_intens: Vec<f32> = Vec::new();

        for (eg, row) in egs.iter().zip(rows) {
            // Reference intensities keyed by fragment label (see fn docs).
            let lookup: std::collections::HashMap<IonAnnot, f32> =
                row.relative_intensities.into_iter().collect();
            let frags: Vec<(IonAnnot, f64)> = eg.iter_fragments().map(|(l, mz)| (*l, mz)).collect();
            for (label, _) in &frags {
                let intensity = lookup.get(label).ok_or_else(|| {
                    LibraryReadingError::SpeclibParse(format!(
                        "fragment {label:?} of precursor {:?} has no reference intensity",
                        row.modified
                    ))
                })?;
                frag_intens.push(*intensity);
            }
            geom.push_row(
                eg.precursor_mz(),
                eg.precursor_charge(),
                eg.rt_seconds(),
                eg.mobility_ook0(),
                &frags,
                &row.stripped,
                &row.modified,
                &[],
                row.is_decoy,
            );
        }

        if frag_intens.len() != geom.frag_labels.len() {
            return Err(LibraryReadingError::SpeclibParse(format!(
                "reference-intensity sidecar ({}) must stay parallel to the fragment-label arena ({})",
                frag_intens.len(),
                geom.frag_labels.len(),
            )));
        }

        geom.seal();
        Ok(LibraryArena::Mzpaf {
            geom,
            frag_intens: Some(frag_intens),
        })
    }

    /// Adapt the legacy [`ElutionGroupCollection`] (produced by the non-speclib
    /// readers) into the arena. `IonAnnot`-labelled groups become
    /// [`LibraryArena::Mzpaf`]; string-labelled groups become
    /// [`LibraryArena::Str`]. Integer-labelled groups have no arena variant (no
    /// live consumer) and are rejected.
    ///
    /// When the reader supplied per-precursor extras (DIA-NN/Skyline/
    /// Spectronaut all carry `relative_intensities`, `modified`/`stripped`
    /// sequences and an `is_decoy` flag), the reference intensities are threaded
    /// into the `frag_intens` sidecar so the timsseek bridge can score against
    /// them; sequences and the decoy flag are threaded into the arena too. When
    /// no extras were supplied (e.g. plain `IonAnnot` JSON, which only exercises
    /// the extraction/geometry path), `frag_intens` stays `None` — matching the
    /// historical behavior where timsseek rejected that shape.
    fn from_elution_groups(egc: ElutionGroupCollection) -> Result<Self, LibraryReadingError> {
        match egc {
            ElutionGroupCollection::MzpafLabels(egs, Some(extras)) => {
                Self::mzpaf_with_intensities(egs, extras)
            }
            ElutionGroupCollection::MzpafLabels(egs, None) => {
                let mut geom =
                    QueryCollection::with_capabilities(LibCapabilities::default_diann_no_decoys());
                for eg in &egs {
                    let frags: Vec<(IonAnnot, f64)> =
                        eg.iter_fragments().map(|(l, mz)| (*l, mz)).collect();
                    geom.push_target(
                        eg.precursor_mz(),
                        eg.precursor_charge(),
                        eg.rt_seconds(),
                        eg.mobility_ook0(),
                        &frags,
                        "",
                        "",
                        &[],
                    );
                }
                geom.seal();
                Ok(LibraryArena::Mzpaf {
                    geom,
                    frag_intens: None,
                })
            }
            ElutionGroupCollection::StringLabels(egs, _) => {
                // String-labelled arenas carry no ion chemistry and ship no
                // decoys: sequence/fragment features unavailable, decoys off.
                let mut geom =
                    QueryCollection::with_capabilities(LibCapabilities::default_unlabeled());
                for eg in &egs {
                    let frags: Vec<(Arc<str>, f64)> = eg
                        .iter_fragments()
                        .map(|(l, mz)| (Arc::<str>::from(l.as_str()), mz))
                        .collect();
                    geom.push_target(
                        eg.precursor_mz(),
                        eg.precursor_charge(),
                        eg.rt_seconds(),
                        eg.mobility_ook0(),
                        &frags,
                        "",
                        "",
                        &[],
                    );
                }
                geom.seal();
                Ok(LibraryArena::Str { geom })
            }
            ElutionGroupCollection::TinyIntLabels(..) | ElutionGroupCollection::IntLabels(..) => {
                warn!("integer-labelled libraries have no LibraryArena variant; rejecting");
                Err(LibraryReadingError::UnableToParseElutionGroups)
            }
        }
    }
}

/// A single spectral-library format reader. Adding a format = one struct + one
/// line in [`registry`], instead of editing an enum, a method, and an ordered
/// try-chain.
pub trait LibraryReader: Send + Sync {
    fn name(&self) -> &'static str;
    /// Cheap probe: header bytes / extension / first data row. Must not read the
    /// whole file.
    fn sniff(&self, path: &Path) -> bool;
    fn read(&self, path: &Path) -> Result<ElutionGroupCollection, LibraryReadingError>;
}

struct DiannParquetReader;
struct DiannTsvReader;
struct SpectronautReader;
struct SkylineReader;
struct JsonReader;

impl LibraryReader for DiannParquetReader {
    fn name(&self) -> &'static str {
        "diann-parquet"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_parquet_library_file(path)
    }

    fn read(&self, path: &Path) -> Result<ElutionGroupCollection, LibraryReadingError> {
        let egs = read_diann_parquet(path).map_err(|e| {
            warn!("Failed to read DIA-NN parquet library file: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        Ok(ElutionGroupCollection::MzpafLabels(
            egs,
            Some(FileReadingExtras::Diann(extras)),
        ))
    }
}

impl LibraryReader for DiannTsvReader {
    fn name(&self) -> &'static str {
        "diann-tsv"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_library_file(path)
    }

    fn read(&self, path: &Path) -> Result<ElutionGroupCollection, LibraryReadingError> {
        let egs = read_diann_tsv(path).map_err(|e| {
            warn!("Failed to read DIA-NN TSV library file: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        Ok(ElutionGroupCollection::MzpafLabels(
            egs,
            Some(FileReadingExtras::Diann(extras)),
        ))
    }
}

impl LibraryReader for SpectronautReader {
    fn name(&self) -> &'static str {
        "spectronaut-tsv"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_spectronaut_library_file(path).is_ok()
    }

    fn read(&self, path: &Path) -> Result<ElutionGroupCollection, LibraryReadingError> {
        let egs = read_spectronaut_tsv(path).map_err(|e| {
            warn!("Failed to read Spectronaut TSV library file: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        Ok(ElutionGroupCollection::MzpafLabels(
            egs,
            Some(FileReadingExtras::Spectronaut(extras)),
        ))
    }
}

impl LibraryReader for SkylineReader {
    fn name(&self) -> &'static str {
        "skyline-csv"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_skyline_library_file(path).is_ok()
    }

    fn read(&self, path: &Path) -> Result<ElutionGroupCollection, LibraryReadingError> {
        let egs = read_skyline_csv(path).map_err(|e| {
            warn!("Failed to read Skyline transition list: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        Ok(ElutionGroupCollection::MzpafLabels(
            egs,
            Some(FileReadingExtras::Skyline(extras)),
        ))
    }
}

impl LibraryReader for JsonReader {
    fn name(&self) -> &'static str {
        "json"
    }

    /// Terminal fallback: always sniffs true so JSON is tried when nothing else
    /// matched. Must be last in the registry.
    fn sniff(&self, _path: &Path) -> bool {
        true
    }

    fn read(&self, path: &Path) -> Result<ElutionGroupCollection, LibraryReadingError> {
        ElutionGroupCollection::try_read_json(path)
    }
}

fn registry() -> &'static [&'static dyn LibraryReader] {
    &[
        &DiannParquetReader,
        &DiannTsvReader,
        &SpectronautReader,
        &SkylineReader,
        &JsonReader,
    ]
}

pub fn read_library_file<T: AsRef<Path>>(path: T) -> Result<LibraryArena, LibraryReadingError> {
    let path = path.as_ref();
    // The DIA-NN `.speclib` reader builds the columnar arena directly (with the
    // reference-intensity sidecar); every other format still produces the legacy
    // `ElutionGroupCollection`, adapted into the arena here. `.speclib` is
    // sniffed first because its `read` path is the only one that can surface an
    // `UnsupportedSpeclibVersion` diagnostic (the sniff has no version gate).
    if sniff_diann_speclib_library_file(path) {
        info!("Dispatching library read to diann-speclib (direct arena build)");
        return read_diann_speclib_library_file(path);
    }
    let mut last_err = None;
    for reader in registry() {
        if reader.sniff(path) {
            info!("Dispatching library read to {}", reader.name());
            match reader.read(path) {
                Ok(egs) => return LibraryArena::from_elution_groups(egs),
                // A sniff can fire on a file the reader then fails to parse
                // (overlapping sniffs). Fall through to the next candidate
                // instead of committing to the first sniff. Keep the FIRST
                // error: readers run specific -> generic, so the earliest
                // sniff-and-fail carries the most useful diagnostic (e.g. a
                // `.speclib` desync), which the always-true JSON fallback's
                // generic error would otherwise clobber.
                Err(e) => {
                    warn!("{} sniffed but failed to read: {:?}", reader.name(), e);
                    last_err.get_or_insert(e);
                }
            }
        }
    }
    // Dead default in practice (JsonReader always sniffs true) — a harmless
    // defensive fallback.
    Err(last_err.unwrap_or(LibraryReadingError::UnableToParseElutionGroups))
}
