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
use super::mzspeclib_io::{
    read_mzspeclib_library_file,
    sniff_mzspeclib_library_file,
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
}

impl ElutionGroupCollection {
    pub fn len(&self) -> usize {
        match self {
            ElutionGroupCollection::StringLabels(egs, _) => egs.len(),
            ElutionGroupCollection::MzpafLabels(egs, _) => egs.len(),
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
        // Next try deserialization via inputed format. Its error is returned
        // as-is: it names the field that is wrong, which
        // `UnableToParseElutionGroups` does not.
        let egs = Self::try_deser_inputed(&file_content)?;
        info!("Successfully deserialized elution groups via inputed format");
        Ok(egs)
    }

    fn try_deser_inputed(content: &str) -> Result<Self, LibraryReadingError> {
        // mzpaf before string: `"y1"` deserializes as either, and only the
        // mzpaf form carries ion chemistry.
        debug!("Attempting deserialization of elution group inputs");
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
        // mzpaf before string, for the same reason as `try_deser_inputed`.
        debug!("Attempting direct deserialization of elution groups");
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
        let mut n_collapsed = 0usize;

        for (eg, row) in egs.iter().zip(rows) {
            // Reference intensities keyed by fragment label (see fn docs).
            let lookup: std::collections::HashMap<IonAnnot, f32> =
                row.relative_intensities.into_iter().collect();

            // Through `FragmentSet` rather than straight into a Vec: two TSV
            // rows with the same series + ordinal + charge for one precursor
            // produce two identical labels, which panics in scoring.
            let mut set = FragmentSet::with_capacity(eg.iter_fragments().count());
            for (label, mz) in eg.iter_fragments() {
                let intensity = lookup.get(label).ok_or_else(|| {
                    LibraryReadingError::SpeclibParse(format!(
                        "fragment {label:?} of precursor {:?} has no reference intensity",
                        row.modified
                    ))
                })?;
                if set.insert(*label, mz, *intensity) == Inserted::Collapsed {
                    n_collapsed += 1;
                }
            }

            set.extend_sidecar(&mut frag_intens);
            geom.push_row(
                eg.precursor_mz(),
                eg.precursor_charge(),
                eg.rt_seconds(),
                eg.mobility_ook0(),
                set.frags(),
                &row.stripped,
                &row.modified,
                &[],
                row.is_decoy,
            );
        }

        if n_collapsed > 0 {
            warn!("{n_collapsed} fragments collapsed onto a duplicate label");
        }

        finish_mzpaf_arena(geom, frag_intens)
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
        }
    }
}

/// What [`FragmentSet::insert`] did with a fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Inserted {
    /// A label not yet in this precursor; stored.
    Added,
    /// The label was already present, so the two peaks were collapsed onto the
    /// more intense one. The caller counts this under its own name.
    Collapsed,
}

/// The fragments of one precursor, with labels unique by construction.
///
/// Fragment labels MUST be unique within a precursor. `linear_get` is
/// first-match, so a duplicate silently shadows one peak — and
/// `ExpectedIntensities::try_from_pairs` rejects duplicates outright, which
/// timsseek's scoring pipeline `.expect()`s. A duplicate reaching the arena is
/// therefore a panic mid-search, per candidate.
///
/// Every reader that builds an mzpaf arena goes through this type, so the
/// invariant holds in one place instead of being re-derived (or, in the TSV
/// readers' case, forgotten) at each site. Collisions collapse onto the more
/// intense peak rather than keeping whichever came first: intensity is the
/// signal being scored, and file order is not meaningful.
#[derive(Debug, Default)]
pub(super) struct FragmentSet {
    frags: Vec<(IonAnnot, f64)>,
    intensities: Vec<f32>,
}

impl FragmentSet {
    pub(super) fn with_capacity(n: usize) -> Self {
        Self {
            frags: Vec::with_capacity(n),
            intensities: Vec::with_capacity(n),
        }
    }

    pub(super) fn insert(&mut self, label: IonAnnot, mz: f64, intensity: f32) -> Inserted {
        if let Some(idx) = self.frags.iter().position(|(l, _)| *l == label) {
            if intensity > self.intensities[idx] {
                self.frags[idx].1 = mz;
                self.intensities[idx] = intensity;
            }
            return Inserted::Collapsed;
        }
        self.frags.push((label, mz));
        self.intensities.push(intensity);
        Inserted::Added
    }

    pub(super) fn is_empty(&self) -> bool {
        self.frags.is_empty()
    }

    /// Empty without releasing the allocation, so one set can be reused across
    /// entries (the DIA-NN `.speclib` reader keeps one per rayon worker).
    pub(super) fn clear(&mut self) {
        self.frags.clear();
        self.intensities.clear();
    }

    #[cfg(test)]
    pub(super) fn len(&self) -> usize {
        self.frags.len()
    }

    pub(super) fn frags(&self) -> &[(IonAnnot, f64)] {
        &self.frags
    }

    /// Append this precursor's intensities to a whole-library sidecar, in the
    /// same order `frags()` will be pushed.
    pub(super) fn extend_sidecar(&self, sidecar: &mut Vec<f32>) {
        sidecar.extend_from_slice(&self.intensities);
    }
}

/// Seal a directly-built mzpaf arena together with its reference-intensity
/// sidecar.
///
/// The sidecar is indexed by the same offsets as `frag_labels`, so a length
/// mismatch means some path pushed a label without an intensity (or the
/// reverse) and every downstream fragment lookup is off by that much. Returned
/// as an error rather than asserted: it is a property of the file being read.
pub(super) fn finish_mzpaf_arena(
    mut geom: QueryCollection<IonAnnot>,
    frag_intens: Vec<f32>,
) -> Result<LibraryArena, LibraryReadingError> {
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

/// A single spectral-library format reader. Adding a format = one struct + one
/// line in [`registry`], instead of editing an enum, a method, and an ordered
/// try-chain.
pub trait LibraryReader: Send + Sync {
    fn name(&self) -> &'static str;
    /// Cheap probe: header bytes / extension / first data row. Must not read the
    /// whole file.
    fn sniff(&self, path: &Path) -> bool;
    /// Read the whole library into the arena.
    ///
    /// Formats that go through [`ElutionGroupCollection`] adapt it with
    /// [`LibraryArena::from_elution_groups`]; the binary `.speclib` and
    /// mzSpecLib readers build the arena directly, because only that path
    /// carries the reference-intensity sidecar.
    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError>;
}

struct MzSpecLibReader;
struct DiannSpeclibReader;
struct DiannParquetReader;
struct DiannTsvReader;
struct SpectronautReader;
struct SkylineReader;
struct JsonReader;

impl LibraryReader for MzSpecLibReader {
    fn name(&self) -> &'static str {
        "mzspeclib"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_mzspeclib_library_file(path)
    }

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        read_mzspeclib_library_file(path)
    }
}

impl LibraryReader for DiannSpeclibReader {
    fn name(&self) -> &'static str {
        "diann-speclib"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_speclib_library_file(path)
    }

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        read_diann_speclib_library_file(path)
    }
}

impl LibraryReader for DiannParquetReader {
    fn name(&self) -> &'static str {
        "diann-parquet"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_parquet_library_file(path)
    }

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        let egs = read_diann_parquet(path).map_err(|e| {
            warn!("Failed to read DIA-NN parquet library file: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        LibraryArena::from_elution_groups(ElutionGroupCollection::MzpafLabels(
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

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        let egs = read_diann_tsv(path).map_err(|e| {
            warn!("Failed to read DIA-NN TSV library file: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        LibraryArena::from_elution_groups(ElutionGroupCollection::MzpafLabels(
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

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        let egs = read_spectronaut_tsv(path).map_err(|e| {
            warn!("Failed to read Spectronaut TSV library file: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        LibraryArena::from_elution_groups(ElutionGroupCollection::MzpafLabels(
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

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        let egs = read_skyline_csv(path).map_err(|e| {
            warn!("Failed to read Skyline transition list: {:?}", e);
            LibraryReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        LibraryArena::from_elution_groups(ElutionGroupCollection::MzpafLabels(
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

    fn read(&self, path: &Path) -> Result<LibraryArena, LibraryReadingError> {
        LibraryArena::from_elution_groups(ElutionGroupCollection::try_read_json(path)?)
    }
}

/// Readers in dispatch order: most specific first, ending with the
/// always-sniffs-true JSON fallback.
///
/// mzSpecLib and `.speclib` lead because their probes are exact (a magic first
/// line and a version-gated header), so they cannot steal another format's
/// file — and `.speclib`'s read is the only one that surfaces an
/// `UnsupportedSpeclibVersion` diagnostic, which a later reader would mask.
fn registry() -> &'static [&'static dyn LibraryReader] {
    &[
        &MzSpecLibReader,
        &DiannSpeclibReader,
        &DiannParquetReader,
        &DiannTsvReader,
        &SpectronautReader,
        &SkylineReader,
        &JsonReader,
    ]
}

pub fn read_library_file<T: AsRef<Path>>(path: T) -> Result<LibraryArena, LibraryReadingError> {
    let path = path.as_ref();
    let mut last_err = None;
    for reader in registry() {
        if reader.sniff(path) {
            info!("Dispatching library read to {}", reader.name());
            match reader.read(path) {
                Ok(arena) => return Ok(arena),
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

#[cfg(test)]
mod tests {
    use super::*;

    fn ion(s: &str) -> IonAnnot {
        IonAnnot::try_from(s).expect("valid annotation")
    }

    /// Fragment labels must be unique within a precursor: `linear_get` is
    /// first-match, and `ExpectedIntensities::try_from_pairs` rejects
    /// duplicates outright — which timsseek's scoring pipeline `.expect()`s.
    /// A duplicate reaching the arena is a panic mid-search, so this is the
    /// invariant that stops it.
    #[test]
    fn duplicate_labels_collapse_onto_the_more_intense_peak() {
        let mut set = FragmentSet::with_capacity(3);
        assert_eq!(set.insert(ion("y1"), 175.1, 0.5), Inserted::Added);
        assert_eq!(set.insert(ion("b3^2"), 200.0, 0.9), Inserted::Added);

        // Weaker duplicate: kept peak is unchanged.
        assert_eq!(set.insert(ion("y1"), 999.9, 0.1), Inserted::Collapsed);
        assert_eq!(set.frags()[0], (ion("y1"), 175.1));

        // Stronger duplicate: takes over both m/z and intensity.
        assert_eq!(set.insert(ion("y1"), 175.2, 0.8), Inserted::Collapsed);
        assert_eq!(set.frags()[0], (ion("y1"), 175.2));

        assert_eq!(set.len(), 2, "a collision must not grow the set");
        let mut sidecar = vec![0.0];
        set.extend_sidecar(&mut sidecar);
        assert_eq!(
            sidecar,
            vec![0.0, 0.8, 0.9],
            "the sidecar stays parallel to frags(), appended in order"
        );

        // Charge is part of the label, so these do not collide.
        assert_eq!(set.insert(ion("y1^2"), 88.0, 0.3), Inserted::Added);
        assert_eq!(set.len(), 3);
    }

    #[test]
    fn clear_keeps_the_set_reusable() {
        let mut set = FragmentSet::with_capacity(2);
        set.insert(ion("y1"), 175.1, 1.0);
        set.clear();
        assert!(set.is_empty());
        assert_eq!(set.insert(ion("y1"), 175.1, 1.0), Inserted::Added);
    }
}
