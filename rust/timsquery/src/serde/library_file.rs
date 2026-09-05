use super::diann_io::{
    read_parquet_library_file as read_diann_parquet,
    read_targets as read_diann_tsv,
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
use super::precursor_extras::PrecursorExtras;
use super::skyline_io::{
    read_targets as read_skyline_csv,
    sniff_skyline_library_file,
};
use super::spectronaut_io::{
    read_targets as read_spectronaut_tsv,
    sniff_spectronaut_library_file,
};
use crate::Target;
use crate::ion::IonAnnot;
use crate::models::capabilities::{
    DecoyPolicy,
    LoadPolicy,
};
use crate::models::{
    Row,
    SourceIdError,
    TargetCapabilities,
    TargetColumns,
    TargetColumnsBuilder,
};
use std::path::Path;
use std::sync::Arc;
use tracing::{
    debug,
    info,
    warn,
};

#[derive(Debug)]
pub enum TargetReadingError {
    IoError(std::io::Error),
    SerdeJsonError(serde_json::Error),
    ElutionGroupInputError(ElutionGroupInputError),
    UnableToParseElutionGroups,
    /// Caller-supplied ids that cannot serve as a result key.
    SourceId(SourceIdError),
    /// A `.speclib` whose version is newer (more negative) than this reader
    /// supports.
    UnsupportedSpeclibVersion(i32),
    /// A `.speclib` structural desync / bad byte count / unexpected value, with
    /// context.
    SpeclibParse(String),
}

impl From<serde_json::Error> for TargetReadingError {
    fn from(err: serde_json::Error) -> Self {
        TargetReadingError::SerdeJsonError(err)
    }
}
impl From<SourceIdError> for TargetReadingError {
    fn from(err: SourceIdError) -> Self {
        TargetReadingError::SourceId(err)
    }
}
impl From<ElutionGroupInputError> for TargetReadingError {
    fn from(err: ElutionGroupInputError) -> Self {
        TargetReadingError::ElutionGroupInputError(err)
    }
}

#[derive(Debug)]
pub enum ElutionGroupCollection {
    StringLabels(Vec<Target<String>>, Option<Vec<PrecursorExtras>>),
    MzpafLabels(Vec<Target<IonAnnot>>, Option<Vec<PrecursorExtras>>),
}

impl ElutionGroupCollection {
    pub fn len(&self) -> usize {
        match self {
            ElutionGroupCollection::StringLabels(egs, _) => egs.len(),
            ElutionGroupCollection::MzpafLabels(egs, _) => egs.len(),
        }
    }

    fn try_read_json(path: &Path) -> Result<Self, TargetReadingError> {
        let file_content = std::fs::read_to_string(path).map_err(TargetReadingError::IoError)?;
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
            Err(_) => Err(TargetReadingError::UnableToParseElutionGroups),
        }
    }

    fn try_deser_inputed(content: &str) -> Result<Self, TargetReadingError> {
        debug!("Attempting deserialization of elution group inputs");
        // Fragments with no labels at all. Minted as unknown-ion annotations
        // (`?1`, `?2`, ...) rather than as integers: the arena has no
        // integer-labelled variant, so integer labels were built and then
        // thrown away, which made this input unloadable.
        //
        // Tried first because `fragment_labels: null` deserializes against
        // every label type, so a later arm would claim it and mint the wrong
        // kind.
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<IonAnnot>>>(content) {
            if eg_inputs.first().is_some_and(|x| x.needs_fragment_labels()) {
                debug!("Filling missing fragment labels with unknown-ion annotations");
                let filled: Result<Vec<_>, _> = eg_inputs
                    .into_iter()
                    .map(|x| x.try_fill_labels_annot())
                    .collect();
                let out: Result<Vec<Target<IonAnnot>>, ElutionGroupInputError> =
                    filled?.into_iter().map(|x| x.try_into()).collect();
                return Ok(ElutionGroupCollection::MzpafLabels(out?, None));
            }
            let out: Result<Vec<Target<IonAnnot>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::MzpafLabels(out?, None));
        }
        debug!("Attempting to deserialize elution group inputs with string labels");
        if let Ok(eg_inputs) = serde_json::from_str::<Vec<ElutionGroupInput<String>>>(content) {
            let out: Result<Vec<Target<String>>, ElutionGroupInputError> =
                eg_inputs.into_iter().map(|x| x.try_into()).collect();
            return Ok(ElutionGroupCollection::StringLabels(out?, None));
        }
        Err(TargetReadingError::UnableToParseElutionGroups)
    }

    fn try_deser_direct(content: &str) -> Result<Self, TargetReadingError> {
        // Ion annotations first, then opaque strings. Integer label types are
        // not tried: the arena has no variant for them, so claiming the input
        // here only turned a loadable file into a rejection.
        debug!("Attempting direct deserialization of elution groups");
        debug!("Attempting to deserialize elution groups with mzpaf labels");
        if let Ok(egs) = serde_json::from_str::<Vec<Target<IonAnnot>>>(content) {
            return Ok(ElutionGroupCollection::MzpafLabels(egs, None));
        }
        debug!("Attempting to deserialize elution groups with string labels");
        if let Ok(egs) = serde_json::from_str::<Vec<Target<String>>>(content) {
            return Ok(ElutionGroupCollection::StringLabels(egs, None));
        }
        Err(TargetReadingError::UnableToParseElutionGroups)
    }
}

/// The label-typed columnar library store returned by [`read_targets`].
///
/// One funnel: every format lands in exactly one variant. DIA-NN family formats
/// (`.speclib`/TSV/parquet) carry ion-chemistry (`IonAnnot`) labels and land in
/// [`TargetTable::Mzpaf`]; string-labelled JSON lands in [`TargetTable::Str`].
///
/// `frag_intens` is the reference-intensity sidecar, parallel to
/// `geom.frag_labels`/`geom.frag_mzs` (same length). The columnar store itself
/// stores query geometry; readers that supply reference intensities populate
/// the sidecar for scoring. Geometry-only `Mzpaf` inputs leave it `None`;
/// the `Str` variant has no intensity sidecar. Extraction ignores intensities.
pub enum TargetTable {
    Mzpaf {
        geom: TargetColumns<IonAnnot>,
        frag_intens: Option<Vec<f32>>,
    },
    Str {
        geom: TargetColumns<Arc<str>>,
    },
}

impl TargetTable {
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
        decoys: DecoyPolicy,
        egs: Vec<Target<IonAnnot>>,
        rows: Vec<PrecursorExtras>,
    ) -> Result<Self, TargetReadingError> {
        if egs.len() != rows.len() {
            return Err(TargetReadingError::SpeclibParse(format!(
                "elution groups ({}) and reader extras ({}) length mismatch",
                egs.len(),
                rows.len()
            )));
        }

        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        let mut frag_intens: Vec<f32> = Vec::new();

        for (eg, row) in egs.iter().zip(rows) {
            // The one filter point for every format that reaches the arena
            // through this adapter. Dropping a shipped decoy here rather than
            // after sealing means its fragments are never parsed and the
            // arena's own `MassShift` downgrade never sees it.
            if !decoys.accepts(row.is_decoy) {
                continue;
            }
            // Reference intensities keyed by fragment label (see fn docs).
            let lookup: std::collections::HashMap<IonAnnot, f32> =
                row.relative_intensities.into_iter().collect();
            let frags: Vec<(IonAnnot, f64)> = eg.iter_fragments().map(|(l, mz)| (*l, mz)).collect();
            for (label, _) in &frags {
                let intensity = lookup.get(label).ok_or_else(|| {
                    TargetReadingError::SpeclibParse(format!(
                        "fragment {label:?} of precursor {:?} has no reference intensity",
                        row.modified_peptide
                    ))
                })?;
                frag_intens.push(*intensity);
            }
            geom.push_row(Row {
                precursor_mz: eg.precursor_mz(),
                charge: eg.precursor_charge(),
                rt_seconds: eg.rt_seconds(),
                mobility: eg.mobility_ook0(),
                frags: &frags,
                seq_strip: &row.stripped_peptide,
                seq_mod: &row.modified_peptide,
                is_decoy: row.is_decoy,
                id: Some(eg.id().to_owned_id()),
                ..Default::default()
            });
        }

        if frag_intens.len() != geom.n_fragments() {
            return Err(TargetReadingError::SpeclibParse(format!(
                "reference-intensity sidecar ({}) must stay parallel to the fragment-label arena ({})",
                frag_intens.len(),
                geom.n_fragments(),
            )));
        }

        let geom = geom.seal(decoys)?;
        Ok(TargetTable::Mzpaf {
            geom,
            frag_intens: Some(frag_intens),
        })
    }

    /// Adapt the legacy [`ElutionGroupCollection`] (produced by the non-speclib
    /// readers) into the arena. `IonAnnot`-labelled groups become
    /// [`TargetTable::Mzpaf`]; string-labelled groups become
    /// [`TargetTable::Str`]. Integer-labelled groups have no arena variant (no
    /// live consumer) and are rejected.
    ///
    /// When the reader supplied per-precursor extras (DIA-NN/Skyline/
    /// Spectronaut all carry `relative_intensities`, `modified`/`stripped`
    /// sequences and an `is_decoy` flag), the reference intensities are threaded
    /// into the `frag_intens` sidecar so the timsseek bridge can score against
    /// them; sequences and the decoy flag are threaded into the arena too. When
    /// no extras were supplied (e.g. plain `IonAnnot` JSON, which only exercises
    /// the extraction/geometry path), `frag_intens` stays `None` -- matching the
    /// historical behavior where timsseek rejected that shape.
    fn from_elution_groups(
        egc: ElutionGroupCollection,
        decoys: DecoyPolicy,
    ) -> Result<Self, TargetReadingError> {
        match egc {
            ElutionGroupCollection::MzpafLabels(egs, Some(extras)) => {
                Self::mzpaf_with_intensities(decoys, egs, extras)
            }
            ElutionGroupCollection::MzpafLabels(egs, None) => {
                let mut geom =
                    TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
                for eg in &egs {
                    let frags: Vec<(IonAnnot, f64)> =
                        eg.iter_fragments().map(|(l, mz)| (*l, mz)).collect();
                    geom.push_row(Row {
                        precursor_mz: eg.precursor_mz(),
                        charge: eg.precursor_charge(),
                        rt_seconds: eg.rt_seconds(),
                        mobility: eg.mobility_ook0(),
                        frags: &frags,
                        id: Some(eg.id().to_owned_id()),
                        ..Default::default()
                    });
                }
                let geom = geom.seal(decoys)?;
                Ok(TargetTable::Mzpaf {
                    geom,
                    frag_intens: None,
                })
            }
            ElutionGroupCollection::StringLabels(egs, _) => {
                // String-labelled arenas carry no ion chemistry and ship no
                // decoys: sequence/fragment features unavailable, decoys off.
                // Sealed `Never` regardless of what the caller asked for -- a
                // derived decoy shifts the precursor but leaves an opaque
                // fragment label untouched (`DecoyShift for Arc<str>` is the
                // identity), so the variants would be scored against the
                // target's own fragments.
                let mut geom = TargetColumnsBuilder::with_capabilities(
                    TargetCapabilities::default_unlabeled(),
                );
                for eg in &egs {
                    let frags: Vec<(Arc<str>, f64)> = eg
                        .iter_fragments()
                        .map(|(l, mz)| (Arc::<str>::from(l.as_str()), mz))
                        .collect();
                    geom.push_row(Row {
                        precursor_mz: eg.precursor_mz(),
                        charge: eg.precursor_charge(),
                        rt_seconds: eg.rt_seconds(),
                        mobility: eg.mobility_ook0(),
                        frags: &frags,
                        id: Some(eg.id().to_owned_id()),
                        ..Default::default()
                    });
                }
                let geom = geom.seal(DecoyPolicy::Never)?;
                Ok(TargetTable::Str { geom })
            }
        }
    }
}

/// A single spectral-library format reader. Adding a format is one struct and
/// one line in [`registry`], rather than an enum, a method and an ordered
/// try-chain.
///
/// `read` returns the arena, not an intermediate. Formats that declare decoy
/// groups have no way to express them through
/// [`ElutionGroupCollection`], so returning that instead forced two of them to
/// bypass the registry through hardcoded `if sniff(..)` blocks; the ones that
/// do go through it call [`TargetTable::from_elution_groups`] themselves.
pub trait LibraryReader: Send + Sync {
    fn name(&self) -> &'static str;
    /// Cheap probe: header bytes / extension / first data row. Must not read the
    /// whole file.
    fn sniff(&self, path: &Path) -> bool;
    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError>;
}

struct DiannParquetReader;
struct DiannTsvReader;
struct SpectronautReader;
struct SkylineReader;

impl LibraryReader for DiannParquetReader {
    fn name(&self) -> &'static str {
        "diann-parquet"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_parquet_library_file(path)
    }

    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError> {
        let egs = read_diann_parquet(path).map_err(|e| {
            warn!("Failed to read DIA-NN parquet library file: {:?}", e);
            TargetReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        TargetTable::from_elution_groups(
            ElutionGroupCollection::MzpafLabels(egs, Some(extras)),
            policy.decoys,
        )
    }
}

impl LibraryReader for DiannTsvReader {
    fn name(&self) -> &'static str {
        "diann-tsv"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_library_file(path)
    }

    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError> {
        let egs = read_diann_tsv(path).map_err(|e| {
            warn!("Failed to read DIA-NN TSV library file: {:?}", e);
            TargetReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        TargetTable::from_elution_groups(
            ElutionGroupCollection::MzpafLabels(egs, Some(extras)),
            policy.decoys,
        )
    }
}

impl LibraryReader for SpectronautReader {
    fn name(&self) -> &'static str {
        "spectronaut-tsv"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_spectronaut_library_file(path).is_ok()
    }

    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError> {
        let egs = read_spectronaut_tsv(path).map_err(|e| {
            warn!("Failed to read Spectronaut TSV library file: {:?}", e);
            TargetReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        TargetTable::from_elution_groups(
            ElutionGroupCollection::MzpafLabels(egs, Some(extras)),
            policy.decoys,
        )
    }
}

impl LibraryReader for SkylineReader {
    fn name(&self) -> &'static str {
        "skyline-csv"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_skyline_library_file(path).is_ok()
    }

    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError> {
        let egs = read_skyline_csv(path).map_err(|e| {
            warn!("Failed to read Skyline transition list: {:?}", e);
            TargetReadingError::UnableToParseElutionGroups
        })?;
        let (egs, extras): (Vec<_>, Vec<_>) = egs.into_iter().unzip();
        TargetTable::from_elution_groups(
            ElutionGroupCollection::MzpafLabels(egs, Some(extras)),
            policy.decoys,
        )
    }
}

struct DiannSpeclibReader;
struct MzSpecLibReader;

impl LibraryReader for DiannSpeclibReader {
    fn name(&self) -> &'static str {
        "diann-speclib"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_diann_speclib_library_file(path)
    }

    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError> {
        read_diann_speclib_library_file(path, policy.decoys)
    }
}

impl LibraryReader for MzSpecLibReader {
    fn name(&self) -> &'static str {
        "mzspeclib"
    }

    fn sniff(&self, path: &Path) -> bool {
        sniff_mzspeclib_library_file(path)
    }

    fn read(&self, path: &Path, policy: LoadPolicy) -> Result<TargetTable, TargetReadingError> {
        read_mzspeclib_library_file(path, policy)
    }
}

// JSON is deliberately not in the registry: it has no sniff, being what a file
// that matched nothing else is tried as. It is the terminal arm of
// `read_targets_with`.

/// Every format, in sniff order.
///
/// Order matters where sniffs overlap, and the array is the only thing that
/// states it. `diann-speclib` is first because its `read` is the one that can
/// surface an `UnsupportedSpeclibVersion`, which its sniff cannot; `mzspeclib`
/// precedes nothing in particular but must follow the binary formats, since its
/// fallback probe reads a line of text.
fn registry() -> &'static [&'static dyn LibraryReader] {
    &[
        &DiannSpeclibReader,
        &DiannParquetReader,
        &MzSpecLibReader,
        &DiannTsvReader,
        &SpectronautReader,
        &SkylineReader,
    ]
}

/// Read a library, keeping whatever decoy rows it ships.
///
/// The common case, and the one every caller wanted before decoy handling
/// became a parameter. Use [`read_targets_with`] to drop them at read time.
pub fn read_targets<T: AsRef<Path>>(path: T) -> Result<TargetTable, TargetReadingError> {
    read_targets_with(
        path,
        LoadPolicy {
            decoys: DecoyPolicy::Never,
            ..Default::default()
        },
    )
}

/// Read a library, deciding up front what a load does with what the file ships.
///
/// [`DecoyPolicy::Force`] drops shipped decoys before they reach the arena, so a
/// caller regenerating its own never pays to parse the file's and the arena
/// never sees a shipped decoy to downgrade `MassShift` over.
pub fn read_targets_with<T: AsRef<Path>>(
    path: T,
    policy: LoadPolicy,
) -> Result<TargetTable, TargetReadingError> {
    let path = path.as_ref();
    let mut last_err = None;
    for reader in registry() {
        if reader.sniff(path) {
            info!("Dispatching library read to {}", reader.name());
            match reader.read(path, policy) {
                Ok(arena) => return Ok(arena),
                // A sniff can fire on a file the reader then fails to parse
                // (overlapping sniffs). Fall through to the next candidate
                // instead of committing to the first sniff. Keep the FIRST
                // error: readers run specific -> generic, so the earliest
                // sniff-and-fail carries the most useful diagnostic, such as a
                // `.speclib` desync, which JSON's generic error would clobber.
                Err(e) => {
                    warn!("{} sniffed but failed to read: {:?}", reader.name(), e);
                    last_err.get_or_insert(e);
                }
            }
        }
    }
    // Nothing sniffed, or everything that sniffed failed to parse. Either way
    // JSON is what is left to try; its own error is reported only if no reader
    // got further, since a `.speclib` desync says more than "not valid JSON".
    info!("Dispatching library read to json (terminal)");
    match ElutionGroupCollection::try_read_json(path) {
        Ok(egs) => TargetTable::from_elution_groups(egs, policy.decoys),
        Err(e) => Err(last_err.unwrap_or(e)),
    }
}

#[cfg(test)]
mod tests {
    use std::path::{
        Path,
        PathBuf,
    };

    /// What a result would actually be keyed by: the ids the sealed arena
    /// holds, read back through the public funnel. Asserting on a `Target`
    /// instead would pass even if sealing dropped the ids.
    fn arena_ids(path: &Path) -> Vec<String> {
        // Both arms read the same way; the label type is what differs.
        fn ids<L: crate::KeyLike>(geom: &crate::models::TargetColumns<L>) -> Vec<String> {
            geom.rows().map(|r| geom.output_id(r).to_string()).collect()
        }
        match super::read_targets(path).expect("fixture loads") {
            super::TargetTable::Mzpaf { geom, .. } => ids(&geom),
            super::TargetTable::Str { geom } => ids(&geom),
        }
    }

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/diann_io_files")
            .join(name)
    }

    /// JSON naming no fragment labels loads, with labels minted.
    ///
    /// The arena has no integer-labelled variant, so minting `u8` labels for
    /// this input built something that was then rejected -- and
    /// `try_fill_labels_annot`, which mints the loadable kind, had no callers at
    /// all. Unknown-ion annotations are what a fragment with no stated identity
    /// gets everywhere else in the tree.
    #[test]
    fn json_without_fragment_labels_loads_with_minted_ones() {
        let arena = super::read_targets(fixture("labelless_targets.json"))
            .expect("label-less JSON is a supported input");
        let super::TargetTable::Mzpaf { geom, .. } = arena else {
            panic!("minted labels carry ion chemistry, so this is the Mzpaf arena");
        };
        assert_eq!(geom.n_rows(), 1);
        let row = geom.rows().next().unwrap();
        let labels: Vec<String> = geom
            .frag_labels(row)
            .iter()
            .map(ToString::to_string)
            .collect();
        assert_eq!(labels, ["?1", "?2"], "one unknown ion per fragment");
    }

    /// Every format goes through `registry()`. Two of them used to be
    /// hardcoded `if sniff(..)` blocks ahead of the loop, because the trait
    /// returned an intermediate they could not express decoy groups through.
    /// A third bypass is the failure mode this guards.
    #[test]
    fn every_format_is_in_the_registry() {
        let names: Vec<_> = super::registry().iter().map(|r| r.name()).collect();
        assert_eq!(
            names,
            [
                "diann-speclib",
                "diann-parquet",
                "mzspeclib",
                "diann-tsv",
                "spectronaut-tsv",
                "skyline-csv",
            ],
            "a format outside this list is dispatched somewhere else"
        );
    }

    /// Sniff order is what the array says. `mzspeclib` falls back to reading a
    /// line of text when the name is unfamiliar, so it has to follow the binary
    /// formats or it probes a parquet file as text.
    #[test]
    fn binary_formats_are_sniffed_before_the_text_probe() {
        let position = |name: &str| {
            super::registry()
                .iter()
                .position(|r| r.name() == name)
                .expect("registered")
        };
        assert!(position("diann-parquet") < position("mzspeclib"));
        assert!(position("diann-speclib") < position("mzspeclib"));
    }

    /// `sample_lib.tsv` carries `transition_group_id`; `sample_lib.txt` does
    /// not, so between them they cover propagation and the minted fallback.
    #[test]
    fn a_diann_tsv_name_survives_into_the_arena() {
        assert_eq!(
            arena_ids(&fixture("sample_lib.tsv")),
            ["AAAAAAALQAK2"],
            "DIA-NN's own name for the precursor, not a counter"
        );
        assert_eq!(arena_ids(&fixture("sample_lib.txt")), ["0", "1"]);
    }

    /// Same, for the parquet variant, where DIA-NN 2.2 spells it
    /// `Precursor.Id`. The Carafe-written fixture has no such column.
    #[test]
    fn a_diann_parquet_name_survives_into_the_arena() {
        assert_eq!(
            arena_ids(&fixture("sample_pq_speclib.parquet"))[0],
            "GREEWESAALQNANTK3"
        );
        assert_eq!(arena_ids(&fixture("carafe_pq_speclib.parquet"))[0], "0");
    }
}
