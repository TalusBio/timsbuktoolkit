use crate::IonAnnot;
use crate::data_sources::reference_library::ReferenceLibrary;
use crate::errors::LibraryReadingError;
use crate::fragment_mass::{
    IsotopeSource,
    isotope_dist_or_averagine,
};
use crate::models::sequence::{
    normalize_to_proforma,
    parse_sequence,
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
use timsquery::models::QueryCollection;
use timsquery::models::capabilities::SeqFeatureState;
use timsquery::serde::read_library_file as read_timsquery_library;
use timsquery::utils::constants::PROTON_MASS;

/// The serializable, on-disk form of a native speclib element. Kept backwards
/// compatible; the load path builds the columnar `ReferenceLibrary` arena
/// directly from these elements (see `Speclib::from_file_with_format`).
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

/// Strip mod annotations — anything inside `(...)` or `[...]` — from a
/// sequence, leaving the bare residue string. The native format ships one
/// (modified) sequence per precursor; the arena's composition-isotope path
/// needs the stripped residues.
fn strip_mods(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth: i32 = 0;
    for c in s.chars() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth = (depth - 1).max(0),
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

/// Summary of a [`finalize_reference_library`] call, for load-time logging.
#[derive(Debug, Clone, Copy)]
pub struct LoadReport {
    /// Physical stored rows (pre decoy expansion), i.e. `QueryCollection::n_rows`.
    pub n_rows: usize,
    pub n_averagine_fallback: usize,
    pub sequence_features: SeqFeatureState,
}

/// Finalize a freshly-narrowed lazy `ReferenceLibrary` arena: apply the decoy
/// strategy, seal, run the whole-library parse gate + averagine tally, and set
/// `caps.sequence_features`. This is the single shared tail of the DEFAULT
/// `.speclib` load (see `speclib_data_flow.md`) — the memory-optimized path
/// that avoids the 9 GB peak RSS of the fully-materialized target+2-decoy
/// expansion.
///
/// `policy` is the raw CLI decoy policy: this is the single place it is resolved
/// (via `map_decoy_strategy`, keyed on whether the arena already ships decoys)
/// and stamped onto `caps.decoys` BEFORE `seal()`, so the seal's
/// `LazyMassShift -> Passthrough` downgrade (the Task-4 gate) sees it. The parse
/// gate walks the MODIFIED sequence blob (the form
/// `RefQuery::materialize_peptide_in_group` parses) and, if any row fails,
/// disables sequence-derived features library-wide. The same pass counts
/// averagine isotope fallbacks for the returned `LoadReport`.
fn finalize_reference_library(
    mut geom: QueryCollection<IonAnnot>,
    frag_intens: Vec<f32>,
    policy: crate::models::DecoyPolicy,
) -> (ReferenceLibrary, LoadReport) {
    let n_stored_decoys = geom.is_decoy.iter().filter(|&&d| d).count();
    geom.caps.decoys = crate::models::map_decoy_strategy(policy, n_stored_decoys > 0);
    geom.seal();

    let n_rows = geom.n_rows();

    // Report the effective decoy strategy (post-seal: `seal()` downgrades
    // LazyMassShift -> Passthrough if the library ships its own decoys). This
    // restores the load-time notice that the search space is expanded by
    // synthetic mass-shift decoys.
    {
        use timsquery::models::capabilities::DecoyStrategy;
        match geom.caps.decoys {
            DecoyStrategy::LazyMassShift { n_decoys, .. } => {
                tracing::warn!(
                    "Library contains no decoys. Generating synthetic ±CH2 mass-shift \
                     decoys: {}x search space ({} targets -> {} scored entries)",
                    n_decoys as usize + 1,
                    n_rows,
                    geom.expanded_len(),
                );
            }
            DecoyStrategy::Passthrough => {
                tracing::info!(
                    "Library ships {} decoys; using them as-is (Passthrough, no synthetic \
                     decoys generated)",
                    n_stored_decoys,
                );
            }
            DecoyStrategy::None => {
                tracing::info!(
                    "Decoy strategy None; scoring {} stored rows as-is (targets + any \
                     shipped decoys)",
                    n_rows
                );
            }
        }
    }
    let mut all_parsable = true;
    let mut n_averagine_fallback = 0usize;
    for tgt in 0..n_rows {
        if all_parsable {
            let modified = &geom.seq_mod_blob[geom.seq_mod_range(tgt)];
            let normalized = normalize_to_proforma(modified);
            if parse_sequence(&normalized).is_none() {
                all_parsable = false;
            }
        }
        let stripped = &geom.seq_strip_blob[geom.seq_strip_range(tgt)];
        let charge = geom.charge[tgt] as f64;
        let neutral_mass = geom.precursor_mz[tgt] * charge - charge * PROTON_MASS;
        let (isotope_src, _envelope) = isotope_dist_or_averagine(stripped, neutral_mass);
        if isotope_src == IsotopeSource::Averagine {
            n_averagine_fallback += 1;
        }
    }

    let sequence_features = if all_parsable {
        SeqFeatureState::Available
    } else {
        SeqFeatureState::Unavailable
    };
    geom.caps.sequence_features = sequence_features;

    if n_averagine_fallback > 0 {
        tracing::warn!(
            "{}/{} library entries used averagine isotope fallback",
            n_averagine_fallback,
            n_rows
        );
    }
    tracing::info!(
        "Lazy ReferenceLibrary arena ready: {} stored rows, sequence_features={:?}",
        n_rows,
        sequence_features
    );

    let report = LoadReport {
        n_rows,
        n_averagine_fallback,
        sequence_features,
    };

    (ReferenceLibrary { geom, frag_intens }, report)
}

/// The spectral library store. Collapsed to the single columnar
/// `ReferenceLibrary` arena representation (the materialized AOS path was
/// deleted in Task 9): both load paths produce a lazy arena, and scoring
/// iterates `RefQuery` flyweights via [`ReferenceLibrary::item_at`].
pub type Speclib = ReferenceLibrary;

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

/// Streams raw `SerSpeclibElement`s out of a native timsseek library file.
///
/// The native path builds the columnar arena directly from these elements (see
/// `Speclib::from_file_with_format`), so the reader stays at the serializable
/// element and does not eagerly build per-row scoring items.
pub struct SpeclibReader<'a> {
    inner: Box<dyn Iterator<Item = Result<SerSpeclibElement, LibraryReadingError>> + Send + 'a>,
}

impl<'a> SpeclibReader<'a> {
    pub fn new<R: Read + Send + 'a>(
        reader: R,
        format: SpeclibFormat,
    ) -> Result<Self, LibraryReadingError> {
        let inner: Box<dyn Iterator<Item = Result<SerSpeclibElement, LibraryReadingError>> + Send> =
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
    type Item = Result<SerSpeclibElement, LibraryReadingError>;

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
    type Item = Result<SerSpeclibElement, LibraryReadingError>;

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

                Some(Ok(elem))
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
    type Item = Result<SerSpeclibElement, LibraryReadingError>;

    fn next(&mut self) -> Option<Self::Item> {
        use serde::Deserialize;

        match SerSpeclibElement::deserialize(&mut self.deserializer) {
            Ok(elem) => Some(Ok(elem)),
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
    /// Whether every sequence in the library parsed (gates sequence-derived
    /// scoring features). Reads the sealed arena's `sequence_features` state.
    pub fn parsable_sequences(&self) -> bool {
        self.geom.caps.sequence_features == SeqFeatureState::Available
    }

    /// Mean number of fragments per entry (0.0 for an empty library).
    ///
    /// Every variant of a target shares the same fragment set, so the
    /// per-target arena length is exact, not a sample.
    pub fn avg_fragments(&self) -> f64 {
        let n_rows = self.geom.n_rows();
        if n_rows == 0 {
            return 0.0;
        }
        self.geom.frag_labels.len() as f64 / n_rows as f64
    }

    pub fn from_file(
        path: &Path,
        decoy_strategy: crate::models::DecoyPolicy,
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
        // `.speclib`/TSV/parquet, Spectronaut, Skyline, JSON), which returns a
        // label-generic `LibraryArena`. One path from here: narrow the arena
        // to the ion-annotated `ReferenceLibrary`, apply the decoy strategy,
        // seal, gate sequence features, and hand back the lazy arena.
        tracing::info!(
            "Loading library via timsquery format detection: {}",
            path.display()
        );
        let arena = read_timsquery_library(path)?;
        let ReferenceLibrary { geom, frag_intens } = ReferenceLibrary::try_from(arena)?;

        // Decoy resolution + seal + parse gate + averagine tally + `LoadReport`
        // all live in the one shared finalize path (see
        // `finalize_reference_library`); the report is logged there, so we drop
        // it here.
        let (lib, _report) = finalize_reference_library(geom, frag_intens, decoy_strategy);
        lib.log_entry_stats();
        Ok(lib)
    }

    pub fn from_file_with_format(
        path: &Path,
        format: SpeclibFormat,
        decoy_strategy: crate::models::DecoyPolicy,
    ) -> Result<Self, LibraryReadingError> {
        let file =
            std::fs::File::open(path).map_err(|e| LibraryReadingError::FileReadingError {
                source: e,
                context: "Error opening speclib file",
                path: PathBuf::from(path),
            })?;

        let reader = SpeclibReader::new(file, format)?;

        // Build the columnar arena directly from the streamed elements (same
        // lazy shape as the `.speclib` path), instead of collecting per-row
        // scoring items into an intermediate Vec. Each element's fragment labels/mzs/
        // intensities are parallel vectors in the native format, so the
        // reference-intensity sidecar is filled in fragment-push order.
        let mut geom =
            QueryCollection::with_capabilities(timsquery::models::LibCapabilities::default_diann());
        let mut frag_intens: Vec<f32> = Vec::new();

        for elem in reader {
            let elem = elem?;
            let eg = &elem.elution_group;
            if eg.fragment_labels.len() != eg.fragment_mzs.len()
                || eg.fragment_labels.len() != eg.fragment_intensities.len()
            {
                return Err(LibraryReadingError::UnsupportedFormat {
                    message: format!(
                        "speclib element {:?}: fragment labels ({}), mzs ({}) and intensities ({}) must be parallel",
                        elem.precursor.sequence,
                        eg.fragment_labels.len(),
                        eg.fragment_mzs.len(),
                        eg.fragment_intensities.len(),
                    ),
                });
            }

            let frags: Vec<(IonAnnot, f64)> = eg
                .fragment_labels
                .iter()
                .cloned()
                .zip(eg.fragment_mzs.iter().cloned())
                .collect();
            frag_intens.extend_from_slice(&eg.fragment_intensities);

            // The native format ships a single (modified) sequence; strip mod
            // annotations for the composition-isotope path.
            let modified = &elem.precursor.sequence;
            let stripped = strip_mods(modified);
            geom.push_row(
                eg.precursor_mz,
                elem.precursor.charge,
                eg.rt_seconds,
                eg.mobility_ook0,
                &frags,
                &stripped,
                modified,
                &[],
                elem.precursor.decoy,
            );
        }

        // Decoy resolution + seal happen inside the one shared finalize path,
        // exactly as the `.speclib` bridge above.
        let (lib, _report) = finalize_reference_library(geom, frag_intens, decoy_strategy);
        lib.log_entry_stats();
        Ok(lib)
    }

    /// Log a one-line summary of the lazy arena's shape at load time.
    fn log_entry_stats(&self) {
        tracing::info!(
            "Speclib stats: lazy arena, {} targets ({} flat scoring entries, {} total fragment slots)",
            self.geom.n_rows(),
            self.len(),
            self.geom.frag_labels.len(),
        );
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
    use crate::data_sources::reference_library::{
        ExpectedIntensity,
        RefQuery,
    };

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

    /// `Speclib` is now a type alias for `ReferenceLibrary` (Task 9 collapsed
    /// the enum), so a loaded library is already the lazy arena. This identity
    /// helper is kept so the fixture assertions below read as
    /// "get the arena" without churning every call site.
    fn expect_lazy(speclib: &Speclib) -> &ReferenceLibrary {
        speclib
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
            .expect("Failed to load DIA-NN TSV library");

        assert!(
            speclib.parsable_sequences(),
            "DIA-NN sample fixture should all parse with modified_peptide"
        );
    }

    #[test]
    fn test_load_skyline_csv_library() {
        use timsquery::traits::QueryGeom;

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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
            .expect("Failed to load Skyline CSV library");

        // Skyline routes through the timsquery bridge (`from_elution_groups`),
        // which now threads the reference intensities through, so it narrows to
        // a lazy `ReferenceLibrary` arena like the DIA-NN formats. No shipped decoys +
        // default IfMissing -> `LazyMassShift`.
        // Fixture has 14 PRTC targets, no decoys -> 14 targets + 28 mass-shift decoys
        assert_eq!(
            speclib.len(),
            42,
            "Expected 42 entries (14 targets + 28 decoys)"
        );

        let lib = expect_lazy(&speclib);
        let n_rows = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();
        assert_eq!(n_rows, 14, "Should have 14 targets");
        assert_eq!(n_decoys, 28, "Should have 28 decoys");

        // Isotope envelope should have been attached (3 isotopes) for every target
        for q in lib.iter().filter(|q| q.geom().variant() == 0) {
            assert_eq!(
                q.expected_precursor_envelope().len(),
                3,
                "Each target should have 3 isotopes in precursor envelope"
            );
            assert!(
                q.fragment_count() > 0,
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
            .expect("Failed to load TXT library");

        // The sample file has 2 unique precursors with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys
        assert_eq!(
            speclib.len(),
            6,
            "Expected 6 entries (2 targets + 4 decoys)"
        );

        let lib = expect_lazy(&speclib);
        let n_rows = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();

        assert_eq!(n_rows, 2, "Should have 2 targets");
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
            .expect("Failed to load Parquet library");

        // The sample parquet file has 3 unique precursors with no decoys
        // Should generate 3x entries: 3 targets + 6 decoys
        assert_eq!(
            speclib.len(),
            9,
            "Expected 9 entries (3 targets + 6 decoys)"
        );

        let lib = expect_lazy(&speclib);
        let n_rows = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();

        assert_eq!(n_rows, 3, "Should have 3 targets");
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
            .expect("Failed to load DIA-NN TSV library");

        // Test file has 2 targets with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys (2 per target)
        assert_eq!(
            speclib.len(),
            6,
            "Should have 6 entries (2 targets + 4 decoys)"
        );

        let lib = expect_lazy(&speclib);
        let n_rows = lib.iter().filter(|q| q.geom().variant() == 0).count();
        let n_decoys = lib.iter().filter(|q| q.geom().variant() != 0).count();

        assert_eq!(n_rows, 2, "Should have 2 target entries");
        assert_eq!(n_decoys, 4, "Should have 4 decoy entries (2 per target)");

        // Each target index (== decoy_group in the old materialized scheme)
        // should have exactly 3 flat variants (1 target + 2 decoys) — this is
        // guaranteed structurally by `ReferenceLibrary::item_at`'s `t,+,-`
        // packing, so assert it directly rather than re-deriving groups.
        let n_target_indices = lib.geom.n_rows();
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
            .expect("Failed to load DIA-NN TSV library");

        let lib = expect_lazy(&speclib);

        // Unified CH2 offset (see `map_decoy_strategy`), replacing the old
        // 12.0 (materialized `IfMissing`) / 14.0 (materialized `Force`) split.
        const CH2_MASS: f64 = 14.0;

        for tgt in 0..lib.geom.n_rows() as u32 {
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

        let speclib = Speclib::from_file(&test_file, crate::models::DecoyPolicy::default())
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

    /// End-to-end `Speclib::from_file` over the real DIA-NN HeLa `.speclib`
    /// fixture (the actual workload path). Proves: the arena narrows to a lazy
    /// library with targets, variant-0 is a target, and the intensity sidecar
    /// threaded through (`iter_expected_fragments` yields at least one pair).
    #[test]
    fn from_file_loads_lazy_library_from_speclib_fixture() {
        use crate::data_sources::reference_library::ScoredIdentity;

        let path = std::path::Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../timsquery/tests/speclib_io_files/diann-hela-diapasef-lib.speclib"
        ));
        assert!(path.exists(), "fixture should exist at {:?}", path);

        let speclib = Speclib::from_file(path, crate::models::DecoyPolicy::default())
            .expect("from_file should load the .speclib fixture");

        let lib = expect_lazy(&speclib);
        assert!(lib.len() > 0, "library should have entries");

        let first = lib.item_at(0);
        assert!(first.is_target(), "flat index 0 must be a target variant");

        let frags: Vec<_> = first.iter_expected_fragments().collect();
        assert!(
            !frags.is_empty(),
            "target 0 must expose at least one (label, intensity) pair \
             (proves the intensity sidecar threaded through)"
        );
    }

    /// Native `SerSpeclibElement` reader (ndjson) builds the lazy arena
    /// directly. The fixture ships one target + one stored decoy, so the
    /// Task-4 seal gate downgrades `LazyMassShift -> Passthrough`: the arena is
    /// 1:1 with the stored rows (no synthetic mass-shift expansion). Proves the
    /// native path produces a lazy `ReferenceLibrary` with the right length, target/
    /// decoy flags, and per-fragment reference intensities.
    #[test]
    fn from_file_with_format_native_ndjson_builds_lazy_arena() {
        use crate::data_sources::reference_library::ScoredIdentity;

        let target = SerSpeclibElement::new(
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
        let decoy = SerSpeclibElement::new(
            PrecursorEntry::new("KEDITPEP".to_string(), 2, true, 0),
            ReferenceEG::new(
                1,
                500.0,
                vec![0, 1, 2],
                vec![300.0, 400.0],
                vec![
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                ],
                vec![1.0, 0.5, 0.2],
                vec![0.6, 0.4],
                0.75,
                120.0,
            ),
        );

        let mut ndjson = String::new();
        ndjson.push_str(&serde_json::to_string(&target).unwrap());
        ndjson.push('\n');
        ndjson.push_str(&serde_json::to_string(&decoy).unwrap());
        ndjson.push('\n');

        let path = std::env::temp_dir().join(format!(
            "timsseek_native_fixture_{}.ndjson",
            std::process::id()
        ));
        std::fs::write(&path, ndjson).unwrap();

        let speclib = Speclib::from_file(&path, crate::models::DecoyPolicy::default())
            .expect("native ndjson should load");
        std::fs::remove_file(&path).ok();

        let lib = expect_lazy(&speclib);
        // Ships a decoy -> Passthrough -> 1 variant/row -> flat len == n_rows.
        assert_eq!(lib.geom.variants_per_row(), 1, "downgraded to Passthrough");
        assert_eq!(lib.len(), 2, "one target + one stored decoy, 1:1");

        assert!(lib.item_at(0).is_target(), "row 0 is the target");
        assert!(!lib.item_at(1).is_target(), "row 1 is the stored decoy");

        let frags: Vec<_> = lib.item_at(0).iter_expected_fragments().collect();
        assert_eq!(frags.len(), 2, "target ships two reference fragments");
        for (_label, intensity) in frags {
            assert!(intensity > 0.0, "reference intensities are positive");
        }
    }

    /// Negative parse-gate coverage on the arena path. `finalize_reference_library`
    /// walks every target's MODIFIED sequence blob; if ANY row fails
    /// `parse_sequence(normalize_to_proforma(..))`, sequence-derived features are
    /// disabled library-wide (`SeqFeatureState::Unavailable`). Here one target
    /// parses (`PEPTIDEK`) and one is poisoned (`GARBAGE!!!`: the `!` bytes are
    /// rejected by both the fast byte-walk parser and the rustyms fallback), so
    /// the gate must report `!parsable_sequences()`. This is the inverse of
    /// `test_diann_tsv_parsable_gate`, and the only test of the OFF branch after
    /// the AOS `test_parse_gate_off_on_poisoned_row` was removed in Task 9.
    #[test]
    fn from_file_native_ndjson_poisoned_row_disables_sequence_features() {
        let good = SerSpeclibElement::new(
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
        // Unparseable modified sequence: `!` is rejected by parse_sequence_fast
        // (`_ => return None`) and by the rustyms pro_forma fallback.
        let poisoned = SerSpeclibElement::new(
            PrecursorEntry::new("GARBAGE!!!".to_string(), 2, false, 1),
            ReferenceEG::new(
                1,
                600.0,
                vec![0, 1, 2],
                vec![300.0, 400.0],
                vec![
                    IonAnnot::try_from("y1").unwrap(),
                    IonAnnot::try_from("y2").unwrap(),
                ],
                vec![1.0, 0.5, 0.2],
                vec![0.7, 0.4],
                0.75,
                120.0,
            ),
        );

        let mut ndjson = String::new();
        ndjson.push_str(&serde_json::to_string(&good).unwrap());
        ndjson.push('\n');
        ndjson.push_str(&serde_json::to_string(&poisoned).unwrap());
        ndjson.push('\n');

        let path = std::env::temp_dir().join(format!(
            "timsseek_poisoned_fixture_{}.ndjson",
            std::process::id()
        ));
        std::fs::write(&path, ndjson).unwrap();

        let speclib = Speclib::from_file(&path, crate::models::DecoyPolicy::default())
            .expect("native ndjson should load even with an unparseable sequence");
        std::fs::remove_file(&path).ok();

        // The poisoned row flips the whole-library gate OFF.
        assert!(
            !speclib.parsable_sequences(),
            "an unparseable modified sequence must disable sequence features library-wide"
        );
    }

    /// The Mzpaf arena WITHOUT the intensity sidecar (`frag_intens: None`) is
    /// the TSV/parquet/Skyline bridge shape — scoring is intensity-driven, so
    /// narrowing it to a `ReferenceLibrary` must be an `Err` (the branch that
    /// causes the disclosed DIA-NN/Skyline regression).
    #[test]
    fn reference_library_rejects_mzpaf_without_intensities() {
        use timsquery::models::QueryCollection;
        use timsquery::models::capabilities::LibCapabilities;
        use timsquery::serde::LibraryArena;

        let mut geom = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        geom.push_target(
            900.4,
            2,
            1.0,
            1.0,
            &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            "PEP",
            "PEP",
            &[],
        );
        geom.seal();
        let arena = LibraryArena::Mzpaf {
            geom,
            frag_intens: None,
        };
        assert!(ReferenceLibrary::try_from(arena).is_err());
    }
}
