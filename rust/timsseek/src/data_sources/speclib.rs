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
    Write,
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
/// directly from these elements (see `Speclib::from_native_file`).
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
}

/// Only the fields the loader reads are in the format.
///
/// `id`, `decoy_group`, `precursor_labels` and `precursor_intensities` used to
/// be written here and dropped at read time — the isotope envelope is
/// recomputed from composition (`IsotopeStrategy::FromComposition`), so storing
/// it cost the writer a full ProForma parse per entry and gave isotopes two
/// sources of truth. Serde ignores unknown fields, so libraries that still
/// carry them load unchanged.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrecursorEntry {
    sequence: String,
    charge: u8,
    decoy: bool,
}

impl PrecursorEntry {
    pub fn new(sequence: String, charge: u8, decoy: bool) -> Self {
        Self {
            sequence,
            charge,
            decoy,
        }
    }
}

/// See [`PrecursorEntry`] for why the precursor-isotope fields are absent.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceEG {
    precursor_mz: f64,
    #[serde(alias = "fragment_mz")]
    fragment_mzs: Vec<f64>,
    fragment_labels: Vec<IonAnnot>,
    fragment_intensities: Vec<f32>,
    #[serde(alias = "mobility")]
    mobility_ook0: f32,
    rt_seconds: f32,
}

impl ReferenceEG {
    pub fn new(
        precursor_mz: f64,
        fragment_mzs: Vec<f64>,
        fragment_labels: Vec<IonAnnot>,
        fragment_intensities: Vec<f32>,
        mobility_ook0: f32,
        rt_seconds: f32,
    ) -> Self {
        Self {
            precursor_mz,
            fragment_mzs,
            fragment_labels,
            fragment_intensities,
            mobility_ook0,
            rt_seconds,
        }
    }
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
/// `.speclib` load — the memory-optimized path that avoids the 9 GB peak RSS
/// of a fully-materialized target+2-decoy expansion.
///
/// `policy` is the raw CLI decoy policy: this is the single place it is resolved
/// (via `map_decoy_strategy`, keyed on whether the arena already ships decoys)
/// and stamped onto `caps.decoys` BEFORE `seal()`, so the seal's
/// `LazyMassShift -> Passthrough` downgrade sees it. The parse
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
    // The first row that fails to parse, kept so the warning can name it. One
    // bad row anywhere disables sequence-derived scoring for the WHOLE library,
    // so "which row" is the only actionable part of that news.
    let mut first_unparsable: Option<String> = None;
    let mut n_averagine_fallback = 0usize;
    for tgt in 0..n_rows {
        if first_unparsable.is_none() {
            let modified = &geom.seq_mod_blob[geom.seq_mod_range(tgt)];
            let normalized = normalize_to_proforma(modified);
            if parse_sequence(&normalized).is_none() {
                first_unparsable = Some(modified.to_string());
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

    let sequence_features = match &first_unparsable {
        None => SeqFeatureState::Available,
        Some(sequence) => {
            tracing::warn!(
                "Sequence-derived scoring features are DISABLED for this entire library: \
                 {sequence:?} could not be parsed. Any non-Unimod modification \
                 (PSI-MOD, RESID, XL-MOD, cross-links) has this effect."
            );
            SeqFeatureState::Unavailable
        }
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
/// since deleted): both load paths produce a lazy arena, and scoring
/// iterates `RefQuery` flyweights via [`ReferenceLibrary::item_at`].
pub type Speclib = ReferenceLibrary;

/// Whether `path` names a native timsseek library, by EXTENSION ONLY.
///
/// This answers *which reader family*, not *which encoding*. Content sniffing
/// cannot answer it: a DIA-NN `.speclib` and a native library are both opaque
/// byte streams, and the point of the extension rule is that a native extension
/// commits to the native reader and surfaces its error rather than falling
/// through the timsquery registry to report some other reader's complaint.
/// Whether the bytes are zstd-wrapped is a separate question, and
/// [`SpeclibReader`] answers that one from the magic number.
fn is_native_extension(path: &Path) -> bool {
    let path_str = path.to_string_lossy().to_lowercase();
    // `.zst` and `.zstd` are both in the wild.
    let stem = path_str
        .strip_suffix(".zst")
        .or_else(|| path_str.strip_suffix(".zstd"))
        .unwrap_or(path_str.as_str());
    stem.ends_with(".ndjson")
}

/// Streams raw `SerSpeclibElement`s out of a native timsseek library file.
///
/// The native path builds the columnar arena directly from these elements (see
/// `Speclib::from_native_file`), so the reader stays at the serializable
/// element and does not eagerly build per-row scoring items.
///
/// The payload is always NDJSON; zstd only adds a decoder underneath, so the
/// boxing is over the byte source rather than over the line parser.
pub struct SpeclibReader<'a> {
    reader: Box<dyn BufRead + Send + 'a>,
}

/// Leading bytes of a zstd frame.
const ZSTD_MAGIC: [u8; 4] = [0x28, 0xB5, 0x2F, 0xFD];

impl<'a> SpeclibReader<'a> {
    /// Compression is detected from the first four bytes, not the file name, so
    /// a mislabelled `.ndjson` that is really zstd (or the reverse) still reads.
    pub fn new<R: Read + Send + 'a>(reader: R) -> Result<Self, LibraryReadingError> {
        let mut buffered = BufReader::new(reader);
        let compressed = buffered
            .fill_buf()
            .map_err(|source| LibraryReadingError::FileReadingError {
                source,
                context: "Error reading the start of the speclib",
                path: None,
            })?
            .starts_with(&ZSTD_MAGIC);

        let reader: Box<dyn BufRead + Send + 'a> = if compressed {
            let decoder = zstd::Decoder::with_buffer(buffered)
                .map_err(|source| LibraryReadingError::Decompression { source })?;
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(buffered)
        };

        Ok(SpeclibReader { reader })
    }
}

impl Iterator for SpeclibReader<'_> {
    type Item = Result<SerSpeclibElement, LibraryReadingError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Looping rather than recursing: a file with a long run of blank lines
        // would otherwise recurse once per line and blow the stack.
        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    if line.trim().is_empty() {
                        continue;
                    }
                    return Some(serde_json::from_str(&line).map_err(|e| {
                        LibraryReadingError::SpeclibParsingError {
                            source: e,
                            context: "Error parsing NDJSON line",
                        }
                    }));
                }
                Err(e) => {
                    return Some(Err(LibraryReadingError::FileReadingError {
                        source: e,
                        context: "Error reading line",
                        path: None,
                    }));
                }
            }
        }
    }
}

/// Writes a native timsseek library: one JSON object per line, zstd-wrapped.
///
/// The exact inverse of [`SpeclibReader`], so
/// what `speclib_build_cli` emits is what `Speclib::from_file` reads back.
pub struct SpeclibWriter<W: Write> {
    encoder: zstd::Encoder<'static, W>,
}

impl<W: Write> SpeclibWriter<W> {
    pub fn new_ndjson_zstd(writer: W) -> Result<Self, std::io::Error> {
        Ok(Self {
            encoder: zstd::Encoder::new(writer, 3)?,
        })
    }

    pub fn append(&mut self, elem: &SerSpeclibElement) -> Result<(), LibraryReadingError> {
        let io_err = |e: std::io::Error| LibraryReadingError::FileReadingError {
            source: e,
            context: "Error writing NDJSON",
            path: None,
        };
        serde_json::to_writer(&mut self.encoder, elem).map_err(|e| {
            LibraryReadingError::SpeclibParsingError {
                source: e,
                context: "Error serializing NDJSON line",
            }
        })?;
        // The newline is the record separator; without it the whole library is
        // one unreadable line.
        self.encoder.write_all(b"\n").map_err(io_err)
    }

    /// Flushes the zstd frame. Skipping this truncates the library.
    pub fn finish(self) -> Result<W, std::io::Error> {
        self.encoder.finish()
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
        decoy_policy: crate::models::DecoyPolicy,
    ) -> Result<Self, LibraryReadingError> {
        // msgpack was removed in the mzcore migration. Without this arm the
        // file reaches the JSON reader and is reported as invalid UTF-8, which
        // names neither the real problem nor the fix.
        let path_str = path.to_string_lossy().to_lowercase();
        if path_str.contains(".msgpack") {
            return Err(LibraryReadingError::UnsupportedFormat {
                message: format!(
                    "{}: msgpack speclibs are no longer supported. Rebuild with \
                     speclib_build_cli, which now emits .ndjson.zst",
                    path.display()
                ),
            });
        }

        // See `is_native_extension`: a native extension commits to the native
        // reader and surfaces its error. A `.speclib` matches no native
        // extension and falls through to the bridge -> timsquery registry.
        if is_native_extension(path) {
            tracing::info!("Loading native speclib from {}", path.display());
            return Self::from_native_file(path, decoy_policy);
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
        let (lib, _report) = finalize_reference_library(geom, frag_intens, decoy_policy);
        lib.log_entry_stats();
        Ok(lib)
    }

    /// Load a native timsseek library (NDJSON, optionally zstd-wrapped —
    /// [`SpeclibReader`] sniffs which).
    fn from_native_file(
        path: &Path,
        decoy_policy: crate::models::DecoyPolicy,
    ) -> Result<Self, LibraryReadingError> {
        let file =
            std::fs::File::open(path).map_err(|e| LibraryReadingError::FileReadingError {
                source: e,
                context: "Error opening speclib file",
                path: Some(PathBuf::from(path)),
            })?;

        let reader = SpeclibReader::new(file)?;

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
            let stripped = timsquery::utils::sequence::strip_mods(modified);
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
        let (lib, _report) = finalize_reference_library(geom, frag_intens, decoy_policy);
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
#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_sources::reference_library::{
        ExpectedIntensity,
        RefQuery,
    };

    /// One native-format element. Fragment m/z are positional stand-ins; no
    /// test here asserts on them, only on labels, intensities and the row
    /// metadata.
    fn element(
        sequence: &str,
        decoy: bool,
        precursor_mz: f64,
        labels: &[&str],
        intensities: &[f32],
    ) -> SerSpeclibElement {
        assert_eq!(labels.len(), intensities.len());
        SerSpeclibElement::new(
            PrecursorEntry::new(sequence.to_string(), 2, decoy),
            ReferenceEG::new(
                precursor_mz,
                (0..labels.len())
                    .map(|i| 300.0 + 100.0 * i as f64)
                    .collect(),
                labels
                    .iter()
                    .map(|l| IonAnnot::try_from(*l).expect("valid annotation"))
                    .collect(),
                intensities.to_vec(),
                0.75,
                120.0,
            ),
        )
    }

    /// `speclib_build_cli` writes with [`SpeclibWriter`] and timsseek reads
    /// with [`SpeclibReader`]; nothing else checks that the two agree, and a
    /// mismatch only shows up as an unreadable library at the end of a long
    /// Koina run.
    #[test]
    fn writer_output_reads_back_through_the_reader() {
        let record = element("PEPTIDEK", false, 450.5, &["y1", "b3^2"], &[0.9, 0.4]);

        let mut writer = SpeclibWriter::new_ndjson_zstd(Vec::new()).expect("encoder");
        writer.append(&record).expect("append");
        // Twice, so the newline separator is exercised rather than the file
        // happening to hold one record.
        writer.append(&record).expect("append");
        let bytes = writer.finish().expect("finish");

        let read: Vec<SerSpeclibElement> = SpeclibReader::new(bytes.as_slice())
            .expect("reader")
            .collect::<Result<_, _>>()
            .expect("every record must parse");

        assert_eq!(read.len(), 2);
        assert_eq!(read[0].precursor.sequence, "PEPTIDEK");
        assert_eq!(
            read[0].elution_group.fragment_mzs,
            element("PEPTIDEK", false, 450.5, &["y1", "b3^2"], &[0.9, 0.4])
                .elution_group
                .fragment_mzs,
            "fragment m/z must survive the round trip"
        );
        let labels: Vec<String> = read[0]
            .elution_group
            .fragment_labels
            .iter()
            .map(|l| l.to_string())
            .collect();
        assert_eq!(labels, vec!["y1".to_string(), "b3^2".to_string()]);
    }

    #[test]
    fn native_extensions_route_to_the_native_reader() {
        use std::path::Path;
        for ext in [
            "lib.ndjson",
            "lib.ndjson.zst",
            "lib.ndjson.zstd",
            "LIB.NDJSON",
        ] {
            assert!(is_native_extension(Path::new(ext)), "{ext} is native");
        }
        // Everything else routes to the timsquery bridge, which sniffs by
        // content. Claiming one of these here would bypass that.
        for ext in ["lib.speclib", "lib.mzSpecLib.txt", "lib.tsv", "lib.zst"] {
            assert!(
                !is_native_extension(Path::new(ext)),
                "{ext} must not be claimed as a native format"
            );
        }
    }

    /// Compression is decided by the magic number, so the two encodings are
    /// interchangeable regardless of what the file is called.
    #[test]
    fn the_reader_sniffs_zstd_rather_than_trusting_the_name() {
        let plain = b"{\"precursor\":{\"sequence\":\"PEPTIDEK\",\"charge\":2,\"decoy\":false,\
                      \"decoy_group\":0},\"elution_group\":{\"id\":0,\"precursor_mz\":500.0,\
                      \"precursor_labels\":[],\"fragment_mzs\":[175.1],\
                      \"fragment_labels\":[\"y1\"],\"precursor_intensities\":[],\
                      \"fragment_intensities\":[1.0],\"mobility_ook0\":0.9,\
                      \"rt_seconds\":10.0}}\n";

        let from_plain: Vec<SerSpeclibElement> = SpeclibReader::new(&plain[..])
            .expect("uncompressed reader")
            .collect::<Result<_, _>>()
            .expect("plain NDJSON parses");

        let compressed = zstd::encode_all(&plain[..], 3).expect("encode");
        assert!(compressed.starts_with(&ZSTD_MAGIC));
        let from_zstd: Vec<SerSpeclibElement> = SpeclibReader::new(compressed.as_slice())
            .expect("compressed reader")
            .collect::<Result<_, _>>()
            .expect("zstd NDJSON parses");

        assert_eq!(from_plain.len(), 1);
        assert_eq!(
            from_plain[0].precursor.sequence,
            from_zstd[0].precursor.sequence
        );
    }

    /// A leftover `.msgpack.zst` must say what happened, not "invalid UTF-8".
    #[test]
    fn msgpack_libraries_report_the_format_removal() {
        let err = Speclib::from_file(
            Path::new("/nonexistent/lib.msgpack.zst"),
            crate::models::DecoyPolicy::default(),
        )
        .expect_err("msgpack is no longer supported");
        match err {
            LibraryReadingError::UnsupportedFormat { message } => {
                assert!(message.contains("msgpack"), "{message}");
                assert!(message.contains("speclib_build_cli"), "{message}");
            }
            other => panic!("expected UnsupportedFormat, got {other:?}"),
        }
    }

    /// `Speclib` is now a type alias for `ReferenceLibrary` (which collapsed
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
        use timsquery::models::capabilities::DECOY_CH2_OFFSET_DA;

        for tgt in 0..lib.geom.n_rows() as u32 {
            let target = RefQuery::new(lib, tgt, 0);
            let plus = RefQuery::new(lib, tgt, 1);
            let minus = RefQuery::new(lib, tgt, 2);

            let target_mz = target.mono_precursor_mz();
            let charge = target.precursor_charge();
            let mz_shift = DECOY_CH2_OFFSET_DA / charge as f64;

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
        assert!(!lib.is_empty(), "library should have entries");

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
    /// seal gate downgrades `LazyMassShift -> Passthrough`: the arena is
    /// 1:1 with the stored rows (no synthetic mass-shift expansion). Proves the
    /// native path produces a lazy `ReferenceLibrary` with the right length, target/
    /// decoy flags, and per-fragment reference intensities.
    #[test]
    fn native_ndjson_load_builds_lazy_arena() {
        use crate::data_sources::reference_library::ScoredIdentity;

        let target = element("PEPTIDEK", false, 500.0, &["y1", "y2"], &[0.8, 0.3]);
        let decoy = element("KEDITPEP", true, 500.0, &["y1", "y2"], &[0.6, 0.4]);

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
    /// rejected by both the fast byte-walk parser and the mzcore fallback), so
    /// the gate must report `!parsable_sequences()`. This is the inverse of
    /// `test_diann_tsv_parsable_gate`, and the only test of the OFF branch after
    /// the materialized `test_parse_gate_off_on_poisoned_row` was removed.
    #[test]
    fn from_file_native_ndjson_poisoned_row_disables_sequence_features() {
        let good = element("PEPTIDEK", false, 500.0, &["y1", "y2"], &[0.8, 0.3]);
        // Unparseable modified sequence: `!` is rejected by parse_sequence_fast
        // (`_ => return None`) and by the mzcore pro_forma fallback.
        let poisoned = element("GARBAGE!!!", false, 600.0, &["y1", "y2"], &[0.7, 0.4]);

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
