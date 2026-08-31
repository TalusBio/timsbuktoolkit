use crate::IonAnnot;
use crate::data_sources::reference_library::ReferenceLibrary;
use crate::errors::TargetReadingError;
use crate::fragment_mass::{
    IsotopeSource,
    isotope_dist_or_averagine,
};
use crate::models::sequence::{
    normalize_to_proforma,
    parse_sequence,
};
use std::path::Path;
use timsquery::models::TargetColumns;
use timsquery::models::capabilities::SeqFeatureState;
use timsquery::serde::read_targets as read_timsquery_library;
use timsquery::utils::constants::PROTON_MASS;

/// Summary of a `finalize_reference_library` call, for load-time logging.
#[derive(Debug, Clone, Copy)]
pub struct LoadReport {
    /// Physical stored rows (pre decoy expansion), i.e. `TargetColumns::n_rows`.
    pub n_rows: usize,
    pub n_averagine_fallback: usize,
    /// Rows whose modified sequence neither parser could turn into a
    /// `ParsedSequence`. If this is nonzero, `sequence_features` is unavailable
    /// for the whole library so target and decoy scores use the same features.
    /// Count every failure for the load report.
    pub n_unparsable_sequences: usize,
    pub sequence_features: SeqFeatureState,
}

/// Finalize a freshly-narrowed lazy `ReferenceLibrary` arena: apply the decoy
/// strategy, seal, run the whole-library parse gate + averagine tally, and set
/// `caps.sequence_features`. This is the single shared tail of the DEFAULT
/// `.speclib` load (see `speclib_data_flow.md`) -- the memory-optimized path
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
    mut geom: TargetColumns<IonAnnot>,
    frag_intens: Vec<f32>,
    policy: crate::models::DecoyPolicy,
) -> Result<(ReferenceLibrary, LoadReport), TargetReadingError> {
    let n_stored_decoys = geom.n_stored_decoys();
    geom.caps.decoys = crate::models::map_decoy_strategy(policy, n_stored_decoys > 0);
    geom.seal()?;

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
            DecoyStrategy::None if n_stored_decoys == 0 => {
                tracing::warn!(
                    "Decoy strategy None and the library ships no decoys; scoring {} \
                     stored rows as-is. FDR will be estimated with nothing to \
                     estimate it from. Use --decoy-strategy if-missing to generate \
                     mass-shift decoys.",
                    n_rows
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
    let mut n_unparsable = 0usize;
    let mut first_unparsable: Option<String> = None;
    let mut n_averagine_fallback = 0usize;
    for tgt in geom.rows() {
        let modified = geom.seq_mod(tgt);
        let normalized = normalize_to_proforma(modified);
        if parse_sequence(&normalized).is_none() {
            n_unparsable += 1;
            first_unparsable.get_or_insert_with(|| modified.to_string());
        }
        let stripped = geom.seq_strip(tgt);
        let charge = geom.charge(tgt) as f64;
        let neutral_mass = geom.precursor_mz(tgt) * charge - charge * PROTON_MASS;
        let (isotope_src, _envelope) = isotope_dist_or_averagine(stripped, neutral_mass);
        if isotope_src == IsotopeSource::Averagine {
            n_averagine_fallback += 1;
        }
    }

    let sequence_features = if n_unparsable == 0 {
        SeqFeatureState::Available
    } else {
        SeqFeatureState::Unavailable
    };
    geom.caps.sequence_features = sequence_features;

    if let Some(example) = &first_unparsable {
        tracing::warn!(
            "{}/{} library entries have an unparsable modified sequence, so \
             sequence features are off for the whole library (first: {:?}). \
             Use ProForma, for example `PEPTC[UNIMOD:4]IDEK`. DIA-NN's \
             `(UniMod:n)` form is converted; other modification spellings may \
             fail.",
            n_unparsable,
            n_rows,
            example
        );
    }

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
        n_unparsable_sequences: n_unparsable,
        sequence_features,
    };

    Ok((ReferenceLibrary { geom, frag_intens }, report))
}

/// The spectral library store. Collapsed to the single columnar
/// `ReferenceLibrary` arena representation (the materialized AOS path was
/// deleted in Task 9): both load paths produce a lazy arena, and scoring
/// iterates `RefQuery` flyweights via [`ReferenceLibrary::item_at`].
pub type Speclib = ReferenceLibrary;

/// The formats this project used to write and no longer reads.
///
/// Nothing produces either one after this change, so a user holding one has a
/// file that can only be rebuilt. Matched by extension and reported by name,
/// because the alternative is a parse error from whichever reader the registry
/// happens to try last, which says nothing about what to do next.
fn retired_format(path: &Path) -> Option<&'static str> {
    let name = path.to_string_lossy().to_lowercase();
    let stem = name
        .strip_suffix(".zst")
        .or_else(|| name.strip_suffix(".zstd"))
        .unwrap_or(&name);
    if stem.ends_with(".msgpack") {
        Some("MessagePack")
    } else if stem.ends_with(".ndjson") {
        Some("the native NDJSON library format")
    } else {
        None
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
        self.geom.n_fragments() as f64 / n_rows as f64
    }

    pub fn from_file(
        path: &Path,
        decoy_policy: crate::models::DecoyPolicy,
    ) -> Result<Self, TargetReadingError> {
        if let Some(format) = retired_format(path) {
            return Err(TargetReadingError::UnsupportedFormat {
                message: format!(
                    "{} is {format}, which timsseek no longer reads. Nothing writes it either; \
                     rebuild the library with `timsseek build-library --fasta <FASTA> \
                     -o <NAME>.mzspeclib.txt.gz`.",
                    path.display(),
                ),
            });
        }

        // Terminal source: bridge to the timsquery reader registry (DIA-NN
        // `.speclib`/TSV/parquet, Spectronaut, Skyline, JSON), which returns a
        // label-generic `TargetTable`. One path from here: narrow the arena
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
        let (lib, _report) = finalize_reference_library(geom, frag_intens, decoy_policy)?;
        lib.log_entry_stats();
        Ok(lib)
    }

    /// Log a one-line summary of the lazy arena's shape at load time.
    fn log_entry_stats(&self) {
        tracing::info!(
            "Speclib stats: lazy arena, {} targets ({} flat scoring entries, {} total fragment slots)",
            self.geom.n_rows(),
            self.len(),
            self.geom.n_fragments(),
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

    /// A user upgrading with a MessagePack library in hand gets told what the
    /// file is and how to replace it, rather than a parse error from whichever
    /// reader the registry tried last.
    #[test]
    fn a_retired_format_says_what_it_is_and_how_to_rebuild() {
        for name in [
            "lib.msgpack",
            "lib.msgpack.zst",
            "lib.msgpack.zstd",
            "lib.ndjson",
            "lib.ndjson.zst",
        ] {
            let err = Speclib::from_file(
                std::path::Path::new(name),
                crate::models::DecoyPolicy::default(),
            )
            .expect_err("a format nothing writes must not be read");
            let msg = format!("{err:?}");
            assert!(msg.contains("no longer reads"), "{name}: {msg}");
            assert!(msg.contains("build-library"), "{name}: {msg}");
        }
    }

    /// Every format the registry still handles has to reach it, so the check
    /// above cannot grow into a gate on things that do load.
    #[test]
    fn a_live_format_is_not_mistaken_for_a_retired_one() {
        for name in [
            "lib.speclib",
            "lib.tsv",
            "lib.parquet",
            "lib.mzspeclib.txt",
            "lib.mzspeclib.txt.gz",
        ] {
            assert_eq!(retired_format(std::path::Path::new(name)), None, "{name}");
        }
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
        // for every flat entry (targets AND decoy variants -- the envelope is
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
        // should have exactly 3 flat variants (1 target + 2 decoys) -- this is
        // guaranteed structurally by `ReferenceLibrary::item_at`'s `t,+,-`
        // packing, so assert it directly rather than re-deriving groups.
        let n_target_indices = lib.geom.n_rows();
        assert_eq!(n_target_indices, 2, "Should have 2 unique targets");
        for tgt in lib.geom.rows() {
            let variants: Vec<u8> = (0..3)
                .map(|v| {
                    RefQuery::new(lib, lib.geom.flat_for(tgt, v))
                        .geom()
                        .variant()
                })
                .collect();
            assert_eq!(
                variants,
                vec![0, 1, 2],
                "each row should have exactly 1 target + 2 decoy variants"
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

        for tgt in lib.geom.rows() {
            let target = RefQuery::new(lib, lib.geom.flat_for(tgt, 0));
            let plus = RefQuery::new(lib, lib.geom.flat_for(tgt, 1));
            let minus = RefQuery::new(lib, lib.geom.flat_for(tgt, 2));

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

        let first = lib.item_at(lib.geom.flats().next().unwrap());
        assert!(first.is_target(), "flat index 0 must be a target variant");

        let frags: Vec<_> = first.iter_expected_fragments().collect();
        assert!(
            !frags.is_empty(),
            "target 0 must expose at least one (label, intensity) pair \
             (proves the intensity sidecar threaded through)"
        );
    }

    /// The OFF branch of the library-scale parse gate, and the only test of it.
    ///
    /// One row parses (`PEPTIDEK`) and one does not (`GARBAGE!!!`, which both
    /// the byte-walk parser and the mzcore fallback reject). Feature
    /// availability is library-scale on purpose, so the one bad row has to
    /// disable sequence features for the good one too: targets and decoys
    /// scored with different features make FDR meaningless.
    ///
    /// Written against `finalize_reference_library` rather than a file, since
    /// the gate is a property of that seam and not of any format.
    #[test]
    fn one_unparsable_sequence_disables_sequence_features_library_wide() {
        use timsquery::models::{
            Row,
            TargetCapabilities,
        };

        let frags = [
            (IonAnnot::try_from("y1").unwrap(), 300.0),
            (IonAnnot::try_from("y2").unwrap(), 400.0),
        ];
        let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        for (mz, sequence) in [(500.0, "PEPTIDEK"), (600.0, "GARBAGE!!!")] {
            geom.push_row(Row {
                precursor_mz: mz,
                charge: 2,
                rt_seconds: 120.0,
                mobility: 0.75,
                frags: &frags,
                seq_strip: sequence,
                seq_mod: sequence,
                ..Default::default()
            });
        }

        let (library, report) = finalize_reference_library(
            geom,
            vec![0.8, 0.3, 0.7, 0.4],
            crate::models::DecoyPolicy::default(),
        )
        .expect("an unparseable sequence degrades the library rather than failing the load");

        assert_eq!(report.n_unparsable_sequences, 1);
        assert!(
            !library.parsable_sequences(),
            "one unparsable row turns the gate off for the whole library"
        );
    }

    /// The Mzpaf arena WITHOUT the intensity sidecar (`frag_intens: None`) is
    /// the TSV/parquet/Skyline bridge shape -- scoring is intensity-driven, so
    /// narrowing it to a `ReferenceLibrary` must be an `Err` (the branch that
    /// causes the disclosed DIA-NN/Skyline regression).
    #[test]
    fn reference_library_rejects_mzpaf_without_intensities() {
        use timsquery::models::capabilities::TargetCapabilities;
        use timsquery::models::{
            Row,
            TargetColumns,
        };
        use timsquery::serde::TargetTable;

        let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        geom.push_row(Row {
            precursor_mz: 900.4,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 1.0,
            frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            ..Default::default()
        });
        geom.seal().expect("fixture ids are usable");
        let arena = TargetTable::Mzpaf {
            geom,
            frag_intens: None,
        };
        assert!(ReferenceLibrary::try_from(arena).is_err());
    }
}
