use smallvec::SmallVec;
use std::path::Path;
use timsquery::IonAnnot;
use timsquery::serde::{
    TargetTable,
    read_targets_with as read_timsquery_library,
};

use crate::errors::TargetReadingError;
use timsquery::models::capabilities::{
    IsotopeStrategy,
    SeqFeatureState,
};
use timsquery::models::{
    FlatIdx,
    GroupCode,
    Query,
    RowIdx,
    TargetColumns,
};
use timsquery::traits::QueryGeom;
use timsquery::utils::constants::PROTON_MASS;

use crate::fragment_mass::{
    IsotopeSource,
    isotope_dist_or_averagine,
};
use crate::models::DecoyMarking;
use crate::models::sequence::{
    Peptide,
    normalize_to_proforma,
    parse_sequence,
};

#[derive(Debug, Clone)]
pub struct ReferenceLibrary {
    pub geom: TargetColumns<IonAnnot>,
    /// Parallel to `geom.frag_labels` / `geom.frag_mzs`; same `frag_off` ranges.
    pub frag_intens: Vec<f32>,
}

pub trait ExpectedIntensity {
    fn iter_expected_fragments(&self) -> impl Iterator<Item = (IonAnnot, f32)>;
    fn expected_precursor_envelope(&self) -> SmallVec<[(i8, f32); 3]>;
}

/// Flyweight over a `ReferenceLibrary`: `QueryGeom` (geometry, delegated to the
/// arena flyweight) + `ExpectedIntensity` (reference intensities + envelope).
#[derive(Clone, Copy)]
pub struct RefQuery<'a> {
    lib: &'a ReferenceLibrary,
    geom: Query<&'a TargetColumns<IonAnnot>, IonAnnot>,
}

impl ReferenceLibrary {
    pub fn len(&self) -> usize {
        self.geom.expanded_len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Maps a scored slot to its `RefQuery`; the flat->(row, variant) transform
    /// belongs to the arena.
    pub fn item_at(&self, flat: FlatIdx) -> RefQuery<'_> {
        RefQuery::new(self, flat)
    }

    pub fn iter(&self) -> impl Iterator<Item = RefQuery<'_>> {
        self.geom.flats().map(move |f| self.item_at(f))
    }

    /// Scored slots in batches of at most `n`; see `TargetColumns::chunks`.
    pub fn chunks(&self, n: usize) -> impl Iterator<Item = Vec<FlatIdx>> + use<> {
        self.geom.chunks(n)
    }

    /// Narrow a label-generic [`TargetTable`] (timsquery's one library funnel)
    /// into the ion-annotated `ReferenceLibrary` timsseek scores against.
    ///
    /// Only the `Mzpaf` arena carries the ion chemistry (`IonAnnot`) AND the
    /// reference-intensity sidecar timsseek needs. Two rejections:
    /// - `Str` (string labels): no ion chemistry, no intensities.
    /// - `Mzpaf` without the sidecar (`frag_intens: None`): the current
    ///   TSV/parquet bridge output. Scoring is intensity-driven, so a lib with
    ///   no reference intensities is unusable; the `.speclib` reader (the
    ///   workload) always populates `Some`.
    pub fn from_arena(arena: TargetTable) -> Result<Self, TargetReadingError> {
        match arena {
            TargetTable::Mzpaf { geom, frag_intens } => {
                let frag_intens =
                    frag_intens.ok_or_else(|| TargetReadingError::UnsupportedFormat {
                        message: "DIA-NN library has no fragment intensities".to_string(),
                    })?;
                Ok(ReferenceLibrary { geom, frag_intens })
            }
            TargetTable::Str { .. } => Err(TargetReadingError::UnsupportedFormat {
                message: "timsseek requires ion-annotated fragments (mzpaf); got string labels"
                    .to_string(),
            }),
        }
    }
}

impl TryFrom<TargetTable> for ReferenceLibrary {
    type Error = TargetReadingError;

    fn try_from(arena: TargetTable) -> Result<Self, Self::Error> {
        Self::from_arena(arena)
    }
}

impl<'a> RefQuery<'a> {
    pub fn new(lib: &'a ReferenceLibrary, flat: FlatIdx) -> Self {
        Self {
            lib,
            geom: Query::new(&lib.geom, flat),
        }
    }

    pub fn geom(&self) -> &Query<&'a TargetColumns<IonAnnot>, IonAnnot> {
        &self.geom
    }
}

impl<'a> ExpectedIntensity for RefQuery<'a> {
    fn iter_expected_fragments(&self) -> impl Iterator<Item = (IonAnnot, f32)> {
        let tgt = self.geom.row();
        let labels = self.lib.geom.frag_labels(tgt);
        let intens = &self.lib.frag_intens[self.lib.geom.frag_range(tgt)];
        debug_assert_eq!(
            labels.len(),
            intens.len(),
            "reference-intensity sidecar desynced from fragment-label arena"
        );
        labels.iter().zip(intens.iter()).map(|(&lab, &i)| (lab, i))
    }

    fn expected_precursor_envelope(&self) -> SmallVec<[(i8, f32); 3]> {
        let tgt = self.geom.row();
        let IsotopeStrategy::FromComposition { n_isotopes } = self.lib.geom.caps.isotopes;
        let seq = self.lib.geom.seq_strip(tgt);
        let charge = self.lib.geom.charge(tgt) as f64;
        let neutral = self.lib.geom.precursor_mz(tgt) * charge - charge * PROTON_MASS;
        let (_src, env) = isotope_dist_or_averagine(seq, neutral);
        (0..n_isotopes as usize)
            .map(|i| (i as i8, env[i]))
            .collect()
    }
}

impl<'a> QueryGeom for RefQuery<'a> {
    type Label = IonAnnot;

    fn source_id(&self) -> Option<timsquery::models::SourceId<'_>> {
        self.geom.source_id()
    }

    fn output_id(&self) -> timsquery::models::SourceId<'_> {
        self.geom.output_id()
    }

    fn mono_precursor_mz(&self) -> f64 {
        self.geom.mono_precursor_mz()
    }

    fn precursor_charge(&self) -> u8 {
        self.geom.precursor_charge()
    }

    fn rt_seconds(&self) -> f32 {
        self.geom.rt_seconds()
    }

    fn mobility_ook0(&self) -> f32 {
        self.geom.mobility_ook0()
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        self.geom.precursor_mz_limits()
    }

    fn precursor_count(&self) -> usize {
        self.geom.precursor_count()
    }

    fn fragment_count(&self) -> usize {
        self.geom.fragment_count()
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        self.geom.iter_precursors()
    }

    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&IonAnnot, f64)> {
        self.geom.iter_fragments_refs()
    }
}

/// Which row a scored result came from, and which group it competes in.
///
/// The two travel together because the scorer holds the run's raw-data index,
/// not the library arena, so it cannot derive one from the other. Both are
/// opaque, so neither can reach an output file or be confused with an id -- the
/// writer resolves them against the arena at the end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RowHandles {
    pub row: RowIdx,
    pub group: GroupCode,
}

/// Arm-neutral identity accessors the scoring loop needs but that are NOT part
/// of `QueryGeom` / `ExpectedIntensity`. Implemented by the `RefQuery`
/// flyweight so the batch scoring loop stays generic (monomorphized,
/// zero-heap) over the concrete type -- see
/// `Scorer::{prescore,score_calibrated}_batch_impl`.
pub trait ScoredIdentity {
    /// Whether this item is a target (vs a decoy variant).
    fn is_target(&self) -> bool;
    /// The arena handles this result carries onward. See [`RowHandles`].
    fn handles(&self) -> RowHandles;
    /// Materialize the output identity `Peptide`.
    fn materialize_peptide(&self) -> Peptide;
}

impl<'a> ScoredIdentity for RefQuery<'a> {
    fn handles(&self) -> RowHandles {
        let row = self.geom().row();
        RowHandles {
            row,
            group: self.lib.geom.decoy_group_code(row),
        }
    }

    fn is_target(&self) -> bool {
        // Honor stored decoys uniformly (correct under `Stored`): a row is a
        // target only when it is not a stored decoy AND is the variant-0 slot.
        let tgt = self.geom().row();
        !self.lib.geom.is_decoy(tgt) && self.geom().variant() == 0
    }

    /// `raw` is the modified-sequence blob slice; parsing is deferred to
    /// `Peptide::parse` and gated on the whole-library parse check. The modified
    /// (not stripped) form is what the `n_mods` feature reads. Lazy decoys are
    /// mass-shift decoys, so any non-target variant is `MassShiftedDecoy`.
    fn materialize_peptide(&self) -> Peptide {
        let tgt = self.geom.row();
        let coll = &self.lib.geom;
        Peptide {
            raw: coll.seq_mod(tgt).into(),
            decoy: if self.geom.variant() == 0 {
                DecoyMarking::Target
            } else {
                DecoyMarking::MassShiftedDecoy
            },
            sequence_features: coll.caps.sequence_features == SeqFeatureState::Available,
        }
    }
}

/// MessagePack and the native NDJSON format, which nothing reads and nothing
/// writes.
///
/// A user holding one has a file that can only be rebuilt, so it is matched by
/// extension and reported by name. The alternative is a parse error from
/// whichever reader the registry tried last, which says nothing about what to
/// do next.
///
/// Removable once no released version wrote either, which is any release after
/// 0.34.
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

/// Loading: the one path from a path on disk to a scored-against arena.
impl ReferenceLibrary {
    /// Read a library of any supported format.
    ///
    /// Every load ends here. The arena arrives from timsquery already sealed,
    /// with its decoy policy resolved against the rows the file turned out to
    /// ship, so this adds only what timsseek needs on top: the whole-library
    /// parse gate and the averagine tally.
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

        // Bridge to the timsquery reader registry (DIA-NN `.speclib`/TSV/parquet,
        // mzSpecLib, Spectronaut, Skyline, JSON), which returns a label-generic
        // `TargetTable`; narrow it to the ion-annotated arena timsseek scores.
        tracing::info!(
            "Loading library via timsquery format detection: {}",
            path.display()
        );
        let arena = read_timsquery_library(path, decoy_policy)?;
        let mut lib = Self::try_from(arena)?;
        lib.report_decoys();
        lib.gate_sequence_features();
        lib.log_entry_stats();
        Ok(lib)
    }

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

    /// What will actually be scored on the decoy side of the FDR estimate.
    ///
    /// Read off the counts rather than off `caps.decoys`: what matters at load
    /// time is whether anything is there, and the one case worth a warning --
    /// nothing derived and nothing shipped -- is not a strategy.
    fn report_decoys(&self) {
        let n_rows = self.geom.n_rows();
        let n_stored_decoys = self.geom.n_stored_decoys();
        let expanded = self.geom.expanded_len();
        if expanded > n_rows {
            // The arena holds no decoy rows here by construction: `MassShift`
            // survives seal only when nothing was shipped, and `Force` got the
            // reader to drop what was.
            tracing::warn!(
                "Deriving synthetic ±CH2 mass-shift decoys: {n_rows} stored rows -> {expanded} \
                 scored entries",
            );
        } else if n_stored_decoys == 0 {
            tracing::warn!(
                "Library ships no decoys and none will be derived; scoring {n_rows} stored \
                 rows as-is. FDR would be estimated with nothing to estimate it from. Use \
                 --decoy-strategy if-missing to derive mass-shift decoys.",
            );
        } else {
            tracing::info!(
                "Library ships {n_stored_decoys} decoys; scoring {n_rows} stored rows as-is, \
                 deriving none",
            );
        }
    }

    /// Walk every row's modified sequence and set `caps.sequence_features`,
    /// counting averagine isotope fallbacks on the same pass.
    ///
    /// Library-scale and not per-row on purpose: per-row would mean targets and
    /// decoys scored with different features, and then FDR means nothing. The
    /// blob walked here is the one `RefQuery::materialize_peptide_in_group`
    /// parses, so a row that fails here is a row that would fail there.
    fn gate_sequence_features(&mut self) {
        let n_rows = self.geom.n_rows();
        let mut n_unparsable = 0usize;
        let mut first_unparsable: Option<String> = None;
        let mut n_averagine_fallback = 0usize;
        for tgt in self.geom.rows() {
            let modified = self.geom.seq_mod(tgt);
            let normalized = normalize_to_proforma(modified);
            if parse_sequence(&normalized).is_none() {
                n_unparsable += 1;
                first_unparsable.get_or_insert_with(|| modified.to_string());
            }
            let stripped = self.geom.seq_strip(tgt);
            let charge = self.geom.charge(tgt) as f64;
            let neutral_mass = self.geom.precursor_mz(tgt) * charge - charge * PROTON_MASS;
            let (isotope_src, _envelope) = isotope_dist_or_averagine(stripped, neutral_mass);
            if isotope_src == IsotopeSource::Averagine {
                n_averagine_fallback += 1;
            }
        }

        self.geom.caps.sequence_features = if n_unparsable == 0 {
            SeqFeatureState::Available
        } else {
            SeqFeatureState::Unavailable
        };

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
    }

    /// Log a one-line summary of the arena's shape at load time.
    fn log_entry_stats(&self) {
        tracing::info!(
            "Library ready: {} stored rows ({} flat scoring entries, {} fragment slots), \
             sequence_features={:?}",
            self.geom.n_rows(),
            self.len(),
            self.geom.n_fragments(),
            self.geom.caps.sequence_features,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::IonAnnot;
    use timsquery::models::Row;

    /// Indices come from the arena; there is no constructor from an integer.
    fn row(lib: &ReferenceLibrary, i: usize) -> RowIdx {
        lib.geom.rows().nth(i).unwrap()
    }

    fn flat(lib: &ReferenceLibrary, i: usize) -> FlatIdx {
        lib.geom.flats().nth(i).unwrap()
    }

    use timsquery::models::TargetColumns;
    use timsquery::models::capabilities::*;

    fn tiny_ref_lib() -> ReferenceLibrary {
        let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        geom.push_row(Row {
            precursor_mz: 900.4,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 1.0,
            frags: &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y8").unwrap(), 800.0),
            ],
            seq_strip: "PEPTIDEK",
            seq_mod: "PEPTIDEK",
            ..Default::default()
        });
        // `Force` because these tests exercise the decoy variants; timsquery's
        // own readers never derive them.
        let geom = geom
            .seal(crate::models::DecoyPolicy::Force)
            .expect("fixture ids are usable");
        ReferenceLibrary {
            geom,
            frag_intens: vec![1.0, 0.5],
        }
    }

    #[test]
    fn target_table_narrows_to_reference_library() {
        use timsquery::models::TargetColumns;
        use timsquery::models::capabilities::TargetCapabilities;
        use timsquery::serde::TargetTable;
        let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        geom.push_row(Row {
            precursor_mz: 900.4,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 1.0,
            frags: &[(timsquery::IonAnnot::try_from("y3").unwrap(), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            ..Default::default()
        });
        let geom = geom
            .seal(crate::models::DecoyPolicy::Never)
            .expect("fixture ids are usable");
        let arena = TargetTable::Mzpaf {
            geom,
            frag_intens: Some(vec![1.0]),
        };
        let lib: ReferenceLibrary = arena.try_into().unwrap();
        assert_eq!(lib.frag_intens.len(), 1);

        let sgeom: TargetColumns<std::sync::Arc<str>> =
            TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        let sgeom = sgeom
            .seal(crate::models::DecoyPolicy::Never)
            .expect("an empty arena seals");
        let s = TargetTable::Str { geom: sgeom };
        assert!(ReferenceLibrary::try_from(s).is_err());
    }

    #[test]
    fn expected_fragments_pair_labels_with_intensities() {
        let lib = tiny_ref_lib();
        let q = RefQuery::new(&lib, lib.geom.flat_for(row(&lib, 0), 0));
        let pairs: Vec<_> = q.iter_expected_fragments().collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, IonAnnot::try_from("y3").unwrap());
        assert!((pairs[0].1 - 1.0).abs() < 1e-6);
        assert!((pairs[1].1 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn precursor_envelope_is_max_normalized_three_peaks() {
        let lib = tiny_ref_lib();
        let q = RefQuery::new(&lib, lib.geom.flat_for(row(&lib, 0), 0));
        let env = q.expected_precursor_envelope();
        assert_eq!(env.len(), 3);
        // Envelopes are MAX-normalized (base peak = 1.0), matching the
        // composition path (`peptide_isotopes`) so both isotope sources are on
        // one scale. Do NOT assert sum == 1.0.
        let max = env.iter().map(|(_, v)| *v).fold(0.0f32, f32::max);
        assert!(
            (max - 1.0).abs() < 1e-4,
            "base peak must be 1.0, got env {env:?}"
        );
        assert!(env.iter().all(|(_, v)| *v >= 0.0 && *v <= 1.0));
        assert_eq!(env[0].0, 0i8); // isotope indices 0,1,2
    }

    #[test]
    fn decoy_variant_reuses_target_intensities() {
        let lib = tiny_ref_lib();
        let t: Vec<_> = RefQuery::new(&lib, lib.geom.flat_for(row(&lib, 0), 0))
            .iter_expected_fragments()
            .collect();
        let d: Vec<_> = RefQuery::new(&lib, lib.geom.flat_for(row(&lib, 0), 1))
            .iter_expected_fragments()
            .collect();
        assert_eq!(t, d, "intensities are variant-independent");
    }

    #[test]
    fn item_at_scores_reference_library() {
        // Scoring reads `RefQuery` flyweights via `item_at`. Variant 0 is the
        // target.
        let lib = tiny_ref_lib();
        let scored = lib.item_at(flat(&lib, 0));
        assert!(scored.is_target());
        // The id is not on the result; it resolves from the arena via the row
        // the handles carry.
        assert_eq!(
            lib.geom.output_id(scored.handles().row),
            timsquery::models::SourceId::Numeric(0)
        );
        let got_frags: Vec<(IonAnnot, f32)> = scored.iter_expected_fragments().collect();
        assert_eq!(got_frags.len(), 2);
        assert!((got_frags[0].1 - 1.0).abs() < 1e-6);
        assert!((got_frags[1].1 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn flat_index_maps_to_target_and_variant_in_tpm_order() {
        let lib = tiny_ref_lib(); // 1 target, n_decoys=2 -> len 3
        assert_eq!(lib.len(), 3);
        assert_eq!(lib.item_at(flat(&lib, 0)).geom().variant(), 0);
        assert_eq!(lib.item_at(flat(&lib, 1)).geom().variant(), 1);
        assert_eq!(lib.item_at(flat(&lib, 2)).geom().variant(), 2);
        assert_eq!(lib.item_at(flat(&lib, 1)).geom().row(), row(&lib, 0));
        let all: Vec<_> = lib.iter().map(|q| q.geom().variant()).collect();
        assert_eq!(all, vec![0, 1, 2]);
    }
}

#[cfg(test)]
mod load_tests {
    use super::*;
    use crate::data_sources::reference_library::{
        ExpectedIntensity,
        RefQuery,
    };

    /// A retired format is recognised by name and rejected, whichever way its
    /// compression is spelled.
    ///
    /// Asserted on the variant and the format name rather than on the error's
    /// wording, so the sentence stays editable. Note these paths do not exist:
    /// the check runs before any filesystem access, which is what makes the
    /// diagnostic reachable for a file the user has but the test does not.
    #[test]
    fn a_retired_format_is_rejected_by_name() {
        for (name, expected) in [
            ("lib.msgpack", "MessagePack"),
            ("lib.msgpack.zst", "MessagePack"),
            ("lib.msgpack.zstd", "MessagePack"),
            ("lib.ndjson", "the native NDJSON library format"),
            ("lib.ndjson.zst", "the native NDJSON library format"),
        ] {
            assert_eq!(
                retired_format(std::path::Path::new(name)),
                Some(expected),
                "{name}"
            );
            let err = ReferenceLibrary::from_file(
                std::path::Path::new(name),
                crate::models::DecoyPolicy::default(),
            )
            .expect_err("a format nothing writes must not be read");
            assert!(
                matches!(err, TargetReadingError::UnsupportedFormat { .. }),
                "{name}: got {err:?}"
            );
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load DIA-NN TSV library");

        // The sample file has 2 unique precursors with no decoys
        // Should generate 3x flat entries: 2 targets + 4 mass-shift decoys
        assert_eq!(
            speclib.len(),
            6,
            "Expected 6 entries (2 targets + 4 decoys)"
        );

        let lib = &speclib;

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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load Skyline CSV library");

        // Skyline routes through the timsquery bridge (`from_elution_groups`),
        // which now threads the reference intensities through, so it narrows to
        // a lazy `ReferenceLibrary` arena like the DIA-NN formats. No shipped decoys +
        // default IfMissing -> `MassShift`.
        // Fixture has 14 PRTC targets, no decoys -> 14 targets + 28 mass-shift decoys
        assert_eq!(
            speclib.len(),
            42,
            "Expected 42 entries (14 targets + 28 decoys)"
        );

        let lib = &speclib;
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load TXT library");

        // The sample file has 2 unique precursors with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys
        assert_eq!(
            speclib.len(),
            6,
            "Expected 6 entries (2 targets + 4 decoys)"
        );

        let lib = &speclib;
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load Parquet library");

        // The sample parquet file has 3 unique precursors with no decoys
        // Should generate 3x entries: 3 targets + 6 decoys
        assert_eq!(
            speclib.len(),
            9,
            "Expected 9 entries (3 targets + 6 decoys)"
        );

        let lib = &speclib;
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load DIA-NN TSV library");

        let lib = &speclib;

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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load DIA-NN TSV library");

        // Test file has 2 targets with no decoys
        // Should generate 3x entries: 2 targets + 4 decoys (2 per target)
        assert_eq!(
            speclib.len(),
            6,
            "Should have 6 entries (2 targets + 4 decoys)"
        );

        let lib = &speclib;
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load DIA-NN TSV library");

        let lib = &speclib;

        // The offset the arena shifts by, rather than a literal, so this cannot
        // drift from `DecoyPolicy::strategy`.
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

        let speclib =
            ReferenceLibrary::from_file(&test_file, crate::models::DecoyPolicy::default())
                .expect("Failed to load DIA-NN TSV library");

        let lib = &speclib;
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

    /// End-to-end `ReferenceLibrary::from_file` over the real DIA-NN HeLa `.speclib`
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

        let speclib = ReferenceLibrary::from_file(path, crate::models::DecoyPolicy::default())
            .expect("from_file should load the .speclib fixture");

        let lib = &speclib;
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
    /// Written against the gate itself rather than a file, since it is a
    /// property of the load path and not of any format.
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

        let geom = geom
            .seal(crate::models::DecoyPolicy::default())
            .expect("fixture ids are usable");
        let mut library = ReferenceLibrary {
            geom,
            frag_intens: vec![0.8, 0.3, 0.7, 0.4],
        };
        assert!(
            library.parsable_sequences(),
            "the arena starts with the reader's optimistic default, so the \
             assertion below is the gate's doing and not the default's"
        );

        library.gate_sequence_features();

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
        let geom = geom
            .seal(crate::models::DecoyPolicy::Never)
            .expect("fixture ids are usable");
        let arena = TargetTable::Mzpaf {
            geom,
            frag_intens: None,
        };
        assert!(ReferenceLibrary::try_from(arena).is_err());
    }
}
