use rand::SeedableRng;
use rand::seq::SliceRandom;
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
    sequence_features: SeqFeatureState,
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

    /// Every scored slot in a deterministic pseudorandom order.
    ///
    /// Phase-1 calibration may stop after a prefix, so its traversal cannot
    /// inherit ordering from the source library (RT, sequence, target/decoy
    /// variant, or otherwise). The arena remains the authority that mints the
    /// opaque flat indices; this method only permutes them.
    pub fn shuffled_flats(&self, seed: u64) -> Vec<FlatIdx> {
        let mut flats: Vec<_> = self.geom.flats().collect();
        flats.shuffle(&mut rand::rngs::StdRng::seed_from_u64(seed));
        flats
    }

    /// Finite library-RT extent across every scored variant.
    pub fn rt_range(&self) -> Option<(f32, f32)> {
        self.iter()
            .map(|q| q.rt_seconds())
            .filter(|rt| rt.is_finite())
            .fold(None, |range, rt| {
                Some(match range {
                    Some((lo, hi)) => (lo.min(rt), hi.max(rt)),
                    None => (rt, rt),
                })
            })
    }

    /// Narrow a label-generic [`TargetTable`] (timsquery's one library funnel)
    /// into the ion-annotated `ReferenceLibrary` timsseek scores against.
    ///
    /// Scoring requires an `Mzpaf` arena with reference fragment intensities.
    /// Rejects string labels, a missing intensity sidecar (for example, plain
    /// target JSON), and nonempty libraries with no annotated fragments.
    /// Readers that supply reference intensities preserve them in the sidecar;
    /// this includes the DIA-NN TSV/Parquet adapters.
    ///
    /// Narrowing only, and crate-internal for that reason: none of the load-time
    /// finalize runs, so `caps.sequence_features` keeps the reader's optimistic
    /// default and sequence-derived features get claimed for rows that will not
    /// parse. [`Self::from_sealed_arena`] is the way in from an arena.
    pub(crate) fn from_arena(arena: TargetTable) -> Result<Self, TargetReadingError> {
        match arena {
            TargetTable::Mzpaf { geom, frag_intens } => {
                let frag_intens =
                    frag_intens.ok_or_else(|| TargetReadingError::UnsupportedFormat {
                        message: "DIA-NN library has no fragment intensities".to_string(),
                    })?;
                // An arena with rows but no fragments is a legitimate arena and
                // an unusable library: every score is driven by fragment
                // intensity, so it would extract nothing for every row and
                // report zero results without ever failing. The likely cause is
                // a file whose peaks carry no annotations at all -- a
                // small-molecule library, say -- read under
                // `UnannotatedPeaks::Skip`, which drops every one of them.
                if geom.n_rows() > 0 && geom.n_fragments() == 0 {
                    return Err(TargetReadingError::UnsupportedFormat {
                        message: format!(
                            "library has {} entries and not one annotated fragment; timsseek \
                             scores annotated peptide fragments, so there is nothing here to \
                             score against",
                            geom.n_rows()
                        ),
                    });
                }
                let sequence_features = geom.capabilities().sequence_features;
                Ok(ReferenceLibrary {
                    geom,
                    frag_intens,
                    sequence_features,
                })
            }
            TargetTable::Str { .. } => Err(TargetReadingError::UnsupportedFormat {
                message: "timsseek requires ion-annotated fragments (mzpaf); got string labels"
                    .to_string(),
            }),
        }
    }
}

/// The public conversion spelling; finalize included, so no `.try_into()`
/// yields a library whose parse gate never ran.
impl TryFrom<TargetTable> for ReferenceLibrary {
    type Error = TargetReadingError;

    fn try_from(arena: TargetTable) -> Result<Self, Self::Error> {
        Self::from_sealed_arena(arena)
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

    /// Which of the three kinds of row this is.
    ///
    /// Target-or-not is the arena's question and is asked there
    /// ([`TargetColumns::is_target_slot`]); every caller that answered it from
    /// one of the two halves alone -- the decoy column or the variant index --
    /// got it wrong for the other half's rows and left the FDR estimated from
    /// nothing. What stays here is only *which kind* of decoy, which timsquery
    /// has no vocabulary for: a non-target in the variant-0 slot is one the
    /// library shipped, and anything else is a variant derived over a target
    /// row.
    fn decoy_marking(&self) -> DecoyMarking {
        let (row, variant) = (self.geom.row(), self.geom.variant());
        if self.lib.geom.is_target_slot(row, variant) {
            DecoyMarking::Target
        } else if variant == 0 {
            DecoyMarking::ReversedDecoy
        } else {
            DecoyMarking::MassShiftedDecoy
        }
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
        let IsotopeStrategy::FromComposition { n_isotopes } = self.lib.geom.capabilities().isotopes;
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

/// Identity accessors used by generic scoring loops alongside `QueryGeom`
/// and `ExpectedIntensity`. `RefQuery` resolves handles without copying row
/// data; `materialize_peptide` creates owned peptide metadata.
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
        self.decoy_marking().is_target()
    }

    /// `raw` is the modified-sequence blob slice; parsing is deferred to
    /// `Peptide::parse` and gated on the whole-library parse check. The modified
    /// (not stripped) form is what the `n_mods` feature reads.
    fn materialize_peptide(&self) -> Peptide {
        let coll = &self.lib.geom;
        Peptide {
            raw: coll.seq_mod(self.geom.row()).into(),
            decoy: self.decoy_marking(),
            sequence_features: self.lib.sequence_features == SeqFeatureState::Available,
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
    /// Narrow a sealed [`TargetTable`] and finish it: decoy reporting, the
    /// whole-library parse gate, the averagine tally.
    ///
    /// The one definition of a finished library. `TargetTable`'s variants and
    /// fields are public, so a caller outside timsseek can assemble an arena
    /// itself, and routing it here is what keeps its `caps.sequence_features`
    /// from standing at the reader's optimistic default -- claiming
    /// sequence-derived features for rows that will not parse. Every route in,
    /// [`Self::from_file`] included, ends here.
    pub(crate) fn from_sealed_arena(arena: TargetTable) -> Result<Self, TargetReadingError> {
        let mut lib = Self::from_arena(arena)?;
        lib.report_decoys();
        lib.gate_sequence_features();
        lib.log_entry_stats();
        Ok(lib)
    }

    /// Read a library of any supported format.
    ///
    /// Every load from disk ends here. The arena arrives from timsquery already
    /// sealed, with the policy's decoy half resolved against the rows the file
    /// turned out to ship, so what is left is the finalize every arena goes
    /// through.
    pub fn from_file(
        path: &Path,
        policy: crate::models::LoadPolicy,
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
        let arena = read_timsquery_library(path, policy)?;
        Self::from_sealed_arena(arena)
    }

    /// Whether every sequence in the library parsed (gates sequence-derived
    /// scoring features). Reads the sealed arena's `sequence_features` state.
    pub fn parsable_sequences(&self) -> bool {
        self.sequence_features == SeqFeatureState::Available
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

        self.sequence_features = if n_unparsable == 0 {
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
            self.sequence_features,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::IonAnnot;
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };

    /// Indices come from the arena; there is no constructor from an integer.
    fn row(lib: &ReferenceLibrary, i: usize) -> RowIdx {
        lib.geom.rows().nth(i).unwrap()
    }

    fn flat(lib: &ReferenceLibrary, i: usize) -> FlatIdx {
        lib.geom.flats().nth(i).unwrap()
    }

    use timsquery::models::capabilities::*;

    fn tiny_ref_lib() -> ReferenceLibrary {
        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
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
            sequence_features: SeqFeatureState::Available,
        }
    }

    #[test]
    fn target_table_narrows_to_reference_library() {
        use timsquery::models::capabilities::TargetCapabilities;
        use timsquery::serde::TargetTable;
        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
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

        let sgeom: TargetColumnsBuilder<std::sync::Arc<str>> =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
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

    #[test]
    fn shuffled_flats_are_a_reproducible_permutation() {
        let lib = tiny_ref_lib();
        let native: Vec<_> = lib.geom.flats().collect();
        let shuffled = lib.shuffled_flats(0x5EED);

        assert_eq!(shuffled, lib.shuffled_flats(0x5EED));
        assert_ne!(shuffled, native, "the test seed must exercise reordering");

        let mut sorted = shuffled;
        sorted.sort_unstable();
        assert_eq!(sorted, native, "shuffle must neither lose nor repeat slots");
    }

    #[test]
    fn rt_range_covers_every_variant_without_minting_indices() {
        let lib = tiny_ref_lib();
        let expected = lib
            .iter()
            .map(|q| q.rt_seconds())
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), rt| {
                (lo.min(rt), hi.max(rt))
            });
        assert_eq!(lib.rt_range(), Some(expected));
    }
}

#[cfg(test)]
mod load_tests {
    use super::*;
    use crate::data_sources::reference_library::{
        ExpectedIntensity,
        RefQuery,
    };
    use timsquery::models::{
        Row,
        TargetCapabilities,
        TargetColumnsBuilder,
    };

    /// A retired format is recognised by name and rejected, whichever way its
    /// compression is spelled.
    ///
    /// Asserted on the variant and the format name rather than on the error's
    /// wording, so the sentence stays editable. Note these paths do not exist:
    /// the check runs before any filesystem access, which is what makes the
    /// diagnostic reachable for a file the user has but the test does not.
    /// A fixture under timsquery's `tests/`, which is where every library
    /// fixture lives: the readers are timsquery's, so the files belong beside
    /// them rather than being copied per consumer.
    fn timsquery_fixture(rel: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("crate dir has a parent")
            .join("timsquery")
            .join("tests")
            .join(rel)
    }

    /// A load that decides only about decoys, leaving everything else default.
    fn deciding_decoys(decoys: crate::models::DecoyPolicy) -> crate::models::LoadPolicy {
        crate::models::LoadPolicy {
            decoys,
            ..Default::default()
        }
    }

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
                crate::models::LoadPolicy::default(),
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
    /// Every format the registry reads, loaded through the one public entry
    /// point, asserting the shape invariant once.
    ///
    /// This replaced seven near-duplicate tests, four of which loaded the same
    /// fixture and asserted the same `len() == 6`. A table makes the differences
    /// between formats the visible part, and a new format one row.
    ///
    /// `flat` is `rows * variants_per_row`, so it also pins the decoy expansion:
    /// a file that ships no decoys resolves to `MassShift` under the default
    /// policy and triples, and one that ships its own stays 1:1.
    #[test]
    fn every_format_loads_to_the_same_shape() {
        use timsquery::traits::QueryGeom;

        struct Case {
            fixture: &'static str,
            rows: usize,
            flat: usize,
            frags_first: usize,
            stored_decoys: usize,
        }

        // Measured, not predicted. A changed reader shows up as a diff here.
        let cases = [
            Case {
                fixture: "diann_io_files/sample_lib.txt",
                rows: 2,
                flat: 6,
                frags_first: 5,
                stored_decoys: 0,
            },
            // The fixture a comment used to exclude as "Skyline format, won't
            // load as DIA-NN". Its header is canonical DIA-NN, timsquery reads
            // it as DIA-NN, and it is the only fixture here carrying a
            // per-precursor `transition_group_id`.
            Case {
                fixture: "diann_io_files/sample_lib.tsv",
                rows: 1,
                flat: 3,
                frags_first: 9,
                stored_decoys: 0,
            },
            Case {
                fixture: "diann_io_files/sample_pq_speclib.parquet",
                rows: 3,
                flat: 9,
                frags_first: 4,
                stored_decoys: 0,
            },
            Case {
                fixture: "skyline_io_files/sample_transition_list.csv",
                rows: 14,
                flat: 42,
                frags_first: 4,
                stored_decoys: 0,
            },
            Case {
                fixture: "mzspeclib_files/diann_export.mzspeclib.txt",
                rows: 9,
                flat: 27,
                frags_first: 20,
                stored_decoys: 0,
            },
            // The only fixture that ships its own decoys, so the only one that
            // stays 1:1 rather than tripling.
            Case {
                fixture: "mzspeclib_files/target_decoy_attribute_set.mzspeclib.txt",
                rows: 10,
                flat: 10,
                frags_first: 82,
                stored_decoys: 5,
            },
        ];

        for case in cases {
            let path = timsquery_fixture(case.fixture);
            assert!(path.exists(), "missing fixture {}", case.fixture);
            let lib = ReferenceLibrary::from_file(&path, crate::models::LoadPolicy::default())
                .unwrap_or_else(|e| panic!("{} failed to load: {e:?}", case.fixture));

            let what = case.fixture;
            assert_eq!(lib.geom.n_rows(), case.rows, "{what}: stored rows");
            assert_eq!(lib.len(), case.flat, "{what}: scored slots");
            assert_eq!(
                lib.geom.n_stored_decoys(),
                case.stored_decoys,
                "{what}: shipped decoys"
            );

            let first = lib
                .iter()
                .find(|q| q.geom().variant() == 0)
                .expect("a library with rows has a target");
            assert_eq!(
                first.fragment_count(),
                case.frags_first,
                "{what}: fragments on the first target"
            );
            assert_eq!(
                first.expected_precursor_envelope().len(),
                3,
                "{what}: three precursor isotopes"
            );
            assert!(
                !first
                    .iter_expected_fragments()
                    .collect::<Vec<_>>()
                    .is_empty(),
                "{what}: reference intensities threaded through"
            );
            // Every fixture's sequences parse. The mzSpecLib ones spell their
            // modifications by name (`C[U:Carbamidomethyl]`), which the
            // normalizer used to expand into a form mzcore rejects -- turning
            // sequence features off for the whole library.
            assert!(
                lib.parsable_sequences(),
                "{what}: sequence features disabled by an unparsable row"
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

        let speclib = ReferenceLibrary::from_file(&test_file, crate::models::LoadPolicy::default())
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

        let speclib = ReferenceLibrary::from_file(path, crate::models::LoadPolicy::default())
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

    /// One row of a hand-built arena, in the terms a test states it: a sequence,
    /// whether the library shipped it as a decoy, and which competition group it
    /// declares.
    ///
    /// Written as `RowSpec::target(..)` / `RowSpec::decoy(..)` so the decoy side
    /// of a fixture is legible at the call site rather than being a positional
    /// `true`.
    #[derive(Clone, Copy)]
    struct RowSpec {
        sequence: &'static str,
        is_decoy: bool,
        /// `None` leaves the row undeclared, which `decoy_group_code` reads as a
        /// group of one.
        group: Option<&'static str>,
    }

    impl RowSpec {
        fn target(sequence: &'static str) -> Self {
            Self {
                sequence,
                is_decoy: false,
                group: None,
            }
        }

        fn decoy(sequence: &'static str) -> Self {
            Self {
                is_decoy: true,
                ..Self::target(sequence)
            }
        }

        fn in_group(self, group: &'static str) -> Self {
            Self {
                group: Some(group),
                ..self
            }
        }
    }

    /// An arena assembled by hand, one row per spec, two fragments each, sealed
    /// under `policy`.
    ///
    /// This is the shape a caller outside timsseek can build for itself --
    /// `TargetTable`'s variants and fields are public -- so it is what the seam
    /// has to finalize. The policy is a parameter because seal resolves it
    /// against the rows: the same specs under two policies are two different
    /// libraries, which is most of what there is to test here.
    fn hand_assembled_arena(policy: crate::models::DecoyPolicy, rows: &[RowSpec]) -> TargetTable {
        let frags = [
            (IonAnnot::try_from("y1").unwrap(), 300.0),
            (IonAnnot::try_from("y2").unwrap(), 400.0),
        ];
        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        for (i, spec) in rows.iter().enumerate() {
            geom.push_row(Row {
                precursor_mz: 500.0 + 100.0 * i as f64,
                charge: 2,
                rt_seconds: 120.0,
                mobility: 0.75,
                frags: &frags,
                seq_strip: spec.sequence,
                seq_mod: spec.sequence,
                is_decoy: spec.is_decoy,
                decoy_group: spec.group.map(Into::into),
                ..Default::default()
            });
        }

        let geom = geom.seal(policy).expect("fixture ids are usable");
        TargetTable::Mzpaf {
            geom,
            frag_intens: Some(vec![0.8; frags.len() * rows.len()]),
        }
    }

    /// A target and the decoy the library shipped for it, competing.
    fn a_shipped_pair() -> [RowSpec; 2] {
        [
            RowSpec::target("PEPTIDEK").in_group("pair-1"),
            RowSpec::decoy("KEDITPEP").in_group("pair-1"),
        ]
    }

    /// Two targets and nothing else, so a policy that derives decoys has
    /// something to derive them from.
    fn two_targets() -> [RowSpec; 2] {
        [RowSpec::target("PEPTIDEK"), RowSpec::target("PEPTIDER")]
    }

    /// What scoring reads off every slot: the mark and the target flag, which
    /// derive from one another and so must never disagree.
    fn marks_of(lib: &ReferenceLibrary) -> Vec<DecoyMarking> {
        lib.iter()
            .map(|query| {
                let mark = query.materialize_peptide().decoy;
                assert_eq!(
                    ScoredIdentity::is_target(&query),
                    mark.is_target(),
                    "is_target disagrees with the mark it derives from"
                );
                mark
            })
            .collect()
    }

    /// Every decoy state a hand-built arena can reach, and what scoring sees in
    /// it.
    ///
    /// The whole point of the table is that the two inputs are independent:
    /// neither the policy nor the rows the file shipped decides this alone. A
    /// marking rule that reads only the mass-shift variant index calls every
    /// shipped decoy a target; one that reads only the decoy column calls every
    /// derived decoy a target. Both leave the FDR estimated from nothing while
    /// every row still looks plausible, so the marks are spelled out per slot
    /// rather than counted.
    ///
    /// `Force` over rows that ship decoys is the one state missing here: the
    /// drop is `DecoyPolicy::accepts` applied by a reader at push time, not by
    /// `seal`, so it is only reachable through a file --
    /// `forcing_mass_shift_decoys_replaces_the_ones_a_library_shipped` below.
    #[test]
    fn what_scoring_sees_is_the_policy_and_the_shipped_rows_together() {
        use crate::models::DecoyPolicy;
        use DecoyMarking::{
            MassShiftedDecoy as Shifted,
            ReversedDecoy as Shipped,
            Target,
        };

        struct Case {
            what: &'static str,
            rows: [RowSpec; 2],
            policy: DecoyPolicy,
            variants_per_row: usize,
            stored_decoys: usize,
            marks: Vec<DecoyMarking>,
        }

        // Measured, not predicted.
        let cases = [
            Case {
                what: "no shipped decoys, derive if missing",
                rows: two_targets(),
                policy: DecoyPolicy::IfMissing,
                variants_per_row: 3,
                stored_decoys: 0,
                marks: vec![Target, Shifted, Shifted, Target, Shifted, Shifted],
            },
            Case {
                what: "no shipped decoys, derive regardless",
                rows: two_targets(),
                policy: DecoyPolicy::Force,
                variants_per_row: 3,
                stored_decoys: 0,
                marks: vec![Target, Shifted, Shifted, Target, Shifted, Shifted],
            },
            // The state `report_decoys` warns about: nothing shipped, nothing
            // derived, so the FDR has nothing to estimate itself from. Pinned
            // here so a later change has to be a deliberate one -- the failure
            // mode of "fixing" it quietly is q-values computed off an empty
            // decoy side, which look like results.
            Case {
                what: "no shipped decoys, derive nothing",
                rows: two_targets(),
                policy: DecoyPolicy::Never,
                variants_per_row: 1,
                stored_decoys: 0,
                marks: vec![Target, Target],
            },
            Case {
                what: "shipped decoys, derive if missing",
                rows: a_shipped_pair(),
                policy: DecoyPolicy::IfMissing,
                variants_per_row: 1,
                stored_decoys: 1,
                marks: vec![Target, Shipped],
            },
            Case {
                what: "shipped decoys, derive nothing",
                rows: a_shipped_pair(),
                policy: DecoyPolicy::Never,
                variants_per_row: 1,
                stored_decoys: 1,
                marks: vec![Target, Shipped],
            },
        ];

        for case in cases {
            let lib =
                ReferenceLibrary::from_sealed_arena(hand_assembled_arena(case.policy, &case.rows))
                    .expect("an mzpaf arena carrying intensities narrows");

            let what = case.what;
            assert_eq!(lib.geom.n_rows(), 2, "{what}: stored rows");
            assert_eq!(
                lib.geom.variants_per_row(),
                case.variants_per_row,
                "{what}: variants per row"
            );
            assert_eq!(
                lib.geom.n_stored_decoys(),
                case.stored_decoys,
                "{what}: shipped decoys"
            );
            assert_eq!(marks_of(&lib), case.marks, "{what}: what scoring sees");
        }
    }

    /// A shipped decoy that never competes fails exactly like a shipped decoy
    /// that reached scoring as a target: the group is what pairs it with its
    /// target, and a decoy alone in its group is a decoy the FDR cannot use.
    #[test]
    fn a_shipped_decoy_competes_in_the_group_of_the_target_it_was_built_from() {
        let lib = ReferenceLibrary::from_sealed_arena(hand_assembled_arena(
            crate::models::DecoyPolicy::IfMissing,
            &[
                RowSpec::target("PEPTIDEK").in_group("pair-1"),
                RowSpec::decoy("KEDITPEP").in_group("pair-1"),
                RowSpec::target("SAMPLERK").in_group("pair-2"),
                RowSpec::decoy("KRELPMAS").in_group("pair-2"),
            ],
        ))
        .expect("an mzpaf arena carrying intensities narrows");

        let groups: Vec<GroupCode> = lib.iter().map(|q| q.handles().group).collect();
        assert_eq!(groups.len(), 4, "1:1 with a library that ships its decoys");
        assert_eq!(groups[0], groups[1], "a pair competes");
        assert_eq!(groups[2], groups[3], "a pair competes");
        assert_ne!(groups[0], groups[2], "separate pairs do not");
    }

    /// A derived decoy has to compete with the row it was derived from, or the
    /// mass shift buys nothing: the group comes from the row, so every variant
    /// of a row shares it.
    #[test]
    fn a_mass_shift_decoy_competes_in_the_group_of_the_row_it_came_from() {
        let lib = ReferenceLibrary::from_sealed_arena(hand_assembled_arena(
            crate::models::DecoyPolicy::IfMissing,
            &two_targets(),
        ))
        .expect("an mzpaf arena carrying intensities narrows");

        let groups: Vec<GroupCode> = lib.iter().map(|q| q.handles().group).collect();
        assert_eq!(groups.len(), 6, "two rows, three variants each");
        assert!(
            groups[..3].iter().all(|g| *g == groups[0]),
            "a row's variants compete"
        );
        assert!(
            groups[3..].iter().all(|g| *g == groups[3]),
            "a row's variants compete"
        );
        assert_ne!(groups[0], groups[3], "separate rows do not");
    }

    /// `Force` over a library that ships its own decoys, the one state a
    /// hand-built arena cannot reach.
    ///
    /// The drop is `DecoyPolicy::accepts`, applied by the reader as it pushes
    /// rows, so `seal` never sees the decoys at all; a fixture that dropped them
    /// itself would be asserting on the test's own arithmetic. Read through
    /// `from_file` for that reason. timsquery's
    /// `skipping_shipped_decoys_leaves_only_targets` covers the arena side of
    /// the same load; this is what scoring then reads off it.
    #[test]
    fn forcing_mass_shift_decoys_replaces_the_ones_a_library_shipped() {
        let path = timsquery_fixture("mzspeclib_files/target_decoy_attribute_set.mzspeclib.txt");
        let lib =
            ReferenceLibrary::from_file(&path, deciding_decoys(crate::models::DecoyPolicy::Force))
                .expect("the fixture loads");

        assert_eq!(
            lib.geom.n_rows(),
            5,
            "the five targets, without their decoys"
        );
        assert_eq!(lib.geom.n_stored_decoys(), 0, "the shipped decoys are gone");
        assert_eq!(lib.geom.variants_per_row(), 3);

        let marks = marks_of(&lib);
        assert_eq!(marks.len(), 15, "five rows, three variants each");
        assert!(
            marks.chunks(3).all(|row| row
                == [
                    DecoyMarking::Target,
                    DecoyMarking::MassShiftedDecoy,
                    DecoyMarking::MassShiftedDecoy
                ]),
            "every row is a target with two derived decoys, got {marks:?}"
        );
    }

    /// The exclusion `decoy_marking` reads as "variant 0 and not a target means
    /// the library shipped it": a library never has both stored decoys and
    /// mass-shift variants, so no shipped decoy is ever also a derived one.
    ///
    /// Asserted over files rather than hand-built arenas because that is where
    /// the exclusion is actually established -- see
    /// `seal_alone_keeps_shipped_decoys_under_force` below.
    #[test]
    fn a_library_read_from_a_file_never_has_both_stored_decoys_and_derived_ones() {
        use crate::models::DecoyPolicy;

        for fixture in [
            "mzspeclib_files/target_decoy_attribute_set.mzspeclib.txt",
            "mzspeclib_files/diann_export.mzspeclib.txt",
        ] {
            for policy in [
                DecoyPolicy::IfMissing,
                DecoyPolicy::Force,
                DecoyPolicy::Never,
            ] {
                let lib = ReferenceLibrary::from_file(
                    &timsquery_fixture(fixture),
                    deciding_decoys(policy),
                )
                .unwrap_or_else(|e| panic!("{fixture} under {policy}: {e:?}"));
                if lib.geom.variants_per_row() > 1 {
                    assert_eq!(
                        lib.geom.n_stored_decoys(),
                        0,
                        "{fixture} under {policy}: derived decoys over a shipped decoy"
                    );
                }
            }
        }
    }

    /// Where the exclusion above is enforced, stated where it can fail.
    ///
    /// The sealed boundary rejects a reader that forgot to apply `Force` while
    /// its intensity sidecar was still aligned with rows.
    #[test]
    fn seal_rejects_shipped_decoys_under_force() {
        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        for spec in a_shipped_pair() {
            geom.push_row(Row {
                precursor_mz: 500.0,
                charge: 2,
                frags: &[(IonAnnot::try_from("y2").unwrap(), 200.0)],
                seq_strip: spec.sequence,
                seq_mod: spec.sequence,
                is_decoy: spec.is_decoy,
                decoy_group: spec.group.map(Into::into),
                ..Default::default()
            });
        }
        let err = geom
            .seal(crate::models::DecoyPolicy::Force)
            .expect_err("derived decoys over stored decoys are invalid");
        assert!(matches!(
            err,
            timsquery::models::SourceIdError::ForceWithStoredDecoys { count: 1 }
        ));
    }

    /// The OFF branch of the library-scale parse gate, and the only test of it.
    ///
    /// One row parses (`PEPTIDEK`) and one does not (`GARBAGE!!!`, which both
    /// the byte-walk parser and the mzcore fallback reject). Feature
    /// availability is library-scale on purpose, so the one bad row has to
    /// disable sequence features for the good one too: targets and decoys
    /// scored with different features make FDR meaningless.
    ///
    /// Driven through the public seam rather than a file, since the gate is a
    /// property of the load path and not of any format, and the seam is the only
    /// way an out-of-crate caller reaches it.
    #[test]
    fn one_unparsable_sequence_disables_sequence_features_library_wide() {
        let library = ReferenceLibrary::from_sealed_arena(hand_assembled_arena(
            crate::models::DecoyPolicy::default(),
            &[RowSpec::target("PEPTIDEK"), RowSpec::target("GARBAGE!!!")],
        ))
        .expect("an mzpaf arena carrying intensities narrows");

        assert!(
            !library.parsable_sequences(),
            "one unparsable row turns the gate off for the whole library"
        );
    }

    /// The ON branch, so the assertion above pins the gate's reading of the rows
    /// and not a constant the seam always writes.
    #[test]
    fn an_arena_whose_sequences_all_parse_keeps_sequence_features_on() {
        let library = ReferenceLibrary::from_sealed_arena(hand_assembled_arena(
            crate::models::DecoyPolicy::default(),
            &two_targets(),
        ))
        .expect("an mzpaf arena carrying intensities narrows");

        assert!(
            library.parsable_sequences(),
            "every sequence parses, so the whole library keeps its features"
        );
    }

    /// Narrowing alone runs no gate: the arena keeps the reader's optimistic
    /// default, unparsable rows and all.
    ///
    /// The difference between the two constructors, stated where it can fail
    /// rather than only in a doc comment. Also why the narrowing one is
    /// crate-internal.
    #[test]
    fn narrowing_an_arena_leaves_sequence_features_at_the_readers_default() {
        let library = ReferenceLibrary::from_arena(hand_assembled_arena(
            crate::models::DecoyPolicy::default(),
            &[RowSpec::target("PEPTIDEK"), RowSpec::target("GARBAGE!!!")],
        ))
        .expect("an mzpaf arena carrying intensities narrows");

        assert!(
            library.parsable_sequences(),
            "the optimistic default stands until the gate reads the rows"
        );
    }

    /// A library nothing can be scored against is refused at load rather than
    /// searched to zero results. The mzSpecLib specification's own small-molecule
    /// example is the shape that gets here: its peaks carry no annotations, so
    /// `UnannotatedPeaks::Skip` drops every one and the rows arrive with no
    /// fragments at all.
    ///
    /// `Skip` is stated rather than left to the default, because it is the only
    /// policy that produces this arena; see the sibling below.
    ///
    /// Distinct from the missing-sidecar rejection below, which is a library that
    /// has fragments and no intensities for them.
    #[test]
    fn a_library_with_no_annotated_fragments_is_refused() {
        let err = ReferenceLibrary::from_file(
            &timsquery_fixture("mzspeclib_files/small_molecule.mzspeclib.txt"),
            crate::models::LoadPolicy {
                decoys: crate::models::DecoyPolicy::IfMissing,
                unannotated: crate::models::UnannotatedPeaks::Skip,
            },
        )
        .expect_err("a library with nothing to score against is not searchable");

        let TargetReadingError::UnsupportedFormat { message } = err else {
            panic!("expected an unsupported-format refusal, got {err:?}");
        };
        assert!(
            message.contains("not one annotated fragment"),
            "the refusal has to say what is missing: {message}"
        );
    }

    /// The same library under the default policy, which keeps an unannotated
    /// peak at the m/z the file measured: the arena has fragments, so the guard
    /// above has nothing to refuse and the library is searchable.
    ///
    /// Five peaks, one row, and no sequence -- a small molecule has none, so the
    /// parse gate closes over the sequence-derived scores and leaves the rest.
    #[test]
    fn the_default_policy_makes_that_same_library_searchable() {
        let lib = ReferenceLibrary::from_file(
            &timsquery_fixture("mzspeclib_files/small_molecule.mzspeclib.txt"),
            crate::models::LoadPolicy::default(),
        )
        .expect("a library whose peaks were kept has something to score against");

        assert_eq!(lib.geom.n_rows(), 1);
        assert_eq!(lib.geom.n_fragments(), 5);
        assert_eq!(lib.frag_intens.len(), 5, "the sidecar stays parallel");
        assert!(
            !lib.parsable_sequences(),
            "a small molecule has no sequence"
        );
    }

    /// Geometry-only arenas can serve extraction, but scoring requires the
    /// reference-intensity sidecar. Reject a missing sidecar at narrowing.
    #[test]
    fn reference_library_rejects_mzpaf_without_intensities() {
        use timsquery::serde::TargetTable;

        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
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
