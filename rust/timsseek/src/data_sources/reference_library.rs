use smallvec::SmallVec;
use std::sync::Arc;
use timsquery::IonAnnot;
use timsquery::serde::LibraryArena;

use crate::errors::LibraryReadingError;
use timsquery::models::capabilities::{
    IsotopeStrategy,
    SeqFeatureState,
};
use timsquery::models::{
    Query,
    QueryCollection,
};
use timsquery::traits::QueryGeom;
use timsquery::utils::constants::PROTON_MASS;

use crate::fragment_mass::isotope_dist_or_averagine;
use crate::models::DecoyMarking;
use crate::models::sequence::{
    Peptide,
    normalize_to_proforma,
    parse_sequence,
};

#[derive(Debug, Clone)]
pub struct ReferenceLibrary {
    pub geom: QueryCollection<IonAnnot>,
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
    geom: Query<&'a QueryCollection<IonAnnot>, IonAnnot>,
}

impl ReferenceLibrary {
    pub fn len(&self) -> usize {
        self.geom.expanded_len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Maps a flat `0..len()` index to a `(row, variant)` `RefQuery`, delegating
    /// the flat->(row,variant) math to the arena's `split_flat` transform.
    pub fn item_at(&self, flat_idx: usize) -> RefQuery<'_> {
        let (row, variant) = self.geom.split_flat(flat_idx);
        RefQuery::new(self, row, variant)
    }

    pub fn iter(&self) -> impl Iterator<Item = RefQuery<'_>> {
        (0..self.len()).map(move |i| self.item_at(i))
    }

    /// Narrow a label-generic [`LibraryArena`] (timsquery's one library funnel)
    /// into the ion-annotated `ReferenceLibrary` timsseek scores against.
    ///
    /// Only the `Mzpaf` arena carries the ion chemistry (`IonAnnot`) AND the
    /// reference-intensity sidecar timsseek needs. Two rejections:
    /// - `Str` (string labels): no ion chemistry, no intensities.
    /// - `Mzpaf` without the sidecar (`frag_intens: None`): the current
    ///   TSV/parquet bridge output. Scoring is intensity-driven, so a lib with
    ///   no reference intensities is unusable; the `.speclib` reader (the
    ///   workload) always populates `Some`.
    pub fn from_arena(arena: LibraryArena) -> Result<Self, LibraryReadingError> {
        match arena {
            LibraryArena::Mzpaf { geom, frag_intens } => {
                let frag_intens =
                    frag_intens.ok_or_else(|| LibraryReadingError::UnsupportedFormat {
                        message: "DIA-NN library has no fragment intensities".to_string(),
                    })?;
                Ok(ReferenceLibrary { geom, frag_intens })
            }
            LibraryArena::Str { .. } => Err(LibraryReadingError::UnsupportedFormat {
                message: "timsseek requires ion-annotated fragments (mzpaf); got string labels"
                    .to_string(),
            }),
        }
    }
}

impl TryFrom<LibraryArena> for ReferenceLibrary {
    type Error = LibraryReadingError;

    fn try_from(arena: LibraryArena) -> Result<Self, Self::Error> {
        Self::from_arena(arena)
    }
}

impl<'a> RefQuery<'a> {
    pub fn new(lib: &'a ReferenceLibrary, tgt: u32, variant: u8) -> Self {
        Self {
            lib,
            geom: Query::new(&lib.geom, tgt, variant),
        }
    }

    pub fn geom(&self) -> &Query<&'a QueryCollection<IonAnnot>, IonAnnot> {
        &self.geom
    }

    /// Materialize the output identity `Peptide` for this flyweight.
    ///
    /// `raw` is the modified-sequence blob slice. `parsed` is filled by
    /// normalizing the modified sequence to ProForma and parsing it — but ONLY when sequence
    /// features are `Available` (the whole-library parse gate passed at build
    /// time), else `None`. Parsing the modified (not stripped) form preserves
    /// the mod set the `n_mods` feature reads. Lazy decoys are mass-shift
    /// decoys, so any non-target variant is `MassShiftedDecoy`.
    pub fn materialize_peptide_in_group(&self, decoy_group: u32) -> Peptide {
        let tgt = self.geom.target_idx();
        let coll = &self.lib.geom;
        let raw: Arc<str> = coll.seq_mod_blob[coll.seq_mod_range(tgt)].into();
        let parsed = if coll.caps.sequence_features == SeqFeatureState::Available {
            let normalized = normalize_to_proforma(&raw);
            parse_sequence(&normalized)
        } else {
            None
        };
        let decoy = if self.geom.variant() == 0 {
            DecoyMarking::Target
        } else {
            DecoyMarking::MassShiftedDecoy
        };
        Peptide {
            raw,
            parsed,
            decoy,
            decoy_group,
        }
    }
}

impl<'a> ExpectedIntensity for RefQuery<'a> {
    fn iter_expected_fragments(&self) -> impl Iterator<Item = (IonAnnot, f32)> {
        let tgt = self.geom.target_idx();
        let r = self.lib.geom.frag_range(tgt);
        let labels = &self.lib.geom.frag_labels[r.clone()];
        let intens = &self.lib.frag_intens[r];
        labels.iter().zip(intens.iter()).map(|(&lab, &i)| (lab, i))
    }

    fn expected_precursor_envelope(&self) -> SmallVec<[(i8, f32); 3]> {
        let tgt = self.geom.target_idx();
        let IsotopeStrategy::FromComposition { n_isotopes } = self.lib.geom.caps.isotopes;
        let seq = &self.lib.geom.seq_strip_blob[self.lib.geom.seq_strip_range(tgt)];
        let charge = self.lib.geom.charge[tgt] as f64;
        let neutral = self.lib.geom.precursor_mz[tgt] * charge - charge * PROTON_MASS;
        let (_src, env) = isotope_dist_or_averagine(seq, neutral);
        (0..n_isotopes as usize)
            .map(|i| (i as i8, env[i]))
            .collect()
    }
}

impl<'a> QueryGeom for RefQuery<'a> {
    type Label = IonAnnot;

    fn id(&self) -> u32 {
        self.geom.id()
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

    fn get_precursor_mz_limits(&self) -> (f64, f64) {
        self.geom.get_precursor_mz_limits()
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

/// Arm-neutral identity accessors the scoring loop needs but that are NOT part
/// of `QueryGeom` / `ExpectedIntensity`. Implemented by the `RefQuery`
/// flyweight so the batch scoring loop stays generic (monomorphized,
/// zero-heap) over the concrete type — see
/// `Scorer::{prescore,score_calibrated}_batch_impl`.
pub trait ScoredIdentity {
    /// Positional library id: lazy target index, or the materialized eg's id.
    fn library_id(&self) -> u32;
    /// Whether this item is a target (vs a decoy variant).
    fn is_target(&self) -> bool;
    /// Target-decoy competition group id.
    fn decoy_group(&self) -> u32;
    /// Materialize the output identity `Peptide`.
    fn materialize_peptide(&self) -> Peptide;
}

impl<'a> ScoredIdentity for RefQuery<'a> {
    fn library_id(&self) -> u32 {
        self.geom().id()
    }

    fn is_target(&self) -> bool {
        // Honor stored decoys uniformly (correct under Passthrough): a row is a
        // target only when it is not a stored decoy AND is the variant-0 slot.
        let tgt = self.geom().target_idx();
        !self.lib.geom.is_decoy[tgt] && self.geom().variant() == 0
    }

    fn decoy_group(&self) -> u32 {
        u32::try_from(self.geom().target_idx())
            .expect("target index exceeds u32::MAX (library_id/decoy_group are u32)")
    }

    fn materialize_peptide(&self) -> Peptide {
        let dg = self.decoy_group();
        RefQuery::materialize_peptide_in_group(self, dg)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::IonAnnot;
    use timsquery::models::QueryCollection;
    use timsquery::models::capabilities::*;

    fn tiny_ref_lib() -> ReferenceLibrary {
        let mut geom = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        geom.push_target(
            900.4,
            2,
            1.0,
            1.0,
            &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y8").unwrap(), 800.0),
            ],
            "PEPTIDEK",
            "PEPTIDEK",
            &[],
        );
        geom.seal();
        ReferenceLibrary {
            geom,
            frag_intens: vec![1.0, 0.5],
        }
    }

    #[test]
    fn library_arena_narrows_to_reference_library() {
        use timsquery::models::QueryCollection;
        use timsquery::models::capabilities::LibCapabilities;
        use timsquery::serde::LibraryArena;
        let mut geom = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        geom.push_target(
            900.4,
            2,
            1.0,
            1.0,
            &[(timsquery::IonAnnot::try_from("y3").unwrap(), 300.0)],
            "PEP",
            "PEP",
            &[],
        );
        geom.seal();
        let arena = LibraryArena::Mzpaf {
            geom,
            frag_intens: Some(vec![1.0]),
        };
        let lib: ReferenceLibrary = arena.try_into().unwrap();
        assert_eq!(lib.frag_intens.len(), 1);

        let mut sgeom: QueryCollection<std::sync::Arc<str>> =
            QueryCollection::with_capabilities(LibCapabilities::default_diann());
        sgeom.seal();
        let s = LibraryArena::Str { geom: sgeom };
        assert!(ReferenceLibrary::try_from(s).is_err());
    }

    #[test]
    fn expected_fragments_pair_labels_with_intensities() {
        let lib = tiny_ref_lib();
        let q = RefQuery::new(&lib, 0, 0);
        let pairs: Vec<_> = q.iter_expected_fragments().collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, IonAnnot::try_from("y3").unwrap());
        assert!((pairs[0].1 - 1.0).abs() < 1e-6);
        assert!((pairs[1].1 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn precursor_envelope_is_max_normalized_three_peaks() {
        let lib = tiny_ref_lib();
        let q = RefQuery::new(&lib, 0, 0);
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
        let t: Vec<_> = RefQuery::new(&lib, 0, 0)
            .iter_expected_fragments()
            .collect();
        let d: Vec<_> = RefQuery::new(&lib, 0, 1)
            .iter_expected_fragments()
            .collect();
        assert_eq!(t, d, "intensities are variant-independent");
    }

    #[test]
    fn item_at_scores_reference_library() {
        // Task 9 collapsed `Speclib` to the single `ReferenceLibrary` arena;
        // scoring reads `RefQuery` flyweights via `item_at` (no materialized
        // arm). Variant 0 is the target.
        let lib = tiny_ref_lib();
        let scored = lib.item_at(0);
        assert!(scored.is_target());
        assert_eq!(scored.library_id(), QueryGeom::id(&scored));
        assert_eq!(scored.decoy_group(), 0);
        let got_frags: Vec<(IonAnnot, f32)> = scored.iter_expected_fragments().collect();
        assert_eq!(got_frags.len(), 2);
        assert!((got_frags[0].1 - 1.0).abs() < 1e-6);
        assert!((got_frags[1].1 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn flat_index_maps_to_target_and_variant_in_tpm_order() {
        let lib = tiny_ref_lib(); // 1 target, n_decoys=2 -> len 3
        assert_eq!(lib.len(), 3);
        assert_eq!(lib.item_at(0).geom().variant(), 0);
        assert_eq!(lib.item_at(1).geom().variant(), 1);
        assert_eq!(lib.item_at(2).geom().variant(), 2);
        assert_eq!(lib.item_at(1).geom().target_idx(), 0);
        let all: Vec<_> = lib.iter().map(|q| q.geom().variant()).collect();
        assert_eq!(all, vec![0, 1, 2]);
    }
}
