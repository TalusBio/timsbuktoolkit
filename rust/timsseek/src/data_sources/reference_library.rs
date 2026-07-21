use smallvec::SmallVec;
use std::sync::Arc;
use timsquery::IonAnnot;
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
use crate::models::query_item::QueryItemToScore;
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
    fn variants_per_target(&self) -> usize {
        match self.geom.caps.decoys {
            timsquery::models::capabilities::DecoyStrategy::LazyMassShift { n_decoys, .. } => {
                n_decoys as usize + 1
            }
            other => panic!("ReferenceLibrary is lazy-only; got {other:?}"),
        }
    }

    pub fn len(&self) -> usize {
        self.geom.n_targets() * self.variants_per_target()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Maps a flat `0..len()` index to a `(target, variant)` `RefQuery`, in
    /// `t,+,-` order (variant 0,1,2 within each target).
    pub fn item_at(&self, flat_idx: usize) -> RefQuery<'_> {
        let vpt = self.variants_per_target();
        let tgt = (flat_idx / vpt) as u32;
        let variant = (flat_idx % vpt) as u8;
        RefQuery::new(self, tgt, variant)
    }

    pub fn iter(&self) -> impl Iterator<Item = RefQuery<'_>> {
        (0..self.len()).map(move |i| self.item_at(i))
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
    /// `raw` is the modified-sequence blob slice. `parsed` mirrors the
    /// materialized load path (`convert_diann_to_query_item`): normalize the
    /// modified sequence to ProForma and parse it — but ONLY when sequence
    /// features are `Available` (the whole-library parse gate passed at build
    /// time), else `None`. Parsing the modified (not stripped) form preserves
    /// the mod set the `n_mods` feature reads. Lazy decoys are mass-shift
    /// decoys, so any non-target variant is `MassShiftedDecoy`.
    pub fn materialize_peptide(&self, decoy_group: u32) -> Peptide {
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

/// Zero-cost two-arm iterator: unifies the distinct concrete iterator types of
/// the `Lazy` / `Materialized` arms behind one `impl Iterator` with NO heap
/// allocation and NO dyn dispatch. Replaces the `Box<dyn Iterator>` that the
/// `ScoredItem` enum used to return (a per-item alloc storm on the scoring
/// path). Cost is one branch per `next()`.
pub enum EitherIter<L, R> {
    Left(L),
    Right(R),
}

impl<T, L, R> Iterator for EitherIter<L, R>
where
    L: Iterator<Item = T>,
    R: Iterator<Item = T>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
        match self {
            EitherIter::Left(l) => l.next(),
            EitherIter::Right(r) => r.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            EitherIter::Left(l) => l.size_hint(),
            EitherIter::Right(r) => r.size_hint(),
        }
    }
}

/// Arm-neutral identity accessors the scoring loop needs but that are NOT part
/// of `QueryGeom` / `ExpectedIntensity`. Implemented by BOTH concrete
/// flyweights (`RefQuery`, `MatQuery`) so the batch scoring loop can be generic
/// (monomorphized, zero-heap) over the concrete type — see
/// `Scorer::{prescore,score_calibrated}_batch`, which match the `Speclib` arm
/// ONCE and then run a fully static-dispatched loop.
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
        self.geom().variant() == 0
    }

    fn decoy_group(&self) -> u32 {
        u32::try_from(self.geom().target_idx())
            .expect("target index exceeds u32::MAX (library_id/decoy_group are u32)")
    }

    fn materialize_peptide(&self) -> Peptide {
        let dg = self.decoy_group();
        // Inherent `RefQuery::materialize_peptide(&self, u32)` (arity picks it).
        RefQuery::materialize_peptide(self, dg)
    }
}

/// Flyweight over one materialized `QueryItemToScore`. Mirrors `RefQuery` but
/// reads the stored geometry / expected intensities directly instead of
/// deriving them from the columnar arena.
#[derive(Clone, Copy)]
pub struct MatQuery<'a>(pub &'a QueryItemToScore);

impl<'a> ScoredIdentity for MatQuery<'a> {
    fn library_id(&self) -> u32 {
        // QueryGeom::id performs the checked u64 -> u32 conversion.
        self.id()
    }

    fn is_target(&self) -> bool {
        self.0.digest.decoy.is_target()
    }

    fn decoy_group(&self) -> u32 {
        self.0.digest.decoy_group
    }

    fn materialize_peptide(&self) -> Peptide {
        self.0.digest.clone()
    }
}

impl<'a> QueryGeom for MatQuery<'a> {
    type Label = IonAnnot;

    fn id(&self) -> u32 {
        // Delegate to the stored eg's `QueryGeom` impl, which performs the
        // checked u64 -> u32 conversion (fail-loud, no silent truncation).
        QueryGeom::id(&self.0.query)
    }

    fn mono_precursor_mz(&self) -> f64 {
        self.0.query.mono_precursor_mz()
    }

    fn precursor_charge(&self) -> u8 {
        self.0.query.precursor_charge()
    }

    fn rt_seconds(&self) -> f32 {
        self.0.query.rt_seconds()
    }

    fn mobility_ook0(&self) -> f32 {
        self.0.query.mobility_ook0()
    }

    fn get_precursor_mz_limits(&self) -> (f64, f64) {
        self.0.query.get_precursor_mz_limits()
    }

    fn precursor_count(&self) -> usize {
        self.0.query.precursor_count()
    }

    fn fragment_count(&self) -> usize {
        self.0.query.fragment_count()
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        QueryGeom::iter_precursors(&self.0.query)
    }

    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&IonAnnot, f64)> {
        QueryGeom::iter_fragments_refs(&self.0.query)
    }
}

impl<'a> ExpectedIntensity for MatQuery<'a> {
    fn iter_expected_fragments(&self) -> impl Iterator<Item = (IonAnnot, f32)> {
        self.0
            .expected_intensity
            .fragment_intensities
            .iter()
            .map(|&(lab, i)| (lab, i))
    }

    fn expected_precursor_envelope(&self) -> SmallVec<[(i8, f32); 3]> {
        // Materialized arm keeps the STORED envelope verbatim — do NOT recompute
        // averagine here (the lazy arm derives it, this arm already has it).
        self.0
            .expected_intensity
            .precursor_intensities
            .iter()
            .copied()
            .collect()
    }
}

/// Arm-neutral scoring flyweight: one accessor surface over both the lazy
/// columnar arena (`RefQuery`) and a materialized `QueryItemToScore`
/// (`MatQuery`). `Speclib::item_at` hands this out for BOTH arms so the batch
/// scoring / calibration loops never branch on the storage layout.
#[derive(Clone, Copy)]
pub enum ScoredItem<'a> {
    Lazy(RefQuery<'a>),
    Materialized(MatQuery<'a>),
}

impl<'a> ScoredIdentity for ScoredItem<'a> {
    fn library_id(&self) -> u32 {
        match self {
            ScoredItem::Lazy(q) => q.library_id(),
            ScoredItem::Materialized(m) => m.library_id(),
        }
    }

    fn is_target(&self) -> bool {
        match self {
            ScoredItem::Lazy(q) => q.is_target(),
            ScoredItem::Materialized(m) => m.is_target(),
        }
    }

    fn decoy_group(&self) -> u32 {
        match self {
            ScoredItem::Lazy(q) => q.decoy_group(),
            ScoredItem::Materialized(m) => m.decoy_group(),
        }
    }

    fn materialize_peptide(&self) -> Peptide {
        match self {
            ScoredItem::Lazy(q) => ScoredIdentity::materialize_peptide(q),
            ScoredItem::Materialized(m) => ScoredIdentity::materialize_peptide(m),
        }
    }
}

impl<'a> QueryGeom for ScoredItem<'a> {
    type Label = IonAnnot;

    fn id(&self) -> u32 {
        match self {
            ScoredItem::Lazy(q) => q.id(),
            ScoredItem::Materialized(m) => m.id(),
        }
    }

    fn mono_precursor_mz(&self) -> f64 {
        match self {
            ScoredItem::Lazy(q) => q.mono_precursor_mz(),
            ScoredItem::Materialized(m) => m.mono_precursor_mz(),
        }
    }

    fn precursor_charge(&self) -> u8 {
        match self {
            ScoredItem::Lazy(q) => q.precursor_charge(),
            ScoredItem::Materialized(m) => m.precursor_charge(),
        }
    }

    fn rt_seconds(&self) -> f32 {
        match self {
            ScoredItem::Lazy(q) => q.rt_seconds(),
            ScoredItem::Materialized(m) => m.rt_seconds(),
        }
    }

    fn mobility_ook0(&self) -> f32 {
        match self {
            ScoredItem::Lazy(q) => q.mobility_ook0(),
            ScoredItem::Materialized(m) => m.mobility_ook0(),
        }
    }

    fn get_precursor_mz_limits(&self) -> (f64, f64) {
        match self {
            ScoredItem::Lazy(q) => q.get_precursor_mz_limits(),
            ScoredItem::Materialized(m) => m.get_precursor_mz_limits(),
        }
    }

    fn precursor_count(&self) -> usize {
        match self {
            ScoredItem::Lazy(q) => q.precursor_count(),
            ScoredItem::Materialized(m) => m.precursor_count(),
        }
    }

    fn fragment_count(&self) -> usize {
        match self {
            ScoredItem::Lazy(q) => q.fragment_count(),
            ScoredItem::Materialized(m) => m.fragment_count(),
        }
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        // EitherIter unifies the two arms' distinct iterator types with no heap
        // alloc / no dyn dispatch (was `Box<dyn Iterator>`).
        match self {
            ScoredItem::Lazy(q) => EitherIter::Left(q.iter_precursors()),
            ScoredItem::Materialized(m) => EitherIter::Right(m.iter_precursors()),
        }
    }

    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&IonAnnot, f64)> {
        match self {
            ScoredItem::Lazy(q) => EitherIter::Left(q.iter_fragments_refs()),
            ScoredItem::Materialized(m) => EitherIter::Right(m.iter_fragments_refs()),
        }
    }
}

impl<'a> ExpectedIntensity for ScoredItem<'a> {
    fn iter_expected_fragments(&self) -> impl Iterator<Item = (IonAnnot, f32)> {
        match self {
            ScoredItem::Lazy(q) => EitherIter::Left(q.iter_expected_fragments()),
            ScoredItem::Materialized(m) => EitherIter::Right(m.iter_expected_fragments()),
        }
    }

    fn expected_precursor_envelope(&self) -> SmallVec<[(i8, f32); 3]> {
        match self {
            ScoredItem::Lazy(q) => q.expected_precursor_envelope(),
            ScoredItem::Materialized(m) => m.expected_precursor_envelope(),
        }
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
    fn materialized_item_scores() {
        use crate::data_sources::speclib::Speclib;
        use crate::models::query_item::QueryItemToScore;
        use crate::models::sequence::SpeclibMeta;

        let item = QueryItemToScore::sample();
        let expected_mono = item.query.mono_precursor_mz();
        let expected_is_target = item.digest.decoy.is_target();
        let expected_lib_id = QueryGeom::id(&item.query);
        let expected_frags: Vec<(IonAnnot, f32)> = item
            .expected_intensity
            .fragment_intensities
            .iter()
            .copied()
            .collect();

        let lib = Speclib::Materialized {
            elems: vec![item],
            meta: SpeclibMeta::default(),
        };

        // Previously panicked (`unimplemented!`); must now score via item_at.
        let scored = lib.item_at(0);
        assert!((scored.mono_precursor_mz() - expected_mono).abs() < 1e-9);
        assert_eq!(scored.is_target(), expected_is_target);
        assert_eq!(scored.library_id(), expected_lib_id);
        let got_frags: Vec<(IonAnnot, f32)> = scored.iter_expected_fragments().collect();
        assert_eq!(got_frags, expected_frags);
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
