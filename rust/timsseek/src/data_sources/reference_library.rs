use smallvec::SmallVec;
use timsquery::IonAnnot;
use timsquery::models::capabilities::IsotopeStrategy;
use timsquery::models::{
    Query,
    QueryCollection,
};
use timsquery::traits::QueryGeom;
use timsquery::utils::constants::PROTON_MASS;

use crate::fragment_mass::isotope_dist_or_averagine;

#[derive(Debug, Clone)]
pub struct ReferenceLibrary {
    pub geom: QueryCollection,
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
    geom: Query<&'a QueryCollection>,
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

    pub fn geom(&self) -> &Query<&'a QueryCollection> {
        &self.geom
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
