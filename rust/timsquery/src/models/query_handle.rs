use std::marker::PhantomData;
use std::ops::Deref;

use crate::models::QueryCollection;
use crate::models::capabilities::{
    DecoyStrategy,
    IsotopeStrategy,
};
use crate::traits::{
    DecoyShift,
    KeyLike,
    QueryGeom,
};
use crate::utils::constants::NEUTRON_MASS;

/// Flyweight handle over a `QueryCollection` arena: `lib` borrows (or owns via
/// `Arc`) the arena, `handle` packs the target index and decoy variant. No
/// decoy geometry is stored anywhere; variants 1/2 compute a ±CH2 mass shift
/// on the fly from `LibCapabilities::decoys`.
#[derive(Debug, Clone, Copy)]
pub struct Query<Lib, L> {
    lib: Lib,
    handle: u64,
    _label: PhantomData<L>,
}

pub type QueryRef<'a, L> = Query<&'a QueryCollection<L>, L>;

impl<Lib, L> Query<Lib, L> {
    pub const VARIANT_BITS: u32 = 2;
    const VARIANT_MASK: u64 = 0b11;

    pub fn new(lib: Lib, tgt: u32, variant: u8) -> Self {
        debug_assert!(u64::from(variant) <= Self::VARIANT_MASK);
        Self {
            lib,
            handle: (u64::from(tgt) << Self::VARIANT_BITS) | u64::from(variant),
            _label: PhantomData,
        }
    }

    pub fn target_idx(&self) -> usize {
        (self.handle >> Self::VARIANT_BITS) as usize
    }

    pub fn variant(&self) -> u8 {
        (self.handle & Self::VARIANT_MASK) as u8
    }
}

impl<Lib: Deref<Target = QueryCollection<L>>, L: KeyLike + DecoyShift> Query<Lib, L> {
    pub fn geom(&self) -> &QueryCollection<L> {
        &self.lib
    }

    /// 0.0 for the target; ±offset for the ± decoy variants.
    fn variant_shift(&self) -> f64 {
        let offset = match self.geom().caps.decoys {
            DecoyStrategy::LazyMassShift { offset, .. } => offset,
            _ => 0.0,
        };
        match self.variant() {
            0 => 0.0,
            1 => offset,
            2 => -offset,
            _ => 0.0,
        }
    }

    /// Number of precursor isotopes this item will be scored with, from the
    /// library's isotope strategy.
    fn n_isotopes(&self) -> u8 {
        match self.geom().caps.isotopes {
            IsotopeStrategy::FromComposition { n_isotopes } => n_isotopes,
        }
    }

    /// Spacing between adjacent isotope peaks in m/z for this precursor's charge.
    fn isotope_step(&self) -> f64 {
        let charge = self.geom().charge[self.target_idx()] as f64;
        NEUTRON_MASS / charge
    }
}

impl<Lib: Deref<Target = QueryCollection<L>>, L: KeyLike + DecoyShift> QueryGeom for Query<Lib, L> {
    type Label = L;

    fn id(&self) -> u32 {
        self.target_idx() as u32 // positional library_id
    }

    fn mono_precursor_mz(&self) -> f64 {
        let tgt = self.target_idx();
        let charge = self.geom().charge[tgt] as f64;
        self.geom().precursor_mz[tgt] + self.variant_shift() / charge
    }

    fn precursor_charge(&self) -> u8 {
        self.geom().charge[self.target_idx()]
    }

    fn rt_seconds(&self) -> f32 {
        self.geom().rt_seconds[self.target_idx()]
    }

    fn mobility_ook0(&self) -> f32 {
        self.geom().mobility[self.target_idx()]
    }

    fn get_precursor_mz_limits(&self) -> (f64, f64) {
        // Span the isotope envelope the item is actually scored with, mirroring
        // `TimsElutionGroup::get_precursor_mz_limits`: isotopes are `0..n`, all
        // non-negative, so the min is the mono and the max is the top isotope.
        let mono = self.mono_precursor_mz();
        let top = self.n_isotopes().saturating_sub(1) as f64;
        (mono, mono + top * self.isotope_step())
    }

    fn precursor_count(&self) -> usize {
        self.n_isotopes() as usize
    }

    fn fragment_count(&self) -> usize {
        self.geom().frag_range(self.target_idx()).len()
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        let mono = self.mono_precursor_mz();
        let step = self.isotope_step();
        (0..self.n_isotopes()).map(move |i| (i as i8, mono + i as f64 * step))
    }

    /// Fragment (label, m/z). m/z is the per-variant SHIFTED value: target
    /// (variant 0) returns the stored m/z; a decoy applies `shift/charge` to
    /// ordinal>2 fragments only, reproducing `create_mass_shifted_decoy`.
    /// Returned BY VALUE, so the computed shift needs no backing storage.
    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&L, f64)> {
        let shift = self.variant_shift();
        let r = self.geom().frag_range(self.target_idx());
        let labels = &self.geom().frag_labels[r.clone()];
        let mzs = &self.geom().frag_mzs[r];
        labels
            .iter()
            .zip(mzs.iter())
            .map(move |(lab, &mz)| (lab, lab.decoy_shift_mz(mz, shift)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IonAnnot;
    use crate::models::QueryCollection;
    use crate::models::capabilities::*;
    use crate::traits::QueryGeom;

    fn one_target_lib() -> QueryCollection<IonAnnot> {
        let mut c = QueryCollection::with_capabilities(LibCapabilities {
            sequence_features: SeqFeatureState::Available,
            fragment_features: FragmentFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::LazyMassShift {
                offset: 14.0,
                n_decoys: 2,
            },
        });
        c.push_target(
            654.855,
            2,
            1.0,
            1.0,
            &[
                (IonAnnot::try_from("y1").unwrap(), 200.0), // ordinal 1 -> NOT shifted
                (IonAnnot::try_from("y8").unwrap(), 896.5),
            ], // ordinal 8 -> shifted
            "PEPTIDEK",
            "PEPTIDEK",
            &[],
        );
        c.seal();
        c
    }

    #[test]
    fn target_variant_is_unshifted() {
        let lib = one_target_lib();
        let q = Query::new(&lib, 0, 0);
        assert_eq!(q.id(), 0);
        assert!((q.mono_precursor_mz() - 654.855).abs() < 1e-9);
        let frags: Vec<_> = q.iter_fragments_refs().collect();
        assert!((frags[1].1 - 896.5).abs() < 1e-9);
    }

    #[test]
    fn plus_decoy_shifts_precursor_and_high_ordinal_only() {
        let lib = one_target_lib();
        let q = Query::new(&lib, 0, 1); // +shift
        // precursor: +14.0/2
        assert!((q.mono_precursor_mz() - (654.855 + 14.0 / 2.0)).abs() < 1e-9);
        let frags: Vec<_> = q.iter_fragments_refs().collect();
        // y1 ordinal 1 -> unshifted
        assert!((frags[0].1 - 200.0).abs() < 1e-9);
        // y8 ordinal 8, charge 1 -> +14.0/1
        assert!((frags[1].1 - (896.5 + 14.0 / 1.0)).abs() < 1e-9);
    }

    #[test]
    fn minus_decoy_shifts_negative() {
        let lib = one_target_lib();
        let q = Query::new(&lib, 0, 2); // -shift
        assert!((q.mono_precursor_mz() - (654.855 - 14.0 / 2.0)).abs() < 1e-9);
    }

    #[test]
    fn precursor_mz_limits_span_isotope_envelope() {
        // n_isotopes=3, charge 2, mono 654.855 -> limits span isotopes 0..3,
        // i.e. (mono, mono + 2*NEUTRON_MASS/2).
        let lib = one_target_lib();
        let q = Query::new(&lib, 0, 0);
        let step = NEUTRON_MASS / 2.0;
        let (lo, hi) = q.get_precursor_mz_limits();
        assert!((lo - 654.855).abs() < 1e-9);
        assert!((hi - (654.855 + 2.0 * step)).abs() < 1e-9);

        assert_eq!(q.precursor_count(), 3);

        let precs: Vec<_> = q.iter_precursors().collect();
        assert_eq!(precs.len(), 3);
        assert_eq!(precs[0].0, 0i8);
        assert!((precs[0].1 - 654.855).abs() < 1e-9);
        assert!((precs[2].1 - (654.855 + 2.0 * step)).abs() < 1e-9);
    }

    #[test]
    fn string_flyweight_never_shifts() {
        use std::sync::Arc;
        let mut c: QueryCollection<Arc<str>> =
            QueryCollection::with_capabilities(LibCapabilities {
                sequence_features: SeqFeatureState::Available,
                fragment_features: FragmentFeatureState::Available,
                isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
                decoys: DecoyStrategy::None,
            });
        c.push_row(
            500.0,
            2,
            1.0,
            1.0,
            &[(Arc::<str>::from("f"), 300.0)],
            "PEP",
            "PEP",
            &[],
            false,
        );
        c.seal();
        let q = Query::new(&c, 0, 0);
        let frags: Vec<_> = q.iter_fragments_refs().collect();
        assert!((frags[0].1 - 300.0).abs() < 1e-9);
    }

    #[test]
    fn precursor_mz_limits_shift_with_decoy_mono() {
        // +decoy (variant 1) shifts the mono by offset/charge; the envelope
        // limits shift by the same amount and keep the same width.
        let lib = one_target_lib();
        let q = Query::new(&lib, 0, 1);
        let step = NEUTRON_MASS / 2.0;
        let mono = 654.855 + 14.0 / 2.0;
        let (lo, hi) = q.get_precursor_mz_limits();
        assert!((lo - mono).abs() < 1e-9);
        assert!((hi - (mono + 2.0 * step)).abs() < 1e-9);
    }
}
