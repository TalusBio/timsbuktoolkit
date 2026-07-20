use std::ops::Deref;

use crate::IonAnnot;
use crate::models::QueryCollection;
use crate::models::capabilities::DecoyStrategy;
use crate::traits::QueryGeom;

/// Flyweight handle over a `QueryCollection` arena: `lib` borrows (or owns via
/// `Arc`) the arena, `handle` packs the target index and decoy variant. No
/// decoy geometry is stored anywhere; variants 1/2 compute a ±CH2 mass shift
/// on the fly from `LibCapabilities::decoys`.
#[derive(Debug, Clone, Copy)]
pub struct Query<L> {
    lib: L,
    handle: u64,
}

pub type QueryRef<'a> = Query<&'a QueryCollection>;
pub type QueryOwned = Query<std::sync::Arc<QueryCollection>>;

impl<L> Query<L> {
    pub const VARIANT_BITS: u32 = 2;
    const VARIANT_MASK: u64 = 0b11;

    pub fn new(lib: L, tgt: u32, variant: u8) -> Self {
        debug_assert!(u64::from(variant) <= Self::VARIANT_MASK);
        Self {
            lib,
            handle: (u64::from(tgt) << Self::VARIANT_BITS) | u64::from(variant),
        }
    }

    pub fn target_idx(&self) -> usize {
        (self.handle >> Self::VARIANT_BITS) as usize
    }

    pub fn variant(&self) -> u8 {
        (self.handle & Self::VARIANT_MASK) as u8
    }
}

impl<L: Deref<Target = QueryCollection>> Query<L> {
    pub fn geom(&self) -> &QueryCollection {
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
}

impl<L: Deref<Target = QueryCollection>> QueryGeom for Query<L> {
    type Label = IonAnnot;

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
        // Single mono precursor label for the DIA-NN path; the scoring path
        // supplies its own isotope labels. Mirror TimsElutionGroup semantics:
        // limits over non-negative isotopes -> just the mono here.
        let mz = self.mono_precursor_mz();
        (mz, mz)
    }

    fn precursor_count(&self) -> usize {
        1
    }

    fn fragment_count(&self) -> usize {
        self.geom().frag_range(self.target_idx()).len()
    }

    fn iter_precursors(&self) -> impl Iterator<Item = (i8, f64)> {
        std::iter::once((0i8, self.mono_precursor_mz()))
    }

    /// Fragment (label, m/z). m/z is the per-variant SHIFTED value: target
    /// (variant 0) returns the stored m/z; a decoy applies `shift/charge` to
    /// ordinal>2 fragments only, reproducing `create_mass_shifted_decoy`.
    /// Returned BY VALUE, so the computed shift needs no backing storage.
    fn iter_fragments_refs(&self) -> impl Iterator<Item = (&IonAnnot, f64)> {
        let shift = self.variant_shift();
        let r = self.geom().frag_range(self.target_idx());
        let labels = &self.geom().frag_labels[r.clone()];
        let mzs = &self.geom().frag_mzs[r];
        labels.iter().zip(mzs.iter()).map(move |(lab, &mz)| {
            let do_shift = lab.try_get_ordinal().map_or(true, |o| o > 2);
            let out = if do_shift {
                mz + shift / lab.get_charge() as f64
            } else {
                mz
            };
            (lab, out)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IonAnnot;
    use crate::models::QueryCollection;
    use crate::models::capabilities::*;
    use crate::traits::QueryGeom;

    fn one_target_lib() -> QueryCollection {
        let mut c = QueryCollection::with_capabilities(LibCapabilities {
            sequence_features: SeqFeatureState::Available,
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
}
