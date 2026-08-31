use std::marker::PhantomData;
use std::ops::Deref;

use crate::models::TargetColumns;
use crate::models::capabilities::{
    DecoyStrategy,
    IsotopeStrategy,
};
use crate::models::target_columns::{
    FlatIdx,
    RowIdx,
};
use crate::traits::{
    DecoyShift,
    KeyLike,
    QueryGeom,
};
use crate::utils::constants::C13_C12_MASS_DIFF;

/// Flyweight handle over a `TargetColumns` arena: `lib` borrows (or owns via
/// `Arc`) the arena, and `row`/`variant` say which scored slot this is. No
/// decoy geometry is stored anywhere; variants 1/2 compute a ±CH2 mass shift
/// on the fly from `TargetCapabilities::decoys`.
///
/// Built from a [`FlatIdx`] and nothing else -- that already names a
/// `(row, variant)` pair, and [`TargetColumns::split_flat`] owns the encoding.
/// Split once here so the accessors stay field reads.
#[derive(Debug, Clone, Copy)]
pub struct Query<Lib, L> {
    lib: Lib,
    row: RowIdx,
    variant: u8,
    _label: PhantomData<L>,
}

pub type QueryRef<'a, L> = Query<&'a TargetColumns<L>, L>;

impl<Lib, L> Query<Lib, L> {
    /// The stored row this handle points at.
    pub fn row(&self) -> RowIdx {
        self.row
    }

    pub fn variant(&self) -> u8 {
        self.variant
    }
}

impl<Lib: Deref<Target = TargetColumns<L>>, L: KeyLike + DecoyShift> Query<Lib, L> {
    /// Prefer [`TargetColumns::item_at`]; this exists for the owning `Lib`
    /// flavours, which cannot borrow the arena to call it.
    pub fn new(lib: Lib, flat: FlatIdx) -> Self {
        let (row, variant) = lib.split_flat(flat);
        Self {
            lib,
            row,
            variant,
            _label: PhantomData,
        }
    }

    pub fn geom(&self) -> &TargetColumns<L> {
        &self.lib
    }

    /// 0.0 for the target; ±offset for the ± decoy variants.
    fn variant_shift(&self) -> f64 {
        let offset = match self.geom().caps.decoys {
            DecoyStrategy::MassShift { offset } => offset,
            DecoyStrategy::Stored => 0.0,
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
        let charge = self.geom().charge(self.row()) as f64;
        C13_C12_MASS_DIFF / charge
    }
}

impl<Lib: Deref<Target = TargetColumns<L>>, L: KeyLike + DecoyShift> QueryGeom for Query<Lib, L> {
    type Label = L;

    fn source_id(&self) -> Option<crate::models::SourceId<'_>> {
        self.geom().source_id(self.row())
    }

    fn output_id(&self) -> crate::models::SourceId<'_> {
        self.geom().output_id(self.row())
    }

    fn mono_precursor_mz(&self) -> f64 {
        let tgt = self.row();
        let charge = self.geom().charge(tgt) as f64;
        self.geom().precursor_mz(tgt) + self.variant_shift() / charge
    }

    fn precursor_charge(&self) -> u8 {
        self.geom().charge(self.row())
    }

    fn rt_seconds(&self) -> f32 {
        self.geom().rt_seconds(self.row())
    }

    fn mobility_ook0(&self) -> f32 {
        self.geom().mobility(self.row())
    }

    fn precursor_mz_limits(&self) -> (f64, f64) {
        // Span the isotope envelope the item is actually scored with, mirroring
        // `Target::precursor_mz_limits`: isotopes are `0..n`, all
        // non-negative, so the min is the mono and the max is the top isotope.
        let mono = self.mono_precursor_mz();
        let top = self.n_isotopes().saturating_sub(1) as f64;
        (mono, mono + top * self.isotope_step())
    }

    fn precursor_count(&self) -> usize {
        self.n_isotopes() as usize
    }

    fn fragment_count(&self) -> usize {
        self.geom().frag_range(self.row()).len()
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
        let r = self.geom().frag_range(self.row());
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
    use crate::models::TargetColumns;
    use crate::models::capabilities::*;
    use crate::models::source_id::SourceId;
    use crate::models::target_columns::Row;
    use crate::traits::QueryGeom;

    /// Rows come from the arena; there is no constructor from an integer.
    fn first_row<L: KeyLike>(lib: &TargetColumns<L>) -> RowIdx {
        lib.rows().next().expect("fixture has a row")
    }

    /// The flyweight is `Copy` and built once per scored slot, so its size is a
    /// real cost. Two words: the arena pointer, then the row and variant inside
    /// the pointer's alignment padding. This is what makes storing the handles
    /// free relative to packing them into a `u64`.
    #[test]
    fn the_flyweight_stays_two_words() {
        let word = size_of::<usize>();
        assert_eq!(size_of::<QueryRef<'_, IonAnnot>>(), 2 * word);
        assert_eq!(
            size_of::<Query<std::sync::Arc<TargetColumns<IonAnnot>>, IonAnnot>>(),
            2 * word
        );
    }

    fn one_target_lib() -> TargetColumns<IonAnnot> {
        let mut c = TargetColumns::with_capabilities(TargetCapabilities {
            sequence_features: SeqFeatureState::Available,
            fragment_features: FragmentFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::MassShift {
                offset: crate::models::capabilities::DECOY_CH2_OFFSET_DA,
            },
        });
        c.push_row(Row {
            precursor_mz: 654.855,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 1.0,
            frags: &[
                (IonAnnot::try_from("y1").unwrap(), 200.0), // ordinal 1 -> NOT shifted
                (IonAnnot::try_from("y8").unwrap(), 896.5), // ordinal 8 -> shifted
            ],
            seq_strip: "PEPTIDEK",
            seq_mod: "PEPTIDEK",
            ..Default::default()
        });
        c.seal().expect("fixture ids are usable");
        c
    }

    #[test]
    fn target_variant_is_unshifted() {
        let lib = one_target_lib();
        let q = Query::new(&lib, lib.flat_for(first_row(&lib), 0));
        assert_eq!(q.output_id(), SourceId::Numeric(0));
        assert!((q.mono_precursor_mz() - 654.855).abs() < 1e-9);
        let frags: Vec<_> = q.iter_fragments_refs().collect();
        assert!((frags[1].1 - 896.5).abs() < 1e-9);
    }

    #[test]
    fn plus_decoy_shifts_precursor_and_high_ordinal_only() {
        let lib = one_target_lib();
        let q = Query::new(&lib, lib.flat_for(first_row(&lib), 1)); // +shift
        // precursor: +offset/2
        assert!((q.mono_precursor_mz() - (654.855 + DECOY_CH2_OFFSET_DA / 2.0)).abs() < 1e-9);
        let frags: Vec<_> = q.iter_fragments_refs().collect();
        // y1 ordinal 1 -> unshifted
        assert!((frags[0].1 - 200.0).abs() < 1e-9);
        // y8 ordinal 8, charge 1 -> +offset/1
        assert!((frags[1].1 - (896.5 + DECOY_CH2_OFFSET_DA / 1.0)).abs() < 1e-9);
    }

    /// Every slot `MASS_SHIFT_VARIANTS` mints must land somewhere `variant_shift`
    /// actually shifts to. A count above 3 mints slots at variant >= 3, which
    /// falls through to `0.0` -- the target's own m/z, scored as a decoy.
    #[test]
    fn every_minted_variant_slot_gets_its_own_m_over_z() {
        let lib = one_target_lib();
        let row = first_row(&lib);
        assert_eq!(lib.variants_per_row(), MASS_SHIFT_VARIANTS);
        let mzs: Vec<f64> = (0..MASS_SHIFT_VARIANTS as u8)
            .map(|v| Query::new(&lib, lib.flat_for(row, v)).mono_precursor_mz())
            .collect();
        for (i, a) in mzs.iter().enumerate() {
            for b in &mzs[i + 1..] {
                assert!(
                    (a - b).abs() > 1e-9,
                    "two scored slots share m/z {a}, so one is the target scored as a decoy"
                );
            }
        }
    }

    #[test]
    fn minus_decoy_shifts_negative() {
        let lib = one_target_lib();
        let q = Query::new(&lib, lib.flat_for(first_row(&lib), 2)); // -shift
        assert!((q.mono_precursor_mz() - (654.855 - DECOY_CH2_OFFSET_DA / 2.0)).abs() < 1e-9);
    }

    #[test]
    fn precursor_mz_limits_span_isotope_envelope() {
        // n_isotopes=3, charge 2, mono 654.855 -> limits span isotopes 0..3,
        // i.e. (mono, mono + 2*C13_C12_MASS_DIFF/2).
        let lib = one_target_lib();
        let q = Query::new(&lib, lib.flat_for(first_row(&lib), 0));
        let step = C13_C12_MASS_DIFF / 2.0;
        let (lo, hi) = q.precursor_mz_limits();
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
        let mut c: TargetColumns<Arc<str>> = TargetColumns::with_capabilities(TargetCapabilities {
            sequence_features: SeqFeatureState::Available,
            fragment_features: FragmentFeatureState::Available,
            isotopes: IsotopeStrategy::FromComposition { n_isotopes: 3 },
            decoys: DecoyStrategy::Stored,
        });
        c.push_row(Row {
            precursor_mz: 500.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 1.0,
            frags: &[(Arc::<str>::from("f"), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            ..Default::default()
        });
        c.seal().expect("fixture ids are usable");
        let q = Query::new(&c, c.flat_for(first_row(&c), 0));
        let frags: Vec<_> = q.iter_fragments_refs().collect();
        assert!((frags[0].1 - 300.0).abs() < 1e-9);
    }

    #[test]
    fn precursor_mz_limits_shift_with_decoy_mono() {
        // +decoy (variant 1) shifts the mono by offset/charge; the envelope
        // limits shift by the same amount and keep the same width.
        let lib = one_target_lib();
        let q = Query::new(&lib, lib.flat_for(first_row(&lib), 1));
        let step = C13_C12_MASS_DIFF / 2.0;
        let mono = 654.855 + DECOY_CH2_OFFSET_DA / 2.0;
        let (lo, hi) = q.precursor_mz_limits();
        assert!((lo - mono).abs() < 1e-9);
        assert!((hi - (mono + 2.0 * step)).abs() < 1e-9);
    }
}
