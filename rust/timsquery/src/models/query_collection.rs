use crate::KeyLike;
use crate::models::capabilities::{
    DecoyStrategy,
    LibCapabilities,
};
use crate::models::query_handle::QueryRef;
use crate::traits::DecoyShift;

#[derive(Debug, Clone)]
pub struct ModDefinition {
    pub token: String, // verbatim, e.g. "[UNIMOD:4]"
    pub mono_delta: f64,
    pub cs_delta: (i16, i16),
}

#[derive(Debug, Clone)]
pub struct QueryCollection<L: KeyLike> {
    pub caps: LibCapabilities,
    // per-target scalars, len = n_rows (no `id` column - library_id = target index)
    pub precursor_mz: Vec<f64>,
    pub charge: Vec<u8>,
    pub rt_seconds: Vec<f32>,
    pub mobility: Vec<f32>,
    // per-row decoy flag (len = n_rows)
    pub is_decoy: Vec<bool>,
    // CSR prefix offsets (n+1)
    pub frag_off: Vec<u32>,
    pub seq_strip_off: Vec<u32>,
    pub seq_mod_off: Vec<u32>,
    pub mod_off: Vec<u32>,
    // fragment arenas (len = total fragments)
    pub frag_labels: Vec<L>,
    pub frag_mzs: Vec<f64>,
    // sequences
    pub seq_strip_blob: String,
    pub seq_mod_blob: String,
    // structured mods
    pub mods: Vec<(u8, u16)>,
    pub mod_registry: Vec<ModDefinition>,
}

impl<L: KeyLike> QueryCollection<L> {
    pub fn with_capabilities(caps: LibCapabilities) -> Self {
        Self {
            caps,
            precursor_mz: Vec::new(),
            charge: Vec::new(),
            rt_seconds: Vec::new(),
            mobility: Vec::new(),
            is_decoy: Vec::new(),
            frag_off: vec![0],
            seq_strip_off: vec![0],
            seq_mod_off: vec![0],
            mod_off: vec![0],
            frag_labels: Vec::new(),
            frag_mzs: Vec::new(),
            seq_strip_blob: String::new(),
            seq_mod_blob: String::new(),
            mods: Vec::new(),
            mod_registry: Vec::new(),
        }
    }

    /// Append one row. `mods` are (position, registry_idx) pairs; the caller
    /// is responsible for having registered the mod tokens in `mod_registry`.
    #[allow(clippy::too_many_arguments)]
    pub fn push_row(
        &mut self,
        precursor_mz: f64,
        charge: u8,
        rt_seconds: f32,
        mobility: f32,
        frags: &[(L, f64)],
        seq_strip: &str,
        seq_mod: &str,
        mods: &[(u8, u16)],
        is_decoy: bool,
    ) {
        self.precursor_mz.push(precursor_mz);
        self.charge.push(charge);
        self.rt_seconds.push(rt_seconds);
        self.mobility.push(mobility);
        for (lab, mz) in frags {
            self.frag_labels.push(lab.clone());
            self.frag_mzs.push(*mz);
        }
        // CSR offsets are u32 (documented ~40x headroom over real arena sizes).
        // If a monster library ever exceeds that, fail loud instead of wrapping
        // an offset and silently corrupting every range past the overflow.
        self.frag_off.push(
            u32::try_from(self.frag_labels.len()).expect("fragment arena exceeds u32 offset range"),
        );
        self.seq_strip_blob.push_str(seq_strip);
        self.seq_strip_off.push(
            u32::try_from(self.seq_strip_blob.len())
                .expect("stripped-seq blob exceeds u32 offset range"),
        );
        self.seq_mod_blob.push_str(seq_mod);
        self.seq_mod_off.push(
            u32::try_from(self.seq_mod_blob.len())
                .expect("modified-seq blob exceeds u32 offset range"),
        );
        self.mods.extend_from_slice(mods);
        self.mod_off
            .push(u32::try_from(self.mods.len()).expect("mods arena exceeds u32 offset range"));
        self.is_decoy.push(is_decoy);
    }

    /// Append one target (defaults `is_decoy = false`).
    #[allow(clippy::too_many_arguments)]
    pub fn push_target(
        &mut self,
        precursor_mz: f64,
        charge: u8,
        rt_seconds: f32,
        mobility: f32,
        frags: &[(L, f64)],
        seq_strip: &str,
        seq_mod: &str,
        mods: &[(u8, u16)],
    ) {
        self.push_row(
            precursor_mz,
            charge,
            rt_seconds,
            mobility,
            frags,
            seq_strip,
            seq_mod,
            mods,
            false,
        );
    }

    /// Number of *physical stored rows*, i.e. how many analytes are held in
    /// memory before decoy expansion. Under `LazyMassShift` every row is a
    /// target; under `Passthrough` the count includes any stored decoy rows.
    /// This is the base for `expanded_len` (the logical, iterator-length count).
    pub fn n_rows(&self) -> usize {
        self.charge.len()
    }

    pub fn frag_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.frag_off[tgt] as usize..self.frag_off[tgt + 1] as usize
    }

    pub fn seq_strip_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.seq_strip_off[tgt] as usize..self.seq_strip_off[tgt + 1] as usize
    }

    pub fn seq_mod_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.seq_mod_off[tgt] as usize..self.seq_mod_off[tgt + 1] as usize
    }

    pub fn mod_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.mod_off[tgt] as usize..self.mod_off[tgt + 1] as usize
    }

    /// Seal after build: enforce the decoy-strategy invariant, then release
    /// excess capacity on every arena.
    ///
    /// `LazyMassShift` requires an all-targets arena (decoys are expressed as an
    /// on-the-fly ±CH2 index transform, never stored). If the library shipped
    /// materialized decoys, downgrade to `Passthrough` so the stored rows are
    /// honored 1:1 instead of being silently re-decoyed.
    pub fn seal(&mut self) {
        if matches!(self.caps.decoys, DecoyStrategy::LazyMassShift { .. })
            && self.is_decoy.iter().any(|&d| d)
        {
            tracing::warn!("library ships decoys; downgrading LazyMassShift -> Passthrough");
            self.caps.decoys = DecoyStrategy::Passthrough;
        }
        // The flyweight packs the decoy variant into `VARIANT_BITS` bits. Fail
        // loud here (where the invariant is established) rather than silently
        // corrupting the packed handle in a release build. `variants_per_row`
        // is inlined because `seal` is not bounded on `DecoyShift`.
        let vpr = match self.caps.decoys {
            DecoyStrategy::LazyMassShift { n_decoys, .. } => n_decoys as usize + 1,
            DecoyStrategy::Passthrough | DecoyStrategy::None => 1,
        };
        assert!(
            vpr <= (1usize << QueryRef::<'_, L>::VARIANT_BITS),
            "variants_per_row {vpr} exceeds the {}-bit decoy-variant packing budget",
            QueryRef::<'_, L>::VARIANT_BITS,
        );
        self.precursor_mz.shrink_to_fit();
        self.charge.shrink_to_fit();
        self.rt_seconds.shrink_to_fit();
        self.mobility.shrink_to_fit();
        self.is_decoy.shrink_to_fit();
        self.frag_off.shrink_to_fit();
        self.seq_strip_off.shrink_to_fit();
        self.seq_mod_off.shrink_to_fit();
        self.mod_off.shrink_to_fit();
        self.frag_labels.shrink_to_fit();
        self.frag_mzs.shrink_to_fit();
        self.seq_strip_blob.shrink_to_fit();
        self.seq_mod_blob.shrink_to_fit();
        self.mods.shrink_to_fit();
        self.mod_registry.shrink_to_fit();
    }
}

/// Decoys as an index transform: the arena stores only targets (under
/// `LazyMassShift`), and expanded/flat indices fan each row out into its
/// target + decoy variants. Bounded on `DecoyShift` so `item_at` can hand out a
/// `QueryRef` (the flyweight that computes decoy geometry on the fly).
impl<L: KeyLike + DecoyShift> QueryCollection<L> {
    /// Variants each stored row expands into: `LazyMassShift` adds `n_decoys`
    /// mass-shifted variants (+1 for the target); `Passthrough`/`None` are 1:1.
    pub fn variants_per_row(&self) -> usize {
        match self.caps.decoys {
            DecoyStrategy::LazyMassShift { n_decoys, .. } => n_decoys as usize + 1,
            DecoyStrategy::Passthrough | DecoyStrategy::None => 1,
        }
    }

    /// Logical count after decoy expansion: the flat iterator length, i.e. how
    /// many analytes are *scored* (`n_rows` * `variants_per_row`). Distinct from
    /// `n_rows`, which is the physical in-memory row count.
    pub fn expanded_len(&self) -> usize {
        self.n_rows() * self.variants_per_row()
    }

    /// Flat `0..expanded_len()` index -> `(row, variant)` flyweight.
    pub fn item_at(&self, flat: usize) -> QueryRef<'_, L> {
        let vpr = self.variants_per_row();
        QueryRef::new(self, (flat / vpr) as u32, (flat % vpr) as u8)
    }

    /// A flat index is a target iff its stored row is a target AND it is the
    /// variant-0 (unshifted) slot.
    pub fn is_target(&self, flat: usize) -> bool {
        let vpr = self.variants_per_row();
        let row = flat / vpr;
        let variant = flat % vpr;
        !self.is_decoy[row] && variant == 0
    }

    /// Target-decoy competition group id: the stored row index.
    pub fn decoy_group(&self, flat: usize) -> u32 {
        let row = flat / self.variants_per_row();
        u32::try_from(row).expect("row index exceeds u32::MAX (library_id/decoy_group are u32)")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IonAnnot;
    use crate::models::capabilities::DecoyStrategy;

    #[test]
    fn lazy_massshift_expands_len_and_flags_targets() {
        let mut c = QueryCollection::with_capabilities(LibCapabilities::default_diann()); // LazyMassShift n=2
        c.push_target(
            500.0,
            2,
            1.0,
            0.8,
            &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            "PEP",
            "PEP",
            &[],
        );
        c.seal();
        assert_eq!(c.variants_per_row(), 3);
        assert_eq!(c.expanded_len(), 3);
        assert!(c.is_target(0)); // variant 0
        assert!(!c.is_target(1)); // +decoy
        assert!(!c.is_target(2)); // -decoy
        assert_eq!(c.decoy_group(1), 0);
    }

    #[test]
    fn passthrough_is_one_variant_per_row_honoring_is_decoy() {
        let mut caps = LibCapabilities::default_diann();
        caps.decoys = DecoyStrategy::Passthrough;
        let mut c = QueryCollection::with_capabilities(caps);
        c.push_row(
            500.0,
            2,
            1.0,
            0.8,
            &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            "PEP",
            "PEP",
            &[],
            false,
        );
        c.push_row(
            510.0,
            2,
            1.0,
            0.8,
            &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            "PEP",
            "PEP",
            &[],
            true,
        );
        c.seal();
        assert_eq!(c.variants_per_row(), 1);
        assert_eq!(c.expanded_len(), 2);
        assert!(c.is_target(0));
        assert!(!c.is_target(1)); // stored decoy row
    }

    #[test]
    fn csr_ranges_recover_per_target_fragments() {
        let mut c = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        // target 0: 2 frags; target 1: 3 frags
        c.push_target(
            // precursor_mz
            500.0,
            // charge
            2,
            // rt
            1.0,
            // mob
            0.8,
            // frags
            &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y4").unwrap(), 400.0),
            ],
            // strip
            "PEPTIDEK",
            // modified
            "PEPTIDEK",
            // mods
            &[],
        );
        c.push_target(
            600.0,
            3,
            2.0,
            0.9,
            &[
                (IonAnnot::try_from("y2").unwrap(), 200.0),
                (IonAnnot::try_from("y5").unwrap(), 500.0),
                (IonAnnot::try_from("y6").unwrap(), 600.0),
            ],
            "SAMPLERK",
            "SAMPLERK",
            &[],
        );
        c.seal();
        assert_eq!(c.n_rows(), 2);
        assert_eq!(c.frag_range(0), 0..2);
        assert_eq!(c.frag_range(1), 2..5);
        assert_eq!(c.frag_mzs[c.frag_range(1)][2], 600.0);
        assert_eq!(&c.seq_strip_blob[c.seq_strip_range(0)], "PEPTIDEK");
        assert_eq!(&c.seq_strip_blob[c.seq_strip_range(1)], "SAMPLERK");
    }

    #[test]
    fn string_labeled_collection_and_is_decoy() {
        use std::sync::Arc;
        let mut c: QueryCollection<Arc<str>> =
            QueryCollection::with_capabilities(LibCapabilities::default_diann());
        c.push_row(
            500.0,
            2,
            1.0,
            0.8,
            &[(Arc::<str>::from("frag_a"), 300.0)],
            "PEP",
            "PEP",
            &[],
            // is_decoy
            false,
        );
        c.push_row(
            600.0,
            2,
            1.0,
            0.8,
            &[(Arc::<str>::from("frag_b"), 400.0)],
            "TIDE",
            "TIDE",
            &[],
            // is_decoy
            true,
        );
        c.seal();
        assert_eq!(c.n_rows(), 2);
        assert_eq!(c.is_decoy, vec![false, true]);
        assert_eq!(
            c.frag_labels[c.frag_range(1)][0],
            Arc::<str>::from("frag_b")
        );
    }
}
