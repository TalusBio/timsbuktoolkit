use crate::KeyLike;
use crate::models::capabilities::{
    DecoyStrategy,
    TargetCapabilities,
};
use crate::models::query_handle::QueryRef;
use crate::models::source_id::{
    OwnedSourceId,
    SourceId,
    SourceIdError,
    SourceIds,
};
use crate::traits::DecoyShift;

pub use index::{
    FlatIdx,
    GroupCode,
    RowIdx,
};

/// The handles the arena hands out: two ways of addressing its memory, and one
/// way of naming a competition group.
///
/// All three are opaque: no constructor from an integer, no accessor yielding
/// one, no `Display`, no `Serialize`. So none of them can be invented by a
/// caller, confused with an id, or written to an output file — a handle can
/// only be obtained from the arena and handed straight back to it.
///
/// Construction is `pub(super)`, i.e. this file, so the arena is the only thing
/// that mints one. `RowIdx` is `pub(in crate::models)` only because the
/// flyweight in `query_handle` packs and unpacks it; unpacking is the last
/// place outside here that builds one.
mod index {
    /// A stored row: `0..n_rows()`.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct RowIdx(u32);

    /// A scored slot after decoy expansion: `0..expanded_len()`. Distinct from
    /// [`RowIdx`] because `variants_per_row` slots map onto one row, so using
    /// one where the other is meant reads a real but wrong row — in range, no
    /// panic, plausible data. The type is what stops that.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct FlatIdx(u32);

    /// Which competition group a row belongs to, as a handle rather than a
    /// value. Rows that compete share one, and that is all a consumer needs:
    /// grouping sorts by it and compares it, and never reads it. Opaque for the
    /// same reason as the indices above — a handle that could be printed would
    /// end up in an output file as though it meant something.
    ///
    /// `Ord` because grouping works by sorting competitors adjacent; any total
    /// order does, since only equality carries meaning.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct GroupCode(u32);

    impl RowIdx {
        pub(in crate::models) fn new(row: u32) -> Self {
            Self(row)
        }

        pub(in crate::models) fn get(self) -> usize {
            self.0 as usize
        }
    }

    impl GroupCode {
        pub(super) fn new(code: u32) -> Self {
            Self(code)
        }

        pub(super) fn get(self) -> usize {
            self.0 as usize
        }
    }

    impl FlatIdx {
        pub(super) fn new(flat: u32) -> Self {
            Self(flat)
        }

        pub(super) fn get(self) -> usize {
            self.0 as usize
        }
    }

    /// `u32::MAX`, not 0. A defaulted index exists only as a placeholder — a
    /// `#[serde(skip)]` field, a test fixture — and 0 is a valid row in every
    /// non-empty arena, so a placeholder that leaked would read row 0 and look
    /// right. `u32::MAX` panics on access instead.
    impl Default for RowIdx {
        fn default() -> Self {
            Self(u32::MAX)
        }
    }

    impl Default for FlatIdx {
        fn default() -> Self {
            Self(u32::MAX)
        }
    }

    /// Same reasoning, and the same placeholder value: a defaulted group is one
    /// nothing has assigned yet. `u32::MAX` makes such rows compete with each
    /// other and with nothing real, which is visible; 0 would silently fold
    /// them into the arena's first group.
    impl Default for GroupCode {
        fn default() -> Self {
            Self(u32::MAX)
        }
    }
}

/// Handles for tests in downstream crates, which score against a real arena in
/// production but assemble results directly in tests.
///
/// Feature-gated, and the feature is only ever enabled by a `dev-dependency`,
/// so "only the arena mints a handle" still holds in every shipped build. The
/// values are dense from 0 because the callers need them *distinct* (a sort key,
/// a competition group) and never resolve them against an arena.
#[cfg(feature = "test-support")]
pub mod test_handles {
    use super::index::{
        GroupCode,
        RowIdx,
    };

    pub fn row(i: u32) -> RowIdx {
        RowIdx::new(i)
    }

    pub fn group(i: u32) -> GroupCode {
        GroupCode::new(i)
    }
}

/// Interned competition groups: `codes[row]` points at `labels`, so rows that
/// compete share one label rather than each storing a copy.
#[derive(Debug, Clone)]
pub struct DecoyGroups {
    codes: Vec<GroupCode>,
    labels: SourceIds,
}

#[derive(Debug, Clone)]
pub struct ModDefinition {
    pub token: String, // verbatim, e.g. "[UNIMOD:4]"
    pub mono_delta: f64,
    pub cs_delta: (i16, i16),
}

#[derive(Debug, Clone)]
pub struct TargetColumns<L: KeyLike> {
    pub caps: TargetCapabilities,
    // per-target scalars, len = n_rows; addressed by `RowIdx`, never by a
    // caller-supplied id (those live in `source_ids`)
    pub(crate) precursor_mz: Vec<f64>,
    pub(crate) charge: Vec<u8>,
    pub(crate) rt_seconds: Vec<f32>,
    pub(crate) mobility: Vec<f32>,
    // per-row decoy flag (len = n_rows)
    pub(crate) is_decoy: Vec<bool>,
    /// What the source file called each row, when it said. Parallel to the row
    /// columns. Minted at [`Self::seal`] for formats that carry no ids, so a
    /// sealed arena always has one per row.
    pub(crate) source_ids: SourceIds,
    /// The competition groups the input declared, interned: a code per row,
    /// plus the labels those codes point at. Rows that compete share a code, so
    /// the label is stored once per group rather than once per row.
    ///
    /// `None` when the input declared none, which is every format today. A
    /// group is then a singleton -- the row competes with its own decoy
    /// variants and nothing else -- so both the code and the label are
    /// derivable from the row and nothing is stored.
    pub(crate) decoy_groups: Option<DecoyGroups>,
    // CSR prefix offsets (n+1)
    pub(crate) frag_off: Vec<u32>,
    pub(crate) seq_strip_off: Vec<u32>,
    pub(crate) seq_mod_off: Vec<u32>,
    pub(crate) mod_off: Vec<u32>,
    // fragment arenas (len = total fragments)
    pub(crate) frag_labels: Vec<L>,
    pub(crate) frag_mzs: Vec<f64>,
    // sequences
    pub(crate) seq_strip_blob: String,
    pub(crate) seq_mod_blob: String,
    // structured mods
    pub(crate) mods: Vec<(u8, u16)>,
    pub(crate) mod_registry: Vec<ModDefinition>,
}

impl<L: KeyLike> TargetColumns<L> {
    pub fn with_capabilities(caps: TargetCapabilities) -> Self {
        Self {
            caps,
            precursor_mz: Vec::new(),
            charge: Vec::new(),
            rt_seconds: Vec::new(),
            mobility: Vec::new(),
            is_decoy: Vec::new(),
            source_ids: SourceIds::default(),
            decoy_groups: None,
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

    /// Attach the ids the file gave its rows, in whichever shape it used.
    /// Call after every `push_row`, before `seal`.
    pub fn set_source_ids<I, S>(&mut self, ids: I) -> Result<(), SourceIdError>
    where
        I: IntoIterator<Item = S>,
        S: Into<OwnedSourceId>,
    {
        let ids: Vec<OwnedSourceId> = ids.into_iter().map(Into::into).collect();
        self.source_ids = SourceIds::owned(ids, self.n_rows())?;
        Ok(())
    }

    /// The rows of this arena, in storage order.
    pub fn rows(&self) -> impl Iterator<Item = RowIdx> + use<L> {
        (0..self.n_rows() as u32).map(RowIdx::new)
    }

    /// Attach the competition groups a file declared, interning them. Rows
    /// sharing a value compete, so unlike source ids these are expected to
    /// repeat -- which is why they are stored as a code per row against a
    /// deduplicated label set rather than a label per row.
    ///
    /// Without this, a group is a singleton and is derived from the row.
    pub fn set_decoy_groups<I, S>(&mut self, groups: I) -> Result<(), SourceIdError>
    where
        I: IntoIterator<Item = S>,
        S: Into<OwnedSourceId>,
    {
        let groups: Vec<OwnedSourceId> = groups.into_iter().map(Into::into).collect();
        if groups.len() != self.n_rows() {
            return Err(SourceIdError::LengthMismatch {
                ids: groups.len(),
                rows: self.n_rows(),
            });
        }

        let mut seen: std::collections::HashMap<OwnedSourceId, GroupCode> =
            std::collections::HashMap::new();
        let mut labels: Vec<OwnedSourceId> = Vec::new();
        let mut codes = Vec::with_capacity(groups.len());
        for group in groups {
            let code = *seen.entry(group.clone()).or_insert_with(|| {
                labels.push(group);
                GroupCode::new((labels.len() - 1) as u32)
            });
            codes.push(code);
        }

        let n_labels = labels.len();
        self.decoy_groups = Some(DecoyGroups {
            codes,
            labels: SourceIds::owned(labels, n_labels)?,
        });
        Ok(())
    }

    /// Which competition group this row is in, as an opaque handle. Rows that
    /// compete share one; nothing else about it is meaningful.
    ///
    /// Derived from the row when the input declared no groups -- a row then
    /// competes only with its own decoy variants, so it is its own group and
    /// there is nothing to store.
    pub fn decoy_group_code(&self, tgt: RowIdx) -> GroupCode {
        match &self.decoy_groups {
            Some(g) => g.codes[tgt.get()],
            None => GroupCode::new(tgt.get() as u32),
        }
    }

    /// The group's id, for output. Resolved here rather than carried alongside
    /// every result, which is the whole point of the code above.
    pub fn decoy_group(&self, tgt: RowIdx) -> SourceId<'_> {
        match &self.decoy_groups {
            Some(g) => g
                .labels
                .get(g.codes[tgt.get()].get())
                .expect("every code indexes a label"),
            None => self.output_id(tgt),
        }
    }

    pub fn source_id(&self, tgt: RowIdx) -> Option<SourceId<'_>> {
        self.source_ids.get(tgt.get())
    }

    /// The id a result for row `tgt` carries. Always a source id: [`Self::seal`]
    /// mints them for formats that carry none, so a row position has no route
    /// into output.
    pub fn output_id(&self, tgt: RowIdx) -> SourceId<'_> {
        self.source_id(tgt)
            .expect("sealed targets have source ids; seal() mints any that are missing")
    }

    pub fn charge(&self, tgt: RowIdx) -> u8 {
        self.charge[tgt.get()]
    }

    pub fn precursor_mz(&self, tgt: RowIdx) -> f64 {
        self.precursor_mz[tgt.get()]
    }

    pub fn rt_seconds(&self, tgt: RowIdx) -> f32 {
        self.rt_seconds[tgt.get()]
    }

    pub fn mobility(&self, tgt: RowIdx) -> f32 {
        self.mobility[tgt.get()]
    }

    pub fn is_decoy(&self, tgt: RowIdx) -> bool {
        self.is_decoy[tgt.get()]
    }

    /// This row's fragment labels and m/z values, positionally paired.
    pub fn frag_labels(&self, tgt: RowIdx) -> &[L] {
        &self.frag_labels[self.frag_range(tgt)]
    }

    pub fn frag_mzs(&self, tgt: RowIdx) -> &[f64] {
        &self.frag_mzs[self.frag_range(tgt)]
    }

    /// This row's stripped (unmodified) sequence.
    pub fn seq_strip(&self, tgt: RowIdx) -> &str {
        &self.seq_strip_blob[self.seq_strip_range(tgt)]
    }

    /// This row's modified sequence, the form sequence features are parsed from.
    pub fn seq_mod(&self, tgt: RowIdx) -> &str {
        &self.seq_mod_blob[self.seq_mod_range(tgt)]
    }

    /// Fragments across every row, for whole-arena statistics.
    pub fn n_fragments(&self) -> usize {
        self.frag_labels.len()
    }

    /// Rows the input itself marked as decoys, as opposed to variants this
    /// arena would generate.
    pub fn n_stored_decoys(&self) -> usize {
        self.is_decoy.iter().filter(|&&d| d).count()
    }

    pub fn frag_range(&self, tgt: RowIdx) -> std::ops::Range<usize> {
        let tgt = tgt.get();
        self.frag_off[tgt] as usize..self.frag_off[tgt + 1] as usize
    }

    pub fn seq_strip_range(&self, tgt: RowIdx) -> std::ops::Range<usize> {
        let tgt = tgt.get();
        self.seq_strip_off[tgt] as usize..self.seq_strip_off[tgt + 1] as usize
    }

    pub fn seq_mod_range(&self, tgt: RowIdx) -> std::ops::Range<usize> {
        let tgt = tgt.get();
        self.seq_mod_off[tgt] as usize..self.seq_mod_off[tgt + 1] as usize
    }

    pub fn mod_range(&self, tgt: RowIdx) -> std::ops::Range<usize> {
        let tgt = tgt.get();
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
        if matches!(self.source_ids, SourceIds::Absent) {
            let n = self.n_rows();
            tracing::warn!(
                "input carries no per-target id; minting {n} self-incremental ids (0..{n}). \
                 Result ids are ours, not the input file's."
            );
            self.source_ids = SourceIds::minted(n);
        }
        // Nothing to build when no groups were declared: a row is then its own
        // group, which `decoy_group_code` derives. Only the case that actually
        // loses information is worth a word.
        if self.decoy_groups.is_none() && matches!(self.caps.decoys, DecoyStrategy::Passthrough) {
            tracing::warn!(
                "library ships its own decoys but declares no competition groups, so each \
                 stored decoy competes alone rather than against its target"
            );
        }
        if matches!(self.caps.decoys, DecoyStrategy::LazyMassShift { .. })
            && self.is_decoy.iter().any(|&d| d)
        {
            tracing::warn!("library ships decoys; downgrading LazyMassShift -> Passthrough");
            self.caps.decoys = DecoyStrategy::Passthrough;
        }
        // The flyweight packs the decoy variant into `VARIANT_BITS` bits. Fail
        // loud here (where the invariant is established) rather than silently
        // corrupting the packed handle in a release build. The free fn is used
        // because `seal` is not bounded on `DecoyShift`.
        let vpr = variants_per_row_for(self.caps.decoys);
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
/// Variants each stored row expands into, from the decoy strategy alone:
/// `LazyMassShift` adds `n_decoys` mass-shifted variants (+1 for the target);
/// `Passthrough`/`None` are 1:1. Free fn (no `DecoyShift` bound) so both the
/// bounded `variants_per_row` and the unbounded `seal` share one definition.
pub fn variants_per_row_for(decoys: DecoyStrategy) -> usize {
    match decoys {
        DecoyStrategy::LazyMassShift { n_decoys, .. } => n_decoys as usize + 1,
        DecoyStrategy::Passthrough | DecoyStrategy::None => 1,
    }
}

impl<L: KeyLike + DecoyShift> TargetColumns<L> {
    /// Variants each stored row expands into (see [`variants_per_row_for`]).
    pub fn variants_per_row(&self) -> usize {
        variants_per_row_for(self.caps.decoys)
    }

    /// Logical count after decoy expansion: the flat iterator length, i.e. how
    /// many analytes are *scored* (`n_rows` * `variants_per_row`). Distinct from
    /// `n_rows`, which is the physical in-memory row count.
    pub fn expanded_len(&self) -> usize {
        self.n_rows() * self.variants_per_row()
    }

    /// Every scored slot, in flat order.
    pub fn flats(&self) -> impl Iterator<Item = FlatIdx> + use<L> {
        (0..self.expanded_len() as u32).map(FlatIdx::new)
    }

    /// The scored slots in batches of at most `n`. The arena owns the split,
    /// so a caller batches work by stating a size rather than by computing
    /// positions — no integer names a slot outside this crate.
    pub fn chunks(&self, n: usize) -> impl Iterator<Item = Vec<FlatIdx>> + use<L> {
        assert!(n > 0, "chunk size must be non-zero");
        let len = self.expanded_len();
        (0..len).step_by(n).map(move |start| {
            (start as u32..(start + n).min(len) as u32)
                .map(FlatIdx::new)
                .collect()
        })
    }

    /// Decompose a scored slot into its stored row and decoy variant. The
    /// single authority for the decoy index-transform encoding, and the only
    /// route from a [`FlatIdx`] to a [`RowIdx`] — which is what keeps the two
    /// index spaces from being used interchangeably.
    pub fn split_flat(&self, flat: FlatIdx) -> (RowIdx, u8) {
        let vpr = self.variants_per_row();
        let flat = flat.get();
        let row = u32::try_from(flat / vpr).expect("row index exceeds u32::MAX");
        (RowIdx::new(row), (flat % vpr) as u8)
    }

    /// Scored slot -> `(row, variant)` flyweight.
    pub fn item_at(&self, flat: FlatIdx) -> QueryRef<'_, L> {
        let (row, variant) = self.split_flat(flat);
        QueryRef::new(self, row, variant)
    }

    /// A slot is a target iff its stored row is a target AND it is the
    /// variant-0 (unshifted) slot.
    pub fn is_target(&self, flat: FlatIdx) -> bool {
        let (row, variant) = self.split_flat(flat);
        !self.is_decoy[row.get()] && variant == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IonAnnot;
    use crate::models::capabilities::DecoyStrategy;

    fn two_rows() -> TargetColumns<IonAnnot> {
        let mut c = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        for _ in 0..2 {
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
        }
        c
    }

    /// Declared groups are interned, so two rows in one group share a code and
    /// the label is stored once rather than per row.
    #[test]
    fn rows_in_one_declared_group_share_a_code() {
        let mut c = two_rows();
        c.set_decoy_groups(["g1", "g1"]).unwrap();
        c.seal();

        let (a, b) = (c.rows().next().unwrap(), c.rows().nth(1).unwrap());
        assert_eq!(c.decoy_group_code(a), c.decoy_group_code(b));
        assert_eq!(c.decoy_group(a), SourceId::Text("g1"));
    }

    /// Without declared groups a row is its own group, derived rather than
    /// stored -- so nothing is allocated and the label is the row's own id.
    #[test]
    fn undeclared_groups_are_derived_not_stored() {
        let mut c = two_rows();
        c.seal();

        let (a, b) = (c.rows().next().unwrap(), c.rows().nth(1).unwrap());
        assert!(c.decoy_groups.is_none(), "nothing stored for minted groups");
        assert_ne!(
            c.decoy_group_code(a),
            c.decoy_group_code(b),
            "each row competes alone"
        );
        assert_eq!(c.decoy_group(a), c.output_id(a));
    }

    /// `output_id` expects a source id on every row; this is what makes that
    /// safe for the formats that carry none.
    #[test]
    fn seal_mints_ids_when_the_input_carried_none() {
        let mut c = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
        for _ in 0..2 {
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
        }
        assert_eq!(c.source_id(RowIdx::new(0)), None);
        c.seal();
        assert_eq!(c.output_id(RowIdx::new(0)), SourceId::Numeric(0));
        assert_eq!(c.output_id(RowIdx::new(1)), SourceId::Numeric(1));
    }

    #[test]
    fn lazy_massshift_expands_len_and_flags_targets() {
        let mut c = TargetColumns::with_capabilities(TargetCapabilities::test_lazy_decoys());
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
        assert!(c.is_target(FlatIdx::new(0))); // variant 0
        assert!(!c.is_target(FlatIdx::new(1))); // +decoy
        assert!(!c.is_target(FlatIdx::new(2))); // -decoy
        assert_eq!(c.decoy_group(RowIdx::new(0)), SourceId::Numeric(0));
    }

    #[test]
    fn passthrough_is_one_variant_per_row_honoring_is_decoy() {
        let mut caps = TargetCapabilities::default_diann();
        caps.decoys = DecoyStrategy::Passthrough;
        let mut c = TargetColumns::with_capabilities(caps);
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
        assert!(c.is_target(FlatIdx::new(0)));
        assert!(!c.is_target(FlatIdx::new(1))); // stored decoy row
    }

    #[test]
    fn csr_ranges_recover_per_target_fragments() {
        let mut c = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
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
        assert_eq!(c.frag_range(RowIdx::new(0)), 0..2);
        assert_eq!(c.frag_range(RowIdx::new(1)), 2..5);
        assert_eq!(c.frag_mzs[c.frag_range(RowIdx::new(1))][2], 600.0);
        assert_eq!(
            &c.seq_strip_blob[c.seq_strip_range(RowIdx::new(0))],
            "PEPTIDEK"
        );
        assert_eq!(
            &c.seq_strip_blob[c.seq_strip_range(RowIdx::new(1))],
            "SAMPLERK"
        );
    }

    #[test]
    fn string_labeled_collection_and_is_decoy() {
        use std::sync::Arc;
        let mut c: TargetColumns<Arc<str>> =
            TargetColumns::with_capabilities(TargetCapabilities::default_diann());
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
            c.frag_labels[c.frag_range(RowIdx::new(1))][0],
            Arc::<str>::from("frag_b")
        );
    }
}
