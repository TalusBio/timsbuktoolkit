use crate::KeyLike;
use crate::models::capabilities::{
    DecoyPolicy,
    DecoyStrategy,
    MASS_SHIFT_VARIANTS,
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
/// caller, confused with an id, or written to an output file -- a handle can
/// only be obtained from the arena and handed straight back to it.
///
/// Construction is `pub(super)`, i.e. this file, so the arena is the only thing
/// that mints one.
mod index {
    /// A stored row: `0..n_rows()`.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct RowIdx(u32);

    /// A scored slot after decoy expansion: `0..expanded_len()`. Distinct from
    /// [`RowIdx`] because `variants_per_row` slots map onto one row, so using
    /// one where the other is meant reads a real but wrong row -- in range, no
    /// panic, plausible data. The type is what stops that.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct FlatIdx(u32);

    /// Which competition group a row belongs to, as a handle rather than a
    /// value. Rows that compete share one, and that is all a consumer needs:
    /// grouping sorts by it and compares it, and never reads it. Opaque for the
    /// same reason as the indices above -- a handle that could be printed would
    /// end up in an output file as though it meant something.
    ///
    /// `Ord` because grouping works by sorting competitors adjacent; any total
    /// order does, since only equality carries meaning.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct GroupCode(u32);

    impl RowIdx {
        pub(super) fn new(row: u32) -> Self {
            Self(row)
        }

        pub(super) fn get(self) -> usize {
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

    /// `u32::MAX`, not 0, because 0 is a valid row in every non-empty arena: a
    /// placeholder that leaked would read row 0 and look right, where `u32::MAX`
    /// panics on access. A defaulted group likewise competes with the other
    /// placeholders and with nothing real, rather than folding into the arena's
    /// first group.
    ///
    /// These exist for exactly one caller: `Identity::sample_default()`, which
    /// `#[derive(ScoreBlock)]` deliberately leaves un-gated so other crates'
    /// tests can build a fixture. Nothing in production defaults a handle. There
    /// is no `Default for FlatIdx` because nothing ever wanted one.
    impl Default for RowIdx {
        fn default() -> Self {
            Self(u32::MAX)
        }
    }

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

/// One row's worth of input to [`TargetColumnsBuilder::push_row`].
///
/// Carries the row's `id`, so a row cannot be stored apart from the name its
/// file gave it. `Default` covers the optional half, so a caller names only
/// what it has:
///
/// ```ignore
/// geom.push_row(Row {
///     precursor_mz: 900.4,
///     charge: 2,
///     rt_seconds: 1.0,
///     mobility: 0.8,
///     frags: &frags,
///     seq_mod: "PEPTIDEK",
///     ..Default::default()
/// });
/// ```
#[derive(Debug, Clone)]
pub struct Row<'a, L: KeyLike> {
    pub precursor_mz: f64,
    pub charge: u8,
    pub rt_seconds: f32,
    pub mobility: f32,
    pub frags: &'a [(L, f64)],
    pub seq_strip: &'a str,
    pub seq_mod: &'a str,
    pub mods: &'a [(u8, u16)],
    /// A decoy the library shipped, as opposed to one the arena derives.
    pub is_decoy: bool,
    /// What the source file called this row, in the shape it used. `None` for a
    /// format that names nothing; [`TargetColumnsBuilder::seal`] mints ids only when
    /// *every* row arrived `None`.
    pub id: Option<OwnedSourceId>,
    /// Which competition group the file put this row in. Rows sharing one
    /// compete, so a shipped decoy has to name the same group as its target --
    /// which means both halves have to be named in one namespace, and a reader
    /// that has only the decoy's side of the pair has to name the target's group
    /// on the target's row too.
    ///
    /// `None` means "the file did not say", and a row the file did not place
    /// competes alone. [`TargetColumnsBuilder::seal`] stores nothing at all when that
    /// is true of every row, and when the groups it did get are all singletons
    /// -- both cases are what deriving the group from the row already does.
    ///
    /// One reader supplies this, and that is a fact about the formats rather
    /// than about this field: mzSpecLib's `related spectrum keys` is the only
    /// place a library states that one entry is the counterpart of another. A
    /// DIA-NN/Spectronaut/Skyline `transition_group_id` names the row it is on,
    /// not a partner, so those readers have nothing to pair by -- a decoy there
    /// is a row flagged as a decoy and nothing more.
    pub decoy_group: Option<OwnedSourceId>,
}

impl<L: KeyLike> Default for Row<'_, L> {
    fn default() -> Self {
        Self {
            precursor_mz: 0.0,
            charge: 0,
            rt_seconds: 0.0,
            mobility: 0.0,
            frags: &[],
            seq_strip: "",
            seq_mod: "",
            mods: &[],
            is_decoy: false,
            id: None,
            decoy_group: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ModDefinition {
    pub token: String, // verbatim, e.g. "[UNIMOD:4]"
    pub mono_delta: f64,
    pub cs_delta: (i16, i16),
}

#[derive(Debug, Clone)]
pub struct TargetColumns<L: KeyLike> {
    pub(crate) caps: TargetCapabilities,
    // per-target scalars, len = n_rows; addressed by `RowIdx`, never by a
    // caller-supplied id (those live in `source_ids`)
    pub(crate) precursor_mz: Vec<f64>,
    pub(crate) charge: Vec<u8>,
    pub(crate) rt_seconds: Vec<f32>,
    pub(crate) mobility: Vec<f32>,
    // per-row decoy flag (len = n_rows)
    pub(crate) is_decoy: Vec<bool>,
    /// Each row's id as pushed. Parallel to the row columns, so a shard merge
    /// cannot separate a row from its name. Drained by [`Self::seal`].
    pub(crate) pending_ids: Vec<Option<OwnedSourceId>>,
    /// Each row's declared competition group as pushed, `None` where the file
    /// declared none. Drained by [`Self::seal`], which is where the interned
    /// [`DecoyGroups`] is built or dropped.
    pub(crate) pending_groups: Vec<Option<OwnedSourceId>>,
    /// Built at [`Self::seal`], minted there for formats that carry no ids, so
    /// a sealed arena always has one per row.
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

/// The only mutable construction state for a [`TargetColumns`] arena.
///
/// Keeping this distinct from the sealed arena makes ids, competition groups,
/// and the decoy strategy mandatory before query methods become available.
/// Consuming [`Self::seal`] is the sole public transition into a scorable arena.
#[derive(Debug, Clone)]
pub struct TargetColumnsBuilder<L: KeyLike> {
    pub(crate) inner: TargetColumns<L>,
}

impl<L: KeyLike> TargetColumnsBuilder<L> {
    pub fn with_capabilities(caps: TargetCapabilities) -> Self {
        Self {
            inner: TargetColumns::empty(caps),
        }
    }

    pub fn push_row(&mut self, row: Row<'_, L>) {
        self.inner.push_row(row);
    }

    pub fn n_rows(&self) -> usize {
        self.inner.n_rows()
    }

    pub fn n_fragments(&self) -> usize {
        self.inner.n_fragments()
    }

    pub fn seal(self, decoys: DecoyPolicy) -> Result<TargetColumns<L>, SourceIdError> {
        self.inner.seal(decoys)
    }
}

impl<L: KeyLike> TargetColumns<L> {
    fn empty(caps: TargetCapabilities) -> Self {
        Self {
            caps,
            precursor_mz: Vec::new(),
            charge: Vec::new(),
            rt_seconds: Vec::new(),
            mobility: Vec::new(),
            is_decoy: Vec::new(),
            pending_ids: Vec::new(),
            pending_groups: Vec::new(),
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
    fn push_row(&mut self, row: Row<'_, L>) {
        let Row {
            precursor_mz,
            charge,
            rt_seconds,
            mobility,
            frags,
            seq_strip,
            seq_mod,
            mods,
            is_decoy,
            id,
            decoy_group,
        } = row;
        self.precursor_mz.push(precursor_mz);
        self.charge.push(charge);
        self.rt_seconds.push(rt_seconds);
        self.mobility.push(mobility);
        self.pending_ids.push(id);
        self.pending_groups.push(decoy_group);
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
    /// Number of *physical stored rows*, i.e. how many analytes are held in
    /// memory before decoy expansion. Under `MassShift` every row is a target;
    /// under `Stored` the count includes any decoy rows the file shipped.
    /// This is the base for `expanded_len` (the logical, iterator-length count).
    pub fn n_rows(&self) -> usize {
        self.charge.len()
    }

    /// The rows of this arena, in storage order.
    pub fn rows(&self) -> impl Iterator<Item = RowIdx> + use<L> {
        (0..self.n_rows() as u32).map(RowIdx::new)
    }

    /// Intern the competition groups the rows declared, or store nothing.
    ///
    /// Rows sharing a group compete, so unlike source ids these are expected to
    /// repeat -- hence a code per row against a deduplicated label set rather
    /// than a label per row.
    ///
    /// Nothing is stored in the two cases where the stored column would say no
    /// more than deriving the group from the row already does: no row declared
    /// one, or every group turned out to be a singleton. A reader can therefore
    /// name a group on every row it pushes without making a library that
    /// contains no pairs pay for a label column.
    ///
    /// A row that declared nothing while others did keeps its own id as its
    /// group, so it competes alone -- the group column has to cover every row
    /// once it exists, since a code indexes it directly.
    ///
    /// Refuses a group that declares two targets at one charge, which is a
    /// group the arena cannot compete: one result survives per group and
    /// charge, so the second target would be dropped. The check belongs here
    /// rather than in each reader because the arena is what holds `is_decoy`
    /// and `charge` beside the labels, so this is the one place that can answer
    /// it for every format at once. It rides the interning loop: one extra hash
    /// probe per target row, and a map bounded by the `(group, charge)` pairs
    /// that carry one.
    fn build_decoy_groups(&mut self) -> Result<(), SourceIdError> {
        let n = self.n_rows();
        debug_assert_eq!(self.pending_groups.len(), n, "a group slot per row");
        if self.pending_groups.iter().all(Option::is_none) {
            self.pending_groups.clear();
            self.pending_groups.shrink_to_fit();
            return Ok(());
        }

        let mut seen: std::collections::HashMap<OwnedSourceId, GroupCode> =
            std::collections::HashMap::new();
        let mut labels: Vec<OwnedSourceId> = Vec::new();
        let mut codes = Vec::with_capacity(n);
        // The target already claiming each `(group, charge)`, so a second one
        // can name the row it collides with.
        let mut claimed: std::collections::HashMap<(GroupCode, u8), usize> =
            std::collections::HashMap::new();
        for (row, group) in self.pending_groups.drain(..).enumerate() {
            let group = match group {
                Some(g) => g,
                // No id column yet -- `seal` builds it after this -- so mirror
                // what `output_id` will report for an unnamed row.
                None => match self.pending_ids.get(row) {
                    Some(Some(id)) => id.clone(),
                    _ => OwnedSourceId::Numeric(row as u64),
                },
            };
            let code = *seen.entry(group.clone()).or_insert_with(|| {
                labels.push(group);
                GroupCode::new((labels.len() - 1) as u32)
            });
            if !self.is_decoy[row] {
                let charge = self.charge[row];
                match claimed.entry((code, charge)) {
                    std::collections::hash_map::Entry::Occupied(first) => {
                        return Err(SourceIdError::GroupHasTwoTargets {
                            group: labels[code.get()].to_string(),
                            first_row: *first.get(),
                            row,
                            charge,
                        });
                    }
                    std::collections::hash_map::Entry::Vacant(slot) => {
                        slot.insert(row);
                    }
                }
            }
            codes.push(code);
        }
        self.pending_groups.shrink_to_fit();

        if labels.len() == n {
            // Every row alone in its own group, which is what `decoy_group_code`
            // derives for a groupless arena.
            return Ok(());
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

    /// The id a result for row `tgt` carries. Always a source id: [`TargetColumnsBuilder::seal`]
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

    /// Seal after build: resolve the decoy policy against the rows that
    /// actually arrived, build the id column, then release excess capacity on
    /// every arena.
    ///
    /// Consuming, and it takes the policy, so an arena is sealed exactly once
    /// and `caps.decoys` is written exactly there. That is what makes the
    /// invariant checkable rather than conventional: `MassShift` needs an
    /// all-targets arena (decoys are an on-the-fly ±CH2 index transform, never
    /// stored), and since the resolution happens after the rows are counted, a
    /// file that shipped decoys resolves to `Stored` instead of being silently
    /// re-decoyed. There is no later pass that could set it back.
    ///
    /// Fails when the pushed ids cannot make a column: two rows sharing an id
    /// (a caller keys results by it, so a repeat hides a row), or a library
    /// mixing numeric and text ids (storing both would coerce the numbers to
    /// strings, and `7` and `"7"` would become two ids in one column).
    ///
    /// Fails too when a declared group holds two targets at one charge, which is
    /// a group that cannot be competed: one result survives per group and
    /// charge, so loading it would silently discard a real identification. That
    /// is refused rather than warned about because the loss leaves no trace in
    /// the output -- the run reports a lower target count and nothing else --
    /// whereas a refusal names the file at the moment it is read.
    fn seal(mut self, decoys: DecoyPolicy) -> Result<Self, SourceIdError> {
        let stored_decoys = self.is_decoy.iter().filter(|&&is_decoy| is_decoy).count();
        if decoys == DecoyPolicy::Force && stored_decoys > 0 {
            return Err(SourceIdError::ForceWithStoredDecoys {
                count: stored_decoys,
            });
        }
        // Groups first: an undeclared row falls back to its own id, which is
        // still in `pending_ids` at this point.
        self.build_decoy_groups()?;
        self.build_source_ids()?;
        let ships_decoys = self.is_decoy.iter().any(|&d| d);
        self.caps.decoys = decoys.strategy(ships_decoys);
        // Nothing to build when no groups were declared: a row is then its own
        // group, which `decoy_group_code` derives. Only the case that actually
        // loses information is worth a word.
        if ships_decoys && self.decoy_groups.is_none() {
            tracing::warn!(
                "library ships its own decoys but declares no competition groups, so each \
                 stored decoy competes alone rather than against its target"
            );
        }
        if ships_decoys && decoys == DecoyPolicy::IfMissing {
            tracing::info!(
                "library ships {} of its own decoys, so none are derived",
                self.n_stored_decoys(),
            );
        }
        self.shrink();
        Ok(self)
    }

    pub fn capabilities(&self) -> &TargetCapabilities {
        &self.caps
    }

    fn build_source_ids(&mut self) -> Result<(), SourceIdError> {
        let n = self.n_rows();
        debug_assert_eq!(self.pending_ids.len(), n, "an id per row, or none");
        let named = self.pending_ids.iter().filter(|id| id.is_some()).count();
        self.source_ids = if named == n {
            SourceIds::owned(self.pending_ids.drain(..).map(Option::unwrap).collect(), n)?
        } else {
            if named > 0 {
                // A library names all its rows or none. One blank cell costs the
                // whole library its names, so say which row did it -- minting
                // overwrites every id, leaving nothing to trace it from.
                let row = self
                    .pending_ids
                    .iter()
                    .position(Option::is_none)
                    .expect("named < n, so some row is unnamed");
                tracing::warn!(
                    "{named} of {n} rows carry an id and {} do not, so the whole library falls \
                     back to minted ids: results will be keyed by row number instead of the \
                     names the file gave. First unnamed row is {row}.",
                    n - named,
                );
            } else {
                tracing::warn!(
                    "input carries no per-target id; minting {n} self-incremental ids (0..{n}). \
                     Result ids are ours, not the input file's."
                );
            }
            self.pending_ids.clear();
            SourceIds::minted(n)
        };
        self.pending_ids.shrink_to_fit();
        Ok(())
    }

    fn shrink(&mut self) {
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
/// `MassShift`), and expanded/flat indices fan each row out into its
/// target + decoy variants. Bounded on `DecoyShift` so `item_at` can hand out a
/// `QueryRef` (the flyweight that computes decoy geometry on the fly).
impl<L: KeyLike + DecoyShift> TargetColumns<L> {
    /// Variants each stored row expands into, from the decoy strategy alone:
    /// `MassShift` adds the ± pair (see [`MASS_SHIFT_VARIANTS`]); `Stored` is
    /// 1:1.
    ///
    /// Both values fit the `u8` that [`Self::split_flat`] casts a variant index
    /// to.
    pub fn variants_per_row(&self) -> usize {
        match self.caps.decoys {
            DecoyStrategy::MassShift { .. } => MASS_SHIFT_VARIANTS,
            DecoyStrategy::Stored => 1,
        }
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
    /// positions -- no integer names a slot outside this crate.
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
    /// route from a [`FlatIdx`] to a [`RowIdx`] -- which is what keeps the two
    /// index spaces from being used interchangeably.
    pub fn split_flat(&self, flat: FlatIdx) -> (RowIdx, u8) {
        let vpr = self.variants_per_row();
        let flat = flat.get();
        let row = u32::try_from(flat / vpr).expect("row index exceeds u32::MAX");
        (RowIdx::new(row), (flat % vpr) as u8)
    }

    /// The inverse of [`Self::split_flat`]: name one of a row's scored slots.
    ///
    /// Both directions live here so the encoding has one owner. This is also the
    /// only way to pair a row with a variant, and it checks the pairing -- a
    /// caller that could build the pair itself could name a variant the library
    /// does not expand into, which would read a real but wrong slot.
    pub fn flat_for(&self, row: RowIdx, variant: u8) -> FlatIdx {
        let vpr = self.variants_per_row();
        assert!(
            (variant as usize) < vpr,
            "variant {variant} does not exist in a library with {vpr} variants per row"
        );
        assert!(row.get() < self.n_rows(), "row is out of range");
        FlatIdx::new((row.get() * vpr + variant as usize) as u32)
    }

    /// Scored slot -> `(row, variant)` flyweight.
    pub fn item_at(&self, flat: FlatIdx) -> QueryRef<'_, L> {
        QueryRef::new(self, flat)
    }

    /// A slot is a target iff its stored row is a target AND it is the
    /// variant-0 (unshifted) slot.
    pub fn is_target(&self, flat: FlatIdx) -> bool {
        let (row, variant) = self.split_flat(flat);
        self.is_target_slot(row, variant)
    }

    /// [`Self::is_target`] for a caller that already holds the pair.
    ///
    /// The flyweights split their [`FlatIdx`] once and keep `(row, variant)`, so
    /// asking them to rebuild a flat index to get this answer would undo that.
    /// Both halves are needed: a stored decoy sits in the variant-0 slot exactly
    /// as a target does, and a mass-shift variant hangs off a row that is itself
    /// a target.
    pub fn is_target_slot(&self, row: RowIdx, variant: u8) -> bool {
        !self.is_decoy[row.get()] && variant == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IonAnnot;

    /// One row of a test arena in the terms the group checks read it: which
    /// group it declares, whether the library shipped it as a decoy, and its
    /// charge.
    ///
    /// Written as `Spec::target(..)` / `Spec::decoy(..)` so the decoy side of a
    /// fixture is legible at the call site rather than a positional `true`.
    #[derive(Clone, Copy)]
    struct Spec {
        group: Option<&'static str>,
        is_decoy: bool,
        charge: u8,
    }

    impl Spec {
        fn target(group: &'static str) -> Self {
            Self {
                group: Some(group),
                is_decoy: false,
                charge: 2,
            }
        }

        fn decoy(group: &'static str) -> Self {
            Self {
                is_decoy: true,
                ..Self::target(group)
            }
        }

        /// A row the file placed in no group, which is every format but
        /// mzSpecLib.
        fn undeclared() -> Self {
            Self {
                group: None,
                is_decoy: false,
                charge: 2,
            }
        }

        fn at_charge(self, charge: u8) -> Self {
            Self { charge, ..self }
        }
    }

    fn arena_of(specs: &[Spec]) -> TargetColumnsBuilder<IonAnnot> {
        let mut c = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        for spec in specs {
            c.push_row(Row {
                precursor_mz: 500.0,
                charge: spec.charge,
                rt_seconds: 1.0,
                mobility: 0.8,
                frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
                seq_strip: "PEP",
                seq_mod: "PEP",
                is_decoy: spec.is_decoy,
                decoy_group: spec.group.map(Into::into),
                ..Default::default()
            });
        }
        c
    }

    fn two_rows() -> TargetColumnsBuilder<IonAnnot> {
        arena_of(&[Spec::undeclared(), Spec::undeclared()])
    }

    /// Declared groups are interned, so a target and the decoy the library
    /// shipped for it share a code and the label is stored once rather than per
    /// row. The ordinary declared pair, which has to keep loading.
    #[test]
    fn rows_in_one_declared_group_share_a_code() {
        let c = arena_of(&[Spec::target("g1"), Spec::decoy("g1")])
            .seal(DecoyPolicy::Never)
            .expect("fixture ids are usable");

        let (a, b) = (c.rows().next().unwrap(), c.rows().nth(1).unwrap());
        assert_eq!(c.decoy_group_code(a), c.decoy_group_code(b));
        assert_eq!(c.decoy_group(a), SourceId::Text("g1"));
    }

    /// One survivor per group *and charge*, so a group naming a target at each
    /// of two charges discards nothing -- that is the ordinary multi-charge
    /// library, and refusing it would reject most of them.
    #[test]
    fn one_group_holds_a_target_at_each_charge() {
        let c = arena_of(&[Spec::target("g1"), Spec::target("g1").at_charge(3)])
            .seal(DecoyPolicy::Never)
            .expect("two charges of one precursor compete separately");

        let (a, b) = (c.rows().next().unwrap(), c.rows().nth(1).unwrap());
        assert_eq!(c.decoy_group_code(a), c.decoy_group_code(b));
    }

    /// A group with two targets at one charge is a group the arena cannot
    /// compete: one of the two would be dropped for the other, and the run would
    /// report a lower target count and nothing else. Sealing refuses it, naming
    /// the group and both rows, since the shared label is all a reader of a
    /// million-row library would otherwise have to go on.
    #[test]
    fn a_group_declaring_two_targets_at_one_charge_does_not_seal() {
        let err = arena_of(&[Spec::target("g1"), Spec::target("g1")])
            .seal(DecoyPolicy::Never)
            .expect_err("a group that cannot be competed is not a sealable arena");

        assert_eq!(
            err,
            SourceIdError::GroupHasTwoTargets {
                group: "g1".to_string(),
                first_row: 0,
                row: 1,
                charge: 2,
            }
        );
    }

    /// Declared groups that are all singletons say no more than deriving from
    /// the row does, so the column is dropped rather than stored. This is what
    /// lets a reader name a group on every row without charging a library that
    /// contains no pairs for a label per row.
    #[test]
    fn all_singleton_groups_are_not_stored() {
        let c = arena_of(&[Spec::target("g1"), Spec::target("g2")])
            .seal(DecoyPolicy::Never)
            .expect("fixture ids are usable");

        assert!(c.decoy_groups.is_none());
        let (a, b) = (c.rows().next().unwrap(), c.rows().nth(1).unwrap());
        assert_ne!(c.decoy_group_code(a), c.decoy_group_code(b));
    }

    /// A row the file left out of a group has to get a code anyway once any row
    /// declared one, since a code indexes the label column directly. It competes
    /// alone, under its own id.
    #[test]
    fn a_row_outside_every_declared_group_competes_alone() {
        let mut c = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        // `b` is the decoy of `a`: a group holds one target per charge, so the
        // pair is the only two-row group there is at a single charge.
        for (id, group, is_decoy) in [
            ("a", Some("pair"), false),
            ("b", Some("pair"), true),
            ("c", None, false),
        ] {
            c.push_row(Row {
                precursor_mz: 500.0,
                charge: 2,
                rt_seconds: 1.0,
                mobility: 0.8,
                frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
                seq_strip: "PEP",
                seq_mod: "PEP",
                is_decoy,
                id: Some(id.into()),
                decoy_group: group.map(Into::into),
                ..Default::default()
            });
        }
        let c = c.seal(DecoyPolicy::Never).expect("fixture ids are usable");

        let rows: Vec<_> = c.rows().collect();
        assert_eq!(c.decoy_group_code(rows[0]), c.decoy_group_code(rows[1]));
        assert_ne!(c.decoy_group_code(rows[0]), c.decoy_group_code(rows[2]));
        assert_eq!(c.decoy_group(rows[2]), SourceId::Text("c"));
    }

    /// Without declared groups a row is its own group, derived rather than
    /// stored -- so nothing is allocated and the label is the row's own id.
    #[test]
    fn undeclared_groups_are_derived_not_stored() {
        let c = two_rows();
        let c = c.seal(DecoyPolicy::Never).expect("fixture ids are usable");

        let (a, b) = (c.rows().next().unwrap(), c.rows().nth(1).unwrap());
        assert!(c.decoy_groups.is_none(), "nothing stored for minted groups");
        assert_ne!(
            c.decoy_group_code(a),
            c.decoy_group_code(b),
            "each row competes alone"
        );
        assert_eq!(c.decoy_group(a), SourceId::Numeric(0));
        assert_eq!(c.decoy_group(b), SourceId::Numeric(1));
    }

    /// `output_id` expects a source id on every row; this is what makes that
    /// safe for the formats that carry none.
    #[test]
    fn seal_mints_ids_when_the_input_carried_none() {
        let c = two_rows();
        let c = c.seal(DecoyPolicy::Never).expect("fixture ids are usable");
        assert_eq!(c.output_id(RowIdx::new(0)), SourceId::Numeric(0));
        assert_eq!(c.output_id(RowIdx::new(1)), SourceId::Numeric(1));
    }

    #[test]
    fn lazy_massshift_expands_len_and_flags_targets() {
        let mut c = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        c.push_row(Row {
            precursor_mz: 500.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 0.8,
            frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            ..Default::default()
        });
        let c = c
            .seal(DecoyPolicy::IfMissing)
            .expect("fixture ids are usable");
        assert_eq!(c.variants_per_row(), 3);
        assert_eq!(c.expanded_len(), 3);
        assert!(c.is_target(FlatIdx::new(0))); // variant 0
        assert!(!c.is_target(FlatIdx::new(1))); // +decoy
        assert!(!c.is_target(FlatIdx::new(2))); // -decoy
        assert_eq!(c.decoy_group(RowIdx::new(0)), SourceId::Numeric(0));
    }

    #[test]
    fn stored_is_one_variant_per_row_honoring_is_decoy() {
        let mut c = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        c.push_row(Row {
            precursor_mz: 500.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 0.8,
            frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            ..Default::default()
        });
        c.push_row(Row {
            precursor_mz: 510.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 0.8,
            frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            is_decoy: true,
            ..Default::default()
        });
        // A shipped decoy is what makes `IfMissing` resolve to `Stored`.
        let c = c
            .seal(DecoyPolicy::IfMissing)
            .expect("fixture ids are usable");
        assert_eq!(c.variants_per_row(), 1);
        assert_eq!(c.expanded_len(), 2);
        assert!(c.is_target(FlatIdx::new(0)));
        assert!(!c.is_target(FlatIdx::new(1))); // stored decoy row
    }

    #[test]
    fn csr_ranges_recover_per_target_fragments() {
        let mut c = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        // target 0: 2 frags; target 1: 3 frags
        c.push_row(Row {
            precursor_mz: 500.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 0.8,
            frags: &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y4").unwrap(), 400.0),
            ],
            seq_strip: "PEPTIDEK",
            seq_mod: "PEPTIDEK",
            ..Default::default()
        });
        c.push_row(Row {
            precursor_mz: 600.0,
            charge: 3,
            rt_seconds: 2.0,
            mobility: 0.9,
            frags: &[
                (IonAnnot::try_from("y2").unwrap(), 200.0),
                (IonAnnot::try_from("y5").unwrap(), 500.0),
                (IonAnnot::try_from("y6").unwrap(), 600.0),
            ],
            seq_strip: "SAMPLERK",
            seq_mod: "SAMPLERK",
            ..Default::default()
        });
        let c = c.seal(DecoyPolicy::Never).expect("fixture ids are usable");
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
        let mut c: TargetColumnsBuilder<Arc<str>> =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        c.push_row(Row {
            precursor_mz: 500.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 0.8,
            frags: &[(Arc::<str>::from("frag_a"), 300.0)],
            seq_strip: "PEP",
            seq_mod: "PEP",
            ..Default::default()
        });
        c.push_row(Row {
            precursor_mz: 600.0,
            charge: 2,
            rt_seconds: 1.0,
            mobility: 0.8,
            frags: &[(Arc::<str>::from("frag_b"), 400.0)],
            seq_strip: "TIDE",
            seq_mod: "TIDE",
            is_decoy: true,
            ..Default::default()
        });
        let c = c.seal(DecoyPolicy::Never).expect("fixture ids are usable");
        assert_eq!(c.n_rows(), 2);
        assert_eq!(c.is_decoy, vec![false, true]);
        assert_eq!(
            c.frag_labels[c.frag_range(RowIdx::new(1))][0],
            Arc::<str>::from("frag_b")
        );
    }
}
