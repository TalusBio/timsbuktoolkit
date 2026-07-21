//! Reader for the DIA-NN `.speclib` binary spectral library format.
//!
//! The format was reverse-engineered from `empirical_library.tsv.speclib` and
//! cross-checked against ProteoWizard/BiblioSpec `DiaNNSpecLibReader.cpp`
//! (`Library::read`). See the `reverse_eng_speclib` spec for the byte layout.
//!
//! The file is a sequence of variable-length sections written back-to-back
//! (little-endian, no padding, no section table); it must be parsed strictly in
//! order, so there is no random access across sections.
//!
//! The reader is three layers:
//!   1. Decode ([`Cursor`], [`Fragment`], [`Peptide`]): typed, zero-copy views
//!      over byte ranges of the in-memory file. No domain knowledge.
//!   2. Emit ([`SpecLib`] + [`EntryIter`]): a pull parser that walks the
//!      variable-length envelope and yields one [`EntryView`] per precursor.
//!   3. Map ([`map_entry`]): [`EntryView`] -> rows pushed directly into a
//!      columnar `QueryCollection<IonAnnot>` (plus a parallel reference-intensity
//!      sidecar). This is where the ion-type table, dedup, and drop stats live.
//!
//! The `Fragment`/`Peptide` split falls out of the layout: a `Product` is a
//! fixed 12-byte record (cheap fixed-offset view), and a `Peptide`'s only
//! variable part (its fragment block) sits at the end after a version-gated
//! fixed header (so scalar fields are O(1) offset reads). The `Entry` envelope
//! has interior variable-length data (an optional embedded decoy peptide, then a
//! length-prefixed name) and so is consumed sequentially by the iterator rather
//! than exposed as a random-access view.

use super::library_file::{
    LibraryArena,
    LibraryReadingError,
};
use crate::ion::IonAnnot;
use crate::models::{
    LibCapabilities,
    QueryCollection,
};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use tracing::{
    info,
    warn,
};

/// Newest (most-negative) `.speclib` version this reader understands. Files with
/// a version below this were written by a newer DIA-NN than was reverse-engineered.
const LATEST_SUPPORTED_VERSION: i32 = -3;

/// Bytes of one `Product`/fragment record.
const PRODUCT_SIZE: usize = 12;

/// Smallest plausible `Entry` on disk (peptide header + zero-fragment count +
/// `dc` + entry_flags + proteotypic + pid_index + empty-name count). Used only
/// to bound pre-allocation against a corrupt count prefix — never a parse gate.
const MIN_ENTRY_BYTES: usize = 40;

// --- Layer 1: decode ------------------------------------------------------

/// Read a little-endian `i32` from `buf` at `offset`. The caller guarantees the
/// range is in bounds (view types only expose offsets validated at construction).
fn le_i32(buf: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap())
}

fn le_f32(buf: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap())
}

/// Little-endian cursor over the whole file in memory. Each accessor advances the
/// position and bounds-checks against the buffer end (a short read is a truncated
/// library, reported as [`LibraryReadingError::SpeclibParse`]).
struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    /// Borrow the next `n` bytes and advance.
    fn take(&mut self, n: usize) -> Result<&'a [u8], LibraryReadingError> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.data.len());
        match end {
            Some(end) => {
                let slice = &self.data[self.pos..end];
                self.pos = end;
                Ok(slice)
            }
            None => Err(LibraryReadingError::SpeclibParse(format!(
                "unexpected end of file: wanted {n} bytes at offset {} of {}",
                self.pos,
                self.data.len()
            ))),
        }
    }

    /// Borrow everything consumed since `start` (an earlier [`Cursor::mark`]).
    fn slice_since(&self, start: usize) -> &'a [u8] {
        &self.data[start..self.pos]
    }

    fn mark(&self) -> usize {
        self.pos
    }

    fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    fn read_i32(&mut self) -> Result<i32, LibraryReadingError> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn read_f64(&mut self) -> Result<f64, LibraryReadingError> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }

    /// `i32 count` that must be non-negative (a length/count prefix).
    fn read_count(&mut self) -> Result<usize, LibraryReadingError> {
        let n = self.read_i32()?;
        if n < 0 {
            return Err(LibraryReadingError::SpeclibParse(format!(
                "negative count/length prefix: {n}"
            )));
        }
        Ok(n as usize)
    }

    /// `i32 length` + `length` raw bytes decoded as latin-1 (1 byte per char), no
    /// terminator. Latin-1 (not UTF-8) because high bytes appear in Windows paths
    /// and accented protein names.
    fn read_str(&mut self) -> Result<String, LibraryReadingError> {
        let n = self.read_count()?;
        Ok(self.take(n)?.iter().map(|&b| b as char).collect())
    }

    /// Skip `n` bytes.
    fn skip_bytes(&mut self, n: usize) -> Result<(), LibraryReadingError> {
        self.take(n).map(|_| ())
    }

    /// `vec<i32>`: `i32 count` + `count` contiguous i32 — skip.
    fn skip_vec_i32(&mut self) -> Result<(), LibraryReadingError> {
        let n = self.read_count()?;
        self.skip_bytes(n * 4)
    }

    fn at_eof(&self) -> bool {
        self.pos >= self.data.len()
    }
}

/// Zero-copy view over one 12-byte `Product`/fragment record.
struct Fragment<'a>(&'a [u8]);

impl Fragment<'_> {
    fn mz(&self) -> f32 {
        le_f32(self.0, 0)
    }

    fn height(&self) -> f32 {
        le_f32(self.0, 4)
    }

    fn charge(&self) -> u8 {
        self.0[8]
    }

    fn typ(&self) -> u8 {
        self.0[9]
    }

    fn index(&self) -> u8 {
        self.0[10]
    }

    fn loss(&self) -> u8 {
        self.0[11]
    }
}

/// Zero-copy view over one `Peptide` record (target or embedded decoy). `buf` is
/// exactly the peptide's byte span; scalar fields are fixed-offset reads and the
/// trailing fragment block is a lazy [`Fragment`] iterator.
///
/// The header length is version-gated: `v <= -2` inserts `lib_qvalue`, `i_im`,
/// `s_im` after `s_rt`. Keeping the version here means the gating lives in one
/// place (offset arithmetic) rather than scattered through a read routine.
struct Peptide<'a> {
    buf: &'a [u8],
    version: i32,
}

impl<'a> Peptide<'a> {
    /// Bytes of the fixed header before the `nfrag` count.
    fn header_len(version: i32) -> usize {
        // index, charge, length, mz, i_rt, s_rt (6 * 4 = 24); v <= -2 adds
        // lib_qvalue, i_im, s_im (3 * 4 = 12).
        if version <= -2 { 36 } else { 24 }
    }

    fn charge(&self) -> i32 {
        le_i32(self.buf, 4)
    }

    fn length(&self) -> i32 {
        le_i32(self.buf, 8)
    }

    fn mz(&self) -> f32 {
        le_f32(self.buf, 12)
    }

    fn i_rt(&self) -> f32 {
        le_f32(self.buf, 16)
    }

    /// Ion mobility (1/K0). Absent (returns 0.0) before version -2.
    fn i_im(&self) -> f32 {
        if self.version <= -2 {
            le_f32(self.buf, 28)
        } else {
            0.0
        }
    }

    fn fragments(&self) -> impl Iterator<Item = Fragment<'a>> {
        // The fragment block starts after the fixed header and the `nfrag`
        // count; its length was validated when the view was constructed.
        let start = Peptide::header_len(self.version) + 4;
        self.buf[start..].chunks_exact(PRODUCT_SIZE).map(Fragment)
    }
}

/// `Peptide::read` — consume one peptide from `c` and return a view over its
/// span. Works for both the target and (when `dc != 0`) the embedded decoy.
fn read_peptide<'a>(c: &mut Cursor<'a>, version: i32) -> Result<Peptide<'a>, LibraryReadingError> {
    let start = c.mark();
    c.skip_bytes(Peptide::header_len(version))?;
    let nfrag = c.read_count()?;
    c.skip_bytes(nfrag * PRODUCT_SIZE)?;
    Ok(Peptide {
        buf: c.slice_since(start),
        version,
    })
}

/// A decoded `Entry` envelope: the target peptide view plus the Entry-level
/// fields needed for mapping. `had_file_decoy` records that a `dc != 0` embedded
/// decoy was read-to-sync and discarded (counted by the mapping layer).
pub(crate) struct EntryView<'a> {
    peptide: Peptide<'a>,
    name: String,
    protein_id: String,
    had_file_decoy: bool,
}

/// `Entry::read` — a target `Peptide`, an optional embedded decoy, then the
/// Entry-level fields.
fn read_entry<'a>(
    c: &mut Cursor<'a>,
    version: i32,
    pg_ids: &[String],
) -> Result<EntryView<'a>, LibraryReadingError> {
    let peptide = read_peptide(c, version)?;

    let dc = c.read_i32()?;
    let had_file_decoy = dc != 0;
    if had_file_decoy {
        // A whole second Peptide is embedded between `dc` and `entry_flags`.
        // Read it (mandatory to stay byte-synced) then discard: a lone file
        // decoy would land in its own singleton FDR group and always "win",
        // silently breaking FDR. Mass-shift decoys are generated downstream.
        let _decoy = read_peptide(c, version)?;
    }

    let _entry_flags = c.read_i32()?;
    let _proteotypic = c.read_i32()?;
    let pid_index = c.read_i32()?;
    let name = c.read_str()?;
    if version <= -3 {
        // pg_qvalue, ptm_qvalue, site_conf (3 * f32).
        c.skip_bytes(12)?;
    }

    let protein_id = if pid_index >= 0 {
        pg_ids.get(pid_index as usize).cloned().unwrap_or_default()
    } else {
        String::new()
    };

    Ok(EntryView {
        peptide,
        name,
        protein_id,
        had_file_decoy,
    })
}

/// Advance `c` past exactly one entry without allocating. Mirrors [`read_entry`]'s
/// cursor walk (skipping the name rather than decoding it) so the serial offset
/// scan stays cheap before the parallel map re-parses each entry from its offset.
fn skip_entry(c: &mut Cursor, version: i32) -> Result<(), LibraryReadingError> {
    let _ = read_peptide(c, version)?; // target peptide
    if c.read_i32()? != 0 {
        let _ = read_peptide(c, version)?; // embedded decoy, kept only to stay synced
    }
    let _entry_flags = c.read_i32()?;
    let _proteotypic = c.read_i32()?;
    let _pid_index = c.read_i32()?;
    let name_len = c.read_count()?;
    c.skip_bytes(name_len)?; // name bytes — skipped, decoded later in parallel
    if version <= -3 {
        c.skip_bytes(12)?; // pg_qvalue, ptm_qvalue, site_conf
    }
    Ok(())
}

// --- Layer 2: emit --------------------------------------------------------

/// Backing storage for the file bytes: an owned buffer (from a `Read` stream) or
/// a memory map (from a path). Decode views borrow `&[u8]` from either, so the
/// rest of the reader is agnostic to which.
enum Backing {
    Owned(Vec<u8>),
    Mapped(memmap2::Mmap),
}

impl Backing {
    fn bytes(&self) -> &[u8] {
        match self {
            Backing::Owned(v) => v,
            Backing::Mapped(m) => m,
        }
    }
}

/// A parsed `.speclib` held in memory. The header and PG `ids` are parsed
/// eagerly on open; entries are parsed and mapped by [`SpecLib::parse_parallel`]
/// (serial offset scan, then parallel per-entry map).
pub(crate) struct SpecLib {
    data: Backing,
    version: i32,
    pg_ids: Vec<String>,
    /// Byte offset of the entries count prefix (section 7).
    entries_start: usize,
    n_entries: usize,
}

impl SpecLib {
    /// Read a whole stream into an owned buffer, then parse the header. Kept for
    /// streams and tests; the path-based file reader prefers [`SpecLib::open_mmap`].
    pub(crate) fn open<R: Read>(mut reader: R) -> Result<Self, LibraryReadingError> {
        let mut data = Vec::new();
        reader
            .read_to_end(&mut data)
            .map_err(LibraryReadingError::IoError)?;
        Self::from_backing(Backing::Owned(data))
    }

    /// Memory-map a `.speclib` file, then parse the header. Avoids a ~file-sized
    /// resident buffer — pages fault in on demand and are OS-reclaimable.
    pub(crate) fn open_mmap(path: &Path) -> Result<Self, LibraryReadingError> {
        let file = File::open(path).map_err(LibraryReadingError::IoError)?;
        // SAFETY: read-only ingest. The library file is not written during load;
        // an mmap of a file mutated/truncated concurrently would be UB.
        let mmap = unsafe { memmap2::Mmap::map(&file) }.map_err(LibraryReadingError::IoError)?;
        Self::from_backing(Backing::Mapped(mmap))
    }

    /// Parse everything up to (not including) the entries section from `data`.
    fn from_backing(data: Backing) -> Result<Self, LibraryReadingError> {
        let mut c = Cursor::new(data.bytes());

        // --- Section 0: header ---
        // A leading i32 >= 0 means v0 with no version word (the value is already
        // gen_decoys); a negative leading word is the version.
        let first = c.read_i32()?;
        let version = if first >= 0 {
            0
        } else {
            let _gen_decoys = c.read_i32()?;
            first
        };
        if version < LATEST_SUPPORTED_VERSION {
            return Err(LibraryReadingError::UnsupportedSpeclibVersion(version));
        }
        let _gen_charges = c.read_i32()?;
        let _infer_proteotypicity = c.read_i32()?;
        let _name = c.read_str()?;
        let _fasta_names = c.read_str()?;

        // --- Section 1: proteins (Isoform[]) — walk & discard ---
        let n_isoforms = c.read_count()?;
        for _ in 0..n_isoforms {
            skip_isoform(&mut c)?;
        }

        // --- Section 2: protein_ids (PG[]) — keep only `ids` per group ---
        let n_pg = c.read_count()?;
        let mut pg_ids = Vec::with_capacity(bounded_capacity(n_pg, c.remaining(), MIN_ENTRY_BYTES));
        for _ in 0..n_pg {
            pg_ids.push(read_pg_keep_ids(&mut c)?);
        }

        // --- Sections 3-5: precursor / name / gene string dictionaries — discard ---
        for _ in 0..3 {
            skip_strings(&mut c)?;
        }

        // --- Section 6: iRT range (2x f64) — discard (per-entry iRT present) ---
        let _irt_min = c.read_f64()?;
        let _irt_max = c.read_f64()?;

        // --- Section 7: entries — recorded, parsed lazily ---
        let n_entries = c.read_count()?;
        let entries_start = c.mark();

        Ok(SpecLib {
            data,
            version,
            pg_ids,
            entries_start,
            n_entries,
        })
    }

    /// Bound pre-allocation so a corrupt entry count can't reserve gigabytes.
    fn capacity_hint(&self) -> usize {
        bounded_capacity(
            self.n_entries,
            self.data.bytes().len().saturating_sub(self.entries_start),
            MIN_ENTRY_BYTES,
        )
    }

    /// Parse and map every entry, doing the per-entry work in parallel.
    ///
    /// Phase A (serial, allocation-free): walk the variable-length envelope
    /// recording each entry's start offset. Boundary finding is inherently
    /// sequential — an entry's position depends on all prior entries' lengths.
    /// This is the ~1 core the format forces; it does no allocation or mapping.
    ///
    /// Phase B (parallel): each entry is re-parsed and mapped independently from
    /// its offset. All per-entry work (`read_entry`, `map_entry`) is pure and
    /// only shares immutable borrows of the buffer + PG ids, so it fans across
    /// cores with no synchronization. Order and drop stats are preserved via an
    /// ordered fold/reduce.
    ///
    /// Returns the mapped entries in file order, merged drop stats, and whether
    /// parsing landed exactly on EOF (a `false` signals a structural desync).
    pub(crate) fn parse_parallel(&self) -> Result<ParsedSpeclib, LibraryReadingError> {
        use rayon::prelude::*;

        // Phase A: serial offset scan, then the trailing section 8 / EOF check.
        let mut c = Cursor::new(self.data.bytes());
        c.pos = self.entries_start;
        let mut offsets = Vec::with_capacity(self.capacity_hint());
        for _ in 0..self.n_entries {
            offsets.push(c.mark());
            skip_entry(&mut c, self.version)?;
        }
        if self.version <= -1 && !c.at_eof() {
            c.skip_vec_i32()?; // section 8: elution_groups (discard)
        }
        let at_eof = c.at_eof();

        // Phase B: parallel parse + map straight into per-worker arena shards,
        // concatenated back into file order. `try_fold`/`try_reduce`
        // short-circuit on the first error and, on an indexed iterator, preserve
        // file order. Each worker owns a `(QueryCollection, frag_intens)` shard
        // it pushes rows into; the reduce step concatenates shards with
        // `append_arena` (which rebases every CSR offset). The 4th accumulator
        // slot is a per-worker dedup scratch buffer reused across every entry
        // that worker maps (see `map_entry`), so the hot loop allocates only the
        // persistent output, not a fresh dedup Vec per entry.
        let (geom, frag_intens, stats, _scratch) = offsets
            .par_iter()
            .try_fold(
                || {
                    (
                        QueryCollection::with_capabilities(
                            LibCapabilities::default_diann_no_decoys(),
                        ),
                        Vec::<f32>::new(),
                        SpeclibDecodeStats::default(),
                        Vec::new(),
                    )
                },
                |(mut geom, mut frag_intens, mut stats, mut scratch),
                 &off|
                 -> Result<_, LibraryReadingError> {
                    let mut cur = Cursor::new(self.data.bytes());
                    cur.pos = off;
                    let view = read_entry(&mut cur, self.version, &self.pg_ids)?;
                    map_entry(view, &mut geom, &mut frag_intens, &mut stats, &mut scratch)?;
                    Ok((geom, frag_intens, stats, scratch))
                },
            )
            .try_reduce(
                || {
                    (
                        QueryCollection::with_capabilities(
                            LibCapabilities::default_diann_no_decoys(),
                        ),
                        Vec::<f32>::new(),
                        SpeclibDecodeStats::default(),
                        Vec::new(),
                    )
                },
                |(mut a_geom, mut a_int, sa, scratch), (b_geom, b_int, sb, _)| {
                    append_arena(&mut a_geom, b_geom);
                    a_int.extend(b_int);
                    let mut merged = sa;
                    merged.merge(&sb);
                    Ok((a_geom, a_int, merged, scratch))
                },
            )?;

        Ok((geom, frag_intens, stats, at_eof))
    }
}

/// `Isoform::read` — read & discard.
fn skip_isoform(c: &mut Cursor) -> Result<(), LibraryReadingError> {
    let _sp = c.read_i32()?;
    let size = c.read_count()?;
    let _id = c.read_str()?;
    let _name = c.read_str()?;
    let _gene = c.read_str()?;
    let _name_index = c.read_i32()?;
    let _gene_index = c.read_i32()?;
    c.skip_bytes(size * 4) // prec[size]
}

/// `PG::read` — keep only the `;`-joined `ids` string.
fn read_pg_keep_ids(c: &mut Cursor) -> Result<String, LibraryReadingError> {
    let size_p = c.read_count()?;
    let ids = c.read_str()?;
    let _names = c.read_str()?;
    let _genes = c.read_str()?;
    c.skip_vec_i32()?; // name_indices
    c.skip_vec_i32()?; // gene_indices
    c.skip_vec_i32()?; // precursors
    c.skip_bytes(size_p * 4)?; // proteins[size_p]
    Ok(ids)
}

/// `read_strings` — read & discard.
fn skip_strings(c: &mut Cursor) -> Result<(), LibraryReadingError> {
    let n = c.read_count()?;
    for _ in 0..n {
        let _ = c.read_str()?;
    }
    Ok(())
}

/// Cap a count-derived capacity at what the remaining bytes could actually hold,
/// so a corrupt/garbage count prefix can't trigger a multi-gigabyte reservation
/// before the per-record bounds checks ever run.
fn bounded_capacity(count: usize, remaining_bytes: usize, min_record_bytes: usize) -> usize {
    count.min(remaining_bytes / min_record_bytes + 1)
}

// --- Layer 3: map ---------------------------------------------------------

/// Counts of fragments/records dropped during a decode — surfaced via `warn!`
/// and returned for callers to assert on.
#[derive(Debug, Default, Clone, Copy)]
pub struct SpeclibDecodeStats {
    /// Fragments flagged `ExcludeFromAssay` (`type & 0x80`). Kept (see
    /// [`map_entry`]); counted for reporting only.
    pub exclude_flagged: usize,
    /// Fragments with a neutral loss (`loss != 0`) — `IonAnnot` cannot hold loss.
    pub loss_dropped: usize,
    /// Fragments whose ion-type code or recovered series was unusable.
    pub unknown_ion_dropped: usize,
    /// Fragments collapsed into an earlier one sharing the same `IonAnnot`.
    pub dedup_dropped: usize,
    /// Embedded (`dc != 0`) file decoys read-to-sync then discarded.
    pub decoys_dropped: usize,
}

impl SpeclibDecodeStats {
    /// Accumulate another partial count set (per-worker merge in the parallel map).
    fn merge(&mut self, o: &SpeclibDecodeStats) {
        self.exclude_flagged += o.exclude_flagged;
        self.loss_dropped += o.loss_dropped;
        self.unknown_ion_dropped += o.unknown_ion_dropped;
        self.dedup_dropped += o.dedup_dropped;
        self.decoys_dropped += o.decoys_dropped;
    }

    fn any_dropped(&self) -> bool {
        self.loss_dropped > 0
            || self.unknown_ion_dropped > 0
            || self.dedup_dropped > 0
            || self.decoys_dropped > 0
    }
}

/// Strip DIA-NN mod annotations — anything inside `(...)` or `[...]` — leaving
/// the bare residue string.
fn strip_mods(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth: i32 = 0;
    for c in s.chars() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth = (depth - 1).max(0),
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

/// Residue count of a mod-stripped sequence.
fn residue_count(stripped: &str) -> usize {
    stripped.chars().filter(|c| c.is_ascii_alphabetic()).count()
}

/// Concatenate `src` onto `dst` in place, rebasing every CSR offset by `dst`'s
/// current arena lengths and each structured-mod registry index by `dst`'s mod
/// registry length. The empty arena is the identity, so this is the associative
/// reduce operator that merges the per-worker shards back into file order.
///
/// `dst.caps` is preserved (all shards share `default_diann_no_decoys`), so merging an
/// empty identity in either position is a no-op on capabilities.
fn append_arena(dst: &mut QueryCollection<IonAnnot>, mut src: QueryCollection<IonAnnot>) {
    // Bases captured BEFORE the backing arenas are appended.
    let frag_base = dst.frag_labels.len();
    let strip_base = dst.seq_strip_blob.len();
    let mod_seq_base = dst.seq_mod_blob.len();
    let mods_base = dst.mods.len();
    let reg_base = dst.mod_registry.len();

    dst.precursor_mz.append(&mut src.precursor_mz);
    dst.charge.append(&mut src.charge);
    dst.rt_seconds.append(&mut src.rt_seconds);
    dst.mobility.append(&mut src.mobility);
    dst.is_decoy.append(&mut src.is_decoy);
    dst.frag_labels.append(&mut src.frag_labels);
    dst.frag_mzs.append(&mut src.frag_mzs);
    dst.seq_strip_blob.push_str(&src.seq_strip_blob);
    dst.seq_mod_blob.push_str(&src.seq_mod_blob);

    for &(pos, reg) in &src.mods {
        let rebased =
            u16::try_from(reg as usize + reg_base).expect("mod registry index exceeds u16 range");
        dst.mods.push((pos, rebased));
    }
    dst.mod_registry.append(&mut src.mod_registry);

    // CSR offset arrays carry a leading 0; skip it and rebase the remainder onto
    // the running arena length. `try_from` mirrors the checked pushes in
    // `QueryCollection` — an overflow fails loud rather than wrapping an offset.
    dst.frag_off.extend(src.frag_off[1..].iter().map(|&o| {
        u32::try_from(o as usize + frag_base).expect("fragment arena exceeds u32 offset range")
    }));
    dst.seq_strip_off
        .extend(src.seq_strip_off[1..].iter().map(|&o| {
            u32::try_from(o as usize + strip_base)
                .expect("stripped-seq blob exceeds u32 offset range")
        }));
    dst.seq_mod_off
        .extend(src.seq_mod_off[1..].iter().map(|&o| {
            u32::try_from(o as usize + mod_seq_base)
                .expect("modified-seq blob exceeds u32 offset range")
        }));
    dst.mod_off.extend(src.mod_off[1..].iter().map(|&o| {
        u32::try_from(o as usize + mods_base).expect("mods arena exceeds u32 offset range")
    }));
}

/// Map a decoded target `EntryView` directly into the columnar arena, folding
/// drop counters into `stats` and appending each kept fragment's reference
/// intensity to `frag_intens` IN THE SAME ORDER as the pushed fragment labels
/// (so the sidecar stays parallel to `geom.frag_labels`).
///
/// `scratch` is a caller-owned dedup buffer reused across entries (one per rayon
/// worker) so the per-entry dedup pass allocates nothing; it is cleared on entry.
fn map_entry(
    entry: EntryView,
    geom: &mut QueryCollection<IonAnnot>,
    frag_intens: &mut Vec<f32>,
    stats: &mut SpeclibDecodeStats,
    scratch: &mut Vec<(IonAnnot, f64, f32)>,
) -> Result<(), LibraryReadingError> {
    let pep = &entry.peptide;
    let name = entry.name;
    if entry.had_file_decoy {
        stats.decoys_dropped += 1;
    }

    // `name` (== transition_group_id) is `modified_peptide + charge`. Strip the
    // exact known-charge decimal suffix — not "trailing digits" — so a
    // C-terminal mod ending in a digit is left intact.
    let charge = pep.charge();
    if charge < 1 || charge > u8::MAX as i32 {
        return Err(LibraryReadingError::SpeclibParse(format!(
            "entry {name:?}: precursor charge {charge} out of range 1..=255 (corrupt record?)"
        )));
    }
    let charge_suffix = charge.to_string();
    let modified_peptide = match name.strip_suffix(&charge_suffix) {
        Some(prefix) => prefix.to_string(),
        None => {
            warn!(
                "speclib entry name {:?} does not end with charge suffix {:?}; keeping whole name",
                name, charge_suffix
            );
            name.clone()
        }
    };
    let stripped_peptide = strip_mods(&modified_peptide);

    // `Peptide.length` may be 0 in some libraries; recover it from the sequence
    // when needed (b/a/c series don't need it, y/x/z series do).
    let resolved_len = if pep.length() > 0 {
        pep.length() as usize
    } else {
        residue_count(&stripped_peptide)
    };

    // (IonAnnot, fragment mz as f64, height). Dedup by IonAnnot keeping max
    // height so a duplicate label can't fail the whole load via
    // `ExpectedIntensities::try_from_pairs`. Reuses the caller's scratch buffer.
    scratch.clear();
    let kept = scratch;

    for f in pep.fragments() {
        // `type & 0x80` => ExcludeFromAssay: DIA-NN's non-quantifier transitions.
        // They carry valid m/z + intensity and are kept — dropping them leaves
        // ~3 fragments/precursor on real libraries, which starves cross-library
        // calibration matching (needs >=5 shared fragments). Counted only.
        if f.typ() & 0x80 != 0 {
            stats.exclude_flagged += 1;
        }
        // IonAnnot cannot represent neutral loss; drop lossy fragments.
        if f.loss() != 0 {
            stats.loss_dropped += 1;
            continue;
        }

        let type_char = match f.typ() & 0x7F {
            1 => 'b',
            2 => 'y',
            3 => 'a',
            4 => 'x',
            5 => 'c',
            6 => 'z',
            other => {
                warn!(
                    "speclib entry {:?}: unknown ion-type code {}; dropping fragment",
                    name, other
                );
                stats.unknown_ion_dropped += 1;
                continue;
            }
        };

        // N-terminal (b/a/c): series = index; C-terminal (y/x/z): series =
        // length - index (isotope arg is always 0 — a library fragment).
        let is_nterm = matches!(type_char, 'b' | 'a' | 'c');
        let series = if is_nterm {
            f.index() as i64
        } else {
            resolved_len as i64 - f.index() as i64
        };
        if series <= 0 || series > u8::MAX as i64 {
            warn!(
                "speclib entry {:?}: fragment series {} out of range (len {}, idx {}); dropping",
                name,
                series,
                resolved_len,
                f.index()
            );
            stats.unknown_ion_dropped += 1;
            continue;
        }

        let ion = match IonAnnot::try_new(type_char, Some(series as u8), f.charge() as i8, 0) {
            Ok(ion) => ion,
            Err(e) => {
                warn!(
                    "speclib entry {:?}: failed to build IonAnnot ({:?}); dropping fragment",
                    name, e
                );
                stats.unknown_ion_dropped += 1;
                continue;
            }
        };

        if let Some(slot) = kept.iter_mut().find(|(k, _, _)| *k == ion) {
            stats.dedup_dropped += 1;
            if f.height() > slot.2 {
                slot.1 = f.mz() as f64;
                slot.2 = f.height();
            }
        } else {
            kept.push((ion, f.mz() as f64, f.height()));
        }
    }

    // (label, mz) pairs for the arena, with the parallel reference-intensity
    // sidecar filled in the same order — dedup already collapsed duplicate
    // labels, so this order is what lands in `geom.frag_labels`.
    let mut frags: Vec<(IonAnnot, f64)> = Vec::with_capacity(kept.len());
    for &(ion, mz, height) in kept.iter() {
        frags.push((ion, mz));
        frag_intens.push(height);
    }

    // Record charge is i32; `push_target` wants u8. Charge was range-checked to
    // 1..=255 above, so the cast is safe.
    geom.push_target(
        pep.mz() as f64,
        charge as u8,
        // Library iRT is dimensionless here; keep it raw (no minute->second
        // scaling) — Phase 1 RT tolerance is unrestricted.
        pep.i_rt(),
        pep.i_im(),
        &frags,
        &stripped_peptide,
        &modified_peptide,
        &[],
    );
    // `entry.protein_id` is intentionally dropped: the columnar arena has no
    // protein column.
    let _ = &entry.protein_id;
    Ok(())
}

// --- Public entry points --------------------------------------------------

/// The built arena, its parallel reference-intensity sidecar, drop statistics,
/// and whether parsing landed exactly on EOF.
type ParsedSpeclib = (
    QueryCollection<IonAnnot>,
    Vec<f32>,
    SpeclibDecodeStats,
    bool,
);

/// Parse a whole `.speclib`. A `false` in the returned EOF flag signals a
/// structural desync.
///
/// Reads the file, parses the header, then parses+maps entries in parallel
/// ([`SpecLib::open`] + [`SpecLib::parse_parallel`]).
pub fn parse_speclib_reader<R: Read>(reader: R) -> Result<ParsedSpeclib, LibraryReadingError> {
    SpecLib::open(reader)?.parse_parallel()
}

/// Cheap sniff: extension is `.speclib` and the first 4 bytes are not the
/// parquet magic `PAR1` (a `*.parquet.*.speclib` clash is handled by the parquet
/// reader). No version check here — a too-new file must still sniff true so
/// `read` returns a real `UnsupportedSpeclibVersion` diagnostic.
pub fn sniff_diann_speclib_library_file<T: AsRef<Path>>(path: T) -> bool {
    let path = path.as_ref();
    let is_speclib = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("speclib"))
        .unwrap_or(false);
    if !is_speclib {
        return false;
    }

    match File::open(path) {
        Ok(mut f) => {
            let mut magic = [0u8; 4];
            match f.read_exact(&mut magic) {
                Ok(()) => &magic != b"PAR1",
                Err(_) => false, // too short to be a real library
            }
        }
        Err(_) => false,
    }
}

/// Read a DIA-NN `.speclib` binary library directly into the columnar arena,
/// returning [`LibraryArena::Mzpaf`] with the reference-intensity sidecar.
pub fn read_diann_speclib_library_file<T: AsRef<Path>>(
    path: T,
) -> Result<LibraryArena, LibraryReadingError> {
    let path = path.as_ref();
    info!("Reading DIA-NN .speclib binary from {}", path.display());

    // mmap the file (not an owned read) — pages fault in on demand, no
    // file-sized resident buffer.
    let (mut geom, frag_intens, stats, at_eof) = SpecLib::open_mmap(path)?.parse_parallel()?;

    if !at_eof {
        // A parse that doesn't land on EOF means the entries were misaligned
        // (wrong field sizes / version gating) — the decoded library would be
        // silently corrupt. Fail loudly rather than return garbage.
        return Err(LibraryReadingError::SpeclibParse(format!(
            "parse did not land on EOF for {} — structural desync",
            path.display()
        )));
    }
    if stats.any_dropped() {
        warn!(
            "DIA-NN .speclib {}: dropped fragments (neutral-loss={}, unknown-ion={}, duplicate-label={}); discarded {} embedded decoys",
            path.display(),
            stats.loss_dropped,
            stats.unknown_ion_dropped,
            stats.dedup_dropped,
            stats.decoys_dropped,
        );
    }
    if stats.exclude_flagged > 0 {
        info!(
            "DIA-NN .speclib {}: kept {} ExcludeFromAssay-flagged fragments",
            path.display(),
            stats.exclude_flagged,
        );
    }
    info!(
        "Parsed {} precursors from DIA-NN .speclib {}",
        geom.n_rows(),
        path.display()
    );

    assert_eq!(
        frag_intens.len(),
        geom.frag_labels.len(),
        "reference-intensity sidecar must stay parallel to the fragment-label arena"
    );
    geom.seal();
    Ok(LibraryArena::Mzpaf {
        geom,
        frag_intens: Some(frag_intens),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("speclib_io_files")
            .join("diann-hela-diapasef-lib.speclib")
    }

    // All expected values below are sourced from `reference_parser.py` (an
    // independent decoder), not from this reader — so the test is not circular.

    #[test]
    fn test_sniff_diann_speclib() {
        assert!(
            sniff_diann_speclib_library_file(fixture_path()),
            "fixture should sniff as a DIA-NN .speclib"
        );
        // A non-.speclib extension must not sniff true.
        let cargo_toml = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
        assert!(!sniff_diann_speclib_library_file(cargo_toml));
    }

    #[test]
    fn test_read_fixture_lands_on_eof_with_expected_count() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (geom, frag_intens, _stats, at_eof) =
            parse_speclib_reader(std::io::BufReader::new(file)).unwrap();
        assert!(at_eof, "parse must land exactly on EOF");
        assert_eq!(geom.n_rows(), 1064, "precursor count from reference parser");
        assert_eq!(
            frag_intens.len(),
            geom.frag_labels.len(),
            "intensity sidecar stays parallel to the fragment arena"
        );
    }

    #[test]
    fn test_first_entry_fields_and_fragments() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (geom, frag_intens, _stats, _eof) =
            parse_speclib_reader(std::io::BufReader::new(file)).unwrap();

        assert_eq!(&geom.seq_mod_blob[geom.seq_mod_range(0)], "AAAGAAATHLEVAR");
        assert_eq!(
            &geom.seq_strip_blob[geom.seq_strip_range(0)],
            "AAAGAAATHLEVAR"
        );
        assert_eq!(geom.charge[0], 2);
        assert!((geom.precursor_mz[0] - 654.85541).abs() < 1e-3);
        // Library iRT kept raw (dimensionless), not scaled to seconds.
        assert!((geom.rt_seconds[0] - (-3.8674114)).abs() < 1e-3);
        assert!((geom.mobility[0] - 1.0254545).abs() < 1e-4);

        // The fixture has Peptide.length == 0, so y-series is recovered from the
        // sequence length (14). ExcludeFromAssay fragments are kept; only
        // neutral-loss/dup drop, so all 12 remain.
        let range = geom.frag_range(0);
        assert_eq!(range.len(), 12, "all non-loss fragments kept");

        // First few fragments in file order, sourced independently. y9 carries
        // the ExcludeFromAssay flag but is kept. `frag_intens` is aligned by
        // index with `frag_labels`/`frag_mzs`.
        let expected = [
            (
                IonAnnot::try_new('y', Some(8), 1, 0).unwrap(),
                896.49591_f64,
                1.0_f32,
            ),
            (
                IonAnnot::try_new('y', Some(9), 1, 0).unwrap(),
                967.53308_f64,
                0.9636101_f32,
            ),
            (
                IonAnnot::try_new('y', Some(7), 1, 0).unwrap(),
                825.45880_f64,
                0.8792037_f32,
            ),
        ];
        let labels = &geom.frag_labels[range.clone()];
        let mzs = &geom.frag_mzs[range.clone()];
        let intens = &frag_intens[range];
        for (i, (exp_ion, exp_mz, exp_int)) in expected.iter().enumerate() {
            assert_eq!(labels[i], *exp_ion, "ion label");
            assert!(
                (mzs[i] - exp_mz).abs() < 1e-2,
                "frag mz {} vs {}",
                exp_mz,
                mzs[i]
            );
            assert!((intens[i] - exp_int).abs() < 1e-4, "frag intensity");
        }
    }

    #[test]
    fn test_mid_entry_fields() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (geom, _intens, _stats, _eof) =
            parse_speclib_reader(std::io::BufReader::new(file)).unwrap();

        assert_eq!(
            &geom.seq_mod_blob[geom.seq_mod_range(532)],
            "LEGNSPQGSNQGVK"
        );
        assert_eq!(geom.charge[532], 2);
        assert!((geom.precursor_mz[532] - 707.85052).abs() < 1e-3);
        assert!((geom.rt_seconds[532] - (-31.935_83)).abs() < 1e-3);
        assert!((geom.mobility[532] - 0.9704546).abs() < 1e-4);
        assert_eq!(geom.frag_range(532).len(), 6);
    }

    #[test]
    fn test_read_library_file_yields_mzpaf_arena_with_parallel_intensities() {
        use crate::serde::read_library_file;
        let arena = read_library_file(fixture_path()).expect("read .speclib as a LibraryArena");
        match arena {
            LibraryArena::Mzpaf { geom, frag_intens } => {
                assert!(geom.n_rows() > 0, "arena must hold precursors");
                assert_eq!(
                    frag_intens.as_ref().unwrap().len(),
                    geom.frag_labels.len(),
                    "frag_intens must stay parallel to the fragment-label arena"
                );
                assert!(
                    !geom.is_decoy.iter().any(|&d| d),
                    "the reader stores targets only"
                );
            }
            LibraryArena::Str { .. } => panic!("DIA-NN .speclib must map to LibraryArena::Mzpaf"),
        }
    }

    #[test]
    fn test_exclude_flagged_kept_and_loss_dropped() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (geom, _intens, stats, _eof) =
            parse_speclib_reader(std::io::BufReader::new(file)).unwrap();
        // Reference parser: 8384 ExcludeFromAssay-flagged (kept, counted only),
        // 152 neutral-loss dropped, no dc != 0 entries.
        assert_eq!(stats.exclude_flagged, 8384);
        assert_eq!(stats.loss_dropped, 152);
        assert_eq!(stats.decoys_dropped, 0);
        // Exclude-flagged fragments are present in the output, not dropped:
        // entry[0] keeps its flagged y9.
        assert!(
            geom.frag_labels[geom.frag_range(0)]
                .iter()
                .any(|l| *l == IonAnnot::try_new('y', Some(9), 1, 0).unwrap()),
            "flagged y9 must be kept"
        );
    }
}
