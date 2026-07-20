//! Reader for the DIA-NN `.speclib` binary spectral library format.
//!
//! The format was reverse-engineered from `empirical_library.tsv.speclib` and
//! cross-checked against ProteoWizard/BiblioSpec `DiaNNSpecLibReader.cpp`
//! (`Library::read`). See the `reverse_eng_speclib` spec for the byte layout.
//!
//! The file is a sequence of variable-length sections written back-to-back
//! (little-endian, no padding, no section table); it must be parsed strictly in
//! order. This reader walks all nine top-level sections but keeps only the
//! per-precursor `Entry` payload (section 7) plus the PG `ids` strings needed to
//! resolve `protein_id`, mapping them onto the
//! [`DiannPrecursorExtras`] + [`TimsElutionGroup<IonAnnot>`] shape shared by the
//! DIA-NN TSV/parquet paths.

use super::diann_io::DiannPrecursorExtras;
use super::library_file::{
    ElutionGroupCollection,
    FileReadingExtras,
    LibraryReadingError,
};
use crate::TimsElutionGroup;
use crate::ion::IonAnnot;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use tinyvec::tiny_vec;
use tracing::{
    info,
    warn,
};

/// Newest (most-negative) `.speclib` version this reader understands. Files with
/// a version below this were written by a newer DIA-NN than was reverse-engineered.
const LATEST_SUPPORTED_VERSION: i32 = -3;

/// Bytes of one `Product`/fragment record.
const PRODUCT_SIZE: usize = 12;

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

    fn read_i32(&mut self) -> Result<i32, LibraryReadingError> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn read_f32(&mut self) -> Result<f32, LibraryReadingError> {
        Ok(f32::from_le_bytes(self.take(4)?.try_into().unwrap()))
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

/// Counts of fragments/records dropped during a decode — surfaced via `warn!`
/// and returned for the probe example to assert on.
#[derive(Debug, Default, Clone, Copy)]
pub struct SpeclibDecodeStats {
    /// Fragments flagged `ExcludeFromAssay` (`type & 0x80`). Kept (see
    /// `build_entry`); counted for reporting only.
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
    fn any_dropped(&self) -> bool {
        self.loss_dropped > 0
            || self.unknown_ion_dropped > 0
            || self.dedup_dropped > 0
            || self.decoys_dropped > 0
    }
}

/// A raw 12-byte `Product`/fragment record, decoded but not yet mapped.
struct RawFragment {
    mz: f32,
    height: f32,
    charge: u8,
    typ: u8,
    index: u8,
    loss: u8,
}

/// A decoded `Peptide` record (target or decoy). Version-gated fields default to
/// their absent value.
struct PeptideRecord {
    index: i32,
    charge: i32,
    length: i32,
    mz: f32,
    i_rt: f32,
    i_im: f32,
    frags: Vec<RawFragment>,
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

/// `Peptide::read` — the shared layout used by both the target and (when
/// `dc != 0`) the embedded decoy peptide.
fn read_peptide(c: &mut Cursor, version: i32) -> Result<PeptideRecord, LibraryReadingError> {
    let index = c.read_i32()?;
    let charge = c.read_i32()?;
    let length = c.read_i32()?;
    let mz = c.read_f32()?;
    let i_rt = c.read_f32()?;
    let _s_rt = c.read_f32()?;

    let mut i_im = 0.0f32;
    if version <= -2 {
        let _lib_qvalue = c.read_f32()?;
        i_im = c.read_f32()?;
        let _s_im = c.read_f32()?;
    }

    // The fragment block is `nfrag` fixed 12-byte records — read it in one go
    // and split into records, rather than per-field.
    let nfrag = c.read_count()?;
    let block = c.take(nfrag * PRODUCT_SIZE)?;
    let frags = block
        .chunks_exact(PRODUCT_SIZE)
        .map(|p| RawFragment {
            mz: f32::from_le_bytes(p[0..4].try_into().unwrap()),
            height: f32::from_le_bytes(p[4..8].try_into().unwrap()),
            charge: p[8],
            typ: p[9],
            index: p[10],
            loss: p[11],
        })
        .collect();

    Ok(PeptideRecord {
        index,
        charge,
        length,
        mz,
        i_rt,
        i_im,
        frags,
    })
}

/// `Entry::read` — a target `Peptide`, an optional embedded decoy, then the
/// Entry-level fields. Returns the mapped `(eg, extras)` pair for the target.
fn read_entry(
    c: &mut Cursor,
    version: i32,
    pg_ids: &[String],
    stats: &mut SpeclibDecodeStats,
) -> Result<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras), LibraryReadingError> {
    let pep = read_peptide(c, version)?;

    let dc = c.read_i32()?;
    if dc != 0 {
        // A whole second Peptide is embedded between `dc` and `entry_flags`.
        // Read it (mandatory to stay byte-synced) then discard: a lone file
        // decoy would land in its own singleton FDR group and always "win",
        // silently breaking FDR. Mass-shift decoys are generated downstream.
        let _decoy = read_peptide(c, version)?;
        stats.decoys_dropped += 1;
    }

    let _entry_flags = c.read_i32()?;
    let _proteotypic = c.read_i32()?;
    let pid_index = c.read_i32()?;
    let name = c.read_str()?;
    if version <= -3 {
        let _pg_qvalue = c.read_f32()?;
        let _ptm_qvalue = c.read_f32()?;
        let _site_conf = c.read_f32()?;
    }

    let protein_id = if pid_index >= 0 {
        pg_ids.get(pid_index as usize).cloned().unwrap_or_default()
    } else {
        String::new()
    };

    build_entry(pep, name, protein_id, stats)
}

/// Map a decoded target `PeptideRecord` + its `name` onto the pipeline shape.
fn build_entry(
    pep: PeptideRecord,
    name: String,
    protein_id: String,
    stats: &mut SpeclibDecodeStats,
) -> Result<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras), LibraryReadingError> {
    // `name` (== transition_group_id) is `modified_peptide + charge`. Strip the
    // exact known-charge decimal suffix — not "trailing digits" — so a
    // C-terminal mod ending in a digit is left intact.
    let charge_suffix = pep.charge.to_string();
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
    let resolved_len = if pep.length > 0 {
        pep.length as usize
    } else {
        residue_count(&stripped_peptide)
    };

    // (IonAnnot, fragment mz as f64, height). Dedup by IonAnnot keeping max
    // height so a duplicate label can't fail the whole load via
    // `ExpectedIntensities::try_from_pairs`.
    let mut kept: Vec<(IonAnnot, f64, f32)> = Vec::with_capacity(pep.frags.len());

    for f in &pep.frags {
        // `type & 0x80` => ExcludeFromAssay: DIA-NN's non-quantifier transitions.
        // They carry valid m/z + intensity and are kept — dropping them leaves
        // ~3 fragments/precursor on real libraries, which starves cross-library
        // calibration matching (needs >=5 shared fragments). Counted only.
        if f.typ & 0x80 != 0 {
            stats.exclude_flagged += 1;
        }
        // IonAnnot cannot represent neutral loss; drop lossy fragments.
        if f.loss != 0 {
            stats.loss_dropped += 1;
            continue;
        }

        let type_char = match f.typ & 0x7F {
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
            f.index as i64
        } else {
            resolved_len as i64 - f.index as i64
        };
        if series <= 0 || series > u8::MAX as i64 {
            warn!(
                "speclib entry {:?}: fragment series {} out of range (len {}, idx {}); dropping",
                name, series, resolved_len, f.index
            );
            stats.unknown_ion_dropped += 1;
            continue;
        }

        let ion = match IonAnnot::try_new(type_char, Some(series as u8), f.charge as i8, 0) {
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
            if f.height > slot.2 {
                slot.1 = f.mz as f64;
                slot.2 = f.height;
            }
        } else {
            kept.push((ion, f.mz as f64, f.height));
        }
    }

    let mut fragment_labels = Vec::with_capacity(kept.len());
    let mut fragment_mzs = Vec::with_capacity(kept.len());
    let mut relative_intensities = Vec::with_capacity(kept.len());
    for (ion, mz, height) in kept {
        fragment_labels.push(ion);
        fragment_mzs.push(mz);
        relative_intensities.push((ion, height));
    }

    let eg = TimsElutionGroup::builder()
        .id(pep.index as u64)
        .mobility_ook0(pep.i_im)
        // Library iRT is dimensionless here; keep it raw (no minute->second
        // scaling) — Phase 1 RT tolerance is unrestricted.
        .rt_seconds(pep.i_rt)
        .fragment_labels(fragment_labels.as_slice().into())
        .fragment_mzs(fragment_mzs)
        .precursor_labels(tiny_vec![0])
        // Record charge is i32; `precursor` wants u8. Charge is realistically
        // 1-10, so the cast is safe.
        .precursor(pep.mz as f64, pep.charge as u8)
        .try_build()
        .map_err(|e| {
            LibraryReadingError::SpeclibParse(format!(
                "failed to build elution group for {name:?}: {e:?}"
            ))
        })?;

    let extras = DiannPrecursorExtras {
        modified_peptide,
        stripped_peptide,
        protein_id,
        is_decoy: false,
        relative_intensities,
    };

    Ok((eg, extras))
}

/// Parse a whole `.speclib`. Returns the mapped entries, drop statistics, and
/// whether parsing landed exactly on EOF (a `false` here signals a structural
/// desync).
///
/// Public so the `probe_speclib` example can drive it over a synthetic buffer.
pub fn parse_speclib_reader<R: Read>(
    mut reader: R,
) -> Result<
    (
        Vec<(TimsElutionGroup<IonAnnot>, DiannPrecursorExtras)>,
        SpeclibDecodeStats,
        bool,
    ),
    LibraryReadingError,
> {
    let mut data = Vec::new();
    reader
        .read_to_end(&mut data)
        .map_err(LibraryReadingError::IoError)?;
    let mut c = Cursor::new(&data);

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
    let mut pg_ids = Vec::with_capacity(n_pg);
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

    // --- Section 7: entries — KEEP ---
    let n_entries = c.read_count()?;
    let mut stats = SpeclibDecodeStats::default();
    let mut entries = Vec::with_capacity(n_entries);
    for _ in 0..n_entries {
        entries.push(read_entry(&mut c, version, &pg_ids, &mut stats)?);
    }

    // --- Section 8: elution_groups (v <= -1, only if bytes remain) — discard ---
    if version <= -1 && !c.at_eof() {
        c.skip_vec_i32()?;
    }

    Ok((entries, stats, c.at_eof()))
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

/// Read a DIA-NN `.speclib` binary library into the shared collection type,
/// emitting `FileReadingExtras::Diann` so the DIA-NN converter is reused.
pub fn read_diann_speclib_library_file<T: AsRef<Path>>(
    path: T,
) -> Result<ElutionGroupCollection, LibraryReadingError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(LibraryReadingError::IoError)?;
    info!("Reading DIA-NN .speclib binary from {}", path.display());

    let (entries, stats, at_eof) = parse_speclib_reader(file)?;

    if !at_eof {
        warn!(
            "DIA-NN .speclib parse did not land on EOF for {} — possible structural desync",
            path.display()
        );
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
        entries.len(),
        path.display()
    );

    let (egs, extras): (Vec<_>, Vec<_>) = entries.into_iter().unzip();
    Ok(ElutionGroupCollection::MzpafLabels(
        egs,
        Some(FileReadingExtras::Diann(extras)),
    ))
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
        let (entries, _stats, at_eof) =
            parse_speclib_reader(std::io::BufReader::new(file)).unwrap();
        assert!(at_eof, "parse must land exactly on EOF");
        assert_eq!(entries.len(), 1064, "precursor count from reference parser");
    }

    #[test]
    fn test_first_entry_fields_and_fragments() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (entries, _stats, _eof) = parse_speclib_reader(std::io::BufReader::new(file)).unwrap();

        let (eg, extras) = &entries[0];
        assert_eq!(extras.modified_peptide, "AAAGAAATHLEVAR");
        assert_eq!(extras.stripped_peptide, "AAAGAAATHLEVAR");
        assert_eq!(eg.precursor_charge(), 2);
        assert!((eg.precursor_mz() - 654.85541).abs() < 1e-3);
        // Library iRT kept raw (dimensionless), not scaled to seconds.
        assert!((eg.rt_seconds() - (-3.8674114)).abs() < 1e-3);
        assert!((eg.mobility_ook0() - 1.0254545).abs() < 1e-4);

        // The fixture has Peptide.length == 0, so y-series is recovered from the
        // sequence length (14). ExcludeFromAssay fragments are kept; only
        // neutral-loss/dup drop, so all 12 remain.
        assert_eq!(eg.fragment_count(), 12, "all non-loss fragments kept");

        // First few fragments in file order, sourced independently. y9 carries
        // the ExcludeFromAssay flag but is kept.
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
        let got: Vec<(IonAnnot, f64)> = eg.iter_fragments().map(|(l, mz)| (*l, mz)).collect();
        for (exp, act) in expected.iter().zip(got.iter()) {
            assert_eq!(exp.0, act.0, "ion label");
            assert!(
                (exp.1 - act.1).abs() < 1e-2,
                "frag mz {} vs {}",
                exp.1,
                act.1
            );
        }
        // Intensities are keyed by label in relative_intensities.
        for (exp_ion, _mz, exp_int) in expected.iter() {
            let found = extras
                .relative_intensities
                .iter()
                .find(|(ion, _)| ion == exp_ion)
                .expect("fragment present in relative_intensities");
            assert!((found.1 - exp_int).abs() < 1e-4);
        }
    }

    #[test]
    fn test_mid_entry_fields() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (entries, _stats, _eof) = parse_speclib_reader(std::io::BufReader::new(file)).unwrap();

        let (eg, extras) = &entries[532];
        assert_eq!(extras.modified_peptide, "LEGNSPQGSNQGVK");
        assert_eq!(eg.precursor_charge(), 2);
        assert!((eg.precursor_mz() - 707.85052).abs() < 1e-3);
        assert!((eg.rt_seconds() - (-31.935829)).abs() < 1e-3);
        assert!((eg.mobility_ook0() - 0.9704546).abs() < 1e-4);
        assert_eq!(eg.fragment_count(), 6);
    }

    #[test]
    fn test_exclude_flagged_kept_and_loss_dropped() {
        let file = std::fs::File::open(fixture_path()).unwrap();
        let (entries, stats, _eof) = parse_speclib_reader(std::io::BufReader::new(file)).unwrap();
        // Reference parser: 8384 ExcludeFromAssay-flagged (kept, counted only),
        // 152 neutral-loss dropped, no dc != 0 entries.
        assert_eq!(stats.exclude_flagged, 8384);
        assert_eq!(stats.loss_dropped, 152);
        assert_eq!(stats.decoys_dropped, 0);
        // Exclude-flagged fragments are present in the output, not dropped:
        // entry[0] keeps its flagged y9.
        let (eg0, _) = &entries[0];
        assert!(
            eg0.iter_fragments()
                .any(|(l, _)| *l == IonAnnot::try_new('y', Some(9), 1, 0).unwrap()),
            "flagged y9 must be kept"
        );
    }
}
