use crate::backend::{
    PerRunTempdir,
    StagedDotD,
};
use crate::common::{
    REQUIRED_DOTD_FILES as REQUIRED,
    make_bar,
    transport_err,
};
use crate::error::StageError;
use crate::resolve::SourceSpec;
use bytes::Bytes;
use std::collections::BTreeMap;
use std::io::{
    Read,
    Seek,
    SeekFrom,
};
use std::path::Path;
use timscentroid::{
    StorageLocation,
    StorageProvider,
};

const BLOCK: usize = 512;
const PREFETCH: usize = 64 * 1024;

/// Abstraction over "something we can pull byte ranges out of".
///
/// `read_range` is for small header reads (up to a few KiB). For large
/// payload copies (often GBs), use `copy_range_to_file` which streams — the
/// S3 impl routes that through one range-GET whose response body writes
/// directly to disk. Never loop `read_range` for big payloads: each call
/// is a separate HTTP GET on S3.
pub(crate) trait TarReader {
    fn read_range(&mut self, offset: u64, len: usize) -> Result<Bytes, StageError>;
    fn total_len(&self) -> u64;
    fn copy_range_to_file(
        &mut self,
        offset: u64,
        size: u64,
        dst: &Path,
        bar: &indicatif::ProgressBar,
    ) -> Result<(), StageError>;
}

/// Range-GET backed reader over an S3 object. Prefetches the first 64 KiB at
/// construction so walks that fit in the prefetch window complete in one
/// round-trip.
pub(crate) struct S3TarReader {
    provider: StorageProvider,
    key: String,
    size: u64,
    /// Prefix buffer covering `[0, prefix.len())` bytes of the object.
    prefix: Bytes,
}

impl S3TarReader {
    pub(crate) fn new(loc: &StorageLocation, key: &str) -> Result<Self, StageError> {
        let uri_for_err = format!("{loc:?}/{key}");
        let provider = StorageProvider::open(loc.clone()).map_err(transport_err(&uri_for_err))?;
        let meta = provider.head(key).map_err(transport_err(&uri_for_err))?;
        let size = meta.size;
        let window = std::cmp::min(PREFETCH as u64, size);
        let prefix = if window == 0 {
            Bytes::new()
        } else {
            provider
                .range_get(key, 0..window)
                .map_err(transport_err(&uri_for_err))?
        };
        Ok(Self {
            provider,
            key: key.to_string(),
            size,
            prefix,
        })
    }
}

impl TarReader for S3TarReader {
    fn read_range(&mut self, offset: u64, len: usize) -> Result<Bytes, StageError> {
        let end = offset
            .checked_add(len as u64)
            .ok_or_else(|| StageError::ShapeMismatch("offset overflow".into()))?;
        if end <= self.prefix.len() as u64 {
            let start = offset as usize;
            return Ok(self.prefix.slice(start..start + len));
        }
        let got = self
            .provider
            .range_get(&self.key, offset..end)
            .map_err(transport_err(&self.key))?;
        if got.len() < len {
            return Err(StageError::ShortRead {
                expected: len as u64,
                actual: got.len() as u64,
            });
        }
        Ok(got)
    }

    fn total_len(&self) -> u64 {
        self.size
    }

    fn copy_range_to_file(
        &mut self,
        offset: u64,
        size: u64,
        dst: &Path,
        bar: &indicatif::ProgressBar,
    ) -> Result<(), StageError> {
        let end = offset
            .checked_add(size)
            .ok_or_else(|| StageError::ShapeMismatch("offset+size overflow".into()))?;
        self.provider
            .range_get_to_file(&self.key, offset..end, dst, bar)
            .map_err(transport_err(&self.key))
    }
}

/// File-backed reader over a local `.tar`.
pub(crate) struct LocalTarReader {
    file: std::fs::File,
    size: u64,
}

impl LocalTarReader {
    pub(crate) fn new(path: &Path) -> Result<Self, StageError> {
        let file = std::fs::File::open(path).map_err(StageError::Io)?;
        let size = file.metadata().map_err(StageError::Io)?.len();
        Ok(Self { file, size })
    }
}

impl TarReader for LocalTarReader {
    fn read_range(&mut self, offset: u64, len: usize) -> Result<Bytes, StageError> {
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(StageError::Io)?;
        let mut buf = vec![0u8; len];
        self.file.read_exact(&mut buf).map_err(StageError::Io)?;
        Ok(Bytes::from(buf))
    }

    fn total_len(&self) -> u64 {
        self.size
    }

    fn copy_range_to_file(
        &mut self,
        offset: u64,
        size: u64,
        dst: &Path,
        bar: &indicatif::ProgressBar,
    ) -> Result<(), StageError> {
        use std::io::Write;
        if let Some(parent) = dst.parent() {
            std::fs::create_dir_all(parent).map_err(StageError::Io)?;
        }
        let mut out = std::fs::File::create(dst).map_err(StageError::Io)?;
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(StageError::Io)?;
        let mut remaining = size;
        let mut buf = vec![0u8; 4 * 1024 * 1024];
        while remaining > 0 {
            let take = std::cmp::min(buf.len() as u64, remaining) as usize;
            self.file
                .read_exact(&mut buf[..take])
                .map_err(StageError::Io)?;
            out.write_all(&buf[..take]).map_err(StageError::Io)?;
            bar.inc(take as u64);
            remaining -= take as u64;
        }
        Ok(())
    }
}

// -------- ustar header parsing --------

#[derive(Debug, Clone)]
struct TarHeader {
    name: String,
    size: u64,
    typeflag: u8,
}

fn parse_header(block: &[u8]) -> Result<TarHeader, StageError> {
    if block.len() != BLOCK {
        return Err(StageError::ShapeMismatch(format!(
            "expected 512 B header, got {}",
            block.len()
        )));
    }
    let name = {
        let raw = &block[0..100];
        let end = raw.iter().position(|&b| b == 0).unwrap_or(raw.len());
        std::str::from_utf8(&raw[..end])
            .map_err(|e| StageError::ShapeMismatch(format!("header name not UTF-8: {e}")))?
            .to_string()
    };
    let size = parse_octal(&block[124..136])?;
    let typeflag = block[156];
    Ok(TarHeader {
        name,
        size,
        typeflag,
    })
}

fn parse_octal(raw: &[u8]) -> Result<u64, StageError> {
    let s = std::str::from_utf8(raw)
        .map_err(|e| StageError::ShapeMismatch(format!("octal field not ASCII: {e}")))?
        .trim_end_matches(|c: char| c == '\0' || c == ' ')
        .trim_start_matches(' ');
    if s.is_empty() {
        return Ok(0);
    }
    u64::from_str_radix(s, 8)
        .map_err(|e| StageError::ShapeMismatch(format!("bad octal `{s}`: {e}")))
}

fn is_zero_block(b: &[u8]) -> bool {
    b.iter().all(|&x| x == 0)
}

fn basename_of(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

/// Walk headers until all required basenames are located, or EOA is reached.
/// Returns a map `basename -> (payload_offset, payload_size)`.
///
/// Typeflag handling:
///   - `'0'` / `0`  : regular file — record if basename is required.
///   - `'x'` / `'g'`: pax extended header — skip payload.
///   - `'L'` / `'K'`: GNU long-name record — hard error.
///   - anything else: unknown — skip payload, don't record.
///
/// Permissive-by-default; the corpus-inventory step will tell us which
/// branches can be pruned later.
fn walk_tar_index(reader: &mut dyn TarReader) -> Result<BTreeMap<String, (u64, u64)>, StageError> {
    let mut out: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    let total = reader.total_len();
    let mut offset: u64 = 0;
    let mut zero_run = 0;
    while offset + BLOCK as u64 <= total {
        let block = reader.read_range(offset, BLOCK)?;
        if is_zero_block(&block) {
            zero_run += 1;
            if zero_run >= 2 {
                break;
            }
            offset += BLOCK as u64;
            continue;
        }
        zero_run = 0;
        let h = parse_header(&block)?;
        let payload_offset = offset + BLOCK as u64;
        let padded = ((h.size + BLOCK as u64 - 1) / BLOCK as u64) * BLOCK as u64;
        match h.typeflag {
            b'0' | 0 => {
                let bn = basename_of(&h.name).to_string();
                if REQUIRED.contains(&bn.as_str()) {
                    out.insert(bn, (payload_offset, h.size));
                }
            }
            b'x' | b'g' => { /* pax extended header — skip payload */ }
            b'L' | b'K' => return Err(StageError::UnsupportedTarFeature),
            _ => { /* unknown typeflag — skip payload */ }
        }
        offset = payload_offset + padded;
        if REQUIRED.iter().all(|r| out.contains_key(*r)) {
            break;
        }
    }
    let missing: Vec<String> = REQUIRED
        .iter()
        .filter(|r| !out.contains_key(**r))
        .map(|s| s.to_string())
        .collect();
    if !missing.is_empty() {
        return Err(StageError::MissingRequiredFiles {
            missing,
            found: out.keys().cloned().collect(),
        });
    }
    Ok(out)
}

// -------- stage helpers --------

pub(crate) fn stage_s3_tar(
    backend: &PerRunTempdir,
    spec: &SourceSpec,
) -> Result<StagedDotD, StageError> {
    let SourceSpec::S3Tar { loc, key } = spec else {
        unreachable!()
    };
    let step = timscentroid::TimedStep::begin("Staging .d (tar)");
    let mut reader = S3TarReader::new(loc, key)?;
    let index = walk_tar_index(&mut reader)?;
    let staged = materialize(backend, &mut reader, &index)?;
    step.finish();
    Ok(staged)
}

pub(crate) fn stage_local_tar(
    backend: &PerRunTempdir,
    spec: &SourceSpec,
) -> Result<StagedDotD, StageError> {
    let SourceSpec::LocalTar { path } = spec else {
        unreachable!()
    };
    let step = timscentroid::TimedStep::begin("Staging .d (local tar)");
    let mut reader = LocalTarReader::new(path)?;
    let index = walk_tar_index(&mut reader)?;
    let staged = materialize(backend, &mut reader, &index)?;
    step.finish();
    Ok(staged)
}

fn materialize(
    backend: &PerRunTempdir,
    reader: &mut dyn TarReader,
    index: &BTreeMap<String, (u64, u64)>,
) -> Result<StagedDotD, StageError> {
    let tempdir = backend.new_tempdir()?;
    let dotd = tempdir.path().join("sample.d");
    std::fs::create_dir_all(&dotd).map_err(StageError::Io)?;
    std::fs::File::create(tempdir.path().join(".lock")).map_err(StageError::Io)?;
    for bn in REQUIRED {
        let (offset, size) = *index.get(*bn).expect("preflight guarantees presence");
        let bar = make_bar(size, bn);
        reader.copy_range_to_file(offset, size, &dotd.join(bn), &bar)?;
        bar.finish_and_clear();
    }
    Ok(StagedDotD::owned(tempdir, dotd))
}

// -------- walker unit tests --------

#[cfg(test)]
mod walker_tests {
    use super::*;

    fn build_ustar_header(name: &str, size: u64, typeflag: u8) -> Vec<u8> {
        let mut buf = vec![0u8; BLOCK];
        let nb = name.as_bytes();
        buf[..nb.len()].copy_from_slice(nb);
        let s = format!("{:011o}", size);
        buf[124..124 + s.len()].copy_from_slice(s.as_bytes());
        buf[135] = 0;
        buf[156] = typeflag;
        for i in 148..156 {
            buf[i] = b' ';
        }
        buf
    }

    fn pad_to_block(payload: &mut Vec<u8>) {
        let pad = (BLOCK - payload.len() % BLOCK) % BLOCK;
        payload.extend(std::iter::repeat(0u8).take(pad));
    }

    struct InMemTar {
        data: Vec<u8>,
    }

    impl TarReader for InMemTar {
        fn read_range(&mut self, offset: u64, len: usize) -> Result<Bytes, StageError> {
            let o = offset as usize;
            if o + len > self.data.len() {
                return Err(StageError::ShortRead {
                    expected: len as u64,
                    actual: (self.data.len() - o) as u64,
                });
            }
            Ok(Bytes::copy_from_slice(&self.data[o..o + len]))
        }

        fn total_len(&self) -> u64 {
            self.data.len() as u64
        }

        fn copy_range_to_file(
            &mut self,
            offset: u64,
            size: u64,
            dst: &Path,
            bar: &indicatif::ProgressBar,
        ) -> Result<(), StageError> {
            let bytes = self.read_range(offset, size as usize)?;
            std::fs::write(dst, &bytes).map_err(StageError::Io)?;
            bar.inc(bytes.len() as u64);
            Ok(())
        }
    }

    fn minimal_tar(entries: &[(&str, Vec<u8>, u8)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (name, payload, typeflag) in entries {
            out.extend_from_slice(&build_ustar_header(name, payload.len() as u64, *typeflag));
            out.extend_from_slice(payload);
            pad_to_block(&mut out);
        }
        out.extend(std::iter::repeat(0u8).take(BLOCK * 2));
        out
    }

    #[test]
    fn walker_finds_required_basenames_flat_tar() {
        let data = minimal_tar(&[
            ("analysis.tdf", vec![1u8; 10], b'0'),
            ("analysis.tdf_bin", vec![2u8; 4096], b'0'),
        ]);
        let mut r = InMemTar { data };
        let idx = walk_tar_index(&mut r).unwrap();
        assert!(idx.contains_key("analysis.tdf"));
        assert!(idx.contains_key("analysis.tdf_bin"));
    }

    #[test]
    fn walker_finds_basenames_nested_in_sample_dot_d() {
        let data = minimal_tar(&[
            ("sample.d/analysis.tdf", vec![1u8; 10], b'0'),
            ("sample.d/analysis.tdf_bin", vec![2u8; 4096], b'0'),
        ]);
        let mut r = InMemTar { data };
        let idx = walk_tar_index(&mut r).unwrap();
        assert_eq!(idx.len(), 2);
    }

    #[test]
    fn walker_skips_pax_x_header() {
        let data = minimal_tar(&[
            ("PaxHeader/extra", vec![1u8; 32], b'x'),
            ("analysis.tdf", vec![1u8; 10], b'0'),
            ("analysis.tdf_bin", vec![2u8; 4096], b'0'),
        ]);
        let mut r = InMemTar { data };
        let idx = walk_tar_index(&mut r).unwrap();
        assert_eq!(idx.len(), 2);
    }

    #[test]
    fn walker_errors_on_gnu_longname() {
        let data = minimal_tar(&[
            ("LongNameBlock", vec![1u8; 32], b'L'),
            ("analysis.tdf", vec![1u8; 10], b'0'),
            ("analysis.tdf_bin", vec![2u8; 4096], b'0'),
        ]);
        let mut r = InMemTar { data };
        let err = walk_tar_index(&mut r).unwrap_err();
        assert!(matches!(err, StageError::UnsupportedTarFeature));
    }

    #[test]
    fn walker_errors_on_missing_basenames() {
        let data = minimal_tar(&[("some_other_file", vec![1u8; 10], b'0')]);
        let mut r = InMemTar { data };
        let err = walk_tar_index(&mut r).unwrap_err();
        assert!(matches!(err, StageError::MissingRequiredFiles { .. }));
    }
}
