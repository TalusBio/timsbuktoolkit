//! TarReader trait + ustar walker + stage_{s3,local}_tar helpers.

use crate::backend::{
    PerRunTempdir,
    StagedDotD,
};
use crate::error::{
    StageError,
    redact_uri,
};
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
const REQUIRED: &[&str] = &["analysis.tdf", "analysis.tdf_bin"];

/// Abstraction over "something we can pull byte ranges out of". S3 impl uses
/// range-GETs; local impl uses `File::seek + read_exact`.
pub trait TarReader {
    fn read_range(&mut self, offset: u64, len: usize) -> Result<Bytes, StageError>;
    fn total_len(&self) -> u64;
}

/// Range-GET backed reader over an S3 object. Prefetches the first 64 KiB at
/// construction so walks that fit in the prefetch window complete in one
/// round-trip.
pub struct S3TarReader {
    provider: StorageProvider,
    key: String,
    size: u64,
    /// Prefix buffer covering `[0, prefix.len())` bytes of the object.
    prefix: Bytes,
}

impl S3TarReader {
    pub fn new(loc: &StorageLocation, key: &str) -> Result<Self, StageError> {
        let uri_for_err = format!("{loc:?}/{key}");
        let provider = StorageProvider::open(loc.clone()).map_err(|e| StageError::Transport {
            uri: redact_uri(&uri_for_err),
            source: e,
        })?;
        let meta = provider.head(key).map_err(|e| StageError::Transport {
            uri: redact_uri(&uri_for_err),
            source: e,
        })?;
        let size = meta.size;
        let window = std::cmp::min(PREFETCH as u64, size);
        let prefix = if window == 0 {
            Bytes::new()
        } else {
            provider
                .range_get(key, 0..window)
                .map_err(|e| StageError::Transport {
                    uri: redact_uri(&uri_for_err),
                    source: e,
                })?
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
            .map_err(|e| StageError::Transport {
                uri: redact_uri(&self.key),
                source: e,
            })?;
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
}

/// File-backed reader over a local `.tar`.
pub struct LocalTarReader {
    file: std::fs::File,
    size: u64,
}

impl LocalTarReader {
    pub fn new(path: &Path) -> Result<Self, StageError> {
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
        write_range_to_file(reader, offset, size, &dotd.join(bn), &bar)?;
        bar.finish_and_clear();
    }
    Ok(StagedDotD::owned(tempdir, dotd))
}

fn write_range_to_file(
    reader: &mut dyn TarReader,
    offset: u64,
    size: u64,
    dst: &Path,
    bar: &indicatif::ProgressBar,
) -> Result<(), StageError> {
    use std::io::Write;
    let mut f = std::fs::File::create(dst).map_err(StageError::Io)?;
    const CHUNK: usize = 4 * 1024 * 1024;
    let mut written = 0u64;
    while written < size {
        let take = std::cmp::min(CHUNK as u64, size - written) as usize;
        let chunk = reader.read_range(offset + written, take)?;
        f.write_all(&chunk).map_err(StageError::Io)?;
        bar.inc(chunk.len() as u64);
        written += chunk.len() as u64;
    }
    Ok(())
}

fn make_bar(size: u64, label: &str) -> indicatif::ProgressBar {
    use std::io::IsTerminal;
    if !std::io::stderr().is_terminal() {
        return indicatif::ProgressBar::hidden();
    }
    let bar = indicatif::ProgressBar::new(size);
    bar.set_style(
        indicatif::ProgressStyle::with_template(
            "{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})",
        )
        .unwrap(),
    );
    bar.set_message(label.to_string());
    bar
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
