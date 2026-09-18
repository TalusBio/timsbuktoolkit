//! Location-derived identity for a single sample. No content access or hashing.

/// Construct only from a location; neither the ID nor its paired display name
/// can be assigned or deserialized independently.
///
/// ```compile_fail
/// use timsseek::sample_identity::SampleIdentity;
/// let forged = SampleIdentity { sample_id: "handmade".into(), sample_name: "run".into() };
/// ```
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct SampleIdentity {
    sample_id: String,
    sample_name: String,
}

impl SampleIdentity {
    /// `location` is an absolute local path (forward slashes on Windows), or a
    /// supported remote URI. The CLI expands local paths before calling this.
    /// No symlink resolution, URI decoding, or content equivalence is attempted.
    pub fn from_location(location: &str) -> std::io::Result<Self> {
        let invalid = || {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Unable to derive sample identity from input {location:?}"),
            )
        };
        let location = location.trim_end_matches('/');
        if let Some((scheme, rest)) = location.split_once("://") {
            if !matches!(scheme, "s3" | "gs" | "az") {
                return Err(invalid());
            }
            let (bucket, key) = rest.split_once('/').ok_or_else(invalid)?;
            if bucket.is_empty() || key.is_empty() {
                return Err(invalid());
            }
        } else if !std::path::Path::new(location).is_absolute() {
            return Err(invalid());
        }
        let (parent, name) = location.rsplit_once('/').ok_or_else(invalid)?;
        let sample_name = sample_name(name).ok_or_else(invalid)?;
        let sample_id = format!(
            "{:016x}-{sample_name}",
            fnv1a64(format!("{parent}/").as_bytes())
        );
        Ok(Self {
            sample_id,
            sample_name,
        })
    }

    pub fn sample_id(&self) -> &str {
        &self.sample_id
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }
}

// Fixed FNV-1a 64-bit over UTF-8 bytes. Not DefaultHasher: its implementation
// is not a persistence contract. Parent text includes its trailing slash.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf29ce484222325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x100000001b3)
    })
}

fn sample_name(name: &str) -> Option<String> {
    // Remote keys can contain backslashes; never turn those into nested output
    // paths on Windows. Control characters are not useful display names either.
    if name.contains('\\') || name.chars().any(char::is_control) {
        return None;
    }
    let mut stem = name;
    loop {
        let before = stem;
        for ext in [".idx", ".tar", ".gz", ".d", ".raw", ".mzML", ".mzml"] {
            if let Some(s) = stem.strip_suffix(ext) {
                stem = s;
            }
        }
        if stem == before {
            break;
        }
    }
    (!matches!(stem, "" | "." | "..")).then(|| stem.to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_hash_vectors() {
        assert_eq!(fnv1a64(b""), 0xcbf29ce484222325);
        assert_eq!(fnv1a64(b"hello"), 0xa430d84680aabd0b);
        let sample = SampleIdentity::from_location("s3://bucket/rerun_1/my-run.d").unwrap();
        assert_eq!(sample.sample_id(), "0965dff92bab1aaf-my-run");
        let json = serde_json::to_value(&sample).unwrap();
        assert_eq!(json["sample_id"], sample.sample_id());
        assert_eq!(json["sample_name"], "my-run");
    }

    #[test]
    fn supported_storage_suffixes() {
        for name in ["run.d", "run.d.tar", "run.d.idx", "run.mzML.gz", "run.raw"] {
            assert_eq!(sample_name(name).as_deref(), Some("run"));
        }
        assert_eq!(sample_name("my-run.v2.d").as_deref(), Some("my-run.v2"));
        assert!(sample_name("run\\..\\other.d").is_none());
        assert!(sample_name("run\0.d").is_none());
    }
}
