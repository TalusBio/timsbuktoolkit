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
        let source = tims_stage::load::PreparedSource::new(location)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
        Self::from_source(&source)
    }

    /// Use the name from the reader selected for loading this source.
    pub fn from_source(source: &tims_stage::load::PreparedSource) -> std::io::Result<Self> {
        let location = source.uri();
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
        let (parent, _) = location.rsplit_once('/').ok_or_else(invalid)?;
        let sample_name = source.sample_name().to_owned();
        if matches!(sample_name.as_str(), "" | "." | "..")
            || sample_name.contains('\\')
            || sample_name.chars().any(char::is_control)
        {
            return Err(invalid());
        }
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
        for name in ["run.d", "run.d.tar", "run.d.idx", "run.idx", "run.tar"] {
            let identity = SampleIdentity::from_location(&format!("s3://bucket/{name}")).unwrap();
            assert_eq!(identity.sample_name(), "run");
        }
        for name in [
            "run\\..\\other.d",
            "run\0.d",
            ".d",
            ".idx",
            "run.raw",
            "run.mzML.gz",
        ] {
            assert!(SampleIdentity::from_location(&format!("s3://bucket/{name}")).is_err());
        }
    }

    #[test]
    fn suffix_case_is_ignored_without_changing_stem_or_parent_case() {
        let expected = SampleIdentity::from_location("s3://bucket/Batch/My-Run.d").unwrap();
        for suffix in [".D", ".D.TaR", ".D.IdX/"] {
            let location = format!("s3://bucket/Batch/My-Run{suffix}");
            assert_eq!(SampleIdentity::from_location(&location).unwrap(), expected);
        }
        assert_eq!(expected.sample_name(), "My-Run");
        for location in ["s3://bucket/batch/My-Run.d", "s3://bucket/Batch/my-run.d"] {
            assert_ne!(SampleIdentity::from_location(location).unwrap(), expected);
        }
        assert_eq!(
            SampleIdentity::from_location(
                std::env::current_dir()
                    .unwrap()
                    .join("Échantillon.D")
                    .to_str()
                    .unwrap()
            )
            .unwrap()
            .sample_name(),
            "Échantillon"
        );
    }
}
