//! Stable sample identity from the input location, never from its contents.

use std::collections::HashMap;

use crate::errors::CliError;
use timsseek::sample_identity::SampleIdentity;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SampleInput {
    uri: String,
    identity: SampleIdentity,
}

impl SampleInput {
    pub(crate) fn uri(&self) -> &str {
        &self.uri
    }

    pub(crate) fn identity(&self) -> &SampleIdentity {
        &self.identity
    }

    pub(crate) fn sample_id(&self) -> &str {
        self.identity.sample_id()
    }

    pub(crate) fn sample_name(&self) -> &str {
        self.identity.sample_name()
    }
}

/// Resolve once, before staging changes an input's location. Local paths become
/// absolute without resolving symlinks; remote URI text is preserved. The parent
/// includes its trailing slash. Supported storage suffixes are not identity.
pub(crate) fn resolve_samples(uris: &[String]) -> Result<Vec<SampleInput>, CliError> {
    let mut seen = HashMap::new();
    let mut samples = Vec::with_capacity(uris.len());
    for uri in uris {
        let expanded = tims_stage::expand_local_uri(uri);
        let invalid = || CliError::Config {
            source: format!("Unable to derive sample identity from input {uri:?}"),
        };
        if expanded
            .trim_end_matches(['/', '\\'])
            .rsplit(['/', '\\'])
            .next()
            .is_none_or(|name| matches!(name, "" | "." | ".."))
        {
            return Err(invalid());
        }
        let location = if tims_stage::is_remote_uri(&expanded) {
            expanded
        } else {
            let absolute = std::path::absolute(&expanded).map_err(|e| CliError::Io {
                source: format!("Resolving sample input: {e}"),
                path: Some(uri.clone()),
            })?;
            let text = absolute.to_str().ok_or_else(invalid)?;
            if cfg!(windows) {
                text.replace('\\', "/")
            } else {
                text.to_owned()
            }
        };
        let identity = SampleIdentity::from_location(&location).map_err(|e| CliError::Config {
            source: e.to_string(),
        })?;
        let sample_id = identity.sample_id();
        if let Some(previous) = seen.insert(sample_id.to_owned(), uri) {
            return Err(CliError::Config {
                source: format!(
                    "Duplicate sample_id {sample_id:?} for inputs {previous:?} and {uri:?}; \
                     duplicate identities are unsupported, including with --overwrite"
                ),
            });
        }
        samples.push(SampleInput {
            uri: uri.clone(),
            identity,
        });
    }
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn resolve(uris: &[&str]) -> Result<Vec<SampleInput>, CliError> {
        resolve_samples(&uris.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }

    #[test]
    fn independent_runs_and_reordering_keep_identity() {
        let uris = [
            "s3://bucket/rerun_1/my-run.d",
            "s3://bucket/rerun_2/my-run.d",
        ];
        let samples = resolve(&uris).unwrap();
        assert_ne!(samples[0].sample_id(), samples[1].sample_id());
        assert_eq!(samples[0].sample_name(), "my-run");
        assert_eq!(samples[0], resolve(&uris[..1]).unwrap()[0]);
        assert_eq!(samples[0], resolve(&[uris[1], uris[0]]).unwrap()[1]);
    }

    #[test]
    fn storage_forms_share_identity_and_collide() {
        let reference = resolve(&["s3://bucket/run.d"]).unwrap().remove(0);
        for suffix in [
            ".d/", ".d.tar", ".d.idx/", ".raw", ".mzML.gz", ".D.IDX/", ".RaW", ".MzMl.GZ",
        ] {
            let uri = format!("s3://bucket/run{suffix}");
            assert_eq!(
                resolve(&[&uri]).unwrap()[0].sample_id(),
                reference.sample_id()
            );
            let error = resolve(&[&reference.uri, &uri]).unwrap_err().to_string();
            assert!(error.contains(&reference.uri));
            assert!(error.contains(&uri));
            assert!(error.contains("Duplicate sample_id"));
        }
    }

    #[test]
    fn relative_and_absolute_local_paths_agree_without_opening_input() {
        let relative = "a1-nonexistent/run.d";
        let absolute = std::env::current_dir().unwrap().join(relative);
        assert_eq!(
            resolve(&[relative]).unwrap()[0].sample_id(),
            resolve(&[absolute.to_str().unwrap()]).unwrap()[0].sample_id()
        );
        assert!(resolve(&[relative, absolute.to_str().unwrap()]).is_err());
        let others = resolve(&["a1-nonexistent/run.d", "another-nonexistent/run.d"]).unwrap();
        assert_ne!(others[0].sample_id(), others[1].sample_id());
    }

    #[test]
    fn missing_or_unsafe_names_are_rejected() {
        for uri in [
            "",
            "/",
            ".",
            "..",
            "s3://bucket/",
            "s3:///run.d",
            "s3://bucket/.d",
            "s3://bucket/..",
        ] {
            assert!(resolve(&[uri]).is_err(), "{uri}");
        }
    }
}
