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
use timscentroid::StorageProvider;

pub(crate) fn stage_s3_prefix(
    backend: &PerRunTempdir,
    spec: &SourceSpec,
) -> Result<StagedDotD, StageError> {
    let SourceSpec::S3Prefix { loc, prefix } = spec else {
        unreachable!()
    };
    let step = timscentroid::TimedStep::begin("Staging .d (prefix)");
    let uri_for_err = format!("{loc:?}/{prefix}");

    // Read-only: do not create local dirs.
    let provider = StorageProvider::open(loc.clone()).map_err(transport_err(&uri_for_err))?;

    let listing = provider
        .list_capped(prefix, backend.config().max_prefix_keys)
        .map_err(|e| match e {
            timscentroid::serialization::SerializationError::PrefixCapExceeded { cap, prefix } => {
                StageError::PrefixCapExceeded { cap, prefix }
            }
            other => StageError::Transport {
                uri: crate::error::redact_uri(&uri_for_err),
                source: other,
            },
        })?;

    // Filter to required basenames.
    let mut by_name: std::collections::HashMap<&'static str, _> = Default::default();
    for meta in &listing {
        let full = meta.location.to_string();
        let bn = full.rsplit('/').next().unwrap_or(full.as_str()).to_string();
        if let Some(r) = REQUIRED.iter().find(|r| ***r == bn) {
            by_name.insert(*r, meta);
        }
    }
    let missing: Vec<String> = REQUIRED
        .iter()
        .filter(|r| !by_name.contains_key(*r))
        .map(|s| s.to_string())
        .collect();
    if !missing.is_empty() {
        return Err(StageError::MissingRequiredFiles {
            missing,
            found: listing.iter().map(|m| m.location.to_string()).collect(),
        });
    }

    let tempdir = backend.new_tempdir()?;
    let dotd = tempdir.path().join("sample.d");
    std::fs::create_dir_all(&dotd).map_err(StageError::Io)?;
    std::fs::File::create(tempdir.path().join(".lock")).map_err(StageError::Io)?;

    for bn in REQUIRED {
        let meta = by_name.get(bn).expect("preflight guarantees presence");
        let bar = make_bar(meta.size, bn);
        let full_key = meta.location.to_string();
        provider
            .get_to_file(&full_key, &dotd.join(bn), &bar)
            .map_err(transport_err(&full_key))?;
        bar.finish_and_clear();
    }
    step.finish();
    Ok(StagedDotD::owned(tempdir, dotd))
}
