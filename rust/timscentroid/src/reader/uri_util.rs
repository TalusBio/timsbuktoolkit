//! Small path-edit helpers for `http::Uri`.
//!
//! `http::Uri` is RFC-3986 and has no path-mutation API, so building sidecar /
//! child URIs from a primary requires rebuilding through `Uri::into_parts`.
//! `child` populates directory manifests (Bruker `analysis.tdf`); `with_suffix`
//! and `sibling` build vendor sidecars (SCIEX `.wiff.scan`, `.timeseries.data`).

use http::Uri;
use http::uri::PathAndQuery;

/// Rebuild `base` with its path replaced by `new_path`, preserving scheme +
/// authority. `new_path` must be a valid URI path (percent-encoded, no spaces).
fn with_path(base: &Uri, new_path: &str) -> Uri {
    let mut parts = base.clone().into_parts();
    parts.path_and_query = Some(
        new_path
            .parse::<PathAndQuery>()
            .unwrap_or_else(|e| panic!("invalid derived URI path {new_path:?}: {e}")),
    );
    Uri::from_parts(parts).expect("reassembled URI is valid")
}

/// Append `/segment` to the URI path — a child within a directory URI.
/// `s3://bkt/sample.d` + `analysis.tdf` → `s3://bkt/sample.d/analysis.tdf`.
pub fn child(base: &Uri, segment: &str) -> Uri {
    let path = base.path().trim_end_matches('/');
    with_path(base, &format!("{path}/{segment}"))
}

/// Append `suffix` directly to the URI path (no separator).
/// `.../foo.wiff` + `.scan` → `.../foo.wiff.scan`; + `2` → `.../foo.wiff2`.
pub fn with_suffix(base: &Uri, suffix: &str) -> Uri {
    let path = base.path().trim_end_matches('/');
    with_path(base, &format!("{path}{suffix}"))
}

/// Replace the last path segment with `name` — a co-located sibling artifact.
/// `.../a/foo.wiff` + `bar.data` → `.../a/bar.data`.
pub fn sibling(base: &Uri, name: &str) -> Uri {
    let path = base.path().trim_end_matches('/');
    let parent = path.rsplit_once('/').map(|(p, _)| p).unwrap_or("");
    with_path(base, &format!("{parent}/{name}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn child_appends_segment_local() {
        let u: Uri = "/data/sample.d".parse().unwrap();
        assert_eq!(
            child(&u, "analysis.tdf").path(),
            "/data/sample.d/analysis.tdf"
        );
    }

    #[test]
    fn child_appends_segment_s3() {
        let u: Uri = "s3://bkt/sample.d".parse().unwrap();
        let c = child(&u, "analysis.tdf_bin");
        assert_eq!(c.scheme_str(), Some("s3"));
        assert_eq!(c.authority().unwrap().as_str(), "bkt");
        assert_eq!(c.path(), "/sample.d/analysis.tdf_bin");
    }

    #[test]
    fn with_suffix_appends_no_separator() {
        let u: Uri = "/data/run.wiff".parse().unwrap();
        assert_eq!(with_suffix(&u, ".scan").path(), "/data/run.wiff.scan");
        assert_eq!(with_suffix(&u, "2").path(), "/data/run.wiff2");
    }

    #[test]
    fn sibling_replaces_last_segment() {
        let u: Uri = "/data/a/run.wiff".parse().unwrap();
        assert_eq!(
            sibling(&u, "run.timeseries.data").path(),
            "/data/a/run.timeseries.data"
        );
    }
}
