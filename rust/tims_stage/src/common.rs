//! Shared internal constants + helpers for the staging impls.

use crate::error::{
    StageError,
    redact_uri,
};
use timscentroid::serialization::SerializationError;

/// Filenames required inside a `.d` for it to be openable by timsrust.
pub(crate) const REQUIRED_DOTD_FILES: &[&str] = &["analysis.tdf", "analysis.tdf_bin"];

/// TTY-aware progress bar. Returns a hidden bar in non-TTY environments so
/// stderr isn't spammed when piped to a file.
pub(crate) fn make_bar(size: u64, label: &str) -> indicatif::ProgressBar {
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

/// Wrap a `SerializationError` into a `StageError::Transport` with a
/// query-stripped URI so presigned credentials never hit the error chain.
pub(crate) fn transport_err(uri: &str) -> impl FnOnce(SerializationError) -> StageError + '_ {
    move |e| StageError::Transport {
        uri: redact_uri(uri),
        source: e,
    }
}
