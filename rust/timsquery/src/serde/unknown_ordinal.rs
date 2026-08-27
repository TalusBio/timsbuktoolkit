//! Shared counter for `?`-labelled unknown ions.

/// Advance a per-precursor `?`-ordinal counter, refusing to wrap.
///
/// `IonAnnot` ordinals are `u8`, so a precursor can carry at most
/// [`u8::MAX`] distinguishable unknown ions. Past that there is no way to
/// keep labels unique, and duplicates are the one thing the arena cannot
/// tolerate — see [`ExpectedIntensities::try_from_pairs`] for why.
///
/// `None` means "cannot allocate another": callers turn it into their own
/// over-capacity error while the row index is still in hand, rather than
/// letting the duplicate surface downstream.
///
/// [`ExpectedIntensities::try_from_pairs`]: https://docs.rs/timsseek
pub(crate) fn next_unknown_ordinal(current: u8) -> Option<u8> {
    current.checked_add(1)
}
