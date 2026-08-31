//! Source ids, as opposed to arena positions.
//!
//! The arena position is self-incremental and feeds decoy grouping and the
//! q-value determinism key, so it stays private to the arena. A source id is
//! what the source file called the row, carried through verbatim.
//!
//! Ids are polymorphic because formats disagree: JSON targets carry integers,
//! DIA-NN carries `transition_group_id` / the `.speclib` entry name, which are
//! strings.
//!
//! The shape is preserved where a consumer distinguishes the two, which today
//! is the JSON result path only: Carafe reads `id` with fastjson and keys
//! results by it, so `"id": 7` must stay a bare number. The parquet writer does
//! not distinguish -- both id columns are `Utf8` and a numeric id is written as
//! its digits -- because a search result is read by column type, not by the
//! shape of one value.

use serde::{
    Deserialize,
    Serialize,
};
use thiserror::Error;

/// One row's id, borrowed from the arena that stores it.
///
/// Deliberately not `Serialize`: the result path writes an [`OwnedSourceId`]
/// (`ChromatogramOutput.id`), and having both be serializable meant the wire
/// contract could be pinned against whichever one a test reached for. Call
/// [`Self::to_owned_id`] to serialize.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SourceId<'a> {
    Numeric(u64),
    Text(&'a str),
}

impl SourceId<'_> {
    pub fn to_owned_id(self) -> OwnedSourceId {
        match self {
            Self::Numeric(n) => OwnedSourceId::Numeric(n),
            Self::Text(s) => OwnedSourceId::Text(s.to_string()),
        }
    }
}

impl std::fmt::Display for SourceId<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Numeric(n) => write!(f, "{n}"),
            Self::Text(s) => f.write_str(s),
        }
    }
}

/// [`SourceId`] detached from the arena, for result structs that outlive the
/// borrow.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(untagged)]
pub enum OwnedSourceId {
    Numeric(u64),
    Text(String),
}

impl OwnedSourceId {
    pub fn as_ref(&self) -> SourceId<'_> {
        match self {
            Self::Numeric(n) => SourceId::Numeric(*n),
            Self::Text(s) => SourceId::Text(s),
        }
    }

    /// Overwrite in place, reusing the existing `String` capacity when both
    /// sides are text. `Target::reset_from` runs per query on the scoring hot
    /// path, so a fresh allocation there is a per-query cost.
    pub fn set_from(&mut self, src: SourceId<'_>) {
        match (&mut *self, src) {
            (Self::Text(dst), SourceId::Text(s)) => {
                dst.clear();
                dst.push_str(s);
            }
            (_, other) => *self = other.to_owned_id(),
        }
    }
}

impl From<u64> for OwnedSourceId {
    fn from(n: u64) -> Self {
        Self::Numeric(n)
    }
}

impl From<String> for OwnedSourceId {
    fn from(s: String) -> Self {
        Self::Text(s)
    }
}

impl From<&str> for OwnedSourceId {
    fn from(s: &str) -> Self {
        Self::Text(s.to_string())
    }
}

impl OwnedSourceId {
    /// The id of a scratch buffer that has not been filled yet (see
    /// [`crate::Target::empty_like`]).
    ///
    /// Allocates nothing: an empty `String` owns no buffer until something is
    /// pushed into it, and [`Self::set_from`] then grows it once -- exactly what
    /// a fresh `String` would cost. A named placeholder like `"<unset>"` would
    /// allocate on every scratch buffer to store bytes nothing reads.
    ///
    /// Deliberately not `Default`, and deliberately not `Numeric(0)`: every
    /// `u64` is a valid id, so a numeric placeholder that leaked into a result
    /// would look exactly like a row keyed on 0 -- the same reason
    /// `RowIdx::default()` is `u32::MAX` rather than 0. An empty id is visibly
    /// unset.
    pub fn placeholder() -> Self {
        Self::Text(String::new())
    }
}

impl std::fmt::Display for OwnedSourceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

/// Which kind of id a library carries is a property of the format, not of a
/// row: JSON targets carry integers, DIA-NN carries strings
/// (`transition_group_id`, or the `.speclib` entry name), Spectronaut/Skyline
/// carry nothing.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum SourceIds {
    #[default]
    Absent,
    Numeric(Vec<u64>),
    /// CSR over one blob, like `seq_strip_blob`: `offsets` has `n_rows + 1`
    /// entries.
    Text {
        blob: String,
        offsets: Vec<u32>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SourceIdError {
    #[error("{ids} source ids for {rows} rows")]
    LengthMismatch { ids: usize, rows: usize },
    /// Callers key results by this id, so a repeat makes one row unreachable.
    #[error("duplicate source id {id}; ids must be unique per row")]
    Duplicate { id: String },
    #[error("source id blob exceeds u32 offset range")]
    BlobTooLarge,
    /// A library names its rows one way or the other. Storing both shapes
    /// together would mean rendering the numbers as strings, which is the
    /// coercion this type exists to avoid, and would make `7` and `"7"` two
    /// different ids in one column.
    ///
    /// Names the offending row: the check runs over the whole library, so
    /// without it a caller has no way to find the row that broke it.
    #[error(
        "row {row} has a text id ({value:?}) but row {first_numeric_row} has a numeric one; \
         a library names its rows one way or the other"
    )]
    MixedShapes {
        row: usize,
        value: String,
        first_numeric_row: usize,
    },
}

impl SourceIds {
    /// All-numeric. Kept dense rather than rendered into the blob, so the common
    /// JSON path does not pay for one.
    fn numeric(ids: Vec<u64>) -> Result<Self, SourceIdError> {
        let mut seen = std::collections::HashSet::with_capacity(ids.len());
        for id in &ids {
            if !seen.insert(*id) {
                return Err(SourceIdError::Duplicate { id: id.to_string() });
            }
        }
        Ok(Self::Numeric(ids))
    }

    /// Ids the file spelled as text, carried through verbatim.
    fn text<I, S>(ids: I) -> Result<Self, SourceIdError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut blob = String::new();
        let mut offsets = vec![0u32];
        for id in ids {
            blob.push_str(id.as_ref());
            offsets.push(u32::try_from(blob.len()).map_err(|_| SourceIdError::BlobTooLarge)?);
        }
        // Dedup over the finished blob, so a row's id is borrowed rather than
        // cloned into the set.
        let mut seen = std::collections::HashSet::with_capacity(offsets.len() - 1);
        for w in offsets.windows(2) {
            let id = &blob[w[0] as usize..w[1] as usize];
            if !seen.insert(id) {
                return Err(SourceIdError::Duplicate { id: id.to_string() });
            }
        }
        Ok(Self::Text { blob, offsets })
    }

    /// Ids of whichever shape the file used, kept as they came. The only way in;
    /// the length check lives here, once.
    pub fn owned(ids: Vec<OwnedSourceId>, n_rows: usize) -> Result<Self, SourceIdError> {
        if ids.len() != n_rows {
            return Err(SourceIdError::LengthMismatch {
                ids: ids.len(),
                rows: n_rows,
            });
        }
        let first_numeric = ids
            .iter()
            .position(|id| matches!(id, OwnedSourceId::Numeric(_)));
        let first_text = ids
            .iter()
            .position(|id| matches!(id, OwnedSourceId::Text(_)));

        match (first_numeric, first_text) {
            // No text at all, including the empty library.
            (_, None) => Self::numeric(
                ids.iter()
                    .map(|id| match id {
                        OwnedSourceId::Numeric(n) => *n,
                        OwnedSourceId::Text(_) => unreachable!("no row holds text"),
                    })
                    .collect(),
            ),
            (None, Some(_)) => Self::text(ids.iter().map(|id| match id {
                OwnedSourceId::Text(s) => s.as_str(),
                OwnedSourceId::Numeric(_) => unreachable!("no row holds a number"),
            })),
            (Some(first_numeric_row), Some(row)) => Err(SourceIdError::MixedShapes {
                row,
                value: ids[row].to_string(),
                first_numeric_row,
            }),
        }
    }

    /// Self-incremental ids for a format that carries none, so results still
    /// have something to be keyed by. Numerically equal to the row positions
    /// they are minted from, but they are ids from here on: stored in the
    /// column, reported in output, and no longer tied to where the row sits.
    pub fn minted(n_rows: usize) -> Self {
        Self::Numeric((0..n_rows as u64).collect())
    }

    pub fn get(&self, row: usize) -> Option<SourceId<'_>> {
        match self {
            Self::Absent => None,
            Self::Numeric(ids) => ids.get(row).copied().map(SourceId::Numeric),
            Self::Text { blob, offsets } => {
                let start = *offsets.get(row)? as usize;
                let end = *offsets.get(row + 1)? as usize;
                Some(SourceId::Text(&blob[start..end]))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn owned(ids: Vec<OwnedSourceId>) -> Result<SourceIds, SourceIdError> {
        let n = ids.len();
        SourceIds::owned(ids, n)
    }

    /// A caller keys results by this id, so a repeat makes one row unreachable.
    #[test]
    fn duplicate_ids_are_rejected() {
        assert_eq!(
            owned(vec![OwnedSourceId::Numeric(7), OwnedSourceId::Numeric(7)]),
            Err(SourceIdError::Duplicate {
                id: "7".to_string()
            })
        );
        assert_eq!(
            owned(vec!["AAAK2".into(), "AAAK2".into()]),
            Err(SourceIdError::Duplicate {
                id: "AAAK2".to_string()
            })
        );
    }

    #[test]
    fn a_row_count_that_disagrees_with_the_ids_is_rejected() {
        assert_eq!(
            SourceIds::owned(vec![OwnedSourceId::Numeric(7)], 2),
            Err(SourceIdError::LengthMismatch { ids: 1, rows: 2 })
        );
    }

    /// Carafe reads `id` with fastjson and keys results by it, so a numeric id
    /// must stay a JSON number even though the type can also hold a string.
    #[test]
    fn each_id_serializes_as_its_own_shape() {
        assert_eq!(
            serde_json::to_string(&OwnedSourceId::Numeric(7)).unwrap(),
            "7"
        );
        assert_eq!(
            serde_json::to_string(&OwnedSourceId::Text("AAAK2".into())).unwrap(),
            r#""AAAK2""#
        );
    }

    /// A library names its rows or it does not. Storing a mixed set would mean
    /// rendering `Numeric(7)` as `"7"`, which is the coercion this type exists
    /// to avoid -- and would then collide with a genuine name of `"7"`.
    #[test]
    fn mixed_shapes_are_refused_rather_than_coerced() {
        let mixed = vec![
            OwnedSourceId::Numeric(7),
            OwnedSourceId::Text("AAAK2".into()),
        ];
        assert_eq!(
            SourceIds::owned(mixed, 2),
            Err(SourceIdError::MixedShapes {
                row: 1,
                value: "AAAK2".to_string(),
                first_numeric_row: 0,
            }),
            "the error has to name the row that broke it"
        );

        // Either shape on its own is fine, and numeric keeps the dense column.
        assert!(matches!(
            SourceIds::owned(vec![OwnedSourceId::Numeric(7)], 1),
            Ok(SourceIds::Numeric(_))
        ));
        assert!(matches!(
            SourceIds::owned(vec![OwnedSourceId::Text("AAAK2".into())], 1),
            Ok(SourceIds::Text { .. })
        ));
    }

    #[test]
    fn text_ids_round_trip_through_the_blob() {
        let ids = owned(vec!["AAAK2".into(), "PEPTIDEK3".into()]).unwrap();
        assert_eq!(ids.get(0), Some(SourceId::Text("AAAK2")));
        assert_eq!(ids.get(1), Some(SourceId::Text("PEPTIDEK3")));
        assert_eq!(ids.get(2), None);
    }
}
