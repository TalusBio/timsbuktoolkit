//! Source ids, as opposed to arena positions.
//!
//! The arena position is self-incremental and feeds decoy grouping and the
//! q-value determinism key, so it stays private to the arena. A source id is
//! what the source file called the row, carried through verbatim.
//!
//! Ids are polymorphic because formats disagree: JSON targets carry integers,
//! DIA-NN carries `transition_group_id` / the `.speclib` entry name, which are
//! strings. Neither is coerced into the other's shape — a numeric id stays a
//! JSON number so Carafe's `"id": 7` is unchanged, and a text id stays a
//! string.

use serde::{
    Deserialize,
    Serialize,
};
use thiserror::Error;

/// `#[serde(transparent)]` keeps the wire format a bare integer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct LibraryId(u64);

impl LibraryId {
    pub fn new(id: u64) -> Self {
        Self(id)
    }

    pub fn get(self) -> u64 {
        self.0
    }
}

/// One row's id, borrowed from the arena that stores it.
///
/// `#[serde(untagged)]` so it serializes as the shape it is — a bare number or
/// a bare string — rather than as a tagged union.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
#[serde(untagged)]
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
}

impl OwnedSourceId {
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

impl Default for OwnedSourceId {
    fn default() -> Self {
        Self::Numeric(0)
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
    Numeric(Vec<LibraryId>),
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
}

impl SourceIds {
    pub fn numeric(ids: Vec<LibraryId>, n_rows: usize) -> Result<Self, SourceIdError> {
        if ids.len() != n_rows {
            return Err(SourceIdError::LengthMismatch {
                ids: ids.len(),
                rows: n_rows,
            });
        }
        let mut seen = std::collections::HashSet::with_capacity(ids.len());
        for id in &ids {
            if !seen.insert(*id) {
                return Err(SourceIdError::Duplicate {
                    id: id.get().to_string(),
                });
            }
        }
        Ok(Self::Numeric(ids))
    }

    /// Ids the file spelled as text, carried through verbatim.
    pub fn text<I, S>(ids: I, n_rows: usize) -> Result<Self, SourceIdError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut blob = String::new();
        let mut offsets = vec![0u32];
        let mut seen = std::collections::HashSet::with_capacity(n_rows);
        for id in ids {
            let id = id.as_ref();
            if !seen.insert(id.to_string()) {
                return Err(SourceIdError::Duplicate { id: id.to_string() });
            }
            blob.push_str(id);
            offsets.push(u32::try_from(blob.len()).map_err(|_| SourceIdError::BlobTooLarge)?);
        }
        if offsets.len() - 1 != n_rows {
            return Err(SourceIdError::LengthMismatch {
                ids: offsets.len() - 1,
                rows: n_rows,
            });
        }
        Ok(Self::Text { blob, offsets })
    }

    /// Ids of whichever shape the file used, kept as they came.
    pub fn owned(ids: Vec<OwnedSourceId>, n_rows: usize) -> Result<Self, SourceIdError> {
        if ids.len() != n_rows {
            return Err(SourceIdError::LengthMismatch {
                ids: ids.len(),
                rows: n_rows,
            });
        }
        // All-numeric stays numeric so the common JSON path keeps its dense
        // integer column instead of paying for a blob.
        if ids.iter().all(|id| matches!(id, OwnedSourceId::Numeric(_))) {
            let nums = ids
                .iter()
                .map(|id| match id {
                    OwnedSourceId::Numeric(n) => LibraryId::new(*n),
                    OwnedSourceId::Text(_) => unreachable!("checked above"),
                })
                .collect();
            return Self::numeric(nums, n_rows);
        }
        Self::text(ids.iter().map(|id| id.to_string()), n_rows)
    }

    /// Self-incremental ids for a format that carries none, so results still
    /// have something to be keyed by. Numerically equal to the row positions
    /// they are minted from, but they are ids from here on: stored in the
    /// column, reported in output, and no longer tied to where the row sits.
    pub fn minted(n_rows: usize) -> Self {
        Self::Numeric((0..n_rows as u64).map(LibraryId::new).collect())
    }

    pub fn get(&self, row: usize) -> Option<SourceId<'_>> {
        match self {
            Self::Absent => None,
            Self::Numeric(ids) => ids.get(row).map(|id| SourceId::Numeric(id.get())),
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

    #[test]
    fn duplicate_ids_are_rejected() {
        let ids = vec![LibraryId::new(7), LibraryId::new(7)];
        assert_eq!(
            SourceIds::numeric(ids, 2),
            Err(SourceIdError::Duplicate {
                id: "7".to_string()
            })
        );
        assert_eq!(
            SourceIds::text(["AAAK2", "AAAK2"], 2),
            Err(SourceIdError::Duplicate {
                id: "AAAK2".to_string()
            })
        );
    }

    /// Carafe reads `id` with fastjson and keys results by it, so a numeric id
    /// must stay a JSON number even though the type can also hold a string.
    #[test]
    fn each_id_serializes_as_its_own_shape() {
        assert_eq!(serde_json::to_string(&SourceId::Numeric(7)).unwrap(), "7");
        assert_eq!(
            serde_json::to_string(&SourceId::Text("AAAK2")).unwrap(),
            r#""AAAK2""#
        );
    }

    #[test]
    fn text_ids_round_trip_through_the_blob() {
        let ids = SourceIds::text(["AAAK2", "PEPTIDEK3"], 2).unwrap();
        assert_eq!(ids.get(0), Some(SourceId::Text("AAAK2")));
        assert_eq!(ids.get(1), Some(SourceId::Text("PEPTIDEK3")));
        assert_eq!(ids.get(2), None);
    }
}
