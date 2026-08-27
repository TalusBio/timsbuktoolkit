//! Source ids, as opposed to arena positions.
//!
//! The arena position is self-incremental and feeds decoy grouping and the
//! q-value determinism key, so it stays private to the arena. [`LibraryId`] is
//! what the source file called the row. Separate types so user-supplied values
//! cannot reach the machinery that needs positions.

use serde::{
    Deserialize,
    Serialize,
};
use thiserror::Error;

/// `#[serde(transparent)]` keeps the wire format a bare integer, so Carafe's
/// `"id": 7` is unchanged.
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

/// Which kind of id a library carries is a property of the format, not of a
/// row: JSON targets and the native speclib carry integers, the DIA-NN tabular
/// formats carry strings (`transition_group_id`), Spectronaut/Skyline carry
/// nothing. `Text` arrives with the reader work that needs it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum SourceIds {
    #[default]
    Absent,
    Numeric(Vec<LibraryId>),
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SourceIdError {
    #[error("{ids} source ids for {rows} rows")]
    LengthMismatch { ids: usize, rows: usize },
    /// Callers key results by this id, so a repeat makes one row unreachable.
    #[error("duplicate source id {id}; ids must be unique per row")]
    Duplicate { id: u64 },
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
                return Err(SourceIdError::Duplicate { id: id.get() });
            }
        }
        Ok(Self::Numeric(ids))
    }

    pub fn get(&self, row: usize) -> Option<LibraryId> {
        match self {
            Self::Absent => None,
            Self::Numeric(ids) => ids.get(row).copied(),
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
            Err(SourceIdError::Duplicate { id: 7 })
        );
    }

    #[test]
    fn library_id_serializes_transparently() {
        assert_eq!(serde_json::to_string(&LibraryId::new(7)).unwrap(), "7");
    }
}
