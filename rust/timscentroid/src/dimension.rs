//! Axis-kind carriers for the run index and speclib.
//!
//! `MobilityKind` records the one fact the query engine needs about the mobility
//! axis: whether it is a searchable TIMS 1/K0 axis (`Ook0`) or not. Non-TIMS
//! sources (mzML, no-IM libraries, FAIMS) all collapse to "not searchable" and
//! take the identical mobility-less fallback path.

/// The kind of mobility axis a run index or speclib carries.
///
/// Only `Ook0` (TIMS inverse-reduced mobility) is searchable. `Absent` (no IM
/// axis at all) and `Unsupported` (an axis present but not modeled, e.g. FAIMS
/// compensation voltage) both mean "the stored mobility scalar is a sentinel"
/// and take the mobility-less fallback path.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum MobilityKind {
    /// TIMS 1/K0 — the only searchable kind (current behavior).
    Ook0,
    /// No IM axis (mzML, no-IM library); the stored scalar is a sentinel.
    Absent,
    /// Axis present but not modeled (e.g. FAIMS CV); the String labels the warning.
    Unsupported(String),
}

/// Defaults to `Ook0` so existing serialized indices (which predate this field)
/// deserialize to the current TIMS behavior via `#[serde(default)]`.
impl Default for MobilityKind {
    fn default() -> Self {
        MobilityKind::Ook0
    }
}

impl MobilityKind {
    /// Whether the mobility axis can be used to filter peaks in the query engine.
    pub fn is_filterable(&self) -> bool {
        matches!(self, MobilityKind::Ook0)
    }

    /// Whether the mobility axis can contribute a scoring feature / calibration.
    pub fn is_scoreable(&self) -> bool {
        matches!(self, MobilityKind::Ook0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ook0_is_filterable_and_scoreable() {
        assert!(MobilityKind::Ook0.is_filterable());
        assert!(MobilityKind::Ook0.is_scoreable());
    }

    #[test]
    fn absent_is_neither_filterable_nor_scoreable() {
        assert!(!MobilityKind::Absent.is_filterable());
        assert!(!MobilityKind::Absent.is_scoreable());
    }

    #[test]
    fn unsupported_is_neither_filterable_nor_scoreable() {
        let k = MobilityKind::Unsupported("FAIMS compensation voltage".to_string());
        assert!(!k.is_filterable());
        assert!(!k.is_scoreable());
    }
}
