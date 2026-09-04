//! Per-precursor fields the arena needs but the geometry does not carry.

use crate::ion::IonAnnot;

/// What a text-format reader knows about a precursor beyond its coordinates.
///
/// One type for every reader. The fields are the ones
/// [`TargetTable::mzpaf_with_intensities`](super::TargetTable) consumes, so a
/// new format is wired in by producing this and nothing else.
#[derive(Debug, Clone)]
pub struct PrecursorExtras {
    /// The peptidoform as the file spelled it, modifications included. The
    /// timsseek parse gate reads this one.
    pub modified_peptide: String,
    /// The bare residues, for the composition-isotope path.
    pub stripped_peptide: String,
    pub protein_id: String,
    /// A decoy the library shipped, as opposed to one the arena derives.
    pub is_decoy: bool,
    /// Reference intensities keyed by fragment label. Matched by label rather
    /// than by position, so a reader that emits fragments in a different order
    /// than it lists intensities cannot silently mispair them.
    pub relative_intensities: Vec<(IonAnnot, f32)>,
}
