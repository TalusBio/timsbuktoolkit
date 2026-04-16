use timsseek::IonAnnot;
use timsseek::data_sources::speclib::{PrecursorEntry, ReferenceEG, SerSpeclibElement};
use timsseek::fragment_mass::elution_group_converter::{
    count_carbon_sulphur_in_sequence, supersimpleprediction,
};
use timsseek::isotopes::peptide_isotopes;

use crate::koina::models::{FragmentPrediction, RtPrediction};

// ── Filters ──────────────────────────────────────────────────────────────────

pub struct EntryFilters {
    pub max_ions: usize,
    pub min_mz: f32,
    pub max_mz: f32,
    pub min_ion_mz: f32,
    pub max_ion_mz: f32,
    pub min_ions: usize,
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Strip bracket-enclosed modifications from a sequence.
///
/// "PEPTC[U:4]IDEK" → "PEPTCIDEK"
pub fn strip_mods(seq: &str) -> String {
    let mut out = String::with_capacity(seq.len());
    let mut depth = 0usize;
    for ch in seq.chars() {
        match ch {
            '[' => depth += 1,
            ']' => {
                depth = depth.saturating_sub(1);
            }
            _ if depth == 0 => out.push(ch),
            _ => {}
        }
    }
    out
}

/// Convert our short mod notation to ProForma for rustyms.
/// `C[U:4]` → `C[UNIMOD:4]`, mass-shift mods like `[+57.02]` pass through.
/// Note: rustyms parses `[U:N]` as element Uranium, not UNIMOD — must expand.
fn to_proforma(seq: &str) -> String {
    seq.replace("[U:", "[UNIMOD:")
}

/// Compute the monoisotopic precursor m/z using rustyms.
/// Input should be the modified sequence (mods included in mass).
fn compute_precursor_mz(modified_seq: &str, charge: u8) -> Option<f64> {
    use rustyms::prelude::*;
    let proforma = to_proforma(modified_seq);
    let peptide = Peptidoform::pro_forma(&proforma, None).ok()?;
    let linear = peptide.as_linear()?;
    let formulas = linear.formulas();
    if formulas.is_empty() {
        return None;
    }
    let mass = formulas[0].monoisotopic_mass().value;
    let proton_mass = 1.007_276_f64;
    Some((mass + proton_mass * charge as f64) / charge as f64)
}

/// Count carbon and sulphur from modified sequence (mods affect formula).
fn count_cs_modified(modified_seq: &str) -> Option<(u16, u16)> {
    let proforma = to_proforma(modified_seq);
    count_carbon_sulphur_in_sequence(&proforma).ok()
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Convert Koina predictions + metadata into a [`SerSpeclibElement`].
///
/// Returns `None` when:
/// - Precursor m/z falls outside `[filters.min_mz, filters.max_mz]`.
/// - Fewer than `filters.min_ions` fragment ions survive the ion m/z filter.
/// - Mass computation fails (malformed sequence).
pub fn build_entry(
    sequence: &str,
    charge: u8,
    decoy: bool,
    decoy_group: u32,
    id: u32,
    fragment: &FragmentPrediction,
    rt: &RtPrediction,
    filters: &EntryFilters,
) -> Option<SerSpeclibElement> {
    // 1. Carbon / sulphur count from modified sequence (includes mod contributions).
    let (ncarbon, nsulphur) = count_cs_modified(sequence)?;
    let iso = peptide_isotopes(ncarbon, nsulphur);

    // 2. Precursor m/z from modified sequence (includes mod masses).
    let precursor_mz = compute_precursor_mz(sequence, charge)?;

    // 4. Filter by precursor m/z range.
    if precursor_mz < filters.min_mz as f64 || precursor_mz > filters.max_mz as f64 {
        return None;
    }

    // 5. Ion mobility prediction.
    let mobility = supersimpleprediction(precursor_mz, charge as i32) as f32;

    // 6. Filter fragments: keep those within ion m/z bounds.
    let min_ion = filters.min_ion_mz as f64;
    let max_ion = filters.max_ion_mz as f64;

    // Collect (mz, intensity, annotation_str) triples that pass the m/z filter.
    let mut passing: Vec<(f64, f32, &str)> = fragment
        .mzs
        .iter()
        .zip(fragment.intensities.iter())
        .zip(fragment.annotations.iter())
        .filter(|((mz, _), _)| **mz >= min_ion && **mz <= max_ion)
        .map(|((mz, intensity), ann)| (*mz, *intensity, ann.as_str()))
        .collect();

    // 7. Sort by intensity descending, truncate to max_ions.
    passing.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    passing.truncate(filters.max_ions);

    // 8. Require minimum number of surviving fragments.
    if passing.len() < filters.min_ions {
        return None;
    }

    // 9. Parse ion annotation strings; skip unparseable.
    let mut fragment_mzs: Vec<f64> = Vec::with_capacity(passing.len());
    let mut fragment_intensities: Vec<f32> = Vec::with_capacity(passing.len());
    let mut fragment_labels: Vec<IonAnnot> = Vec::with_capacity(passing.len());

    for (mz, intensity, ann_str) in passing {
        match IonAnnot::try_from(ann_str) {
            Ok(ion) => {
                fragment_mzs.push(mz);
                fragment_intensities.push(intensity);
                fragment_labels.push(ion);
            }
            Err(_) => {
                tracing::debug!("Skipping unparseable ion annotation: {ann_str}");
            }
        }
    }

    // After dropping unparseable annotations the count may fall below min_ions.
    if fragment_labels.len() < filters.min_ions {
        return None;
    }

    // 10. Build precursor labels and intensities from the isotope distribution.
    let precursor_labels: Vec<i8> = vec![0i8, 1i8, 2i8];
    let precursor_intensities: Vec<f32> = vec![iso[0], iso[1], iso[2]];

    // 11. Assemble the element.
    let precursor = PrecursorEntry::new(sequence.to_owned(), charge, decoy, decoy_group);
    let elution_group = ReferenceEG::new(
        id,
        precursor_mz,
        precursor_labels,
        fragment_mzs,
        fragment_labels,
        precursor_intensities,
        fragment_intensities,
        mobility,
        rt.irt,
    );

    Some(SerSpeclibElement::new(precursor, elution_group))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::koina::models::{FragmentPrediction, RtPrediction};

    fn make_filters(min_ions: usize) -> EntryFilters {
        EntryFilters {
            max_ions: 10,
            min_mz: 400.0,
            max_mz: 2000.0,
            min_ion_mz: 150.0,
            max_ion_mz: 2000.0,
            min_ions,
        }
    }

    fn four_ion_fragment() -> FragmentPrediction {
        FragmentPrediction {
            annotations: vec![
                "y3^1".to_string(),
                "y4^1".to_string(),
                "b3^1".to_string(),
                "b4^1".to_string(),
            ],
            mzs: vec![400.0, 500.0, 300.0, 350.0],
            intensities: vec![0.9, 0.8, 0.7, 0.6],
        }
    }

    #[test]
    fn test_strip_mods() {
        assert_eq!(strip_mods("PEPTC[U:4]IDEK"), "PEPTCIDEK");
        assert_eq!(strip_mods("PEPTM[+15.995]IDEK"), "PEPTMIDEK");
        assert_eq!(strip_mods("PEPTIDEK"), "PEPTIDEK");
    }

    #[test]
    fn test_compute_precursor_mz() {
        let mz = compute_precursor_mz("PEPTIDEK", 2).unwrap();
        assert!(mz > 400.0 && mz < 600.0, "Expected reasonable m/z, got {mz}");
    }

    #[test]
    fn test_precursor_mz_includes_mod_mass() {
        // to_proforma converts [U:4] → [UNIMOD:4] before rustyms
        let mz_unmod = compute_precursor_mz("PEPTCIDEK", 2).unwrap();
        let mz_mod = compute_precursor_mz("PEPTC[U:4]IDEK", 2).unwrap();
        let diff = mz_mod - mz_unmod;
        // Cys carbamidomethyl = +57.02 Da, at charge 2 = +28.51 m/z
        assert!(
            (diff - 28.51).abs() < 0.1,
            "Mod should add ~28.51 m/z at z=2, got diff={diff}"
        );
    }


    #[test]
    fn test_build_entry_basic() {
        let fragment = four_ion_fragment();
        let rt = RtPrediction { irt: 30.0 };
        let filters = make_filters(3);

        let result = build_entry(
            "PEPTIDEK",
            2,
            false,
            0,
            42,
            &fragment,
            &rt,
            &filters,
        );

        assert!(result.is_some(), "Expected Some but got None");
    }

    #[test]
    fn test_build_entry_skipped_too_few_ions() {
        // Only 1 ion — min_ions = 3 → should return None.
        let fragment = FragmentPrediction {
            annotations: vec!["y3^1".to_string()],
            mzs: vec![400.0],
            intensities: vec![0.9],
        };
        let rt = RtPrediction { irt: 30.0 };
        let filters = make_filters(3);

        let result = build_entry(
            "PEPTIDEK",
            2,
            false,
            0,
            1,
            &fragment,
            &rt,
            &filters,
        );

        assert!(result.is_none(), "Expected None but got Some");
    }

    #[test]
    fn test_build_entry_precursor_mz_filter() {
        let fragment = four_ion_fragment();
        let rt = RtPrediction { irt: 30.0 };
        // Tight precursor window that excludes PEPTIDEK/2 (~476 Da)
        let filters = EntryFilters {
            max_ions: 10,
            min_mz: 1800.0,
            max_mz: 2000.0,
            min_ion_mz: 150.0,
            max_ion_mz: 2000.0,
            min_ions: 3,
        };

        let result = build_entry(
            "PEPTIDEK",
            2,
            false,
            0,
            2,
            &fragment,
            &rt,
            &filters,
        );

        assert!(result.is_none(), "Expected None due to precursor m/z filter");
    }

    #[test]
    fn test_build_entry_decoy() {
        let fragment = four_ion_fragment();
        let rt = RtPrediction { irt: 30.0 };
        let filters = make_filters(3);

        let result = build_entry(
            "KEDITREP",
            2,
            true,
            99,
            7,
            &fragment,
            &rt,
            &filters,
        );

        // Decoy peptide should still build an entry if the mz/ions pass filters
        // (KEDITREP ~476 should pass default 400–2000 window).
        assert!(result.is_some(), "Expected Some for decoy entry");
    }
}
