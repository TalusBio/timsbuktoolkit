//! One C/S envelope method for the entire library, resolved before extraction.
//!
//! The existing three-bin C/S approximation is retained; this is not a full
//! elemental isotope calculation. Modifications contribute signed C/S deltas.
//! A mass-only annotation cannot establish those deltas. Generated mass shifts
//! reuse their stored parent's envelope, including when counts come from mass.
use std::collections::{
    BTreeMap,
    HashMap,
};

use mzcore::prelude::{
    AmbiguousMolecule,
    AminoAcid,
    Element,
    MolecularFormula,
    Molecule,
};
use mzcore::sequence::SimpleModificationInner;
use serde::Serialize;
use timsquery::IonAnnot;
use timsquery::chemistry::analyte::{
    AnalyteRef,
    FormulaBasis,
    PeptideRef,
};
use timsquery::chemistry::ontologies;
use timsquery::models::target_columns::RowValues;
use timsquery::models::{
    RowIdx,
    TargetColumns,
};
use timsquery::utils::constants::PROTON_MASS;

use super::averagine::isotope_dist_from_mass;
use crate::isotopes::peptide_isotopes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum IsotopeMethod {
    CompositionCs,
    MassEstimatedCs,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum UnavailableReason {
    MissingChemistry,
    IncompleteModifications,
    UnresolvedModification,
    UnresolvedResidues,
    UnsupportedFormula,
    InvalidCounts,
    ModelRange,
}

type Counts = (i64, i64);
type Resolution = Result<Counts, UnavailableReason>;

/// Reported with the scoring plan; cached envelopes never enter feature metadata.
#[derive(Debug, Clone, Serialize)]
pub struct IsotopePlan {
    pub method: IsotopeMethod,
    pub composition_rows: usize,
    pub total_rows: usize,
    pub unavailable: BTreeMap<UnavailableReason, usize>,
    #[serde(skip)]
    envelopes: RowValues<[f32; 3]>,
}

impl IsotopePlan {
    pub(crate) fn resolve(geom: &TargetColumns<IonAnnot>) -> Result<Self, String> {
        let mut resolver = Resolver::default();
        let mut unavailable = BTreeMap::new();
        let mut composition_rows = 0;
        let mut conflict = None;
        let counts = geom.map_rows(|row| {
            let result = resolver.analyte(geom.analyte(row));
            match result {
                Ok(Ok(cs)) => {
                    composition_rows += 1;
                    Some(cs)
                }
                Ok(Err(reason)) => {
                    *unavailable.entry(reason).or_insert(0) += 1;
                    None
                }
                Err(message) => {
                    conflict = Some(format!("entry {:?}: {message}", geom.output_id(row)));
                    None
                }
            }
        });
        if let Some(message) = conflict {
            return Err(message);
        }
        let total_rows = geom.n_rows();
        let method = if composition_rows == total_rows {
            IsotopeMethod::CompositionCs
        } else {
            IsotopeMethod::MassEstimatedCs
        };
        let mut invalid_mass = None;
        let envelopes = geom.map_rows(|row| match method {
            IsotopeMethod::CompositionCs => {
                let (c, s) = counts[row].expect("library plan promised C/S counts");
                peptide_isotopes(c as u16, s as u16)
            }
            IsotopeMethod::MassEstimatedCs => {
                // Stored geometry, never the synthetic variant's shifted mass.
                let z = f64::from(geom.charge(row));
                let envelope = isotope_dist_from_mass(geom.precursor_mz(row) * z - z * PROTON_MASS);
                if !envelope.iter().all(|v| v.is_finite()) {
                    invalid_mass = Some(format!(
                        "entry {}: precursor mass exceeds the C/S isotope model's numerical range",
                        geom.output_id(row)
                    ));
                }
                envelope
            }
        });
        if let Some(message) = invalid_mass {
            return Err(message);
        }
        Ok(Self {
            method,
            composition_rows,
            total_rows,
            unavailable,
            envelopes,
        })
    }

    pub(crate) fn envelope(&self, row: RowIdx) -> &[f32; 3] {
        &self.envelopes[row]
    }
}

impl std::fmt::Display for IsotopePlan {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let method = match self.method {
            IsotopeMethod::CompositionCs => "composition-derived C/S",
            IsotopeMethod::MassEstimatedCs => "mass-estimated C/S (averagine)",
        };
        write!(
            f,
            "Precursor isotopes: {method} for all {} entries ({}/{} with usable composition counts)",
            self.total_rows, self.composition_rows, self.total_rows
        )?;
        for (reason, count) in &self.unavailable {
            let reason = match reason {
                UnavailableReason::MissingChemistry => "missing or unresolved chemistry",
                UnavailableReason::IncompleteModifications => "incomplete modifications",
                UnavailableReason::UnresolvedModification => "unresolved modification composition",
                UnavailableReason::UnresolvedResidues => "unresolved residue C/S counts",
                UnavailableReason::UnsupportedFormula => {
                    "unsupported formula basis or labelled C/S"
                }
                UnavailableReason::InvalidCounts => "invalid C/S counts",
                UnavailableReason::ModelRange => "outside the isotope model's numerical range",
            };
            write!(f, "; {count} {reason}")?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct Resolver {
    modifications: HashMap<String, Resolution>,
}

impl Resolver {
    // Conflicting comparable declarations are invalid input, not absent chemistry.
    fn analyte(&mut self, analyte: AnalyteRef<'_>) -> Result<Resolution, String> {
        let peptide = analyte.peptide.known().map(|p| self.peptide(p));
        let formula = analyte.formula.known().map(|f| {
            if f.basis == FormulaBasis::Unspecified {
                return Err(UnavailableReason::UnsupportedFormula);
            }
            elements_cs(&f.elements).and_then(valid_counts)
        });
        if let (Some(Ok(a)), Some(Ok(b))) = (peptide, formula)
            && analyte
                .formula
                .known()
                .is_some_and(|f| f.basis == FormulaBasis::NeutralMolecule)
            && a != b
        {
            return Err(
                "declared neutral formula disagrees with peptide/modification C/S counts".into(),
            );
        }
        // An explicitly supplied formula can establish composition without sequence.
        // An observed-ion formula already includes adduct atoms; don't add them again.
        if let Some(result) = formula {
            return Ok(result);
        }
        Ok(peptide.unwrap_or(Err(UnavailableReason::MissingChemistry)))
    }

    fn peptide(&mut self, peptide: PeptideRef<'_>) -> Resolution {
        let modifications = peptide
            .modifications
            .known()
            .ok_or(UnavailableReason::IncompleteModifications)?;
        let mut counts = (0i64, 0i64);
        if let Some((c, s)) = super::elution_group_converter::count_cs_fast(peptide.residues) {
            counts = (i64::from(c), i64::from(s));
        } else {
            for residue in peptide.residues.bytes() {
                // X's placeholder zero formula is not evidence of zero C/S.
                if residue == b'X' {
                    return Err(UnavailableReason::UnresolvedResidues);
                }
                let aa = AminoAcid::try_from(residue as char)
                    .map_err(|_| UnavailableReason::UnresolvedResidues)?;
                let formulas = aa.formulas();
                let mut alternatives = formulas.iter().map(formula_cs);
                let cs = alternatives
                    .next()
                    .ok_or(UnavailableReason::UnresolvedResidues)??;
                if !alternatives.all(|v| v == Ok(cs)) {
                    return Err(UnavailableReason::UnresolvedResidues);
                }
                counts.0 += cs.0;
                counts.1 += cs.1;
            }
        }
        for (_, modification) in modifications.iter() {
            let annotation = modification.annotation();
            let delta = if let Some(&delta) = self.modifications.get(annotation) {
                delta
            } else {
                let delta = modification_cs(annotation);
                self.modifications.insert(annotation.to_owned(), delta);
                delta
            };
            let (c, s) = delta?;
            counts.0 += c;
            counts.1 += s;
        }
        valid_counts(counts)
    }
}

fn modification_cs(annotation: &str) -> Resolution {
    let ((parsed, _), _) = SimpleModificationInner::pro_forma(
        annotation,
        &mut Default::default(),
        &mut Default::default(),
        ontologies(),
    )
    .map_err(|_| UnavailableReason::UnresolvedModification)?;
    let modification = parsed
        .defined()
        .ok_or(UnavailableReason::UnresolvedModification)?;
    if matches!(
        modification.as_ref(),
        SimpleModificationInner::Mass(..) | SimpleModificationInner::Info(..)
    ) {
        return Err(UnavailableReason::UnresolvedModification);
    }
    formula_cs(&modification.formula())
}

fn formula_cs(formula: &MolecularFormula) -> Resolution {
    if *formula.additional_mass() != 0.0 {
        return Err(UnavailableReason::UnresolvedModification);
    }
    elements_cs(formula.elements())
}

fn elements_cs(elements: &[(Element, Option<std::num::NonZeroU16>, i32)]) -> Resolution {
    let mut counts = (0, 0);
    for &(element, isotope, count) in elements {
        if matches!(element, Element::C | Element::S) && isotope.is_some() {
            return Err(UnavailableReason::UnsupportedFormula);
        }
        match element {
            Element::C => counts.0 += i64::from(count),
            Element::S => counts.1 += i64::from(count),
            _ => {}
        }
    }
    Ok(counts)
}

fn valid_counts(cs: Counts) -> Resolution {
    if !(0..=i64::from(u16::MAX)).contains(&cs.0) || !(0..=i64::from(u16::MAX)).contains(&cs.1) {
        Err(UnavailableReason::InvalidCounts)
    } else if !peptide_isotopes(cs.0 as u16, cs.1 as u16)
        .iter()
        .all(|v| v.is_finite())
    {
        Err(UnavailableReason::ModelRange)
    } else {
        Ok(cs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_sources::reference_library::{
        ExpectedIntensity,
        ReferenceLibrary,
    };
    use timsquery::chemistry::analyte::{
        Analyte,
        Formula,
        Property,
    };
    use timsquery::models::capabilities::{
        DecoyPolicy,
        TargetCapabilities,
    };
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };
    use timsquery::serde::TargetTable;

    fn library(
        analytes: &[Analyte],
        shipped_decoy: bool,
    ) -> Result<ReferenceLibrary, crate::errors::TargetReadingError> {
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        for (i, analyte) in analytes.iter().enumerate() {
            builder.push_row(Row {
                precursor_mz: 500.0 + i as f64 * 100.0,
                charge: 2,
                frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
                analyte: analyte.as_input(),
                is_decoy: shipped_decoy && i == 1,
                ..Default::default()
            });
        }
        ReferenceLibrary::try_from(TargetTable::Mzpaf {
            geom: builder
                .seal(if shipped_decoy {
                    DecoyPolicy::Never
                } else {
                    DecoyPolicy::Force
                })
                .unwrap(),
            frag_intens: Some(vec![1.0; analytes.len()]),
        })
    }

    fn envelope(lib: &ReferenceLibrary, row: usize, variant: u8) -> Vec<f32> {
        let geom = lib.geometry();
        let row = geom.rows().nth(row).unwrap();
        lib.item_at(geom.flat_for(row, variant))
            .expected_precursor_envelope()
            .iter()
            .map(|(_, i)| *i)
            .collect()
    }

    #[test]
    fn modifications_contribute_signed_cs_and_other_atoms_leave_model_unchanged() {
        // PEPTCIDEK has C43 S1; carbamidomethyl adds C2, sulfur adds S1.
        for (sequence, expected) in [
            ("PEPTCIDEK", (43, 1)),
            ("PEPTC[UNIMOD:4]IDEK", (45, 1)),
            ("PEPTC[Formula:C2H3NO]IDEK", (45, 1)),
            ("PEPTS[UNIMOD:21]IDE", (37, 0)),
            ("PEPTS[MOD:00046]IDE", (37, 0)),
            ("PEPBK", (25, 0)), // D/N alternatives agree on C/S.
            ("PEPDK", (25, 0)),
            ("[UNIMOD:1]-PEPTCIDEK", (45, 1)),
            ("PEPTC[Formula:S-1]IDEK", (43, 0)),
            ("PEPTC[Formula:S]IDEK", (43, 2)),
            ("PEPTC[Formula:C-1]IDEK", (42, 1)),
            ("PEPTC[Formula:O]IDEK", (43, 1)),
        ] {
            let lib = library(&[Analyte::from_sequence(sequence)], false).unwrap();
            assert_eq!(
                lib.scoring_plan().isotopes().method,
                IsotopeMethod::CompositionCs,
                "{sequence}"
            );
            assert_eq!(
                envelope(&lib, 0, 0),
                peptide_isotopes(expected.0, expected.1)
            );
        }
    }

    #[test]
    fn one_mass_only_entry_selects_mass_for_every_row_and_variant() {
        for shipped in [false, true] {
            let lib = library(
                &[
                    Analyte::from_sequence("PEPTC[UNIMOD:4]IDEK"),
                    Analyte::from_sequence("PEPTC[+57.021]IDEK"),
                ],
                shipped,
            )
            .unwrap();
            let plan = lib.scoring_plan().isotopes();
            assert_eq!(plan.method, IsotopeMethod::MassEstimatedCs);
            assert_eq!(plan.composition_rows, 1);
            assert_eq!(
                plan.unavailable[&UnavailableReason::UnresolvedModification],
                1
            );
            for row in 0..2 {
                let expected =
                    isotope_dist_from_mass((500.0 + row as f64 * 100.0) * 2.0 - 2.0 * PROTON_MASS);
                for variant in 0..if shipped { 1 } else { 3 } {
                    assert_eq!(envelope(&lib, row, variant), expected);
                }
            }
        }
    }

    #[test]
    fn formula_only_neutral_and_ion_inputs_need_no_sequence() {
        for basis in [FormulaBasis::NeutralMolecule, FormulaBasis::ObservedIon] {
            let analyte = Analyte::from_formula(&mzcore::molecular_formula!(C 12 H 20 S 2), basis);
            let lib = library(&[analyte], false).unwrap();
            assert_eq!(
                lib.scoring_plan().isotopes().method,
                IsotopeMethod::CompositionCs
            );
            assert_eq!(envelope(&lib, 0, 0), peptide_isotopes(12, 2));
            assert_eq!(lib.scoring_plan().width(), 0);
        }
    }

    #[test]
    fn invalid_unknown_and_labelled_composition_select_mass() {
        for analyte in [
            Analyte::default(),
            Analyte::from_sequence("PEPXIDEK"),
            Analyte::from_sequence("<13C>PEPTIDEK"),
            Analyte::from_sequence("PEPTC[Formula:C-100]IDEK"),
            Analyte::from_sequence("PEPTC[Formula:[13C2]]IDEK"),
            Analyte::from_formula(
                &mzcore::molecular_formula!(C 10 H 20),
                FormulaBasis::Unspecified,
            ),
            Analyte {
                formula: Property::Known(Formula {
                    elements: vec![(Element::C, std::num::NonZeroU16::new(13), 2)],
                    basis: FormulaBasis::NeutralMolecule,
                }),
                ..Default::default()
            },
        ] {
            let lib = library(std::slice::from_ref(&analyte), false).unwrap();
            assert_eq!(
                lib.scoring_plan().isotopes().method,
                IsotopeMethod::MassEstimatedCs,
                "{analyte:?}"
            );
        }
    }

    #[test]
    fn comparable_formula_conflicts_are_rejected_and_agreement_is_not_double_counted() {
        let mut analyte = Analyte::from_sequence("PEPTC[UNIMOD:4]IDEK");
        for (carbon, valid) in [(45, true), (44, false)] {
            analyte.formula = Property::Known(Formula {
                elements: vec![(Element::C, None, carbon), (Element::S, None, 1)],
                basis: FormulaBasis::NeutralMolecule,
            });
            let result = library(std::slice::from_ref(&analyte), false);
            if valid {
                assert_eq!(envelope(&result.unwrap(), 0, 0), peptide_isotopes(45, 1));
            } else {
                assert!(format!("{:?}", result.unwrap_err()).contains("disagrees"));
            }
        }
    }

    #[test]
    fn numerical_model_limits_are_checked() {
        assert_eq!(
            valid_counts((65536, 0)),
            Err(UnavailableReason::InvalidCounts)
        );
        assert_eq!(valid_counts((65535, 0)), Err(UnavailableReason::ModelRange));
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        builder.push_row(Row {
            precursor_mz: 1e9,
            charge: 2,
            frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
            ..Default::default()
        });
        let result = ReferenceLibrary::try_from(TargetTable::Mzpaf {
            geom: builder.seal(DecoyPolicy::Never).unwrap(),
            frag_intens: Some(vec![1.0]),
        });
        assert!(
            matches!(result, Err(crate::errors::TargetReadingError::InvalidLibrary { message }) if message.contains("numerical range"))
        );
    }

    #[test]
    fn diann_and_mzspeclib_modification_annotations_reach_the_same_envelope() {
        let diann = concat!(
            "ModifiedPeptide\tStrippedPeptide\tPrecursorMz\tPrecursorCharge\tTr_recalibrated\tIonMobility\t",
            "ProteinID\tDecoy\tFragmentMz\tFragmentType\tFragmentNumber\tFragmentCharge\tFragmentLossType\tRelativeIntensity\n",
            "PEPTC(UniMod:4)IDEK\tPEPTCIDEK\t500\t2\t1\t1\tP1\t0\t300\ty\t3\t1\tnoloss\t1\n",
        );
        let mzspeclib = concat!(
            "<mzSpecLib>\nMS:1003186|library format version=1.0\nMS:1003188|library name=isotopes\n",
            "<Spectrum=1>\nMS:1003061|library spectrum name=independent label\n",
            "MS:1000744|selected ion m/z=500\nMS:1000041|charge state=2\n",
            "MS:1003059|number of peaks=1\n<Analyte=1>\n",
            "MS:1003270|proforma peptidoform ion notation=PEPTC[U:Carbamidomethyl]IDEK/2\n",
            "<Peaks>\n300\t1\ty3/0.0\n",
        );
        for (suffix, text) in [(".tsv", diann), (".mzspeclib.txt", mzspeclib)] {
            let file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
            std::fs::write(file.path(), text).unwrap();
            let lib = ReferenceLibrary::from_file(file.path(), Default::default())
                .unwrap_or_else(|e| panic!("{suffix}: {e:?}"));
            assert_eq!(
                lib.scoring_plan().isotopes().method,
                IsotopeMethod::CompositionCs
            );
            for query in lib.iter() {
                let values: Vec<_> = query
                    .expected_precursor_envelope()
                    .iter()
                    .map(|(_, i)| *i)
                    .collect();
                assert_eq!(values, peptide_isotopes(45, 1), "{suffix}");
            }
        }
    }
}
