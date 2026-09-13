//! Library-wide operation selection, shared by feature names and computation.
//!
//! All stored targets and decoys must supply an operation's required facts.
//! Generated mass-shift variants reuse the same residues/modification count;
//! their geometry changes do not establish different chemistry requirements.
//! A partial peptide may supply complete residues, but a partial modification
//! list does not establish its total count. Formula availability is independent.
//!
//! Each enabled operation contributes its derive-generated names and values in
//! one order. Disabled operations contribute neither values nor missingness
//! indicators. Existing observation-derived blocks retain their fixed widths;
//! their all-NaN checks remain data-dependent guards. Sequence operations have
//! no Parquet score columns; their decisions and coverage are file metadata.
use super::blocks::lazy::FragmentIsotopeScores;
use super::blocks::sequence_counts::{
    ModificationCounts,
    ResidueCounts,
};
use super::blocks::{
    NameSink,
    ScoreBlock,
};
use super::results::ScoringFields;
use crate::fragment_mass::isotope_plan::IsotopePlan;
use serde::Serialize;
use std::sync::Arc;
use timsquery::chemistry::analyte::{
    AnalyteRef,
    PeptideRef,
    PropertyRef,
};
use timsquery::ion::{
    IonAnnot,
    IonSeriesOrdinal,
};
use timsquery::models::TargetColumns;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Requirement {
    ResidueSequence,
    ModificationCount,
}

// One registration supplies operation identity, metadata and dispatch. Requirements
// describe input facts; several operations may depend on the same requirement.
macro_rules! operations {
    ($($operation:ident => $block:ty, $label:literal);+ $(;)?) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
        #[serde(rename_all = "snake_case")]
        pub enum Operation { $($operation),+ }
        impl Operation {
            const ALL: &'static [Self] = &[$(Self::$operation),+];
            fn requirement(self) -> Requirement {
                match self { $(Self::$operation => <$block as ScoreBlock>::requirement().expect("operation must declare its requirement")),+ }
            }
            fn names(self) -> Vec<Arc<str>> {
                let mut names = NameSink::new();
                match self { $(Self::$operation => <$block as ScoreBlock>::nonlinear_feature_names(&mut names)),+ }
                names.into_names()
            }
            fn label(self) -> &'static str {
                match self { $(Self::$operation => $label),+ }
            }
            fn project(self, peptide: PeptideRef<'_>, emit: &mut impl FnMut(&[f64])) {
                match self { $(Self::$operation => emit(&<$block>::compute(peptide).nonlinear_feature_array())),+ }
            }
        }
    };
}
operations! {
    ResidueCounts => ResidueCounts, "Residue counts";
    ModificationCounts => ModificationCounts, "Modification counts";
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct Coverage {
    pub known: usize,
    pub recovered: usize,
    pub missing: usize,
    pub not_applicable: usize,
    pub unresolved: usize,
}

impl Coverage {
    fn record<T>(&mut self, value: PropertyRef<'_, T>, accept_recovered: bool) {
        match value {
            PropertyRef::Known(_) => self.known += 1,
            PropertyRef::Missing => self.missing += 1,
            PropertyRef::NotApplicable => self.not_applicable += 1,
            PropertyRef::Unresolved {
                recovered: Some(_), ..
            } if accept_recovered => self.recovered += 1,
            PropertyRef::Unresolved { .. } => self.unresolved += 1,
        }
    }

    fn available(&self, rows: usize) -> bool {
        rows > 0 && self.known + self.recovered == rows
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct OperationDecision {
    pub operation: Operation,
    pub requirement: Requirement,
    pub enabled: bool,
    pub coverage: Coverage,
    pub columns: Vec<Arc<str>>,
}

impl std::fmt::Display for OperationDecision {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let label = self.operation.label();
        let c = &self.coverage;
        write!(
            f,
            "{label}: {} ({} known, {} recovered; {} missing, {} not applicable, {} unresolved)",
            if self.enabled { "enabled" } else { "disabled" },
            c.known,
            c.recovered,
            c.missing,
            c.not_applicable,
            c.unresolved
        )
    }
}

/// Resolved over every retained fragment of every stored target and decoy.
#[derive(Debug, Clone, Serialize)]
pub struct FragmentIsotopeDecision {
    pub enabled: bool,
    pub usable_fragments: usize,
    pub total_fragments: usize,
    pub columns: Vec<Arc<str>>,
}

/// Owned by its reference library; callers cannot install a plan from another library.
#[derive(Debug, Clone, Serialize)]
pub struct ScoringPlan {
    decoys: timsquery::models::capabilities::DecoyResolution,
    isotopes: IsotopePlan,
    fragment_isotopes: FragmentIsotopeDecision,
    #[serde(skip)]
    linear_indices: Vec<usize>,
    rows: usize,
    unmodified_rows: usize,
    operations: Vec<OperationDecision>,
}
impl ScoringPlan {
    pub(crate) fn resolve(geom: &TargetColumns<IonAnnot>) -> Result<Self, String> {
        let isotopes = IsotopePlan::resolve(geom)?;
        let total_fragments = geom.n_fragments();
        let usable_fragments = geom
            .rows()
            .flat_map(|row| geom.frag_labels(row))
            .filter(|label| {
                !matches!(label.series_ordinal(), IonSeriesOrdinal::unknown { .. })
                    && label.get_charge() > 0
                    && label.try_with_offset_neutrons(1).is_ok()
            })
            .count();
        let enabled = total_fragments > 0 && usable_fragments == total_fragments;
        let mut isotope_names = NameSink::new();
        FragmentIsotopeScores::linear_feature_names(&mut isotope_names);
        let isotope_names = isotope_names.into_names();
        let mut linear_names = NameSink::new();
        ScoringFields::linear_feature_names(&mut linear_names);
        let linear_indices = linear_names
            .into_names()
            .iter()
            .enumerate()
            .filter_map(|(i, name)| (enabled || !isotope_names.contains(name)).then_some(i))
            .collect();
        let fragment_isotopes = FragmentIsotopeDecision {
            enabled,
            usable_fragments,
            total_fragments,
            columns: if enabled { isotope_names } else { Vec::new() },
        };
        let mut unmodified_rows = 0;
        let mut residues = Coverage::default();
        let mut modifications = Coverage::default();
        for row in geom.rows() {
            let analyte = geom.analyte(row);
            residues.record(analyte.peptide, true);
            if let Some(peptide) = analyte.peptide.recovered() {
                // A recovered modification list need not be complete.
                if peptide
                    .modifications
                    .known()
                    .is_some_and(|mods| mods.is_empty())
                {
                    unmodified_rows += 1;
                }
                modifications.record(peptide.modifications, false);
            } else {
                modifications.record(analyte.peptide, false);
            }
        }
        let rows = geom.n_rows();
        let operations = Operation::ALL
            .iter()
            .map(|&operation| {
                let requirement = operation.requirement();
                let coverage = match requirement {
                    Requirement::ResidueSequence => residues.clone(),
                    Requirement::ModificationCount => modifications.clone(),
                };
                let enabled = coverage.available(rows);
                OperationDecision {
                    operation,
                    requirement,
                    enabled,
                    coverage,
                    columns: if enabled {
                        operation.names()
                    } else {
                        Vec::new()
                    },
                }
            })
            .collect();
        Ok(Self {
            decoys: *geom.decoy_resolution(),
            isotopes,
            fragment_isotopes,
            linear_indices,
            rows,
            unmodified_rows,
            operations,
        })
    }

    /// Plan-level counts, shared by CLI and viewer reporting.
    pub fn summary(&self) -> String {
        format!(
            "{} library entries; {} unmodified (known empty modification list); {}; fragment isotope scoring: {} ({}/{} usable fragments)",
            self.rows,
            self.unmodified_rows,
            self.isotopes,
            if self.fragment_isotopes.enabled {
                "enabled"
            } else {
                "disabled"
            },
            self.fragment_isotopes.usable_fragments,
            self.fragment_isotopes.total_fragments
        )
    }

    pub fn fragment_isotopes(&self) -> &FragmentIsotopeDecision {
        &self.fragment_isotopes
    }

    /// Positions in the derive-generated ScoringFields linear lane. Names and
    /// values use the same selection; disabled scores never enter a model.
    pub(crate) fn linear_indices(&self) -> &[usize] {
        &self.linear_indices
    }

    pub fn isotopes(&self) -> &IsotopePlan {
        &self.isotopes
    }

    pub fn unmodified_rows(&self) -> usize {
        self.unmodified_rows
    }

    pub fn operations(&self) -> &[OperationDecision] {
        &self.operations
    }

    pub fn enabled(&self, operation: Operation) -> bool {
        self.operations
            .iter()
            .any(|o| o.operation == operation && o.enabled)
    }

    pub fn width(&self) -> usize {
        self.operations.iter().map(|o| o.columns.len()).sum()
    }

    pub fn names(&self) -> impl Iterator<Item = &Arc<str>> {
        self.operations.iter().flat_map(|o| &o.columns)
    }

    /// Disabled operations perform no chemistry access, computation, or projection.
    pub(crate) fn project<'a>(
        &self,
        analyte: impl FnOnce() -> AnalyteRef<'a>,
        mut emit: impl FnMut(&[f64]),
    ) {
        if self.width() == 0 {
            return;
        }
        let peptide = analyte()
            .peptide
            .recovered()
            .expect("plan promised residues");
        for operation in &self.operations {
            if !operation.enabled {
                continue;
            }
            operation.operation.project(peptide, &mut emit);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsquery::chemistry::analyte::{
        Analyte,
        Property,
    };
    use timsquery::models::capabilities::DecoyPolicy;
    use timsquery::models::{
        Row,
        TargetCapabilities,
        TargetColumnsBuilder,
    };

    fn arena(last: &Analyte) -> TargetColumns<IonAnnot> {
        let normal = Analyte::from_sequence("PEPTIDE");
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        for i in 0..100 {
            builder.push_row(Row {
                analyte: if i == 99 { last } else { &normal }.as_input(),
                entry_name: Some("label unrelated to chemistry"),
                ..Default::default()
            });
        }
        builder.seal(DecoyPolicy::Force).unwrap()
    }

    #[test]
    fn complete_modification_counts_do_not_require_a_supported_formula_or_ontology() {
        for sequence in [
            "PEPT[+79.966]IDE",
            "PEPTS[MOD:00046]IDE",
            "PEPTN[Glycan:HexNAc]IDE",
            &"A".repeat(300),
        ] {
            let geom = arena(&Analyte::from_sequence(sequence));
            let plan = ScoringPlan::resolve(&geom).unwrap();
            assert!(plan.enabled(Operation::ResidueCounts), "{sequence}");
            assert!(plan.enabled(Operation::ModificationCounts), "{sequence}");
            assert_eq!(plan.width(), 22);
            assert_eq!(geom.n_rows(), 100);
            assert_eq!(geom.flats().count(), 300);
        }
    }

    #[test]
    fn partial_modifications_disable_only_modification_counts_for_every_variant() {
        for analyte in [
            Analyte::from_sequence_fields("", "PEPTIDE").unwrap(),
            Analyte::from_sequence_fields("PEP[unresolved]TIDE", "PEPTIDE").unwrap(),
        ] {
            let geom = arena(&analyte);
            let plan = ScoringPlan::resolve(&geom).unwrap();
            assert!(plan.enabled(Operation::ResidueCounts));
            assert!(!plan.enabled(Operation::ModificationCounts));
            assert_eq!(plan.width(), 21);
            assert!(!plan.names().any(|n| n.contains("peptide_n_mods")));
            for flat in geom.flats() {
                let row = geom.item_at(flat).row();
                let mut values = Vec::new();
                plan.project(|| geom.analyte(row), |v| values.extend_from_slice(v));
                assert_eq!(values.len(), 21);
                assert_eq!(values[0], 7.0);
                assert_eq!(geom.entry_name(row), Some("label unrelated to chemistry"));
            }
        }
    }

    #[test]
    fn no_sequence_disables_all_sequence_access_and_preserves_coverage_reasons() {
        for peptide in [
            Property::Missing,
            Property::NotApplicable,
            Property::Unresolved {
                recovered: None,
                annotation: "unreadable".into(),
            },
        ] {
            let geom = arena(&Analyte {
                peptide,
                ..Default::default()
            });
            let plan = ScoringPlan::resolve(&geom).unwrap();
            assert_eq!(plan.width(), 0);
            assert_eq!(plan.unmodified_rows, 99);
            let coverage = &plan.operations()[0].coverage;
            assert_eq!(coverage.known, 99);
            assert_eq!(
                coverage.missing + coverage.not_applicable + coverage.unresolved,
                1
            );
            plan.project(
                || panic!("disabled operations must not look up chemistry"),
                |_| panic!("disabled operations must not emit values"),
            );
        }
    }

    #[test]
    fn stored_decoys_participate_in_coverage() {
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        let normal = Analyte::from_sequence("PEPTIDE");
        builder.push_row(Row {
            analyte: normal.as_input(),
            ..Default::default()
        });
        builder.push_row(Row {
            is_decoy: true,
            ..Default::default()
        });
        let geom = builder.seal(DecoyPolicy::Never).unwrap();
        let plan = ScoringPlan::resolve(&geom).unwrap();
        assert_eq!(plan.width(), 0);
        assert_eq!(plan.operations()[0].coverage.missing, 1);
    }

    #[test]
    fn recovered_residues_do_not_make_a_partial_modification_list_complete() {
        let mut peptide = Analyte::from_sequence("PEPTIDE")
            .peptide
            .known()
            .unwrap()
            .clone();
        peptide.modifications = Property::Unresolved {
            recovered: Some(vec![]),
            annotation: "partial list".into(),
        };
        let geom = arena(&Analyte {
            peptide: Property::Unresolved {
                recovered: Some(peptide),
                annotation: "partial chemistry".into(),
            },
            ..Default::default()
        });
        let plan = ScoringPlan::resolve(&geom).unwrap();
        assert_eq!(plan.operations()[0].coverage.recovered, 1);
        assert!(
            plan.operations()[0]
                .to_string()
                .contains("99 known, 1 recovered")
        );
        assert_eq!(plan.unmodified_rows(), 99);
        assert!(plan.summary().contains("99 unmodified"));
        let metadata = serde_json::to_value(&plan).unwrap();
        assert_eq!(metadata["operations"][0]["operation"], "residue_counts");
        assert_eq!(metadata["operations"][0]["requirement"], "residue_sequence");

        assert_eq!(plan.operations()[1].coverage.unresolved, 1);
        assert_eq!(plan.width(), 21);
    }

    #[test]
    #[should_panic(expected = "plan promised residues")]
    fn broken_plan_promise_is_an_invariant_failure() {
        let geom = arena(&Analyte::from_sequence("PEPTIDE"));
        ScoringPlan::resolve(&geom).unwrap().project(
            || AnalyteRef {
                peptide: PropertyRef::Missing,
                formula: PropertyRef::Missing,
            },
            |_| {},
        );
    }
}
