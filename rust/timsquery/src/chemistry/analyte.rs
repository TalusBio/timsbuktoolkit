//! # Library analyte metadata
//!
//! `Row.analyte` supplies independent optional peptide and molecular-formula facts.
//! `Row.entry_name` is a display label; `Row.id` is the external source key. Neither
//! is interpreted as chemistry. `RowIdx` addresses storage in one owning library;
//! it is not an external ID or a target/decoy relationship.
//!
//! `TargetColumns::analyte(row)` returns borrowed `AnalyteRef` properties:
//!
//! - `Missing`: no fact supplied.
//! - `NotApplicable`: explicitly inapplicable; absence alone does not establish this.
//! - `Known`: complete at that level. Known residues can have unresolved modifications.
//! - `Unresolved`: preserved chemical annotation, optionally with partial recovery.
//!
//! Readers convert explicit sequence fields once. The arena packs residues and
//! located modifications, interns definitions, and remaps them when reader shards
//! merge. Original sequence spelling is discarded where structure is represented;
//! unresolved chemical content remains available through the property. Names and
//! source IDs retain their original text.
//!
//! | Reader | Chemistry source | Name |
//! |---|---|---|
//! | DIA-NN TSV/Parquet | Modified and stripped sequence fields, checked for disagreement | `transition_group_id` / `Precursor.Id`, when supplied |
//! | DIA-NN binary | Format-defined modified-peptide-plus-charge field | Original field, unchanged |
//! | Spectronaut / Skyline | Explicit peptide fields | No separate entry-name field mapped |
//! | mzSpecLib | Already-parsed analyte; peptide or formula | Library spectrum name, when supplied |
//! | Prediction sink | Explicit ProForma, independently of generated label | Same label as prediction-file output |
//! | Target / ElutionGroupInput JSON arrays | Current schemas supply no analyte facts | Source ID remains available; no name inferred |
//!
//! An mzSpecLib spectrum declaring multiple analytes is rejected rather than
//! silently selecting one. Unsupported peptide structures preserve their annotation.
//! Molecular formulas retain their declared basis (or `Unspecified`) and signed
//! electron counts. Peptide and formula facts can coexist: sealing rejects conflicting
//! neutral formulas for unmodified canonical peptides. Scoring finalization also
//! checks comparable neutral-formula C/S counts against modified peptides. Full
//! formula comparison for modified/ambiguous structures and ion-basis formulas
//! remains unsupported; coexistence alone does not certify chemical consistency.
//!
//! Search candidates carry a row and competition metadata, not copied sequences.
//! Rescorers and the dashboard receive the owning `ReferenceLibrary`. Its library-wide
//! plan enables residue and modification-count operations independently. Unknown
//! modification composition does not prevent counting a complete modification set;
//! one row lacking residues disables residue features for all targets and decoys.
//! Disabled operations have no ML projections. Global/labile/ambiguous modifications
//! remain unresolved where their complete set cannot be represented. The isotope
//! model in timsseek includes modification C/S deltas or explicitly based formulas.
//! It selects composition-derived counts only when every stored target and decoy
//! supports them; otherwise all entries use mass-estimated C/S. Both paths retain
//! the same approximate C/S calculator. Synthetic variants reuse parent envelopes.
//! Labelled C/S isotopes and unspecified formula bases are unresolved for this model;
//! other elements do not contribute. This is not a full elemental isotope model.
//!
//! Peptide deduplication compares complete structural keys plus charge and m/z,
//! then chooses the highest score, preferring a target on a tie. Equivalent
//! modification spellings therefore collapse. Incomparable entries remain separate;
//! formula equality and labels do not establish chemical identity or competition.
//!
//! Results format version 5 resolves metadata from the library. `sequence` is nullable
//! and canonically formatted; `entry_name`, `molecular_formula`, and `formula_basis`
//! are nullable columns. Unformattable properties produce null. The viewer's Analyte
//! column displays sequence or formula. `Analyte` serialization preserves property
//! states and supports programmatic metadata round-trips; it does not change the
//! existing geometry-only target-list JSON schemas.
use super::{
    CANONICAL_AA_LETTERS,
    normalize_to_proforma,
    ontologies,
};
use mzcore::prelude::IsAminoAcid;
use serde::{
    Deserialize,
    Serialize,
};
use std::collections::HashMap;
use std::ops::Range;

/// Evidence state, independent of whether an operation is enabled library-wide.
/// `Unresolved.recovered` is partial; `known()` deliberately excludes it.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub enum Property<T> {
    #[default]
    Missing,
    NotApplicable,
    Known(T),
    Unresolved {
        recovered: Option<T>,
        annotation: String,
    },
}
impl<T> Property<T> {
    fn map<U>(self, f: impl FnOnce(T) -> U) -> Property<U> {
        match self {
            Self::Missing => Property::Missing,
            Self::NotApplicable => Property::NotApplicable,
            Self::Known(value) => Property::Known(f(value)),
            Self::Unresolved {
                recovered,
                annotation,
            } => Property::Unresolved {
                recovered: recovered.map(f),
                annotation,
            },
        }
    }

    pub fn known(&self) -> Option<&T> {
        if let Self::Known(v) = self {
            Some(v)
        } else {
            None
        }
    }

    fn view<'a, U>(&'a self, f: impl FnOnce(&'a T) -> U) -> PropertyRef<'a, U> {
        match self {
            Self::Missing => PropertyRef::Missing,
            Self::NotApplicable => PropertyRef::NotApplicable,
            Self::Known(v) => PropertyRef::Known(f(v)),
            Self::Unresolved {
                recovered,
                annotation,
            } => PropertyRef::Unresolved {
                recovered: recovered.as_ref().map(f),
                annotation,
            },
        }
    }
}
#[derive(Debug, Clone, Copy, Default, Serialize)]
pub enum PropertyRef<'a, T> {
    #[default]
    Missing,
    NotApplicable,
    Known(T),
    Unresolved {
        recovered: Option<T>,
        annotation: &'a str,
    },
}
impl<T> PropertyRef<'_, T> {
    pub fn known(self) -> Option<T> {
        if let Self::Known(v) = self {
            Some(v)
        } else {
            None
        }
    }

    /// Return the value from `Known`, or the partial value from `Unresolved`.
    /// Unlike `known()`, this does not promise completeness. For example, a
    /// recovered peptide may have known residues but an incomplete modification list.
    ///
    /// ```
    /// use timsquery::chemistry::analyte::PropertyRef;
    /// let partial = PropertyRef::Unresolved {
    ///     recovered: Some("PEPTIDE"), annotation: "unresolved chemistry",
    /// };
    /// assert_eq!(partial.known(), None);
    /// assert_eq!(partial.recovered(), Some("PEPTIDE"));
    /// ```
    pub fn recovered(self) -> Option<T> {
        match self {
            Self::Known(v) => Some(v),
            Self::Unresolved { recovered, .. } => recovered,
            _ => None,
        }
    }

    fn owned<U>(self, f: impl FnOnce(T) -> U) -> Property<U> {
        match self {
            Self::Missing => Property::Missing,
            Self::NotApplicable => Property::NotApplicable,
            Self::Known(v) => Property::Known(f(v)),
            Self::Unresolved {
                recovered,
                annotation,
            } => Property::Unresolved {
                recovered: recovered.map(f),
                annotation: annotation.into(),
            },
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum ModificationSite {
    NTerm,
    Residue(u32),
    CTerm,
}
/// A located occurrence; registry indices never cross the builder boundary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Modification {
    pub site: ModificationSite,
    pub definition: ModificationDefinition,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ModificationKind {
    Unimod(u32),
    Mass(f64),
    Other,
}
/// Canonical or unresolved modification annotation, independent of its site.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ModificationDefinition {
    annotation: String,
    kind: ModificationKind,
}
impl ModificationDefinition {
    pub fn unimod(id: u32) -> Self {
        Self {
            annotation: format!("UNIMOD:{id}"),
            kind: ModificationKind::Unimod(id),
        }
    }

    pub fn mass(delta: f64) -> Self {
        Self {
            annotation: format!("{delta:+}"),
            kind: ModificationKind::Mass(delta),
        }
    }

    pub fn other(annotation: String) -> Self {
        Self {
            annotation,
            kind: ModificationKind::Other,
        }
    }

    pub fn annotation(&self) -> &str {
        &self.annotation
    }

    fn intern_key(&self) -> String {
        format!("{:?}:{}", self.kind, self.annotation)
    }

    pub fn kind(&self) -> &ModificationKind {
        &self.kind
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Peptide {
    pub residues: String,
    pub modifications: Property<Vec<Modification>>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum FormulaBasis {
    NeutralMolecule,
    ObservedIon,
    #[default]
    Unspecified,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Formula {
    pub elements: Vec<(
        mzcore::chemistry::Element,
        Option<std::num::NonZeroU16>,
        i32,
    )>,
    pub basis: FormulaBasis,
}
impl Formula {
    pub fn notation(&self) -> Option<String> {
        mzcore::chemistry::MolecularFormula::new(&self.elements, &[]).map(|f| f.hill_notation())
    }
}
impl FormulaBasis {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NeutralMolecule => "neutral_molecule",
            Self::ObservedIon => "observed_ion",
            Self::Unspecified => "unspecified",
        }
    }
}

/// Owned reader staging and serialization. No identity or charge is inferred here.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Analyte {
    pub peptide: Property<Peptide>,
    pub formula: Property<Formula>,
}
/// Independent optional facts supplied to the common row builder.
#[derive(Debug, Clone, Copy, Default)]
pub struct AnalyteInput<'a> {
    pub peptide: PropertyRef<'a, PeptideInput<'a>>,
    pub formula: PropertyRef<'a, &'a Formula>,
}
#[derive(Debug, Clone, Copy)]
pub struct PeptideInput<'a> {
    pub residues: &'a str,
    pub modifications: PropertyRef<'a, &'a [Modification]>,
}
/// Borrowed facts from one row of the owning arena.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct AnalyteRef<'a> {
    pub peptide: PropertyRef<'a, PeptideRef<'a>>,
    pub formula: PropertyRef<'a, &'a Formula>,
}
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PeptideRef<'a> {
    pub residues: &'a str,
    pub modifications: PropertyRef<'a, ModificationListRef<'a>>,
}
#[derive(Debug, Clone, Copy)]
pub struct ModificationListRef<'a> {
    entries: &'a [(ModificationSite, usize)],
    registry: &'a [ModificationDefinition],
}
impl Serialize for ModificationListRef<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.collect_seq(self.iter().map(|(site, definition)| Modification {
            site,
            definition: definition.clone(),
        }))
    }
}
impl<'a> ModificationListRef<'a> {
    pub fn iter(self) -> impl Iterator<Item = (ModificationSite, &'a ModificationDefinition)> {
        self.entries
            .iter()
            .map(|&(site, idx)| (site, &self.registry[idx]))
    }

    pub fn len(self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(self) -> bool {
        self.entries.is_empty()
    }
}
#[derive(Debug, thiserror::Error)]
pub enum SequenceFormatError {
    #[error("modifications are missing or unresolved")]
    IncompleteModifications,
    #[error(transparent)]
    Write(#[from] std::fmt::Error),
}

impl PeptideRef<'_> {
    pub fn write_sequence(self, out: &mut impl std::fmt::Write) -> Result<(), SequenceFormatError> {
        let mods = self
            .modifications
            .known()
            .ok_or(SequenceFormatError::IncompleteModifications)?;
        for (site, m) in mods.iter() {
            if site == ModificationSite::NTerm {
                write!(out, "[{}]", m.annotation)?;
            }
        }
        if mods.iter().any(|(site, _)| site == ModificationSite::NTerm) {
            out.write_char('-')?;
        }
        for (i, residue) in self.residues.chars().enumerate() {
            out.write_char(residue)?;
            for (site, m) in mods.iter() {
                if site == ModificationSite::Residue(i as u32) {
                    write!(out, "[{}]", m.annotation)?;
                }
            }
        }
        if mods.iter().any(|(site, _)| site == ModificationSite::CTerm) {
            out.write_char('-')?;
        }
        for (site, m) in mods.iter() {
            if site == ModificationSite::CTerm {
                write!(out, "[{}]", m.annotation)?;
            }
        }
        Ok(())
    }

    pub fn sequence(self) -> Option<String> {
        let mut s = String::new();
        self.write_sequence(&mut s).ok()?;
        Some(s)
    }
}
impl Analyte {
    pub fn as_input(&self) -> AnalyteInput<'_> {
        AnalyteInput {
            peptide: self.peptide.view(|p| PeptideInput {
                residues: &p.residues,
                modifications: p.modifications.view(|m| m.as_slice()),
            }),
            formula: self.formula.view(|f| f),
        }
    }

    /// Check independent stripped/modified fields at the reader boundary.
    pub fn from_sequence_fields(sequence: &str, stripped: &str) -> Result<Self, String> {
        let mut analyte = Self::from_sequence(sequence);
        if let Some(peptide) = analyte.peptide.known()
            && !stripped.is_empty()
            && peptide.residues != stripped
        {
            return Err(format!(
                "modified and stripped sequence disagree: {sequence:?}, {stripped:?}"
            ));
        }
        if !stripped.is_empty() && matches!(analyte.peptide, Property::Missing) {
            analyte.peptide = Property::Known(Peptide {
                residues: stripped.into(),
                modifications: Property::Missing,
            });
        } else if !stripped.is_empty() && matches!(analyte.peptide, Property::Unresolved { .. }) {
            analyte.peptide = Property::Known(Peptide {
                residues: stripped.into(),
                modifications: Property::Unresolved {
                    recovered: None,
                    annotation: sequence.into(),
                },
            });
        }
        Ok(analyte)
    }

    pub fn from_sequence(sequence: &str) -> Self {
        if sequence.is_empty() {
            return Self::default();
        }
        let normalized = normalize_to_proforma(sequence);
        if let Some(peptide) = parse_simple(&normalized) {
            return Self {
                peptide: Property::Known(peptide),
                ..Self::default()
            };
        }
        match mzcore::sequence::Peptidoform::pro_forma(&normalized, ontologies()) {
            Ok((p, _)) => Self::from_peptidoform(&p),
            Err(_) => Self {
                peptide: Property::Unresolved {
                    recovered: None,
                    annotation: sequence.into(),
                },
                ..Self::default()
            },
        }
    }

    pub fn from_peptidoform(p: &mzcore::sequence::Peptidoform<mzcore::sequence::Linked>) -> Self {
        if p.sequence()
            .iter()
            .any(|e| e.aminoacid.aminoacid().one_letter_code().is_none())
        {
            return Self {
                peptide: Property::Unresolved {
                    recovered: None,
                    annotation: p.to_string(),
                },
                ..Default::default()
            };
        }
        let Some(p) = p.as_simple_linear().filter(|p| p.is_semi_ambiguous()) else {
            let residues = p
                .sequence()
                .iter()
                .map(|e| e.aminoacid.aminoacid().one_letter_code().unwrap_or('X'))
                .collect();
            return Self {
                peptide: Property::Known(Peptide {
                    residues,
                    modifications: Property::Unresolved {
                        recovered: None,
                        annotation: p.to_string(),
                    },
                }),
                ..Self::default()
            };
        };
        let residues: String = p
            .sequence()
            .iter()
            .map(|e| e.aminoacid.aminoacid().one_letter_code().unwrap_or('X'))
            .collect();
        let mut mods = Vec::new();
        for (i, e) in p.sequence().iter().enumerate() {
            for m in &e.modifications {
                mods.push(convert_mod(ModificationSite::Residue(i as u32), m));
            }
        }
        for m in p.get_n_term() {
            mods.push(convert_mod(ModificationSite::NTerm, m));
        }
        for m in p.get_c_term() {
            mods.push(convert_mod(ModificationSite::CTerm, m));
        }
        mods.sort_by(|a, b| {
            (a.site, &a.definition.annotation).cmp(&(b.site, &b.definition.annotation))
        });
        Self {
            peptide: Property::Known(Peptide {
                residues,
                modifications: Property::Known(mods),
            }),
            ..Self::default()
        }
    }

    pub fn from_formula(
        formula: &mzcore::chemistry::MolecularFormula,
        basis: FormulaBasis,
    ) -> Self {
        let recovered = Formula {
            elements: formula.elements().to_vec(),
            basis,
        };
        Self {
            peptide: Property::Missing,
            formula: if *formula.additional_mass() == 0.0 {
                Property::Known(recovered)
            } else {
                Property::Unresolved {
                    recovered: Some(recovered),
                    annotation: formula.hill_notation(),
                }
            },
        }
    }
}
fn convert_mod(site: ModificationSite, m: &mzcore::sequence::Modification) -> Modification {
    use mzcore::sequence::{
        Modification as M,
        SimpleModificationInner as S,
    };
    let definition = match m {
        M::Simple(s) => match s.as_ref() {
            S::Database { id, .. } if id.ontology == mzcore::ontology::Ontology::Unimod => {
                match id.id() {
                    mzcv::AccessionCode::Numeric(n) => ModificationDefinition::unimod(n),
                    _ => ModificationDefinition::other(m.to_string()),
                }
            }
            S::Mass(_, mass, _) => ModificationDefinition::mass(mass.value),
            _ => ModificationDefinition::other(m.to_string()),
        },
        _ => ModificationDefinition::other(m.to_string()),
    };
    Modification { site, definition }
}

#[derive(Debug, Clone)]
struct StoredPeptide {
    residues: Range<usize>,
    modifications: Property<Range<usize>>,
}
#[derive(Debug, Clone)]
struct StoredAnalyte {
    peptide: Property<StoredPeptide>,
    formula: Property<Formula>,
}
/// Packed row metadata, with one shared registry across the arena.
#[derive(Debug, Clone, Default)]
pub(crate) struct AnalyteColumns {
    rows: Vec<StoredAnalyte>,
    residues: String,
    modifications: Vec<(ModificationSite, usize)>,
    registry: Vec<ModificationDefinition>,
    intern: HashMap<String, usize>,
}
impl AnalyteColumns {
    pub fn push(&mut self, input: AnalyteInput<'_>) {
        let peptide = input.peptide.owned(|p| {
            let start = self.residues.len();
            self.residues.push_str(p.residues);
            let modifications = p.modifications.owned(|mods| {
                let start = self.modifications.len();
                let mut mods: Vec<_> = mods.iter().collect();
                mods.sort_by(|a, b| {
                    (a.site, &a.definition.annotation).cmp(&(b.site, &b.definition.annotation))
                });
                for m in mods {
                    // Site belongs to the occurrence, not its interned definition.
                    let key = m.definition.intern_key();
                    let idx = *self.intern.entry(key).or_insert_with(|| {
                        let idx = self.registry.len();
                        self.registry.push(m.definition.clone());
                        idx
                    });
                    self.modifications.push((m.site, idx));
                }
                start..self.modifications.len()
            });
            StoredPeptide {
                residues: start..self.residues.len(),
                modifications,
            }
        });
        self.rows.push(StoredAnalyte {
            peptide,
            formula: input.formula.owned(Clone::clone),
        });
    }

    pub fn get(&self, row: usize) -> AnalyteRef<'_> {
        let r = &self.rows[row];
        AnalyteRef {
            peptide: r.peptide.view(|p| PeptideRef {
                residues: &self.residues[p.residues.clone()],
                modifications: p.modifications.view(|range| ModificationListRef {
                    entries: &self.modifications[range.clone()],
                    registry: &self.registry,
                }),
            }),
            formula: r.formula.view(|f| f),
        }
    }

    pub fn append(&mut self, other: Self) {
        let residue_base = self.residues.len();
        let modification_base = self.modifications.len();
        let remap: Vec<_> = other
            .registry
            .into_iter()
            .map(|definition| {
                *self
                    .intern
                    .entry(definition.intern_key())
                    .or_insert_with(|| {
                        let index = self.registry.len();
                        self.registry.push(definition);
                        index
                    })
            })
            .collect();
        self.residues.push_str(&other.residues);
        self.modifications.extend(
            other
                .modifications
                .into_iter()
                .map(|(site, index)| (site, remap[index])),
        );
        self.rows
            .extend(other.rows.into_iter().map(|row| StoredAnalyte {
                peptide: row.peptide.map(|peptide| StoredPeptide {
                    residues: peptide.residues.start + residue_base
                        ..peptide.residues.end + residue_base,
                    modifications: peptide.modifications.map(|range| {
                        range.start + modification_base..range.end + modification_base
                    }),
                }),
                formula: row.formula,
            }));
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn finish(&mut self) {
        self.intern = HashMap::new();
        self.rows.shrink_to_fit();
        self.residues.shrink_to_fit();
        self.modifications.shrink_to_fit();
        self.registry.shrink_to_fit();
    }

    /// Check stored sites, formulas, and comparable peptide/formula declarations.
    /// Public builders run this when sealing:
    ///
    /// ```
    /// use timsquery::chemistry::analyte::{Analyte, FormulaBasis};
    /// use timsquery::models::{Row, TargetColumnsBuilder};
    /// use timsquery::models::capabilities::DecoyPolicy;
    ///
    /// let mut analyte = Analyte::from_sequence("A");
    /// for (formula, valid) in [
    ///     (mzcore::molecular_formula!(C 3 H 7 N 1 O 2), true),
    ///     (mzcore::molecular_formula!(C 2 H 6 O 1), false),
    /// ] {
    ///     analyte.formula = Analyte::from_formula(&formula, FormulaBasis::NeutralMolecule).formula;
    ///     let mut builder = TargetColumnsBuilder::<timsquery::ion::IonAnnot>::with_capabilities(
    ///         timsquery::models::TargetCapabilities::default_diann());
    ///     builder.push_row(Row { analyte: analyte.as_input(), ..Default::default() });
    ///     assert_eq!(builder.seal(DecoyPolicy::Never).is_ok(), valid);
    /// }
    /// ```
    pub fn validate(&self) -> Result<(), String> {
        for i in 0..self.rows.len() {
            let analyte = self.get(i);
            if let Some(formula) = analyte.formula.recovered() {
                let mut species = std::collections::HashSet::new();
                for &(element, isotope, count) in &formula.elements {
                    if (element != mzcore::chemistry::Element::Electron && count < 0)
                        || !element.is_valid(isotope)
                        || !species.insert((element, isotope))
                    {
                        return Err(format!(
                            "row {i}: molecular formula has negative atom counts, invalid isotopes or repeated element species"
                        ));
                    }
                }
            }

            if let (Some(peptide), Some(formula)) =
                (analyte.peptide.known(), analyte.formula.known())
                && formula.basis == FormulaBasis::NeutralMolecule
                && let Some(expected) = unmodified_formula(peptide)
            {
                let supplied = mzcore::chemistry::MolecularFormula::new(&formula.elements, &[]);
                if supplied
                    .as_ref()
                    .is_some_and(|supplied| supplied != &expected)
                {
                    return Err(format!(
                        "row {i}: declared neutral formula disagrees with unmodified peptide"
                    ));
                }
            }

            if let Some(p) = self.get(i).peptide.recovered() {
                if p.residues.is_empty() || !p.residues.bytes().all(|c| c.is_ascii_uppercase()) {
                    return Err(format!("row {i}: residues must be uppercase ASCII"));
                }
                if let Some(mods) = p.modifications.recovered() {
                    for (site, _) in mods.iter() {
                        if let ModificationSite::Residue(pos) = site
                            && pos as usize >= p.residues.len()
                        {
                            return Err(format!(
                                "row {i}: modification site {pos} exceeds residue count"
                            ));
                        }
                    }
                }
            }
        }
        for definition in &self.registry {
            let canonical = match definition.kind {
                ModificationKind::Unimod(id) => ModificationDefinition::unimod(id).annotation,
                ModificationKind::Mass(m) => ModificationDefinition::mass(m).annotation,
                ModificationKind::Other => definition.annotation.clone(),
            };
            if canonical != definition.annotation
                || definition.annotation.is_empty()
                || matches!(definition.kind, ModificationKind::Mass(m) if !m.is_finite())
            {
                return Err("empty modification annotation or non-finite mass".into());
            }
        }
        Ok(())
    }
}
impl AnalyteRef<'_> {
    pub fn to_owned(self) -> Analyte {
        Analyte {
            peptide: self.peptide.owned(|p| Peptide {
                residues: p.residues.into(),
                modifications: p.modifications.owned(|m| {
                    m.iter()
                        .map(|(site, d)| Modification {
                            site,
                            definition: d.clone(),
                        })
                        .collect()
                }),
            }),
            formula: self.formula.owned(Clone::clone),
        }
    }
}

/// Common explicit sequence grammar. Avoid ontology lookup for already numeric
/// accessions/mass deltas; unfamiliar syntax goes through mzcore intact.
///
/// ```
/// use timsquery::chemistry::analyte::Analyte;
///
/// let numeric = Analyte::from_sequence("PEPM[UNIMOD:35]AS[MOD:00046]");
/// let named = Analyte::from_sequence("PEPM[U:Oxidation]AS[M:O-phospho-L-serine]");
/// assert_eq!(numeric, named);
/// let peptide = named.peptide.known().unwrap();
/// assert_eq!(peptide.residues, "PEPMAS");
/// assert_eq!(peptide.modifications.known().unwrap().len(), 2);
/// ```
fn parse_simple(s: &str) -> Option<Peptide> {
    fn definition(s: &str) -> Option<ModificationDefinition> {
        if let Some(id) = s.strip_prefix("UNIMOD:") {
            return id
                .trim()
                .parse::<u16>()
                .ok()
                .map(|id| ModificationDefinition::unimod(id.into()));
        }
        if s.starts_with(['+', '-']) {
            return s
                .parse::<f64>()
                .ok()
                .filter(|m| m.is_finite())
                .map(ModificationDefinition::mass);
        }
        None
    }
    let bytes = s.as_bytes();
    let mut residues = String::new();
    let mut modifications = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i].is_ascii_uppercase() {
            residues.push(bytes[i] as char);
            i += 1;
            continue;
        }
        let terminal = bytes[i] == b'-';
        if terminal {
            i += 1;
        }
        if bytes.get(i) != Some(&b'[') {
            return None;
        }
        let close = i + 1 + bytes[i + 1..].iter().position(|&b| b == b']')?;
        let definition = definition(s[i + 1..close].trim())?;
        let site = if residues.is_empty() {
            if terminal || bytes.get(close + 1) != Some(&b'-') {
                return None;
            }
            ModificationSite::NTerm
        } else if terminal {
            if close + 1 != bytes.len() {
                return None;
            }
            ModificationSite::CTerm
        } else {
            ModificationSite::Residue((residues.len() - 1) as u32)
        };
        i = close + 1;
        if site == ModificationSite::NTerm {
            i += 1;
        }
        modifications.push(Modification { site, definition });
    }
    if residues.is_empty() {
        return None;
    }
    Some(Peptide {
        residues,
        modifications: Property::Known(modifications),
    })
}

/// Compare only declarations whose relationship this slice establishes.
/// Modified/ambiguous peptides and ion-basis formulas remain independent;
/// modification-aware composition is a separate operation.
fn unmodified_formula(peptide: PeptideRef<'_>) -> Option<mzcore::chemistry::MolecularFormula> {
    use mzcore::chemistry::AmbiguousMolecule;
    if !peptide.modifications.known()?.is_empty() {
        return None;
    }
    let mut formula = mzcore::molecular_formula!(H 2 O 1);
    for residue in peptide.residues.bytes() {
        if !CANONICAL_AA_LETTERS.contains(&residue) {
            return None;
        }
        let aa = mzcore::sequence::AminoAcid::try_from(residue).ok()?;
        let formulas = aa.formulas();
        let part = formulas.iter().next()?;
        formula.ref_mut_checked_add(part)?;
    }
    Some(formula)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::capabilities::DecoyPolicy;
    use crate::models::{
        Query,
        Row,
        TargetCapabilities,
        TargetColumnsBuilder,
    };
    use crate::traits::Target;

    #[test]
    fn geometry_and_formula_rows_need_no_peptide() {
        use mzcore::chemistry::Element::{
            C,
            H,
            O,
        };
        let formula = Formula {
            elements: vec![(H, None, 6), (C, None, 2), (O, None, 1)],
            basis: FormulaBasis::NeutralMolecule,
        };
        let analyte = Analyte {
            formula: Property::Known(formula),
            ..Default::default()
        };
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        let frags = [(std::sync::Arc::<str>::from("fragment-label"), 42.0)];
        for input in [AnalyteInput::default(), analyte.as_input()] {
            builder.push_row(Row {
                precursor_mz: 47.0,
                charge: 1,
                frags: &frags,
                analyte: input,
                entry_name: Some("label—not a sequence/2"),
                ..Default::default()
            });
        }
        let geom = builder.seal(DecoyPolicy::Never).unwrap();
        for row in geom.rows() {
            assert!(matches!(geom.analyte(row).peptide, PropertyRef::Missing));
            assert_eq!(geom.entry_name(row), Some("label—not a sequence/2"));
            let query = Query::new(&geom, geom.flat_for(row, 0));
            assert_eq!(query.mono_precursor_mz(), 47.0);
            assert_eq!(query.iter_fragments_refs().count(), 1);
        }
        let row = geom.rows().nth(1).unwrap();
        let json = serde_json::to_string(&geom.analyte(row)).unwrap();
        let restored: Analyte = serde_json::from_str(&json).unwrap();
        assert_eq!(restored, analyte);
        assert_eq!(
            restored.formula.known().unwrap().notation().as_deref(),
            Some("C2H6O1")
        );
    }

    #[test]
    fn definitions_are_canonical_and_remapped_when_shards_merge() {
        let mut left = AnalyteColumns::default();
        let mut right = AnalyteColumns::default();
        let first = Analyte::from_sequence("AC(UniMod:4)M(UniMod:35)");
        let second = Analyte::from_sequence("M[Oxidation]C[Carbamidomethyl]");
        left.push(first.as_input());
        right.push(second.as_input());
        left.append(right);
        left.validate().unwrap();
        assert_eq!(left.registry.len(), 2);
        assert_eq!(left.get(0).to_owned(), first);
        assert_eq!(left.get(1).to_owned(), second);
        assert_eq!(
            left.get(1).peptide.known().unwrap().sequence().as_deref(),
            Some("M[UNIMOD:35]C[UNIMOD:4]")
        );
        let json = serde_json::to_string(&left.get(1)).unwrap();
        let restored: Analyte = serde_json::from_str(&json).unwrap();
        assert_eq!(
            restored, second,
            "serialization must not include other rows' registry entries"
        );
    }

    #[test]
    fn terminal_and_mass_modifications_survive_formatting() {
        for sequence in [
            "[UNIMOD:1]-AC[+57.021464]D-[UNIMOD:2]",
            "M[Oxidation]PEPTIDE",
        ] {
            let input = Analyte::from_sequence(sequence);
            let mut columns = AnalyteColumns::default();
            columns.push(input.as_input());
            columns.validate().unwrap();
            let output = columns.get(0).peptide.known().unwrap().sequence().unwrap();
            assert_eq!(Analyte::from_sequence(&output), input);
        }
    }

    #[test]
    fn proforma_named_ontology_prefixes_are_preserved_and_resolved() {
        // ProForma 2.1 §7.8: one-letter prefixes qualify modification names.
        let sequence = "PEPM[U:Oxidation]AS[M:O-phospho-L-serine]";
        assert_eq!(normalize_to_proforma(sequence), sequence);
        let analyte = Analyte::from_sequence(sequence);
        let peptide = analyte.peptide.known().expect("valid named modifications");
        assert_eq!(peptide.residues, "PEPMAS");
        assert_eq!(peptide.modifications.known().unwrap().len(), 2);
        assert_eq!(
            analyte,
            Analyte::from_sequence("PEPM[UNIMOD:35]AS[MOD:00046]")
        );
    }

    #[test]
    fn stripped_residues_survive_missing_modified_sequence() {
        let analyte = Analyte::from_sequence_fields("", "PEPTIDE").unwrap();
        let peptide = analyte.peptide.known().unwrap();
        assert_eq!(peptide.residues, "PEPTIDE");
        assert!(matches!(peptide.modifications, Property::Missing));
        let mut columns = AnalyteColumns::default();
        columns.push(analyte.as_input());
        columns.validate().unwrap();
        assert!(columns.get(0).peptide.known().unwrap().sequence().is_none());
        assert_eq!(columns.get(0).to_owned(), analyte);
        assert!(matches!(
            Analyte::from_sequence_fields("", "").unwrap().peptide,
            Property::Missing
        ));
    }

    #[test]
    fn unresolved_chemistry_is_preserved_at_the_narrowest_known_level() {
        let analyte = Analyte::from_sequence_fields("PEP[unknown-label]TIDE", "PEPTIDE").unwrap();
        let peptide = analyte.peptide.known().unwrap();
        assert_eq!(peptide.residues, "PEPTIDE");
        assert!(
            matches!(&peptide.modifications, Property::Unresolved { annotation, .. } if annotation == "PEP[unknown-label]TIDE")
        );
        let mut columns = AnalyteColumns::default();
        columns.push(analyte.as_input());
        assert!(columns.get(0).peptide.known().unwrap().sequence().is_none());
        assert_eq!(columns.get(0).to_owned(), analyte);
        assert!(Analyte::from_sequence_fields("PEPTIDE", "DIFFERENT").is_err());
    }

    #[test]
    fn unsupported_global_chemistry_is_not_silently_dropped() {
        let analyte = Analyte::from_sequence("<13C>PEPTIDE");
        let mut columns = AnalyteColumns::default();
        columns.push(analyte.as_input());
        columns.validate().unwrap();
        let peptide = columns.get(0).peptide.known().expect("residues retained");
        assert_eq!(peptide.residues, "PEPTIDE");
        assert!(
            matches!(peptide.modifications, PropertyRef::Unresolved { annotation, .. } if annotation.contains("13C"))
        );
        assert!(peptide.sequence().is_none());
    }

    #[test]
    fn explicit_inapplicability_and_unresolved_annotations_round_trip() {
        let analyte = Analyte {
            peptide: Property::NotApplicable,
            formula: Property::Unresolved {
                recovered: None,
                annotation: "formula supplied in unsupported notation".into(),
            },
        };
        let mut columns = AnalyteColumns::default();
        columns.push(analyte.as_input());
        columns.validate().unwrap();
        let restored: Analyte =
            serde_json::from_str(&serde_json::to_string(&columns.get(0)).unwrap()).unwrap();
        assert_eq!(restored, analyte);
    }

    #[test]
    fn validation_checks_recovered_modification_sites() {
        let analyte = Analyte {
            peptide: Property::Unresolved {
                annotation: "partial".into(),
                recovered: Some(Peptide {
                    residues: "PEP".into(),
                    modifications: Property::Unresolved {
                        annotation: "partial mods".into(),
                        recovered: Some(vec![Modification {
                            site: ModificationSite::Residue(3),
                            definition: ModificationDefinition::unimod(4),
                        }]),
                    },
                }),
            },
            ..Default::default()
        };
        let mut columns = AnalyteColumns::default();
        columns.push(analyte.as_input());
        assert!(columns.validate().unwrap_err().contains("site 3"));
    }
}

#[cfg(test)]
#[path = "analyte_parser_tests.rs"]
mod parser_tests;
