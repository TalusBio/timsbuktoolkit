# Library analyte metadata

`Row.analyte` supplies independent optional peptide and molecular-formula facts.
`Row.entry_name` is a display label; `Row.id` is the external source key. Neither
is interpreted as chemistry. `RowIdx` addresses storage in one owning library;
it is not an external ID or a target/decoy relationship.

`TargetColumns::analyte(row)` returns borrowed `AnalyteRef` properties:

- `Missing`: no fact supplied.
- `NotApplicable`: explicitly inapplicable; absence alone does not establish this.
- `Known`: complete at that level. Known residues can have unresolved modifications.
- `Unresolved`: preserved chemical annotation, optionally with partial recovery.

Readers convert explicit sequence fields once. The arena packs residues and
located modifications, interns definitions, and remaps them when reader shards
merge. Original sequence spelling is discarded where structure is represented;
unresolved chemical content remains available through the property. Names and
source IDs retain their original text.

| Reader | Chemistry source | Name |
|---|---|---|
| DIA-NN TSV/Parquet | Modified and stripped sequence fields, checked for disagreement | `transition_group_id` / `Precursor.Id`, when supplied |
| DIA-NN binary | Format-defined modified-peptide-plus-charge field | Original field, unchanged |
| Spectronaut / Skyline | Explicit peptide fields | No separate entry-name field mapped |
| mzSpecLib | Already-parsed analyte; peptide or formula | Library spectrum name, when supplied |
| Prediction sink | Explicit ProForma, independently of generated label | Same label as prediction-file output |
| Target / ElutionGroupInput JSON arrays | Current schemas supply no analyte facts | Source ID remains available; no name inferred |

An mzSpecLib spectrum declaring multiple analytes is rejected rather than
silently selecting one. Unsupported peptide structures preserve their annotation.
Molecular formulas retain their declared basis (or `Unspecified`) and signed
electron counts. Peptide and formula facts can coexist: sealing rejects conflicting
neutral formulas for unmodified canonical peptides. Comparison of modified or
ambiguous structures and ion-basis formulas is deferred to composition support;
coexistence alone does not certify chemical consistency.

Search candidates carry a row and competition metadata, not copied sequences.
Rescorers and the dashboard receive the owning `ReferenceLibrary`. The existing
library-wide sequence gate reads stored structure: supported Unimod/mass
modifications and at most 254 residues. Missing, partial or unsupported structure
disables sequence features for the entire library. Global/labile/ambiguous
modifications are preserved as unresolved instead of silently omitted. The isotope
model still uses residues and its existing averagine fallback, not modification
or declared-formula composition.

Peptide deduplication compares complete structural keys plus charge and m/z,
then chooses the highest score, preferring a target on a tie. Equivalent
modification spellings therefore collapse. Incomparable entries remain separate;
formula equality and labels do not establish chemical identity or competition.

Results format version 4 resolves metadata from the library. `sequence` is nullable
and canonically formatted; `entry_name`, `molecular_formula`, and `formula_basis`
are nullable columns. Unformattable properties produce null. The viewer's Analyte
column displays sequence or formula. `Analyte` serialization preserves property
states and supports programmatic metadata round-trips; it does not change the
existing geometry-only target-list JSON schemas.
