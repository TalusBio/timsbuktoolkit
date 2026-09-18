use arrow::array::{
    ArrayRef,
    new_empty_array,
};
use arrow::datatypes::*;
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use timsquery::IonAnnot;
use timsquery::models::{
    RowIdx,
    SourceId,
    TargetColumns,
};
use tracing::debug;

use super::blocks::{
    ColSink,
    NameSink,
    SchemaSink,
    ScoreBlock,
};
use super::results::FinalResult;
use crate::sample_identity::SampleIdentity;

/// Bumped when a column's meaning or type changes, so a reader can tell a new
/// file from an old one rather than silently misreading it.
///
/// - 2: `library_id` and `decoy_group_id` widened to UInt64. `library_id` was
///   already the caller's id but was truncated to 32 bits; `decoy_group_id`
///   used to be the arena row position and is now an id declared by the input
///   or minted at load.
/// - 3: both became Utf8, because an id keeps the shape its source used and
///   DIA-NN names its precursors with a string (`transition_group_id`). A
///   numeric id is written as its digits.
/// - 4: sequence is canonical and nullable, resolved from stored analyte facts;
///   entry_name, molecular_formula and formula_basis are separate nullable columns.
///
/// `decoy_group_id` equals `library_id` on a row whose library declared no
/// competition group: the row competes only with its own mass-shift variants
/// and is its own group. mzSpecLib can declare a group (`related spectrum
/// keys`), and then a target and its shipped decoy share one `decoy_group_id`
/// across two `library_id`s. Other formats declare nothing, so for them the
/// duplication is expected, not a bug.
pub const RESULTS_FORMAT_VERSION: u32 = 4;

// ---------------------------------------------------------------------------
// Build a RecordBatch from a slice of FinalResult
// ---------------------------------------------------------------------------

/// The ids a reader joins on, resolved from the arena.
///
/// A block rather than two hand-written `o.str` calls in `emit_row` plus two in
/// `FinalResult::column_schema`: that split is exactly the schema/data drift
/// `#[derive(ScoreBlock)]` exists to prevent. A block keeps the pair adjacent,
/// so the two bodies cannot disagree about names, types, or order.
///
/// Hand-written, like `Identity`, because the derive works on owned scalar
/// fields and these are borrowed from the arena.
pub(crate) struct Ids<'a> {
    library_id: SourceId<'a>,
    decoy_group_id: SourceId<'a>,
}

impl<'a> Ids<'a> {
    fn for_row(geom: &'a TargetColumns<IonAnnot>, row: RowIdx) -> Self {
        Self {
            library_id: geom.output_id(row),
            decoy_group_id: geom.decoy_group(row),
        }
    }

    /// A text id goes in borrowed; only a numeric one is rendered.
    fn emit(o: &mut ColSink, name: &str, id: SourceId<'_>) {
        match id {
            SourceId::Text(s) => o.str(name, s),
            SourceId::Numeric(n) => o.str(name, &n.to_string()),
        }
    }
}

impl ScoreBlock for Ids<'_> {
    fn columns(&self, o: &mut ColSink) {
        Self::emit(o, "library_id", self.library_id);
        Self::emit(o, "decoy_group_id", self.decoy_group_id);
    }

    fn column_schema(o: &mut SchemaSink) {
        o.str("library_id");
        o.str("decoy_group_id");
    }

    fn nonlinear_feature_names(_: &mut NameSink) {}
}

/// Context-dependent metadata; schema and values stay together like score blocks.
pub(crate) struct AnalyteColumns<'a> {
    analyte: timsquery::chemistry::analyte::AnalyteRef<'a>,
    entry_name: Option<&'a str>,
}
impl ScoreBlock for AnalyteColumns<'_> {
    fn columns(&self, sink: &mut ColSink) {
        let sequence = self.analyte.peptide.known().and_then(|p| p.sequence());
        let formula = self.analyte.formula.known().and_then(|f| f.notation());
        sink.optional_str("sequence", sequence.as_deref());
        sink.optional_str("entry_name", self.entry_name);
        sink.optional_str("molecular_formula", formula.as_deref());
        sink.optional_str(
            "formula_basis",
            self.analyte.formula.known().map(|f| f.basis.as_str()),
        );
    }

    fn column_schema(sink: &mut SchemaSink) {
        for name in [
            "sequence",
            "entry_name",
            "molecular_formula",
            "formula_basis",
        ] {
            sink.optional_str(name);
        }
    }

    fn nonlinear_feature_names(_: &mut super::blocks::NameSink) {}
}

/// Emit identifiers and analyte metadata first, then scoring and result metadata.
///
/// Resolve `library_id` and `decoy_group_id` through the result's arena row.
/// The owned source ID retained for rescoring order does not supply these columns.
fn emit_row(r: &FinalResult, geom: &TargetColumns<IonAnnot>, sink: &mut ColSink) {
    Ids::for_row(geom, r.scoring.identity.row).columns(sink);
    AnalyteColumns {
        analyte: geom.analyte(r.scoring.identity.row),
        entry_name: geom.entry_name(r.scoring.identity.row),
    }
    .columns(sink);
    r.scoring.columns(sink);
    r.result_meta().columns(sink);
    sink.end_row();
}

/// Zero-row RecordBatch built directly from the block schema -- no fake row.
fn empty_batch() -> std::io::Result<RecordBatch> {
    let mut ss = super::blocks::SchemaSink::new();
    FinalResult::column_schema(&mut ss);
    let fields = ss.into_fields();
    let arrays: Vec<ArrayRef> = fields
        .iter()
        .map(|f| new_empty_array(f.data_type()))
        .collect();
    let schema = Arc::new(Schema::new(fields));
    RecordBatch::try_new(schema, arrays)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

/// Build a `RecordBatch` from a slice of `FinalResult`.
///
/// **COMPILE-TIME SAFETY:** each block's `columns()` is generated by walking
/// the same field list that defines the block struct (see
/// `#[derive(ScoreBlock)]` in `timsseek_macros`), so a new field cannot
/// silently skip the Parquet projection -- there is no hand-written destructure
/// to fall out of sync.
pub fn build_record_batch(
    results: &[FinalResult],
    geom: &TargetColumns<IonAnnot>,
) -> std::io::Result<RecordBatch> {
    if results.is_empty() {
        return empty_batch();
    }
    let mut sink = ColSink::new();
    for r in results {
        emit_row(r, geom, &mut sink);
    }

    let (fields, arrays) = sink.finish();
    let schema = Arc::new(Schema::new(fields));
    RecordBatch::try_new(schema, arrays)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

fn output_batch(
    results: &[FinalResult],
    geom: &TargetColumns<IonAnnot>,
    raw: bool,
) -> std::io::Result<RecordBatch> {
    let batch = build_record_batch(results, geom)?;
    if !raw {
        return Ok(batch);
    }
    let mut schema = SchemaSink::new();
    super::blocks::result_meta::ResultMeta::column_schema(&mut schema);
    let omitted = schema.into_fields();
    let keep: Vec<_> = batch
        .schema()
        .fields()
        .iter()
        .enumerate()
        .filter_map(|(i, field)| {
            (!omitted.iter().any(|omit| omit.name() == field.name())).then_some(i)
        })
        .collect();
    batch.project(&keep).map_err(std::io::Error::other)
}

// ---------------------------------------------------------------------------
// Buffered Parquet writer
// ---------------------------------------------------------------------------

/// Borrows the arena the results were scored against to resolve each result's
/// row handle into the external library and competition-group IDs.
pub struct ResultParquetWriter<'a> {
    writer: ArrowWriter<File>,
    buffer: Vec<FinalResult>,
    row_group_size: usize,
    geom: &'a TargetColumns<IonAnnot>,
    raw: bool,
}

impl<'a> ResultParquetWriter<'a> {
    /// Identity is required at construction, including for zero-row artifacts.
    pub fn new(
        path: impl AsRef<Path>,
        row_group_size: usize,
        library: &'a crate::data_sources::reference_library::ReferenceLibrary,
        sample: &SampleIdentity,
    ) -> std::io::Result<Self> {
        Self::with_mode(path, row_group_size, library, sample, false)
    }

    /// Common scores only: no competition, discriminant score or q-value columns.
    pub fn raw(
        path: impl AsRef<Path>,
        row_group_size: usize,
        library: &'a crate::ReferenceLibrary,
        sample: &SampleIdentity,
    ) -> std::io::Result<Self> {
        Self::with_mode(path, row_group_size, library, sample, true)
    }

    fn with_mode(
        path: impl AsRef<Path>,
        row_group_size: usize,
        library: &'a crate::ReferenceLibrary,
        sample: &SampleIdentity,
        raw: bool,
    ) -> std::io::Result<Self> {
        let geom = library.geometry();
        let file = match File::create_new(path.as_ref()) {
            Ok(f) => f,
            Err(err) => {
                tracing::error!("Failed to create file {:?}: {}", path.as_ref(), err);
                return Err(err);
            }
        };

        // Build schema from a zero-row batch
        let empty_batch = output_batch(&[], geom, raw)?;
        let schema = empty_batch.schema();

        let kv = vec![
            KeyValue {
                key: "sample_id".into(),
                value: Some(sample.sample_id().to_owned()),
            },
            KeyValue {
                key: "sample_name".into(),
                value: Some(sample.sample_name().to_owned()),
            },
            KeyValue {
                key: "result_mode".into(),
                value: Some(if raw { "raw" } else { "rescored" }.into()),
            },
            KeyValue {
                key: "scoring_plan".into(),
                value: Some(
                    serde_json::to_string(library.scoring_plan()).map_err(std::io::Error::other)?,
                ),
            },
            KeyValue {
                key: "all_sequence_counts_enabled".to_string(),
                value: Some(library.all_sequence_counts_enabled().to_string()),
            },
            KeyValue {
                key: "results_format_version".to_string(),
                value: Some(RESULTS_FORMAT_VERSION.to_string()),
            },
        ];
        let props = WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .set_key_value_metadata(Some(kv))
            .build();

        let writer =
            ArrowWriter::try_new(file, schema, Some(props)).map_err(std::io::Error::other)?;

        Ok(Self {
            writer,
            buffer: Vec::with_capacity(row_group_size),
            row_group_size,
            geom,
            raw,
        })
    }

    pub fn add_raw(&mut self, result: super::results::ScoredCandidate) -> std::io::Result<()> {
        assert!(self.raw, "raw scores require the raw output schema");
        self.add(FinalResult {
            scoring: result.scoring,
            delta_group_ln1p_diff: f32::NAN,
            delta_group_ln1p_ratio: f32::NAN,
            discriminant_score: f32::NAN,
            qvalue: f32::NAN,
        })
    }

    pub fn add(&mut self, result: FinalResult) -> std::io::Result<()> {
        self.buffer.push(result);
        if self.buffer.len() >= self.row_group_size {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        debug!("Flushing {} results to parquet", self.buffer.len());
        let batch = output_batch(&self.buffer, self.geom, self.raw)?;
        self.writer.write(&batch).map_err(std::io::Error::other)?;
        self.buffer.clear();
        Ok(())
    }

    pub fn close(mut self) -> std::io::Result<()> {
        self.flush()?;
        self.writer.close().map_err(std::io::Error::other)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use parquet::file::reader::{
        FileReader,
        SerializedFileReader,
    };
    use std::fs::File;
    use timsquery::models::capabilities::TargetCapabilities;
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };

    #[test]
    fn sample_metadata_survives_empty_and_nonempty_files_in_both_modes() {
        let geom = one_row_arena();
        let library = crate::ReferenceLibrary::try_from(timsquery::serde::TargetTable::Mzpaf {
            frag_intens: Some(vec![1.0; geom.n_fragments()]),
            geom,
        })
        .unwrap();
        let sample =
            crate::sample_identity::SampleIdentity::from_location("s3://bucket/run.d").unwrap();
        let dir = tempfile::tempdir().unwrap();
        for raw in [false, true] {
            for nrows in [0, 1] {
                let path = dir.path().join(format!("{raw}-{nrows}.parquet"));
                let mut writer = if raw {
                    ResultParquetWriter::raw(&path, 1, &library, &sample)
                } else {
                    ResultParquetWriter::new(&path, 1, &library, &sample)
                }
                .unwrap();
                if nrows == 1 {
                    writer.add(sample_in(library.geometry())).unwrap();
                }
                writer.close().unwrap();
                let reader = SerializedFileReader::new(File::open(path).unwrap()).unwrap();
                let meta = reader.metadata().file_metadata();
                assert_eq!(meta.num_rows(), nrows);
                for (key, value) in [
                    ("sample_id", "07cff6d98863b0e4-run"),
                    ("sample_name", "run"),
                ] {
                    let entries: Vec<_> = meta
                        .key_value_metadata()
                        .unwrap()
                        .iter()
                        .filter(|kv| kv.key == key)
                        .collect();
                    assert_eq!(entries.len(), 1);
                    assert_eq!(entries[0].value.as_deref(), Some(value));
                }
            }
        }
    }

    /// A sealed arena with one row per `(sequence, id)`, for the writer to
    /// resolve ids against. `None` for an id leaves the row to be minted.
    #[test]
    fn raw_output_omits_fdr_and_competition_columns_including_empty_files() {
        let geom = one_row_arena();
        let row = sample_in(&geom);
        for rows in [&[][..], std::slice::from_ref(&row)] {
            let batch = output_batch(rows, &geom, true).unwrap();
            assert_eq!(batch.num_rows(), rows.len());
            for name in [
                "qvalue",
                "discriminant_score",
                "delta_group_ln1p_diff",
                "delta_group_ln1p_ratio",
            ] {
                assert!(batch.schema().field_with_name(name).is_err());
            }
            assert!(batch.schema().field_with_name("main_score").is_ok());
            assert!(batch.schema().field_with_name("library_id").is_ok());
        }
    }

    fn arena_of(rows: &[(&str, Option<&str>)]) -> TargetColumns<IonAnnot> {
        let mut geom = TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        for (seq, id) in rows {
            geom.push_row(Row {
                precursor_mz: 900.4,
                charge: 2,
                rt_seconds: 1.0,
                mobility: 1.0,
                frags: &[(IonAnnot::try_from("y3").unwrap(), 300.0)],
                analyte: timsquery::chemistry::analyte::Analyte::from_sequence(seq).as_input(),
                id: id.map(Into::into),
                ..Default::default()
            });
        }
        geom.seal(crate::models::DecoyPolicy::Never)
            .expect("fixture ids are usable")
    }

    fn one_row_arena() -> TargetColumns<IonAnnot> {
        arena_of(&[("PEPTIDEK", None)])
    }

    /// `FinalResult::sample()` carries a placeholder row, which by design reads
    /// no arena. Point it at a real one before writing.
    fn sample_in(geom: &TargetColumns<IonAnnot>) -> FinalResult {
        let mut r = FinalResult::sample();
        r.scoring.identity.row = geom.rows().next().expect("one row");
        r
    }

    /// Golden byte-compat contract for the Parquet schema.
    ///
    /// The maintained contract: a SUBSET of columns downstream consumers rely on
    /// to join a row and threshold it. Each listed `column_name -> (dtype,
    /// nullable)` must be present and match. Everything else the writer emits
    /// (all the feature/diagnostic columns) is free to change or disappear --
    /// this test does NOT reject extra columns, and column ORDER is not asserted
    /// (downstream reads by name). Add a column here only when you commit to
    /// keeping it stable.
    const GOLDEN_SCHEMA: &[(&str, &str, bool)] = &[
        ("sequence", "Utf8", true),
        ("entry_name", "Utf8", true),
        ("molecular_formula", "Utf8", true),
        ("formula_basis", "Utf8", true),
        ("library_id", "Utf8", false),
        ("decoy_group_id", "Utf8", false),
        ("is_target", "Boolean", false),
        ("precursor_mz", "Float64", false),
        ("precursor_charge", "UInt8", false),
        ("main_score", "Float32", false),
        ("discriminant_score", "Float32", false),
        ("qvalue", "Float32", false),
    ];

    #[test]
    fn formula_only_output_keeps_label_and_has_null_sequence() {
        use arrow::array::{
            Array,
            StringArray,
        };
        use mzcore::chemistry::Element::{
            C,
            H,
            O,
        };
        use timsquery::chemistry::analyte::{
            Analyte,
            Formula,
            FormulaBasis,
            Property,
        };
        let analyte = Analyte {
            formula: Property::Known(Formula {
                elements: vec![(H, None, 6), (C, None, 2), (O, None, 1)],
                basis: FormulaBasis::NeutralMolecule,
            }),
            ..Default::default()
        };
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        builder.push_row(Row {
            analyte: analyte.as_input(),
            id: Some("source/123".into()),
            entry_name: Some("Compound label / not sequence"),
            ..Default::default()
        });
        let geom = builder.seal(crate::models::DecoyPolicy::Never).unwrap();
        let batch = build_record_batch(&[sample_in(&geom)], &geom).unwrap();
        let column = |name: &str| {
            batch
                .column_by_name(name)
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
        };
        assert!(column("sequence").is_null(0));
        assert_eq!(column("library_id").value(0), "source/123");
        assert_eq!(
            column("entry_name").value(0),
            "Compound label / not sequence"
        );
        assert_eq!(column("molecular_formula").value(0), "C2H6O1");
        assert_eq!(column("formula_basis").value(0), "neutral_molecule");
        assert_eq!(
            batch.schema(),
            build_record_batch(&[], &geom).unwrap().schema()
        );
    }

    #[test]
    fn schema_first_matches_populated_data_path() {
        // The writer's schema comes from the empty (SchemaSink) batch; real data
        // goes through the ColSink path. They MUST have identical schemas (name,
        // dtype, nullability, ORDER) or the writer rejects populated batches.
        let geom = one_row_arena();
        let empty = build_record_batch(&[], &geom).expect("empty");
        let populated = build_record_batch(&[sample_in(&geom)], &geom).expect("populated");
        assert_eq!(
            empty.schema(),
            populated.schema(),
            "schema-first empty batch must byte-match the ColSink data-path schema"
        );
    }

    /// Output IDs must come from the arena row referenced by the result.
    #[test]
    fn ids_resolve_from_the_arena_row() {
        use arrow::array::StringArray;

        // The names DIA-NN would give these rows, which is the case that must
        // not come back as digits.
        let geom = arena_of(&[
            ("PEPTIDEK", Some("PEPTIDEK2")),
            ("AAAAAAALQAK", Some("AAAAAAALQAK2")),
        ]);

        let second = geom.rows().nth(1).expect("two rows");
        let mut result = FinalResult::sample();
        result.scoring.identity.row = second;

        let batch = build_record_batch(&[result], &geom).expect("batch");
        let ids = batch
            .column_by_name("library_id")
            .expect("library_id column")
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("library_id is Utf8");
        assert_eq!(ids.value(0), "AAAAAAALQAK2");

        // No declared groups, so the row is its own group and the group id is
        // the id the row reports.
        let groups = batch
            .column_by_name("decoy_group_id")
            .expect("decoy_group_id column")
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("decoy_group_id is Utf8");
        assert_eq!(groups.value(0), "AAAAAAALQAK2");
    }

    #[test]
    fn parquet_schema_contains_maintained_subset() {
        use std::collections::BTreeMap;
        let batch = build_record_batch(&[], &one_row_arena()).expect("schema");
        let got: BTreeMap<String, (String, bool)> = batch
            .schema()
            .fields()
            .iter()
            .map(|f| {
                (
                    f.name().clone(),
                    (format!("{:?}", f.data_type()), f.is_nullable()),
                )
            })
            .collect();

        for (name, dtype, nullable) in GOLDEN_SCHEMA {
            let got_ty = got
                .get(*name)
                .unwrap_or_else(|| panic!("maintained column missing: `{name}`"));
            assert_eq!(
                got_ty,
                &(dtype.to_string(), *nullable),
                "dtype/nullable mismatch for `{name}`"
            );
        }
    }

    #[test]
    fn delta_columns_name_their_ln1p_formulas() {
        let batch = build_record_batch(&[], &one_row_arena()).expect("schema");
        let schema = batch.schema();

        assert!(schema.index_of("delta_group_ln1p_diff").is_ok());
        assert!(schema.index_of("delta_group_ln1p_ratio").is_ok());
        assert!(schema.index_of("delta_group").is_err());
        assert!(schema.index_of("delta_group_ratio").is_err());
    }

    #[test]
    fn all_sequence_counts_enabled_key_in_parquet_metadata() {
        let tmp = tempfile::NamedTempFile::new().expect("tmpfile");
        let path = tmp.path().to_path_buf();
        // Close the NamedTempFile so ResultParquetWriter can create the file
        // (File::create_new fails if the path already exists on some platforms,
        //  but NamedTempFile holds it open; drop it first).
        drop(tmp);
        {
            let geom = one_row_arena();
            let library = crate::data_sources::reference_library::ReferenceLibrary::try_from(
                timsquery::serde::TargetTable::Mzpaf {
                    frag_intens: Some(vec![1.0; geom.n_fragments()]),
                    geom,
                },
            )
            .unwrap();
            let sample = SampleIdentity::from_location("s3://bucket/run.d").unwrap();
            let writer =
                ResultParquetWriter::new(&path, 1024, &library, &sample).expect("create writer");
            writer.close().expect("close");
        }
        let file = File::open(&path).expect("open");
        let reader = SerializedFileReader::new(file).expect("reader");
        let meta = reader.metadata().file_metadata();
        let kv_list = meta.key_value_metadata().expect("kv metadata present");
        let plan = kv_list
            .iter()
            .find(|k| k.key == "scoring_plan")
            .expect("plan metadata");
        let plan: serde_json::Value = serde_json::from_str(plan.value.as_deref().unwrap()).unwrap();
        assert_eq!(plan["rows"], 1);
        assert_eq!(plan["isotopes"]["method"], "composition_cs");
        assert_eq!(plan["isotopes"]["composition_rows"], 1);
        assert!(plan["isotopes"].get("envelopes").is_none());
        assert_eq!(plan["unmodified_rows"], 1);
        assert_eq!(plan["operations"][0]["requirement"], "residue_sequence");
        assert_eq!(
            plan["operations"][0]["columns"].as_array().unwrap().len(),
            21
        );
        assert_eq!(plan["operations"][1]["enabled"], true);
        assert!(!kv_list.iter().any(|k| k.key == "parsable_sequences"));
        let found: Vec<_> = kv_list
            .iter()
            .filter(|k| k.key == "all_sequence_counts_enabled")
            .collect();
        assert_eq!(
            found.len(),
            1,
            "expected exactly one all_sequence_counts_enabled key"
        );
        assert_eq!(
            found[0].value.as_deref(),
            Some("true"),
            "all_sequence_counts_enabled value should be 'true'"
        );
    }
}
