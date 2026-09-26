//! A library predicted straight into the arena timsseek scores, with no file in
//! between.
//!
//! msspeculator hands every precursor to a [`LibrarySink`]; this one spells the
//! rows into [`timsquery::models::TargetColumns`] instead of into mzSpecLib or DIA-NN text. The
//! result has to agree with loading the same library off disk on every value a
//! scorer reads, so every spelling decision here mirrors the mzSpecLib reader in
//! `timsquery::serde::mzspeclib_io` -- the two routes agree on the numbers they
//! store, not on the strings they pass through.
//!
//! The sink is moved onto msspeculator's writer thread and never handed back, so
//! the finished rows come out through shared state: [`sink`] returns a handle
//! that stays with the caller and the sink that goes into `stream_library`.

use std::sync::{
    Arc,
    Mutex,
};

use anyhow::{
    Context,
    Result,
    bail,
};
use msspeculator_inference::{
    LibraryProvenance,
    LibrarySink,
    LibraryStats,
    Peak,
    SpectrumRow,
};
use timsquery::IonAnnot;
use timsquery::ion::{
    IonSeriesOrdinal,
    NeutralLoss,
};
use timsquery::models::{
    Row,
    TargetCapabilities,
    TargetColumnsBuilder,
};
use timsquery::serde::{
    TargetReadingError,
    TargetTable,
};
use timsseek::DecoyPolicy;
use timsseek::data_sources::reference_library::ReferenceLibrary;

use crate::errors::CliError;

/// Convert predicted gradient minutes to acquisition seconds; indices remain unchanged.
const SECONDS_PER_MINUTE: f32 = 60.0;

/// A library that was predicted rather than read, plus what produced it.
pub(crate) struct PredictedLibrary {
    pub library: ReferenceLibrary,
    pub provenance: LibraryProvenance,
}

/// The caller's half of [`sink`]: what the prediction produced, once it has.
pub(crate) struct PredictedLibraryHandle {
    shared: Arc<Mutex<Handoff>>,
}

/// The [`LibrarySink`] half of [`sink`], for `stream_library` to consume.
pub(crate) struct PredictedLibrarySink {
    shared: Arc<Mutex<Handoff>>,
    state: SinkState,
}

enum SinkState {
    AwaitingHeader {
        arena: ArenaRows,
    },
    Building {
        arena: ArenaRows,
        provenance: Box<LibraryProvenance>,
        normalized_axis: timsquery::RtAxis,
    },
    Finished,
}

/// What crosses from the writer thread back to the caller.
///
/// `completed` is `Some` only once `finish` ran, which is not the same as the stream
/// having succeeded: msspeculator's `run_library` joins its inference workers
/// before its writer, so a worker that panics drops the result channel, the
/// writer's loop over it ends cleanly, and `finish` publishes whatever arrived
/// before the panic. `stream_library` still returns the error, which is why
/// [`PredictedLibraryHandle::into_library`] takes the [`LibraryStats`] only a
/// successful stream produces and checks its count against what arrived.
#[derive(Default)]
struct Handoff {
    completed: Option<Completed>,
}

struct Completed {
    provenance: Box<LibraryProvenance>,
    arena: ArenaRows,
}

/// Build a sink and the handle that collects from it.
pub(crate) fn sink(decoys: DecoyPolicy) -> (PredictedLibraryHandle, PredictedLibrarySink) {
    let shared = Arc::new(Mutex::new(Handoff::default()));
    (
        PredictedLibraryHandle {
            shared: Arc::clone(&shared),
        },
        PredictedLibrarySink {
            shared,
            state: SinkState::AwaitingHeader {
                arena: ArenaRows::new(decoys),
            },
        },
    )
}

impl PredictedLibraryHandle {
    /// Seal what the prediction handed over into a scorable library.
    ///
    /// `stats` is the witness that the stream succeeded: it is what
    /// `stream_library` returns, and it returns nothing on the error path. Taking
    /// it by reference is what stops the sequence a `finish`-ran flag cannot
    /// catch -- log the error from `stream_library`, then seal the rows that
    /// arrived before it -- because there is no `LibraryStats` to pass.
    ///
    /// `stats.precursors` is then checked against the rows that arrived.
    /// msspeculator increments it once per `spectrum` call, immediately before
    /// making it, so the two are the same count of the same thing: what the sink
    /// was handed, ahead of the decoy policy dropping any of it.
    pub(crate) fn into_library(self, stats: &LibraryStats) -> Result<PredictedLibrary, CliError> {
        let completed = self
            .shared
            .lock()
            .expect("handoff mutex poisoned")
            .completed
            .take();
        let Some(Completed {
            provenance,
            arena: arena_rows,
        }) = completed
        else {
            return Err(CliError::LibraryBuild {
                source: "the prediction stream handed over no rows, so it did not finish"
                    .to_string(),
            });
        };
        if arena_rows.received != stats.precursors {
            return Err(CliError::LibraryBuild {
                source: format!(
                    "the prediction reported {} precursors but the sink received {}, so what it \
                     handed over is a prefix of the library that was asked for",
                    stats.precursors, arena_rows.received,
                ),
            });
        }
        let arena = arena_rows.seal().map_err(|e| CliError::LibraryBuild {
            source: format!("assembling the predicted library: {e:?}"),
        })?;
        let library = ReferenceLibrary::try_from(arena).map_err(|e| CliError::LibraryBuild {
            source: format!("finalizing the predicted library: {e:?}"),
        })?;
        Ok(PredictedLibrary {
            library,
            provenance: *provenance,
        })
    }
}

impl LibrarySink for PredictedLibrarySink {
    fn header(&mut self, provenance: &LibraryProvenance) -> Result<()> {
        match std::mem::replace(&mut self.state, SinkState::Finished) {
            SinkState::AwaitingHeader { arena } => {
                self.state = SinkState::Building {
                    arena,
                    provenance: Box::new(provenance.clone()),
                    normalized_axis: timsquery::RtAxis::NormalizedIndex {
                        scale: Some(provenance.retention.normalized.scale.to_owned()),
                    },
                };
                Ok(())
            }
            _ => bail!("prediction header called more than once or after finish"),
        }
    }

    fn spectrum(&mut self, row: &SpectrumRow<'_>) -> Result<()> {
        match &mut self.state {
            SinkState::Building {
                arena,
                normalized_axis,
                ..
            } => {
                arena.push(PredictedRow::from_spectrum(row)?, normalized_axis);
                Ok(())
            }
            SinkState::AwaitingHeader { .. } => bail!("prediction spectrum before header"),
            SinkState::Finished => bail!("prediction spectrum after finish"),
        }
    }

    fn finish(&mut self) -> Result<()> {
        match std::mem::replace(&mut self.state, SinkState::Finished) {
            SinkState::Building {
                arena, provenance, ..
            } => {
                self.shared
                    .lock()
                    .expect("handoff mutex poisoned")
                    .completed = Some(Completed { provenance, arena });
                Ok(())
            }
            SinkState::AwaitingHeader { .. } => bail!("prediction finished before header"),
            SinkState::Finished => bail!("prediction finished more than once"),
        }
    }
}

/// One prediction, owned and already spelled the way the arena stores it.
#[derive(Debug)]
struct PredictedRow {
    /// `{proforma}/{charge}`, which is the name the mzSpecLib writer gives the
    /// same precursor under `MS:1003061|library spectrum name`, so a row keeps
    /// one name across both routes. Stored independently of chemistry.
    id: String,
    /// The source peptide's pair id, `None` when decoys are off. All modified
    /// forms and all charges of one peptide share it, so a group is a peptide
    /// rather than a target/decoy couple.
    group: Option<String>,
    analyte: timsquery::chemistry::analyte::Analyte,
    precursor_mz: f64,
    charge: u8,
    rt_value: f32,
    rt_axis: timsquery::RtAxis,
    mobility: f32,
    frags: Vec<(IonAnnot, f64)>,
    /// Parallel to `frags`, kept apart because the arena stores the two in
    /// different places: labels and m/z in the columns, intensities in a sidecar.
    intensities: Vec<f32>,
    is_decoy: bool,
}

impl PredictedRow {
    fn from_spectrum(row: &SpectrumRow<'_>) -> Result<Self> {
        // Preserve the writer's label; chemistry comes from the explicit proforma.
        let id = format!("{}/{}", row.proforma, row.charge);
        let charge = u8::try_from(row.charge)
            .with_context(|| format!("precursor charge {} of {id}", row.charge))?;

        let mut frags = Vec::with_capacity(row.peaks.len());
        let mut intensities = Vec::with_capacity(row.peaks.len());
        for peak in &row.peaks {
            frags.push((
                ion_annot(peak).with_context(|| format!("fragment of {id}"))?,
                peak.mz,
            ));
            intensities.push(peak.intensity as f32);
        }

        Ok(Self {
            // `Display` reads a decoy's residues with its interior reversed, so
            // this is the decoy's own sequence rather than its target's.
            analyte: timsquery::chemistry::analyte::Analyte::from_sequence(row.proforma),
            id,
            group: row.decoy_pair_id.map(|pair| pair.to_string()),
            precursor_mz: row.precursor_mz,
            charge,
            // `rt` is the quantity the mzSpecLib writer puts under the term the
            // file route reads back as `start_time`: gradient minutes with a
            // chromatography context, the normalized index without one. `irt`
            // carries the index alongside a gradient time, under a term the
            // reader keeps as a scan param and never scores on.
            rt_value: if row.irt.is_some() {
                row.rt * SECONDS_PER_MINUTE
            } else {
                row.rt
            },
            rt_axis: if row.irt.is_some() {
                timsquery::RtAxis::Seconds
            } else {
                timsquery::RtAxis::NormalizedIndex { scale: None }
            },
            mobility: row.mobility as f32,
            frags,
            intensities,
            is_decoy: row.decoy,
        })
    }
}

/// The packed label for a predicted peak.
///
/// A failure is an error rather than the `?N` placeholder the file route mints:
/// a library on disk can carry an annotation this build cannot spell, while a
/// prediction carries a `b` or a `y` with no loss and no isotope, so anything
/// else here means the peak is not what the row says it is. Silently labelling
/// it unknown would keep its intensity in the sidecar under a label nothing can
/// match.
fn ion_annot(peak: &Peak<'_>) -> Result<IonAnnot> {
    let mut letters = peak.ion.chars();
    let (Some(series), None) = (letters.next(), letters.next()) else {
        bail!(
            "ion series {:?} is not a single mzPAF series letter",
            peak.ion
        );
    };
    let ordinal =
        u8::try_from(peak.ordinal).with_context(|| format!("ordinal {}", peak.ordinal))?;
    let charge = i8::try_from(peak.charge).with_context(|| format!("charge {}", peak.charge))?;
    let series = IonSeriesOrdinal::from_series_char(series, Some(ordinal))
        .map_err(|e| anyhow::anyhow!("series {series}{ordinal}: {e}"))?;
    IonAnnot::new(series, NeutralLoss::None, charge, 0)
        .map_err(|e| anyhow::anyhow!("packing {series}^{charge}: {e}"))
}

/// Writer-owned columns. Each row enters the builder during `spectrum`.
struct ArenaRows {
    geom: TargetColumnsBuilder<IonAnnot>,
    frag_intens: Vec<f32>,
    received: usize,
    decoys: DecoyPolicy,
}

impl ArenaRows {
    fn new(decoys: DecoyPolicy) -> Self {
        Self {
            geom: TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann()),
            frag_intens: Vec::new(),
            received: 0,
            decoys,
        }
    }

    fn push(&mut self, row: PredictedRow, normalized_axis: &timsquery::RtAxis) {
        self.received += 1;
        // Count every prediction, including decoys excluded by the policy.
        if !self.decoys.accepts(row.is_decoy) {
            return;
        }
        let group = row.group.clone().unwrap_or_else(|| row.id.clone());
        self.frag_intens.extend_from_slice(&row.intensities);
        self.geom.push_row(Row {
            precursor_mz: row.precursor_mz,
            charge: row.charge,
            rt: Some(timsquery::RtCoordinate {
                value: timsquery::LibraryRT(row.rt_value),
                axis: if matches!(row.rt_axis, timsquery::RtAxis::NormalizedIndex { .. }) {
                    normalized_axis
                } else {
                    &row.rt_axis
                },
            }),
            mobility: row.mobility,
            frags: &row.frags,
            analyte: row.analyte.as_input(),
            entry_name: Some(&row.id),
            is_decoy: row.is_decoy,
            id: Some(row.id.clone().into()),
            decoy_group: Some(group.into()),
        });
    }

    fn seal(self) -> Result<TargetTable, TargetReadingError> {
        let Self {
            geom,
            frag_intens,
            received,
            decoys,
        } = self;
        // A zero-row arena seals: ids and groups are both vacuously consistent, the
        // parse gate has nothing to reject, and the result searches to zero results
        // without ever reporting why. msspeculator only refuses an empty digest, so
        // a FASTA whose peptides all fall outside the windows gets this far.
        if received == 0 {
            return Err(TargetReadingError::SpeclibParse(
            "the prediction produced no precursors; check the charge range and the length range \
             against the peptides the FASTA actually digests to"
                .to_string(),
        ));
        }
        // The same empty arena the guard above refuses, reached the other way round:
        // every row that arrived was a shipped decoy and the policy dropped all of
        // them, leaving nothing to derive decoys against.
        if geom.n_rows() == 0 {
            return Err(TargetReadingError::SpeclibParse(format!(
                "all {} predicted rows were shipped decoys, which {decoys:?} drops, so the library \
             holds no targets",
                received,
            )));
        }

        let geom = geom.seal(decoys)?;
        if frag_intens.len() != geom.n_fragments() {
            return Err(TargetReadingError::SpeclibParse(format!(
                "reference-intensity sidecar ({}) must stay parallel to the fragment-label arena ({})",
                frag_intens.len(),
                geom.n_fragments(),
            )));
        }

        Ok(TargetTable::Mzpaf {
            geom,
            frag_intens: Some(frag_intens),
        })
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::sync::LazyLock;

    use msspeculator_core::peptide::Peptide;
    use msspeculator_inference::{
        ProteinGroup,
        Residues,
    };
    use timsquery::models::{
        DecoyStrategy,
        RowIdx,
        TargetColumns,
    };

    use super::*;

    static TEST_PROVENANCE: LazyLock<LibraryProvenance> = LazyLock::new(|| {
        let fasta =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test_data/tiny.fasta");
        let mut prediction = crate::build_library::resolve_prediction(
            fasta,
            &crate::config::LibraryConfig::default(),
        );
        prediction.decoys = false;
        prediction.max_fragments = Some(4);
        crate::build_library::predict_in_memory(&prediction, DecoyPolicy::Never)
            .expect("fixture provenance")
            .provenance
    });

    fn test_provenance(scale: Option<&'static str>) -> LibraryProvenance {
        let mut provenance = (*TEST_PROVENANCE).clone();
        if let Some(scale) = scale {
            provenance.retention.normalized.scale = scale;
        }
        provenance
    }

    impl PredictedLibrarySink {
        fn record_test_provenance(&mut self) {
            self.header(&test_provenance(None)).unwrap();
        }
    }

    /// A prediction row assembled by hand.
    ///
    /// `SpectrumRow`'s fields are public and `Residues`/`ProteinGroup` have
    /// constructors for a caller building one itself, so these tests drive the
    /// trait method rather than a stand-in for it. The one field they cannot
    /// reach is `header`'s `LibraryProvenance`, which only msspeculator can
    /// build.
    struct Fixture {
        peptide: Peptide,
        proteins: Vec<String>,
        members: Vec<u32>,
        proforma: String,
        stripped: String,
    }

    impl Fixture {
        fn new(proforma: &str, stripped: &str) -> Self {
            Self {
                peptide: Peptide::new(stripped.to_string(), Vec::new()),
                proteins: vec!["sp|P00001|TEST".to_string()],
                members: vec![0],
                proforma: proforma.to_string(),
                stripped: stripped.to_string(),
            }
        }

        fn row<'a>(
            &'a self,
            charge: i64,
            decoy: bool,
            pair: Option<usize>,
            peaks: Vec<Peak<'a>>,
        ) -> SpectrumRow<'a> {
            SpectrumRow {
                stripped: Residues::target(&self.stripped),
                proteins: ProteinGroup::new(&self.proteins, &self.members, decoy),
                peptide: &self.peptide,
                proforma: &self.proforma,
                decoy,
                decoy_pair_id: pair,
                charge,
                precursor_mz: 500.25 + charge as f64,
                neutral_mass: 998.5,
                // The fixture value a context-free prediction writes: an index,
                // not a duration.
                rt: 1.559414,
                irt: None,
                mobility: 0.85,
                peaks,
            }
        }
    }

    fn peak(
        ion: &'static str,
        ordinal: i64,
        charge: i64,
        mz: f64,
        intensity: f64,
    ) -> Peak<'static> {
        Peak {
            mz,
            intensity,
            ion,
            ordinal,
            charge,
        }
    }

    fn peaks(n: usize) -> Vec<Peak<'static>> {
        (1..=n)
            .map(|i| peak("y", i as i64, 1, 200.0 + i as f64, 1.0 / i as f64))
            .collect()
    }

    /// The stats a stream that handed over `precursors` rows would return, which
    /// is the only field [`PredictedLibraryHandle::into_library`] reads.
    fn stats(precursors: usize) -> LibraryStats {
        LibraryStats {
            precursors,
            ..LibraryStats::default()
        }
    }

    /// Drive the sink with typed provenance from one real tiny prediction.
    ///
    /// Borrows the rows rather than taking them so one set can be put through
    /// this route and the file route both.
    fn build(rows: &[SpectrumRow<'_>], decoys: DecoyPolicy) -> PredictedLibrary {
        build_with_scale(rows, decoys, None)
    }

    fn build_with_scale(
        rows: &[SpectrumRow<'_>],
        decoys: DecoyPolicy,
        scale: Option<&'static str>,
    ) -> PredictedLibrary {
        let (handle, mut collector) = sink(decoys);
        collector
            .header(&test_provenance(scale))
            .expect("header converts");
        for row in rows {
            collector.spectrum(row).expect("row converts");
        }
        collector.finish().expect("stream finishes");
        handle
            .into_library(&stats(rows.len()))
            .expect("library seals")
    }

    /// Everything a scorer can observe about a library, so two of them can be
    /// compared without reaching into the arena's private columns.
    fn projection(lib: &PredictedLibrary) -> Vec<String> {
        let mut rows: Vec<_> = lib
            .library
            .geometry()
            .rows()
            .map(|tgt| {
                let geom = lib.library.geometry();
                format!(
                    "{:?} {:?} {} {} {} {} {:?} {:?} {:?}",
                    geom.output_id(tgt),
                    geom.decoy_group(tgt),
                    geom.charge(tgt),
                    geom.precursor_mz(tgt),
                    geom.library_rt(tgt).map(|rt| rt.0).unwrap(),
                    geom.analyte(tgt)
                        .peptide
                        .known()
                        .and_then(|p| p.sequence())
                        .unwrap_or_default(),
                    geom.frag_labels(tgt),
                    geom.frag_mzs(tgt),
                    &lib.library.fragment_intensities()[geom.frag_range(tgt)],
                )
            })
            .collect();
        rows.sort_unstable();
        rows
    }

    #[test]
    fn arrival_order_does_not_change_library_contents() {
        let first = Fixture::new("AAAPEPTIDEK", "AAAPEPTIDEK");
        let second = Fixture::new("MMMPEPTIDER", "MMMPEPTIDER");

        let ordered = build(
            &[
                first.row(2, false, None, peaks(3)),
                first.row(3, false, None, peaks(2)),
                second.row(2, false, None, peaks(4)),
            ],
            DecoyPolicy::Never,
        );
        let shuffled = build(
            &[
                second.row(2, false, None, peaks(4)),
                first.row(3, false, None, peaks(2)),
                first.row(2, false, None, peaks(3)),
            ],
            DecoyPolicy::Never,
        );

        assert_eq!(projection(&ordered), projection(&shuffled));
    }

    #[test]
    fn a_target_and_its_decoy_that_share_a_pair_id_land_in_one_competition_group() {
        let target = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let decoy = Fixture::new("PDITPEEK", "PDITPEEK");
        let lib = build(
            &[
                target.row(2, false, Some(7), peaks(3)),
                decoy.row(2, true, Some(7), peaks(3)),
            ],
            DecoyPolicy::Never,
        );

        let geom = lib.library.geometry();
        let groups: Vec<String> = geom
            .rows()
            .map(|tgt| format!("{:?}", geom.decoy_group(tgt)))
            .collect();
        let ids: Vec<String> = geom
            .rows()
            .map(|tgt| format!("{:?}", geom.output_id(tgt)))
            .collect();
        assert_eq!(groups[0], groups[1]);
        assert_ne!(
            groups[0], ids[0],
            "a declared pair is a group of its own, not the target's id",
        );
    }

    /// Only what `decoy_group` reports, which is all this can see: whether the
    /// label column was stored or dropped is not observable from outside
    /// timsquery, and both states give a groupless row its own id back.
    #[test]
    fn a_row_that_declares_no_pair_reports_its_own_id_as_its_group() {
        let first = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let second = Fixture::new("PEPTIDER", "PEPTIDER");
        let lib = build(
            &[
                first.row(2, false, None, peaks(2)),
                second.row(2, false, None, peaks(2)),
            ],
            DecoyPolicy::Never,
        );

        let geom = lib.library.geometry();
        for tgt in geom.rows() {
            assert_eq!(
                format!("{:?}", geom.decoy_group(tgt)),
                format!("{:?}", geom.output_id(tgt)),
                "a row that competes alone reports its own id as its group",
            );
        }
    }

    #[test]
    fn forcing_derived_decoys_drops_a_shipped_decoy_and_leaves_its_intensities_out_of_the_sidecar()
    {
        let target = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let decoy = Fixture::new("PDITPEEK", "PDITPEEK");
        let lib = build(
            &[
                target.row(2, false, Some(1), peaks(3)),
                decoy.row(2, true, Some(1), peaks(5)),
            ],
            DecoyPolicy::Force,
        );

        assert_eq!(lib.library.geometry().n_rows(), 1);
        assert_eq!(lib.library.fragment_intensities().len(), 3);
        assert!(matches!(
            lib.library.geometry().capabilities().decoys,
            DecoyStrategy::MassShift { .. }
        ));
    }

    #[test]
    fn spectrum_pushes_into_columns_before_finish() {
        let target = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let decoy = Fixture::new("PDITPEEK", "PDITPEEK");
        let (handle, mut collector) = sink(DecoyPolicy::Force);
        collector.record_test_provenance();
        collector
            .spectrum(&target.row(2, false, Some(1), peaks(3)))
            .unwrap();
        collector
            .spectrum(&decoy.row(2, true, Some(1), peaks(5)))
            .unwrap();
        let SinkState::Building { arena, .. } = &collector.state else {
            panic!("header should start building");
        };
        assert_eq!(arena.received, 2);
        assert_eq!(arena.geom.n_rows(), 1);
        assert_eq!(arena.frag_intens.len(), 3);
        collector.finish().unwrap();
        assert_eq!(
            handle
                .into_library(&stats(2))
                .unwrap()
                .library
                .geometry()
                .n_rows(),
            1
        );
    }

    #[test]
    fn the_intensity_sidecar_is_as_long_as_the_fragment_label_arena() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let lib = build(
            &[
                fixture.row(2, false, None, peaks(6)),
                fixture.row(3, false, None, peaks(4)),
            ],
            DecoyPolicy::Never,
        );

        assert_eq!(lib.library.fragment_intensities().len(), 6 + 4);
    }

    #[test]
    fn a_predicted_library_keeps_sequence_features_for_sequences_that_parse() {
        let fixture = Fixture::new("PEPC[UNIMOD:4]IDEK", "PEPCIDEK");
        let lib = build(&[fixture.row(2, false, None, peaks(3))], DecoyPolicy::Never);

        let tgt = lib.library.geometry().rows().next().unwrap();
        assert_eq!(
            lib.library
                .geometry()
                .analyte(tgt)
                .peptide
                .known()
                .and_then(|p| p.sequence())
                .unwrap_or_default(),
            "PEPC[UNIMOD:4]IDEK"
        );
        assert_eq!(
            lib.library
                .geometry()
                .analyte(tgt)
                .peptide
                .known()
                .map_or("", |p| p.residues),
            "PEPCIDEK"
        );
        assert!(lib.library.all_sequence_counts_enabled());
    }

    #[test]
    fn a_predicted_peak_gets_the_label_the_file_route_parses_from_the_same_annotation() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let lib = build(
            &[fixture.row(
                2,
                false,
                None,
                vec![
                    peak("y", 7, 1, 800.4, 1.0),
                    peak("y", 7, 2, 400.7, 0.5),
                    peak("b", 3, 1, 324.1, 0.25),
                ],
            )],
            DecoyPolicy::Never,
        );

        let tgt = lib.library.geometry().rows().next().unwrap();
        assert_eq!(
            lib.library.geometry().frag_labels(tgt),
            [
                IonAnnot::try_from("y7").unwrap(),
                IonAnnot::try_from("y7^2").unwrap(),
                IonAnnot::try_from("b3").unwrap(),
            ],
        );
    }

    #[test]
    fn declared_scale_survives_prediction_and_file_loading_only_for_normalized_rows() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        for scale in [None, Some("reference anchors")] {
            let declared_scale = scale.unwrap_or(TEST_PROVENANCE.retention.normalized.scale);
            for irt in [None, Some(37.75)] {
                let rows = [SpectrumRow {
                    irt,
                    rt: 2.0,
                    ..fixture.row(2, false, None, peaks(2))
                }];
                let predicted = build_with_scale(&rows, DecoyPolicy::Never, scale);
                let TargetTable::Mzpaf {
                    geom: from_file, ..
                } = via_file_table_with_scale(&rows, Some(declared_scale))
                else {
                    panic!("expected mzpaf")
                };
                let expected = if irt.is_some() {
                    timsquery::RtAxis::Seconds
                } else {
                    timsquery::RtAxis::NormalizedIndex {
                        scale: Some(declared_scale.to_owned()),
                    }
                };
                let sunk = predicted.library.geometry();
                assert_eq!(sunk.rt_axis(), &expected);
                assert_eq!(from_file.rt_axis(), &expected);
                assert_eq!(
                    sunk.library_rt(sunk.rows().next().unwrap()).map(|rt| rt.0),
                    from_file
                        .library_rt(from_file.rows().next().unwrap())
                        .map(|rt| rt.0)
                );
            }
        }
    }

    #[test]
    fn retention_is_stored_on_the_scale_the_file_route_reads_the_same_library_back_on() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let lib = build(&[fixture.row(2, false, None, peaks(2))], DecoyPolicy::Never);

        let tgt = lib.library.geometry().rows().next().unwrap();
        assert!(
            (lib.library
                .geometry()
                .library_rt(tgt)
                .map(|rt| rt.0)
                .unwrap()
                - 1.559414)
                .abs()
                < 1e-3
        );
    }

    #[test]
    fn a_row_predicted_against_a_gradient_stores_that_gradient_time_and_not_its_index() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        // A chromatography context puts gradient minutes in `rt` and the index
        // in `irt`, and the mzSpecLib reader reads the gradient time back, so
        // the index is the number this route must not store.
        let contextual = SpectrumRow {
            rt: 31.5,
            irt: Some(37.75),
            ..fixture.row(2, false, None, peaks(2))
        };
        let lib = build(&[contextual], DecoyPolicy::Never);

        let tgt = lib.library.geometry().rows().next().unwrap();
        assert!(
            (lib.library
                .geometry()
                .library_rt(tgt)
                .map(|rt| rt.0)
                .unwrap()
                - 1890.0)
                .abs()
                < 1e-3
        );
    }

    #[test]
    fn a_stream_that_never_finished_hands_over_no_library() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (handle, mut collector) = sink(DecoyPolicy::Never);
        collector.record_test_provenance();
        collector
            .spectrum(&fixture.row(2, false, None, peaks(3)))
            .expect("row converts");

        let Err(error) = handle.into_library(&stats(1)) else {
            panic!("a library that was never finished is not a library");
        };
        assert!(
            format!("{error}").contains("did not finish"),
            "unexpected error: {error}",
        );
    }

    /// The header `MzSpecLibSink::header` would write, minus the provenance
    /// pairs, written by hand.
    ///
    /// `header` takes a `LibraryProvenance`, which is `#[non_exhaustive]` with no
    /// constructor, so no test outside msspeculator can call it. Its `spectrum`
    /// does not depend on it having run, so every per-row value below still comes
    /// out of the real writer; what is reproduced here is only the framing the
    /// reader needs to parse one -- the format version, and the two attribute
    /// sets that carry `spectrum origin type`, which is where `is_decoy` reads
    /// from. Only the normalized-scale provenance is included when requested.
    fn write_mzspeclib_header(out: &mut impl std::io::Write, scale: Option<&str>) {
        for line in [
            "<mzSpecLib>",
            "MS:1003186|library format version=1.0",
            "MS:1003188|library name=both_routes",
            "<AttributeSet Spectrum=all>",
            "MS:1000511|ms level=2",
            "MS:1003072|spectrum origin type=MS:1003074|predicted spectrum",
            "MS:1003065|spectrum aggregation type=MS:1003074|predicted spectrum",
            "<AttributeSet Spectrum=Decoy>",
            "MS:1003072|spectrum origin type=MS:1003195|unnatural peptidoform decoy spectrum",
        ] {
            writeln!(out, "{line}").expect("header writes");
            if line == "MS:1003188|library name=both_routes"
                && let Some(scale) = scale
            {
                writeln!(
                    out,
                    "[1]MS:1003275|other attribute name=msspeculator:retention.normalized.scale"
                )
                .unwrap();
                writeln!(out, "[1]MS:1003276|other attribute value={scale}").unwrap();
            }
        }
    }

    /// The other route: msspeculator's own mzSpecLib writer to a file, then this
    /// project's reader back off it, which is what a `build-library` followed by
    /// a `search` does.
    fn via_file_table(rows: &[SpectrumRow<'_>]) -> TargetTable {
        via_file_table_with_scale(rows, Some(TEST_PROVENANCE.retention.normalized.scale))
    }

    fn via_file_table_with_scale(rows: &[SpectrumRow<'_>], scale: Option<&str>) -> TargetTable {
        let dir = tempfile::tempdir().expect("temp dir");
        // The name the sniffer dispatches on: it takes `.mzspeclib.` anywhere in
        // the file name, and reads plain text for anything not ending `.gz`.
        let path = dir.path().join("both_routes.mzspeclib.txt");
        {
            let mut file = std::fs::File::create(&path).expect("library file");
            write_mzspeclib_header(&mut file, scale);
            // The sink appends to the same handle at the cursor the header left,
            // and closes it when this scope drops it.
            let mut writer = msspeculator_inference::mzspeclib::MzSpecLibSink::new(file, &path);
            for row in rows {
                writer.spectrum(row).expect("row writes");
            }
            writer.finish().expect("file finishes");
        }
        timsquery::serde::read_targets_with(
            &path,
            timsseek::LoadPolicy {
                decoys: DecoyPolicy::Never,
                ..Default::default()
            },
        )
        .expect("the file this project writes reads back")
    }

    fn via_file(rows: &[SpectrumRow<'_>]) -> TargetColumns<IonAnnot> {
        let TargetTable::Mzpaf { geom, .. } = via_file_table(rows) else {
            panic!("an mzSpecLib library is mzpaf-labelled");
        };
        geom
    }

    #[test]
    fn modified_prediction_and_reload_share_composition_envelopes() {
        use timsseek::data_sources::reference_library::ExpectedIntensity;
        use timsseek::fragment_mass::isotope_plan::IsotopeMethod;
        let peptide = Fixture::new("PEPTC[UNIMOD:4]IDEK", "PEPTCIDEK");
        let rows = [peptide.row(2, false, None, peaks(4))];
        let predicted = build(&rows, DecoyPolicy::Never);
        let reloaded = ReferenceLibrary::try_from(via_file_table(&rows)).unwrap();
        for lib in [&predicted.library, &reloaded] {
            assert_eq!(
                lib.scoring_plan().isotopes().method,
                IsotopeMethod::CompositionCs
            );
            let envelope: Vec<_> = lib
                .iter()
                .next()
                .unwrap()
                .expected_precursor_envelope()
                .iter()
                .map(|(_, value)| *value)
                .collect();
            assert_eq!(envelope, timsseek::isotopes::peptide_isotopes(45, 1));
        }
    }

    /// This row's fragments as label to m/z, which is how the two routes can be
    /// compared at all: the writer sorts peaks by m/z and the sink keeps
    /// `(position, ion type)`, so the sequences differ where the mapping does not.
    fn fragments(geom: &TargetColumns<IonAnnot>, at: RowIdx) -> BTreeMap<String, f64> {
        geom.frag_labels(at)
            .iter()
            .map(|label| label.to_string())
            .zip(geom.frag_mzs(at).iter().copied())
            .collect()
    }

    /// Find the row carrying `id`; persistence is free to choose storage order.
    fn row_named(geom: &TargetColumns<IonAnnot>, id: &str) -> RowIdx {
        geom.rows()
            .find(|row| geom.output_id(*row).to_string() == id)
            .unwrap_or_else(|| panic!("no row named {id}"))
    }

    /// The premise the rest of this module rests on: one prediction through both
    /// routes agrees on every value a scorer reads.
    ///
    /// Every such quantity rather than a sample of them, because the two routes
    /// spell each one independently and any of them could drift. Two charge
    /// states of one peptide, which msspeculator gives one `decoy_pair_id`, so
    /// the group column is stored rather than derived and the comparison of it
    /// means something.
    ///
    /// Fragment m/z is compared to 1e-5 rather than exactly: the writer prints
    /// six decimal places, so the file route reads back a rounded value of what
    /// the sink keeps whole.
    #[test]
    fn a_prediction_lands_on_the_same_values_through_this_sink_as_through_a_file() {
        let peptide = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let rows = [
            peptide.row(
                2,
                false,
                Some(4),
                vec![
                    peak("y", 7, 1, 800.4, 1.0),
                    peak("b", 3, 1, 324.1, 0.25),
                    peak("y", 7, 2, 400.7, 0.5),
                ],
            ),
            peptide.row(3, false, Some(4), peaks(4)),
        ];

        let predicted = build(&rows, DecoyPolicy::Never);
        let sunk = predicted.library.geometry();
        let from_file = via_file(&rows);
        assert_eq!(sunk.n_rows(), from_file.n_rows());

        for row in sunk.rows() {
            let id = sunk.output_id(row).to_string();
            let mirror = row_named(&from_file, &id);
            assert_eq!(
                sunk.analyte(row).to_owned(),
                from_file.analyte(mirror).to_owned(),
                "chemistry of {id}"
            );
            assert_eq!(
                sunk.entry_name(row),
                from_file.entry_name(mirror),
                "label of {id}"
            );

            assert_eq!(sunk.charge(row), from_file.charge(mirror), "charge of {id}");
            assert_eq!(
                sunk.decoy_group(row).to_string(),
                from_file.decoy_group(mirror).to_string(),
                "competition group of {id}",
            );
            assert_ne!(
                sunk.decoy_group(row).to_string(),
                id,
                "a declared pair is a group of its own, so this comparison is not the id twice",
            );
            assert!(
                (sunk.precursor_mz(row) - from_file.precursor_mz(mirror)).abs() < 1e-6,
                "precursor m/z of {id}",
            );
            assert!(
                (sunk.library_rt(row).map(|rt| rt.0).unwrap()
                    - from_file.library_rt(mirror).map(|rt| rt.0).unwrap())
                .abs()
                    < 1e-3,
                "retention of {id}: {} against {}",
                sunk.library_rt(row).map(|rt| rt.0).unwrap(),
                from_file.library_rt(mirror).map(|rt| rt.0).unwrap(),
            );
            assert!(
                (sunk.mobility(row) - from_file.mobility(mirror)).abs() < 1e-6,
                "mobility of {id}",
            );

            assert_eq!(sunk.rt_axis(), from_file.rt_axis());
            let (sunk_frags, file_frags) = (fragments(sunk, row), fragments(&from_file, mirror));
            assert_eq!(
                sunk_frags.keys().collect::<Vec<_>>(),
                file_frags.keys().collect::<Vec<_>>(),
                "fragment labels of {id}",
            );
            for (label, mz) in &sunk_frags {
                assert!(
                    (mz - file_frags[label]).abs() < 1e-5,
                    "m/z of {label} on {id}: {mz} against {}",
                    file_frags[label],
                );
            }
        }
    }

    /// A shipped decoy is the row where the two routes are most easily made to
    /// disagree, and it agrees on everything a scorer reads.
    ///
    /// msspeculator's writer claims `MS:1003195|unnatural peptidoform decoy
    /// spectrum` for a decoy, which replaces the `MS:1003074|predicted spectrum`
    /// its `Spectrum=all` set gives everything else. Only the latter says the
    /// peak m/z values are calculated, and a decoy with no fragments is an FDR
    /// estimated against nothing, so the decoy declaration must cost the entry
    /// none of its peaks.
    #[test]
    fn a_shipped_decoy_keeps_its_fragments_through_this_sink_and_through_a_file() {
        let target = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let decoy = Fixture::new("PDITPEEK", "PDITPEEK");
        let rows = [
            target.row(2, false, Some(4), peaks(3)),
            decoy.row(2, true, Some(4), peaks(4)),
        ];

        let predicted = build(&rows, DecoyPolicy::Never);
        let sunk = predicted.library.geometry();
        let from_file = via_file(&rows);

        let sunk_decoy = row_named(sunk, "PDITPEEK/2");
        let file_decoy = row_named(&from_file, "PDITPEEK/2");
        assert!(sunk.is_decoy(sunk_decoy) && from_file.is_decoy(file_decoy));

        assert_eq!(
            sunk.decoy_group(sunk_decoy).to_string(),
            from_file.decoy_group(file_decoy).to_string(),
        );
        assert_eq!(sunk.charge(sunk_decoy), from_file.charge(file_decoy));

        let (sunk_frags, file_frags) = (
            fragments(sunk, sunk_decoy),
            fragments(&from_file, file_decoy),
        );
        assert_eq!(sunk_frags.len(), 4);
        assert_eq!(
            sunk_frags.keys().collect::<Vec<_>>(),
            file_frags.keys().collect::<Vec<_>>(),
        );
        for (label, mz) in &sunk_frags {
            assert!(
                (mz - file_frags[label]).abs() < 1e-5,
                "m/z of {label} on the decoy: {mz} against {}",
                file_frags[label],
            );
        }
        // The target alongside it keeps its own count, so the decoy's peaks were
        // not read off the wrong entry.
        assert_eq!(
            fragments(&from_file, row_named(&from_file, "PEPTIDEK/2")).len(),
            3,
        );
    }

    /// The shape a panicked inference worker leaves behind: `finish` ran and
    /// published, so the rows are there, but fewer of them than the stream
    /// counted. Sealing that would put an FDR on a fraction of the proteome.
    #[test]
    fn a_stream_that_published_fewer_rows_than_it_counted_hands_over_no_library() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (handle, mut collector) = sink(DecoyPolicy::Never);
        collector.record_test_provenance();
        for charge in 2..=3 {
            collector
                .spectrum(&fixture.row(charge, false, None, peaks(3)))
                .expect("row converts");
        }
        collector.finish().expect("stream finishes");

        let Err(error) = handle.into_library(&stats(5)) else {
            panic!("a prefix of a library is not a library");
        };
        let error = format!("{error}");
        for number in ["5", "2"] {
            assert!(
                error.contains(number),
                "the error names both counts, got: {error}",
            );
        }
    }

    #[test]
    fn a_prediction_that_produced_no_precursors_at_all_is_an_error_and_not_an_empty_library() {
        let (handle, mut collector) = sink(DecoyPolicy::Never);
        collector.record_test_provenance();
        collector.finish().expect("stream finishes");

        let Err(error) = handle.into_library(&stats(0)) else {
            panic!("a library with nothing in it searches to nothing and reports no reason");
        };
        assert!(
            format!("{error}").contains("charge range"),
            "the error says what to check, got: {error}",
        );
    }

    #[test]
    fn a_stream_that_published_rows_without_a_header_hands_over_no_library() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (handle, mut collector) = sink(DecoyPolicy::Never);
        let error = collector
            .spectrum(&fixture.row(2, false, None, peaks(3)))
            .expect_err("spectrum before header must fail");
        assert!(format!("{error}").contains("before header"));
        assert!(format!("{}", collector.finish().unwrap_err()).contains("before header"));

        let Err(error) = handle.into_library(&stats(1)) else {
            panic!("a library whose provenance is unknown is not one this can record");
        };
        assert!(
            format!("{error}").contains("did not finish"),
            "unexpected error: {error}",
        );
    }

    #[test]
    fn callbacks_after_finish_are_rejected() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (_handle, mut collector) = sink(DecoyPolicy::Never);
        collector.record_test_provenance();
        collector.finish().unwrap();
        assert!(format!("{}", collector.finish().unwrap_err()).contains("more than once"));
        assert!(
            format!(
                "{}",
                collector
                    .spectrum(&fixture.row(2, false, None, peaks(3)))
                    .unwrap_err()
            )
            .contains("after finish")
        );
    }

    #[test]
    fn a_peak_whose_series_is_not_an_mzpaf_letter_fails_the_row_rather_than_losing_the_peak() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (handle, mut collector) = sink(DecoyPolicy::Never);
        collector.record_test_provenance();
        let error = collector
            .spectrum(&fixture.row(2, false, None, vec![peak("Q", 3, 1, 300.0, 1.0)]))
            .expect_err("Q is not a series letter");
        assert!(
            format!("{error:#}").contains("PEPTIDEK/2"),
            "unexpected error: {error:#}",
        );

        // Nothing partial survived the rejected row: the peak it could not label
        // took the whole row with it, so there is no library left to seal.
        collector.finish().expect("stream finishes");
        let Err(error) = handle.into_library(&stats(0)) else {
            panic!("a rejected row left a library behind");
        };
        assert!(
            format!("{error}").contains("no precursors"),
            "unexpected error: {error}",
        );
    }
}
