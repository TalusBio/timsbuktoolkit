//! A library predicted straight into the arena timsseek scores, with no file in
//! between.
//!
//! msspeculator hands every precursor to a [`LibrarySink`]; this one spells the
//! rows into [`TargetColumns`] instead of into mzSpecLib or DIA-NN text. The
//! result has to agree with loading the same library off disk on every value a
//! scorer reads, so every spelling decision here mirrors the mzSpecLib reader in
//! `timsquery::serde::mzspeclib_io` -- the two routes agree on the numbers they
//! store, not on the strings they pass through. Row order is the exception: the
//! file route stores prediction-arrival order, this one sorts (see
//! [`build_arena`]).
//!
//! The sink is moved onto msspeculator's writer thread and never handed back, so
//! the finished rows come out through shared state: [`sink`] returns a handle
//! that stays with the caller and the sink that goes into `stream_library`.

// Reached only by its own tests until `search --fasta` calls it; the allow goes
// with that command.
#![allow(dead_code)]

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
    TargetColumns,
};
use timsquery::serde::{
    TargetReadingError,
    TargetTable,
};
use timsseek::DecoyPolicy;
use timsseek::data_sources::reference_library::ReferenceLibrary;

use crate::errors::CliError;

/// The scale the file route puts retention on: sixty times what msspeculator
/// reported, whatever that number means.
///
/// With no chromatography context it is a dimensionless PROCAL-anchored index,
/// which the mzSpecLib writer has to declare in `minute` because the vocabulary
/// has no unit for an index, and the mzSpecLib reader multiplies by sixty
/// regardless -- which recovers the declared minutes only because the pinned
/// mzannotate leaves a `minute`-declared value alone (a release divides it by
/// sixty). Scoring only requires `rt_seconds` to be ordered, so scaling one
/// route and not the other gives two libraries 60x apart that each work alone
/// and disagree the moment a run holds one of each --
/// `check_rt_scale_compatibility` then reports a pair with no overlapping RT
/// range. That the unit is wrong for the context-free case is
/// https://github.com/TalusBio/timsbuktoolkit/issues/115, which is a different
/// question from the two routes agreeing.
const SECONDS_PER_MINUTE: f32 = 60.0;

/// A library that was predicted rather than read, plus what produced it.
pub(crate) struct PredictedLibrary {
    pub library: ReferenceLibrary,
    /// msspeculator's provenance, as the JSON its own sidecar carries, for a
    /// caller recording what a run used.
    pub provenance: serde_json::Value,
}

/// The caller's half of [`sink`]: what the prediction produced, once it has.
pub(crate) struct PredictedLibraryHandle {
    shared: Arc<Mutex<Handoff>>,
}

/// The [`LibrarySink`] half of [`sink`], for `stream_library` to consume.
pub(crate) struct PredictedLibrarySink {
    shared: Arc<Mutex<Handoff>>,
    /// Every row so far, owned. A [`SpectrumRow`] borrows from the prediction it
    /// came out of, so nothing kept past `spectrum` can borrow; and the rows are
    /// held rather than pushed as they arrive because [`build_arena`] needs the
    /// whole set to order it.
    rows: Vec<PredictedRow>,
}

/// What crosses from the writer thread back to the caller.
///
/// `rows` is `Some` only once `finish` ran, which is not the same as the stream
/// having succeeded: msspeculator's `run_library` joins its inference workers
/// before its writer, so a worker that panics drops the result channel, the
/// writer's loop over it ends cleanly, and `finish` publishes whatever arrived
/// before the panic. `stream_library` still returns the error, which is why
/// [`PredictedLibraryHandle::into_library`] takes the [`LibraryStats`] only a
/// successful stream produces and checks its count against what arrived.
#[derive(Default)]
struct Handoff {
    provenance: Option<serde_json::Value>,
    rows: Option<Vec<PredictedRow>>,
}

/// Build a sink and the handle that collects from it.
pub(crate) fn sink() -> (PredictedLibraryHandle, PredictedLibrarySink) {
    let shared = Arc::new(Mutex::new(Handoff::default()));
    (
        PredictedLibraryHandle {
            shared: Arc::clone(&shared),
        },
        PredictedLibrarySink {
            shared,
            rows: Vec::new(),
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
    ///
    /// The decoy policy arrives here rather than at [`sink`] because it decides
    /// which rows reach the arena, and no row reaches it until every row has
    /// arrived.
    pub(crate) fn into_library(
        self,
        stats: &LibraryStats,
        decoys: DecoyPolicy,
    ) -> Result<PredictedLibrary, CliError> {
        let handoff = std::mem::take(&mut *self.shared.lock().expect("handoff mutex poisoned"));
        let Some(rows) = handoff.rows else {
            return Err(CliError::LibraryBuild {
                source: "the prediction stream handed over no rows, so it did not finish"
                    .to_string(),
            });
        };
        if rows.len() != stats.precursors {
            return Err(CliError::LibraryBuild {
                source: format!(
                    "the prediction reported {} precursors but the sink received {}, so what it \
                     handed over is a prefix of the library that was asked for",
                    stats.precursors,
                    rows.len(),
                ),
            });
        }
        let Some(provenance) = handoff.provenance else {
            return Err(CliError::LibraryBuild {
                source: "the prediction stream handed over rows without a header, so its \
                         provenance is unknown"
                    .to_string(),
            });
        };
        let arena = build_arena(rows, decoys).map_err(|e| CliError::LibraryBuild {
            source: format!("assembling the predicted library: {e:?}"),
        })?;
        let library =
            ReferenceLibrary::from_sealed_arena(arena).map_err(|e| CliError::LibraryBuild {
                source: format!("finalizing the predicted library: {e:?}"),
            })?;
        Ok(PredictedLibrary {
            library,
            provenance,
        })
    }
}

impl PredictedLibrarySink {
    fn record_provenance(&mut self, provenance: serde_json::Value) {
        self.shared
            .lock()
            .expect("handoff mutex poisoned")
            .provenance = Some(provenance);
    }

    /// Stand in for [`LibrarySink::header`], which no test can call:
    /// `LibraryProvenance` is `#[non_exhaustive]` and has no constructor, so only
    /// msspeculator can build one.
    #[cfg(test)]
    fn record_test_provenance(&mut self) {
        self.record_provenance(serde_json::json!({ "generator": { "tool": "test" } }));
    }
}

impl LibrarySink for PredictedLibrarySink {
    fn header(&mut self, provenance: &LibraryProvenance) -> Result<()> {
        // Kept as the JSON `to_json` builds, unflattened: this is a record of
        // what produced the library, and nothing here reads it.
        self.record_provenance(provenance.to_json());
        Ok(())
    }

    fn spectrum(&mut self, row: &SpectrumRow<'_>) -> Result<()> {
        self.rows.push(PredictedRow::from_spectrum(row)?);
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.shared.lock().expect("handoff mutex poisoned").rows =
            Some(std::mem::take(&mut self.rows));
        Ok(())
    }
}

/// One prediction, owned and already spelled the way the arena stores it.
#[derive(Debug)]
struct PredictedRow {
    /// `{proforma}/{charge}`, which is the name the mzSpecLib writer gives the
    /// same precursor under `MS:1003061|library spectrum name`, so a row keeps
    /// one name across both routes. Pushed as the arena's `seq_mod` too, which
    /// wants the same string: a proforma peptidoform ion with its charge.
    id: String,
    /// The source peptide's pair id, `None` when decoys are off. All modified
    /// forms and all charges of one peptide share it, so a group is a peptide
    /// rather than a target/decoy couple.
    group: Option<String>,
    seq_strip: String,
    precursor_mz: f64,
    charge: u8,
    rt_seconds: f32,
    mobility: f32,
    frags: Vec<(IonAnnot, f64)>,
    /// Parallel to `frags`, kept apart because the arena stores the two in
    /// different places: labels and m/z in the columns, intensities in a sidecar.
    intensities: Vec<f32>,
    is_decoy: bool,
}

impl PredictedRow {
    fn from_spectrum(row: &SpectrumRow<'_>) -> Result<Self> {
        // The charge suffix is part of the stored sequence, because the file
        // route stores mzcore's `PeptidoformIon` spelling and that renders its
        // charge carriers. The modification spelling still differs -- ours is
        // `C[UNIMOD:4]`, mzcore's is `C[U:Carbamidomethyl]` -- and both parse,
        // so the parse gate reaches the same verdict from either.
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
            seq_strip: row.stripped.to_string(),
            id,
            group: row.decoy_pair_id.map(|pair| pair.to_string()),
            precursor_mz: row.precursor_mz,
            charge,
            // `rt` is the quantity the mzSpecLib writer puts under the term the
            // file route reads back as `start_time`: gradient minutes with a
            // chromatography context, the normalized index without one. `irt`
            // carries the index alongside a gradient time, under a term the
            // reader keeps as a scan param and never scores on.
            rt_seconds: row.rt * SECONDS_PER_MINUTE,
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

/// Push every kept row into a sealed arena, with the intensity sidecar beside it.
///
/// Rows go in sorted rather than in arrival order. msspeculator calls `spectrum`
/// as predictions come back from a pool of workers, so the same FASTA predicted
/// twice arrives in a different order; `RowIdx` is then the canonicalization key
/// for rescoring, where `canonicalize_and_shuffle` sorts by
/// `(identity.row, precursor_charge)` and cuts folds positionally. Arrival order
/// would put a peptide in a different fold on every run and give one library two
/// sets of q-values. The buffer this drains exists for that and not for speed.
///
/// `{proforma}/{charge}` is a total order over the rows and is also the source
/// id, so two rows can only tie by being one duplicated id, which the seal
/// rejects.
fn build_arena(
    mut rows: Vec<PredictedRow>,
    decoys: DecoyPolicy,
) -> Result<TargetTable, TargetReadingError> {
    // A zero-row arena seals: ids and groups are both vacuously consistent, the
    // parse gate has nothing to reject, and the result searches to zero results
    // without ever reporting why. msspeculator only refuses an empty digest, so
    // a FASTA whose peptides all fall outside the windows gets this far.
    if rows.is_empty() {
        return Err(TargetReadingError::SpeclibParse(
            "the prediction produced no precursors; check the charge range and the length range \
             against the peptides the FASTA actually digests to"
                .to_string(),
        ));
    }
    rows.sort_unstable_by(|a, b| a.id.cmp(&b.id));

    let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
    let mut frag_intens: Vec<f32> = Vec::new();

    for row in &rows {
        // Dropped before the push, which is what makes `Force` mean anything:
        // neither the row nor its intensities reach the arena, so the seal sees
        // an all-target library and resolves to the derived decoys the policy
        // asked for.
        if !decoys.accepts(row.is_decoy) {
            continue;
        }
        // Ids and group labels are one namespace, because a row that names no
        // group falls back to its own id. A peptide with no pair id is therefore
        // its own singleton group, which is what the seal drops the column for,
        // and text on both sides so no row can mix an id shape with a group's.
        let group = row.group.clone().unwrap_or_else(|| row.id.clone());
        frag_intens.extend_from_slice(&row.intensities);
        geom.push_row(Row {
            precursor_mz: row.precursor_mz,
            charge: row.charge,
            rt_seconds: row.rt_seconds,
            mobility: row.mobility,
            frags: &row.frags,
            seq_strip: &row.seq_strip,
            seq_mod: &row.id,
            is_decoy: row.is_decoy,
            id: Some(row.id.clone().into()),
            decoy_group: Some(group.into()),
            // No structured mods: the token registry is timsquery-internal, so
            // no reader outside it can name a modification, and the modified
            // sequence carries the same information.
            ..Default::default()
        });
    }

    // The same empty arena the guard above refuses, reached the other way round:
    // every row that arrived was a shipped decoy and the policy dropped all of
    // them, leaving nothing to derive decoys against.
    if geom.n_rows() == 0 {
        return Err(TargetReadingError::SpeclibParse(format!(
            "all {} predicted rows were shipped decoys, which {decoys:?} drops, so the library \
             holds no targets",
            rows.len(),
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

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use msspeculator_core::peptide::Peptide;
    use msspeculator_inference::{
        ProteinGroup,
        Residues,
    };
    use timsquery::models::{
        DecoyStrategy,
        RowIdx,
        SeqFeatureState,
    };

    use super::*;

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

    /// Drive the sink the way `stream_library` does, minus the header, whose
    /// `LibraryProvenance` cannot be built here.
    ///
    /// Borrows the rows rather than taking them so one set can be put through
    /// this route and the file route both.
    fn build(rows: &[SpectrumRow<'_>], decoys: DecoyPolicy) -> PredictedLibrary {
        let (handle, mut collector) = sink();
        collector.record_test_provenance();
        for row in rows {
            collector.spectrum(row).expect("row converts");
        }
        collector.finish().expect("stream finishes");
        handle
            .into_library(&stats(rows.len()), decoys)
            .expect("library seals")
    }

    /// Everything a scorer can observe about a library, so two of them can be
    /// compared without reaching into the arena's private columns.
    fn projection(lib: &PredictedLibrary) -> Vec<String> {
        lib.library
            .geom
            .rows()
            .map(|tgt| {
                let geom = &lib.library.geom;
                format!(
                    "{:?} {:?} {} {} {} {} {:?} {:?} {:?}",
                    geom.output_id(tgt),
                    geom.decoy_group(tgt),
                    geom.charge(tgt),
                    geom.precursor_mz(tgt),
                    geom.rt_seconds(tgt),
                    geom.seq_mod(tgt),
                    geom.frag_labels(tgt),
                    geom.frag_mzs(tgt),
                    &lib.library.frag_intens[geom.frag_range(tgt)],
                )
            })
            .collect()
    }

    #[test]
    fn rows_handed_over_out_of_order_produce_the_same_library_as_rows_handed_over_in_order() {
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
        assert_eq!(ordered.library.frag_intens, shuffled.library.frag_intens);
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

        let geom = &lib.library.geom;
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

        let geom = &lib.library.geom;
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

        assert_eq!(lib.library.geom.n_rows(), 1);
        assert_eq!(lib.library.frag_intens.len(), 3);
        assert!(matches!(
            lib.library.geom.caps.decoys,
            DecoyStrategy::MassShift { .. }
        ));
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

        assert_eq!(lib.library.frag_intens.len(), 6 + 4);
    }

    #[test]
    fn a_predicted_library_keeps_sequence_features_for_sequences_that_parse() {
        let fixture = Fixture::new("PEPC[UNIMOD:4]IDEK", "PEPCIDEK");
        let lib = build(&[fixture.row(2, false, None, peaks(3))], DecoyPolicy::Never);

        let tgt = lib.library.geom.rows().next().unwrap();
        assert_eq!(lib.library.geom.seq_mod(tgt), "PEPC[UNIMOD:4]IDEK/2");
        assert_eq!(lib.library.geom.seq_strip(tgt), "PEPCIDEK");
        assert_eq!(
            lib.library.geom.caps.sequence_features,
            SeqFeatureState::Available,
        );
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

        let tgt = lib.library.geom.rows().next().unwrap();
        assert_eq!(
            lib.library.geom.frag_labels(tgt),
            [
                IonAnnot::try_from("y7").unwrap(),
                IonAnnot::try_from("y7^2").unwrap(),
                IonAnnot::try_from("b3").unwrap(),
            ],
        );
    }

    #[test]
    fn retention_is_stored_on_the_scale_the_file_route_reads_the_same_library_back_on() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let lib = build(&[fixture.row(2, false, None, peaks(2))], DecoyPolicy::Never);

        let tgt = lib.library.geom.rows().next().unwrap();
        assert!((lib.library.geom.rt_seconds(tgt) - 93.56484).abs() < 1e-3);
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

        let tgt = lib.library.geom.rows().next().unwrap();
        assert!((lib.library.geom.rt_seconds(tgt) - 1890.0).abs() < 1e-3);
    }

    #[test]
    fn a_stream_that_never_finished_hands_over_no_library() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (handle, mut collector) = sink();
        collector
            .spectrum(&fixture.row(2, false, None, peaks(3)))
            .expect("row converts");

        let Err(error) = handle.into_library(&stats(1), DecoyPolicy::Never) else {
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
    /// from. The provenance pairs are dropped because no per-row value comes from
    /// them.
    fn write_mzspeclib_header(out: &mut impl std::io::Write) {
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
        }
    }

    /// The other route: msspeculator's own mzSpecLib writer to a file, then this
    /// project's reader back off it, which is what a `build-library` followed by
    /// a `search` does.
    fn via_file(rows: &[SpectrumRow<'_>]) -> TargetColumns<IonAnnot> {
        let dir = tempfile::tempdir().expect("temp dir");
        // The name the sniffer dispatches on: it takes `.mzspeclib.` anywhere in
        // the file name, and reads plain text for anything not ending `.gz`.
        let path = dir.path().join("both_routes.mzspeclib.txt");
        {
            let mut file = std::fs::File::create(&path).expect("library file");
            write_mzspeclib_header(&mut file);
            // The sink appends to the same handle at the cursor the header left,
            // and closes it when this scope drops it.
            let mut writer = msspeculator_inference::mzspeclib::MzSpecLibSink::new(file, &path);
            for row in rows {
                writer.spectrum(row).expect("row writes");
            }
            writer.finish().expect("file finishes");
        }
        let TargetTable::Mzpaf { geom, .. } =
            timsquery::serde::read_targets_with(&path, DecoyPolicy::Never)
                .expect("the file this project writes reads back")
        else {
            panic!("an mzSpecLib library is mzpaf-labelled");
        };
        geom
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

    /// Find the row carrying `id`, since row order is the one thing the two
    /// routes are not expected to share (see [`build_arena`]).
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
        let sunk = &predicted.library.geom;
        let from_file = via_file(&rows);
        assert_eq!(sunk.n_rows(), from_file.n_rows());

        for row in sunk.rows() {
            let id = sunk.output_id(row).to_string();
            let mirror = row_named(&from_file, &id);

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
                (sunk.rt_seconds(row) - from_file.rt_seconds(mirror)).abs() < 1e-3,
                "retention of {id}: {} against {}",
                sunk.rt_seconds(row),
                from_file.rt_seconds(mirror),
            );
            assert!(
                (sunk.mobility(row) - from_file.mobility(mirror)).abs() < 1e-6,
                "mobility of {id}",
            );

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

    /// A decoy is the one row the two routes do not agree on, and it loses
    /// everything a scorer would match it by.
    ///
    /// msspeculator's writer claims `MS:1003195|unnatural peptidoform decoy
    /// spectrum` for a decoy, which replaces the `MS:1003074|predicted spectrum`
    /// its `Spectrum=all` set gives everything else. Only the latter is in the
    /// closure of origin types that declare their peak m/z calculated, so the
    /// reader stops treating a decoy's stated m/z as theoretical, finds neither a
    /// composition nor a mass-error suffix to recover one from, and drops every
    /// peak. The same thing happens to the committed
    /// `msspeculator_built.mzspeclib.txt`, whose two decoys load with no
    /// fragments at all.
    ///
    /// Asserted rather than left as a surprise so that the sink is not read as
    /// matching the file route on a shipped decoy: on this it is the file route
    /// that is wrong, and a fix there should fail this test.
    #[test]
    fn a_shipped_decoy_keeps_its_fragments_through_this_sink_and_loses_them_through_a_file() {
        let target = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let decoy = Fixture::new("PDITPEEK", "PDITPEEK");
        let rows = [
            target.row(2, false, Some(4), peaks(3)),
            decoy.row(2, true, Some(4), peaks(4)),
        ];

        let predicted = build(&rows, DecoyPolicy::Never);
        let sunk = &predicted.library.geom;
        let from_file = via_file(&rows);

        let sunk_decoy = row_named(sunk, "PDITPEEK/2");
        let file_decoy = row_named(&from_file, "PDITPEEK/2");
        assert!(sunk.is_decoy(sunk_decoy) && from_file.is_decoy(file_decoy));

        // Everything but the peaks still agrees, so the divergence is the origin
        // type's effect on the peak list and nothing wider.
        assert_eq!(
            sunk.decoy_group(sunk_decoy).to_string(),
            from_file.decoy_group(file_decoy).to_string(),
        );
        assert_eq!(sunk.charge(sunk_decoy), from_file.charge(file_decoy));

        assert_eq!(fragments(sunk, sunk_decoy).len(), 4);
        assert_eq!(
            fragments(&from_file, file_decoy).len(),
            0,
            "the file route drops a decoy's peaks; see this test's doc comment",
        );
        // The target alongside it keeps its own, so this is the decoy's origin
        // type and not the file.
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
        let (handle, mut collector) = sink();
        collector.record_test_provenance();
        for charge in 2..=3 {
            collector
                .spectrum(&fixture.row(charge, false, None, peaks(3)))
                .expect("row converts");
        }
        collector.finish().expect("stream finishes");

        let Err(error) = handle.into_library(&stats(5), DecoyPolicy::Never) else {
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
        let (handle, mut collector) = sink();
        collector.record_test_provenance();
        collector.finish().expect("stream finishes");

        let Err(error) = handle.into_library(&stats(0), DecoyPolicy::Never) else {
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
        let (handle, mut collector) = sink();
        collector
            .spectrum(&fixture.row(2, false, None, peaks(3)))
            .expect("row converts");
        collector.finish().expect("stream finishes");

        let Err(error) = handle.into_library(&stats(1), DecoyPolicy::Never) else {
            panic!("a library whose provenance is unknown is not one this can record");
        };
        assert!(
            format!("{error}").contains("provenance"),
            "unexpected error: {error}",
        );
    }

    #[test]
    fn a_peak_whose_series_is_not_an_mzpaf_letter_fails_the_row_rather_than_losing_the_peak() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (handle, mut collector) = sink();
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
        let Err(error) = handle.into_library(&stats(0), DecoyPolicy::Never) else {
            panic!("a rejected row left a library behind");
        };
        assert!(
            format!("{error}").contains("no precursors"),
            "unexpected error: {error}",
        );
    }
}
