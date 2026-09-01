//! A library predicted straight into the arena timsseek scores, with no file in
//! between.
//!
//! msspeculator hands every precursor to a [`LibrarySink`]; this one spells the
//! rows into [`TargetColumns`] instead of into mzSpecLib or DIA-NN text. The
//! result has to be indistinguishable from loading the same library off disk, so
//! every spelling decision here mirrors the mzSpecLib reader in
//! `timsquery::serde::mzspeclib_io` -- the two routes are compared by the numbers
//! they produce, not by the strings they pass through.
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
/// regardless. Scoring only requires `rt_seconds` to be ordered, so scaling one
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
    /// caller recording what a run used. `None` when the stream never reached
    /// its header, which is also the only way to hold this type without one:
    /// `LibraryProvenance` cannot be constructed outside msspeculator.
    pub provenance: Option<serde_json::Value>,
}

/// The caller's half of [`sink`]: what the prediction produced, once it has.
pub(crate) struct PredictedLibraryHandle {
    shared: Arc<Mutex<Handoff>>,
}

/// The [`LibrarySink`] half of [`sink`], for `stream_library` to consume.
pub(crate) struct PredictedLibrarySink {
    shared: Arc<Mutex<Handoff>>,
    /// Every row so far, owned. A [`SpectrumRow`] borrows from the prediction it
    /// came out of, so nothing kept past `spectrum` can borrow, and the arena
    /// cannot be built row by row anyway (see [`build_arena`]).
    rows: Vec<PredictedRow>,
}

/// What crosses from the writer thread back to the caller.
///
/// `rows` is `Some` only once `finish` ran, which msspeculator calls on the
/// success path alone. A stream that failed therefore leaves nothing here and
/// [`PredictedLibraryHandle::into_library`] refuses, rather than handing back the
/// prefix of a library that arrived before the error.
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
    /// The decoy policy arrives here rather than at [`sink`] because it decides
    /// which rows reach the arena, and no row reaches it until every row has
    /// arrived.
    pub(crate) fn into_library(self, decoys: DecoyPolicy) -> Result<PredictedLibrary, CliError> {
        let handoff = std::mem::take(&mut *self.shared.lock().expect("handoff mutex poisoned"));
        let Some(rows) = handoff.rows else {
            return Err(CliError::LibraryBuild {
                source: "the prediction stream handed over no rows, so it did not finish"
                    .to_string(),
            });
        };
        let arena = build_arena(rows, decoys).map_err(|e| CliError::LibraryBuild {
            source: format!("assembling the predicted library: {e:?}"),
        })?;
        // The one finalize every library goes through, file-loaded or predicted:
        // decoy reporting, the whole-library parse gate, the entry tally.
        let library =
            ReferenceLibrary::from_sealed_arena(arena).map_err(|e| CliError::LibraryBuild {
                source: format!("finalizing the predicted library: {e:?}"),
            })?;
        Ok(PredictedLibrary {
            library,
            provenance: handoff.provenance,
        })
    }
}

impl LibrarySink for PredictedLibrarySink {
    fn header(&mut self, provenance: &LibraryProvenance) -> Result<()> {
        // Kept as the JSON `to_json` builds, unflattened: this is a record of
        // what produced the library, and nothing here reads it.
        self.shared
            .lock()
            .expect("handoff mutex poisoned")
            .provenance = Some(provenance.to_json());
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
#[derive(Debug, Clone, PartialEq)]
struct PredictedRow {
    /// `{proforma}/{charge}`, which is the name the mzSpecLib writer gives the
    /// same precursor under `MS:1003061|library spectrum name`, so a row keeps
    /// one name across both routes.
    id: String,
    /// The source peptide's pair id, `None` when decoys are off. All modified
    /// forms and all charges of one peptide share it, so a group is a peptide
    /// rather than a target/decoy couple.
    group: Option<String>,
    seq_strip: String,
    seq_mod: String,
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
            seq_mod: id.clone(),
            id,
            group: row.decoy_pair_id.map(|pair| pair.to_string()),
            precursor_mz: row.precursor_mz,
            charge,
            // `rt` is the quantity the mzSpecLib writer puts under the term the
            // file route reads back: gradient minutes with a chromatography
            // context, the normalized index without one. `irt` is the second
            // copy that exists only in the first case, and the file route never
            // sees it, so neither does this.
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
            seq_mod: &row.seq_mod,
            is_decoy: row.is_decoy,
            id: Some(row.id.clone().into()),
            decoy_group: Some(group.into()),
            // No structured mods: the token registry is timsquery-internal, so
            // no reader outside it can name a modification, and the modified
            // sequence carries the same information.
            ..Default::default()
        });
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
    use msspeculator_core::peptide::Peptide;
    use msspeculator_inference::{
        ProteinGroup,
        Residues,
    };
    use timsquery::models::{
        DecoyStrategy,
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

    /// Drive the sink the way `stream_library` does, minus the header it cannot
    /// be handed here.
    fn build(rows: Vec<SpectrumRow<'_>>, decoys: DecoyPolicy) -> PredictedLibrary {
        let (handle, mut collector) = sink();
        for row in &rows {
            collector.spectrum(row).expect("row converts");
        }
        collector.finish().expect("stream finishes");
        handle.into_library(decoys).expect("library seals")
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
            vec![
                first.row(2, false, None, peaks(3)),
                first.row(3, false, None, peaks(2)),
                second.row(2, false, None, peaks(4)),
            ],
            DecoyPolicy::Never,
        );
        let shuffled = build(
            vec![
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
            vec![
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

    #[test]
    fn a_library_whose_rows_declare_no_pair_id_keeps_no_group_column() {
        let first = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let second = Fixture::new("PEPTIDER", "PEPTIDER");
        let lib = build(
            vec![
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
            vec![
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
            vec![
                fixture.row(2, false, None, peaks(6)),
                fixture.row(3, false, None, peaks(4)),
            ],
            DecoyPolicy::Never,
        );

        assert_eq!(
            lib.library.frag_intens.len(),
            lib.library.geom.n_fragments()
        );
        assert_eq!(lib.library.frag_intens.len(), 10);
    }

    #[test]
    fn a_predicted_library_keeps_sequence_features_for_sequences_that_parse() {
        let fixture = Fixture::new("PEPC[UNIMOD:4]IDEK", "PEPCIDEK");
        let lib = build(
            vec![fixture.row(2, false, None, peaks(3))],
            DecoyPolicy::Never,
        );

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
            vec![fixture.row(
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
        let lib = build(
            vec![fixture.row(2, false, None, peaks(2))],
            DecoyPolicy::Never,
        );

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
        let lib = build(vec![contextual], DecoyPolicy::Never);

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

        let Err(error) = handle.into_library(DecoyPolicy::Never) else {
            panic!("a library that was never finished is not a library");
        };
        assert!(
            format!("{error}").contains("did not finish"),
            "unexpected error: {error}",
        );
    }

    #[test]
    fn a_peak_whose_series_is_not_an_mzpaf_letter_fails_the_row_rather_than_losing_the_peak() {
        let fixture = Fixture::new("PEPTIDEK", "PEPTIDEK");
        let (_handle, mut collector) = sink();
        let error = collector
            .spectrum(&fixture.row(2, false, None, vec![peak("Q", 3, 1, 300.0, 1.0)]))
            .expect_err("Q is not a series letter");
        assert!(
            format!("{error:#}").contains("PEPTIDEK/2"),
            "unexpected error: {error:#}",
        );
    }
}
