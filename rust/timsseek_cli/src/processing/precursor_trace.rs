use std::collections::HashMap;
use std::path::Path;
use std::sync::Mutex;

use serde::Serialize;
use timsquery::Target;
use timsquery::models::{
    FlatIdx,
    GroupCode,
    RowIdx,
};
use timsseek::ReferenceLibrary;
use timsseek::scoring::pipeline::{
    Phase3Stage,
    Phase3Trace,
};
use timsseek::scoring::results::ScoringFields;

use crate::artifacts::PRECURSOR_TRACE_TSV;

type CompetitionKey = (GroupCode, u8);
type ScoredKey = (RowIdx, bool, u64);

struct Request {
    modified: String,
    charge: u8,
    matches: Vec<FlatIdx>,
}

/// Resolve the requested keys to library slots once, then reuse the selection
/// for each sample. Only matching rows and their decoys are kept.
pub(crate) struct TraceSelection {
    requests: Vec<Request>,
    keys: HashMap<FlatIdx, (ScoredKey, CompetitionKey)>,
    decoys: HashMap<CompetitionKey, Vec<FlatIdx>>,
}

impl TraceSelection {
    pub(crate) fn read(path: &Path, library: &ReferenceLibrary) -> std::io::Result<Self> {
        Self::from_tsv(std::fs::File::open(path)?, library)
    }

    fn from_tsv(input: impl std::io::Read, library: &ReferenceLibrary) -> std::io::Result<Self> {
        let mut input = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .from_reader(input);
        let headers = input.headers()?;
        let column = |name| {
            headers
                .iter()
                .position(|header| header == name)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("trace TSV is missing the {name} column"),
                    )
                })
        };
        let modified_col = column("modified")?;
        let charge_col = column("charge")?;
        let mut requests = Vec::new();
        let mut lookup = HashMap::<u8, HashMap<String, usize>>::new();
        let mut max_len = 0;
        for row in input.records() {
            let row = row?;
            let modified = row.get(modified_col).unwrap_or_default();
            let charge = row
                .get(charge_col)
                .and_then(|charge| charge.parse::<u8>().ok())
                .filter(|&charge| charge > 0);
            let Some(charge) = charge.filter(|_| !modified.is_empty()) else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "expected a modified sequence and positive charge on trace TSV line {}",
                        row.position().expect("CSV record position").line()
                    ),
                ));
            };
            let by_sequence = lookup.entry(charge).or_default();
            if by_sequence.contains_key(modified) {
                continue;
            }
            by_sequence.insert(modified.to_owned(), requests.len());
            max_len = max_len.max(modified.len());
            requests.push(Request {
                modified: modified.to_owned(),
                charge,
                matches: Vec::new(),
            });
        }

        let geom = library.geometry();
        let mut selection = Self {
            requests,
            keys: HashMap::new(),
            decoys: HashMap::new(),
        };
        let mut sequence = SequenceBuffer {
            text: String::with_capacity(max_len),
            max_len,
        };
        for row in geom.rows() {
            let Some(by_sequence) = lookup.get(&geom.charge(row)) else {
                continue;
            };
            let Some(peptide) = geom.analyte(row).peptide.known() else {
                continue;
            };
            let Some(mods) = peptide.modifications.known() else {
                continue;
            };
            let modified = if mods.is_empty() {
                peptide.residues
            } else {
                sequence.text.clear();
                if peptide.write_sequence(&mut sequence).is_err() {
                    continue;
                }
                &sequence.text
            };
            if let Some(&request) = by_sequence.get(modified) {
                let flat = geom.flat_for(row, 0);
                let (key, group) = selection.insert_slot(library, flat);
                selection.requests[request].matches.push(flat);
                if key.1 {
                    selection.decoys.entry(group).or_default();
                    for variant in 1..geom.variants_per_row() {
                        let variant = u8::try_from(variant).expect("library variant fits u8");
                        let decoy = geom.flat_for(row, variant);
                        selection.insert_slot(library, decoy);
                        selection.decoys.get_mut(&group).unwrap().push(decoy);
                    }
                }
            }
        }
        if geom.n_stored_decoys() > 0 && !selection.decoys.is_empty() {
            for row in geom.rows().filter(|&row| geom.is_decoy(row)) {
                let group = (geom.decoy_group_code(row), geom.charge(row));
                if let Some(decoys) = selection.decoys.get_mut(&group) {
                    let decoy = geom.flat_for(row, 0);
                    decoys.push(decoy);
                    selection.insert_slot(library, decoy);
                }
            }
        }
        Ok(selection)
    }

    fn insert_slot(
        &mut self,
        library: &ReferenceLibrary,
        flat: FlatIdx,
    ) -> (ScoredKey, CompetitionKey) {
        let geom = library.geometry();
        let (row, _) = geom.split_flat(flat);
        let key = (
            row,
            geom.is_target(flat),
            library.item_at(flat).mono_precursor_mz().to_bits(),
        );
        let group = (geom.decoy_group_code(row), geom.charge(row));
        self.keys.insert(flat, (key, group));
        (key, group)
    }
}

/// Use the same sequence formatter as Parquet. A string longer than every
/// requested key cannot match.
struct SequenceBuffer {
    text: String,
    max_len: usize,
}

impl std::fmt::Write for SequenceBuffer {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        if text.len() > self.max_len - self.text.len() {
            return Err(std::fmt::Error);
        }
        self.text.push_str(text);
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Winner {
    row: RowIdx,
    is_target: bool,
    mz: f64,
    score: f32,
}

#[derive(Default)]
enum Phase4 {
    #[default]
    NotReached,
    Deduplicated(Winner),
    Lost(Winner),
    Won,
}

#[derive(Default)]
struct Downstream {
    phase4: Phase4,
    qvalue: Option<f32>,
    written: bool,
}

impl Downstream {
    fn progress(&self) -> u8 {
        if self.written {
            5
        } else if self.qvalue.is_some() {
            4
        } else {
            match self.phase4 {
                Phase4::Won => 3,
                Phase4::Lost(_) => 2,
                Phase4::Deduplicated(_) => 1,
                Phase4::NotReached => 0,
            }
        }
    }
}

pub(super) struct PrecursorTrace<'a> {
    selection: &'a TraceSelection,
    pub(super) phase3: HashMap<FlatIdx, Mutex<Phase3Trace>>,
    downstream: HashMap<ScoredKey, Downstream>,
}

impl<'a> PrecursorTrace<'a> {
    pub(super) fn new(selection: &'a TraceSelection) -> Self {
        let phase3 = selection
            .keys
            .keys()
            .map(|&flat| (flat, Mutex::default()))
            .collect();
        let downstream: HashMap<_, _> = selection
            .keys
            .values()
            .map(|&(key, _)| (key, Downstream::default()))
            .collect();
        assert_eq!(
            selection.keys.len(),
            downstream.len(),
            "ambiguous traced scoring slot"
        );
        Self {
            selection,
            phase3,
            downstream,
        }
    }

    fn downstream(&mut self, scoring: &ScoringFields) -> Option<&mut Downstream> {
        let id = &scoring.identity;
        self.downstream
            .get_mut(&(id.row, id.is_target, id.precursor_mz.to_bits()))
    }

    pub(super) fn dropped(
        &mut self,
        loser: &ScoringFields,
        winner: &ScoringFields,
        deduplicated: bool,
    ) {
        if let Some(row) = self.downstream(loser) {
            let winner = Winner {
                row: winner.identity.row,
                is_target: winner.identity.is_target,
                mz: winner.identity.precursor_mz,
                score: winner.primary.main_score,
            };
            row.phase4 = if deduplicated {
                Phase4::Deduplicated(winner)
            } else {
                Phase4::Lost(winner)
            };
        }
    }

    pub(super) fn won(&mut self, winner: &ScoringFields) {
        if let Some(row) = self.downstream(winner) {
            row.phase4 = Phase4::Won;
        }
    }

    pub(super) fn final_q(&mut self, scoring: &ScoringFields, qvalue: f32, written: bool) {
        if let Some(row) = self.downstream(scoring) {
            row.qvalue = Some(qvalue);
            row.written = written;
        }
    }

    pub(super) fn write(
        &mut self,
        library: &ReferenceLibrary,
        parquet_path: &Path,
        raw: bool,
    ) -> Result<(), timsseek::errors::TimsSeekError> {
        let path = parquet_path.with_file_name(PRECURSOR_TRACE_TSV);
        let result =
            std::fs::File::create_new(&path).and_then(|file| self.write_tsv(library, file, raw));
        result.map_err(|source| timsseek::errors::TimsSeekError::Io {
            source,
            path: Some(path),
        })
    }

    fn write_tsv(
        &mut self,
        library: &ReferenceLibrary,
        output: impl std::io::Write,
        raw: bool,
    ) -> std::io::Result<()> {
        let geom = library.geometry();
        let mut output = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(output);
        if self.selection.requests.is_empty() {
            output.write_record(TRACE_COLUMNS)?;
        }
        for request in &self.selection.requests {
            // A peptide/charge may match multiple rows. Prefer targets, then
            // the furthest progress, then the highest main score.
            let rank = |flat| {
                let trace = self.phase3[&flat].lock().expect("trace slot poisoned");
                let (key, _) = self.selection.keys[&flat];
                (
                    (key.1, self.downstream[&key].progress(), trace.stage),
                    trace.main_score.unwrap_or(f32::NEG_INFINITY),
                )
            };
            let best = request.matches.iter().copied().max_by(|&a, &b| {
                let (a_stage, a_score) = rank(a);
                let (b_stage, b_score) = rank(b);
                a_stage
                    .cmp(&b_stage)
                    .then_with(|| a_score.total_cmp(&b_score))
            });

            let base = TraceRow {
                modified: &request.modified,
                charge: request.charge,
                ..Default::default()
            };
            let Some(best) = best else {
                output.serialize(TraceRow {
                    in_library: "n",
                    reached_stages: "library".into(),
                    dropped_stage: Some("library"),
                    reason: Some("not in library"),
                    ..base
                })?;
                continue;
            };
            let (key, group) = self.selection.keys[&best];
            let paired = self.selection.decoys.get(&group).filter(|_| key.1);
            for flat in std::iter::once(&best).chain(paired.into_iter().flatten()) {
                let (key, _) = self.selection.keys[flat];
                let (row, is_target, mz_bits) = key;
                let (_, variant) = geom.split_flat(*flat);
                let phase3 = self
                    .phase3
                    .get_mut(flat)
                    .unwrap()
                    .get_mut()
                    .expect("trace slot poisoned");
                let end = &self.downstream[&key];
                let (last_stage, dropped_stage, reason, lost_competition, winner) =
                    fate(phase3, end, raw);
                let rt = phase3.rt_window;
                output.serialize(TraceRow {
                    in_library: "y",
                    sequence: geom.analyte(row).peptide.known().and_then(|p| p.sequence()),
                    library_id: Some(geom.output_id(row).to_string()),
                    is_target: Some(if is_target { "y" } else { "n" }),
                    variant: Some(variant),
                    precursor_mz: Some(f64::from_bits(mz_bits)),
                    reached_stages: reached_stages(phase3.stage, last_stage, raw),
                    dropped_stage,
                    reason,
                    main_score: phase3.main_score,
                    apex_rt_seconds: phase3.apex_rt_seconds,
                    library_rt: geom.library_rt(row).map(|rt| rt.0),
                    calibrated_rt_seconds: rt.and_then(|rt| rt.center()).map(|center| center.0),
                    rt_window_start_seconds: rt
                        .map(|rt| f64::from(rt.range_millis().start()) / 1000.0),
                    rt_window_end_seconds: rt.map(|rt| f64::from(rt.range_millis().end()) / 1000.0),
                    matched_fragments: phase3.matched_fragments,
                    lost_competition,
                    winner_library_id: winner.map(|w| geom.output_id(w.row).to_string()),
                    winner_sequence: winner
                        .and_then(|w| geom.analyte(w.row).peptide.known())
                        .and_then(|p| p.sequence()),
                    winner_is_target: winner.map(|w| if w.is_target { "y" } else { "n" }),
                    winner_precursor_mz: winner.map(|w| w.mz),
                    winner_main_score: winner.map(|w| w.score),
                    qvalue: end.qvalue,
                    ..base
                })?;
            }
        }
        output.flush()
    }
}

type Fate = (
    &'static str,
    Option<&'static str>,
    Option<&'static str>,
    Option<&'static str>,
    Option<Winner>,
);

const PHASE3_STAGES: [(Phase3Stage, &str); 6] = [
    (Phase3Stage::PrecursorRange, "phase3_precursor_range"),
    (Phase3Stage::ExpectedFragments, "phase3_expected_fragments"),
    (Phase3Stage::Extraction, "phase3_extraction"),
    (Phase3Stage::Apex, "phase3_apex"),
    (Phase3Stage::Finalize, "phase3_finalize"),
    (Phase3Stage::Scored, "phase3_scored"),
];
const SUPERVISED_STAGES: [&str; 4] = [
    "phase4_deduplication",
    "phase4_competition",
    "phase5_rescore",
    "phase6_qvalue_filter",
];

fn reached_stages(phase3: Phase3Stage, last_stage: &str, raw: bool) -> String {
    let mut stages = String::from("library");
    let later = (!raw).then_some(SUPERVISED_STAGES).into_iter().flatten();
    for stage in PHASE3_STAGES
        .iter()
        .take_while(|&&(stage, _)| stage <= phase3)
        .map(|&(_, name)| name)
        .chain(later)
        .chain(["phase6_written"])
    {
        stages.push('|');
        stages.push_str(stage);
        if stage == last_stage {
            break;
        }
    }
    stages
}

fn fate(phase3: &Phase3Trace, end: &Downstream, raw: bool) -> Fate {
    let stage = PHASE3_STAGES
        .iter()
        .find(|&&(stage, _)| stage == phase3.stage)
        .expect("Phase 3 stage has a trace label")
        .1;
    if let Some(skip) = phase3.skip {
        return (stage, Some(stage), Some(skip.as_str()), None, None);
    }
    match end.phase4 {
        Phase4::Deduplicated(winner) => (
            "phase4_deduplication",
            Some("phase4_deduplication"),
            Some("duplicate sequence, charge and precursor m/z"),
            None,
            Some(winner),
        ),
        Phase4::Lost(winner) => (
            "phase4_competition",
            Some("phase4_competition"),
            Some("lost target-decoy competition"),
            Some("y"),
            Some(winner),
        ),
        Phase4::Won if end.written => ("phase6_written", None, None, Some("n"), None),
        Phase4::Won if end.qvalue.is_some() => (
            "phase6_qvalue_filter",
            Some("phase6_qvalue_filter"),
            Some("did not pass output q-value threshold"),
            Some("n"),
            None,
        ),
        Phase4::Won => (
            "phase5_rescore",
            Some("phase5_rescore"),
            Some("no final q-value"),
            Some("n"),
            None,
        ),
        Phase4::NotReached if raw => ("phase6_written", None, None, None, None),
        Phase4::NotReached => (stage, None, None, None, None),
    }
}

const TRACE_COLUMNS: &[&str] = &[
    "modified",
    "charge",
    "in_library",
    "sequence",
    "library_id",
    "is_target",
    "variant",
    "precursor_mz",
    "reached_stages",
    "dropped_stage",
    "reason",
    "main_score",
    "apex_rt_seconds",
    "library_rt",
    "calibrated_rt_seconds",
    "rt_window_start_seconds",
    "rt_window_end_seconds",
    "matched_fragments",
    "lost_competition",
    "winner_library_id",
    "winner_sequence",
    "winner_is_target",
    "winner_precursor_mz",
    "winner_main_score",
    "qvalue",
];

#[derive(Default, Serialize)]
struct TraceRow<'a> {
    modified: &'a str,
    charge: u8,
    in_library: &'static str,
    sequence: Option<String>,
    library_id: Option<String>,
    is_target: Option<&'static str>,
    variant: Option<u8>,
    precursor_mz: Option<f64>,
    reached_stages: String,
    dropped_stage: Option<&'static str>,
    reason: Option<&'static str>,
    main_score: Option<f32>,
    apex_rt_seconds: Option<f32>,
    library_rt: Option<f32>,
    calibrated_rt_seconds: Option<f32>,
    rt_window_start_seconds: Option<f64>,
    rt_window_end_seconds: Option<f64>,
    matched_fragments: Option<u8>,
    lost_competition: Option<&'static str>,
    winner_library_id: Option<String>,
    winner_sequence: Option<String>,
    winner_is_target: Option<&'static str>,
    winner_precursor_mz: Option<f64>,
    winner_main_score: Option<f32>,
    qvalue: Option<f32>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use timscentroid::rt_mapping::{
        CycleToRTMapping,
        MS1CycleIndex,
    };
    use timsquery::chemistry::analyte::Analyte;
    use timsquery::models::{
        Row,
        TargetColumnsBuilder,
    };
    use timsquery::{
        IonAnnot,
        PeakTolerance,
        QueriableData,
    };
    use timsseek::Scorer;
    use timsseek::traits::MappableRTCycles;

    struct EmptyIndex(CycleToRTMapping<MS1CycleIndex>);

    impl MappableRTCycles for EmptyIndex {
        fn ms1_cycle_mapping(&self) -> &CycleToRTMapping<MS1CycleIndex> {
            &self.0
        }

        fn mobility_kind(&self) -> &timscentroid::MobilityKind {
            &timscentroid::MobilityKind::Ook0
        }
    }

    impl<A: Send + Sync> QueriableData<A> for EmptyIndex {
        fn add_query(&self, _: &mut A, _: &PeakTolerance) {}
    }

    fn library_with_mass_shift_decoys() -> ReferenceLibrary {
        let mut builder = TargetColumnsBuilder::with_capabilities(
            timsquery::models::TargetCapabilities::default_diann(),
        );
        let analyte = Analyte::from_sequence("SC[UNIMOD:4]LYVK");
        builder.push_row(Row {
            analyte: analyte.as_input(),
            charge: 2,
            precursor_mz: 500.0,
            frags: &[(IonAnnot::try_from("y1").unwrap(), 300.0)],
            ..Default::default()
        });
        ReferenceLibrary::try_from(timsquery::serde::TargetTable::Mzpaf {
            geom: builder.seal(timsseek::DecoyPolicy::Force).unwrap(),
            frag_intens: Some(vec![1.0]),
        })
        .unwrap()
    }

    #[test]
    fn traces_phase3_skip_and_missing_library_key() {
        let library = library_with_mass_shift_decoys();
        let tsv =
            b"modified\tcharge\tbucket\nSC[UNIMOD:4]LYVK\t2\tabsent\nSC[UNIMOD:4]LYVK\t3\tnear\n";
        let selection = TraceSelection::from_tsv(tsv.as_slice(), &library).unwrap();
        let mut trace = PrecursorTrace::new(&selection);
        let scorer = Scorer {
            index: EmptyIndex(CycleToRTMapping::new(vec![0, 1_000])),
            broad_tolerance: Default::default(),
            fragmented_range: (400.0, 600.0).try_into().unwrap(),
        };
        let batch: Vec<_> = library.geometry().flats().collect();
        let (scored, _, skips) = scorer.score_traced_batch(&library, &batch, None, &trace.phase3);
        assert!(scored.is_empty());
        assert_eq!(skips.fragments_outside_scan_range, 3);

        let mut tsv = Vec::new();
        trace.write_tsv(&library, &mut tsv, false).unwrap();
        let mut output = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(tsv.as_slice());
        let columns = output.headers().unwrap().clone();
        assert!(columns.iter().eq(TRACE_COLUMNS.iter().copied()));
        let column = |name| columns.iter().position(|field| field == name).unwrap();
        let rows: Vec<_> = output.records().map(Result::unwrap).collect();
        assert_eq!(rows.len(), 4);
        for (i, row) in rows[..3].iter().enumerate() {
            assert_eq!(&row[column("charge")], "2");
            assert_eq!(&row[column("in_library")], "y");
            assert_eq!(&row[column("is_target")], if i == 0 { "y" } else { "n" });
            assert_eq!(&row[column("dropped_stage")], "phase3_extraction");
            assert_eq!(
                &row[column("reached_stages")],
                "library|phase3_precursor_range|phase3_expected_fragments|phase3_extraction"
            );
            assert_eq!(&row[column("reason")], "fragments outside scan range");
            assert_eq!(&row[column("rt_window_start_seconds")], "0.0");
            assert_eq!(&row[column("rt_window_end_seconds")], "1.0");
        }
        assert_eq!(&rows[3][column("charge")], "3");
        assert_eq!(&rows[3][column("in_library")], "n");
        assert_eq!(&rows[3][column("reason")], "not in library");
    }

    #[test]
    fn traces_competition_winner_and_filtered_qvalue() {
        use timsseek::ScoredCandidate;
        use timsseek::data_sources::reference_library::ScoredIdentity;

        let library = library_with_mass_shift_decoys();
        let selection = TraceSelection::from_tsv(
            b"modified\tcharge\nSC[UNIMOD:4]LYVK\t2\n".as_slice(),
            &library,
        )
        .unwrap();
        let mut trace = PrecursorTrace::new(&selection);
        let candidates = library
            .geometry()
            .flats()
            .zip([9.0, 7.0, 10.0])
            .map(|(flat, score)| {
                let query = library.item_at(flat);
                let mut scoring = ScoringFields::sample_default();
                scoring.identity.row = query.handles().row;
                scoring.identity.group = query.handles().group;
                scoring.identity.is_target = query.is_target();
                scoring.identity.precursor_mz = query.mono_precursor_mz();
                scoring.primary.main_score = score;
                *trace.phase3.get_mut(&flat).unwrap().get_mut().unwrap() = Phase3Trace {
                    stage: Phase3Stage::Scored,
                    main_score: Some(score),
                    ..Default::default()
                };
                ScoredCandidate { scoring }
            })
            .collect();
        let competed = super::super::target_decoy_compete(candidates, &library, Some(&mut trace));
        let [winner] = competed.as_slice() else {
            panic!("expected one competition winner");
        };
        assert_eq!(winner.scoring.primary.main_score, 10.0);
        trace.final_q(&winner.scoring, 0.05, false);

        let mut tsv = Vec::new();
        trace.write_tsv(&library, &mut tsv, false).unwrap();
        let mut output = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(tsv.as_slice());
        let columns = output.headers().unwrap().clone();
        let column = |name| columns.iter().position(|field| field == name).unwrap();
        let rows: Vec<_> = output.records().map(Result::unwrap).collect();
        assert_eq!(&rows[0][column("lost_competition")], "y");
        assert_eq!(&rows[0][column("winner_main_score")], "10.0");
        assert_eq!(
            &rows[0][column("winner_precursor_mz")],
            &rows[2][column("precursor_mz")]
        );
        assert_eq!(&rows[0][column("qvalue")], "");
        assert_eq!(&rows[2][column("lost_competition")], "n");
        assert_eq!(&rows[2][column("qvalue")], "0.05");
        assert_eq!(&rows[2][column("dropped_stage")], "phase6_qvalue_filter");
        assert!(rows[2][column("reached_stages")].ends_with(
            "|phase4_deduplication|phase4_competition|phase5_rescore|phase6_qvalue_filter"
        ));
    }
}
