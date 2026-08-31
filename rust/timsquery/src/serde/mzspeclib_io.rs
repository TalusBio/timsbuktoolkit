//! mzSpecLib, read through [`mzannotate`].
//!
//! The parse is not ours. What is ours is the mapping onto the arena, and that
//! mapping is where the formats differ: the exports we have seen disagree on
//! which term carries the precursor m/z, whether a retention time is declared
//! at all, and which mobility term they use. Every field below is a ladder over
//! those spellings rather than a single lookup.
//!
//! Peaks split three ways. A peak whose annotation resolves to one representable
//! ion is stored at its theoretical m/z; a peak that resolves but cannot be
//! spelled as an [`IonAnnot`] gets an unknown label at the same theoretical m/z;
//! a peak that is unannotated, or annotated several ways at once, is skipped.
//! The third case is skipped rather than stored at its observed m/z, because an
//! arena mixing observed and theoretical masses is invisible downstream.
//!
//! Finding that theoretical m/z is itself a ladder, because the exports differ
//! again: the composition when the annotation carries one, the mzPAF mass-error
//! suffix when it does not, and the peak's own m/z when the spectrum declares
//! its masses calculated rather than measured. See [`theoretical_mz`].
//!
//! What a library loses on the way in is counted rather than logged per row,
//! and reported once. A field being absent and a field being unreadable are
//! different, and the counts are what tell them apart; see [`Degradation`].

use std::io::BufRead;
use std::path::Path;

use micromzpaf::{
    IonAnnot,
    IonSeriesOrdinal,
    NeutralLoss,
    UnknownIonCounter,
};
use mzannotate::fragment::{
    Fragment,
    FragmentType,
};
use mzannotate::mzdata::params::{
    ControlledVocabulary,
    ParamValue,
};
use mzannotate::mzdata::prelude::PeakCollection;
use mzannotate::mzspeclib::{
    AnalyteTarget,
    AttributeValue,
    MzSpecLibTextParser,
};
use mzannotate::spectrum::AnnotatedSpectrum;
use mzcore::chemistry::{
    AmbiguousMolecule,
    MassMode,
    MolecularCharge,
    Molecule,
};
use mzcore::prelude::IsAminoAcid;
use tracing::warn;

use super::library_file::{
    TargetReadingError,
    TargetTable,
};
use super::psims_origin_type;
use crate::chemistry::ontologies;
use crate::models::{
    Row,
    TargetCapabilities,
    TargetColumns,
};

// ── The controlled vocabulary ────────────────────────────────────────────────
//
// Every bare `u32` accession in this file is a PSI-MS accession, and comparing
// them as integers is only correct because the namespace is established first.
// There are exactly two places that establish it, and every integer downstream
// of them has already passed through one:
//
//   - `is_ms_term`, for values mzannotate flattened into an mzdata `Param`,
//     which splits the CURIE into a number and a namespace.
//   - the `MS:` prefix check in `origin_types`, for a term arriving as text.
//
// It matters because PSI-MS and the Unit Ontology both allocate 7-digit
// accessions in overlapping ranges: UO already reaches 1010060, above every
// accession named here. No colliding term is allocated today, so this guards a
// range overlap rather than a live bug.
//
// Where the whole `Curie` survives, `mzcv::curie!` compares both halves at once
// and no constant is needed.

/// `MS:1003072|spectrum origin type`.
const SPECTRUM_ORIGIN_TYPE: u32 = 1_003_072;
/// `MS:1002815|inverse reduced ion mobility`, the timsTOF mobility axis.
const INVERSE_REDUCED_MOBILITY: u32 = 1_002_815;
/// `MS:1002476|ion mobility drift time`, the drift-tube spelling.
const ION_MOBILITY_DRIFT_TIME: u32 = 1_002_476;
/// msspeculator's decoy-to-target link, carried as a project-defined name/value
/// pair because no PSI-MS term identifies the target a decoy was derived from.
const MSSPECULATOR_PAIR_ID: &str = "msspeculator:decoy_pair_id";

/// Whether `accession` names a decoy spectrum.
///
/// `None` for a term outside [`psims_origin_type::ALL`], which is a term this
/// build has never heard of. The caller reports those rather than reading them
/// as targets: guessing wrong on an unknown decoy subtype yields an FDR that is
/// wrong instead of absent.
///
/// The `is_a` closure was resolved against the published vocabulary by
/// `scripts/gen_psims_origin_type.py`, so this is a lookup rather than a walk.
fn subsumes_decoy(accession: u32) -> Option<bool> {
    psims_origin_type::ALL
        .binary_search(&accession)
        .ok()
        .map(|_| psims_origin_type::DECOY.binary_search(&accession).is_ok())
}

/// Whether this spectrum's peak m/z values are calculated rather than measured.
///
/// It decides what an annotation with no composition and no mass-error suffix
/// means. On a predicted spectrum the peak m/z already *is* the theoretical
/// mass, so the peak is usable; on an observed one there is no theoretical mass
/// to be had and the peak is skipped rather than mixed in at its observed m/z.
fn declares_theoretical_mz(accession: u32) -> bool {
    psims_origin_type::THEORETICAL_MZ
        .binary_search(&accession)
        .is_ok()
}

// ── Reading a file ───────────────────────────────────────────────────────────

/// True for a file this reader should try, by extension or by first line.
///
/// The extension is spelled several ways in the wild (`.mzSpecLib.txt` in the
/// specification's own examples, `.mzspeclib.txt` elsewhere) and gzipped
/// libraries add `.gz`, so the header check is what actually decides for
/// anything unusual.
pub fn sniff_mzspeclib_library_file(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if name.contains(".mzspeclib.") {
        return true;
    }
    // gzip is opaque to a text probe, so an unrecognised name plus a `.gz`
    // suffix is not this format.
    if name.ends_with(".gz") {
        return false;
    }
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let mut first = String::new();
    std::io::BufReader::new(file).read_line(&mut first).is_ok() && first.trim_end() == "<mzSpecLib>"
}

/// Open `path`, transparently decompressing a gzipped library.
///
/// `.gz` is the recommended spelling for a stored library, so reading it is not
/// optional. A proteome-scale library is about 1.2 GB gzipped against 6.1 GB
/// plain.
fn open_reader(path: &Path) -> Result<Box<dyn BufRead>, TargetReadingError> {
    let file = std::fs::File::open(path).map_err(TargetReadingError::IoError)?;
    let gzipped = path
        .file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.to_ascii_lowercase().ends_with(".gz"));
    Ok(if gzipped {
        Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    })
}

/// Read a library into a sealed arena.
///
/// Builds the arena directly rather than going through
/// [`ElutionGroupCollection`](super::library_file::ElutionGroupCollection),
/// because the decoy groups a library declares have no representation there.
pub fn read_mzspeclib_library_file(path: &Path) -> Result<TargetTable, TargetReadingError> {
    let reader = open_reader(path)?;
    let parser = MzSpecLibTextParser::open(reader, Some(path.to_path_buf()), ontologies())
        .map_err(|e| TargetReadingError::SpeclibParse(format!("mzSpecLib header: {e}")))?;

    let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
    let mut frag_intens: Vec<f32> = Vec::new();
    let mut groups: Vec<String> = Vec::new();
    let mut any_group_declared = false;
    let mut degradation = Degradation::default();

    for (index, spectrum) in parser.enumerate() {
        let spectrum = spectrum.map_err(|e| {
            TargetReadingError::SpeclibParse(format!("mzSpecLib spectrum {}: {e}", index + 1))
        })?;
        let row = SpectrumRow::extract(&spectrum, &mut degradation)?;
        // Both halves of a pair have to land in one namespace. A decoy names
        // its target by `<Spectrum=N>` key, so an unpaired row is its own key
        // rather than its source id, which is a different string entirely.
        any_group_declared |= row.declared_group.is_some();
        groups.push(row.declared_group.unwrap_or_else(|| row.key.to_string()));

        frag_intens.extend(row.frags.iter().map(|(_, _, intensity)| *intensity));
        let frags: Vec<(IonAnnot, f64)> = row
            .frags
            .iter()
            .map(|(label, mz, _)| (*label, *mz))
            .collect();
        geom.push_row(Row {
            precursor_mz: row.precursor_mz,
            charge: row.charge,
            rt_seconds: row.rt_seconds,
            mobility: row.mobility,
            frags: &frags,
            seq_strip: &row.seq_strip,
            seq_mod: &row.seq_mod,
            is_decoy: row.is_decoy,
            id: Some(row.id.into()),
            ..Default::default()
        });
    }

    // A file that declares no pairing gets no group column: a row is then its
    // own group and nothing is allocated.
    if any_group_declared {
        geom.set_decoy_groups(groups)?;
    }
    geom.seal()?;
    degradation.report(geom.n_rows());

    Ok(TargetTable::Mzpaf {
        geom,
        frag_intens: Some(frag_intens),
    })
}

/// What a library lost on the way into the arena, counted rather than logged
/// per row: a library with a million unannotated peaks should produce one line,
/// not a million.
#[derive(Default)]
struct Degradation {
    unannotated_peaks: usize,
    ambiguous_peaks: usize,
    unrepresentable_labels: usize,
    peaks_without_theoretical_mz: usize,
    rows_without_sequence: usize,
    // No count of entries missing a retention time. Zero is an ordinary value
    // on a normalized scale rather than an absence -- it is an anchor point, and
    // it is what msspeculator's own libraries write -- so there is nothing to
    // test for. See the CONTEXT.md entry on retention time.
    /// A set, not a list: this one is pushed per spectrum rather than per file,
    /// so a five-million-row library with one unrecognised subtype would
    /// otherwise hold a 20 MB `Vec` and print five million identical lines.
    unrecognised_origin_types: std::collections::BTreeSet<u32>,
}

impl Degradation {
    fn report(&self, rows: usize) {
        if self.unannotated_peaks > 0 || self.ambiguous_peaks > 0 {
            warn!(
                "mzSpecLib: skipped {} unannotated and {} ambiguously annotated peak(s) across {rows} entries; \
                 an arena mixing observed and theoretical m/z would be invisible downstream",
                self.unannotated_peaks, self.ambiguous_peaks
            );
        }
        if self.unrepresentable_labels > 0 {
            warn!(
                "mzSpecLib: {} peak(s) carry an ion this build cannot spell and were labelled unknown",
                self.unrepresentable_labels
            );
        }
        if self.peaks_without_theoretical_mz > 0 {
            warn!(
                "mzSpecLib: {} peak(s) have neither a composition nor a mass-error suffix and were skipped",
                self.peaks_without_theoretical_mz
            );
        }
        if self.rows_without_sequence > 0 {
            warn!(
                "mzSpecLib: {} of {rows} entries carry no peptidoform; sequence features are unavailable for them",
                self.rows_without_sequence
            );
        }
        for accession in &self.unrecognised_origin_types {
            warn!(
                "mzSpecLib: spectrum origin type MS:{accession} is outside this build's copy of \
                 the vocabulary (psi-ms {}) and was read as a target. Re-run \
                 scripts/gen_psims_origin_type.py if the term is newer than that.",
                psims_origin_type::DATA_VERSION,
            );
        }
    }
}

/// One spectrum, reduced to what the arena stores.
struct SpectrumRow {
    /// The `<Spectrum=N>` key, which is what a declared pair refers to.
    key: u32,
    id: String,
    precursor_mz: f64,
    charge: u8,
    rt_seconds: f32,
    mobility: f32,
    is_decoy: bool,
    declared_group: Option<String>,
    seq_strip: String,
    seq_mod: String,
    /// `(label, theoretical m/z, intensity)`.
    frags: Vec<(IonAnnot, f64, f32)>,
}

impl SpectrumRow {
    fn extract<M: mzcore::chemistry::MassOutputMode>(
        spectrum: &AnnotatedSpectrum<M>,
        degradation: &mut Degradation,
    ) -> Result<Self, TargetReadingError> {
        reject_unresolved_attribute_sets(spectrum)?;

        let analyte = spectrum.analytes.first();
        let peptidoform = analyte.and_then(|a| match &a.target {
            AnalyteTarget::PeptidoformIon(ion) => Some(ion),
            _ => None,
        });

        let (seq_strip, seq_mod) = match peptidoform {
            Some(ion) => (stripped_sequence(ion), ion.to_string()),
            None => {
                degradation.rows_without_sequence += 1;
                (String::new(), String::new())
            }
        };

        let peak_mz_is_theoretical = origin_types(spectrum).any(declares_theoretical_mz);
        let mut counter = UnknownIonCounter::default();
        let mut frags = Vec::with_capacity(spectrum.peaks.len());
        for peak in spectrum.peaks.iter() {
            let observed = peak.mz.value;
            match peak.annotations.as_slice() {
                [] => {
                    degradation.unannotated_peaks += 1;
                    continue;
                }
                [annotation] => {
                    let Some(mz) = theoretical_mz(annotation, observed, peak_mz_is_theoretical)
                    else {
                        degradation.peaks_without_theoretical_mz += 1;
                        continue;
                    };
                    let label = match to_ion_annot(annotation) {
                        Some(label) => label,
                        None => {
                            degradation.unrepresentable_labels += 1;
                            counter
                                .next_unknown(annotation.charge.value.max(1) as i8)
                                .map_err(|e| {
                                    TargetReadingError::SpeclibParse(format!(
                                        "ran out of unknown-ion labels: {e}"
                                    ))
                                })?
                        }
                    };
                    frags.push((label, mz, peak.intensity));
                }
                // Several identities for one peak: no single theoretical mass
                // exists, so there is nothing honest to store.
                _ => degradation.ambiguous_peaks += 1,
            }
        }

        Ok(Self {
            key: spectrum.key,
            id: source_id(spectrum),
            precursor_mz: precursor_mz(spectrum, peptidoform),
            charge: charge(spectrum, peptidoform),
            rt_seconds: rt_seconds(spectrum),
            mobility: mobility(spectrum),
            is_decoy: is_decoy(spectrum, degradation),
            declared_group: declared_group(spectrum),
            seq_strip,
            seq_mod,
            frags,
        })
    }
}

/// A claimed attribute set the header never defined is an error.
///
/// Reading it as a target produces an FDR number that is wrong rather than
/// absent, and the decoy declaration is exactly what such a set usually carries.
/// mzannotate drops an unresolvable claim silently but leaves the claim itself
/// in the spectrum's leftover attributes, which is what this sees.
fn reject_unresolved_attribute_sets<M: mzcore::chemistry::MassOutputMode>(
    spectrum: &AnnotatedSpectrum<M>,
) -> Result<(), TargetReadingError> {
    // A resolved claim contributes its own attributes, so the origin type is
    // present; an unresolved one leaves only the claim.
    let origin_type_present = spectrum
        .description
        .params
        .iter()
        .any(|p| is_ms_term(p, SPECTRUM_ORIGIN_TYPE));
    for attribute in spectrum.attributes.iter().flatten() {
        if attribute.name.accession != mzcv::curie!(MS:1003212) {
            continue;
        }
        if !origin_type_present {
            return Err(TargetReadingError::SpeclibParse(format!(
                "spectrum {} claims attribute set {} which the header does not define; \
                 refusing to read its entries as targets",
                spectrum.key, attribute.value
            )));
        }
    }
    Ok(())
}

/// The unmodified residues, one letter each. `X` for a residue with no
/// one-letter code, matching what mzannotate's own writer emits.
fn stripped_sequence(ion: &mzcore::sequence::PeptidoformIon) -> String {
    ion.peptidoforms()
        .first()
        .map(|peptide| {
            peptide
                .sequence()
                .iter()
                .map(|s| s.aminoacid.aminoacid().one_letter_code().unwrap_or('X'))
                .collect()
        })
        .unwrap_or_default()
}

/// What the file called this precursor. `MS:1003061|library spectrum name` when
/// it has one, otherwise the `<Spectrum=N>` key.
fn source_id<M: mzcore::chemistry::MassOutputMode>(spectrum: &AnnotatedSpectrum<M>) -> String {
    if spectrum.description.id.is_empty() {
        spectrum.key.to_string()
    } else {
        spectrum.description.id.clone()
    }
}

/// `MS:1000744|selected ion m/z` (DIA-NN), then
/// `MS:1003208|experimental precursor monoisotopic m/z` (Spectronaut, which
/// mzannotate routes to the isolation window), then the analyte's own mass.
fn precursor_mz<M: mzcore::chemistry::MassOutputMode>(
    spectrum: &AnnotatedSpectrum<M>,
    peptidoform: Option<&mzcore::sequence::PeptidoformIon>,
) -> f64 {
    if let Some(precursor) = spectrum.description.precursor.first() {
        if let Some(mz) = precursor.ions.first().map(|i| i.mz).filter(|mz| *mz > 0.0) {
            return mz;
        }
        let target = f64::from(precursor.isolation_window.target);
        if target > 0.0 {
            return target;
        }
    }
    // Last resort: the analyte's own composition, the way mzannotate's writer
    // derives `MS:1003053|theoretical monoisotopic m/z` from it.
    peptidoform
        .and_then(|ion| {
            let carriers = ion.get_charge_carriers()?;
            let charge = carriers.charge().value as f64;
            if charge == 0.0 {
                return None;
            }
            let formula = (ion.formulas() + carriers.formula()).to_vec().pop()?;
            Some(formula.monoisotopic_mass().value / charge)
        })
        .unwrap_or(0.0)
}

/// `MS:1000041|charge state`, then the peptidoform's own charge carriers.
fn charge<M: mzcore::chemistry::MassOutputMode>(
    spectrum: &AnnotatedSpectrum<M>,
    peptidoform: Option<&mzcore::sequence::PeptidoformIon>,
) -> u8 {
    spectrum
        .description
        .precursor
        .first()
        .and_then(|p| p.ions.first())
        .and_then(|i| i.charge)
        .or_else(|| {
            peptidoform
                .and_then(|ion| ion.get_charge_carriers().map(MolecularCharge::charge))
                .map(|c| c.value as i32)
        })
        .and_then(|c| u8::try_from(c.abs()).ok())
        .unwrap_or(0)
}

/// Whatever the file's retention time is, scaled from the minutes mzannotate
/// normalises every declared unit into.
///
/// Not necessarily a duration. `MS:1000896|normalized retention time` is an
/// index on a reference scale, which is what msspeculator writes and what the
/// Spectronaut export carries, so this can hold a dimensionless value where
/// zero and negative numbers are ordinary. Calibration fits a monotone path
/// from it to observed time, so ordering is what has to survive here, not
/// units. See the CONTEXT.md entry on retention time.
fn rt_seconds<M: mzcore::chemistry::MassOutputMode>(spectrum: &AnnotatedSpectrum<M>) -> f32 {
    const SECONDS_PER_MINUTE: f64 = 60.0;
    spectrum
        .description
        .acquisition
        .scans
        .first()
        .map_or(0.0, |scan| (scan.start_time * SECONDS_PER_MINUTE) as f32)
}

/// `MS:1002815|inverse reduced ion mobility` first, because that is the axis a
/// timsTOF search filters on. `MS:1002476|ion mobility drift time` is the
/// drift-tube spelling both vendor exports happen to use.
fn mobility<M: mzcore::chemistry::MassOutputMode>(spectrum: &AnnotatedSpectrum<M>) -> f32 {
    let Some(scan) = spectrum.description.acquisition.scans.first() else {
        return 0.0;
    };
    let read = |accession: u32| {
        scan.params
            .iter()
            .flat_map(|p| p.iter())
            .find(|p| is_ms_term(p, accession))
            .and_then(|p| p.value.to_f32().ok())
    };
    read(INVERSE_REDUCED_MOBILITY)
        .or_else(|| read(ION_MOBILITY_DRIFT_TIME))
        .unwrap_or(0.0)
}

/// Whether the library shipped this entry as a decoy.
///
/// A library naming `MS:1003195` reads as a decoy whether the writer meant
/// SpectraST's old spelling of shuffle-and-reposition or the current
/// unnatural-peptidoform one, because membership comes from the `is_a` closure
/// rather than from the term's name.
fn is_decoy<M: mzcore::chemistry::MassOutputMode>(
    spectrum: &AnnotatedSpectrum<M>,
    degradation: &mut Degradation,
) -> bool {
    for accession in origin_types(spectrum) {
        match subsumes_decoy(accession) {
            Some(decoy) => return decoy,
            None => {
                degradation.unrecognised_origin_types.insert(accession);
            }
        }
    }
    false
}

/// The `spectrum origin type` accessions on this entry.
///
/// The term is not repeatable, so more than one is a malformed file rather than
/// a case to reconcile; the callers take the first one they can resolve.
///
/// The value arrives as the stringified term, `MS:1003195|name`, so only the
/// CURIE half is parsed. A value in another namespace yields nothing rather than
/// a bare number: this is one of the two namespace boundaries described at the
/// top of this file.
fn origin_types<M: mzcore::chemistry::MassOutputMode>(
    spectrum: &AnnotatedSpectrum<M>,
) -> impl Iterator<Item = u32> + '_ {
    spectrum
        .description
        .params
        .iter()
        .filter(|param| is_ms_term(param, SPECTRUM_ORIGIN_TYPE))
        .filter_map(|param| {
            let value = param.value.to_string();
            let (curie, _name) = value.split_once('|')?;
            curie.strip_prefix("MS:")?.parse().ok()
        })
}

/// Whether an mzdata param is the named PSI-MS term.
///
/// mzannotate flattens a recognised attribute into an mzdata `Param`, which
/// splits the CURIE into `Option<u32>` plus `Option<ControlledVocabulary>`, so
/// there is no `Curie` left to compare in one go.
fn is_ms_term(param: &mzannotate::mzdata::params::Param, accession: u32) -> bool {
    param.controlled_vocabulary == Some(ControlledVocabulary::MS)
        && param.accession == Some(accession)
}

/// Which entries compete. msspeculator's project-defined pair id first, then
/// SpectraST's `related spectrum keys`.
///
/// Both name the target a decoy was derived from, so the target and its decoy
/// land in one group. `None` leaves the row as its own group.
fn declared_group<M: mzcore::chemistry::MassOutputMode>(
    spectrum: &AnnotatedSpectrum<M>,
) -> Option<String> {
    if let Some(param) = spectrum
        .description
        .params
        .iter()
        .find(|p| p.name == MSSPECULATOR_PAIR_ID)
    {
        return Some(param.value.to_string());
    }
    spectrum
        .attributes
        .iter()
        .flatten()
        .find(|a| a.name.accession == mzcv::curie!(MS:1003259))
        .map(|a| match &a.value {
            AttributeValue::Scalar(v) => v.to_string(),
            other => other.to_string(),
        })
}

/// The theoretical m/z for an annotated peak.
///
/// From the composition when the annotation carries one. When it does not, the
/// mzPAF mass-error suffix says how far the observed peak sits from the
/// theoretical mass, so subtracting it recovers the same number. Failing both,
/// a spectrum that declares its m/z values calculated is already reporting the
/// theoretical mass, which is the case for every predicted library including
/// the ones this project writes. An observed spectrum with neither has no
/// theoretical mass to offer and the peak is skipped.
fn theoretical_mz<M: mzcore::chemistry::MassOutputMode>(
    annotation: &Fragment<M>,
    observed: f64,
    peak_mz_is_theoretical: bool,
) -> Option<f64> {
    if let Some(mz) = annotation.mz(MassMode::Monoisotopic) {
        return Some(mz.value);
    }
    if let Some(deviation) = annotation.deviation.as_ref() {
        return Some(match deviation {
            mzcore::quantities::Tolerance::Absolute(offset) => observed - offset.value,
            mzcore::quantities::Tolerance::Relative(ppm) => observed / (1.0 + ppm.value / 1e6),
        });
    }
    peak_mz_is_theoretical.then_some(observed)
}

/// Map an annotation onto the packed label, or `None` when this build cannot
/// spell it.
///
/// Structural rather than a string round trip: the series letter, ordinal,
/// charge and isotope come off the typed annotation. Only the neutral loss goes
/// through a formula, and that lookup is keyed by composition, so `CH3SOH` and
/// `CH4OS` resolve to one loss rather than two labels for one ion.
fn to_ion_annot<M: mzcore::chemistry::MassOutputMode>(
    annotation: &Fragment<M>,
) -> Option<IonAnnot> {
    let series = match &annotation.ion {
        FragmentType::a(pos, 0) => ('a', pos.series_number),
        FragmentType::b(pos, 0) => ('b', pos.series_number),
        FragmentType::c(pos, 0) => ('c', pos.series_number),
        FragmentType::x(pos, 0) => ('x', pos.series_number),
        FragmentType::y(pos, 0) => ('y', pos.series_number),
        FragmentType::z(pos, 0) => ('z', pos.series_number),
        FragmentType::Immonium(_, residue) => {
            let ordinal =
                IonSeriesOrdinal::try_immonium(residue.aminoacid.aminoacid().one_letter_code()?)
                    .ok()?;
            return IonAnnot::new(
                ordinal,
                neutral_loss(annotation)?,
                charge_of(annotation)?,
                isotope_of(annotation)?,
            )
            .ok();
        }
        FragmentType::Precursor => {
            return IonAnnot::new(
                IonSeriesOrdinal::precursor,
                neutral_loss(annotation)?,
                charge_of(annotation)?,
                isotope_of(annotation)?,
            )
            .ok();
        }
        _ => return None,
    };
    let ordinal = u8::try_from(series.1).ok()?;
    IonAnnot::new(
        IonSeriesOrdinal::from_series_char(series.0, Some(ordinal)).ok()?,
        neutral_loss(annotation)?,
        charge_of(annotation)?,
        isotope_of(annotation)?,
    )
    .ok()
}

/// The annotation's loss, or `None` when it carries one this build has no
/// discriminant for. `None` routes the peak to an unknown label rather than
/// silently dropping the loss, which would put a lossy peak's m/z on the
/// unlossy label.
fn neutral_loss<M: mzcore::chemistry::MassOutputMode>(
    annotation: &Fragment<M>,
) -> Option<NeutralLoss> {
    match annotation.neutral_loss.as_slice() {
        [] => Some(NeutralLoss::None),
        [loss] => {
            let formula = loss.to_string();
            let formula = formula.strip_prefix('-').unwrap_or(&formula);
            NeutralLoss::from_expression(formula).ok().flatten()
        }
        _ => None,
    }
}

fn charge_of<M: mzcore::chemistry::MassOutputMode>(annotation: &Fragment<M>) -> Option<i8> {
    i8::try_from(annotation.charge.value)
        .ok()
        .filter(|c| *c != 0)
}

fn isotope_of<M: mzcore::chemistry::MassOutputMode>(annotation: &Fragment<M>) -> Option<i8> {
    match annotation.isotope.as_slice() {
        [] => Some(0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::TargetColumns;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/mzspeclib_files")
            .join(name)
    }

    /// Read through the public funnel, so anything `seal` drops fails here the
    /// way it would for a caller.
    fn arena(name: &str) -> TargetColumns<IonAnnot> {
        match super::super::library_file::read_targets(fixture(name)).expect("fixture loads") {
            TargetTable::Mzpaf { geom, .. } => geom,
            TargetTable::Str { .. } => panic!("mzSpecLib carries ion chemistry"),
        }
    }

    /// Five targets and five decoys, where the only thing marking a decoy is
    /// `MS:1003212|library attribute set name=DECOY` pointing at a header set
    /// that carries the origin type. Before the upstream attribute-set fix this
    /// file read as ten targets, so the count is the assertion that matters.
    #[test]
    fn decoys_declared_only_through_an_attribute_set_are_decoys() {
        let geom = arena("target_decoy_attribute_set.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 10);
        let decoys = geom.rows().filter(|r| geom.is_decoy(*r)).count();
        assert_eq!(decoys, 5, "five shipped decoys, declared only by set name");
    }

    /// `MS:1003259|related spectrum keys` names the target a decoy came from,
    /// so the two share a group and compete. Five pairs, ten rows.
    #[test]
    fn a_declared_pair_becomes_one_decoy_group() {
        let geom = arena("target_decoy_attribute_set.mzspeclib.txt");
        let groups: std::collections::HashSet<_> =
            geom.rows().map(|r| geom.decoy_group_code(r)).collect();
        assert_eq!(groups.len(), 5, "each decoy competes with its own target");
    }

    /// The two exports disagree on nearly every term, so reading both is what
    /// exercises the ladders rather than one file's habits.
    #[test]
    fn the_two_vendor_exports_disagree_and_both_load() {
        // DIA-NN writes `MS:1000744|selected ion m/z` and no retention time.
        let diann = arena("diann_export.mzspeclib.txt");
        assert!(diann.n_rows() > 0);
        let first = diann.rows().next().unwrap();
        assert!((diann.precursor_mz(first) - 778.412_96).abs() < 1e-4);
        assert_eq!(diann.charge(first), 2);
        assert_eq!(diann.rt_seconds(first), 0.0, "DIA-NN declares no RT");
        assert_eq!(diann.output_id(first).to_string(), "AAAAAAAAAAAAAAAASAGGK2");

        // Spectronaut writes `MS:1003208|experimental precursor monoisotopic
        // m/z`, which mzannotate routes to the isolation window rather than to
        // the selected ion, so the DIA-NN field reads zero here.
        let spectronaut = arena("spectronaut_export.mzspeclib.txt");
        let first = spectronaut.rows().next().unwrap();
        assert!((spectronaut.precursor_mz(first) - 405.763_44).abs() < 1e-3);
        assert_eq!(spectronaut.charge(first), 2);
        assert!(
            (spectronaut.mobility(first) - 0.762_965_5).abs() < 1e-6,
            "ion mobility drift time, the spelling both exports use"
        );
        assert!(
            spectronaut.rt_seconds(first) > 0.0,
            "Spectronaut declares RT"
        );
    }

    /// Fragment labels carry ion chemistry rather than an opaque name, and the
    /// m/z stored beside them is theoretical.
    #[test]
    fn fragment_labels_carry_ion_chemistry() {
        let geom = arena("diann_export.mzspeclib.txt");
        let first = geom.rows().next().unwrap();
        let labels: Vec<String> = geom
            .frag_labels(first)
            .iter()
            .map(ToString::to_string)
            .collect();
        assert!(labels.contains(&"b6".to_string()), "got {labels:?}");
        assert!(labels.contains(&"y5".to_string()), "got {labels:?}");

        let b6 = geom
            .frag_labels(first)
            .iter()
            .position(|l| l.to_string() == "b6")
            .unwrap();
        assert!(
            (geom.frag_mzs(first)[b6] - 427.229_95).abs() < 1e-3,
            "b6 at its theoretical m/z"
        );
    }

    /// Compression is the form a stored library is meant to keep, so a gzipped
    /// library has to be indistinguishable from its plain twin.
    #[test]
    fn a_gzipped_library_produces_the_same_arena() {
        let plain = arena("diann_export.mzspeclib.txt");
        let gzipped = arena("diann_export.mzspeclib.txt.gz");
        assert_eq!(plain.n_rows(), gzipped.n_rows());
        for (a, b) in plain.rows().zip(gzipped.rows()) {
            assert_eq!(
                plain.output_id(a).to_string(),
                gzipped.output_id(b).to_string()
            );
            assert_eq!(plain.precursor_mz(a), gzipped.precursor_mz(b));
            assert_eq!(plain.charge(a), gzipped.charge(b));
            assert_eq!(plain.rt_seconds(a), gzipped.rt_seconds(b));
            assert_eq!(plain.mobility(a), gzipped.mobility(b));
            assert_eq!(plain.frag_mzs(a), gzipped.frag_mzs(b));
            assert_eq!(plain.frag_labels(a), gzipped.frag_labels(b));
        }
    }

    /// Reading an undefined set as all-targets produces an FDR number that is
    /// wrong rather than absent, which is the worse failure.
    #[test]
    fn an_undefined_attribute_set_is_an_error() {
        let Err(err) = super::super::library_file::read_targets(fixture(
            "undefined_attribute_set.mzspeclib.txt",
        )) else {
            panic!("an unresolvable set claim must not load");
        };
        let msg = format!("{err:?}");
        assert!(
            msg.contains("DECOY"),
            "names the set it could not resolve: {msg}"
        );
    }

    /// The format this project writes, read back by the reader that has to read
    /// it. Built by `timsseek build-library` from a one-protein FASTA, so it
    /// exercises every term msspeculator actually emits rather than the terms a
    /// vendor happens to use: `MS:1002815|inverse reduced ion mobility` for the
    /// timsTOF axis, `MS:1000896` for retention, and the project-defined
    /// `msspeculator:decoy_pair_id` for the target/decoy pairing.
    #[test]
    fn the_format_this_project_writes_reads_back() {
        let geom = arena("msspeculator_built.mzspeclib.txt.gz");
        assert_eq!(geom.n_rows(), 4, "two targets and their two decoys");
        assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 2);

        let groups: std::collections::HashSet<_> =
            geom.rows().map(|r| geom.decoy_group_code(r)).collect();
        assert_eq!(
            groups.len(),
            2,
            "`decoy_pair_id` puts each decoy in its target's group"
        );

        let first = geom.rows().next().unwrap();
        assert_eq!(geom.output_id(first).to_string(), "VLSAAKPEDR/2");
        assert_eq!(geom.charge(first), 2);
        assert!((geom.precursor_mz(first) - 543.301_11).abs() < 1e-4);
        assert!(
            (geom.mobility(first) - 1.182_379_6).abs() < 1e-6,
            "the timsTOF mobility axis, which neither vendor export uses"
        );
        assert_eq!(geom.seq_strip(first), "VLSAAKPEDR");
        assert_eq!(geom.frag_labels(first).len(), 4);
    }

    /// A third-party library carries no `msspeculator:` attributes and no
    /// pairing. It loads, with every row its own group.
    #[test]
    fn a_library_without_provenance_loads() {
        let geom = arena("spectronaut_export.mzspeclib.txt");
        assert!(geom.n_rows() > 0);
        let groups: std::collections::HashSet<_> =
            geom.rows().map(|r| geom.decoy_group_code(r)).collect();
        assert_eq!(
            groups.len(),
            geom.n_rows(),
            "no pairing declared, so a row competes only with its own variants"
        );
    }

    /// Every entry in this fixture declares `MS:1000894|retention time` in
    /// `UO:0000010|second`, so every entry has to come back with the seconds it
    /// declared.
    ///
    /// Six of ten used to read `0.0`. mzannotate gated its retention-time branch
    /// on the attribute group holding exactly two members, and SpectraST groups a
    /// consensus retention time with `MS:1003174|attribute maximum` and
    /// `MS:1003175|attribute minimum`, so whether an entry kept its retention
    /// time depended on what else its writer put in the group. Fixed upstream.
    ///
    /// Asserting the declared values rather than a distinctness property is what
    /// makes that failure impossible to pass: six entries collapsing onto one
    /// value satisfies "more than one distinct value".
    #[test]
    fn every_declared_retention_time_survives_the_load() {
        let geom = arena("target_decoy_attribute_set.mzspeclib.txt");
        let rts: Vec<f32> = geom.rows().map(|r| geom.rt_seconds(r)).collect();

        // Read straight off the fixture, target and decoy of each pair sharing
        // one value.
        let declared = [
            1189.6, 1189.6, 2622.1, 2622.1, 4164.1, 4164.1, 1566.4, 1566.4, 1175.3, 1175.3,
        ];
        assert_eq!(rts.len(), declared.len());
        for (row, (got, want)) in rts.iter().zip(declared).enumerate() {
            assert!(
                (got - want).abs() < 0.1,
                "row {row} declares {want} s, got {got}"
            );
        }
    }

    /// `binary_search` is only correct on a sorted slice, and the generator is
    /// what sorts them. If a hand edit or a generator change breaks the order,
    /// lookups start missing terms that are present, which is silent.
    #[test]
    fn the_generated_lists_are_sorted() {
        for (name, list) in [
            ("ALL", psims_origin_type::ALL),
            ("DECOY", psims_origin_type::DECOY),
            ("THEORETICAL_MZ", psims_origin_type::THEORETICAL_MZ),
        ] {
            assert!(list.is_sorted(), "{name} is not sorted: {list:?}");
            assert!(!list.is_empty(), "{name} is empty");
        }
    }

    /// Every subset has to be a subset. A term in `DECOY` but not in `ALL`
    /// would be unreachable, because `subsumes_decoy` gates on `ALL` first.
    #[test]
    fn every_subset_is_contained_in_the_whole_subtree() {
        for accession in psims_origin_type::DECOY
            .iter()
            .chain(psims_origin_type::THEORETICAL_MZ)
        {
            assert!(
                psims_origin_type::ALL.binary_search(accession).is_ok(),
                "MS:{accession} is in a subset but not in ALL"
            );
        }
    }

    /// Membership comes from the `is_a` closure the generator resolved, so a
    /// subtype is a decoy without being named here, and a sibling of `decoy
    /// spectrum` is not.
    #[test]
    fn decoy_subtypes_subsume_and_siblings_do_not() {
        assert_eq!(subsumes_decoy(1_003_192), Some(true), "decoy spectrum");
        assert_eq!(subsumes_decoy(1_003_193), Some(true), "a subtype");
        assert_eq!(subsumes_decoy(1_003_195), Some(true), "another subtype");
        assert_eq!(subsumes_decoy(1_003_073), Some(false), "observed");
        assert_eq!(subsumes_decoy(1_003_074), Some(false), "predicted");
        assert_eq!(
            subsumes_decoy(1_003_424),
            Some(false),
            "Spectronaut's origin type"
        );
    }

    /// A term outside the subtree is unrecognised, not a target. The caller
    /// reports it rather than deciding silently.
    #[test]
    fn a_term_outside_the_subtree_does_not_resolve() {
        assert_eq!(subsumes_decoy(9_999_999), None);
        assert_eq!(subsumes_decoy(0), None);
    }

    /// Both terms whose definitions state the m/z is calculated, and one whose
    /// definition states the opposite.
    #[test]
    fn only_calculated_mz_origin_types_declare_theoretical_mz() {
        assert!(declares_theoretical_mz(1_003_074), "predicted spectrum");
        assert!(declares_theoretical_mz(1_003_424), "Spectronaut's");
        assert!(!declares_theoretical_mz(1_003_073), "observed spectrum");
    }

    /// The namespace boundary. Comparing accessions as integers is correct only
    /// after the vocabulary is established, so the same number in another
    /// vocabulary must not match. UO reaches 1010060, above every accession
    /// this file names, so the ranges genuinely overlap.
    #[test]
    fn an_accession_in_another_vocabulary_is_not_an_ms_term() {
        use mzannotate::mzdata::params::{
            Param,
            Value,
        };

        let term = |cv: Option<ControlledVocabulary>| Param {
            name: "spectrum origin type".to_string(),
            value: Value::String("MS:1003192|decoy spectrum".to_string()),
            accession: Some(SPECTRUM_ORIGIN_TYPE),
            controlled_vocabulary: cv,
            unit: Default::default(),
        };

        assert!(is_ms_term(
            &term(Some(ControlledVocabulary::MS)),
            SPECTRUM_ORIGIN_TYPE
        ));
        assert!(
            !is_ms_term(&term(Some(ControlledVocabulary::UO)), SPECTRUM_ORIGIN_TYPE),
            "the same number in the Unit Ontology is a different term"
        );
        assert!(
            !is_ms_term(&term(None), SPECTRUM_ORIGIN_TYPE),
            "a param with no vocabulary names no controlled term"
        );
    }
}
