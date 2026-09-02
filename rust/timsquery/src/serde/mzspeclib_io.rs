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
//! a peak annotated several ways at once is skipped, having no one theoretical
//! mass to be stored at.
//!
//! An unannotated peak is the case the caller decides, through
//! [`UnannotatedPeaks`]: an arena mixing observed and theoretical masses is
//! invisible downstream, so storing one at its observed m/z is a decision stated
//! before the read rather than a fallback taken during it. The decision is only
//! ever reached for a library that annotates nothing, because mzPAF's `?` -- the
//! writer marking a peak noise -- arrives as the same empty annotation list a
//! file with no annotation column produces. See [`AnnotationCensus`].
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
    Attribute,
    AttributeValue,
    EntryType,
    LibraryHeader,
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
use crate::models::capabilities::{
    LoadPolicy,
    UnannotatedPeaks,
};
use crate::models::{
    Row,
    TargetCapabilities,
    TargetColumns,
};

/// What the parser yields. `MzSpecLibTextParser`'s `Iterator::Item` fixes the
/// mass-output mode, so nothing downstream of it is generic over one.
type Spectrum = AnnotatedSpectrum<mzcore::chemistry::OutputMolecularFormula>;

/// One peak annotation, in the same fixed mode.
type Annotation = Fragment<mzcore::chemistry::OutputMolecularFormula>;

// ── The controlled vocabulary ────────────────────────────────────────────────
//
// Every bare `u32` accession in this file is a PSI-MS accession, and comparing
// them as integers is only correct because the namespace is established first.
// There are exactly two places that establish it, and every integer downstream
// of them has already passed through one:
//
//   - `is_ms_term`, for values mzannotate flattened into an mzdata `Param`,
//     which splits the CURIE into a number and a namespace.
//   - the `MS:` prefix check in `ms_accession`, for a term arriving as text.
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
/// `MS:1003065|spectrum aggregation type`, the term's other parent.
const SPECTRUM_AGGREGATION_TYPE: u32 = 1_003_065;
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
    psims_origin_type::lookup(accession).map(|term| term.decoy)
}

/// Whether this spectrum's peak m/z values are calculated rather than measured.
///
/// It decides what an annotation with no composition and no mass-error suffix
/// means. On a predicted spectrum the peak m/z already *is* the theoretical
/// mass, so the peak is usable; on an observed one there is no theoretical mass
/// to be had and the peak is skipped rather than mixed in at its observed m/z.
fn declares_theoretical_mz(accession: u32) -> bool {
    psims_origin_type::lookup(accession).is_some_and(|term| term.theoretical_mz)
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
pub fn read_mzspeclib_library_file(
    path: &Path,
    policy: LoadPolicy,
) -> Result<TargetTable, TargetReadingError> {
    read_counting_degradation(path, policy).map(|(table, _)| table)
}

/// The read itself, with the counts still in hand.
///
/// Split out so a test can assert on them: a count that is only ever warned
/// about is a count nothing checks.
fn read_counting_degradation(
    path: &Path,
    policy: LoadPolicy,
) -> Result<(TargetTable, Degradation), TargetReadingError> {
    let reader = open_reader(path)?;
    let parser = MzSpecLibTextParser::open(reader, Some(path.to_path_buf()), ontologies())
        .map_err(|e| TargetReadingError::SpeclibParse(format!("mzSpecLib header: {e}")))?;

    let mut geom = TargetColumns::with_capabilities(TargetCapabilities::default_diann());
    let mut frag_intens: Vec<f32> = Vec::new();
    let mut degradation = Degradation::default();
    // Read off the header before the iterator takes the parser, which is also
    // the only place it can be read: it is what an entry declaring nothing about
    // its mass provenance falls back to.
    let library_theoretical_mz = library_declares_theoretical_mz(parser.header());
    let mut census = AnnotationCensus::new(policy.unannotated);

    for (index, spectrum) in parser.enumerate() {
        let spectrum = spectrum.map_err(|e| {
            TargetReadingError::SpeclibParse(format!("mzSpecLib spectrum {}: {e}", index + 1))
        })?;
        let unannotated = census.plan_for(&spectrum)?;
        let row = SpectrumRow::extract(
            &spectrum,
            library_theoretical_mz,
            unannotated,
            &mut degradation,
        )?;
        // Dropped here, so the decoy's peaks are never pushed and its group
        // never interns.
        if !policy.decoys.accepts(row.is_decoy) {
            continue;
        }
        // Counted past that filter rather than inside the extract, because the
        // warning divides by the sealed row count and a dropped decoy is in
        // neither set. An entry the file gave peaks for and the arena got none
        // of is scored against nothing, which the per-peak counts do not say on
        // their own: they cannot distinguish a library losing a peak here and
        // there from one losing whole entries.
        if row.frags.is_empty() && !spectrum.peaks.is_empty() {
            degradation.rows_without_fragments += 1;
        }

        // Both halves of a pair have to land in one namespace. A decoy names its
        // target by `<Spectrum=N>` key, so every row names its own key and a
        // decoy overrides that with its target's -- not its source id, which is
        // a different string entirely. A library that declares no pair therefore
        // ends up with one key per row, all singletons, and `seal` drops the
        // column rather than storing a group per row.
        let group = row
            .declared_group
            .unwrap_or_else(|| row.key.to_string())
            .into();

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
            id: Some(row.id.into()),
            decoy_group: Some(group),
            ..Default::default()
        });
    }

    let geom = geom.seal(policy.decoys)?;
    degradation.report(geom.n_rows());

    Ok((
        TargetTable::Mzpaf {
            geom,
            frag_intens: Some(frag_intens),
        },
        degradation,
    ))
}

/// Whether this library annotates its peaks at all, decided from the first
/// entry that has any.
///
/// The decision is the library's rather than the peak's because the peak cannot
/// carry it: mzannotate collapses mzPAF's `?` onto the same empty annotation
/// list a bare peak list produces, so a noise-marked peak and an unannotated one
/// are one shape here. An entry that annotates something has marked the rest of
/// its peaks noise; a library that annotates nothing is the library the caller's
/// [`UnannotatedPeaks`] was stated for.
///
/// A leading entry with no peaks decides nothing, or a library whose first entry
/// happens to declare none would have that entry answer for all the rest.
struct AnnotationCensus {
    /// The caller's decision, reached only once the library turns out to
    /// annotate nothing.
    policy: UnannotatedPeaks,
    annotates: Option<bool>,
}

impl AnnotationCensus {
    fn new(policy: UnannotatedPeaks) -> Self {
        Self {
            policy,
            annotates: None,
        }
    }

    /// What this entry's unannotated peaks get.
    ///
    /// [`UnannotatedPeaks::KeepAll`] discards the annotations an entry does
    /// carry, so it never asks the question and the census stays undecided
    /// under it.
    fn plan_for(&mut self, spectrum: &Spectrum) -> Result<UnannotatedPeaks, TargetReadingError> {
        if self.policy == UnannotatedPeaks::KeepAll {
            return Ok(UnannotatedPeaks::KeepAll);
        }
        let annotates = spectrum
            .peaks
            .iter()
            .any(|peak| !peak.annotations.is_empty());
        match self.annotates {
            // Nothing in the corpus does this. It is refused rather than
            // resolved because the two answers are not reconcilable: the
            // earlier entries' unannotated peaks were read as the only peaks
            // the library has, and this entry says they were noise.
            Some(false) if annotates => {
                return Err(TargetReadingError::SpeclibParse(format!(
                    "spectrum {} ({}) carries an annotated peak, in a library whose earlier \
                     entries annotate none; its unannotated peaks have already been read as \
                     the only peaks it has",
                    spectrum.key,
                    source_id(spectrum),
                )));
            }
            None if !spectrum.peaks.is_empty() => self.annotates = Some(annotates),
            _ => {}
        }
        if self.annotates == Some(true) {
            // An annotated peak anywhere in the library makes an unannotated
            // one mzPAF's `?`, which is the writer calling it noise.
            return Ok(UnannotatedPeaks::Skip);
        }
        Ok(self.policy)
    }
}

/// What a library lost on the way into the arena, counted rather than logged
/// per row: a library with a million unannotated peaks should produce one line,
/// not a million.
#[derive(Default)]
struct Degradation {
    /// Every peak the file declared, so the peak-level counts below have a
    /// denominator the way the row-level ones have the sealed row count.
    total_peaks: usize,
    unannotated_peaks: usize,
    ambiguous_peaks: usize,
    unrepresentable_labels: usize,
    peaks_without_theoretical_mz: usize,
    /// Peaks stored under a placeholder label at the m/z the file measured,
    /// which is the one way an observed mass reaches the arena.
    peaks_kept_at_observed_mz: usize,
    /// Peaks past the unknown labels an entry has to number them with, dropped
    /// quietest first. See [`loudest_within_the_label_ceiling`].
    peaks_over_the_label_ceiling: usize,
    rows_without_sequence: usize,
    rows_without_fragments: usize,
    // No count of entries missing a retention time. Zero is an ordinary value
    // on a normalized scale rather than an absence -- it is an anchor point, and
    // it is what msspeculator's own libraries write -- so there is nothing to
    // test for. See the CONTEXT.md entry on retention time.
    /// A set, because this is reached per spectrum rather than per file and a
    /// library repeats one unrecognised term on every row.
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
        // Both reasons land on a `?N` label and are indistinguishable
        // downstream, so the split is stated here or nowhere: one is a peak at
        // its theoretical m/z whose ion this build cannot spell, the other a
        // peak at the m/z the file measured.
        let unknown_labels = self.unrepresentable_labels + self.peaks_kept_at_observed_mz;
        if unknown_labels > 0 {
            warn!(
                "mzSpecLib: {unknown_labels} of {} peak(s) reached the arena under an unknown \
                 (`?N`) label -- {} unannotated, {} carrying an ion this build cannot spell",
                self.total_peaks, self.peaks_kept_at_observed_mz, self.unrepresentable_labels,
            );
        }
        if self.peaks_over_the_label_ceiling > 0 {
            warn!(
                "mzSpecLib: {} of {} peak(s) are past the {UNKNOWN_LABELS_PER_ENTRY} unknown \
                 labels one entry can carry and were dropped, quietest first",
                self.peaks_over_the_label_ceiling, self.total_peaks,
            );
        }
        if self.peaks_kept_at_observed_mz > 0 {
            warn!(
                "mzSpecLib: {} of {} peak(s) are stored at the m/z the file measured rather than \
                 a theoretical one; calibration fits its m/z offsets assuming a library mass is \
                 theoretical, so for these it fits the library's own error along with the run's",
                self.peaks_kept_at_observed_mz, self.total_peaks,
            );
        }
        if self.rows_without_sequence > 0 {
            warn!(
                "mzSpecLib: {} of {rows} entries carry no peptidoform; sequence features are unavailable for them",
                self.rows_without_sequence
            );
        }
        if self.rows_without_fragments > 0 {
            warn!(
                "mzSpecLib: {} of {rows} entries declare peaks none of which could be read; they score against nothing",
                self.rows_without_fragments
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
    /// Parallel to `intensities`, and pushed straight into the arena.
    frags: Vec<(IonAnnot, f64)>,
    /// The reference-intensity sidecar for `frags`, kept separate because the
    /// arena stores the two in different columns.
    intensities: Vec<f32>,
}

impl SpectrumRow {
    fn extract(
        spectrum: &Spectrum,
        library_declares_theoretical_mz: bool,
        unannotated: UnannotatedPeaks,
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

        let peak_mz_is_theoretical =
            peak_mz_is_theoretical(spectrum, library_declares_theoretical_mz);
        let mut counter = UnknownIonCounter::default();
        let mut frags = Vec::with_capacity(spectrum.peaks.len());
        let mut intensities = Vec::with_capacity(spectrum.peaks.len());
        degradation.total_peaks += spectrum.peaks.len();
        let under_ceiling = loudest_within_the_label_ceiling(spectrum, unannotated);
        for (index, peak) in spectrum.peaks.iter().enumerate() {
            if under_ceiling.as_ref().is_some_and(|keep| !keep[index]) {
                degradation.peaks_over_the_label_ceiling += 1;
                continue;
            }
            let observed = peak.mz.value;
            // `KeepAll` is uniform by construction: it drops the annotations an
            // entry does carry rather than storing some of its peaks at a
            // theoretical m/z and the rest at a measured one.
            let annotations: &[Annotation] = if unannotated == UnannotatedPeaks::KeepAll {
                &[]
            } else {
                peak.annotations.as_slice()
            };
            match annotations {
                [] => match unannotated {
                    UnannotatedPeaks::Skip => degradation.unannotated_peaks += 1,
                    UnannotatedPeaks::Fail => {
                        return Err(TargetReadingError::SpeclibParse(format!(
                            "spectrum {} ({}) has an unannotated peak at m/z {observed} and this \
                             library annotates none, so nothing in it can be given a theoretical \
                             mass",
                            spectrum.key,
                            source_id(spectrum),
                        )));
                    }
                    UnannotatedPeaks::Keep | UnannotatedPeaks::KeepAll => {
                        degradation.peaks_kept_at_observed_mz += 1;
                        // Charge 1: the charge would have come from the
                        // annotation, and this peak has none to read it off.
                        frags.push((mint_unknown(&mut counter, 1)?, observed));
                        intensities.push(peak.intensity);
                    }
                },
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
                            mint_unknown(&mut counter, annotation.charge.value.max(1) as i8)?
                        }
                    };
                    frags.push((label, mz));
                    intensities.push(peak.intensity);
                }
                // Several identities for one peak: no single theoretical mass
                // exists, so there is nothing to store.
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
            intensities,
        })
    }
}

/// A claimed attribute set the header never defined is an error.
///
/// Reading it as a target produces an FDR number that is wrong rather than
/// absent, and the decoy declaration is exactly what such a set usually carries.
/// mzannotate drops an unresolvable claim silently but leaves the claim itself
/// in the spectrum's leftover attributes, which is what this sees.
fn reject_unresolved_attribute_sets(spectrum: &Spectrum) -> Result<(), TargetReadingError> {
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
fn source_id(spectrum: &Spectrum) -> String {
    if spectrum.description.id.is_empty() {
        spectrum.key.to_string()
    } else {
        spectrum.description.id.clone()
    }
}

/// `MS:1000744|selected ion m/z` (DIA-NN), then
/// `MS:1003208|experimental precursor monoisotopic m/z` (Spectronaut, which
/// mzannotate routes to the isolation window), then the analyte's own mass.
fn precursor_mz(
    spectrum: &Spectrum,
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
fn charge(spectrum: &Spectrum, peptidoform: Option<&mzcore::sequence::PeptidoformIon>) -> u8 {
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
fn rt_seconds(spectrum: &Spectrum) -> f32 {
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
fn mobility(spectrum: &Spectrum) -> f32 {
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
///
/// One decoy subtype anywhere among the entry's origin types is enough, for the
/// reason [`origin_types`] gives: a further origin type naming how the masses
/// were produced does not withdraw the decoy claim.
fn is_decoy(spectrum: &Spectrum, degradation: &mut Degradation) -> bool {
    let mut decoy = false;
    for accession in origin_types(spectrum) {
        match subsumes_decoy(accession) {
            Some(subsumes) => decoy |= subsumes,
            None => {
                degradation.unrecognised_origin_types.insert(accession);
            }
        }
    }
    decoy
}

/// The `spectrum origin type` accessions on this entry.
///
/// Several of them are several assertions rather than a contradiction to
/// resolve. mzSpecLib v1.0 §4.1.11 gives the last instance precedence only where
/// the values differ *and cannot be reconciled*, and it permits an attribute to
/// be repeated deliberately to encode multiple values. Here they reconcile:
/// predicted-ness and decoy-ness are orthogonal facts, so an entry naming both
/// `MS:1003074|predicted spectrum` and `MS:1003195` is asserting both and
/// neither instance displaces the other. Both callers read the whole set
/// accordingly -- a decoy subtype anywhere in it makes the entry a decoy, a
/// calculated-m/z term anywhere in it makes the peak m/z theoretical.
fn origin_types(spectrum: &Spectrum) -> impl Iterator<Item = u32> + '_ {
    terms_valued_by(spectrum, SPECTRUM_ORIGIN_TYPE)
}

/// The PSI-MS accessions this entry gives as the value of `accession`.
///
/// Both containers a declaration can land in. mzannotate flattens an ungrouped
/// attribute it recognises into `description.params` but leaves a grouped one it
/// does not special-case in `attributes`, and mzSpecLib v1.0 §4.1.6 writes the
/// aggregation type grouped with its replicate counts. Reading one container
/// would resolve a term for one of its two legal spellings.
fn terms_valued_by(spectrum: &Spectrum, accession: u32) -> impl Iterator<Item = u32> + '_ {
    let ungrouped = spectrum
        .description
        .params
        .iter()
        .filter(move |param| is_ms_term(param, accession))
        .map(|param| param.value.to_string());
    let grouped = spectrum
        .attributes
        .iter()
        .flatten()
        .filter(move |attribute| is_ms_attribute(attribute, accession))
        .map(|attribute| attribute.value.to_string());
    ungrouped
        .chain(grouped)
        .filter_map(|value| ms_accession(&value))
}

/// The accession in a stringified term, `MS:1003195|name`, so only the CURIE
/// half is parsed.
///
/// A value in another namespace yields nothing rather than a bare number: this
/// is one of the two namespace boundaries described at the top of this file.
fn ms_accession(value: &str) -> Option<u32> {
    let (curie, _name) = value.split_once('|')?;
    curie.strip_prefix("MS:")?.parse().ok()
}

/// Whether this entry's peak m/z values are calculated rather than measured.
///
/// `MS:1003074|predicted spectrum` is `is_a` both `MS:1003072|spectrum origin
/// type` and `MS:1003065|spectrum aggregation type`, so both axes are consulted:
/// reading one of a term's two parents and not the other would be arbitrary.
///
/// The origin-type slot carries two orthogonal facts, how the masses were
/// produced and whether the analyte is a decoy, because the decoy subtree the
/// format's own examples use hangs off the spectrum's origin type and
/// `MS:1002217|decoy peptide` is a peptide-identification attribute rather than
/// a library-entry one. Mass provenance is not read off that subtree, which does
/// not determine one: `MS:1003195` is defined as "a decoy spectrum that is
/// either a real spectrum of an unnatural peptidoform (e.g. a synthetic peptide
/// that cannot be found in nature), or an artificial spectrum predicted for such
/// unnatural peptidoform", so the vocabulary itself has the term spanning
/// measured and predicted. An entry whose own origin types are all decoy
/// subtypes -- or that declares none, the empty case the `all` below accepts --
/// has left the question to the library's `Spectrum=all` set. mzSpecLib v1.0
/// §4.1.4 and §4.1.11 resolve attribute sets per accession with the later set
/// replacing the earlier, and §4.1.12 Example 1 shows exactly this override, so
/// the entry is not contradicting the library when it names a decoy subtype.
///
/// Anything else on the origin-type axis is an override meant as one:
/// `MS:1003073|observed spectrum` on an entry of a predicted library is that
/// entry saying its masses were measured, and it is believed.
fn peak_mz_is_theoretical(spectrum: &Spectrum, library_declares_theoretical_mz: bool) -> bool {
    if origin_types(spectrum).any(declares_theoretical_mz)
        || terms_valued_by(spectrum, SPECTRUM_AGGREGATION_TYPE).any(declares_theoretical_mz)
    {
        return true;
    }
    origin_types(spectrum).all(|accession| subsumes_decoy(accession) == Some(true))
        && library_declares_theoretical_mz
}

/// Whether the library's `Spectrum=all` set declares its entries' m/z values
/// calculated, on either axis that can carry the claim.
///
/// The header is parsed by `MzSpecLibTextParser::open`, so this is answerable
/// before the first entry is yielded.
fn library_declares_theoretical_mz(header: &LibraryHeader) -> bool {
    header
        .attribute_classes
        .get(&EntryType::Spectrum)
        .into_iter()
        .flatten()
        .filter(|set| set.id == "all")
        .flat_map(|set| set.attributes.iter().flatten())
        .filter(|attribute| {
            attribute.name.accession == mzcv::curie!(MS:1003072)
                || attribute.name.accession == mzcv::curie!(MS:1003065)
        })
        .filter_map(|attribute| ms_accession(&attribute.value.to_string()))
        .any(declares_theoretical_mz)
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

/// Whether an mzSpecLib attribute is the named PSI-MS term.
///
/// An attribute mzannotate left grouped still carries its whole `Curie`, so the
/// namespace and the accession are compared in one go rather than separately as
/// [`is_ms_term`] has to.
fn is_ms_attribute(attribute: &Attribute, accession: u32) -> bool {
    attribute.name.accession
        == mzcv::Curie {
            cv: mzcv::ControlledVocabulary::MS,
            accession: mzcv::AccessionCode::Numeric(accession),
        }
}

/// Which entries compete. msspeculator's project-defined pair id first, then
/// SpectraST's `related spectrum keys`.
///
/// Both name the target a decoy was derived from, so the target and its decoy
/// land in one group. `None` leaves the row as its own group.
fn declared_group(spectrum: &Spectrum) -> Option<String> {
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
fn theoretical_mz(
    annotation: &Annotation,
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

/// Unknown labels one entry can carry. `UnknownIonCounter` numbers them in a
/// `u8` from 1, so the 256th mint is an error rather than a dropped peak.
const UNKNOWN_LABELS_PER_ENTRY: usize = u8::MAX as usize;

/// Which of this entry's peaks fit under [`UNKNOWN_LABELS_PER_ENTRY`], or `None`
/// when every one of them does.
///
/// Only a plan that mints a label per peak can reach the ceiling, and only on a
/// spectrum with more peaks than there are labels -- nothing in the corpus comes
/// within a hundred of it, but a DDA consensus spectrum read under
/// [`UnannotatedPeaks::KeepAll`] runs to thousands. Which peaks survive is not
/// arbitrary: every score is driven by fragment intensity, so the loudest are
/// the ones worth the labels.
///
/// Sorted rather than kept in a bounded heap. An entry carries peaks in the tens
/// to hundreds, where an `O(k log k)` sort per row costs nothing and a heap is
/// more code to be wrong in.
fn loudest_within_the_label_ceiling(
    spectrum: &Spectrum,
    unannotated: UnannotatedPeaks,
) -> Option<Vec<bool>> {
    let mints_a_label_per_peak = matches!(
        unannotated,
        UnannotatedPeaks::Keep | UnannotatedPeaks::KeepAll
    );
    if !mints_a_label_per_peak || spectrum.peaks.len() <= UNKNOWN_LABELS_PER_ENTRY {
        return None;
    }
    let intensity = |index: usize| spectrum.peaks[index].intensity;
    let mut by_intensity: Vec<usize> = (0..spectrum.peaks.len()).collect();
    by_intensity.sort_unstable_by(|a, b| intensity(*b).total_cmp(&intensity(*a)));

    let mut keep = vec![false; spectrum.peaks.len()];
    for index in by_intensity.into_iter().take(UNKNOWN_LABELS_PER_ENTRY) {
        keep[index] = true;
    }
    Some(keep)
}

/// The next `?N` label for this entry.
///
/// One counter per entry, so `N` numbers the unknowns within a spectrum rather
/// than within the library, and there are 255 of them.
fn mint_unknown(
    counter: &mut UnknownIonCounter,
    charge: i8,
) -> Result<IonAnnot, TargetReadingError> {
    counter.next_unknown(charge).map_err(|e| {
        TargetReadingError::SpeclibParse(format!("ran out of unknown-ion labels: {e}"))
    })
}

/// Map an annotation onto the packed label, or `None` when this build cannot
/// spell it.
///
/// Structural rather than a string round trip: the series letter, ordinal,
/// charge and isotope come off the typed annotation. Only the neutral loss goes
/// through a formula, and that lookup is keyed by composition, so `CH3SOH` and
/// `CH4OS` resolve to one loss rather than two labels for one ion.
fn to_ion_annot(annotation: &Annotation) -> Option<IonAnnot> {
    // The match yields only the kind; the annotation is built once below.
    let backbone = |letter: char, position: &mzcore::sequence::PeptidePosition| {
        IonSeriesOrdinal::from_series_char(letter, Some(u8::try_from(position.series_number).ok()?))
            .ok()
    };
    let series = match &annotation.ion {
        FragmentType::a(pos, 0) => backbone('a', pos)?,
        FragmentType::b(pos, 0) => backbone('b', pos)?,
        FragmentType::c(pos, 0) => backbone('c', pos)?,
        FragmentType::x(pos, 0) => backbone('x', pos)?,
        FragmentType::y(pos, 0) => backbone('y', pos)?,
        FragmentType::z(pos, 0) => backbone('z', pos)?,
        FragmentType::Immonium(_, residue) => {
            IonSeriesOrdinal::try_immonium(residue.aminoacid.aminoacid().one_letter_code()?).ok()?
        }
        FragmentType::Precursor => IonSeriesOrdinal::precursor,
        _ => return None,
    };
    IonAnnot::new(
        series,
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
fn neutral_loss(annotation: &Annotation) -> Option<NeutralLoss> {
    match annotation.neutral_loss.as_slice() {
        [] => Some(NeutralLoss::None),
        // A loss, and only a loss. `Gain` and `SideChainLoss` are outside what
        // this build can spell, and are rejected on the variant rather than on
        // their spelling.
        //
        // The multiplier stays in the string: Hill notation renders
        // `Loss(2, H2O)` as `-2H2O`, which the composition lookup resolves to
        // the two-water entry. Matching only `Loss(1, ..)` would drop it.
        [loss @ mzcore::chemistry::NeutralLoss::Loss(..)] => {
            NeutralLoss::lookup_formula(&loss.to_string())
        }
        _ => None,
    }
}

fn charge_of(annotation: &Annotation) -> Option<i8> {
    i8::try_from(annotation.charge.value)
        .ok()
        .filter(|c| *c != 0)
}

fn isotope_of(annotation: &Annotation) -> Option<i8> {
    match annotation.isotope.as_slice() {
        [] => Some(0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::TargetColumns;
    use crate::models::capabilities::DecoyPolicy;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/mzspeclib_files")
            .join(name)
    }

    /// A load that decides only about decoys, leaving everything else default.
    fn deciding_decoys(decoys: DecoyPolicy) -> LoadPolicy {
        LoadPolicy {
            decoys,
            ..Default::default()
        }
    }

    /// The other half of the policy, alone.
    fn deciding_unannotated(unannotated: UnannotatedPeaks) -> LoadPolicy {
        LoadPolicy {
            unannotated,
            ..Default::default()
        }
    }

    /// Read through the public funnel, so anything `seal` drops fails here the
    /// way it would for a caller.
    fn arena(name: &str) -> TargetColumns<IonAnnot> {
        match super::super::library_file::read_targets(fixture(name)).expect("fixture loads") {
            TargetTable::Mzpaf { geom, .. } => geom,
            TargetTable::Str { .. } => panic!("mzSpecLib carries ion chemistry"),
        }
    }

    /// A library of small molecules is still a library to this reader: the rows
    /// arrive, and everything peptide-shaped about them is absent rather than
    /// wrong. No peak in it carries an annotation, so under the default every
    /// one is kept at the m/z the file measured under a `?N` label, and an
    /// analyte identified by `MS:1000866|molecular formula` yields no sequence.
    ///
    /// The five m/z values, labels and intensities are the fixture's own,
    /// unaltered: an observed m/z reaching the arena unchanged is the whole
    /// point of keeping the peak.
    #[test]
    fn a_library_that_annotates_nothing_keeps_its_peaks_at_the_observed_mz() {
        let geom = arena("small_molecule.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 1);

        let row = geom.rows().next().expect("one row");
        assert_eq!(geom.output_id(row).to_string(), "JWH-250-5OH");
        assert_eq!(geom.seq_strip(row), "");
        assert_eq!(geom.seq_mod(row), "");
        assert!(!geom.is_decoy(row), "no origin type is not a decoy claim");

        let labels: Vec<String> = geom
            .frag_labels(row)
            .iter()
            .map(ToString::to_string)
            .collect();
        assert_eq!(labels, ["?1", "?2", "?3", "?4", "?5"]);
        for label in geom.frag_labels(row) {
            assert_eq!(
                label.get_charge(),
                1,
                "an unannotated peak has no charge to read off"
            );
        }
        assert_eq!(
            geom.frag_mzs(row),
            [51.0236, 52.3338, 63.9969, 65.0390, 121.0649]
        );

        let TargetTable::Mzpaf { frag_intens, .. } = read_mzspeclib_library_file(
            &fixture("small_molecule.mzspeclib.txt"),
            LoadPolicy::default(),
        )
        .expect("fixture loads") else {
            panic!("mzSpecLib carries ion chemistry");
        };
        // The sidecar is `f32`, so the fixture's own digits are narrowed rather
        // than rewritten to what an `f32` happens to hold.
        let declared: [f64; 5] = [41921.14, 18399.55, 22028.22, 263790.31, 80814008.00];
        assert_eq!(
            frag_intens.expect("mzSpecLib fills the sidecar"),
            declared.map(|i| i as f32)
        );
    }

    /// `Skip` is the caller saying an arena of observed masses is not what they
    /// want, and it leaves this library with nothing in it. timsseek's
    /// `a_library_with_no_annotated_fragments_is_refused` is what turns that
    /// into a refusal; here it is only the count.
    #[test]
    fn skipping_unannotated_peaks_leaves_that_library_no_fragments() {
        let TargetTable::Mzpaf { geom, .. } = read_mzspeclib_library_file(
            &fixture("small_molecule.mzspeclib.txt"),
            deciding_unannotated(UnannotatedPeaks::Skip),
        )
        .expect("fixture loads") else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_rows(), 1, "the row arrives either way");
        assert_eq!(geom.n_fragments(), 0);
    }

    /// `Fail` is for a caller who would rather hear about it than search a
    /// library of observed masses, so the refusal has to say which entry it
    /// found one in.
    #[test]
    fn failing_on_an_unannotated_peak_names_the_entry_it_was_in() {
        let Err(err) = read_mzspeclib_library_file(
            &fixture("small_molecule.mzspeclib.txt"),
            deciding_unannotated(UnannotatedPeaks::Fail),
        ) else {
            panic!("a library of unannotated peaks must not load under Fail");
        };
        let msg = format!("{err:?}");
        assert!(msg.contains("JWH-250-5OH"), "names the entry: {msg}");
    }

    /// The mixed library the specification itself ships: 820 of its 1242 peaks
    /// are annotated and the other 422 are mzPAF `?`. Mixing is what a
    /// consensus library looks like, so the policy must not reach it -- the
    /// annotated peaks have already said the `?` ones are noise, and every
    /// policy but `KeepAll` has to leave the arena exactly as `Skip` does.
    #[test]
    fn a_library_that_annotates_anything_reads_the_same_under_every_policy() {
        let read = |unannotated| {
            let TargetTable::Mzpaf { geom, .. } = read_mzspeclib_library_file(
                &fixture("target_decoy_attribute_set.mzspeclib.txt"),
                deciding_unannotated(unannotated),
            )
            .expect("fixture loads") else {
                panic!("mzSpecLib carries ion chemistry");
            };
            (geom.n_rows(), geom.n_fragments())
        };
        for unannotated in [
            UnannotatedPeaks::Skip,
            UnannotatedPeaks::Fail,
            UnannotatedPeaks::Keep,
        ] {
            assert_eq!(read(unannotated), (10, 736), "under {unannotated:?}");
        }
    }

    /// `KeepAll` is the one policy that asks the library nothing: it discards
    /// the annotations every peak has, so the arena it produces is uniform
    /// whatever the file annotated. The mixed library is 1242 peaks and 1242
    /// fragments under it, against the 736 every other policy comes back with.
    #[test]
    fn keep_all_stores_every_peak_and_discards_the_annotations() {
        let TargetTable::Mzpaf { geom, frag_intens } = read_mzspeclib_library_file(
            &fixture("target_decoy_attribute_set.mzspeclib.txt"),
            deciding_unannotated(UnannotatedPeaks::KeepAll),
        )
        .expect("fixture loads") else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_rows(), 10);
        assert_eq!(geom.n_fragments(), 1242, "every peak the file declares");
        assert_eq!(
            frag_intens.expect("mzSpecLib fills the sidecar").len(),
            1242
        );
        for row in geom.rows() {
            for label in geom.frag_labels(row) {
                assert_eq!(
                    label.to_string().chars().next(),
                    Some('?'),
                    "an annotation this policy discarded left a label behind"
                );
            }
        }

        // Nothing consults the library's own annotations under this policy, so
        // the entry-wise contradiction the other policies refuse is not one
        // here: there is no decision left for an entry to contradict.
        let TargetTable::Mzpaf { geom, .. } = read_mzspeclib_library_file(
            &fixture("entry_annotates_after_one_that_does_not.mzspeclib.txt"),
            deciding_unannotated(UnannotatedPeaks::KeepAll),
        )
        .expect("the mixture is not a contradiction under KeepAll") else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_fragments(), 4, "both peaks of both entries");
    }

    /// An entry with more peaks than it has unknown labels to number them with.
    /// A DDA consensus spectrum is that shape; no committed fixture is, because
    /// 300 peak lines is not a file anyone reads in a review, so this one is
    /// written here instead.
    ///
    /// Intensity rises with m/z, which makes the 255 survivors the top 255 by
    /// m/z as well and the count assertion a claim about which peaks were kept.
    #[test]
    fn an_entry_over_the_label_ceiling_keeps_only_its_loudest_peaks() {
        const PEAKS: usize = 300;
        let mut library = String::from(
            "<mzSpecLib>\n\
             MS:1003186|library format version=1.0\n\
             MS:1003188|library name=over_the_label_ceiling\n\
             <Spectrum=1>\n\
             MS:1003061|library spectrum name=CONSENSUS\n\
             MS:1000744|selected ion m/z=700.0\n",
        );
        library.push_str(&format!("MS:1003059|number of peaks={PEAKS}\n"));
        library.push_str("<Analyte=1>\nMS:1000041|charge state=1\n<Peaks>\n");
        for peak in 0..PEAKS {
            library.push_str(&format!("{}\t{}.0\n", 100.0 + peak as f64, peak + 1));
        }
        let dir = tempfile::tempdir().expect("a temp dir");
        let path = dir.path().join("over_the_label_ceiling.mzspeclib.txt");
        std::fs::write(&path, library).expect("the generated library is written");

        // The ceiling belongs to the label, not to the policy: `Keep` mints one
        // label per peak for a library that annotates nothing, so a fully
        // unannotated entry reaches it exactly as `KeepAll` does.
        for unannotated in [UnannotatedPeaks::Keep, UnannotatedPeaks::KeepAll] {
            let (table, degradation) =
                read_counting_degradation(&path, deciding_unannotated(unannotated))
                    .expect("an entry over the ceiling loads rather than failing");
            let TargetTable::Mzpaf { geom, frag_intens } = table else {
                panic!("mzSpecLib carries ion chemistry");
            };
            let row = geom.rows().next().expect("one row");
            assert_eq!(
                geom.frag_labels(row).len(),
                UNKNOWN_LABELS_PER_ENTRY,
                "under {unannotated:?}"
            );
            assert_eq!(degradation.total_peaks, PEAKS);
            assert_eq!(
                degradation.peaks_over_the_label_ceiling,
                PEAKS - UNKNOWN_LABELS_PER_ENTRY
            );

            let quietest_kept = 1.0 + (PEAKS - UNKNOWN_LABELS_PER_ENTRY) as f32;
            let intens = frag_intens.expect("mzSpecLib fills the sidecar");
            assert_eq!(
                intens.iter().copied().fold(f32::INFINITY, f32::min),
                quietest_kept,
                "the 45 quietest peaks are the ones dropped"
            );
            assert_eq!(
                geom.frag_mzs(row).first().copied(),
                Some(100.0 + (PEAKS - UNKNOWN_LABELS_PER_ENTRY) as f64),
                "the survivors stay in the order the file wrote them"
            );
        }
    }

    /// The library-level decision waits for an entry that has peaks. This
    /// library's first entry declares none, and its second annotates both of
    /// its own: an entry with nothing to say deciding for the library would
    /// make the second entry a contradiction and refuse the file.
    #[test]
    fn an_entry_with_no_peaks_decides_nothing_about_the_library() {
        let geom = arena("first_entry_declares_no_peaks.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 2);
        let (empty, annotated) = (geom.rows().next().unwrap(), geom.rows().nth(1).unwrap());
        assert_eq!(geom.frag_labels(empty).len(), 0, "it declared no peaks");
        assert_eq!(geom.frag_labels(annotated).len(), 2);
    }

    /// Nothing in the corpus mixes the two entry-wise, and the reader would
    /// have read one library two ways if it did: the first entry's peaks stored
    /// as the only peaks it has, the second's marked noise and dropped.
    #[test]
    fn an_annotated_entry_after_an_unannotated_one_is_refused() {
        let Err(err) = read_mzspeclib_library_file(
            &fixture("entry_annotates_after_one_that_does_not.mzspeclib.txt"),
            LoadPolicy::default(),
        ) else {
            panic!("a library that annotates only some of its entries must not load");
        };
        let msg = format!("{err:?}");
        assert!(msg.contains("SAMPLERK/2"), "names the entry: {msg}");
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

    /// `DecoyPolicy::Force` drops shipped decoys at the row level, which is
    /// what makes `--decoy-strategy force` mean anything: with none in the
    /// arena, `seal()` has no reason to rewrite `MassShift` to `Stored`, so the
    /// mass-shift decoys the user asked for are the ones
    /// scored.
    #[test]
    fn skipping_shipped_decoys_leaves_only_targets() {
        let path = fixture("target_decoy_attribute_set.mzspeclib.txt");
        let TargetTable::Mzpaf { geom, .. } =
            read_mzspeclib_library_file(&path, deciding_decoys(DecoyPolicy::Force))
                .expect("fixture loads")
        else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_rows(), 5, "the five targets, without their decoys");
        assert_eq!(geom.n_stored_decoys(), 0);

        // The intensity sidecar has to shrink with the rows, or it stops being
        // parallel to the fragment arena.
        let TargetTable::Mzpaf { frag_intens, .. } =
            read_mzspeclib_library_file(&path, deciding_decoys(DecoyPolicy::Force))
                .expect("fixture loads")
        else {
            unreachable!()
        };
        assert_eq!(
            frag_intens.expect("mzSpecLib fills the sidecar").len(),
            geom.n_fragments()
        );
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

    /// Two predicted entries, alike but for the second claiming the `Decoy` set.
    /// A decoy declaration replaces the origin type without saying anything
    /// about how the masses were produced, so the decoy has to keep the peaks
    /// its target keeps.
    #[test]
    fn a_decoy_declaration_costs_an_entry_none_of_its_fragments() {
        let geom = arena("minimal_predicted_decoy.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 2);
        let (target, decoy) = (geom.rows().next().unwrap(), geom.rows().nth(1).unwrap());
        assert!(!geom.is_decoy(target) && geom.is_decoy(decoy));

        let labels = |row| {
            geom.frag_labels(row)
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
        };
        assert_eq!(labels(decoy).len(), 2, "both peaks, as the file declares");
        assert_eq!(labels(decoy), labels(target));
    }

    /// The `Spectrum=all` fallback is the library's claim, not an override of the
    /// entry's own: an entry naming `MS:1003073|observed spectrum` has said its
    /// m/z values are measured and its unresolvable peaks are dropped, while the
    /// decoy beside it -- which said nothing about mass provenance -- keeps its.
    #[test]
    fn an_entry_declaring_its_masses_observed_is_taken_at_its_word() {
        let geom = arena("entry_overrides_library_origin_type.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 3);
        let row = |name: &str| {
            geom.rows()
                .find(|r| geom.output_id(*r).to_string() == name)
                .expect("fixture names every entry")
        };
        assert_eq!(geom.frag_labels(row("PEPTIDEK/2")).len(), 2, "predicted");
        assert_eq!(geom.frag_labels(row("PDITPEEK/2")).len(), 2, "the decoy");
        assert_eq!(
            geom.frag_labels(row("SAMPLERK/2")).len(),
            0,
            "an observed spectrum offers no theoretical m/z to store"
        );
    }

    /// A whole entry lost is not the same shape of loss as a peak lost here and
    /// there, so it gets its own count. Both entries in this library declare two
    /// peaks and neither peak can be given a theoretical m/z.
    #[test]
    fn an_entry_whose_every_peak_was_dropped_is_counted() {
        let (table, degradation) = read_counting_degradation(
            &fixture("observed_without_annotation_masses.mzspeclib.txt"),
            deciding_decoys(DecoyPolicy::Never),
        )
        .expect("fixture loads");
        let TargetTable::Mzpaf { geom, .. } = table else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_rows(), 2);
        assert_eq!(geom.n_fragments(), 0, "nothing reached the arena");
        assert_eq!(degradation.rows_without_fragments, 2);
        assert_eq!(degradation.peaks_without_theoretical_mz, 4);

        // An entry that declares peaks and keeps them is not counted, and
        // neither is one that declares none.
        let (_, degradation) = read_counting_degradation(
            &fixture("minimal_predicted_decoy.mzspeclib.txt"),
            deciding_decoys(DecoyPolicy::Never),
        )
        .expect("fixture loads");
        assert_eq!(degradation.rows_without_fragments, 0);
    }

    /// mzSpecLib v1.0 §4.1.6 writes `MS:1003065|spectrum aggregation type`
    /// grouped with its replicate count, and mzannotate keeps a group it does
    /// not special-case out of the flattened params. Nothing else in this
    /// library says how the masses were produced -- there is no `MS:1003072`
    /// anywhere in it, so no library fallback either -- which leaves the grouped
    /// declaration as the only thing that can give a bare `b2` a theoretical
    /// m/z.
    #[test]
    fn a_grouped_aggregation_type_declares_an_entry_predicted() {
        let geom = arena("grouped_predicted_aggregation_type.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 2);
        let (target, decoy) = (geom.rows().next().unwrap(), geom.rows().nth(1).unwrap());
        assert!(!geom.is_decoy(target) && geom.is_decoy(decoy));
        assert_eq!(geom.frag_labels(target).len(), 2, "the target's two peaks");
        assert_eq!(geom.frag_labels(decoy).len(), 2, "and the decoy's");
    }

    /// An entry naming a decoy subtype and `MS:1003074|predicted spectrum` is
    /// asserting both, which mzSpecLib v1.0 §4.1.11 leaves room for: repeated
    /// instances whose values reconcile are not a case for last-wins. So it is a
    /// decoy, and its m/z values are calculated -- not whichever of the two the
    /// file happened to write first.
    #[test]
    fn an_entry_naming_two_origin_types_asserts_both_of_them() {
        let geom = arena("entry_declares_predicted_and_decoy.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 2);
        let (target, both) = (geom.rows().next().unwrap(), geom.rows().nth(1).unwrap());
        assert!(!geom.is_decoy(target));
        assert!(
            geom.is_decoy(both),
            "a second origin type does not withdraw the decoy claim"
        );
        assert_eq!(
            geom.frag_labels(both).len(),
            2,
            "and the decoy claim does not withdraw the predicted one"
        );
    }

    /// A shipped decoy that `--decoy-strategy force` drops is not an entry
    /// scoring against nothing, because it does not score. The count and the row
    /// total it is reported against have to measure the same entries, or the
    /// warning asserts a consequence for rows that are not in the arena and can
    /// exceed its own denominator.
    ///
    /// This also pins the library fallback on the library: the `Spectrum=all`
    /// set here declares its entries observed, so the decoy that overrode the
    /// origin type with `MS:1003195` has nothing to fall back to and its bare
    /// annotations reach no theoretical m/z.
    #[test]
    fn a_dropped_decoy_is_not_reported_as_scoring_against_nothing() {
        let path = fixture("shipped_decoy_without_readable_peaks.mzspeclib.txt");

        let (table, degradation) =
            read_counting_degradation(&path, deciding_decoys(DecoyPolicy::Never))
                .expect("fixture loads");
        let TargetTable::Mzpaf { geom, .. } = table else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_rows(), 2, "the target and its shipped decoy");
        assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 1);
        assert_eq!(
            geom.n_fragments(),
            0,
            "an observed library offers neither entry a theoretical m/z"
        );
        assert_eq!(degradation.rows_without_fragments, 2);

        let (table, degradation) =
            read_counting_degradation(&path, deciding_decoys(DecoyPolicy::Force))
                .expect("fixture loads");
        let TargetTable::Mzpaf { geom, .. } = table else {
            panic!("mzSpecLib carries ion chemistry");
        };
        assert_eq!(geom.n_rows(), 1, "the decoy never reached the arena");
        assert_eq!(
            degradation.rows_without_fragments, 1,
            "the target lost its peaks; the dropped decoy lost nothing"
        );
        assert!(degradation.rows_without_fragments <= geom.n_rows());
    }

    /// The format this project writes, read back by the reader that has to read
    /// it. It exercises every term msspeculator actually emits rather than the
    /// terms a vendor happens to use: `MS:1002815|inverse reduced ion mobility`
    /// for the timsTOF axis, `MS:1000896` for retention, and the
    /// project-defined `msspeculator:decoy_pair_id` for the target/decoy
    /// pairing.
    ///
    /// Stored as plain text, not gzipped, so it can be read in a review; the
    /// gzip path is covered by `diann_export.mzspeclib.txt.gz`. It was written
    /// by `timsseek build-library` and records how, in its own header: every
    /// `msspeculator:` attribute in the `<mzSpecLib>` block names one setting,
    /// including the FASTA's path and BLAKE2b digest, the model, and the
    /// digestion bounds. Regenerating it means reading those back off the
    /// fixture rather than guessing at flags.
    #[test]
    fn the_format_this_project_writes_reads_back() {
        let geom = arena("msspeculator_built.mzspeclib.txt");
        assert_eq!(geom.n_rows(), 4, "two targets and their two decoys");
        assert_eq!(geom.rows().filter(|r| geom.is_decoy(*r)).count(), 2);

        let groups: std::collections::HashSet<_> =
            geom.rows().map(|r| geom.decoy_group_code(r)).collect();
        assert_eq!(
            groups.len(),
            2,
            "`decoy_pair_id` puts each decoy in its target's group"
        );

        // Four peaks per entry, decoys included: the writer's `Decoy` set
        // replaces the origin type and leaves the aggregation type alone, and a
        // decoy with no fragments would be an FDR estimated against nothing.
        assert_eq!(geom.n_fragments(), 16);
        for row in geom.rows() {
            assert_eq!(
                geom.frag_labels(row).len(),
                4,
                "{} carries the four peaks it declares",
                geom.output_id(row),
            );
        }

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
        // The only fixture declaring `MS:1000896` in `minute`, and so the only
        // one that pins the minute path: the entry says 1.559414, and 60x that
        // is what a reader has to come back with. Which mzannotate build is in
        // the graph decides it -- the pinned one leaves a `minute`-declared
        // value alone for the `* 60` here to scale, a release divides it by
        // sixty first and hands back a sixtieth of the number.
        assert!(
            (geom.rt_seconds(first) - 93.564_84).abs() < 1e-3,
            "declared 1.559414 minute, got {} s",
            geom.rt_seconds(first),
        );
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
    /// The declared values, not a distinctness property. Several rows collapsing
    /// onto one value satisfies "more than one distinct value", and that is the
    /// shape a dropped retention time takes: six of these ten are grouped with
    /// `MS:1003174|attribute maximum` and `MS:1003175|attribute minimum`, which
    /// is where a reader is liable to lose them.
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
