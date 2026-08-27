//! Reader for the mzSpecLib text format (HUPO-PSI).
//!
//! # Why every field needs a fallback ladder
//!
//! mzSpecLib says *how* to write a controlled-vocabulary term, not *which*
//! term a writer must use. The two reference exports disagree on nearly
//! everything this reader needs:
//!
//! | field | DIA-NN writes | Spectronaut writes |
//! |---|---|---|
//! | precursor m/z | `MS:1000744` selected ion m/z | `MS:1003208` experimental precursor monoisotopic m/z |
//! | retention time | *nothing at all* | `MS:1000896` normalized retention time, in minutes |
//! | ion mobility | `MS:1002476` (written as `0.0`) | `MS:1002476` |
//!
//! So each field is resolved by trying terms in priority order, and RT carries
//! a unit that must be honoured rather than assumed.
//!
//! # Peak resolution
//!
//! The arena wants *theoretical* m/z. A peak list may carry either: whether it
//! does is declared by `MS:1003072|spectrum origin type`, which this reader
//! does not consult (see below). What it uses instead is the annotation's
//! mass-error suffix, `theoretical = observed - error` — correct for an
//! observed list, and a no-op on a theoretical one, where the error is `0.0`.
//!
//! So a theoretical mass exists only once a single identity is pinned:
//!
//! | annotation | kept | m/z |
//! |---|---|---|
//! | resolved, representable | real label | theoretical |
//! | resolved, not representable (`y1-HCOOH`) | unknown label | theoretical |
//! | resolved, no `/error` suffix | as above | observed, and counted |
//! | unannotated (`?`), tied ambiguity, malformed suffix | no | — |
//!
//! A peak with no single identity is skipped rather than stored at observed
//! m/z: an arena mixing observed and theoretical masses would be invisible
//! downstream. A known-but-unspellable identity is kept because it still has an
//! exact mass — only the label is lost. The suffix is optional in mzPAF, so the
//! third row is possible and is the one case where the mixture does happen;
//! `kept_at_observed_mz` is how it shows up.
//!
//! Ambiguous (comma-separated) annotations take the alternative with the
//! smallest absolute mass error. If that one is unrepresentable the peak gets an
//! unknown label rather than falling back to a worse-matching representable
//! alternative, which would assign both a wrong identity and a wrong mass.
//!
//! # Not implemented
//!
//! - **`MS:1003072|spectrum origin type`.** It distinguishes `MS:1003073`
//!   (observed), `MS:1003074` (predicted) and `MS:1003424` (theoretical m/z,
//!   observed intensity), i.e. exactly whether the subtraction above is needed.
//!   Both vendored fixtures write `/0.0` throughout, so the subtraction is a
//!   no-op on them either way.
//! - **Decoys.** `<AttributeSet Spectrum=DECOY>` + `MS:1003212` marks decoy
//!   spectra in SpectraST exports. Every row here is pushed as a target; see
//!   `ignored_attribute_set_entries`.

use crate::ion::{
    IonAnnot,
    UnknownIonCounter,
};
use crate::models::{
    LibCapabilities,
    QueryCollection,
};
use crate::serde::library_file::{
    FragmentSet,
    Inserted,
    LibraryArena,
    LibraryReadingError,
    finish_mzpaf_arena,
};
use micromzpaf::{
    MassError,
    split_mass_error,
};
use std::io::{
    BufRead,
    BufReader,
};
use std::path::Path;
use tracing::{
    info,
    warn,
};

/// First non-empty line of an mzSpecLib text file.
const MAGIC: &str = "<mzSpecLib>";

// CV term ladders, ordered most- to least-specific; the first present wins.

/// Precursor m/z. Experimental monoisotopic is preferred over `selected ion
/// m/z`, which on a quadrupole instrument is the isolation-window centre and
/// need not be the monoisotopic peak.
const PRECURSOR_MZ_TERMS: &[&str] = &[
    "MS:1003208", // experimental precursor monoisotopic m/z
    "MS:1003053", // theoretical monoisotopic m/z
    "MS:1000744", // selected ion m/z
];
/// Retention time. `normalized retention time` is an iRT-style index rather
/// than a clock reading, but it is what Spectronaut exports and the only RT
/// signal available in those files.
const RT_TERMS: &[&str] = &[
    "MS:1000894", // retention time
    "MS:1000896", // normalized retention time
];
/// Ion mobility. Note these are different quantities, not spellings of one:
/// `MS:1002815` is inverse reduced ion mobility (1/K0), which is what the
/// arena wants, while `MS:1002476` is a drift time. They are tried in that
/// order and a drift-time-only library is counted, since treating a drift time
/// as 1/K0 is only valid for instruments that report it that way.
const MOBILITY_INVERSE_REDUCED: &str = "MS:1002815";
const MOBILITY_DRIFT_TIME: &str = "MS:1002476";

const CHARGE_TERM: &str = "MS:1000041";
const STRIPPED_SEQ_TERM: &str = "MS:1000888";
const PROFORMA_TERM: &str = "MS:1003270";
const UNIT_TERM: &str = "UO:0000000";

const UNIT_MINUTE: &str = "UO:0000031";
const UNIT_SECOND: &str = "UO:0000010";

/// Declare the per-library tally so that the counters, the "is anything
/// wrong?" test and the log line all come from one list.
///
/// Spelled as a macro because the alternative — a struct plus a hand-written
/// `||` chain plus a hand-written `warn!` — is three places to update per
/// counter and the compiler checks none of them. That already went wrong once.
macro_rules! anomaly_counters {
    ($( $(#[$doc:meta])* $field:ident => $label:literal, )+) => {
        /// Per-library tally of everything that did not land verbatim in the
        /// arena.
        ///
        /// Reported once at the end of a load rather than per row: a consensus
        /// library can carry thousands of unannotated peaks, and a line each
        /// would bury the signal.
        #[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
        pub(crate) struct MzSpecLibStats {
            /// Peaks stored with their parsed annotation. The only counter
            /// here that is not an anomaly.
            pub kept_annotated: usize,
            $( $(#[$doc])* pub $field: usize, )+
        }

        impl MzSpecLibStats {
            /// Every anomaly counter paired with how to say it.
            fn anomalies(&self) -> impl Iterator<Item = (&'static str, usize)> {
                [ $( ($label, self.$field), )+ ].into_iter()
            }

            /// Fold `other` in. Used to hold a spectrum's counts aside until
            /// it is known to be kept, so a dropped spectrum does not also
            /// report the peaks it would have contributed.
            fn merge(&mut self, other: &Self) {
                self.kept_annotated += other.kept_annotated;
                $( self.$field += other.$field; )+
            }
        }
    };
}

anomaly_counters! {
    /// Identity known but not representable (a loss outside the table, a
    /// modified immonium). Stored with an unknown label and an exact mass.
    kept_unknown_label => "kept with an unknown label",
    /// `?` — no annotation, so no mass error, so no theoretical m/z.
    skipped_unannotated => "skipped as unannotated",
    /// Comma-separated alternatives that tied on absolute mass error.
    skipped_ambiguous => "skipped as ambiguous",
    /// An `/error` suffix that would not parse. The peak has no recoverable
    /// mass, so it cannot be stored at any label.
    skipped_malformed_mass_error => "skipped for a malformed mass error",
    /// Peaks kept at their OBSERVED m/z because the annotation carried no
    /// `/error` suffix (which mzPAF makes optional). Everything else in the
    /// arena is theoretical, so this is the one counter that means the arena
    /// mixes the two.
    kept_at_observed_mz => "kept at observed m/z (no mass-error suffix)",
    /// Peaks dropped because their label collided with one already in the
    /// precursor.
    dropped_duplicate_label => "dropped for a duplicate label",
    /// Peaks dropped because the precursor had already spent all 255 unknown
    /// labels, so no distinct one was left.
    dropped_unknown_over_capacity => "dropped with the unknown labels exhausted",
    /// Spectra with no `MS:1000888` stripped sequence, whose bare residues
    /// were derived from the proforma instead. Not required by mzSpecLib.
    stripped_sequence_derived => "spectra with a derived stripped sequence",
    /// Spectra with no retention-time term at all.
    spectra_without_rt => "spectra without an RT",
    /// Spectra whose mobility came from a drift time rather than 1/K0.
    spectra_with_drift_time_mobility => "spectra using a drift time as mobility",
    /// Spectra whose retention time carried a unit this reader does not know.
    spectra_with_unknown_rt_unit => "spectra with an unknown RT unit",
    /// Spectra dropped for missing or unparseable precursor m/z, charge or
    /// sequence.
    dropped_malformed_spectrum => "spectra dropped as malformed",
    /// Precursors dropped for having no usable peak left.
    dropped_empty_precursors => "precursors dropped as empty",
    /// Attributes inside an `<AttributeSet ...>` block, which this reader does
    /// not apply. Both vendored fixtures have these blocks empty; a writer
    /// that hoists an RT unit group or a DECOY marker into one would otherwise
    /// lose it with no trace.
    ignored_attribute_set_entries => "attribute-set entries ignored",
}

/// One spectrum converted into the shape [`QueryCollection::push_row`] takes.
struct ArenaRow {
    precursor_mz: f64,
    charge: u8,
    rt_seconds: f32,
    mobility: f32,
    frags: FragmentSet,
    stripped: String,
    modified: String,
}

impl MzSpecLibStats {
    fn report(&self, path: &Path) {
        let flagged: Vec<String> = self
            .anomalies()
            .filter(|(_, n)| *n > 0)
            .map(|(label, n)| format!("{n} {label}"))
            .collect();
        if flagged.is_empty() {
            info!(
                "mzSpecLib {}: {} peaks, all annotated and representable",
                path.display(),
                self.kept_annotated
            );
        } else {
            warn!(
                "mzSpecLib {}: kept {} annotated peaks; {}",
                path.display(),
                self.kept_annotated,
                flagged.join(", "),
            );
        }
    }
}

/// Which `<...>` block an attribute was written in.
///
/// `[n]` group ids are scoped to their block, so two blocks can both use `[2]`
/// for unrelated groups — `spectronaut.mzSpecLib.txt` already does. Carrying
/// the block makes a group id unique within a spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BlockId {
    Spectrum,
    Analyte(u32),
    Interpretation(u32),
}

impl BlockId {
    /// Parse a `<Analyte=1>`-style header. `None` for a header that opens no
    /// attribute block (`<Peaks>`, `<AttributeSet ...>`).
    fn parse(header: &str) -> Option<Self> {
        let inner = header.strip_prefix('<')?.strip_suffix('>')?;
        if inner == "Spectrum" || inner.starts_with("Spectrum=") {
            return Some(Self::Spectrum);
        }
        let (kind, id) = inner.split_once('=')?;
        let id = id.parse().ok()?;
        match kind {
            "Analyte" => Some(Self::Analyte(id)),
            "Interpretation" => Some(Self::Interpretation(id)),
            _ => None,
        }
    }
}

/// One `ACC|name=value` attribute, with its `[n]` group tag scoped to the
/// block it was written in.
#[derive(Debug, Clone)]
struct Attr {
    group: Option<(BlockId, u32)>,
    accession: String,
    value: String,
}

impl Attr {
    fn parse(line: &str, block: BlockId) -> Option<Attr> {
        let (group, rest) = match line.strip_prefix('[') {
            Some(r) => {
                let (g, r) = r.split_once(']')?;
                (Some((block, g.parse().ok()?)), r)
            }
            None => (None, line),
        };
        let (key, value) = rest.split_once('=')?;
        let accession = key.split('|').next()?.to_string();
        Some(Attr {
            group,
            accession,
            value: value.to_string(),
        })
    }

    /// The accession out of a `ACC|name` *value* (as opposed to a key), for
    /// terms whose value is itself a CV term — e.g. `unit=UO:0000031|minute`.
    fn value_accession(&self) -> &str {
        self.value.split('|').next().unwrap_or(&self.value)
    }
}

/// Attributes collected for one spectrum (its own plus its analyte's).
#[derive(Debug, Default)]
struct AttrBag(Vec<Attr>);

impl AttrBag {
    /// Spectrum-block attributes win over Analyte/Interpretation ones. Without
    /// this the answer would depend on which block the writer emitted first.
    fn find(&self, accession: &str) -> Option<&Attr> {
        let matching = || self.0.iter().filter(|a| a.accession == accession);
        matching()
            .find(|a| a.group.is_none_or(|(b, _)| b == BlockId::Spectrum))
            .or_else(|| matching().next())
    }

    fn first_of(&self, accessions: &[&str]) -> Option<&Attr> {
        accessions.iter().find_map(|a| self.find(a))
    }

    fn f64_of(&self, accessions: &[&str]) -> Option<f64> {
        self.first_of(accessions)?.value.parse().ok()
    }

    /// The unit term attached to `attr` via its `[n]` group, if any.
    ///
    /// The group is matched on `(block, id)`, not `id` alone: an unrelated
    /// `[2]` in the Analyte block must not supply the unit for a `[2]` in the
    /// Spectrum block. Getting that wrong is a silent 60x RT error.
    fn unit_for(&self, attr: &Attr) -> Option<&str> {
        let group = attr.group?;
        self.0
            .iter()
            .find(|a| a.group == Some(group) && a.accession == UNIT_TERM)
            .map(|a| a.value_accession())
    }
}

/// A spectrum accumulated from the text stream, before conversion.
#[derive(Debug, Default)]
struct RawSpectrum {
    attrs: AttrBag,
    /// `(observed mz, intensity, annotation)`
    peaks: Vec<(f64, f32, String)>,
}

/// Why a peak cannot be stored at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SkipReason {
    /// No annotation, so no mass error, so no theoretical m/z.
    Unannotated,
    /// Alternatives that pin no single identity, so likewise no theoretical
    /// m/z. Storing the observed one would mix the two.
    Ambiguous,
    /// The `/error` suffix was present but unparseable. Distinct from
    /// `Ambiguous`: one malformed annotation is malformed, not ambiguous.
    MalformedMassError,
}

/// What resolving one peak's annotation produced.
///
/// The mass error rides along on the two variants that keep the peak, because
/// it is what recovers theoretical m/z from the observed value in the file. It
/// is absent from `Skip` because a skipped peak has no mass to recover.
enum Resolved {
    /// Parsed cleanly; store with this label.
    Annotated(IonAnnot, Option<MassError>),
    /// Identity known, not representable. Store with an unknown label and the
    /// exact mass.
    UnknownLabel(Option<MassError>),
    Skip(SkipReason),
}

/// Resolve one annotation string into a storage decision.
fn resolve_annotation(annotation: &str) -> Resolved {
    let annotation = annotation.trim();
    if annotation.is_empty() || annotation == "?" {
        return Resolved::Skip(SkipReason::Unannotated);
    }

    // Splitting the error off comes first: it works even when the ion will not
    // parse, which is exactly the case that still needs an exact mass. A
    // malformed suffix on ANY alternative makes the whole peak unresolvable —
    // dropping just that one would silently turn an ambiguous peak into an
    // unambiguous one.
    let mut alternatives = Vec::new();
    for alt in annotation.split(',') {
        let Ok(parsed) = split_mass_error(alt.trim()) else {
            return Resolved::Skip(SkipReason::MalformedMassError);
        };
        alternatives.push(parsed);
    }

    let (ion_str, mass_error) = if let [single] = alternatives[..] {
        single
    } else {
        // Closest by absolute mass error. Comparing a Da magnitude against a
        // ppm one would be meaningless, but a library uses one unit
        // throughout. Errors are parsed decimal literals, so equal ones
        // compare exactly and a tie pins no identity.
        let magnitude = |m: &Option<MassError>| match m {
            Some(MassError::Da(v) | MassError::Ppm(v)) => v.abs(),
            None => f64::INFINITY,
        };
        let best = alternatives
            .iter()
            .map(|(_, e)| magnitude(e))
            .fold(f64::INFINITY, f64::min);
        if !best.is_finite() {
            return Resolved::Skip(SkipReason::Ambiguous);
        }
        let mut winners = alternatives.iter().filter(|(_, e)| magnitude(e) == best);
        let winner = *winners.next().expect("the minimum came from this iterator");
        if winners.next().is_some() {
            return Resolved::Skip(SkipReason::Ambiguous);
        }
        winner
    };
    match IonAnnot::try_from(ion_str) {
        Ok(ion) => Resolved::Annotated(ion, mass_error),
        // Keep the peak and its exact mass, lose only the label.
        Err(_) => Resolved::UnknownLabel(mass_error),
    }
}

/// Convert one accumulated spectrum into an arena row, counting whichever way
/// it failed.
///
/// The counting lives here rather than inside [`spectrum_row`] so that every
/// `?` in there lands on a tally. Without it a library that dropped half its
/// spectra for a missing charge still reports "all annotated and
/// representable", which is the one thing this module exists to prevent.
fn convert_spectrum(raw: &RawSpectrum, stats: &mut MzSpecLibStats) -> Option<ArenaRow> {
    // A malformed spectrum bails out of `spectrum_row` via `?`, possibly after
    // incrementing RT/mobility counters on the way. Those are held aside and
    // discarded on the bail-out, so it is reported once as malformed rather
    // than also as "without an RT" and "using a drift time".
    //
    // An *empty* precursor is structurally fine, so its per-peak counts are
    // kept: they are the explanation for why it came out empty.
    let mut local = MzSpecLibStats::default();
    let Some(row) = spectrum_row(raw, &mut local) else {
        stats.dropped_malformed_spectrum += 1;
        return None;
    };
    stats.merge(&local);

    if row.frags.is_empty() {
        stats.dropped_empty_precursors += 1;
        return None;
    }
    Some(row)
}

/// The conversion proper: `None` means the spectrum lacked something
/// structural (precursor m/z, charge, a parseable RT or mobility, a sequence).
fn spectrum_row(raw: &RawSpectrum, stats: &mut MzSpecLibStats) -> Option<ArenaRow> {
    let precursor_mz = raw.attrs.f64_of(PRECURSOR_MZ_TERMS)?;
    let charge: u8 = raw.attrs.find(CHARGE_TERM)?.value.parse().ok()?;

    let rt_seconds = match raw.attrs.first_of(RT_TERMS) {
        Some(attr) => {
            let v: f64 = attr.value.parse().ok()?;
            // Honour the unit rather than assuming: Spectronaut writes minutes.
            match raw.attrs.unit_for(attr) {
                Some(UNIT_SECOND) => v,
                Some(UNIT_MINUTE) | None => v * 60.0,
                // An unrecognized unit is counted, not silently treated as
                // minutes: guessing wrong here is a 60x error in the RT the
                // whole extraction window is built around.
                Some(_) => {
                    stats.spectra_with_unknown_rt_unit += 1;
                    v * 60.0
                }
            }
        }
        None => {
            stats.spectra_without_rt += 1;
            0.0
        }
    };

    let mobility = match raw.attrs.find(MOBILITY_INVERSE_REDUCED) {
        Some(a) => a.value.parse().ok()?,
        None => match raw.attrs.find(MOBILITY_DRIFT_TIME) {
            Some(a) => {
                let drift: f32 = a.value.parse().ok()?;
                // DIA-NN writes `MS:1002476|ion mobility drift time=0.0` on
                // every spectrum, meaning "unset". Counting that as a real
                // drift time warns on every DIA-NN load and buries the case
                // this counter exists for.
                if drift != 0.0 {
                    stats.spectra_with_drift_time_mobility += 1;
                }
                drift
            }
            // Absent is fine — an unset mobility is 0.0, same as DIA-NN writes.
            // A *present but malformed* one drops the spectrum instead.
            None => 0.0,
        },
    };

    // The proforma term carries a trailing `/charge` that is not part of the
    // peptidoform.
    let modified = raw.attrs.find(PROFORMA_TERM).map(|a| {
        a.value
            .rsplit_once('/')
            .map(|(p, _)| p.to_string())
            .unwrap_or_else(|| a.value.clone())
    });

    // `MS:1000888` is not required by mzSpecLib. Deriving the stripped form
    // from the proforma is better than defaulting to "", which reads as a
    // zero-residue peptide and silently routes the row to averagine isotopes
    // — tallied downstream as `n_averagine_fallback`, at a layer that can no
    // longer say the sequence was simply absent.
    let stripped = match raw.attrs.find(STRIPPED_SEQ_TERM) {
        Some(a) => a.value.clone(),
        None => {
            stats.stripped_sequence_derived += 1;
            crate::utils::sequence::strip_mods(modified.as_deref()?)
        }
    };
    let modified = modified.unwrap_or_else(|| stripped.clone());
    if stripped.is_empty() && modified.is_empty() {
        return None;
    }

    let mut frags = FragmentSet::with_capacity(raw.peaks.len());
    let mut unknown_ions = UnknownIonCounter::new();

    for (observed_mz, intensity, annotation) in &raw.peaks {
        let (label, mass_error, annotated) = match resolve_annotation(annotation) {
            Resolved::Annotated(ion, error) => (ion, error, true),
            Resolved::UnknownLabel(error) => match unknown_ions.next(1) {
                Ok(ion) => (ion, error, false),
                Err(_) => {
                    stats.dropped_unknown_over_capacity += 1;
                    continue;
                }
            },
            Resolved::Skip(SkipReason::Unannotated) => {
                stats.skipped_unannotated += 1;
                continue;
            }
            Resolved::Skip(SkipReason::Ambiguous) => {
                stats.skipped_ambiguous += 1;
                continue;
            }
            Resolved::Skip(SkipReason::MalformedMassError) => {
                stats.skipped_malformed_mass_error += 1;
                continue;
            }
        };

        // The `/error` suffix is optional in mzPAF. Without it the observed
        // m/z is the best available value, but it is NOT the theoretical one
        // the rest of the arena holds.
        let mz = match mass_error {
            Some(e) => e.theoretical_from_observed(*observed_mz),
            None => *observed_mz,
        };

        // `FragmentSet` owns the per-precursor uniqueness invariant; a
        // collision collapses onto the more intense peak. Unknown labels come
        // off a monotonic counter and cannot collide, so this only ever fires
        // for annotated ones. Every kept counter is behind this check, so a
        // collided peak is tallied once, as dropped.
        if frags.insert(label, mz, *intensity) == Inserted::Collapsed {
            stats.dropped_duplicate_label += 1;
            continue;
        }
        if mass_error.is_none() {
            stats.kept_at_observed_mz += 1;
        }
        if annotated {
            stats.kept_annotated += 1;
        } else {
            stats.kept_unknown_label += 1;
        }
    }

    Some(ArenaRow {
        precursor_mz,
        charge,
        rt_seconds: rt_seconds as f32,
        mobility,
        frags,
        stripped,
        modified,
    })
}

/// Cheap probe: the format's magic first line.
pub fn sniff_mzspeclib_library_file<T: AsRef<Path>>(path: T) -> bool {
    let Ok(file) = std::fs::File::open(path.as_ref()) else {
        return false;
    };
    // Lazy, so only the first non-empty line is read: this stays O(1) on a
    // multi-gigabyte library.
    BufReader::new(file)
        .lines()
        .map_while(Result::ok)
        .find(|l| !l.trim().is_empty())
        .is_some_and(|l| l.trim() == MAGIC)
}

/// Which part of a spectrum block the line loop is inside.
enum Section {
    Attributes,
    Peaks,
}

/// Convert an accumulated spectrum and append it to the arena.
fn flush(
    current: Option<RawSpectrum>,
    geom: &mut QueryCollection<IonAnnot>,
    frag_intens: &mut Vec<f32>,
    stats: &mut MzSpecLibStats,
) {
    let Some(raw) = current else { return };
    let Some(row) = convert_spectrum(&raw, stats) else {
        return;
    };
    row.frags.extend_sidecar(frag_intens);
    geom.push_row(
        row.precursor_mz,
        row.charge,
        row.rt_seconds,
        row.mobility,
        row.frags.frags(),
        &row.stripped,
        &row.modified,
        &[],
        false,
    );
}

/// Read an mzSpecLib text file into the columnar arena.
pub fn read_mzspeclib_library_file<T: AsRef<Path>>(
    path: T,
) -> Result<LibraryArena, LibraryReadingError> {
    let path = path.as_ref();
    let file = std::fs::File::open(path).map_err(LibraryReadingError::IoError)?;
    let reader = BufReader::new(file);

    let mut geom = QueryCollection::with_capabilities(LibCapabilities::default_diann_no_decoys());
    let mut frag_intens: Vec<f32> = Vec::new();
    let mut stats = MzSpecLibStats::default();

    let mut current: Option<RawSpectrum> = None;
    let mut section = Section::Attributes;
    // Which block the attributes now being read belong to. `[n]` group ids are
    // scoped to it, so it has to be threaded into every `Attr::parse`.
    let mut block = BlockId::Spectrum;
    let mut in_attribute_set = false;

    for line in reader.lines() {
        let line = line.map_err(LibraryReadingError::IoError)?;
        let trimmed = line.trim_end();

        // Any `<...>` header ends the peak list; only `<Peaks>` opens one.
        // `<Analyte=n>` and `<Interpretation=n>` attributes are folded into
        // the spectrum's bag: this reader wants the union, not the hierarchy,
        // but each keeps its own group scope.
        if trimmed.starts_with('<') {
            if trimmed.starts_with("<Spectrum=") {
                flush(current.take(), &mut geom, &mut frag_intens, &mut stats);
                current = Some(RawSpectrum::default());
            }
            section = if trimmed == "<Peaks>" {
                Section::Peaks
            } else {
                Section::Attributes
            };
            // `<AttributeSet ...>` declares defaults applied to whole classes
            // of spectra by name. Honouring them is not implemented, so a
            // non-empty one is counted rather than silently dropped: it can
            // carry the RT unit, or the DECOY marker.
            in_attribute_set = trimmed.starts_with("<AttributeSet");
            block = BlockId::parse(trimmed).unwrap_or(block);
            continue;
        }
        if trimmed.is_empty() {
            section = Section::Attributes;
            continue;
        }

        let Some(spec) = current.as_mut() else {
            if in_attribute_set && Attr::parse(trimmed, BlockId::Spectrum).is_some() {
                stats.ignored_attribute_set_entries += 1;
            }
            continue; // library-level header
        };

        match section {
            Section::Peaks => {
                let mut cols = trimmed.split('\t');
                let (Some(mz), Some(intensity)) = (cols.next(), cols.next()) else {
                    continue;
                };
                let (Ok(mz), Ok(intensity)) =
                    (mz.trim().parse::<f64>(), intensity.trim().parse::<f32>())
                else {
                    continue;
                };
                let annotation = cols.next().unwrap_or("?").to_string();
                spec.peaks.push((mz, intensity, annotation));
            }
            Section::Attributes => {
                if let Some(attr) = Attr::parse(trimmed, block) {
                    spec.attrs.0.push(attr);
                }
            }
        }
    }
    flush(current.take(), &mut geom, &mut frag_intens, &mut stats);

    stats.report(path);

    if geom.n_rows() == 0 {
        return Err(LibraryReadingError::SpeclibParse(format!(
            "mzSpecLib {} yielded no usable precursors",
            path.display()
        )));
    }
    finish_mzpaf_arena(geom, frag_intens)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/mzspeclib_io_files")
            .join(name)
    }

    #[test]
    fn sniffs_only_mzspeclib() {
        assert!(sniff_mzspeclib_library_file(fixture("diann.mzSpecLib.txt")));
        assert!(sniff_mzspeclib_library_file(fixture(
            "spectronaut.mzSpecLib.txt"
        )));
        // A DIA-NN TSV must not be claimed.
        let tsv = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/diann_io_files/sample_lib.tsv");
        assert!(!sniff_mzspeclib_library_file(tsv));
    }

    #[test]
    fn reads_diann_export() {
        let arena = read_mzspeclib_library_file(fixture("diann.mzSpecLib.txt")).unwrap();
        let LibraryArena::Mzpaf { geom, frag_intens } = arena else {
            panic!("mzSpecLib must build an mzpaf arena");
        };
        assert!(geom.n_rows() > 0);
        let intens = frag_intens.expect("reference intensities are populated");
        assert_eq!(intens.len(), geom.frag_labels.len());
    }

    /// DIA-NN's export writes every mass error as exactly `0.0`, so observed
    /// and theoretical coincide and the m/z must pass through untouched.
    #[test]
    fn zero_mass_error_leaves_mz_unchanged() {
        let arena = read_mzspeclib_library_file(fixture("diann.mzSpecLib.txt")).unwrap();
        let LibraryArena::Mzpaf { geom, .. } = arena else {
            unreachable!()
        };
        // First spectrum's first peak in the fixture: 427.22995, annotated b6/0.0
        let first = geom.frag_mzs[0];
        assert!(
            (first - 427.22995).abs() < 1e-6,
            "expected the observed m/z verbatim, got {first}"
        );
    }

    /// Spectronaut's export carries `-H2O`/`-NH3` losses. `IonAnnot`
    /// represents those, so they must keep real labels rather than degrading
    /// to unknown.
    #[test]
    fn spectronaut_losses_keep_real_labels() {
        let arena = read_mzspeclib_library_file(fixture("spectronaut.mzSpecLib.txt")).unwrap();
        let LibraryArena::Mzpaf { geom, .. } = arena else {
            unreachable!()
        };
        let has_loss = geom
            .frag_labels
            .iter()
            .any(|l| l.loss() != micromzpaf::NeutralLoss::None);
        assert!(has_loss, "expected at least one loss-bearing label");
        let unknowns = geom
            .frag_labels
            .iter()
            .filter(|l| l.try_get_ordinal().is_none() && l.loss() == micromzpaf::NeutralLoss::None)
            .count();
        assert_eq!(unknowns, 0, "no peak should need an unknown label here");
    }

    /// RT is unit-tagged; Spectronaut writes minutes and the arena wants
    /// seconds. Getting this wrong is a silent 60x error.
    #[test]
    fn retention_time_honours_its_unit() {
        let arena = read_mzspeclib_library_file(fixture("spectronaut.mzSpecLib.txt")).unwrap();
        let LibraryArena::Mzpaf { geom, .. } = arena else {
            unreachable!()
        };
        // Fixture's first spectrum: normalized retention time = 28.658491 min.
        let rt = geom.rt_seconds[0];
        assert!(
            (rt - 28.658491 * 60.0).abs() < 0.01,
            "expected minutes converted to seconds, got {rt}"
        );
    }

    /// `[n]` group ids are scoped to their block, so an Analyte `[2]` must not
    /// supply the unit for a Spectrum `[2]`.
    ///
    /// `spectronaut.mzSpecLib.txt` already carries both: `[2]` is the RT +
    /// unit pair in its Spectrum block and the NCBI TaxID in its Analyte
    /// block. It survives today only because the Analyte's `[2]` happens to
    /// carry no unit term. If it did, the RT would be read in the wrong unit —
    /// a silent 60x error in the value the extraction window is built on.
    #[test]
    fn group_ids_do_not_leak_between_blocks() {
        let grouped = |block: BlockId, id: u32, accession: &str, value: &str| Attr {
            group: Some((block, id)),
            accession: accession.to_string(),
            value: value.to_string(),
        };
        let rt = grouped(BlockId::Spectrum, 2, RT_TERMS[1], "10.0");
        let bag = AttrBag(vec![
            rt.clone(),
            // Same group id, different block, and it does carry a unit.
            grouped(BlockId::Analyte(1), 2, "MS:1001467", "9606"),
            grouped(BlockId::Analyte(1), 2, UNIT_TERM, UNIT_SECOND),
        ]);

        assert_eq!(
            bag.unit_for(&rt),
            None,
            "the Analyte block's unit must not reach the Spectrum block's group"
        );

        // And when the unit is in the same block, it is found.
        let mut same_block = bag;
        same_block
            .0
            .push(grouped(BlockId::Spectrum, 2, UNIT_TERM, UNIT_SECOND));
        assert_eq!(same_block.unit_for(&rt), Some(UNIT_SECOND));
    }

    /// The `/error` suffix is optional in mzPAF, so a bare `y5` is legal and
    /// its m/z stays observed while the rest of the arena is theoretical.
    /// Both vendored fixtures write `/0.0` everywhere, so nothing else covers
    /// this. A malformed suffix is a different verdict from an ambiguous one.
    #[test]
    fn a_missing_mass_error_is_counted_and_a_malformed_one_is_skipped() {
        assert!(matches!(
            resolve_annotation("y5"),
            Resolved::Annotated(_, None)
        ));
        assert!(matches!(
            resolve_annotation("y5/not-a-number"),
            Resolved::Skip(SkipReason::MalformedMassError)
        ));

        let attr = |accession: &str, value: &str| Attr {
            group: None,
            accession: accession.to_string(),
            value: value.to_string(),
        };
        let raw = RawSpectrum {
            attrs: AttrBag(vec![
                attr(PRECURSOR_MZ_TERMS[0], "500.25"),
                attr(CHARGE_TERM, "2"),
                attr(STRIPPED_SEQ_TERM, "PEPTIDEK"),
                attr(RT_TERMS[0], "10.0"),
            ]),
            peaks: vec![
                (175.1, 1.0, "y1/0.0".to_string()),
                (288.2, 1.0, "y2".to_string()),
            ],
        };

        let mut stats = MzSpecLibStats::default();
        let row = convert_spectrum(&raw, &mut stats).expect("both peaks are usable");
        assert_eq!(row.frags.len(), 2);
        assert_eq!(stats.kept_annotated, 2);
        assert_eq!(
            stats.kept_at_observed_mz, 1,
            "only the suffix-less peak keeps its observed m/z"
        );
    }

    /// A collided peak must be reported once, as dropped — not as kept AND
    /// dropped, which is what counting before the uniqueness check produced.
    #[test]
    fn a_duplicate_label_is_counted_once_as_dropped() {
        let attr = |accession: &str, value: &str| Attr {
            group: None,
            accession: accession.to_string(),
            value: value.to_string(),
        };
        let raw = RawSpectrum {
            attrs: AttrBag(vec![
                attr(PRECURSOR_MZ_TERMS[0], "500.25"),
                attr(CHARGE_TERM, "2"),
                attr(STRIPPED_SEQ_TERM, "PEPTIDEK"),
                attr(RT_TERMS[0], "10.0"),
            ]),
            peaks: vec![
                (175.1, 1.0, "y1/0.0".to_string()),
                (175.2, 2.0, "y1/0.0".to_string()),
            ],
        };

        let mut stats = MzSpecLibStats::default();
        let row = convert_spectrum(&raw, &mut stats).expect("the first peak is usable");
        assert_eq!(row.frags.len(), 1);
        assert_eq!(stats.kept_annotated, 1);
        assert_eq!(stats.dropped_duplicate_label, 1);
    }

    /// DIA-NN writes `MS:1002476|ion mobility drift time=0.0` on every
    /// spectrum, meaning "unset". Counting that warns on every DIA-NN load.
    #[test]
    fn a_zero_drift_time_is_absent_not_a_drift_time() {
        let mut stats = MzSpecLibStats::default();
        let attr = |accession: &str, value: &str| Attr {
            group: None,
            accession: accession.to_string(),
            value: value.to_string(),
        };
        let with_drift = |drift: &str| RawSpectrum {
            attrs: AttrBag(vec![
                attr(PRECURSOR_MZ_TERMS[0], "500.25"),
                attr(CHARGE_TERM, "2"),
                attr(STRIPPED_SEQ_TERM, "PEPTIDEK"),
                attr(RT_TERMS[0], "10.0"),
                attr(MOBILITY_DRIFT_TIME, drift),
            ]),
            peaks: vec![(175.1, 1.0, "y1/0.0".to_string())],
        };

        assert!(convert_spectrum(&with_drift("0.0"), &mut stats).is_some());
        assert_eq!(stats.spectra_with_drift_time_mobility, 0);
        assert!(convert_spectrum(&with_drift("0.85"), &mut stats).is_some());
        assert_eq!(stats.spectra_with_drift_time_mobility, 1);
    }

    /// A spectrum missing a structural field is dropped by `?` deep inside
    /// `spectrum_row`. Nothing there increments a counter, so this is what
    /// stops the load from reporting "all annotated and representable" while
    /// having silently lost half the library.
    #[test]
    fn structurally_broken_spectra_are_counted_not_swallowed() {
        let attr = |accession: &str, value: &str| Attr {
            group: None,
            accession: accession.to_string(),
            value: value.to_string(),
        };
        let good = |charge: &str| RawSpectrum {
            attrs: AttrBag(vec![
                attr(PRECURSOR_MZ_TERMS[0], "500.25"),
                attr(CHARGE_TERM, charge),
                attr(STRIPPED_SEQ_TERM, "PEPTIDEK"),
                attr(RT_TERMS[0], "10.0"),
            ]),
            peaks: vec![(175.1, 1.0, "y1/0.0".to_string())],
        };

        let mut stats = MzSpecLibStats::default();
        assert!(convert_spectrum(&good("2"), &mut stats).is_some());
        // Unparseable charge — the `?` that used to vanish.
        assert!(convert_spectrum(&good("not-a-number"), &mut stats).is_none());
        // No peak survives resolution, which is a different failure.
        let mut empty = good("2");
        empty.peaks = vec![(175.1, 1.0, "?".to_string())];
        assert!(convert_spectrum(&empty, &mut stats).is_none());

        assert_eq!(stats.kept_annotated, 1);
        assert_eq!(stats.dropped_malformed_spectrum, 1);
        assert_eq!(stats.dropped_empty_precursors, 1);
        assert_eq!(stats.skipped_unannotated, 1);
        assert!(
            stats.anomalies().any(|(_, n)| n > 0),
            "the report must not claim a clean load"
        );
    }

    #[test]
    fn resolves_unambiguous_representable() {
        assert!(matches!(
            resolve_annotation("y5/-0.0005"),
            Resolved::Annotated(_, Some(MassError::Da(-0.0005)))
        ));
    }

    /// Known identity, unrepresentable spelling: keep the peak and the exact
    /// mass, erase only the label. The mass error must survive so theoretical
    /// m/z stays exact.
    #[test]
    fn unrepresentable_loss_keeps_peak_with_unknown_label() {
        assert!(matches!(
            resolve_annotation("y1-HCOOH/0.0003"),
            Resolved::UnknownLabel(Some(MassError::Da(0.0003)))
        ));
    }

    #[test]
    fn unannotated_peak_is_skipped() {
        assert!(matches!(
            resolve_annotation("?"),
            Resolved::Skip(SkipReason::Unannotated)
        ));
    }

    /// Closest-by-error wins; if that alternative is unrepresentable the peak
    /// takes an unknown label rather than falling back to the representable
    /// one, which would assign a wrong identity and a wrong mass.
    #[test]
    fn ambiguity_resolves_to_the_closest_not_the_representable() {
        // a2 is representable and further; y2-CO2-NH3 is closer and is not.
        assert!(
            matches!(
                resolve_annotation("a2/-0.0040,y2-CO2-NH3/-0.0001"),
                Resolved::UnknownLabel(Some(MassError::Da(-0.0001)))
            ),
            "the closest alternative wins even when unrepresentable"
        );

        // When the closest one IS representable, it is used.
        assert!(matches!(
            resolve_annotation("a2/-0.0001,y2-CO2-NH3/-0.0040"),
            Resolved::Annotated(..)
        ));
    }

    /// An exact tie pins no identity, so no theoretical m/z exists and the peak
    /// cannot be stored without mixing observed and theoretical masses.
    #[test]
    fn tied_ambiguity_is_skipped() {
        assert!(matches!(
            resolve_annotation("a2/-0.0004,y2-CO2-NH3/-0.0004"),
            Resolved::Skip(SkipReason::Ambiguous)
        ));
    }

    #[test]
    fn mass_error_recovers_theoretical() {
        // Real SpectraST peak: y1 for C-terminal R, observed 175.1184.
        let e = MassError::Da(-0.0005);
        assert!((e.theoretical_from_observed(175.1184) - 175.1189).abs() < 1e-9);
    }

    /// The public entry point must dispatch mzSpecLib itself. It is the
    /// registry's first entry precisely because the JSON reader at the end
    /// accepts anything, so a regression in the sniff falls through to it and
    /// surfaces as a generic parse error rather than as "wrong reader".
    #[test]
    fn public_read_library_file_dispatches_mzspeclib() {
        use crate::serde::read_library_file;
        for name in ["diann.mzSpecLib.txt", "spectronaut.mzSpecLib.txt"] {
            let arena = read_library_file(fixture(name))
                .unwrap_or_else(|e| panic!("{name} must load through the registry: {e:?}"));
            let LibraryArena::Mzpaf { geom, frag_intens } = arena else {
                panic!("{name} must land in the mzpaf arena");
            };
            assert!(geom.n_rows() > 0, "{name} produced no precursors");
            assert!(
                frag_intens.is_some(),
                "{name} must populate the reference-intensity sidecar"
            );
        }
    }
}
