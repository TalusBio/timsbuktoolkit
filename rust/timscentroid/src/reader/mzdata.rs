//! mzML / mzdata ingest → `IndexedTimstofPeaks`.
//!
//! mzML has no ion-mobility axis, so the built index carries a non-`Ook0`
//! [`MobilityKind`] (`Absent`, or `Unsupported` for FAIMS/drift) and the query
//! engine takes the mobility-less path. See `.plans/mzml-support-design.md` §3.
//!
//! Single sequential pass. Centroided input only (profile is error-logged and
//! the file skipped). MS1↔MS2 are reconciled on a shared integer cycle: each
//! MS1 and the MS2 window-sweep that follows it (until the next MS1) share one
//! cycle index, whose RT is the MS1 scan start time.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;

use half::f16;
use http::Uri;
use mzdata::io::mzml::MzMLReader;
use mzdata::prelude::*;
use tracing::warn;

use super::{
    Manifest,
    RawReader,
    ReadError,
    ResolvedSource,
    Sniff,
    path_ends_with,
};
use crate::centroiding::CentroidingConfig;
use crate::dimension::MobilityKind;
use crate::geometry::QuadrupoleIsolationScheme;
use crate::indexing::{
    IndexedPeak,
    IndexedPeakGroup,
    IndexedTimstofPeaks,
};
use crate::rt_mapping::{
    CycleToRTMapping,
    MS1CycleIndex,
    RTIndex,
    WindowCycleIndex,
};

/// CV accessions for the ion-mobility-type scan param.
const ACC_INVERSE_REDUCED_IM: u32 = 1002815; // 1/K0 — searchable
const ACC_FAIMS_CV: u32 = 1001581; // FAIMS compensation voltage
const ACC_DRIFT_TIME: u32 = 1002476; // drift time (not 1/K0)

/// Bucket size for peak groups, matching the TDF path.
const BUCKET_SIZE: usize = 4096;

/// The `1.0` mobility placeholder stored on every peak of a non-`Ook0` run —
/// present only to satisfy `IndexedPeakGroup::new` (non-NaN, ≥0); never used
/// once the mobility filter is unrestricted.
const MOBILITY_SENTINEL: f32 = 1.0;

/// Reader for open mzML. (Native Thermo `.raw` via mzdata's `thermo` feature is
/// deferred — it needs `MZReader::open_path`, not the mzML-only `MzMLReader`
/// used here, and would extend `sniff` to claim `.raw`.)
pub struct MzdataReader;

impl RawReader for MzdataReader {
    fn name(&self) -> &'static str {
        "mzdata"
    }

    fn sniff(&self, uri: &Uri) -> Sniff {
        // `.mzML.gz` is out of scope: it ends with `.gz`, not `.mzml`, so it is
        // correctly NOT claimed here (registry → UnknownFormat, loud).
        if path_ends_with(uri, ".mzml") {
            Sniff::Yes
        } else {
            Sniff::No
        }
    }

    fn manifest(&self, uri: &Uri) -> Manifest {
        Manifest {
            entry: uri.clone(),
            required: vec![uri.clone()],
            optional: vec![],
        }
    }

    fn caches_to_idx(&self) -> bool {
        false
    }

    fn file_extensions(&self) -> &'static [&'static str] {
        &["mzML", "mzml"]
    }

    fn read(
        &self,
        src: &ResolvedSource,
        cfg: &CentroidingConfig,
    ) -> Result<IndexedTimstofPeaks, ReadError> {
        from_mzml_file(&src.entry_path(), cfg)
    }
}

/// Detect the mobility axis from a scan's `ion_mobility_type` param accession.
/// Keys on the ACCESSION, not `has_ion_mobility()`: a real Astral file carries
/// `FAIMS compensation voltage` (value 0), so `has_ion_mobility()` returns true
/// and `ion_mobility()` returns 0.0 — both would mislead.
///
/// Returns `None` if the spectrum carries no scan yet (defer to a later
/// spectrum); `Some(Absent)` when a scan exists but declares no IM-type param.
fn detect_mobility_kind(spec: &impl SpectrumLike) -> Option<MobilityKind> {
    let scan = spec.acquisition().first_scan()?;
    Some(match scan.ion_mobility_type().and_then(|p| p.accession) {
        Some(ACC_INVERSE_REDUCED_IM) => MobilityKind::Ook0,
        Some(ACC_FAIMS_CV) => MobilityKind::Unsupported("FAIMS compensation voltage".to_string()),
        Some(ACC_DRIFT_TIME) => MobilityKind::Unsupported("drift time".to_string()),
        _ => MobilityKind::Absent,
    })
}

/// Absolute `(mz_start, mz_end)` isolation bounds for an MS2 spectrum, or `None`
/// if the window is unusable (skip + warn at the call site).
fn window_bounds(win: &mzdata::spectrum::IsolationWindow) -> Option<(f32, f32)> {
    use mzdata::spectrum::IsolationWindowState as S;
    match win.flags {
        // Complete/Explicit → lower/upper are ABSOLUTE m/z (do NOT do target±).
        S::Complete | S::Explicit => Some((win.lower_bound, win.upper_bound)),
        // Offset → bounds are relative to target.
        S::Offset => Some((win.target - win.lower_bound, win.target + win.upper_bound)),
        S::Unknown => None,
    }
}

/// Accumulates the peaks of one distinct isolation window across all cycles.
struct WindowAccum {
    mz_start: f32,
    mz_end: f32,
    peaks: Vec<IndexedPeak<WindowCycleIndex>>,
}

/// Build an `IndexedTimstofPeaks` from an mzML file. Peer to
/// `IndexedTimstofPeaks::from_timstof_file`.
pub fn from_mzml_file(
    path: &Path,
    _cfg: &CentroidingConfig,
) -> Result<IndexedTimstofPeaks, ReadError> {
    let file = std::fs::File::open(path)?;
    let reader = MzMLReader::new(BufReader::new(file));

    let mut mobility_kind: Option<MobilityKind> = None;
    let mut ms1_peaks: Vec<IndexedPeak<MS1CycleIndex>> = Vec::new();
    let mut cycle_rts_ms: Vec<u32> = Vec::new();
    // window key (rounded lower,upper) → slot in `windows`
    let mut window_index: HashMap<(i64, i64), usize> = HashMap::new();
    let mut windows: Vec<WindowAccum> = Vec::new();
    let mut cycle: Option<u32> = None;
    let mut warned_ms2_before_ms1 = false;
    let mut warned_unknown_window = false;

    for spec in reader {
        if spec.signal_continuity().is_profile() {
            return Err(ReadError::Build(format!(
                "profile-mode spectrum in {path:?}; only centroided mzML is supported \
                 (on-the-fly centroiding is not implemented)"
            )));
        }
        if mobility_kind.is_none() {
            // Latch on the first spectrum that actually carries a scan, so a
            // scanless leading spectrum can't wrongly lock in `Absent`.
            mobility_kind = detect_mobility_kind(&spec);
        }

        if spec.ms_level() == 1 {
            let rt_min = spec.start_time();
            if !rt_min.is_finite() {
                return Err(ReadError::Build(format!(
                    "{path:?}: MS1 scan has a non-finite start time ({rt_min}); cannot build a \
                     cycle→RT mapping"
                )));
            }
            let c = cycle.map_or(0, |c| c + 1);
            cycle = Some(c);
            let rt_ms = (rt_min * 60_000.0).round().max(0.0) as u32;
            cycle_rts_ms.push(rt_ms);
            push_peaks(&spec, MS1CycleIndex::new(c), &mut ms1_peaks);
        } else {
            let Some(c) = cycle else {
                if !warned_ms2_before_ms1 {
                    warn!("mzML: MS2 scan before any MS1; skipping (warm-up/calibration?)");
                    warned_ms2_before_ms1 = true;
                }
                continue;
            };
            let Some(prec) = spec.precursor() else {
                continue;
            };
            let Some((mz_start, mz_end)) = window_bounds(prec.isolation_window()) else {
                if !warned_unknown_window {
                    warn!("mzML: isolation window with Unknown bounds; skipping spectrum");
                    warned_unknown_window = true;
                }
                continue;
            };
            // Dedup windows by their bounds at millidalton (1e-3 m/z) precision —
            // fine enough to keep genuinely distinct DIA windows apart, coarse
            // enough that float jitter maps a recurring window to one bin.
            let key = (
                (mz_start as f64 * 1000.0).round() as i64,
                (mz_end as f64 * 1000.0).round() as i64,
            );
            let slot = *window_index.entry(key).or_insert_with(|| {
                windows.push(WindowAccum {
                    mz_start,
                    mz_end,
                    peaks: Vec::new(),
                });
                windows.len() - 1
            });
            push_peaks(&spec, WindowCycleIndex::new(c), &mut windows[slot].peaks);
        }
    }

    let mobility_kind = mobility_kind.unwrap_or(MobilityKind::Absent);
    if !matches!(mobility_kind, MobilityKind::Ook0) {
        warn!(
            "mzML mobility axis is not TIMS 1/K0 ({mobility_kind:?}); \
             the mobility filter and score feature are disabled for this run"
        );
    }

    if cycle_rts_ms.is_empty() || ms1_peaks.is_empty() {
        return Err(ReadError::Build(format!(
            "{path:?} yielded no MS1 cycles/peaks; nothing to index"
        )));
    }
    // `CycleToRTMapping` binary-searches the RT array, so non-monotonic RTs
    // would silently return wrong cycle indices. Refuse loudly instead.
    if let Some(w) = cycle_rts_ms.windows(2).find(|w| w[1] < w[0]) {
        return Err(ReadError::Build(format!(
            "{path:?}: cycle RTs are not monotonically increasing ({} ms then {} ms); \
             spectra are out of acquisition order",
            w[0], w[1]
        )));
    }
    warn_on_overlapping_windows(&windows);

    let (ms1_group, _stats) = IndexedPeakGroup::new(
        ms1_peaks,
        CycleToRTMapping::new(cycle_rts_ms.clone()),
        BUCKET_SIZE,
    );

    let mut ms2_window_groups = Vec::with_capacity(windows.len());
    for w in windows {
        if w.peaks.is_empty() {
            continue; // IndexedPeakGroup::new cannot handle an empty peak list
        }
        // A single isolation window → one ring, with a non-degenerate IM AABB
        // strictly enclosing the sentinel so `simplify_vw_preserve` cannot
        // collapse coincident IM corners.
        let scheme = QuadrupoleIsolationScheme::from_xxyy(std::iter::once((
            w.mz_start as f64,
            w.mz_end as f64,
            0.0,
            2.0,
        )));
        let (group, _stats) = IndexedPeakGroup::new(
            w.peaks,
            CycleToRTMapping::new(cycle_rts_ms.clone()),
            BUCKET_SIZE,
        );
        ms2_window_groups.push((scheme, group));
    }

    if ms2_window_groups.is_empty() {
        warn!(
            "mzML {path:?} produced no MS2 window groups; only MS1 will be queryable \
             (is this a DIA file?)"
        );
    }

    Ok(IndexedTimstofPeaks {
        ms2_window_groups,
        ms1_peaks: ms1_group,
        mobility_kind,
    })
}

/// Push a spectrum's centroid peaks into `out`, dropping non-finite m/z and
/// clamping negative intensities (both would trip `IndexedPeakGroup::new`).
fn push_peaks<T: RTIndex>(spec: &impl SpectrumLike, cycle_index: T, out: &mut Vec<IndexedPeak<T>>) {
    for p in spec.peaks().iter() {
        let mz = p.mz;
        if !mz.is_finite() {
            continue;
        }
        let intensity = p.intensity.max(0.0);
        out.push(IndexedPeak {
            mz: mz as f32,
            intensity,
            mobility_ook0: f16::from_f32(MOBILITY_SENTINEL),
            cycle_index,
        });
    }
}

/// First strictly-overlapping window pair, if any. Half-open: adjacent DIA
/// windows SHARE an edge (`upper == next.lower`), which is NOT an overlap — only
/// a strict interior intersection (`a.lo < b.hi && b.lo < a.hi`) counts, so a
/// normal contiguous DIA scheme does not trip it.
fn first_overlapping_pair(windows: &[WindowAccum]) -> Option<(usize, usize)> {
    for (i, a) in windows.iter().enumerate() {
        for (j, b) in windows.iter().enumerate().skip(i + 1) {
            if a.mz_start < b.mz_end && b.mz_start < a.mz_end {
                return Some((i, j));
            }
        }
    }
    None
}

fn warn_on_overlapping_windows(windows: &[WindowAccum]) {
    if let Some((i, j)) = first_overlapping_pair(windows) {
        let (a, b) = (&windows[i], &windows[j]);
        warn!(
            "mzML: overlapping isolation windows ([{:.3},{:.3}] vs [{:.3},{:.3}]); \
             staggered/demux acquisition is out of scope — results may be approximate",
            a.mz_start, a.mz_end, b.mz_start, b.mz_end
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mzdata::spectrum::{
        IsolationWindow,
        IsolationWindowState,
    };

    fn win(mz_start: f32, mz_end: f32) -> WindowAccum {
        WindowAccum {
            mz_start,
            mz_end,
            peaks: Vec::new(),
        }
    }

    #[test]
    fn window_bounds_complete_and_explicit_are_absolute() {
        // Thermo DIA lands `Complete`: lower/upper are ABSOLUTE m/z.
        let w = IsolationWindow {
            target: 401.43,
            lower_bound: 400.43,
            upper_bound: 402.43,
            flags: IsolationWindowState::Complete,
        };
        assert_eq!(window_bounds(&w), Some((400.43, 402.43)));
        let w2 = IsolationWindow {
            flags: IsolationWindowState::Explicit,
            ..w
        };
        assert_eq!(window_bounds(&w2), Some((400.43, 402.43)));
    }

    #[test]
    fn window_bounds_offset_is_relative_to_target() {
        let w = IsolationWindow {
            target: 500.0,
            lower_bound: 1.0,
            upper_bound: 1.5,
            flags: IsolationWindowState::Offset,
        };
        assert_eq!(window_bounds(&w), Some((499.0, 501.5)));
    }

    #[test]
    fn window_bounds_unknown_is_none() {
        let w = IsolationWindow {
            target: 500.0,
            lower_bound: 0.0,
            upper_bound: 0.0,
            flags: IsolationWindowState::Unknown,
        };
        assert_eq!(window_bounds(&w), None);
    }

    #[test]
    fn edge_touching_windows_do_not_overlap() {
        // Contiguous DIA: win N upper == win N+1 lower. Must NOT be flagged.
        let ws = vec![win(400.0, 402.0), win(402.0, 404.0), win(404.0, 406.0)];
        assert_eq!(first_overlapping_pair(&ws), None);
    }

    #[test]
    fn interior_overlap_is_detected() {
        let ws = vec![win(400.0, 403.0), win(402.0, 405.0)];
        assert_eq!(first_overlapping_pair(&ws), Some((0, 1)));
    }
}
