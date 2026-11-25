//! Joint subject/rival scoring after all Phase 3 batches and before competition.

use super::window_map::{
    self,
    Groups,
    Intervals,
    MAX_GROUPS_PER_ROW,
    Span,
    WindowMap,
};
use crate::data_sources::reference_library::{
    ExpectedIntensity,
    RefQuery,
    ReferenceLibrary,
};
use crate::errors::TimsSeekError;
use crate::rt_calibration::CalibrationResult;
use crate::scoring::blocks::shadow::ShadowScores;
use crate::scoring::pipeline::{
    Scorer,
    TOP_N_FRAGMENTS,
};
use crate::scoring::results::ScoredCandidate;
use crate::traits::ScorerQueriable;
use half::f16;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use smallvec::SmallVec;
use std::cmp::Ordering;
use std::fmt::Write;
use std::num::NonZeroU16;
use timscentroid::indexing::IndexedPeakGroup;
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
    WindowCycleIndex,
};
use timscentroid::utils::OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use timsquery::chemistry::analyte::PropertyRef;
use timsquery::models::FlatIdx;
use timsquery::models::tolerance::MzTolerance;
use timsquery::{
    PeakTolerance,
    Target,
    TupleRange,
};

const MAX_HALF_WIDTH: u16 = 16;
const MAX_POPPED: u8 = 64;
const MAX_PROBES: u8 = 8;
const MAX_CHANNELS: usize = 2 * TOP_N_FRAGMENTS;
const MAX_BUCKET_HEADS: usize = 2 * MAX_HALF_WIDTH as usize + 1;
const INVALID_GROUP: u16 = u16::MAX;
const _: () = assert!(TOP_N_FRAGMENTS == 8 && MAX_CHANNELS == 16);

#[derive(Debug, Clone, Copy)]
pub(super) struct ShadowContext {
    flat: FlatIdx,
    apex_cycle: u32,
    radius: NonZeroU16,
    group: u16,
    mobility: [f32; 2],
}

const _: () = assert!(size_of::<Option<ShadowContext>>() == 20);

/// Phase 2.5 geometry and calibrated peak tolerances for an eager run index.
pub struct ShadowRun<'a> {
    windows: WindowMap<'a>,
    tolerance: PeakTolerance,
    overlap_tolerance: MzTolerance,
    n_cycles: u32,
    n_offsets: u32,
}

pub(super) fn bounds(message: &'static str) -> TimsSeekError {
    TimsSeekError::ShadowBounds { message }
}

fn cycle_radius(fwhm: f32) -> NonZeroU16 {
    assert!(fwhm.is_finite());
    let radius = fwhm.ceil().clamp(1.0, f32::from(MAX_HALF_WIDTH)) as u16;
    NonZeroU16::new(radius).expect("radius is clamped above zero")
}

fn global_apex(offset: usize, local: usize) -> u32 {
    u32::try_from(offset.checked_add(local).expect("apex cycle addition fits"))
        .expect("apex lies on the u32 cycle axis")
}

impl<I: ScorerQueriable> Scorer<I> {
    pub fn prepare_shadows(
        &self,
        calibration: Option<&CalibrationResult>,
    ) -> Result<Option<ShadowRun<'_>>, TimsSeekError> {
        let Some(index) = self.index.eager_index() else {
            return Ok(None);
        };
        let windows = WindowMap::new(index)?;
        if windows.is_empty() {
            return Ok(None);
        }
        let n_cycles = u32::try_from(index.ms1_cycle_mapping().len())
            .map_err(|_| bounds("too many MS1 cycles"))?;
        let n_groups = u32::try_from(windows.group_count()).expect("checked group count");
        let n_offsets = n_groups
            .checked_mul(n_cycles)
            .and_then(|n| n.checked_add(1))
            .filter(|&n| n < u32::MAX)
            .ok_or_else(|| bounds("too many group/cycle offsets"))?;
        let tolerance = match calibration {
            Some(calibration) => calibration.unrestricted_rt_tolerance().peak_tolerance(),
            None => self.broad_tolerance.peak_tolerance(),
        };
        Ok(Some(ShadowRun {
            windows,
            overlap_tolerance: recentered(tolerance.ms2),
            tolerance,
            n_cycles,
            n_offsets,
        }))
    }
}

impl ShadowRun<'_> {
    pub(super) fn context(
        &self,
        query: &RefQuery<'_>,
        flat: FlatIdx,
        offset: usize,
        local: usize,
        fwhm: f32,
    ) -> Option<ShadowContext> {
        let (group, mobility) = if self.windows.mobility {
            let mobility = self
                .tolerance
                .mobility_range_f16(query.mobility_ook0())
                .map(|range| [range.start().to_f32(), range.end().to_f32()])
                .unwrap_or([f32::NEG_INFINITY, f32::INFINITY]);
            let group = match self
                .windows
                .mobility_groups(self.precursor_band(query), span(mobility))
            {
                Ok(groups) => groups?.0,
                // Phase 3.5 reports this before allocating the membership index.
                Err(_) => INVALID_GROUP,
            };
            (group, mobility)
        } else {
            (
                self.windows.static_group(query.mono_precursor_mz())?,
                [f32::NEG_INFINITY, f32::INFINITY],
            )
        };
        let apex_cycle = global_apex(offset, local);
        assert!(apex_cycle < self.n_cycles);
        Some(ShadowContext {
            flat,
            apex_cycle,
            radius: cycle_radius(fwhm),
            group,
            mobility,
        })
    }

    fn precursor_band(&self, query: &RefQuery<'_>) -> Span {
        let (lo, hi) = query.precursor_mz_limits();
        let mz = self.tolerance.quad_range_f32((lo as f32, hi as f32));
        (f64::from(mz.start()), f64::from(mz.end()))
    }

    fn memberships(&self, context: ShadowContext, library: &ReferenceLibrary) -> Groups {
        if !self.windows.mobility {
            return Groups::from_slice(&[context.group]);
        }
        let query = library.item_at(context.flat);
        let (primary, groups) = self
            .windows
            .mobility_groups(self.precursor_band(&query), span(context.mobility))
            .expect("Phase 3 validated mobility group count")
            .expect("Phase 3 found a mobility group");
        assert_eq!(primary, context.group);
        groups
    }

    fn bucket(&self, group: u16, cycle: u32) -> usize {
        (u32::from(group) * self.n_cycles + cycle) as usize
    }

    /// Score a whole run at once, preserving Phase 3 row order and main scores.
    pub fn score(
        &self,
        candidates: &mut [ScoredCandidate],
        library: &ReferenceLibrary,
    ) -> Result<(), TimsSeekError> {
        let order = RowOrder {
            candidates,
            library,
        };
        let index = self.build_csr(&order)?;
        let targets = candidates
            .iter()
            .filter(|row| row.scoring.identity.is_target)
            .count();
        let decoys = candidates.len() - targets;
        let mut targets_with_rival = 0;
        let mut decoys_with_rival = 0;
        let largest = (0..self.windows.group_count())
            .map(|id| index.group(self.bucket(id as u16, 0), self.n_cycles).len())
            .max()
            .unwrap_or(0);
        let mut scores = Vec::with_capacity(largest);
        #[cfg(not(feature = "rayon"))]
        let mut scratch = SubjectScratch::default();
        for group in 0..self.windows.group_count() {
            let group = u16::try_from(group).expect("checked group count");
            let bucket = self.bucket(group, 0);
            let members = index.group(bucket, self.n_cycles);
            let order = RowOrder {
                candidates,
                library,
            };
            let score = |scratch: &mut SubjectScratch, &row: &u32| {
                self.score_subject(&index, &order, row, group, scratch)
            };
            #[cfg(feature = "rayon")]
            members
                .par_iter()
                .map_init(SubjectScratch::default, score)
                .collect_into_vec(&mut scores);
            #[cfg(not(feature = "rayon"))]
            {
                scores.clear();
                scores.extend(members.iter().map(|row| score(&mut scratch, row)));
            }
            assert_eq!(members.len(), scores.len());
            for (&row, &(score, has_rival)) in members.iter().zip(&scores) {
                let candidate = &mut candidates[row as usize];
                if candidate
                    .shadow_context
                    .is_some_and(|context| context.group == group)
                {
                    candidate.scoring.shadow = score;
                    if has_rival {
                        match candidate.scoring.identity.is_target {
                            true => targets_with_rival += 1,
                            false => decoys_with_rival += 1,
                        }
                    }
                }
            }
        }
        // Missing = NaN features, including rows without a shadow context.
        let (targets_missing, decoys_missing) =
            candidates.iter().fold((0usize, 0usize), |(t, d), row| {
                match (
                    row.scoring.shadow.shadow_cosine.is_nan(),
                    row.scoring.identity.is_target,
                ) {
                    (true, true) => (t + 1, d),
                    (true, false) => (t, d + 1),
                    (false, true | false) => (t, d),
                }
            });
        tracing::info!(
            targets_with_rival,
            targets_without_rival = targets - targets_with_rival,
            decoys_with_rival,
            decoys_without_rival = decoys - decoys_with_rival,
            targets_missing,
            decoys_missing,
            "Phase 3.5 shadow rival counts"
        );
        Ok(())
    }

    fn build_csr(&self, order: &RowOrder<'_>) -> Result<Csr, TimsSeekError> {
        let candidates = order.candidates;
        if candidates.len() >= u32::MAX as usize {
            return Err(bounds("too many scored rows"));
        }
        let mut n_members = 0u32;
        for context in candidates.iter().filter_map(|row| row.shadow_context) {
            if usize::from(context.group) >= self.windows.group_count() {
                return Err(bounds("invalid subject group or too many memberships"));
            }
            if context.apex_cycle >= self.n_cycles {
                return Err(bounds("shadow apex outside the MS1 cycle axis"));
            }
            let groups = self.memberships(context, order.library);
            assert!(groups.len() <= MAX_GROUPS_PER_ROW);
            assert!(
                groups
                    .iter()
                    .enumerate()
                    .all(|(i, id)| !groups[..i].contains(id))
            );
            n_members = n_members
                .checked_add(groups.len() as u32)
                .filter(|&n| n < u32::MAX)
                .ok_or_else(|| bounds("too many shadow memberships"))?;
        }
        // The preflight above checks the exact size before either CSR allocation.
        let mut offsets = vec![0u32; self.n_offsets as usize];
        let mut rows = vec![0u32; n_members as usize];
        for context in candidates.iter().filter_map(|row| row.shadow_context) {
            for group in self.memberships(context, order.library) {
                offsets[self.bucket(group, context.apex_cycle) + 1] += 1;
            }
        }
        for index in 1..offsets.len() {
            offsets[index] += offsets[index - 1];
        }
        assert_eq!(offsets.last().copied(), Some(n_members));
        for (row, candidate) in candidates.iter().enumerate() {
            let Some(context) = candidate.shadow_context else {
                continue;
            };
            for group in self.memberships(context, order.library) {
                let cursor = &mut offsets[self.bucket(group, context.apex_cycle)];
                rows[*cursor as usize] = u32::try_from(row).expect("checked row count");
                *cursor += 1;
            }
        }
        offsets.copy_within(..self.n_offsets as usize - 1, 1);
        offsets[0] = 0;
        assert!(offsets.windows(2).all(|pair| pair[0] <= pair[1]));
        assert_eq!(offsets.last().copied(), Some(n_members));
        let mut remaining = rows.as_mut_slice();
        for pair in offsets.windows(2) {
            let (bucket, rest) = remaining.split_at_mut((pair[1] - pair[0]) as usize);
            bucket.sort_unstable_by(|&left, &right| order.compare(left, right));
            assert!(
                bucket
                    .windows(2)
                    .all(|pair| !order.compare(pair[0], pair[1]).is_gt())
            );
            remaining = rest;
        }
        assert!(remaining.is_empty());
        Ok(Csr { offsets, rows })
    }

    fn score_subject(
        &self,
        index: &Csr,
        order: &RowOrder<'_>,
        row: u32,
        group: u16,
        scratch: &mut SubjectScratch,
    ) -> (ShadowScores, bool) {
        let subject = order.context(row);
        if subject.group != group {
            return (ShadowScores::default(), false);
        }
        let cycle = WindowCycleIndex::from(MS1CycleIndex::new(subject.apex_cycle));
        let peaks = self.windows.peaks(group);
        if !peaks.has_frame(cycle) {
            return (ShadowScores::default(), false);
        }

        let query = order.library.item_at(subject.flat);
        let fragments = top_fragments(&query, self.tolerance.ms2, self.overlap_tolerance);
        if fragments.is_empty() {
            return (ShadowScores::default(), false);
        }

        if self.windows.mobility {
            self.windows.project(
                group,
                self.precursor_band(&query),
                span(subject.mobility),
                &mut scratch.subject_im,
            );
            assert!(!scratch.subject_im.is_empty());
        } else {
            assert!(scratch.subject_im.is_empty());
        }
        let cycles = cycle_bounds(subject.apex_cycle, subject.radius, self.n_cycles);
        let start = self.bucket(group, cycles.0);
        let end = self.bucket(group, cycles.1);
        let rival = find_rival(
            &mut scratch.heads,
            index,
            start..=end,
            row,
            |left, right| order.compare(left, right),
            |rival| {
                let context = order.context(rival);
                assert!(
                    context.apex_cycle.abs_diff(subject.apex_cycle)
                        <= u32::from(subject.radius.get())
                );
                let query = order.library.item_at(context.flat);
                if self.windows.mobility {
                    self.windows.project(
                        group,
                        self.precursor_band(&query),
                        span(context.mobility),
                        &mut scratch.rival_im,
                    );
                    assert!(!scratch.rival_im.is_empty());
                    if !window_map::intersects(&scratch.subject_im, &scratch.rival_im) {
                        return None;
                    }
                } else {
                    assert_eq!(context.group, group);
                }
                Some(top_fragments(
                    &query,
                    self.tolerance.ms2,
                    self.overlap_tolerance,
                ))
            },
            &fragments,
        );
        let scores = score_at_apex(
            peaks,
            cycle,
            &fragments,
            rival.as_deref().unwrap_or(&[]),
            &scratch.subject_im,
            &mut scratch.channels,
        );
        (scores, rival.is_some())
    }
}

fn span([lo, hi]: [f32; 2]) -> Span {
    (f64::from(lo), f64::from(hi))
}

fn cycle_bounds(cycle: u32, radius: NonZeroU16, n_cycles: u32) -> (u32, u32) {
    assert!(cycle < n_cycles);
    assert!(radius.get() <= MAX_HALF_WIDTH);
    let radius = u32::from(radius.get());
    (
        cycle.saturating_sub(radius),
        cycle.saturating_add(radius).min(n_cycles - 1),
    )
}

struct Csr {
    offsets: Vec<u32>,
    rows: Vec<u32>,
}

impl Csr {
    fn group(&self, start: usize, n_cycles: u32) -> &[u32] {
        &self.rows[self.offsets[start] as usize..self.offsets[start + n_cycles as usize] as usize]
    }
}

struct RowOrder<'a> {
    candidates: &'a [ScoredCandidate],
    library: &'a ReferenceLibrary,
}

impl RowOrder<'_> {
    fn context(&self, row: u32) -> ShadowContext {
        self.candidates[row as usize]
            .shadow_context
            .expect("only rows with shadow context enter CSR")
    }

    // Ascending ordering places the highest score first.
    fn compare(&self, left: u32, right: u32) -> Ordering {
        let left_score = self.candidates[left as usize].scoring.primary.main_score;
        let right_score = self.candidates[right as usize].scoring.primary.main_score;
        right_score.total_cmp(&left_score).then_with(|| {
            let left = self.library.item_at(self.context(left).flat);
            let right = self.library.item_at(self.context(right).flat);
            compare_ties(&left, &right)
        })
    }
}

#[derive(Default)]
struct SequenceHash {
    value: u64,
    tail: [u8; 8],
    used: u8,
}

impl Write for SequenceHash {
    fn write_str(&mut self, value: &str) -> std::fmt::Result {
        for byte in value.bytes() {
            self.tail[usize::from(self.used)] = byte;
            self.used += 1;
            if self.used == 8 {
                self.value = splitmix64(self.value ^ u64::from_le_bytes(self.tail));
                self.used = 0;
                self.tail = [0; 8];
            }
        }
        Ok(())
    }
}

fn splitmix64(input: u64) -> u64 {
    let mut z = input.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

impl SequenceHash {
    fn finish(self, charge: u8) -> u64 {
        let hash = if self.used == 0 {
            self.value
        } else {
            splitmix64(self.value ^ u64::from_le_bytes(self.tail))
        };
        splitmix64(hash ^ u64::from(charge))
    }
}

fn sequence_key(query: &RefQuery<'_>) -> (u64, u64) {
    let geom = query.geom();
    let mut hash = SequenceHash::default();
    let annotation = match geom.geom().analyte(geom.row()).peptide {
        PropertyRef::Known(peptide) => match peptide.modifications {
            PropertyRef::Known(_) => {
                assert!(u32::try_from(peptide.residues.len()).is_ok());
                peptide
                    .write_sequence(&mut hash)
                    .expect("known modified sequence");
                None
            }
            PropertyRef::Unresolved { annotation, .. } => Some(annotation),
            PropertyRef::Missing | PropertyRef::NotApplicable => None,
        },
        PropertyRef::Unresolved { annotation, .. } => Some(annotation),
        PropertyRef::Missing | PropertyRef::NotApplicable => None,
    };
    if let Some(annotation) = annotation {
        hash.write_str(annotation).expect("hash sink is infallible");
    }
    (
        hash.finish(query.precursor_charge()),
        query.mono_precursor_mz().to_bits(),
    )
}

fn compare_ties(left: &RefQuery<'_>, right: &RefQuery<'_>) -> Ordering {
    sequence_key(left).cmp(&sequence_key(right)).then_with(|| {
        let left = top_predictions(left);
        let right = top_predictions(right);
        left.iter()
            .map(Prediction::bits)
            .cmp(right.iter().map(Prediction::bits))
    })
}

#[derive(Clone, Copy)]
struct BucketHead {
    cursor: u32,
    end: u32,
}

type BucketHeads = SmallVec<[BucketHead; MAX_BUCKET_HEADS]>;

#[derive(Default)]
struct SubjectScratch {
    heads: BucketHeads,
    channels: ApexChannels,
    subject_im: Intervals,
    rival_im: Intervals,
}

fn sift_down(
    heads: &mut [BucketHead],
    mut parent: usize,
    higher: impl Fn(&BucketHead, &BucketHead) -> bool,
) {
    for _ in 0..MAX_BUCKET_HEADS.ilog2() {
        let left = 2 * parent + 1;
        if left >= heads.len() {
            return;
        }
        let right = left + 1;
        let best = if right < heads.len() && higher(&heads[right], &heads[left]) {
            right
        } else {
            left
        };
        if !higher(&heads[best], &heads[parent]) {
            return;
        }
        heads.swap(parent, best);
        parent = best;
    }
    assert!(2 * parent + 1 >= heads.len());
}

fn find_rival(
    heads: &mut BucketHeads,
    index: &Csr,
    buckets: std::ops::RangeInclusive<usize>,
    subject: u32,
    compare: impl Fn(u32, u32) -> Ordering,
    mut eligible_fragments: impl FnMut(u32) -> Option<Fragments>,
    fragments: &[Fragment],
) -> Option<Fragments> {
    assert!(buckets.clone().count() <= MAX_BUCKET_HEADS);
    heads.clear();
    for bucket in buckets {
        let cursor = index.offsets[bucket];
        let end = index.offsets[bucket + 1];
        if cursor < end {
            assert!(heads.len() < MAX_BUCKET_HEADS);
            heads.push(BucketHead { cursor, end });
        }
    }
    let higher = |left: &BucketHead, right: &BucketHead| {
        compare(
            index.rows[left.cursor as usize],
            index.rows[right.cursor as usize],
        )
        .is_lt()
    };
    for parent in (0..heads.len() / 2).rev() {
        sift_down(heads, parent, higher);
    }
    let mut probes = 0;
    for _ in 0..MAX_POPPED {
        let Some(head) = heads.first_mut() else {
            break;
        };
        let row = index.rows[head.cursor as usize];
        head.cursor += 1;
        if head.cursor == head.end {
            heads.swap_remove(0);
        }
        sift_down(heads, 0, higher);
        if row == subject {
            continue;
        }
        if let Some(rival) = eligible_fragments(row) {
            probes += 1;
            assert!(probes <= MAX_PROBES);
            if has_overlap(fragments, &rival) {
                assert_ne!(row, subject);
                return Some(rival);
            }
            if probes == MAX_PROBES {
                break;
            }
        }
    }
    None
}

fn recentered(tolerance: MzTolerance) -> MzTolerance {
    match tolerance {
        MzTolerance::Absolute((lo, hi)) => {
            let radius = lo.midpoint(hi);
            MzTolerance::Absolute((radius, radius))
        }
        MzTolerance::Ppm((lo, hi)) => {
            let radius = lo.midpoint(hi);
            MzTolerance::Ppm((radius, radius))
        }
    }
}

#[derive(Clone, Copy)]
struct Prediction {
    mz: f32,
    weight: f32,
}

impl Prediction {
    fn bits(&self) -> (u32, u32) {
        (self.mz.to_bits(), self.weight.to_bits())
    }
}

fn top_predictions(query: &RefQuery<'_>) -> SmallVec<[Prediction; TOP_N_FRAGMENTS]> {
    let mut top = SmallVec::<[Prediction; TOP_N_FRAGMENTS]>::new();
    for ((label, mz), (expected_label, weight)) in query
        .iter_fragments_refs()
        .zip(query.iter_expected_fragments())
    {
        debug_assert_eq!(*label, expected_label);
        if !weight.is_finite() || weight <= 0.0 {
            continue;
        }
        let mz = mz as f32;
        assert!(mz.is_finite() && mz > 0.0);
        let position = top.partition_point(|fragment| {
            fragment.weight > weight
                || (fragment.weight == weight && fragment.mz.to_bits() <= mz.to_bits())
        });
        if position < TOP_N_FRAGMENTS {
            if top.len() == TOP_N_FRAGMENTS {
                top.pop();
            }
            top.insert(position, Prediction { mz, weight });
        }
    }
    assert!(top.len() <= TOP_N_FRAGMENTS);
    top
}

#[derive(Clone, Copy)]
struct Fragment {
    mz: f32,
    predicted: f64,
    range: TupleRange<f32>,
    overlap_range: TupleRange<f32>,
}

impl Fragment {
    fn new(prediction: Prediction, tolerance: MzTolerance, overlap_tolerance: MzTolerance) -> Self {
        let Prediction { mz, weight } = prediction;
        Self {
            mz,
            predicted: f64::from(weight),
            range: tolerance.mz_range_f32(mz),
            overlap_range: overlap_tolerance.mz_range_f32(mz),
        }
    }

    fn shares_channel(&self, other: &Self) -> bool {
        ((self.overlap_range.start() <= other.mz) & (other.mz <= self.overlap_range.end()))
            | ((other.overlap_range.start() <= self.mz) & (self.mz <= other.overlap_range.end()))
    }
}

type Fragments = SmallVec<[Fragment; TOP_N_FRAGMENTS]>;

fn has_overlap(subject: &[Fragment], rival: &[Fragment]) -> bool {
    assert!(subject.len() <= TOP_N_FRAGMENTS);
    assert!(rival.len() <= TOP_N_FRAGMENTS);
    let chunks = subject.chunks_exact(4);
    let tail = chunks.remainder();
    for other in rival {
        let mut hit = false;
        for chunk in chunks.clone() {
            hit |= chunk[0].shares_channel(other)
                | chunk[1].shares_channel(other)
                | chunk[2].shares_channel(other)
                | chunk[3].shares_channel(other);
        }
        for fragment in tail {
            hit |= fragment.shares_channel(other);
        }
        if hit {
            return true;
        }
    }
    false
}

fn top_fragments(
    query: &RefQuery<'_>,
    tolerance: MzTolerance,
    overlap_tolerance: MzTolerance,
) -> Fragments {
    top_predictions(query)
        .into_iter()
        .map(|prediction| Fragment::new(prediction, tolerance, overlap_tolerance))
        .collect()
}

struct ApexChannel {
    mz: f32,
    range: TupleRange<f32>,
    predicted: [f64; 2],
}

type ApexChannels = SmallVec<[ApexChannel; MAX_CHANNELS]>;

fn build_channels(fragments: &[Fragment], rival: &[Fragment], channels: &mut ApexChannels) {
    assert!(fragments.len() <= TOP_N_FRAGMENTS);
    assert!(rival.len() <= TOP_N_FRAGMENTS);
    channels.clear();
    channels.extend(fragments.iter().map(|fragment| ApexChannel {
        mz: fragment.mz,
        range: fragment.range,
        predicted: [fragment.predicted, 0.0],
    }));
    for other in rival {
        let shared = fragments
            .iter()
            .enumerate()
            .filter(|(_, fragment)| fragment.shares_channel(other))
            .min_by(|(i, left), (j, right)| {
                distance_order(other.mz, left.mz, right.mz).then(i.cmp(j))
            });
        if let Some((index, _)) = shared {
            let channel = &mut channels[index];
            channel.predicted[1] += other.predicted;
            // Retain the candidate center and cover both extraction ranges.
            channel.range = TupleRange::try_new(
                channel.range.start().min(other.range.start()),
                channel.range.end().max(other.range.end()),
            )
            .expect("union covers two ordered fragment ranges");
        } else {
            assert!(channels.len() < MAX_CHANNELS);
            channels.push(ApexChannel {
                mz: other.mz,
                range: other.range,
                predicted: [0.0, other.predicted],
            });
        }
    }
    assert!((fragments.len()..=fragments.len() + rival.len()).contains(&channels.len()));
}

fn distance_order(mz: f32, left: f32, right: f32) -> Ordering {
    (mz - left).abs().total_cmp(&(mz - right).abs())
}

/// Assign a centroid to the nearest eligible channel across both predictions.
fn channel(channels: &[ApexChannel], mz: f32) -> Option<usize> {
    channels
        .iter()
        .enumerate()
        .filter(|(_, channel)| channel.range.contains(mz))
        .min_by(|(i, left), (j, right)| distance_order(mz, left.mz, right.mz).then(i.cmp(j)))
        .map(|(index, _)| index)
}

fn score_at_apex(
    peaks: &IndexedPeakGroup<WindowCycleIndex>,
    cycle: WindowCycleIndex,
    fragments: &[Fragment],
    rival: &[Fragment],
    mobility: &[Span],
    channels: &mut ApexChannels,
) -> ShadowScores {
    build_channels(fragments, rival, channels);
    let observed = observed_at_apex(peaks, cycle, channels, mobility);
    residual_scores(channels, &observed[..channels.len()], !rival.is_empty())
}

fn observed_at_apex(
    peaks: &IndexedPeakGroup<WindowCycleIndex>,
    cycle: WindowCycleIndex,
    channels: &[ApexChannel],
    mobility: &[Span],
) -> [f64; MAX_CHANNELS] {
    let mut observed = [0.0; MAX_CHANNELS];
    if !peaks.has_frame(cycle) {
        return observed;
    }
    let mut rounded;
    let mobility = if mobility.is_empty() {
        mobility
    } else {
        // The index stores centroid mobilities rounded to f16.
        rounded = Intervals::new();
        assert!(mobility.len() <= rounded.inline_size());
        rounded.extend(
            mobility
                .iter()
                .map(|&(lo, hi)| (f16::from_f64(lo).to_f64(), f16::from_f64(hi).to_f64())),
        );
        rounded.as_slice()
    };
    let cycles = TupleRange::try_new(cycle, cycle).expect("one cycle is an ordered range");
    let mobility_envelope = match (mobility.first(), mobility.last()) {
        (Some(first), Some(last)) => Restricted(
            TupleRange::try_new(f16::from_f64(first.0), f16::from_f64(last.1))
                .expect("mobility projection is sorted"),
        ),
        _ => Unrestricted,
    };
    assert!(channels.len() <= MAX_CHANNELS);
    let mut ranges: SmallVec<[(u8, TupleRange<f32>); MAX_CHANNELS]> = (0u8..)
        .zip(channels.iter().map(|channel| channel.range))
        .collect();
    ranges.sort_unstable_by(|a, b| a.1.start().total_cmp(&b.1.start()));
    let mut scan = |lo, hi, single: Option<u8>| {
        let mz = TupleRange::try_new(lo, hi).expect("merged channel ranges are ordered");
        peaks.for_each_peak(mz, Restricted(cycles), mobility_envelope, |peak| {
            if mobility.is_empty() || window_map::contains(mobility, peak.mobility_ook0.to_f64()) {
                let index = single.map_or_else(
                    || channel(channels, peak.mz).expect("merged range covers the centroid"),
                    usize::from,
                );
                observed[index] += f64::from(peak.intensity);
            }
        });
    };
    if let Some((&(first_id, first), rest)) = ranges.split_first() {
        let (mut lo, mut hi) = (first.start(), first.end());
        let mut single = Some(first_id);
        for &(index, range) in rest {
            if range.start() > hi {
                scan(lo, hi, single);
                lo = range.start();
                hi = range.end();
                single = Some(index);
            } else {
                hi = hi.max(range.end());
                single = None;
            }
        }
        scan(lo, hi, single);
    }
    observed
}

/// Solve the two-variable NNLS over all channels, including unmatched ones.
fn fit_scales(channels: &[ApexChannel], observed: &[f64]) -> [f64; 2] {
    assert_eq!(channels.len(), observed.len());
    let (mut aa, mut ab, mut bb, mut ay, mut by) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (channel, &y) in channels.iter().zip(observed) {
        let [a, b] = channel.predicted;
        aa += a * a;
        ab += a * b;
        bb += b * b;
        ay += a * y;
        by += b * y;
    }
    assert!(aa > 0.0);
    assert!(bb > 0.0);
    let determinant = aa * bb - ab * ab;
    let scales = if determinant > 8.0 * f64::EPSILON * aa * bb {
        let a = (bb * ay - ab * by) / determinant;
        let b = (aa * by - ab * ay) / determinant;
        if a >= 0.0 && b >= 0.0 {
            Some([a, b])
        } else {
            None
        }
    } else {
        None
    }
    .unwrap_or_else(|| {
        // Near-singular systems also use an axis approximation. Prefer the
        // subject when the error reductions tie.
        let ay = ay.max(0.0);
        let by = by.max(0.0);
        if ay * ay / aa >= by * by / bb {
            [ay / aa, 0.0]
        } else {
            [0.0, by / bb]
        }
    });
    assert!(
        scales
            .iter()
            .all(|&scale| scale.is_finite() && scale >= 0.0)
    );
    scales
}

fn residual_scores(channels: &[ApexChannel], observed: &[f64], has_rival: bool) -> ShadowScores {
    assert_eq!(channels.len(), observed.len());
    let rival_scale = if has_rival {
        fit_scales(channels, observed)[1]
    } else {
        0.0
    };
    let mut dot = 0.0;
    let mut norm_predicted = 0.0;
    let mut norm_residual = 0.0;
    let mut intensity_observed = 0.0;
    let mut intensity_residual = 0.0;
    let mut lost = 0u8;
    let mut present = 0u8;
    for (channel, &before) in channels.iter().zip(observed) {
        let [candidate_predicted, rival_predicted] = channel.predicted;
        if candidate_predicted == 0.0 {
            continue;
        }
        let residual = before - rival_scale * rival_predicted;
        // Treat rounding residue from an exactly fitted channel as zero.
        let after = if residual > before * (8.0 * f64::EPSILON) {
            residual
        } else {
            0.0
        };
        assert!(after >= 0.0);
        assert!(after <= before);
        dot += candidate_predicted * after;
        norm_predicted += candidate_predicted * candidate_predicted;
        norm_residual += after * after;
        intensity_observed += before;
        intensity_residual += after;
        lost += u8::from(before > 0.0 && after == 0.0);
        present += u8::from(before > 0.0);
    }
    if intensity_observed == 0.0 {
        assert_eq!(present, 0);
        return ShadowScores::default();
    }
    assert!(lost <= present);
    assert!(present > 0);
    let cosine = if norm_residual > 0.0 {
        (dot / (norm_predicted * norm_residual).sqrt()).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let scores = ShadowScores {
        shadow_cosine: cosine as f32,
        shadow_retained: (intensity_residual / intensity_observed) as f32,
        shadow_lost: f32::from(lost) / f32::from(present),
    };
    assert!(
        [
            scores.shadow_cosine,
            scores.shadow_retained,
            scores.shadow_lost
        ]
        .iter()
        .all(|score| (0.0..=1.0).contains(score))
    );
    scores
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DecoyPolicy;
    use crate::scoring::results::ScoringFields;
    use timscentroid::geometry::QuadrupoleIsolationScheme;
    use timscentroid::indexing::IndexedPeak;
    use timscentroid::rt_mapping::CycleToRTMapping;
    use timscentroid::{
        IndexedTimstofPeaks,
        MobilityKind,
    };
    use timsquery::models::{
        Row,
        TargetCapabilities,
        TargetColumnsBuilder,
    };
    use timsquery::serde::TargetTable;

    const ABS: MzTolerance = MzTolerance::Absolute((0.010, 0.010));

    fn union_channels(fragments: &[Fragment], rival: &[Fragment]) -> ApexChannels {
        let mut channels = ApexChannels::new();
        build_channels(fragments, rival, &mut channels);
        channels
    }

    fn score_at_apex(
        peaks: &IndexedPeakGroup<WindowCycleIndex>,
        cycle: WindowCycleIndex,
        fragments: &[Fragment],
        rival: &[Fragment],
        mobility: &[Span],
    ) -> ShadowScores {
        super::score_at_apex(
            peaks,
            cycle,
            fragments,
            rival,
            mobility,
            &mut ApexChannels::new(),
        )
    }

    fn fragment(mz: f32, weight: f32) -> Fragment {
        Fragment::new(Prediction { mz, weight }, ABS, ABS)
    }

    fn peaks(
        values: &[(f32, f32, f32, u32)],
        frames: &[u32],
    ) -> IndexedPeakGroup<WindowCycleIndex> {
        let peaks = values
            .iter()
            .map(|&(mz, intensity, mobility, cycle)| IndexedPeak {
                mz,
                intensity,
                mobility_ook0: f16::from_f32(mobility),
                cycle_index: WindowCycleIndex::new(cycle),
            })
            .collect();
        IndexedPeakGroup::testing_new(peaks, CycleToRTMapping::new(frames.to_vec()), 4).0
    }

    fn assert_scores(scores: ShadowScores, expected: [f32; 3]) {
        let actual = [
            scores.shadow_cosine,
            scores.shadow_retained,
            scores.shadow_lost,
        ];
        for (actual, expected) in actual.into_iter().zip(expected) {
            if expected.is_nan() {
                assert!(actual.is_nan());
            } else {
                assert!((actual - expected).abs() < 1e-6, "{actual} != {expected}");
            }
        }
    }

    #[test]
    fn joint_fit_leaves_the_known_subject_contribution() {
        let subject = [fragment(300.0, 1.0), fragment(400.0, 1.0)];
        let rival = [fragment(400.005, 1.0), fragment(500.0, 1.0)];
        let channels = union_channels(&subject, &rival);
        let observed = [2.0, 5.0, 3.0];
        assert_eq!(channels.len(), 3);
        assert_eq!(fit_scales(&channels, &observed), [2.0, 3.0]);
        let index = peaks(
            &[
                (300.0, 2.0, 1.0, 0),
                (400.003, 5.0, 1.0, 0),
                (500.0, 3.0, 1.0, 0),
            ],
            &[0],
        );
        assert_scores(
            score_at_apex(&index, WindowCycleIndex::new(0), &subject, &rival, &[]),
            [1.0, 4.0 / 7.0, 0.0],
        );

        // The older rival-only one-axis fit removes four units instead of three.
        let rival_only = (observed[1] + observed[2]) / 2.0;
        let residual = [observed[0], observed[1] - rival_only];
        assert_eq!(rival_only, 4.0);
        assert_eq!(residual, [2.0, 1.0]);
        let cosine = residual.iter().sum::<f64>()
            / (2.0 * residual.iter().map(|value| value * value).sum::<f64>()).sqrt();
        assert_eq!(cosine, 3.0 / 10.0_f64.sqrt());
        assert_eq!(
            residual.iter().sum::<f64>() / (observed[0] + observed[1]),
            3.0 / 7.0
        );
    }

    #[test]
    fn negative_interior_uses_the_rival_axis() {
        let subject = [fragment(300.0, 1.0), fragment(400.0, 1.0)];
        let rival = [fragment(400.0, 1.0), fragment(500.0, 1.0)];
        let channels = union_channels(&subject, &rival);
        let observed = [0.0, 2.0, 6.0];
        assert_eq!((2.0 * 2.0 - 8.0) / 3.0, -4.0 / 3.0);
        assert_eq!(fit_scales(&channels, &observed), [0.0, 4.0]);
        assert_scores(residual_scores(&channels, &observed, true), [0.0, 0.0, 1.0]);
    }

    #[test]
    fn fully_removed_channel_counts_toward_lost() {
        let subject = [fragment(300.0, 1.0), fragment(400.0, 1.0)];
        let rival = [fragment(300.0, 1.0), fragment(500.0, 1.0)];
        let channels = union_channels(&subject, &rival);
        let scores = residual_scores(&channels, &[2.0, 4.0, 8.0], true);
        assert_eq!(scores.shadow_lost, 0.5);
        assert_scores(scores, [1.0 / 2.0_f32.sqrt(), 4.0 / 6.0, 0.5]);
    }

    #[test]
    fn proportional_and_near_singular_ties_prefer_the_subject_axis() {
        let subject = [fragment(300.0, 1.0), fragment(400.0, 2.0)];
        let rival = [fragment(300.0, 2.0), fragment(400.0, 4.0)];
        let channels = union_channels(&subject, &rival);
        let observed = [3.0, 6.0];
        // union_channels follows the supplied subject order.
        assert_eq!(fit_scales(&channels, &observed), [3.0, 0.0]);
        assert_scores(residual_scores(&channels, &observed, true), [1.0, 1.0, 0.0]);

        let subject = [fragment(300.0, 1.0)];
        let epsilon = 2.0_f32.powi(-25);
        let rival = [fragment(300.0, 1.0), fragment(400.0, epsilon)];
        let channels = union_channels(&subject, &rival);
        let observed = [2.0, f64::from(epsilon)];
        assert_eq!(fit_scales(&channels, &observed), [2.0, 0.0]);
        assert_scores(residual_scores(&channels, &observed, true), [1.0, 1.0, 0.0]);
    }

    #[test]
    fn overlap_requires_a_centre_inside_the_recentered_band() {
        let subject = [fragment(300.0, 1.0); 8];
        let band_overlap_only = [fragment(300.015, 1.0)];
        assert!(subject[0].range.intersects(band_overlap_only[0].range));
        assert!(!has_overlap(&subject, &band_overlap_only));
        assert!(has_overlap(&subject, &[fragment(300.008, 1.0)]));
        let ppm = MzTolerance::Ppm((-8.0, 15.0));
        let abs = |mz| Fragment::new(Prediction { mz, weight: 1.0 }, ppm, recentered(ppm));
        let anchor = abs(1000.0);
        assert_eq!(anchor.range.as_tuple(), (1000.008, 1000.015));
        assert!(anchor.shares_channel(&abs(1000.003)));
        assert!(!anchor.shares_channel(&abs(1000.005)));
    }

    #[test]
    fn merged_ranges_assign_each_centroid_once_to_the_nearest_eligible_channel() {
        let subject = [fragment(100.0, 1.0)];
        let rival = [fragment(100.012, 1.0)];
        let channels = union_channels(&subject, &rival);
        assert_eq!(channel(&channels, 100.004), Some(0));
        // At exactly equal distances the first eligible channel wins.
        let delta = 1.0 / 256.0;
        let ties = [
            ApexChannel {
                mz: 100.0 - delta,
                ..channels[0]
            },
            ApexChannel {
                mz: 100.0 + delta,
                ..channels[0]
            },
        ];
        assert_eq!(channel(&ties, 100.0), Some(0));

        for input in [
            [(100.004, 2.0, 1.0, 0), (100.008, 3.0, 1.0, 0)],
            [(100.008, 3.0, 1.0, 0), (100.004, 2.0, 1.0, 0)],
        ] {
            let index = peaks(&input, &[0]);
            let observed = observed_at_apex(&index, WindowCycleIndex::new(0), &channels, &[]);
            assert_eq!(&observed[..2], &[2.0, 3.0]);
            assert_scores(
                score_at_apex(&index, WindowCycleIndex::new(0), &subject, &rival, &[]),
                [1.0, 1.0, 0.0],
            );
        }
    }

    #[test]
    fn channel_distance_is_gated_by_calibrated_ranges() {
        let shifted = MzTolerance::Absolute((-0.008, 0.015));
        let make = |mz| Fragment::new(Prediction { mz, weight: 1.0 }, shifted, recentered(shifted));
        let channels = union_channels(&[make(100.000), make(100.010)], &[]);
        assert_eq!(channel(&channels, 100.011), Some(0));
    }

    #[test]
    fn ceil_radius_and_global_cycle_include_the_right_frames() {
        let radius = cycle_radius(1.9);
        let apex = global_apex(2, 1);
        assert_eq!(radius.get(), 2);
        assert_eq!(apex, 3);
        let (lo, hi) = cycle_bounds(apex, radius, 7);
        assert_eq!((lo, hi), (1, 5));
        let csr = Csr {
            offsets: vec![0, 0, 0, 0, 1, 1, 2, 3],
            rows: vec![0, 1, 2],
        };
        let subject = [fragment(300.0, 1.0)];
        let rival = find_rival(
            &mut BucketHeads::new(),
            &csr,
            lo as usize..=hi as usize,
            0,
            |left, right| right.cmp(&left),
            |row| {
                assert_eq!(row, 1, "row 2 is stronger but three cycles away");
                Some(Fragments::from_slice(&subject))
            },
            &subject,
        );
        assert!(rival.is_some());
        let cycle = WindowCycleIndex::from(MS1CycleIndex::new(apex));
        let channels = union_channels(&[fragment(300.0, 1.0)], &[]);
        let group = peaks(
            &[(300.0, 2.0, 1.0, 3), (300.0, 100.0, 1.0, 4)],
            &[0, 100, 201, 302, 402, 402],
        )
        .with_frame_presence(vec![false, true, false, true, true, false]);
        assert!(group.has_frame(cycle));
        assert_eq!(observed_at_apex(&group, cycle, &channels, &[])[0], 2.0);
        assert_scores(
            score_at_apex(
                &group,
                WindowCycleIndex::new(2),
                &[fragment(300.0, 1.0)],
                &[],
                &[],
            ),
            [f32::NAN; 3],
        );
    }

    #[test]
    fn rival_probe_and_pop_limits_charge_different_rows() {
        let single_bucket = |n| Csr {
            offsets: vec![0, n],
            rows: (0..n).collect(),
        };
        let subject = [fragment(300.0, 1.0)];
        let match_at = |row, target| {
            std::iter::once(if row == target {
                subject[0]
            } else {
                fragment(400.0, 1.0)
            })
            .collect()
        };

        // Self precedes the eight failed eligible probes. The ninth cannot win.
        let index = single_bucket(10);
        let mut tried = Vec::with_capacity(9);
        assert!(
            find_rival(
                &mut BucketHeads::new(),
                &index,
                0..=0,
                0,
                |a, b| a.cmp(&b),
                |row| {
                    tried.push(row);
                    Some(match_at(row, 9))
                },
                &subject
            )
            .is_none()
        );
        assert_eq!(tried, (1..=8).collect::<Vec<_>>());
        assert!(
            find_rival(
                &mut BucketHeads::new(),
                &index,
                0..=0,
                0,
                |a, b| a.cmp(&b),
                |row| { Some(match_at(row, 8)) },
                &subject
            )
            .is_some()
        );

        let subject_im = [(0.92, 0.95), (1.05, 1.08)];
        let ineligible_im = [(0.97, 1.03)];
        assert!(
            find_rival(
                &mut BucketHeads::new(),
                &index,
                0..=0,
                0,
                |a, b| a.cmp(&b),
                |row| {
                    let mobility = if row < 9 {
                        &ineligible_im[..]
                    } else {
                        &subject_im[..]
                    };
                    window_map::intersects(&subject_im, mobility).then(|| match_at(row, 9))
                },
                &subject
            )
            .is_some()
        );

        let index = single_bucket(65);
        for (rival, expected) in [(63, true), (64, false)] {
            assert_eq!(
                find_rival(
                    &mut BucketHeads::new(),
                    &index,
                    0..=0,
                    0,
                    |a, b| a.cmp(&b),
                    |row| { (row == rival).then(|| match_at(row, rival)) },
                    &subject
                )
                .is_some(),
                expected
            );
        }
    }

    fn library(weights: &[[f32; 2]], mobility: &[f32]) -> ReferenceLibrary {
        assert_eq!(weights.len(), mobility.len());
        let mut builder =
            TargetColumnsBuilder::with_capabilities(TargetCapabilities::default_diann());
        let fragments = [
            (timsquery::IonAnnot::try_from("y1").unwrap(), 300.0),
            (timsquery::IonAnnot::try_from("y2").unwrap(), 400.0),
        ];
        for &mobility in mobility {
            builder.push_row(Row {
                precursor_mz: 500.0,
                charge: 2,
                mobility,
                frags: &fragments,
                ..Default::default()
            });
        }
        let geom = builder.seal(DecoyPolicy::Never).unwrap();
        ReferenceLibrary::try_from(TargetTable::Mzpaf {
            geom,
            frag_intens: Some(weights.iter().flatten().copied().collect()),
        })
        .unwrap()
    }

    fn index(
        groups: Vec<(
            QuadrupoleIsolationScheme,
            IndexedPeakGroup<WindowCycleIndex>,
        )>,
        mobility: MobilityKind,
    ) -> IndexedTimstofPeaks {
        let ms1 = [IndexedPeak {
            mz: 500.0,
            intensity: 1.0,
            mobility_ook0: f16::ONE,
            cycle_index: MS1CycleIndex::new(0),
        }];
        let (ms1, _) = IndexedPeakGroup::testing_new(
            ms1.to_vec(),
            CycleToRTMapping::new(vec![0, 100, 200, 300]),
            4,
        );
        IndexedTimstofPeaks::from_parts(ms1, groups).with_mobility_kind(mobility)
    }

    #[test]
    fn score_ties_compare_top_fragment_bits_and_survive_row_permutation() {
        let library = library(&[[1.0, 0.0], [2.0, 0.5], [2.0, 1.0]], &[1.0; 3]);
        let flats: Vec<_> = library.geometry().flats().collect();
        let preferred = library.item_at(flats[1]);
        let alternate = library.item_at(flats[2]);
        assert_eq!(
            sequence_key(&preferred),
            (splitmix64(2), 500.0_f64.to_bits())
        );
        assert_eq!(sequence_key(&preferred), sequence_key(&alternate));
        assert_eq!(compare_ties(&preferred, &alternate), Ordering::Less);

        let group = peaks(&[(300.0, 4.0, 1.0, 1), (400.0, 1.0, 1.0, 1)], &[0, 100]);
        let window = QuadrupoleIsolationScheme::from_xxyy([(490.0, 510.0, 0.0, 2.0)].into_iter());
        let scorer = Scorer {
            index: index(vec![(window, group)], MobilityKind::Absent),
            broad_tolerance: timsquery::Tolerance {
                ms2: ABS,
                ..Default::default()
            },
            fragmented_range: (490.0, 510.0).try_into().unwrap(),
            denoise_traces: false,
        };
        let run = scorer.prepare_shadows(None).unwrap().unwrap();
        let make = |row: usize| {
            let flat = flats[row];
            let mut candidate = ScoredCandidate::new(ScoringFields::sample_default());
            candidate.scoring.primary.main_score = if row == 0 { 1.0 } else { 2.0 };
            // Put tied rivals in different CSR buckets to exercise the heap tie too.
            let local = if row == 2 { 2 } else { 1 };
            candidate.shadow_context = run.context(&library.item_at(flat), flat, 0, local, 1.9);
            candidate
        };
        for permutation in [[0, 1, 2], [0, 2, 1], [2, 0, 1]] {
            let mut rows = permutation.map(make);
            run.score(&mut rows, &library).unwrap();
            let scores = rows
                .iter()
                .find(|row| row.scoring.primary.main_score == 1.0)
                .unwrap();
            assert_scores(scores.scoring.shadow, [0.0, 0.0, 1.0]);
        }

        let mut alternate_only = [make(0), make(2)];
        run.score(&mut alternate_only, &library).unwrap();
        assert_scores(alternate_only[0].scoring.shadow, [1.0, 0.5, 0.0]);
    }

    #[test]
    fn sequence_hash_pads_the_final_little_endian_chunk() {
        let mut hash = SequenceHash::default();
        hash.write_str("ABCDE").unwrap();
        hash.write_str("FGHI").unwrap();
        let expected = splitmix64(u64::from_le_bytes(*b"ABCDEFGH"));
        let expected = splitmix64(expected ^ u64::from(b'I'));
        assert_eq!(hash.finish(3), splitmix64(expected ^ 3));
    }

    #[test]
    fn disjoint_mobility_ranges_mask_centroids_and_other_group_rivals() {
        let group_0 = QuadrupoleIsolationScheme::from_xxyy(
            [(490.0, 510.0, 0.90, 0.95), (490.0, 510.0, 1.05, 1.10)].into_iter(),
        );
        let group_1 =
            QuadrupoleIsolationScheme::from_xxyy([(490.0, 510.0, 0.97, 1.03)].into_iter());
        let mz = (499.0, 501.0);
        let bounds = (0.92, 1.08);
        let mut subject_im = Intervals::new();
        window_map::project(&group_0, mz, bounds, &mut subject_im);
        assert_eq!(subject_im.as_slice(), &[(0.92, 0.95), (1.05, 1.08)]);

        let group = peaks(
            &[
                (300.0, 2.0, 0.94, 0),
                (300.0, 100.0, 1.00, 0),
                (300.0, 3.0, 1.07, 0),
            ],
            &[0],
        );
        let subject = [fragment(300.0, 1.0)];
        let channels = union_channels(&subject, &[]);
        assert_eq!(
            observed_at_apex(&group, WindowCycleIndex::new(0), &channels, &subject_im)[0],
            5.0
        );

        let mut rival_im = Intervals::new();
        window_map::project(&group_0, mz, (0.97, 1.03), &mut rival_im);
        assert!(!window_map::intersects(&subject_im, &rival_im));
        window_map::project(&group_1, mz, (0.97, 1.03), &mut rival_im);
        assert!(!rival_im.is_empty());
        let memberships = Csr {
            offsets: vec![0, 1, 2],
            rows: vec![0, 1],
        };
        assert!(
            find_rival(
                &mut BucketHeads::new(),
                &memberships,
                0..=0,
                0,
                |a, b| a.cmp(&b),
                |_| { panic!("self and other group must not be probed") },
                &subject
            )
            .is_none()
        );
    }

    #[test]
    fn mobility_mask_includes_centroids_on_quantized_edges() {
        let geometry = QuadrupoleIsolationScheme::from_xxyy(
            [(490.0, 510.0, 0.90, 0.95), (490.0, 510.0, 1.05, 1.10)].into_iter(),
        );
        let mut intervals = Intervals::new();
        window_map::project(&geometry, (499.0, 501.0), (0.89, 1.11), &mut intervals);
        assert!(f16::from_f64(0.90).to_f64() < 0.90);
        assert!(f16::from_f64(0.95).to_f64() > 0.95);
        let group = peaks(
            &[
                (300.0, 2.0, 0.90, 0),
                (300.0, 3.0, 0.95, 0),
                (300.0, 100.0, 0.96, 0),
            ],
            &[0],
        );
        let channels = union_channels(&[fragment(300.0, 1.0)], &[]);
        let observed = observed_at_apex(&group, WindowCycleIndex::new(0), &channels, &intervals);
        assert_eq!(observed[0], 5.0);
    }

    #[test]
    fn mobility_run_fits_a_rival_from_its_secondary_group() {
        use timsquery::models::tolerance::MobilityTolerance;

        let library = library(&[[1.0, 0.0], [2.0, 0.5], [2.0, 1.0]], &[1.0, 1.13, 1.45]);
        let group_0 = QuadrupoleIsolationScheme::from_xxyy(
            [
                (490.0, 510.0, 0.90, 0.95),
                (490.0, 510.0, 1.05, 1.10),
                (490.0, 510.0, 1.40, 1.50),
            ]
            .into_iter(),
        );
        let group_1 =
            QuadrupoleIsolationScheme::from_xxyy([(490.0, 510.0, 1.07, 1.20)].into_iter());
        let primary = peaks(
            &[
                (300.0, 4.0, 0.94, 1),
                (300.0, 100.0, 1.00, 1),
                (400.0, 1.0, 1.07, 1),
            ],
            &[0, 100],
        );
        let secondary = peaks(&[(300.0, 10.0, 1.07, 1)], &[0, 100]);
        let scorer = Scorer {
            index: index(
                vec![(group_0, primary), (group_1, secondary)],
                MobilityKind::Ook0,
            ),
            broad_tolerance: timsquery::Tolerance {
                ms2: ABS,
                mobility: MobilityTolerance::Absolute((0.08, 0.08)),
                ..Default::default()
            },
            fragmented_range: (490.0, 510.0).try_into().unwrap(),
            denoise_traces: false,
        };
        let run = scorer.prepare_shadows(None).unwrap().unwrap();
        let mut rows: Vec<_> = library
            .geometry()
            .flats()
            .enumerate()
            .map(|(row, flat)| {
                let mut candidate = ScoredCandidate::new(ScoringFields::sample_default());
                candidate.scoring.primary.main_score = [1.0, 2.0, 3.0][row];
                candidate.shadow_context = run.context(&library.item_at(flat), flat, 0, 1, 1.9);
                candidate
            })
            .collect();
        let groups: Vec<_> = rows
            .iter()
            .map(|row| row.shadow_context.unwrap().group)
            .collect();
        assert_eq!(groups, [0, 1, 0]);
        let rival_groups = run.memberships(rows[1].shadow_context.unwrap(), &library);
        assert_eq!(rival_groups.as_slice(), &[0, 1]);
        run.score(&mut rows, &library).unwrap();
        assert_scores(rows[0].scoring.shadow, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn sloped_projection_excludes_peaks_inside_only_the_bounding_box() {
        let geometry: QuadrupoleIsolationScheme = serde_json::from_value(serde_json::json!([
            {"Trapezoid": {"im": [0.80, 1.20], "mz_lo": [480.0, 520.0], "mz_hi": [505.0, 545.0]}},
            {"Aabb": {"mz": [480.0, 545.0], "im": [1.30, 1.40]}}
        ]))
        .unwrap();
        let mut intervals = Intervals::new();
        window_map::project(&geometry, (510.0, 511.0), (0.80, 1.40), &mut intervals);
        assert_eq!(intervals.len(), 2);
        assert!((intervals[0].0 - 0.85).abs() < 1e-12);
        assert!((intervals[0].1 - 1.11).abs() < 1e-12);
        assert_eq!(intervals[1], (1.30, 1.40));
        let group = peaks(
            &[
                (300.0, 2.0, 0.90, 0),
                (300.0, 100.0, 1.15, 0),
                (300.0, 200.0, 1.25, 0),
            ],
            &[0],
        );
        let channels = union_channels(&[fragment(300.0, 1.0)], &[]);
        assert_eq!(
            observed_at_apex(&group, WindowCycleIndex::new(0), &channels, &intervals)[0],
            2.0
        );
    }
}
