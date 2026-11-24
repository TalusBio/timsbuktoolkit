//! Acquisition cycles of a diaPASEF run, validated once before indexing.
//!
//! A cycle starts at an MS1 frame and runs to the next one. Indexing assumes every
//! cycle holds the MS1 frame and a frame for every MS2 window group, so that cycle
//! indices line up across MS1 and all window groups. Instruments occasionally write
//! an incomplete cycle mid-run (an MS1 frame followed directly by the next MS1).
//! Such a cycle is dropped everywhere and the remaining cycles are renumbered, so
//! the indexer's invariant holds. Dropping leaves a wider RT step between the
//! neighbouring cycles, which is what was measured; filling with empty frames would
//! claim zero intensity where nothing was acquired.

use timsrust::core::MSLevel;
use timsrust::tdf::TdfFrameReader;
use tracing::{
    debug,
    warn,
};

/// More interior cycles than this fraction dropped means the data is broken, not
/// hiccuping: refuse to index it.
pub const DROPPED_FRACTION_MAX: f64 = 0.01;

/// First cycle index handed out. Index 0 is never a cycle (its RT slot stays 0).
const FIRST_CYCLE: u32 = 1;

#[derive(Debug, Clone)]
pub struct CycleTable {
    /// Cycle index per frame index; `None` for frames in a dropped cycle.
    cycle_of_frame: Vec<Option<u32>>,
    /// Incomplete cycles between the first and last complete cycle.
    pub dropped_interior: Vec<DroppedCycle>,
    pub n_complete: u32,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DroppedCycle {
    pub first_frame: usize,
    pub rt_seconds: f64,
    pub has_ms1: bool,
    pub window_groups_present: u32,
    pub window_groups_expected: u32,
}

#[derive(Debug, thiserror::Error)]
pub enum CycleError {
    #[error("no complete acquisition cycle (MS1 + every MS2 window group) in the run")]
    NoCompleteCycles,
    #[error(
        "{dropped} of {total} cycles are incomplete (max {max_percent}%); first at \
         frame {first_frame}, {first_rt_seconds:.2} s"
    )]
    TooManyDropped {
        dropped: usize,
        total: usize,
        max_percent: f64,
        first_frame: usize,
        first_rt_seconds: f64,
    },
    #[error("cannot read frame metadata: {0}")]
    Metadata(String),
}

/// One frame's role in cycle assembly: the facts `CycleTable::build` needs.
#[derive(Debug, Clone, Copy)]
pub enum FrameKind {
    Ms1,
    Ms2 { window_group: u8 },
}

impl CycleTable {
    pub fn from_frame_reader(frame_reader: &TdfFrameReader) -> Result<Self, CycleError> {
        let mut indices: Vec<usize> = frame_reader.iter_indices().collect();
        indices.sort_unstable();
        let mut frames = Vec::with_capacity(indices.len());
        for index in indices {
            let meta = frame_reader
                .get_info(index)
                .map_err(|e| CycleError::Metadata(format!("frame {index}: {e:?}")))?;
            let kind = match meta.ms_level() {
                MSLevel::MS1 => FrameKind::Ms1,
                MSLevel::MS2 => FrameKind::Ms2 {
                    window_group: meta.window_group(),
                },
                MSLevel::Unknown => {
                    return Err(CycleError::Metadata(format!(
                        "frame {index}: unknown MS level"
                    )));
                }
            };
            frames.push((index, kind, meta.rt_in_seconds()));
        }
        Self::build(&frames)
    }

    /// Assign cycle indices to `frames` (frame index, kind, RT seconds), sorted by
    /// frame index.
    pub fn build(frames: &[(usize, FrameKind, f64)]) -> Result<Self, CycleError> {
        assert!(
            frames.windows(2).all(|w| w[0].0 < w[1].0),
            "frames must be sorted"
        );
        let window_groups_expected = {
            let mut groups: Vec<u8> = frames
                .iter()
                .filter_map(|(_, kind, _)| match kind {
                    FrameKind::Ms1 => None,
                    FrameKind::Ms2 { window_group } => Some(*window_group),
                })
                .collect();
            groups.sort_unstable();
            groups.dedup();
            groups
        };
        let raw_cycles = split_at_ms1(frames);
        let complete: Vec<bool> = raw_cycles
            .iter()
            .map(|cycle| is_complete(cycle, &window_groups_expected))
            .collect();
        let Some(first) = complete.iter().position(|c| *c) else {
            return Err(CycleError::NoCompleteCycles);
        };
        let last = complete
            .iter()
            .rposition(|c| *c)
            .expect("a first implies a last");

        let size = frames.last().map_or(0, |(index, _, _)| index + 1);
        let mut cycle_of_frame = vec![None; size];
        let mut dropped_interior = Vec::with_capacity(raw_cycles.len());
        let mut next = FIRST_CYCLE;
        for (position, (cycle, is_complete)) in raw_cycles.iter().zip(&complete).enumerate() {
            if *is_complete {
                for (index, _, _) in *cycle {
                    cycle_of_frame[*index] = Some(next);
                }
                next += 1;
                continue;
            }
            let dropped = describe(cycle, window_groups_expected.len());
            if first < position && position < last {
                warn!(
                    first_frame = dropped.first_frame,
                    rt_seconds = dropped.rt_seconds,
                    has_ms1 = dropped.has_ms1,
                    window_groups_present = dropped.window_groups_present,
                    window_groups_expected = dropped.window_groups_expected,
                    "Dropping incomplete acquisition cycle"
                );
                dropped_interior.push(dropped);
            } else {
                // Before the first or after the last complete cycle: run start/end.
                debug!(?dropped, "Dropping partial cycle at run edge");
            }
        }
        let n_complete = next - FIRST_CYCLE;
        let total = last - first + 1;
        let fraction = dropped_interior.len() as f64 / total as f64;
        if fraction > DROPPED_FRACTION_MAX {
            let first_dropped = &dropped_interior[0];
            return Err(CycleError::TooManyDropped {
                dropped: dropped_interior.len(),
                total,
                max_percent: DROPPED_FRACTION_MAX * 100.0,
                first_frame: first_dropped.first_frame,
                first_rt_seconds: first_dropped.rt_seconds,
            });
        }
        assert_eq!(n_complete as usize + dropped_interior.len(), total);
        Ok(Self {
            cycle_of_frame,
            dropped_interior,
            n_complete,
        })
    }

    /// Cycle index of a frame, or `None` if its cycle was dropped.
    pub fn cycle_of(&self, frame_index: usize) -> Option<u32> {
        self.cycle_of_frame.get(frame_index).copied().flatten()
    }
}

type Frame = (usize, FrameKind, f64);

/// Split frames into raw cycles at each MS1 frame. Frames before the first MS1
/// form a leading cycle without MS1.
fn split_at_ms1(frames: &[Frame]) -> Vec<&[Frame]> {
    let mut cycles = Vec::with_capacity(frames.len());
    let mut start = 0;
    for (position, (_, kind, _)) in frames.iter().enumerate() {
        match kind {
            FrameKind::Ms1 if position > start => {
                cycles.push(&frames[start..position]);
                start = position;
            }
            FrameKind::Ms1 | FrameKind::Ms2 { .. } => {}
        }
    }
    if start < frames.len() {
        cycles.push(&frames[start..]);
    }
    cycles
}

fn is_complete(cycle: &[Frame], window_groups_expected: &[u8]) -> bool {
    let has_ms1 = matches!(cycle.first(), Some((_, FrameKind::Ms1, _)));
    has_ms1
        && window_groups_expected.iter().all(|group| {
            cycle.iter().any(|(_, kind, _)| match kind {
                FrameKind::Ms1 => false,
                FrameKind::Ms2 { window_group } => window_group == group,
            })
        })
}

fn describe(cycle: &[Frame], window_groups_expected: usize) -> DroppedCycle {
    let (first_frame, first_kind, rt_seconds) = cycle[0];
    let mut present: Vec<u8> = cycle
        .iter()
        .filter_map(|(_, kind, _)| match kind {
            FrameKind::Ms1 => None,
            FrameKind::Ms2 { window_group } => Some(*window_group),
        })
        .collect();
    present.sort_unstable();
    present.dedup();
    DroppedCycle {
        first_frame,
        rt_seconds,
        has_ms1: matches!(first_kind, FrameKind::Ms1),
        window_groups_present: u32::try_from(present.len()).expect("window groups fit u32"),
        window_groups_expected: u32::try_from(window_groups_expected)
            .expect("window groups fit u32"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A run of `n` cycles of MS1 + groups 1..=3, with `skip` cycles missing all MS2.
    fn run(n: usize, skip: &[usize]) -> Vec<Frame> {
        let mut frames = Vec::with_capacity(n * 4);
        for cycle in 0..n {
            frames.push((frames.len(), FrameKind::Ms1, cycle as f64));
            if !skip.contains(&cycle) {
                for window_group in 1..=3 {
                    frames.push((frames.len(), FrameKind::Ms2 { window_group }, cycle as f64));
                }
            }
        }
        frames
    }

    #[test]
    fn complete_run_numbers_cycles_from_one() {
        let table = CycleTable::build(&run(5, &[])).unwrap();
        assert_eq!(table.n_complete, 5);
        assert!(table.dropped_interior.is_empty());
        assert_eq!(table.cycle_of(0), Some(1));
        assert_eq!(table.cycle_of(19), Some(5));
    }

    #[test]
    fn interior_ms1_only_cycle_is_dropped_and_later_cycles_shift_down() {
        // Mirrors 261020_HELA_200ng_postcal: one MS1 frame with no MS2 frames.
        let frames = run(200, &[100]);
        let table = CycleTable::build(&frames).unwrap();
        assert_eq!(table.n_complete, 199);
        assert_eq!(table.dropped_interior.len(), 1);
        assert_eq!(table.dropped_interior[0].window_groups_present, 0);
        let orphan_ms1 = 100 * 4;
        assert_eq!(table.cycle_of(orphan_ms1), None);
        assert_eq!(table.cycle_of(orphan_ms1 - 1), Some(100));
        assert_eq!(table.cycle_of(orphan_ms1 + 1), Some(101));
    }

    #[test]
    fn partial_cycles_at_run_edges_are_dropped_without_counting() {
        let mut frames = vec![(0, FrameKind::Ms2 { window_group: 2 }, 0.0)];
        frames.extend(run(3, &[]).into_iter().map(|(i, k, rt)| (i + 1, k, rt)));
        frames.push((13, FrameKind::Ms1, 9.0));
        frames.push((14, FrameKind::Ms2 { window_group: 1 }, 9.0));
        let table = CycleTable::build(&frames).unwrap();
        assert_eq!(table.n_complete, 3);
        assert!(table.dropped_interior.is_empty());
        assert_eq!(table.cycle_of(0), None);
        assert_eq!(table.cycle_of(14), None);
    }

    #[test]
    fn too_many_interior_drops_is_an_error() {
        let result = CycleTable::build(&run(50, &[10, 20]));
        assert!(matches!(
            result,
            Err(CycleError::TooManyDropped {
                dropped: 2,
                total: 50,
                ..
            })
        ));
    }

    #[test]
    fn run_without_a_complete_cycle_is_an_error() {
        let frames = [
            (0, FrameKind::Ms2 { window_group: 1 }, 0.0),
            (1, FrameKind::Ms2 { window_group: 2 }, 0.1),
        ];
        assert!(matches!(
            CycleTable::build(&frames),
            Err(CycleError::NoCompleteCycles)
        ));
    }
}
