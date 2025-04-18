use crate::utils::streaming_calculators::RunningStatsCalculator;
use crate::errors::Result;
use crate::traits::KeyLike;
use nohash_hasher::BuildNoHashHasher;
use std::collections::HashMap;
use std::sync::Arc;


pub struct ImsMzRollingCalculator {
    pub mobility: RunningStatsCalculator,
    pub mz: RunningStatsCalculator,
}

impl ImsMzRollingCalculator {
    fn new(intensity: u64, mobility: f32, mz: f32) -> Self {
        Self {
            mobility: RunningStatsCalculator::new(intensity, mobility),
            mz: RunningStatsCalculator::new(intensity, mz),
        }
    }
    fn add(&mut self, intensity: u64, mobility: f32, mz: f32) {
        self.mobility.add(intensity, mobility);
        self.mz.add(intensity, mz);
    }
}

/// Hashmap that represents a series of scan+tof values
/// that are grouped by retention time represented in ms (u32)
pub struct SparseRTCollection {
    inner:HashMap<u32, ImsMzRollingCalculator, BuildNoHashHasher<u32>>,
} 

impl SparseRTCollection {
    fn add(&mut self, rt_ms: u32, mobility: f32, mz: f32, intensity: u64) {
        self.inner.entry(rt_ms)
            .and_modify(|curr| {
                curr.add(intensity, mobility, mz);
            })
            .or_insert(ImsMzRollingCalculator::new(
                intensity, mobility, mz,
            ));
    }
}

/// Represents "there is a value at each positions in the retention time"
#[derive(Debug, Clone)]
pub struct DenseRTCollection {
    // QUESTION: should I keep this an arc and do the binary search myself or
    // abstract to some `rt_to_index_mapper`?
    // 2025-Apr-16 JSPP
    pub reference_rt_ms: Arc<[u32]>,
    pub scan_tof_calc: Vec<Option<ImsMzRollingCalculator>>,
}

impl DenseRTCollection {
    fn new(reference_rt_ms: Arc<[u32]>) -> Self {
        // Check that its sorted ... let make it panic for now
        // Since I am not using this anywhere where this would be a recoverable
        // error.
        assert!(
            reference_rt_ms.is_sorted(),
            "DenseRTCollection::new reference_rt_ms must be sorted"
        );

        let num_frames = reference_rt_ms.len();
        Self {
            reference_rt_ms,
            scan_tof_calc: vec![None; num_frames],
        }
    }
    fn add(&mut self, rt_ms: u32, mobility: f32, mz: f32, intensity: u64) {
        let mut pos = self.reference_rt_ms.partition_point(|&x| x <= rt_ms);
        if pos == self.reference_rt_ms.len() {
            // If we are less than 1 second above the end, we can just subtract 1.
            let diff = rt_ms - self.reference_rt_ms[pos - 1];
            if diff < 1_000 {
                pos -= 1;
            }
        }
        match self.scan_tof_calc.get_mut(pos) {
            Some(x) => match x {
                Some(x) => {
                    x.add(intensity, mobility, mz);
                }
                None => {
                    self.scan_tof_calc[pos] = Some(ImsMzRollingCalculator::new(
                        intensity, mobility, mz,
                    ))
                }
            },
            None => {
                let max_pos = self.scan_tof_calc.len() - 1;
                println!(
                    "DenseRTCollection::add out of bounds {:?} > {:?}",
                    pos, max_pos
                );
                // TODO: Decide if this is really the behavior I want.
                panic!()
            }
        }
    }
}


/// Represents a series of scan+tof values that are grouped by retention time represented in ms (u32)
/// Each element in the vec is an independent 'track' (each is a precursor/transition)
#[derive(Debug, Clone)]
pub struct ParallelTracks {
    // Sparse(Vec<SparseRTCollection>),
    tracks: Vec<DenseRTCollection>
}

impl ParallelTracks {
    fn with_capacity(num_tracks: usize, reference_rt_ms: Option<Arc<[u32]>>) -> Self {
        match reference_rt_ms {
            Some(rt) => {
                let mut scan_tof_calc = Vec::with_capacity(num_tracks);
                for _ in 0..num_tracks {
                    scan_tof_calc.push(DenseRTCollection::new(rt.clone()));
                }
                Self::Dense(scan_tof_calc)
            }
            None => {
                let mut scan_tof_calc = Vec::with_capacity(num_tracks);
                for _ in 0..num_tracks {
                    scan_tof_calc
                        .push(SparseRTCollection::with_hasher(BuildNoHashHasher::default()));
                }
                Self::Sparse(scan_tof_calc)
            }
        }
    }

    fn add(&mut self, idx: usize, rt_ms: u32, mobility: f32, tof_index: f32, intensity: u64) {
        match self {
            Self::Sparse(sparse) => sparse[idx].add(rt_ms, mobility, tof_index, intensity),
            Self::Dense(dense) => dense[idx].add(rt_ms, mobility, tof_index, intensity),
        }
    }

    fn iter_mut(&mut self) -> impl Iterator<Item = &mut SparseRTCollection> {
        match self {
            Self::Sparse(sparse) => sparse.iter_mut(),
            Self::Dense(dense) => dense.iter_mut(),
        }
    }
}


#[derive(Debug, Clone)]
pub struct ParitionedCMGAggregator<FH: KeyLike> {
    pub scan_tof_calc: ParallelTracks,
    pub keys: Vec<FH>,
    pub context_key_num: usize,
    pub expected_scan_index: usize,
    pub expected_tof_indices: Vec<u32>,
}

impl<FH: KeyLike> ParitionedCMGAggregator<FH>
{
    pub fn new(
        keys: Vec<FH>,
        expected_scan_index: usize,
        expected_tof_indices: Vec<u32>,
        expected_rt_ms: Option<Arc<[u32]>>,
    ) -> Self {
        assert!(expected_tof_indices.len() == keys.len());
        let scan_tof_calc = ParallelTracks::with_capacity(keys.len(), expected_rt_ms);

        Self {
            scan_tof_calc,
            keys,
            context_key_num: 0,
            expected_scan_index,
            expected_tof_indices,
        }
    }

    pub fn set_context(&mut self, context: FH) -> Result<()> {
        // Find the position in the keys.
        let pos = self.keys.iter().position(|x| x == &context);
        if let Some(pos) = pos {
            self.context_key_num = pos;
            Ok(())
        } else {
            let msg = format!(
                "Context Not Found, wante any of {:?}, got {:?}",
                self.keys, context
            );
            Err(crate::TimsqueryError::Other(msg))
        }
    }

    pub fn add(&mut self, rt_ms: u32, scan_index: usize, tof_index: u32, intensity: u64) {
        self.scan_tof_calc.add(
            self.context_key_num,
            rt_ms,
            scan_index,
            tof_index,
            intensity,
        );
    }

}

