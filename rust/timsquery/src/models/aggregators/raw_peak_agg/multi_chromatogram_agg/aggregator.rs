use super::super::chromatogram_agg::{
    ChromatomobilogramStatsArrays,
    ScanTofStatsCalculatorPair,
};
use crate::errors::Result;
use nohash_hasher::BuildNoHashHasher;
use serde::Serialize;
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Arc;

pub type SparseRTCollection = HashMap<u32, ScanTofStatsCalculatorPair, BuildNoHashHasher<u32>>;

// TODO: Consider using a typestate instead of an enum here...
pub enum RTCollection {
    Sparse(SparseRTCollection),
    Dense(DenseRTCollection),
}

#[derive(Debug, Clone)]
pub enum ParallelTracks {
    Sparse(Vec<SparseRTCollection>),
    Dense(Vec<DenseRTCollection>),
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

    fn add(&mut self, idx: usize, rt_ms: u32, scan_index: usize, tof_index: u32, intensity: u64) {
        match self {
            Self::Sparse(scan_tof_calc) => {
                scan_tof_calc[idx]
                    .entry(rt_ms)
                    .and_modify(|curr| {
                        curr.add(intensity, scan_index, tof_index);
                    })
                    .or_insert(ScanTofStatsCalculatorPair::new(
                        intensity, scan_index, tof_index,
                    ));
            }
            Self::Dense(scan_tof_calc) => {
                scan_tof_calc[idx].add(rt_ms, scan_index, tof_index, intensity);
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct DenseRTCollection {
    pub reference_rt_ms: Arc<[u32]>,
    pub scan_tof_calc: Vec<Option<ScanTofStatsCalculatorPair>>,
}

impl DenseRTCollection {
    fn new(reference_rt_ms: Arc<[u32]>) -> Self {
        // Check that its sorted ... let make it panic for now
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

    fn add(&mut self, rt_ms: u32, scan_index: usize, tof_index: u32, intensity: u64) {
        // Not the biggest fan of how branchy this code is ...
        let mut pos = self.reference_rt_ms.partition_point(|&x| x <= rt_ms);
        if pos == self.reference_rt_ms.len() {
            // If we are less than 1 second above the end, we can just subtract 1.
            let diff = rt_ms - self.reference_rt_ms[pos - 1];
            if diff < 1_000 {
                pos -= 1;
            }
        }
        // ~ This version is only ~5% faster (which I personally find surprising),
        // leaving here as a reference. I am keeping the exact solution for now.
        // let stride_ms = (self.reference_rt_ms.last().unwrap() - self.reference_rt_ms.first().unwrap()) / (self.scan_tof_calc.len() as u32 - 1);
        // let pos = (rt_ms - self.reference_rt_ms.first().unwrap()) / stride_ms;
        // let pos = usize::try_from(pos).unwrap().min(self.scan_tof_calc.len() - 1);

        match self.scan_tof_calc.get_mut(pos) {
            Some(x) => match x {
                Some(x) => {
                    x.add(intensity, scan_index, tof_index);
                }
                None => {
                    self.scan_tof_calc[pos] = Some(ScanTofStatsCalculatorPair::new(
                        intensity, scan_index, tof_index,
                    ))
                }
            },
            None => {
                let max_pos = self.scan_tof_calc.len() - 1;
                println!(
                    "DenseRTCollection::add out of bounds {:?} > {:?}",
                    pos, max_pos
                );
                panic!()
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct ParitionedCMGAggregator<
    FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug,
> {
    pub scan_tof_calc: ParallelTracks,
    pub keys: Vec<FH>,
    pub context_key_num: usize,
    pub expected_scan_index: usize,
    pub expected_tof_indices: Vec<u32>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    ParitionedCMGAggregator<FH>
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

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: HashMap<FH, ChromatomobilogramStatsArrays>,
}
