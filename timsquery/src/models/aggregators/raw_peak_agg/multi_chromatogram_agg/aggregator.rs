use super::super::chromatogram_agg::{
    ChromatomobilogramStatsArrays,
    ScanTofStatsCalculatorPair,
};
use crate::errors::Result;
use nohash_hasher::BuildNoHashHasher;
use serde::Serialize;
use std::collections::HashMap;
use std::hash::Hash;

type SparseRTCollection = HashMap<u32, ScanTofStatsCalculatorPair, BuildNoHashHasher<u32>>;

#[derive(Debug, Clone)]
pub struct ParitionedCMGAggregator<
    FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug,
> {
    pub scan_tof_calc: Vec<SparseRTCollection>,
    pub keys: Vec<FH>,
    pub context_key_num: usize,
    pub expected_scan_index: usize,
    pub expected_tof_indices: Vec<u32>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    ParitionedCMGAggregator<FH>
{
    pub fn new(keys: Vec<FH>, expected_scan_index: usize, expected_tof_indices: Vec<u32>) -> Self {
        let mut scan_tof_calc = Vec::with_capacity(keys.len());
        for _ in 0..keys.len() {
            scan_tof_calc.push(SparseRTCollection::with_hasher(BuildNoHashHasher::default()));
        }
        assert!(expected_tof_indices.len() == keys.len());

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
        self.scan_tof_calc[self.context_key_num]
            .entry(rt_ms)
            .and_modify(|curr| {
                curr.add(intensity, scan_index, tof_index);
            })
            .or_insert(ScanTofStatsCalculatorPair::new(
                intensity, scan_index, tof_index,
            ));
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: HashMap<FH, ChromatomobilogramStatsArrays>,
}
