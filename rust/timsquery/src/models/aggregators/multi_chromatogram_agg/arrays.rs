use super::aggregator::{
    DenseRTCollection,
    ParallelTracks,
    ParitionedCMGAggregator,
    SparseRTCollection,
};
use serde::Serialize;
use std::collections::HashMap;
use std::f64;
use std::hash::Hash;
use std::sync::Arc;
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};
use tracing::warn;
use crate::KeyLike;

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrayStats<FH: KeyLike> {
    pub retention_time_miliseconds: Arc<[u32]>,
    pub weighted_ims_mean: Vec<f64>,
    pub ims_means: HashMap<FH, Vec<f64>>,
    pub mz_means: HashMap<FH, Vec<f64>>,
    // TODO consider if I want to add the standard deviations ... RN they dont
    // seem to be that useful.
    pub intensities: HashMap<FH, Vec<u64>>,
}


impl <FH: KeyLike> PartitionedCMGArrayStats<FH> {

}

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: KeyLike> {
    pub transition_stats: Vec<ChromatomobilogramStatsArrays>,
    pub transition_keys: Vec<FH>,
    pub retention_times_ms: Arc<[u32]>,
    pub expected_scan_index: usize,
    pub expected_tof_indices: Vec<u32>,
}

impl<FH: KeyLike> PartitionedCMGArrays<FH> {
    pub fn new_with_sparse(
        collections: Vec<SparseRTCollection>,
        keys: Vec<FH>,
        expected_scan_index: usize,
        expected_tof_indices: Vec<u32>,
    ) -> Self {
        let mut transition_stats = Vec::with_capacity(keys.len());
        let mut uniq_rts: Vec<u32> = collections.iter().flatten().map(|x| *x.0).collect();
        uniq_rts.sort_unstable();
        uniq_rts.dedup();
        let uniq_rts: Arc<[u32]> = uniq_rts.into();

        for (id_ind, _id_key) in keys.iter().enumerate() {
            let mut id_cmgs = ChromatomobilogramStatsArrays::new();
            id_cmgs.retention_time_miliseconds = uniq_rts.clone();
            let local_id_mapping = &collections[id_ind];

            for rt_key in uniq_rts.iter() {
                let scan_tof_mapping = local_id_mapping.get(rt_key);
                if let Some(scan_tof_mapping) = scan_tof_mapping {
                    id_cmgs
                        .scan_index_means
                        .push(scan_tof_mapping.scan.mean().unwrap());
                    // id_cmgs
                    //     .scan_index_sds
                    //     .push(scan_tof_mapping.scan.standard_deviation().unwrap());
                    id_cmgs
                        .tof_index_means
                        .push(scan_tof_mapping.tof.mean().unwrap());
                    // id_cmgs
                    //     .tof_index_sds
                    //     .push(scan_tof_mapping.tof.standard_deviation().unwrap());
                    id_cmgs.intensities.push(scan_tof_mapping.tof.weight());
                }
            }
            transition_stats.push(id_cmgs);
        }

        Self {
            retention_times_ms: uniq_rts,
            transition_stats,
            transition_keys: keys,
            expected_tof_indices,
            expected_scan_index,
        }
    }

    pub fn new_with_dense(
        collections: Vec<DenseRTCollection>,
        keys: Vec<FH>,
        expected_scan_index: usize,
        expected_tof_indices: Vec<u32>,
    ) -> Self {
        let mut transition_stats = Vec::with_capacity(keys.len());
        let uniq_rts = collections.first().unwrap().reference_rt_ms.clone();
        // Q: Does this test pointer equality or values directly?
        // A: Pointer equality now that I changed eq to prt_eq
        assert!(
            collections
                .iter()
                .all(|x| Arc::ptr_eq(&x.reference_rt_ms, &uniq_rts)),
            "Expected all rts to come from the same Arc"
        );

        for (id_ind, _id_key) in keys.iter().enumerate() {
            let mut id_cmgs = ChromatomobilogramStatsArrays::empty_with_rts(
                uniq_rts.clone(),
                expected_tof_indices[id_ind],
                expected_scan_index,
            );
            let local_id_mapping = &collections[id_ind];

            local_id_mapping.scan_tof_calc.iter().enumerate().for_each(
                |(rt_key, scan_tof_mapping)| match scan_tof_mapping {
                    None => {}
                    Some(sts) => {
                        let scan_mean = sts.scan.mean().unwrap();
                        let tof_mean = sts.tof.mean().unwrap();
                        // let scan_sd = sts.scan.standard_deviation().unwrap();
                        // let tof_sd = sts.tof.standard_deviation().unwrap();
                        let inten = sts.tof.weight();

                        id_cmgs.intensities[rt_key] = inten;
                        id_cmgs.tof_index_means[rt_key] = tof_mean;
                        id_cmgs.scan_index_means[rt_key] = scan_mean;
                        // id_cmgs.scan_index_sds[rt_key] = scan_sd;
                        // id_cmgs.tof_index_sds[rt_key] = tof_sd;
                    }
                },
            );

            transition_stats.push(id_cmgs);
        }

        Self {
            retention_times_ms: uniq_rts,
            transition_stats,
            transition_keys: keys,
            expected_tof_indices,
            expected_scan_index,
        }
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    From<ParitionedCMGAggregator<FH>> for PartitionedCMGArrays<FH>
{
    fn from(item: ParitionedCMGAggregator<FH>) -> Self {
        // The main idea here is to "fill" the gaps in all XICs with 0s
        // If a transition was not observed at a specific retention time, it would not have
        // an entry if the storage is sparse.
        let keys = item.keys;
        match item.scan_tof_calc {
            ParallelTracks::Dense(scan_tof_calc) => Self::new_with_dense(
                scan_tof_calc,
                keys,
                item.expected_scan_index,
                item.expected_tof_indices,
            ),
            ParallelTracks::Sparse(scan_tof_calc) => Self::new_with_sparse(
                scan_tof_calc,
                keys,
                item.expected_scan_index,
                item.expected_tof_indices,
            ),
        }
    }
}

impl<FH: KeyLike> PartitionedCMGArrayStats<FH> {
    /// This step in essence converts the tof/scan indices to
    /// mz/ims units.
    pub fn new(
        other: PartitionedCMGArrays<FH>,
        mz_converter: &Tof2MzConverter,
        ims_converter: &Scan2ImConverter,
    ) -> Self {
        // Q: Why am I converting this to a hashmap and then re-converting it
        //    into a vec(array)
        // TODO: make a constructor to make sure everything here
        // Is actually getting added and is the right size/shape.
        let mut ims_means = HashMap::new();
        let mut mz_means = HashMap::new();
        let mut intensities_out = HashMap::new();

        let unique_rts = other.retention_times_ms.clone();

        // Products are used to calculate the weighted mean
        let mut tof_products: Vec<u64> = vec![0; unique_rts.len() + 1];
        let mut scan_products: Vec<u64> = vec![0; unique_rts.len() + 1];
        let mut summed_intensity_vec: Vec<u64> = vec![0; unique_rts.len() + 1];

        for (v, k) in other
            .transition_stats
            .into_iter()
            .zip(other.transition_keys.iter())
        {
            for i in 0..v.retention_time_miliseconds.len() {
                tof_products[i] += v.tof_index_means[i] as u64 * v.intensities[i];
                scan_products[i] += v.scan_index_means[i] as u64 * v.intensities[i];
                summed_intensity_vec[i] += v.intensities[i];
            }

            let intensities = v.intensities;
            let imss = v
                .scan_index_means
                .into_iter()
                .map(|x| ims_converter.convert(x))
                .collect();
            let mzs = v
                .tof_index_means
                .into_iter()
                .map(|x| mz_converter.convert(x))
                .collect();

            ims_means.insert(k.clone(), imss);
            mz_means.insert(k.clone(), mzs);
            intensities_out.insert(k.clone(), intensities);
        }

        let weighted_ims_mean = scan_products
            .into_iter()
            .zip(summed_intensity_vec.iter())
            .map(|(x, y)| {
                if *y < 10 {
                    return f64::NAN;
                }
                let out = (x / y) as f64;
                if !(0.0..=1100.0).contains(&out) {
                    warn!("Bad mobility value: {:?}, input was {:?}", out, x);
                }
                ims_converter.convert(out)
            })
            .collect();

        PartitionedCMGArrayStats {
            retention_time_miliseconds: unique_rts,
            ims_means,
            mz_means,
            weighted_ims_mean,
            intensities: intensities_out,
        }
    }

    pub fn len(&self) -> usize {
        self.retention_time_miliseconds.len()
    }

    // pub fn is_empty(&self) -> bool {
    //     self.retention_time_miliseconds.is_empty()
    // }
}
