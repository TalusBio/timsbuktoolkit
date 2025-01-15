use super::super::chromatogram_agg::ChromatomobilogramStatsArrays;
use super::aggregator::ParitionedCMGAggregator;
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use serde::Serialize;
use std::collections::{
    BTreeMap,
    HashMap,
    HashSet,
};
use std::f64;
use std::hash::Hash;
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};
use tracing::warn;

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrayStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub retention_time_miliseconds: Vec<u32>,
    pub weighted_ims_mean: Vec<f64>,
    pub summed_intensity: Vec<u64>,
    pub ims_means: HashMap<FH, Vec<f64>>,
    pub mz_means: HashMap<FH, Vec<f64>>,
    // TODO consider if I want to add the standard deviations ... RN they dont
    // seem to be that useful.
    pub intensities: HashMap<FH, Vec<u64>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: Vec<ChromatomobilogramStatsArrays>,
    pub transition_keys: Vec<FH>,
    pub expected_scan_index: usize,
    pub expected_tof_indices: Vec<u32>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    From<ParitionedCMGAggregator<FH>> for PartitionedCMGArrays<FH>
{
    fn from(item: ParitionedCMGAggregator<FH>) -> Self {
        let mut transition_stats = Vec::with_capacity(item.keys.len());
        let mut uniq_rts: Vec<u32> = item.scan_tof_calc.iter().flatten().map(|x| *x.0).collect();
        uniq_rts.sort_unstable();
        uniq_rts.dedup();
        let mut out_expected_tof_indices: Vec<u32> = Vec::with_capacity(item.keys.len());

        for (id_ind, _id_key) in item.keys.iter().enumerate() {
            let mut id_cmgs = ChromatomobilogramStatsArrays::new();
            let local_id_mapping = &item.scan_tof_calc[id_ind];
            out_expected_tof_indices.push(item.expected_tof_indices[id_ind]);

            for rt_key in uniq_rts.iter() {
                let scan_tof_mapping = local_id_mapping.get(rt_key);
                if let Some(scan_tof_mapping) = scan_tof_mapping {
                    id_cmgs.retention_time_miliseconds.push(*rt_key);
                    id_cmgs
                        .scan_index_means
                        .push(scan_tof_mapping.scan.mean().unwrap());
                    id_cmgs
                        .scan_index_sds
                        .push(scan_tof_mapping.scan.standard_deviation().unwrap());
                    id_cmgs
                        .tof_index_means
                        .push(scan_tof_mapping.tof.mean().unwrap());
                    id_cmgs
                        .tof_index_sds
                        .push(scan_tof_mapping.tof.standard_deviation().unwrap());
                    id_cmgs.intensities.push(scan_tof_mapping.tof.weight());
                }
            }
            id_cmgs.sort_by_rt();
            transition_stats.push(id_cmgs);
        }

        Self {
            transition_stats,
            transition_keys: item.keys,
            expected_tof_indices: out_expected_tof_indices,
            expected_scan_index: item.expected_scan_index,
        }
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> PartitionedCMGArrayStats<FH> {
    pub fn new(
        other: PartitionedCMGArrays<FH>,
        mz_converter: &Tof2MzConverter,
        ims_converter: &Scan2ImConverter,
    ) -> Self {
        // TODO: make a constructor to make sure everything here
        // Is actually getting added and is the right size/shape.
        let mut ims_means = HashMap::new();
        let mut mz_means = HashMap::new();

        let unique_rts = other
            .transition_stats
            .iter()
            .flat_map(|v| v.retention_time_miliseconds.clone())
            // This is quite a bit of cloning ... I should re-write it to something more
            // performant.
            .collect::<HashSet<u32>>();
        let mut unique_rts = unique_rts.into_iter().collect::<Vec<u32>>();
        unique_rts.sort();

        let mut summed_intensity_vec: Vec<u64> = vec![0; unique_rts.len()];
        let mut intensity_logsums_vec: Vec<f64> = vec![0.0; unique_rts.len()];
        let mut npeaks_tree: Vec<u8> = vec![0; unique_rts.len()];
        let mut weighted_scan_index_mean_vec: Vec<Option<RunningStatsCalculator>> =
            vec![None; unique_rts.len()];

        let mut tmp_transition_intensities: Vec<Vec<u64>> = Vec::new();

        for (v, k) in other
            .transition_stats
            .into_iter()
            .zip(other.transition_keys.iter())
        {
            type ScanTofIntensityTuple = (f64, f64, u64);
            type ScanTofIntensityVecs = (Vec<f64>, (Vec<f64>, Vec<u64>));
            // TODO: replace with a vec ...
            let mut tmp_tree: BTreeMap<u32, ScanTofIntensityTuple> = BTreeMap::new();

            for i in 0..v.retention_time_miliseconds.len() {
                let rt = v.retention_time_miliseconds[i];
                let final_index = unique_rts
                    .binary_search(&rt)
                    .expect("All RTs should be registered above.");

                let scan = v.scan_index_means[i];
                let tof = v.tof_index_means[i];
                let inten = v.intensities[i];

                if inten > 100 {
                    npeaks_tree[final_index] += 1;
                }
                tmp_tree.entry(rt).or_insert((scan, tof, inten));
                summed_intensity_vec[final_index] += inten;
                if inten > 10 {
                    intensity_logsums_vec[final_index] += (inten as f64).ln();
                }

                if let Some(mut x) = weighted_scan_index_mean_vec[final_index] {
                    x.add(scan, inten);
                } else {
                    weighted_scan_index_mean_vec[final_index] =
                        Some(RunningStatsCalculator::new(inten, scan));
                }
            }

            // Now we fill with nans the missing values.
            // Q: Do I really need to do this here?
            for rt in unique_rts.iter() {
                tmp_tree.entry(*rt).or_insert((f64::NAN, f64::NAN, 0));
            }

            let (out_scans, (out_tofs, out_intens)): ScanTofIntensityVecs = tmp_tree
                .into_iter()
                .map(|(_, (scan, tof, inten))| {
                    (
                        ims_converter.convert(scan),
                        (mz_converter.convert(tof), inten),
                    )
                })
                .unzip();
            ims_means.insert(k.clone(), out_scans);
            mz_means.insert(k.clone(), out_tofs);
            tmp_transition_intensities.push(out_intens);
        }

        let retention_time_miliseconds = unique_rts;
        let summed_intensity = summed_intensity_vec;

        let weighted_ims_mean = weighted_scan_index_mean_vec
            .into_iter()
            .map(|x| {
                let out = x
                    .expect("At least one value should be given for every RT")
                    .mean()
                    .expect("At least one value should be present");
                if !(0.0..=1000.0).contains(&out) {
                    warn!("Bad mobility value: {:?}, input was {:?}", out, x);
                }

                ims_converter.convert(out)
            })
            .collect();
        let intensities = other
            .transition_keys
            .into_iter()
            .zip(tmp_transition_intensities)
            .collect();

        PartitionedCMGArrayStats {
            retention_time_miliseconds,
            ims_means,
            mz_means,
            weighted_ims_mean,
            intensities,
            summed_intensity,
        }
    }
}
