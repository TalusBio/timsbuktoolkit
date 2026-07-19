//! Cross-row coincidence-voting weighting.
//!
//! Per fragment row: subtract the row's global median (baseline), clip to >=0,
//! matched-filter with the known gaussian (sigma~1), then convert to a robust
//! soft vote `v = sigmoid((z - tau)/s)` where `z = (m - median(m))/MAD`. Votes
//! are summed across rows (expected-intensity weighted) into a support trace;
//! `profile[t] *= 1 + K * support_norm(t)`. Because interferents are
//! independent across rows, only a genuine apex gets many concurrent votes; a
//! vote saturates so a 10x-tall interferent counts the same as a real peak.
//! Disabled when K<=0.

use super::{
    ApexConfig,
    ApexScratch,
    Rows,
};

/// Gaussian matched-filter kernel (sigma~1), 5-tap.
const G: [f32; 5] = [0.13534, 0.60653, 1.0, 0.60653, 0.13534];

/// Multiply `profile` in place by the coincidence-vote weight. `rows` are the
/// active fragment rows (see [`super::active_rows`]).
pub(crate) fn weight(
    profile: &mut [f32],
    rows: &Rows<'_>,
    cfg: &ApexConfig,
    scratch: &mut ApexScratch,
) {
    let vk = cfg.vote_k;
    if vk <= 0.0 || profile.is_empty() || rows.is_empty() {
        return;
    }
    kernel(
        profile,
        rows,
        cfg.vote_tau,
        cfg.vote_s.max(1e-3),
        vk,
        scratch,
    );
}

/// Pure numeric core, decoupled from the extraction data model for testing.
fn kernel(
    profile: &mut [f32],
    rows: &Rows<'_>,
    tau: f32,
    soft: f32,
    vk: f32,
    scratch: &mut ApexScratch,
) {
    let n = profile.len();
    let gsum: f32 = G.iter().sum();

    // `support` must start zeroed (accumulated across rows); r/m/tmp are fully
    // overwritten per row.
    let ApexScratch {
        support, r, m, tmp, ..
    } = scratch;
    support.clear();
    support.resize(n, 0.0);
    r.resize(n, 0.0);
    m.resize(n, 0.0);
    tmp.clear();
    let mut wsum = 0.0f32;

    for &(chrom, w) in rows.iter() {
        // Baseline = global median of the row. select_nth is O(n) vs
        // O(n log n) for a full sort; only the k-th order statistic is needed.
        tmp.clear();
        tmp.extend_from_slice(chrom);
        let med = *tmp.select_nth_unstable_by(n / 2, f32::total_cmp).1;
        for t in 0..n {
            r[t] = (chrom[t] - med).max(0.0);
        }

        // Matched filter (edge-clamped convolution with G).
        for t in 0..n {
            let mut acc = 0.0f32;
            for (di, &g) in G.iter().enumerate() {
                let idx = t as isize + di as isize - 2;
                if idx >= 0 && (idx as usize) < n {
                    acc += g * r[idx as usize];
                }
            }
            m[t] = acc / gsum;
        }

        // Robust scale: MAD of m (both medians via O(n) selection).
        tmp.clear();
        tmp.extend_from_slice(m.as_slice());
        let mmed = *tmp.select_nth_unstable_by(n / 2, f32::total_cmp).1;
        for t in 0..n {
            tmp[t] = (m[t] - mmed).abs();
        }
        let mad = (*tmp.select_nth_unstable_by(n / 2, f32::total_cmp).1).max(1e-6);

        // Soft vote per cycle.
        for t in 0..n {
            let z = (m[t] - mmed) / mad;
            let v = 1.0 / (1.0 + (-(z - tau) / soft).exp());
            support[t] += w * v;
        }
        wsum += w;
    }

    if wsum <= 0.0 {
        return;
    }
    for t in 0..n {
        profile[t] *= 1.0 + vk * support[t] / wsum;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bump(n: usize, c: usize, height: f32) -> Vec<f32> {
        (0..n)
            .map(|t| {
                let d = t as f32 - c as f32;
                height * (-(d * d) / 2.0).exp()
            })
            .collect()
    }

    #[test]
    fn many_concurrent_votes_beat_one_tall_interferent() {
        let n = 11;
        let cfg = ApexConfig {
            vote_k: 10.0,
            vote_tau: 1.28,
            vote_s: 0.42,
            ..Default::default()
        };
        let mut scratch = ApexScratch::default();

        // Five fragments each with a modest coincident peak at cycle 5.
        let cols: Vec<Vec<f32>> = (0..5).map(|_| bump(n, 5, 1.0)).collect();
        let rows_many: Rows = cols.iter().map(|c| (c.as_slice(), 1.0)).collect();
        let mut p_many = vec![1.0f32; n];
        weight(&mut p_many, &rows_many, &cfg, &mut scratch);

        // One single 10x-tall interferent peak, four flat rows.
        let tall = bump(n, 5, 10.0);
        let flat = vec![0.0f32; n];
        let rows_one: Rows = vec![
            (tall.as_slice(), 1.0),
            (flat.as_slice(), 1.0),
            (flat.as_slice(), 1.0),
            (flat.as_slice(), 1.0),
            (flat.as_slice(), 1.0),
        ];
        let mut p_one = vec![1.0f32; n];
        weight(&mut p_one, &rows_one, &cfg, &mut scratch);

        // Concurrent multi-fragment support wins over a lone tall interferent.
        assert!(
            p_many[5] > 1.0,
            "coincident apex must be up-weighted: {}",
            p_many[5]
        );
        assert!(
            p_many[5] > p_one[5],
            "5 concurrent votes {} should beat 1 tall interferent {}",
            p_many[5],
            p_one[5]
        );
    }

    #[test]
    fn vote_noop_when_disabled() {
        let n = 7;
        let cfg = ApexConfig {
            vote_k: 0.0,
            ..Default::default()
        };
        let mut scratch = ApexScratch::default();
        let c = bump(n, 3, 1.0);
        let rows: Rows = vec![(c.as_slice(), 1.0)];
        let mut p = vec![1.0f32; n];
        weight(&mut p, &rows, &cfg, &mut scratch);
        assert!(p.iter().all(|&x| x == 1.0));
    }
}
