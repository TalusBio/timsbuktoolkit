//! Ratio-free temporal-coelution weighting.
//!
//! For each cycle, over a +/-W window, each active fragment's XIC slice is
//! centered and normalized to a unit vector. The coelution score is the
//! expected-intensity-weighted mean pairwise correlation of those unit
//! vectors. Because slices are centered+normalized, the score is independent
//! of absolute fragment ratios: it rewards genuine co-eluting peaks and
//! suppresses uncorrelated random-interferent pileups.
//! `profile[t] *= 1 + K * max(coel, 0)`. Disabled when K<=0.

use super::{
    ApexConfig,
    ApexScratch,
    Rows,
};

/// Multiply `profile` in place by the coelution weight. `rows` are the active
/// fragment rows (see [`super::active_rows`]).
pub(crate) fn weight(
    profile: &mut [f32],
    rows: &Rows<'_>,
    cfg: &ApexConfig,
    scratch: &mut ApexScratch,
) {
    let k = cfg.coel_k;
    if k <= 0.0 || profile.is_empty() || rows.len() < 2 {
        return;
    }
    let w = cfg.coel_w.max(1.0) as usize;
    kernel(profile, rows, w, k, scratch);
}

/// Pure numeric core, decoupled from the extraction data model for testing.
fn kernel(profile: &mut [f32], rows: &Rows<'_>, w: usize, k: f32, scratch: &mut ApexScratch) {
    let n = profile.len();
    let win = 2 * w + 1;
    let ApexScratch {
        cn, active, wacc, ..
    } = scratch;
    cn.clear();
    cn.resize(rows.len() * win, 0.0);
    wacc.clear();
    wacc.resize(win, 0.0);

    for (t, p) in profile.iter_mut().enumerate() {
        // The weight only ever multiplies profile[t]; a zero cycle stays zero
        // for any coel, so skip its O(rows*win) window work.
        if *p == 0.0 {
            continue;
        }

        let lo = t.saturating_sub(w);
        let hi = (t + w + 1).min(n);
        let wl = hi - lo;
        active.clear();

        for (fi, (chrom, _wt)) in rows.iter().enumerate() {
            let slice = &chrom[lo..hi];
            let mean: f32 = slice.iter().sum::<f32>() / wl as f32;
            let base = fi * win;
            let mut nsq = 0.0f32;
            for j in 0..wl {
                let v = slice[j] - mean;
                cn[base + j] = v;
                nsq += v * v;
            }
            let norm = nsq.sqrt();
            if norm > 1e-12 {
                for j in 0..wl {
                    cn[base + j] /= norm;
                }
                active.push(fi);
            }
        }

        // Expected-intensity-weighted mean pairwise correlation, computed in
        // O(active) rather than O(active^2): for unit vectors u_a,
        //   sum_{a<b} w_a w_b <u_a,u_b> = (||W||^2 - sum w_a^2) / 2,
        //   sum_{a<b} w_a w_b          = ((sum w_a)^2 - sum w_a^2) / 2,
        // where W = sum_a w_a u_a. The /2 cancels in the ratio.
        wacc[..wl].fill(0.0);
        let mut sw = 0.0f32;
        let mut sw2 = 0.0f32;
        for &fi in active.iter() {
            let base = fi * win;
            let wa = rows[fi].1;
            sw += wa;
            sw2 += wa * wa;
            for j in 0..wl {
                wacc[j] += wa * cn[base + j];
            }
        }
        let mut wnorm2 = 0.0f32;
        for &x in &wacc[..wl] {
            wnorm2 += x * x;
        }
        let den = sw * sw - sw2;
        let coel = if den > 0.0 {
            ((wnorm2 - sw2) / den).max(0.0)
        } else {
            0.0
        };
        *p *= 1.0 + k * coel;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Gaussian bump centered at `c` over `n` cycles.
    fn bump(n: usize, c: usize) -> Vec<f32> {
        (0..n)
            .map(|t| {
                let d = t as f32 - c as f32;
                (-(d * d) / 2.0).exp()
            })
            .collect()
    }

    #[test]
    fn coelution_rewards_coeluting_over_scattered() {
        let n = 11;
        let cfg = ApexConfig {
            coel_k: 1.0,
            coel_w: 2.0,
            ..Default::default()
        };
        let mut scratch = ApexScratch::default();

        // Co-eluting: three fragments all peaking at cycle 5.
        let (a, b, c) = (bump(n, 5), bump(n, 5), bump(n, 5));
        let rows_coel: Rows = vec![
            (a.as_slice(), 1.0),
            (b.as_slice(), 1.0),
            (c.as_slice(), 1.0),
        ];
        let mut p_coel = vec![1.0f32; n];
        weight(&mut p_coel, &rows_coel, &cfg, &mut scratch);

        // Scattered: three fragments peaking at different cycles.
        let (d, e, f) = (bump(n, 3), bump(n, 5), bump(n, 7));
        let rows_scat: Rows = vec![
            (d.as_slice(), 1.0),
            (e.as_slice(), 1.0),
            (f.as_slice(), 1.0),
        ];
        let mut p_scat = vec![1.0f32; n];
        weight(&mut p_scat, &rows_scat, &cfg, &mut scratch);

        // The apex of genuinely co-eluting fragments is boosted...
        assert!(
            p_coel[5] > 1.0,
            "co-eluting apex must be up-weighted: {}",
            p_coel[5]
        );
        // ...and more than a scattered pileup at the same cycle.
        assert!(
            p_coel[5] > p_scat[5],
            "co-eluting {} should beat scattered {}",
            p_coel[5],
            p_scat[5]
        );
    }

    #[test]
    fn coelution_noop_when_disabled_or_too_few_rows() {
        let n = 7;
        let mut scratch = ApexScratch::default();
        let a = bump(n, 3);
        let b = bump(n, 3);

        // Disabled (k = 0).
        let cfg_off = ApexConfig {
            coel_k: 0.0,
            ..Default::default()
        };
        let rows: Rows = vec![(a.as_slice(), 1.0), (b.as_slice(), 1.0)];
        let mut p = vec![1.0f32; n];
        weight(&mut p, &rows, &cfg_off, &mut scratch);
        assert!(p.iter().all(|&x| x == 1.0));

        // Fewer than two rows: nothing to correlate.
        let cfg_on = ApexConfig {
            coel_k: 1.0,
            ..Default::default()
        };
        let rows_one: Rows = vec![(a.as_slice(), 1.0)];
        let mut p2 = vec![1.0f32; n];
        weight(&mut p2, &rows_one, &cfg_on, &mut scratch);
        assert!(p2.iter().all(|&x| x == 1.0));
    }
}
