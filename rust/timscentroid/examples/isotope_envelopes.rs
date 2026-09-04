//! Centroiding quality through isotope envelopes.
//!
//! Seeds are confident precursors from a search (`mz,charge,rt_s,mobility` CSV,
//! header row, extra columns ignored). For each seed the MS1 index is queried
//! at M+0 .. M+4 (spacing 1.00335 / z) within +-10 ppm, +-2 % mobility and
//! +-2 cycles of the apex. Two things are measured per isotope:
//!
//! - visible: any centroid in that box
//! - centroids at the apex cycle: how many centroids the isotope became on the
//!   apex frame. One is right. More means the centroider split the ion, along
//!   TOF (quantization jitter) or along scans.
//!
//! The same is measured at M+0.5, M+1.5, ... as a null: what the box finds
//! where no isotope should be.
//!
//! Usage:
//!   isotope_envelopes <file.d> <seeds.csv> [--ms1-mz-bins N] [--ms2-mz-bins N]
//!                     [--ms1-cap N/W] [--ms2-cap N/W] [--seeds N]

use half::f16;
use timscentroid::centroiding::{
    MzTolerance,
    WindowCap,
};
use timscentroid::rt_mapping::{
    MS1CycleIndex,
    RTIndex,
};
use timscentroid::utils::OptionallyRestricted::Restricted;
use timscentroid::utils::TupleRange;
use timscentroid::{
    IndexedTimstofPeaks,
    IndexingCentroidingConfig,
};
use timsrust::TimsTofPath;

const ISOTOPE_DA: f64 = 1.003_354_83;
const N_ISOTOPES: usize = 5;
const PPM: f64 = 10.0;
const MOB_PCT: f32 = 2.0;
const CYCLE_HALO: u32 = 2;

struct Seed {
    mz: f64,
    charge: u8,
    rt_s: f64,
    mobility: f32,
}

#[derive(Default, Clone, Copy)]
struct Box_ {
    visible: bool,
    centroids_at_apex: u32,
    intensity: f64,
}

fn parse_cap(s: &str) -> WindowCap {
    let (n, w) = s.split_once('/').expect("cap as N/W");
    WindowCap {
        max_peaks: n.parse().unwrap(),
        window_da: w.parse().unwrap(),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: isotope_envelopes <file.d> <seeds.csv> [flags]");
    let seeds_path = args.get(2).expect("seeds csv");
    let mut cfg = IndexingCentroidingConfig::default();
    let mut max_seeds = usize::MAX;
    let mut i = 3;
    while i < args.len() {
        let v = args.get(i + 1).expect("flag value");
        match args[i].as_str() {
            "--ms1-mz-bins" => cfg.ms1.mz_tol = MzTolerance::Bins(v.parse().unwrap()),
            "--ms2-mz-bins" => cfg.ms2.mz_tol = MzTolerance::Bins(v.parse().unwrap()),
            "--ms1-cap" => cfg.ms1.window_cap = Some(parse_cap(v)),
            "--ms2-cap" => cfg.ms2.window_cap = Some(parse_cap(v)),
            "--seeds" => max_seeds = v.parse().unwrap(),
            f => panic!("unknown flag {f}"),
        }
        i += 2;
    }

    let seeds: Vec<Seed> = std::fs::read_to_string(seeds_path)
        .unwrap()
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split(',').collect();
            Seed {
                mz: c[0].parse().unwrap(),
                charge: c[1].parse().unwrap(),
                rt_s: c[2].parse().unwrap(),
                mobility: c[3].parse().unwrap(),
            }
        })
        .take(max_seeds)
        .collect();

    println!("file: {path}");
    println!("seeds: {} from {seeds_path}", seeds.len());
    println!("ms1: {:?}", cfg.ms1);
    println!("ms2: {:?}", cfg.ms2);

    let file = TimsTofPath::new(path).unwrap();
    let st = std::time::Instant::now();
    let (index, stats) = IndexedTimstofPeaks::from_timstof_file(&file, cfg);
    println!(
        "index built in {:.1?}: ms1 {} peaks, ms2 {} peaks",
        st.elapsed(),
        stats.ms1_stats.indexing_stats.num_peaks,
        stats.ms2_stats.indexing_stats.num_peaks
    );

    let mut on: Vec<[Box_; N_ISOTOPES]> = Vec::with_capacity(seeds.len());
    let mut off: Vec<[Box_; N_ISOTOPES]> = Vec::with_capacity(seeds.len());
    for s in &seeds {
        let apex = index.rt_ms_to_cycle_index((s.rt_s * 1000.0) as u32);
        let c = apex.as_u32();
        let cycles = TupleRange::try_new(
            MS1CycleIndex::new(c.saturating_sub(CYCLE_HALO)),
            MS1CycleIndex::new(c + CYCLE_HALO),
        )
        .unwrap();
        let mob = TupleRange::try_new(
            f16::from_f32(s.mobility * (1.0 - MOB_PCT / 100.0)),
            f16::from_f32(s.mobility * (1.0 + MOB_PCT / 100.0)),
        )
        .unwrap();
        let measure = |offset_isotopes: f64| -> Box_ {
            let mz = s.mz + offset_isotopes * ISOTOPE_DA / s.charge as f64;
            let mzr = TupleRange::try_new(
                (mz * (1.0 - PPM / 1e6)) as f32,
                (mz * (1.0 + PPM / 1e6)) as f32,
            )
            .unwrap();
            let mut b = Box_::default();
            for p in index.query_peaks_ms1(mzr, Restricted(cycles), Restricted(mob)) {
                b.visible = true;
                b.intensity += p.intensity as f64;
                if p.cycle_index == apex {
                    b.centroids_at_apex += 1;
                }
            }
            b
        };
        let mut row_on = [Box_::default(); N_ISOTOPES];
        let mut row_off = [Box_::default(); N_ISOTOPES];
        for k in 0..N_ISOTOPES {
            row_on[k] = measure(k as f64);
            row_off[k] = measure(k as f64 + 0.5);
        }
        on.push(row_on);
        off.push(row_off);
    }

    report("isotope positions M+0..M+4", &on);
    report("null positions M+0.5..M+4.5", &off);

    // Envelope depth: consecutive visible isotopes from M+0.
    let mut depth_hist = [0usize; N_ISOTOPES + 1];
    for row in &on {
        let d = row.iter().take_while(|b| b.visible).count();
        depth_hist[d] += 1;
    }
    println!("\nenvelope depth (consecutive isotopes visible from M+0):");
    for (d, n) in depth_hist.iter().enumerate() {
        println!(
            "  {d} isotopes: {:>5.1} %",
            100.0 * *n as f64 / on.len().max(1) as f64
        );
    }

    // Averagine: expected envelope for a peptide of this mass, compared to the
    // observed one by cosine. 1.0 is a textbook envelope; a split or merged
    // isotope, or a co-isolated ion, pulls it down. Seeds with nothing at M+0
    // are counted separately rather than scored.
    let mut cos: Vec<f64> = Vec::new();
    let mut no_m0 = 0usize;
    for (s, row) in seeds.iter().zip(&on) {
        if !row[0].visible {
            no_m0 += 1;
            continue;
        }
        let mass = (s.mz - PROTON) * s.charge as f64;
        let expected = averagine_envelope(mass);
        let observed: [f64; N_ISOTOPES] = std::array::from_fn(|k| row[k].intensity);
        cos.push(cosine(&expected, &observed));
    }
    cos.sort_by(|a, b| a.total_cmp(b));
    if !cos.is_empty() {
        let frac =
            |t: f64| 100.0 * cos.iter().filter(|&&c| c >= t).count() as f64 / cos.len() as f64;
        println!(
            "\naveragine cosine (M+0..M+4, n={}, {} seeds had no M+0):",
            cos.len(),
            no_m0
        );
        println!(
            "  p10 {:.3}  p50 {:.3}  p90 {:.3}   >=0.90: {:.1} %   >=0.95: {:.1} %   >=0.98: {:.1} %",
            cos[cos.len() / 10],
            cos[cos.len() / 2],
            cos[cos.len() * 9 / 10],
            frac(0.90),
            frac(0.95),
            frac(0.98)
        );
    }

    // M+1 / M+0 ratio for what is visible on both; a sanity check that the
    // envelope is an envelope and not two unrelated ions.
    let mut ratios: Vec<f64> = on
        .iter()
        .filter(|r| r[0].visible && r[1].visible && r[0].intensity > 0.0)
        .map(|r| r[1].intensity / r[0].intensity)
        .collect();
    ratios.sort_by(|a, b| a.total_cmp(b));
    if !ratios.is_empty() {
        println!(
            "M+1 / M+0 intensity ratio: p10 {:.2}  p50 {:.2}  p90 {:.2}  (n={})",
            ratios[ratios.len() / 10],
            ratios[ratios.len() / 2],
            ratios[ratios.len() * 9 / 10],
            ratios.len()
        );
    }
}

const PROTON: f64 = 1.007_276_47;

/// Expected relative intensities of the first `N_ISOTOPES` aggregated
/// isotopologues for an averagine peptide of the given neutral mass.
///
/// Averagine: C 4.9384, H 7.7583, N 1.3577, O 1.4773, S 0.0417 per 111.1254 Da
/// (Senko 1995). Element counts are rounded, then each element's per-atom
/// isotope distribution is convolved in `count` times, truncated at M+4.
fn averagine_envelope(mass: f64) -> [f64; N_ISOTOPES] {
    let units = mass / 111.1254;
    // (per-atom offset probabilities, atoms)
    let elements: [(&[f64], usize); 5] = [
        (&[0.9893, 0.0107], (4.9384 * units).round() as usize),
        (&[0.999_885, 0.000_115], (7.7583 * units).round() as usize),
        (&[0.996_36, 0.003_64], (1.3577 * units).round() as usize),
        (
            &[0.997_57, 0.000_38, 0.002_05],
            (1.4773 * units).round() as usize,
        ),
        (
            &[0.9499, 0.0075, 0.0425, 0.0, 0.0001],
            (0.0417 * units).round() as usize,
        ),
    ];
    let mut dist = [0.0; N_ISOTOPES];
    dist[0] = 1.0;
    for (atom, count) in elements {
        for _ in 0..count {
            let mut next = [0.0; N_ISOTOPES];
            for (i, &p) in dist.iter().enumerate() {
                if p == 0.0 {
                    continue;
                }
                for (j, &q) in atom.iter().enumerate() {
                    if i + j < N_ISOTOPES {
                        next[i + j] += p * q;
                    }
                }
            }
            dist = next;
        }
    }
    dist
}

fn cosine(a: &[f64; N_ISOTOPES], b: &[f64; N_ISOTOPES]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na == 0.0 || nb == 0.0 {
        0.0
    } else {
        dot / (na * nb)
    }
}

fn report(name: &str, rows: &[[Box_; N_ISOTOPES]]) {
    println!("\n== {name} ==");
    println!(
        "{:>8} {:>10} {:>22} {:>16}",
        "isotope", "visible %", "centroids@apex (mean)", "split (>1) %"
    );
    for k in 0..N_ISOTOPES {
        let vis: Vec<&Box_> = rows.iter().map(|r| &r[k]).filter(|b| b.visible).collect();
        let n = rows.len().max(1) as f64;
        let apex_hits: Vec<u32> = vis
            .iter()
            .map(|b| b.centroids_at_apex)
            .filter(|&c| c > 0)
            .collect();
        let mean = apex_hits.iter().sum::<u32>() as f64 / apex_hits.len().max(1) as f64;
        let split =
            apex_hits.iter().filter(|&&c| c > 1).count() as f64 / apex_hits.len().max(1) as f64;
        println!(
            "{:>8} {:>9.1} % {:>22.2} {:>15.1} %",
            format!("M+{k}"),
            100.0 * vis.len() as f64 / n,
            mean,
            100.0 * split
        );
    }
}
