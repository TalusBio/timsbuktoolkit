//! Centroiding quality through isotope envelopes, on a handful of frames.
//!
//! No index is built. For each seed precursor (`mz,charge,rt_s,mobility,
//! sequence` CSV from a search) the MS1 frames around its apex and the MS2
//! frames of the window group that isolates it are centroided with the given
//! config, and envelopes are read straight off those centroids:
//!
//! - MS1: M+0..M+4 of the precursor at +-10 ppm, +-2 % mobility.
//! - MS2: M+0..M+2 of every b/y fragment (z = 1, 2) computed from the
//!   sequence, at +-15 ppm, +-2 % of the precursor mobility. Sequences with
//!   modifications are skipped for fragments.
//!
//! Per isotope: visible (any centroid in the box), centroids on the apex frame
//! (one is right; more is a split ion), and an averagine cosine over the
//! envelope. The same boxes shifted by half an isotope give the null rate.
//!
//! Frames are centroided once and shared across seeds, so a config change on
//! 400 seeds costs seconds.
//!
//! Usage:
//!   isotope_frames <file.d> <seeds.csv> [--seeds N] [--halo F]
//!       [--rt-center S --rt-window S]   (seeds with apex in center +- window)
//!       [--box-im-pct X] [--box-ms1-ppm X] [--box-ms2-ppm X]   (readout box)
//!       [--ms1-mz-bins N] [--ms2-mz-bins N] [--ms1-ppm X] [--ms2-ppm X]
//!       [--ms1-im-pct X] [--ms2-im-pct X] [--ms1-im-scans N] [--ms2-im-scans N]
//!       [--ms1-cap N/W] [--ms2-cap N/W] [--ms1-early-stop N] [--ms2-early-stop N]
//!       [--ms1-max-peaks N] [--ms2-max-peaks N] [--transitive true|false]
//!
//! The readout box must be wider than anything the config under test can move
//! a centroid by, or the metric measures the box, not the centroiding. With
//! the default 2 % mobility box, configs that merge across 2-3 % of mobility
//! score lower purely because merged centroids land outside it; at 4 % they
//! do not. Compare at two widths before trusting a ranking.
//!
//! "No centroiding" baseline: `--ms1-im-scans 0 --ms2-im-scans 0
//! --ms1-max-peaks 10000000 --ms2-max-peaks 10000000 --ms1-early-stop 0`
//! with the default Bins(0): every raw entry is its own centroid.

use std::collections::HashMap;

use rayon::prelude::*;
use timscentroid::centroiding::{
    MzTolerance,
    PeakCentroider,
    WindowCap,
};
use timscentroid::{
    CentroidingConfig,
    IndexingCentroidingConfig,
};
use timsrust::converters::{
    ConvertableDomain,
    Tof2MzConverter2,
};
use timsrust::readers::TdfBlob;
use timsrust::{
    FramePeaks,
    MSLevel,
    TimsTofPath,
};

const ISOTOPE_DA: f64 = 1.003_354_83;
const PROTON: f64 = 1.007_276_47;
const WATER: f64 = 18.010_565;
const MS1_PPM: f64 = 10.0;
const MS2_PPM: f64 = 15.0;
const MOB_PCT: f64 = 2.0;
const N_MS1: usize = 5;
const N_MS2: usize = 3;

struct Seed {
    mz: f64,
    charge: u8,
    rt_s: f64,
    mobility: f64,
    sequence: Option<String>,
}

#[derive(Clone, Copy)]
struct Cent {
    mz: f64,
    im: f64,
    intensity: f64,
}

#[derive(Default, Clone, Copy)]
struct Box_ {
    visible: bool,
    at_apex: u32,
    intensity: f64,
    /// Spread among the apex-frame centroids in the box, when there is more
    /// than one: which axis the split is on.
    apex_mz_spread_ppm: f64,
    apex_im_spread_pct: f64,
    /// Mobility of the brightest centroid in the box, relative to the seed's
    /// mobility, in percent. Across a peptide's fragments and isotopes this is
    /// the within-peptide mobility spread: the quantity a mobility window has
    /// to cover once the query is centered on the observed mobility.
    best_im_offset_pct: f64,
}

fn residue_mass(c: u8) -> Option<f64> {
    Some(match c {
        b'G' => 57.021_46,
        b'A' => 71.037_11,
        b'S' => 87.032_03,
        b'P' => 97.052_76,
        b'V' => 99.068_41,
        b'T' => 101.047_68,
        b'C' => 103.009_19,
        b'L' | b'I' => 113.084_06,
        b'N' => 114.042_93,
        b'D' => 115.026_94,
        b'Q' => 128.058_58,
        b'K' => 128.094_96,
        b'E' => 129.042_59,
        b'M' => 131.040_49,
        b'H' => 137.058_91,
        b'F' => 147.068_41,
        b'R' => 156.101_11,
        b'Y' => 163.063_33,
        b'W' => 186.079_31,
        _ => return None,
    })
}

/// (m/z, charge, neutral mass) of b and y ions, z = 1 and 2, within 100..1700.
fn fragments(seq: &str) -> Vec<(f64, u8, f64)> {
    let masses: Option<Vec<f64>> = seq.bytes().map(residue_mass).collect();
    let Some(masses) = masses else {
        return Vec::new();
    };
    let n = masses.len();
    let mut out = Vec::new();
    for i in 1..n {
        let b: f64 = masses[..i].iter().sum();
        let y: f64 = masses[n - i..].iter().sum::<f64>() + WATER;
        for neutral in [b, y] {
            for z in 1u8..=2 {
                let mz = (neutral + z as f64 * PROTON) / z as f64;
                if (100.0..=1700.0).contains(&mz) {
                    out.push((mz, z, neutral));
                }
            }
        }
    }
    out
}

fn averagine<const N: usize>(mass: f64) -> [f64; N] {
    let units = mass / 111.1254;
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
    let mut dist = [0.0; N];
    dist[0] = 1.0;
    for (atom, count) in elements {
        for _ in 0..count {
            let mut next = [0.0; N];
            for (i, &p) in dist.iter().enumerate() {
                if p == 0.0 {
                    continue;
                }
                for (j, &q) in atom.iter().enumerate() {
                    if i + j < N {
                        next[i + j] += p * q;
                    }
                }
            }
            dist = next;
        }
    }
    dist
}

fn cosine<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na == 0.0 || nb == 0.0 {
        0.0
    } else {
        dot / (na * nb)
    }
}

fn parse_cap(s: &str) -> WindowCap {
    let (n, w) = s.split_once('/').expect("cap as N/W");
    WindowCap {
        max_peaks: n.parse().unwrap(),
        window_da: w.parse().unwrap(),
    }
}

/// Centroids of `frames` in a box, per frame: (visible, hits on `apex`, sum).
fn measure(
    cents: &HashMap<usize, Vec<Cent>>,
    frames: &[usize],
    apex: usize,
    mz: f64,
    ppm: f64,
    im: f64,
    im_pct: f64,
) -> Box_ {
    let lo = mz * (1.0 - ppm / 1e6);
    let hi = mz * (1.0 + ppm / 1e6);
    let im_lo = im * (1.0 - im_pct / 100.0);
    let im_hi = im * (1.0 + im_pct / 100.0);
    let mut b = Box_ {
        best_im_offset_pct: f64::NAN,
        ..Box_::default()
    };
    let mut best_intensity = 0.0;
    let (mut mz_min, mut mz_max, mut im_min, mut im_max) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &f in frames {
        let Some(v) = cents.get(&f) else { continue };
        let l = v.partition_point(|c| c.mz < lo);
        for c in &v[l..] {
            if c.mz > hi {
                break;
            }
            if c.im >= im_lo && c.im <= im_hi {
                b.visible = true;
                b.intensity += c.intensity;
                if c.intensity > best_intensity {
                    best_intensity = c.intensity;
                    b.best_im_offset_pct = (c.im - im) / im * 100.0;
                }
                if f == apex {
                    b.at_apex += 1;
                    mz_min = mz_min.min(c.mz);
                    mz_max = mz_max.max(c.mz);
                    im_min = im_min.min(c.im);
                    im_max = im_max.max(c.im);
                }
            }
        }
    }
    if b.at_apex > 1 {
        b.apex_mz_spread_ppm = (mz_max - mz_min) / mz * 1e6;
        b.apex_im_spread_pct = (im_max - im_min) / im * 100.0;
    }
    b
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: isotope_frames <file.d> <seeds.csv> [flags]");
    let seeds_path = args.get(2).expect("seeds csv");
    let mut cfg = IndexingCentroidingConfig::default();
    let mut max_seeds = usize::MAX;
    let mut halo: usize = 1;
    // Keep only seeds whose apex is within `rt_window` seconds of `rt_center`.
    // With a window of one or two cycles this is the fastest possible loop:
    // a handful of frames centroided, a config comparison in about a second.
    let mut rt_center: Option<f64> = None;
    let mut rt_window: f64 = 5.0;
    // The box the envelope is read from. Independent of the centroiding
    // config, so it must be wide enough that no config under test pushes a
    // real centroid out of it; check a ranking at two widths before trusting it.
    let mut box_im_pct: f64 = MOB_PCT;
    let mut box_ms1_ppm: f64 = MS1_PPM;
    let mut box_ms2_ppm: f64 = MS2_PPM;
    let mut i = 3;
    while i < args.len() {
        let v = args.get(i + 1).expect("flag value");
        match args[i].as_str() {
            "--ms1-mz-bins" => cfg.ms1.mz_tol = MzTolerance::Bins(v.parse().unwrap()),
            "--ms2-mz-bins" => cfg.ms2.mz_tol = MzTolerance::Bins(v.parse().unwrap()),
            "--ms1-ppm" => cfg.ms1.mz_tol = MzTolerance::Ppm(v.parse().unwrap()),
            "--ms2-ppm" => cfg.ms2.mz_tol = MzTolerance::Ppm(v.parse().unwrap()),
            "--ms1-cap" => cfg.ms1.window_cap = Some(parse_cap(v)),
            "--ms2-cap" => cfg.ms2.window_cap = Some(parse_cap(v)),
            "--ms1-early-stop" => cfg.ms1.early_stop_iterations = v.parse().unwrap(),
            "--ms2-early-stop" => cfg.ms2.early_stop_iterations = v.parse().unwrap(),
            "--ms1-max-peaks" => cfg.ms1.max_peaks = v.parse().unwrap(),
            "--ms2-max-peaks" => cfg.ms2.max_peaks = v.parse().unwrap(),
            "--ms1-im-pct" => {
                cfg.ms1.im_tol = timscentroid::centroiding::ImTolerance::Pct(v.parse().unwrap())
            }
            "--ms2-im-pct" => {
                cfg.ms2.im_tol = timscentroid::centroiding::ImTolerance::Pct(v.parse().unwrap())
            }
            "--transitive" => {
                let t: bool = v.parse().unwrap();
                cfg.ms1.transitive = t;
                cfg.ms2.transitive = t;
            }
            "--seeds" => max_seeds = v.parse().unwrap(),
            "--halo" => halo = v.parse().unwrap(),
            "--rt-center" => rt_center = Some(v.parse().unwrap()),
            "--rt-window" => rt_window = v.parse().unwrap(),
            "--box-im-pct" => box_im_pct = v.parse().unwrap(),
            "--box-ms1-ppm" => box_ms1_ppm = v.parse().unwrap(),
            "--box-ms2-ppm" => box_ms2_ppm = v.parse().unwrap(),
            "--ms1-im-scans" => {
                cfg.ms1.im_tol = timscentroid::centroiding::ImTolerance::Scans(v.parse().unwrap())
            }
            "--ms2-im-scans" => {
                cfg.ms2.im_tol = timscentroid::centroiding::ImTolerance::Scans(v.parse().unwrap())
            }
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
            let seq = c
                .get(4)
                .map(|s| s.split('/').next().unwrap_or(s).to_string());
            Seed {
                mz: c[0].parse().unwrap(),
                charge: c[1].parse().unwrap(),
                rt_s: c[2].parse().unwrap(),
                mobility: c[3].parse().unwrap(),
                sequence: seq.filter(|s| !s.contains('[') && !s.is_empty()),
            }
        })
        .filter(|s| rt_center.is_none_or(|c| (s.rt_s - c).abs() <= rt_window))
        .take(max_seeds)
        .collect();
    if seeds.is_empty() {
        eprintln!("no seeds in the requested RT window");
        std::process::exit(2);
    }

    let file = TimsTofPath::new(path).unwrap();
    let frame_reader = file.load_frame_reader().unwrap();
    let metadata = file.load_metadata().unwrap();
    // Per-frame calibrations for both axes, as the index build uses them. The
    // run-level `metadata.im_converter` is up to 2.7 % off at low scan
    // indices, which would bias every mobility position read here.
    let (mz_calibrations, ims_calibrations) = metadata.get_calibration().unwrap();
    let metas = &frame_reader.frame_metas;

    println!("file: {path}");
    println!(
        "seeds: {} from {seeds_path}, halo {halo} frame(s) each side",
        seeds.len()
    );
    println!("ms1: {:?}", cfg.ms1);
    println!("ms2: {:?}", cfg.ms2);

    // --- pick frames per seed -------------------------------------------
    let ms1_frames: Vec<usize> = metas
        .iter()
        .enumerate()
        .filter(|(_, m)| m.ms_level == MSLevel::MS1)
        .map(|(i, _)| i)
        .collect();
    let nearest = |list: &[usize], rt: f64| -> usize {
        let p = list.partition_point(|&i| metas[i].rt_in_seconds < rt);
        let cands = [p.saturating_sub(1), p.min(list.len() - 1)];
        cands
            .iter()
            .map(|&k| list[k])
            .min_by(|&a, &b| {
                (metas[a].rt_in_seconds - rt)
                    .abs()
                    .total_cmp(&(metas[b].rt_in_seconds - rt).abs())
            })
            .unwrap()
    };
    let with_halo = |list: &[usize], apex: usize| -> Vec<usize> {
        let k = list.iter().position(|&i| i == apex).unwrap();
        let lo = k.saturating_sub(halo);
        let hi = (k + halo).min(list.len() - 1);
        list[lo..=hi].to_vec()
    };

    // MS2 frames grouped by window group index.
    let mut ms2_by_group: HashMap<u8, Vec<usize>> = HashMap::new();
    for (i, m) in metas.iter().enumerate() {
        if let (MSLevel::MS2, Some(wg)) = (m.ms_level, &m.window_group) {
            ms2_by_group.entry(wg.window_group).or_default().push(i);
        }
    }
    // Which group isolates (mz, scan)?
    let group_for = |mz: f64, scan: f64| -> Option<u8> {
        for (g, frames) in &ms2_by_group {
            let q = &metas[frames[0]]
                .window_group
                .as_ref()
                .unwrap()
                .quadrupole_settings;
            for k in 0..q.isolation_mz.len() {
                let half = q.isolation_width[k] / 2.0;
                if (mz - q.isolation_mz[k]).abs() <= half
                    && scan >= q.scan_starts[k] as f64
                    && scan <= q.scan_ends[k] as f64
                {
                    return Some(*g);
                }
            }
        }
        None
    };

    struct Plan {
        ms1: Vec<usize>,
        ms1_apex: usize,
        ms2: Vec<usize>,
        ms2_apex: Option<usize>,
    }
    let plans: Vec<Plan> = seeds
        .iter()
        .map(|s| {
            let a1 = nearest(&ms1_frames, s.rt_s);
            let scan = metadata.im_converter.invert(s.mobility);
            let (ms2, ms2_apex) = match group_for(s.mz, scan) {
                Some(g) => {
                    let list = &ms2_by_group[&g];
                    let a2 = nearest(list, s.rt_s);
                    (with_halo(list, a2), Some(a2))
                }
                None => (Vec::new(), None),
            };
            Plan {
                ms1: with_halo(&ms1_frames, a1),
                ms1_apex: a1,
                ms2,
                ms2_apex,
            }
        })
        .collect();

    // --- centroid every needed frame once ---------------------------------
    let mut needed: Vec<usize> = plans
        .iter()
        .flat_map(|p| p.ms1.iter().chain(p.ms2.iter()).copied())
        .collect();
    needed.sort_unstable();
    needed.dedup();
    let st = std::time::Instant::now();
    let cents: HashMap<usize, Vec<Cent>> = needed
        .par_iter()
        .map_init(
            || {
                let mk = |c: CentroidingConfig| {
                    PeakCentroider::with_capacity(
                        500_000,
                        c,
                        metadata.mz_converter,
                        metadata.im_converter,
                    )
                };
                (
                    mk(cfg.ms1),
                    mk(cfg.ms2),
                    TdfBlob::with_capacity(500_000),
                    timsrust::Frame {
                        meta: Default::default(),
                        peaks: FramePeaks::with_capacity(1000, 500_000),
                    },
                )
            },
            |(c1, c2, blob, frame), &idx| {
                frame_reader.get_buffered(idx, frame, blob).unwrap();
                let cal = mz_calibrations
                    .get_by_id(frame.meta.calibration.calibration_id)
                    .unwrap();
                let im_conv = ims_calibrations
                    .get_by_id(frame.meta.calibration.calibration_id)
                    .unwrap()
                    .get_conversion_function();
                let conv = Tof2MzConverter2::try_from_calibration(
                    cal,
                    frame.meta.calibration.t1,
                    frame.meta.calibration.t2,
                )
                .unwrap();
                let centroider = if frame.meta.ms_level == MSLevel::MS1 {
                    c1
                } else {
                    c2
                };
                let (_, peaks) = centroider.centroid_frame(frame);
                let mut v: Vec<Cent> = peaks
                    .map(|p| Cent {
                        mz: conv.convert(p.tof_index),
                        im: im_conv(p.scan_index as f64),
                        intensity: p.corrected_intensity,
                    })
                    .collect();
                v.sort_by(|a, b| a.mz.total_cmp(&b.mz));
                (idx, v)
            },
        )
        .collect();
    println!(
        "centroided {} frames in {:.1?} ({} MS1-level seeds, {} with an MS2 window)",
        cents.len(),
        st.elapsed(),
        plans.len(),
        plans.iter().filter(|p| p.ms2_apex.is_some()).count()
    );

    // --- MS1 envelopes ------------------------------------------------------
    let mut on1: Vec<[Box_; N_MS1]> = Vec::new();
    let mut off1: Vec<[Box_; N_MS1]> = Vec::new();
    let mut cos1: Vec<f64> = Vec::new();
    for (s, p) in seeds.iter().zip(&plans) {
        let z = s.charge as f64;
        let row: [Box_; N_MS1] = std::array::from_fn(|k| {
            measure(
                &cents,
                &p.ms1,
                p.ms1_apex,
                s.mz + k as f64 * ISOTOPE_DA / z,
                box_ms1_ppm,
                s.mobility,
                box_im_pct,
            )
        });
        let null: [Box_; N_MS1] = std::array::from_fn(|k| {
            measure(
                &cents,
                &p.ms1,
                p.ms1_apex,
                s.mz + (k as f64 + 0.5) * ISOTOPE_DA / z,
                box_ms1_ppm,
                s.mobility,
                box_im_pct,
            )
        });
        if row[0].visible {
            let expected = averagine::<N_MS1>((s.mz - PROTON) * z);
            let observed: [f64; N_MS1] = std::array::from_fn(|k| row[k].intensity);
            cos1.push(cosine(&expected, &observed));
        }
        on1.push(row);
        off1.push(null);
    }
    println!("\n##### MS1 precursor envelopes #####");
    report(&on1, "isotope positions");
    report(&off1, "null positions (+half isotope)");
    cosine_report(&mut cos1, on1.len());

    // --- MS2 fragment envelopes ----------------------------------------------
    let mut on2: Vec<[Box_; N_MS2]> = Vec::new();
    let mut off2: Vec<[Box_; N_MS2]> = Vec::new();
    let mut cos2: Vec<f64> = Vec::new();
    let mut n_frag = 0usize;
    for (s, p) in seeds.iter().zip(&plans) {
        let (Some(seq), Some(apex)) = (&s.sequence, p.ms2_apex) else {
            continue;
        };
        for (fmz, fz, neutral) in fragments(seq) {
            n_frag += 1;
            let z = fz as f64;
            let row: [Box_; N_MS2] = std::array::from_fn(|k| {
                measure(
                    &cents,
                    &p.ms2,
                    apex,
                    fmz + k as f64 * ISOTOPE_DA / z,
                    box_ms2_ppm,
                    s.mobility,
                    box_im_pct,
                )
            });
            if !row[0].visible {
                continue;
            }
            let null: [Box_; N_MS2] = std::array::from_fn(|k| {
                measure(
                    &cents,
                    &p.ms2,
                    apex,
                    fmz + (k as f64 + 0.5) * ISOTOPE_DA / z,
                    box_ms2_ppm,
                    s.mobility,
                    box_im_pct,
                )
            });
            let expected = averagine::<N_MS2>(neutral);
            let observed: [f64; N_MS2] = std::array::from_fn(|k| row[k].intensity);
            cos2.push(cosine(&expected, &observed));
            on2.push(row);
            off2.push(null);
        }
    }
    println!(
        "\n##### MS2 fragment envelopes: {} of {} theoretical b/y ions detected at M+0 #####",
        on2.len(),
        n_frag
    );
    report(&on2, "isotope positions (detected fragments)");
    report(&off2, "null positions (+half isotope)");
    cosine_report(&mut cos2, on2.len());

    // Within-peptide mobility spread: brightest centroid of each detected
    // fragment M+0 (and each visible precursor isotope) relative to the seed's
    // observed mobility. This, not the library-vs-observed offset, is what a
    // mobility window centered on the observed mobility has to cover.
    let mut frag_off: Vec<f64> = on2
        .iter()
        .map(|r| r[0].best_im_offset_pct)
        .filter(|x| x.is_finite())
        .collect();
    let mut iso_off: Vec<f64> = on1
        .iter()
        .flat_map(|r| r.iter().skip(1).map(|b| b.best_im_offset_pct))
        .filter(|x| x.is_finite())
        .collect();
    for (name, v) in [
        ("fragments M+0", &mut frag_off),
        ("precursor M+1..M+4", &mut iso_off),
    ] {
        if v.is_empty() {
            continue;
        }
        v.sort_by(|a, b| a.total_cmp(b));
        let med = v[v.len() / 2];
        let mut dev: Vec<f64> = v.iter().map(|x| (x - med).abs()).collect();
        dev.sort_by(|a, b| a.total_cmp(b));
        let mad = dev[dev.len() / 2];
        let within =
            |t: f64| 100.0 * v.iter().filter(|x| x.abs() <= t).count() as f64 / v.len() as f64;
        println!(
            "within-peptide mobility offset, {name} (n={}): median {:+.2} %  MAD {:.2} %  (3*1.4826*MAD = {:.2} %)   |offset| <= 1 %: {:.0} %  <= 2 %: {:.0} %  <= 3 %: {:.0} %",
            v.len(),
            med,
            mad,
            3.0 * 1.4826 * mad,
            within(1.0),
            within(2.0),
            within(3.0)
        );
    }
}

fn report<const N: usize>(rows: &[[Box_; N]], name: &str) {
    println!("-- {name} (n={}) --", rows.len());
    println!(
        "{:>8} {:>10} {:>16} {:>12} {:>22}",
        "isotope", "visible %", "centroids@apex", "split>1 %", "split spread p50 (mz ppm / im %)"
    );
    for k in 0..N {
        let vis: Vec<&Box_> = rows.iter().map(|r| &r[k]).filter(|b| b.visible).collect();
        let hits: Vec<u32> = vis.iter().map(|b| b.at_apex).filter(|&c| c > 0).collect();
        let mean = hits.iter().sum::<u32>() as f64 / hits.len().max(1) as f64;
        let split = hits.iter().filter(|&&c| c > 1).count() as f64 / hits.len().max(1) as f64;
        let mut mzs: Vec<f64> = vis
            .iter()
            .filter(|b| b.at_apex > 1)
            .map(|b| b.apex_mz_spread_ppm)
            .collect();
        let mut ims: Vec<f64> = vis
            .iter()
            .filter(|b| b.at_apex > 1)
            .map(|b| b.apex_im_spread_pct)
            .collect();
        mzs.sort_by(|a, b| a.total_cmp(b));
        ims.sort_by(|a, b| a.total_cmp(b));
        let med = |v: &[f64]| {
            if v.is_empty() {
                f64::NAN
            } else {
                v[v.len() / 2]
            }
        };
        println!(
            "{:>8} {:>9.1} % {:>16.2} {:>11.1} % {:>12.2} / {:>6.2}",
            format!("M+{k}"),
            100.0 * vis.len() as f64 / rows.len().max(1) as f64,
            mean,
            100.0 * split,
            med(&mzs),
            med(&ims)
        );
    }
}

fn cosine_report(cos: &mut [f64], n_total: usize) {
    if cos.is_empty() {
        return;
    }
    cos.sort_by(|a, b| a.total_cmp(b));
    let frac = |t: f64| 100.0 * cos.iter().filter(|&&c| c >= t).count() as f64 / cos.len() as f64;
    println!(
        "averagine cosine (n={} of {}): p10 {:.3}  p50 {:.3}  p90 {:.3}   >=0.90 {:.1} %   >=0.95 {:.1} %   >=0.98 {:.1} %",
        cos.len(),
        n_total,
        cos[cos.len() / 10],
        cos[cos.len() / 2],
        cos[cos.len() * 9 / 10],
        frac(0.90),
        frac(0.95),
        frac(0.98)
    );
}
