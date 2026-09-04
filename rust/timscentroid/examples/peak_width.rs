//! Raw peak width probe: how wide is a real peak in m/z and in mobility, in
//! the units the centroider's tolerances are expressed in?
//!
//! For each sampled frame, the brightest raw peaks are anchors. Around each
//! anchor the connected run of raw peaks along TOF (same scan, gap <= 1 bin)
//! gives the m/z profile width, and the connected run along scan index (same
//! TOF run, gap <= 2 scans) gives the mobility width. Everything a run touches
//! is claimed so one ion is measured once. Frames come from the middle of the
//! gradient, where they are crowded.
//!
//! The question answered: does the configured `mz_tol` / `im_tol` cover one
//! raw peak, or does it slice through it? Note the raw data is line spectra
//! along TOF (one index per ion per scan), so the "m/z profile" here is the
//! scan-to-scan quantization jitter, not a sampled peak shape; along scans it
//! is a real profile.
//!
//! Usage: cargo run --release --example peak_width -- <file.d> [max_frames]

use rayon::prelude::*;
use timscentroid::IndexingCentroidingConfig;
use timscentroid::centroiding::{
    ImTolerance,
    MzTolerance,
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

const ANCHORS_PER_FRAME: usize = 400;
const GATHER_PPM: f64 = 30.0;
const SCAN_GAP: i64 = 2;
const TOF_GAP: u32 = 1;
const RT_FRACTION: (f64, f64) = (0.3, 0.7);

#[derive(Clone, Copy)]
struct Col {
    anchor_intensity: u32,
    mz: f64,
    mz_width_ppm: f64,
    tof_bins: u32,
    im_width_scans: u32,
    im_half_pct: f64,
    /// TOF drift across the mobility cloud: span of TOF indices (within +-4
    /// bins of the anchor) over the scans of the mobility run, and the
    /// Pearson correlation of TOF index with scan index over those peaks.
    cloud_tof_span: u32,
    cloud_corr: f64,
    cloud_n: u32,
}

#[derive(Clone, Copy)]
struct P {
    scan: u16,
    tof: u32,
    intensity: u32,
}

fn pct(v: &mut [f64], p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.total_cmp(b));
    v[((v.len() - 1) as f64 * p).round() as usize]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: peak_width <file.d> [max_frames]");
    let max_frames: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(150);

    let file = TimsTofPath::new(path).unwrap();
    let frame_reader = file.load_frame_reader().unwrap();
    let metadata = file.load_metadata().unwrap();
    let (mz_calibrations, _ims) = metadata.get_calibration().unwrap();
    let cfg = IndexingCentroidingConfig::default();
    let rt_max = frame_reader
        .frame_metas
        .iter()
        .map(|m| m.rt_in_seconds)
        .fold(0.0, f64::max);
    println!("file: {path}");
    println!(
        "frames sampled from rt {:.0}-{:.0} s of {:.0} s, up to {max_frames} per level",
        RT_FRACTION.0 * rt_max,
        RT_FRACTION.1 * rt_max,
        rt_max
    );

    // TOF bin size in ppm at a few m/z, using the global converter.
    println!("\nTOF bin width (one raw bin, ppm of m/z) with the run's converter:");
    for mz in [300.0, 500.0, 800.0, 1200.0] {
        let t = metadata.mz_converter.invert(mz).round();
        let w =
            (metadata.mz_converter.convert(t + 1.0) - metadata.mz_converter.convert(t)) / mz * 1e6;
        println!("  m/z {mz:>6.0}: {w:.2} ppm/bin");
    }
    println!("mobility per scan at scan 400: {:.5} 1/K0/scan", {
        let a = metadata.im_converter.convert(400.0);
        let b = metadata.im_converter.convert(401.0);
        (a - b).abs()
    });

    for (level, level_cfg) in [(MSLevel::MS1, cfg.ms1), (MSLevel::MS2, cfg.ms2)] {
        let candidates: Vec<usize> = frame_reader
            .frame_metas
            .iter()
            .enumerate()
            .filter(|(_, m)| {
                m.ms_level == level
                    && m.rt_in_seconds >= RT_FRACTION.0 * rt_max
                    && m.rt_in_seconds <= RT_FRACTION.1 * rt_max
            })
            .map(|(i, _)| i)
            .collect();
        let stride = (candidates.len() / max_frames).max(1);
        let indices: Vec<usize> = candidates
            .into_iter()
            .step_by(stride)
            .take(max_frames)
            .collect();

        let results: Vec<(usize, Vec<Col>)> = indices
            .par_iter()
            .map_init(
                || {
                    (
                        TdfBlob::with_capacity(500_000),
                        timsrust::Frame {
                            meta: Default::default(),
                            peaks: FramePeaks::with_capacity(1000, 500_000),
                        },
                    )
                },
                |(blob, frame), &idx| {
                    frame_reader.get_buffered(idx, frame, blob).unwrap();
                    let cal = mz_calibrations
                        .get_by_id(frame.meta.calibration.calibration_id)
                        .unwrap();
                    let conv = Tof2MzConverter2::try_from_calibration(
                        cal,
                        frame.meta.calibration.t1,
                        frame.meta.calibration.t2,
                    )
                    .unwrap();
                    let cols = measure_frame(frame, &conv, &metadata.im_converter);
                    (frame.peaks.len(), cols)
                },
            )
            .collect();

        report(level, level_cfg, &results);
    }
}

fn measure_frame(
    frame: &timsrust::Frame,
    conv: &Tof2MzConverter2,
    im_conv: &impl ConvertableDomain,
) -> Vec<Col> {
    let mut peaks: Vec<P> = frame
        .peaks
        .iter_peaks()
        .map(|p| P {
            scan: p.scan_index,
            tof: p.tof_index,
            intensity: p.intensity,
        })
        .collect();
    peaks.sort_unstable_by_key(|p| p.tof);
    let n = peaks.len();
    let mut by_int: Vec<usize> = (0..n).collect();
    by_int.sort_unstable_by(|&a, &b| peaks[b].intensity.cmp(&peaks[a].intensity));
    let mut claimed = vec![false; n];
    let mut cols = Vec::with_capacity(ANCHORS_PER_FRAME);

    for &a in &by_int {
        if cols.len() >= ANCHORS_PER_FRAME {
            break;
        }
        if claimed[a] {
            continue;
        }
        let anchor = peaks[a];
        let mz_a = conv.convert(anchor.tof);
        let tof_lo = conv.invert(mz_a * (1.0 - GATHER_PPM / 1e6)).floor() as u32;
        let tof_hi = conv.invert(mz_a * (1.0 + GATHER_PPM / 1e6)).ceil() as u32;
        let l = peaks.partition_point(|p| p.tof < tof_lo);
        let r = peaks.partition_point(|p| p.tof <= tof_hi);
        let window = &peaks[l..r];

        // m/z profile: same scan (+-1), connected in TOF around the anchor.
        let mut tofs: Vec<u32> = window
            .iter()
            .filter(|p| (p.scan as i64 - anchor.scan as i64).abs() <= 1)
            .map(|p| p.tof)
            .collect();
        tofs.sort_unstable();
        tofs.dedup();
        let (tof_min, tof_max) = connected_run(&tofs, anchor.tof, TOF_GAP as i64);

        // mobility profile: within that TOF run, connected in scan.
        let mut scans: Vec<u32> = window
            .iter()
            .filter(|p| p.tof >= tof_min && p.tof <= tof_max)
            .map(|p| p.scan as u32)
            .collect();
        scans.sort_unstable();
        scans.dedup();
        let (scan_min, scan_max) = connected_run(&scans, anchor.scan as u32, SCAN_GAP);

        // Drift: does the TOF index move with scan across the cloud?
        let cloud: Vec<(f64, f64)> = window
            .iter()
            .filter(|p| {
                (p.scan as u32) >= scan_min
                    && (p.scan as u32) <= scan_max
                    && (p.tof as i64 - anchor.tof as i64).abs() <= 4
            })
            .map(|p| (p.scan as f64, p.tof as f64))
            .collect();
        let (cloud_tof_span, cloud_corr) = if cloud.len() >= 4 {
            let tmin = cloud.iter().map(|c| c.1).fold(f64::MAX, f64::min);
            let tmax = cloud.iter().map(|c| c.1).fold(f64::MIN, f64::max);
            let n = cloud.len() as f64;
            let ms = cloud.iter().map(|c| c.0).sum::<f64>() / n;
            let mt = cloud.iter().map(|c| c.1).sum::<f64>() / n;
            let cov: f64 = cloud.iter().map(|c| (c.0 - ms) * (c.1 - mt)).sum();
            let vs: f64 = cloud.iter().map(|c| (c.0 - ms).powi(2)).sum();
            let vt: f64 = cloud.iter().map(|c| (c.1 - mt).powi(2)).sum();
            let corr = if vs > 0.0 && vt > 0.0 {
                cov / (vs * vt).sqrt()
            } else {
                0.0
            };
            ((tmax - tmin) as u32, corr)
        } else {
            (0, f64::NAN)
        };

        // Claim everything in the box so the ion is measured once.
        for (i, p) in window.iter().enumerate() {
            if p.tof >= tof_min
                && p.tof <= tof_max
                && (p.scan as u32) >= scan_min
                && (p.scan as u32) <= scan_max
            {
                claimed[l + i] = true;
            }
        }
        claimed[a] = true;

        let mz_width_ppm = (conv.convert(tof_max + 1) - conv.convert(tof_min)) / mz_a * 1e6;
        let im_a = im_conv.convert(anchor.scan as f64);
        let im_lo = im_conv.convert(scan_min as f64);
        let im_hi = im_conv.convert(scan_max as f64 + 1.0);
        let im_half_pct = (im_hi - im_lo).abs() / 2.0 / im_a * 100.0;
        cols.push(Col {
            anchor_intensity: anchor.intensity,
            mz: mz_a,
            mz_width_ppm,
            tof_bins: tof_max - tof_min + 1,
            im_width_scans: scan_max - scan_min + 1,
            im_half_pct,
            cloud_tof_span,
            cloud_corr,
            cloud_n: cloud.len() as u32,
        });
    }
    cols
}

/// Extent of the run of consecutive values (gap <= `gap`) containing `at`,
/// over a sorted, deduplicated slice.
fn connected_run(sorted: &[u32], at: u32, gap: i64) -> (u32, u32) {
    let pos = sorted.partition_point(|&v| v < at);
    let mut lo = pos;
    while lo > 0 && sorted[lo] as i64 - sorted[lo - 1] as i64 <= gap {
        lo -= 1;
    }
    let mut hi = pos;
    while hi + 1 < sorted.len() && sorted[hi + 1] as i64 - sorted[hi] as i64 <= gap {
        hi += 1;
    }
    (sorted[lo], sorted[hi])
}

fn report(level: MSLevel, cfg: timscentroid::CentroidingConfig, results: &[(usize, Vec<Col>)]) {
    let mut raw: Vec<f64> = results.iter().map(|(n, _)| *n as f64).collect();
    let cols: Vec<Col> = results
        .iter()
        .flat_map(|(_, c)| c.iter().copied())
        .collect();
    println!(
        "\n===== {level:?}: {} frames, raw peaks/frame p50 {:.0}, {} anchors =====",
        results.len(),
        pct(&mut raw, 0.5),
        cols.len()
    );
    println!(
        "tolerances in use: mz {:?}, mobility {:?}",
        cfg.mz_tol, cfg.im_tol
    );
    // "Wider than tolerance" is answered in the tolerance's own unit.
    let mz_wider = |c: &Col| match cfg.mz_tol {
        MzTolerance::Ppm(p) => c.mz_width_ppm / 2.0 > p,
        MzTolerance::Bins(b) => (c.tof_bins.saturating_sub(1)) / 2 > b,
    };
    let im_wider = |c: &Col| match cfg.im_tol {
        ImTolerance::Pct(p) => c.im_half_pct > p,
        ImTolerance::Scans(s) => (c.im_width_scans.saturating_sub(1)) / 2 > s as u32,
    };

    // Bright vs dim anchors, so the answer is not dominated by noise anchors.
    let mut sorted = cols.clone();
    sorted.sort_by_key(|c| std::cmp::Reverse(c.anchor_intensity));
    let cut = sorted.len() / 4;
    for (name, slice) in [
        ("all anchors", &sorted[..]),
        ("brightest quarter", &sorted[..cut.max(1)]),
        ("dimmest quarter", &sorted[sorted.len() - cut.max(1)..]),
    ] {
        let mut mzw: Vec<f64> = slice.iter().map(|c| c.mz_width_ppm).collect();
        let mut bins: Vec<f64> = slice.iter().map(|c| c.tof_bins as f64).collect();
        let mut imw: Vec<f64> = slice.iter().map(|c| c.im_width_scans as f64).collect();
        let mut imp: Vec<f64> = slice.iter().map(|c| c.im_half_pct).collect();
        let n = slice.len() as f64;
        let mz_split = slice.iter().filter(|c| mz_wider(c)).count() as f64 / n;
        let im_split = slice.iter().filter(|c| im_wider(c)).count() as f64 / n;
        let single_bin = slice.iter().filter(|c| c.tof_bins == 1).count() as f64 / n;
        let single_scan = slice.iter().filter(|c| c.im_width_scans == 1).count() as f64 / n;
        let mut mzs: Vec<f64> = slice.iter().map(|c| c.mz).collect();
        println!(
            "\n-- {name} (n={}, anchor m/z p50 {:.0}) --",
            slice.len(),
            pct(&mut mzs, 0.5)
        );
        println!(
            "m/z profile width:  p10 {:>5.2}  p50 {:>5.2}  p90 {:>5.2} ppm   ({:.1} / {:.1} / {:.1} TOF bins)",
            pct(&mut mzw, 0.1),
            pct(&mut mzw, 0.5),
            pct(&mut mzw, 0.9),
            pct(&mut bins, 0.1),
            pct(&mut bins, 0.5),
            pct(&mut bins, 0.9)
        );
        println!(
            "mobility width:     p10 {:>5.0}  p50 {:>5.0}  p90 {:>5.0} scans  (half-width p50 {:.2} %, p90 {:.2} %)",
            pct(&mut imw, 0.1),
            pct(&mut imw, 0.5),
            pct(&mut imw, 0.9),
            pct(&mut imp, 0.5),
            pct(&mut imp, 0.9)
        );
        println!(
            "wider than tolerance: m/z {:.0} %   mobility {:.0} %     single TOF bin {:.0} %   single scan {:.0} %",
            mz_split * 100.0,
            im_split * 100.0,
            single_bin * 100.0,
            single_scan * 100.0
        );
        // Drift over clouds with at least 8 raw peaks (enough to see a trend).
        let clouds: Vec<&Col> = slice.iter().filter(|c| c.cloud_n >= 8).collect();
        let mut span: Vec<f64> = clouds.iter().map(|c| c.cloud_tof_span as f64).collect();
        let mut corr: Vec<f64> = clouds.iter().map(|c| c.cloud_corr).collect();
        let strong = clouds.iter().filter(|c| c.cloud_corr.abs() > 0.5).count() as f64
            / clouds.len().max(1) as f64;
        let pos = clouds.iter().filter(|c| c.cloud_corr > 0.5).count() as f64
            / clouds.len().max(1) as f64;
        println!(
            "TOF drift across the mobility cloud (n={} clouds with >=8 peaks): span p50 {:.0}  p90 {:.0} bins;  corr(scan,tof) p50 {:+.2};  |corr|>0.5 in {:.0} % (positive {:.0} %)",
            clouds.len(),
            pct(&mut span, 0.5),
            pct(&mut span, 0.9),
            pct(&mut corr, 0.5),
            strong * 100.0,
            pos * 100.0
        );
    }
}
