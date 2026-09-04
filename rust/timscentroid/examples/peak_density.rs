//! Peak density probe: how many centroided peaks land in each 100 Da m/z
//! window, per frame, per MS level, with the default centroiding config.
//!
//! Answers "would a per-window cap of N peaks / 100 Da bite, and where?"
//! before any cap is implemented. Streams frames, so memory stays flat.
//!
//! Usage: cargo run --release --example peak_density -- <file.d> [frame_stride]

use rayon::prelude::*;
use timscentroid::centroiding::{
    PeakCentroider,
    StoppingReason,
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

const WINDOW_DA: f64 = 100.0;
const MAX_MZ: f64 = 2000.0;
const N_WINDOWS: usize = (MAX_MZ / WINDOW_DA) as usize;
const CAPS: &[usize] = &[50, 100, 200, 500, 1000, 2000];

#[derive(Default, Clone)]
struct FrameStat {
    raw_peaks: usize,
    centroided: usize,
    per_window: [u32; N_WINDOWS],
    reason: Option<StoppingReason>,
}

fn pct(v: &mut [u32], p: f64) -> u32 {
    if v.is_empty() {
        return 0;
    }
    v.sort_unstable();
    v[((v.len() - 1) as f64 * p).round() as usize]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: peak_density <file.d> [frame_stride]");
    let stride: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(1);

    let file = TimsTofPath::new(path).unwrap();
    let frame_reader = file.load_frame_reader().unwrap();
    let metadata = file.load_metadata().unwrap();
    let (mz_calibrations, _ims) = metadata.get_calibration().unwrap();
    // Uncapped on purpose: the point is to see the density a cap would act on.
    let mut cfg = IndexingCentroidingConfig::default();
    cfg.ms1.window_cap = None;
    cfg.ms2.window_cap = None;
    println!("file: {path}");
    println!(
        "frames: {}, stride: {stride}",
        frame_reader.frame_metas.len()
    );
    println!("config ms1: {:?}", cfg.ms1);
    println!("config ms2: {:?}", cfg.ms2);

    for (level, level_cfg) in [(MSLevel::MS1, cfg.ms1), (MSLevel::MS2, cfg.ms2)] {
        let indices: Vec<usize> = frame_reader
            .frame_metas
            .iter()
            .enumerate()
            .filter(|(_, m)| m.ms_level == level)
            .map(|(i, _)| i)
            .step_by(stride)
            .collect();
        let st = std::time::Instant::now();
        let stats: Vec<FrameStat> = indices
            .par_iter()
            .with_min_len(50)
            .map_init(
                || {
                    (
                        PeakCentroider::with_capacity(
                            500_000,
                            level_cfg,
                            metadata.mz_converter,
                            metadata.im_converter,
                        ),
                        TdfBlob::with_capacity(500_000),
                        timsrust::Frame {
                            meta: Default::default(),
                            peaks: FramePeaks::with_capacity(1000, 500_000),
                        },
                    )
                },
                |(centroider, blob, frame), &idx| {
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
                    let raw_peaks = frame.peaks.len();
                    let (summary, peaks) = centroider.centroid_frame(frame);
                    let mut fs = FrameStat {
                        raw_peaks,
                        reason: Some(summary.stopping_reason),
                        ..Default::default()
                    };
                    for p in peaks {
                        let mz = conv.convert(p.tof_index as f64);
                        let w = ((mz / WINDOW_DA) as usize).min(N_WINDOWS - 1);
                        fs.per_window[w] += 1;
                        fs.centroided += 1;
                    }
                    fs
                },
            )
            .collect();
        let elapsed = st.elapsed();
        report(level, level_cfg, &stats, elapsed);
    }
}

fn report(
    level: MSLevel,
    cfg: CentroidingConfig,
    stats: &[FrameStat],
    elapsed: std::time::Duration,
) {
    let n = stats.len();
    println!("\n===== {level:?}: {n} frames centroided in {elapsed:.1?} =====");
    if n == 0 {
        return;
    }
    let mut raw: Vec<u32> = stats.iter().map(|s| s.raw_peaks as u32).collect();
    let mut cen: Vec<u32> = stats.iter().map(|s| s.centroided as u32).collect();
    let total_cen: u64 = cen.iter().map(|&x| x as u64).sum();
    let total_raw: u64 = raw.iter().map(|&x| x as u64).sum();
    println!(
        "raw peaks/frame:        p50 {:>7}  p95 {:>7}  max {:>7}  total {:>12}",
        pct(&mut raw, 0.5),
        pct(&mut raw, 0.95),
        raw.last().unwrap(),
        total_raw
    );
    println!(
        "centroided peaks/frame: p50 {:>7}  p95 {:>7}  max {:>7}  total {:>12}  (max_peaks={})",
        pct(&mut cen, 0.5),
        pct(&mut cen, 0.95),
        cen.last().unwrap(),
        total_cen,
        cfg.max_peaks
    );
    let mut reasons = [0usize; 3];
    for s in stats {
        match s.reason {
            Some(StoppingReason::EarlyStop) => reasons[0] += 1,
            Some(StoppingReason::MaxPeaks) => reasons[1] += 1,
            Some(StoppingReason::AllTaken) => reasons[2] += 1,
            None => {}
        }
    }
    println!(
        "stopping reason: early_stop {}  max_peaks {}  all_taken {}",
        reasons[0], reasons[1], reasons[2]
    );

    println!("\nper 100 Da window, centroided peaks per frame (over frames):");
    println!(
        "{:>10} {:>8} {:>8} {:>8} {:>8}",
        "window", "p50", "p95", "max", "mean"
    );
    for w in 0..N_WINDOWS {
        let mut v: Vec<u32> = stats.iter().map(|s| s.per_window[w]).collect();
        let sum: u64 = v.iter().map(|&x| x as u64).sum();
        if sum == 0 {
            continue;
        }
        let p50 = pct(&mut v, 0.5);
        let p95 = pct(&mut v, 0.95);
        println!(
            "{:>4}-{:<5} {:>8} {:>8} {:>8} {:>8.1}",
            (w as f64 * WINDOW_DA) as u32,
            ((w + 1) as f64 * WINDOW_DA) as u32,
            p50,
            p95,
            v.last().unwrap(),
            sum as f64 / n as f64
        );
    }

    println!("\nretained under a per-window cap (top-N by intensity per 100 Da per frame):");
    println!("{:>6} {:>14} {:>8}", "cap", "peaks", "kept%");
    for &cap in CAPS {
        let kept: u64 = stats
            .iter()
            .map(|s| {
                s.per_window
                    .iter()
                    .map(|&c| (c as usize).min(cap) as u64)
                    .sum::<u64>()
            })
            .sum();
        println!(
            "{:>6} {:>14} {:>7.1}%",
            cap,
            kept,
            100.0 * kept as f64 / total_cen.max(1) as f64
        );
    }
}
