//! Sanity check on TOF bin width: converter-derived bin width at several m/z,
//! the raw TOF index at those m/z, and the spacing between consecutive raw
//! peaks within one scan of a crowded MS1 frame. If the stored data were a
//! sampled profile, consecutive peaks within a scan would sit 1 bin apart most
//! of the time; line spectra sit many bins apart.
//!
//! Usage: cargo run --release --example tof_bins -- <file.d>

use timsrust::converters::ConvertableDomain;
use timsrust::readers::TdfBlob;
use timsrust::{
    FramePeaks,
    MSLevel,
    TimsTofPath,
};

fn main() {
    let path = std::env::args().nth(1).expect("usage: tof_bins <file.d>");
    let file = TimsTofPath::new(&path).unwrap();
    let frame_reader = file.load_frame_reader().unwrap();
    let metadata = file.load_metadata().unwrap();
    let conv = metadata.mz_converter;

    println!("m/z      tof_index    bin_mDa   bin_ppm");
    for mz in [200.0, 300.0, 500.0, 700.0, 1000.0, 1200.0, 1700.0] {
        let t = conv.invert(mz).round();
        let d_mz = conv.convert(t + 1.0) - conv.convert(t);
        println!(
            "{:>6.0}  {:>10.0}  {:>8.3}  {:>8.2}",
            mz,
            t,
            d_mz * 1e3,
            d_mz / mz * 1e6
        );
    }

    // Crowded MS1 frame: the one nearest the middle of the gradient.
    let rt_max = frame_reader
        .frame_metas
        .iter()
        .map(|m| m.rt_in_seconds)
        .fold(0.0, f64::max);
    let idx = frame_reader
        .frame_metas
        .iter()
        .enumerate()
        .filter(|(_, m)| m.ms_level == MSLevel::MS1)
        .min_by(|a, b| {
            (a.1.rt_in_seconds - rt_max / 2.0)
                .abs()
                .total_cmp(&(b.1.rt_in_seconds - rt_max / 2.0).abs())
        })
        .map(|(i, _)| i)
        .unwrap();
    let mut blob = TdfBlob::with_capacity(500_000);
    let mut frame = timsrust::Frame {
        meta: Default::default(),
        peaks: FramePeaks::with_capacity(1000, 500_000),
    };
    frame_reader
        .get_buffered(idx, &mut frame, &mut blob)
        .unwrap();
    let p = &frame.peaks;
    println!(
        "\nMS1 frame {idx} at rt {:.0} s: {} peaks, {} scans, max tof index {}",
        frame.meta.rt_in_seconds,
        p.len(),
        p.scan_offsets.len().saturating_sub(1),
        p.tof_indices.iter().max().unwrap()
    );

    // Two ways to turn a scan index into 1/K0: the run-level converter and the
    // per-frame calibration the index build uses. If they disagree, a probe
    // using one and a search using the other are offset in mobility.
    let (_mz_cal, ims_cal) = metadata.get_calibration().unwrap();
    let per_frame = ims_cal
        .get_by_id(frame.meta.calibration.calibration_id)
        .unwrap()
        .get_conversion_function();
    println!(
        "\nscan -> 1/K0: run-level converter vs this frame's calibration (id {}):",
        frame.meta.calibration.calibration_id
    );
    for scan in [50.0, 200.0, 400.0, 600.0, 700.0] {
        let a = metadata.im_converter.convert(scan);
        let b = per_frame(scan);
        println!(
            "  scan {:>4.0}: run-level {:.4}   per-frame {:.4}   diff {:+.2} %",
            scan,
            a,
            b,
            (a - b) / b * 100.0
        );
    }

    // Consecutive-peak TOF spacing within a scan. Within a scan the TOF
    // indices are stored ascending.
    let mut hist = [0usize; 12];
    let mut n = 0usize;
    for w in p.scan_offsets.windows(2) {
        let (a, b) = (w[0] as usize, w[1] as usize);
        for i in a + 1..b {
            let d = p.tof_indices[i].saturating_sub(p.tof_indices[i - 1]) as usize;
            hist[d.min(11)] += 1;
            n += 1;
        }
    }
    println!("spacing between consecutive raw peaks within a scan (TOF bins):");
    for (d, c) in hist.iter().enumerate() {
        let label = if d == 11 {
            ">=11".to_string()
        } else {
            d.to_string()
        };
        println!(
            "  {:>4}: {:>6.1} %",
            label,
            100.0 * *c as f64 / n.max(1) as f64
        );
    }
    println!("intensity: p50 {} (raw counts)", {
        let mut v = p.intensities.clone();
        v.sort_unstable();
        v[v.len() / 2]
    });
}
