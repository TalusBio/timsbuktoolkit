// Throwaway: probe a real Thermo/Astral mzML to validate the mzdata ingest assumptions
// in .plans/mzml-support-design.md §3. Run:
//   cargo run --release -p timscentroid --features mzdata --example probe_mzml -- <path.mzML>
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;

use mzdata::io::mzml::MzMLReader;
use mzdata::prelude::*;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: probe_mzml <path.mzML>");
    let cap: usize = std::env::var("CAP")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(8000);

    let reader = MzMLReader::new(BufReader::new(File::open(&path).unwrap()));

    let mut n = 0usize;
    let mut ms1 = 0usize;
    let mut ms2 = 0usize;
    let mut cont: BTreeMap<String, usize> = BTreeMap::new();
    let mut windows: BTreeMap<(i64, i64, i64, String), usize> = BTreeMap::new();
    let mut rts: Vec<f64> = Vec::new();
    let mut any_im = false;
    let mut im_example: Option<f64> = None;
    let mut unsorted_ms2 = 0usize;
    let mut peak_total = 0usize;
    let mut ms2_seen = 0usize;

    for spec in reader {
        n += 1;
        let lvl = spec.ms_level();
        let sc = format!("{:?}", spec.signal_continuity());
        *cont.entry(sc).or_default() += 1;

        if spec.has_ion_mobility() {
            any_im = true;
            if im_example.is_none() {
                im_example = spec.ion_mobility();
                // what param is mzdata keying IM off of?
                if let Some(scan) = spec.acquisition().first_scan() {
                    println!("-- IM debug (first spectrum w/ has_ion_mobility) --");
                    println!("  ion_mobility()={:?}", spec.ion_mobility());
                    println!("  ion_mobility_type param={:?}", scan.ion_mobility_type());
                    for p in scan.params() {
                        println!(
                            "  scan.param: acc={:?} name={:?} val={:?} unit={:?}",
                            p.accession, p.name, p.value, p.unit
                        );
                    }
                }
            }
        }

        if lvl == 1 {
            ms1 += 1;
            if rts.len() < 30 {
                rts.push(spec.start_time());
            }
        } else {
            ms2 += 1;
            if let Some(prec) = spec.precursor() {
                let w = prec.isolation_window();
                let key = (
                    (w.target as f64 * 100.0).round() as i64,
                    (w.lower_bound as f64 * 100.0).round() as i64,
                    (w.upper_bound as f64 * 100.0).round() as i64,
                    format!("{:?}", w.flags),
                );
                *windows.entry(key).or_default() += 1;
            }
            // sortedness of the peak m/z stream (first few MS2 only)
            if ms2_seen < 50 {
                ms2_seen += 1;
                let peaks = spec.peaks();
                let mut last = f64::NEG_INFINITY;
                let mut bad = false;
                let mut cnt = 0usize;
                for p in peaks.iter() {
                    if p.mz < last {
                        bad = true;
                    }
                    last = p.mz;
                    cnt += 1;
                }
                peak_total += cnt;
                if bad {
                    unsorted_ms2 += 1;
                }
            }
        }

        if n >= cap {
            break;
        }
    }

    println!("== scanned {n} spectra (cap={cap}) ==");
    println!("ms1={ms1}  ms2={ms2}");
    println!("signal_continuity: {cont:?}");
    println!("has_ion_mobility(any)={any_im}  example={im_example:?}");
    println!("distinct isolation windows: {}", windows.len());
    println!("unsorted_ms2 (of first 50)={unsorted_ms2}  sampled_peaks={peak_total}");

    println!("-- first RT values (start_time, mzdata-normalized) --");
    for (i, rt) in rts.iter().enumerate() {
        let d = if i > 0 { rt - rts[i - 1] } else { 0.0 };
        println!("  rt[{i}]={rt:.6}  delta={d:.6}");
    }

    println!("-- isolation windows (target, lower, upper, flags) x count [first 20 by m/z] --");
    for ((t, lo, hi, fl), c) in windows.iter().take(20) {
        println!(
            "  target={:.2} lower={:.2} upper={:.2} flags={} width={:.2} x{}",
            *t as f64 / 100.0,
            *lo as f64 / 100.0,
            *hi as f64 / 100.0,
            fl,
            (*hi - *lo) as f64 / 100.0,
            c
        );
    }
}
