// Throwaway probe for the DIA-NN `.speclib` binary reader (serde/diann_speclib_io.rs).
// Eyeball the first entries against reference_parser.py / the paired TSV oracle,
// and exercise the embedded-decoy (dc != 0) branch via a synthetic buffer.
//
// Run against a real library (no local-path fallback in code):
//   SPECLIB_PATH=reverse_eng_speclib/empirical_library.tsv.speclib \
//     cargo run -p timsquery --example probe_speclib
// Print more entries with N=... :
//   SPECLIB_PATH=... N=10 cargo run -p timsquery --example probe_speclib
use std::io::Cursor;

use timsquery::serde::diann_speclib_io::parse_speclib_reader;

fn main() {
    // The synthetic dc=1 check needs no external file — always run it first.
    check_embedded_decoy_branch();

    let Ok(path) = std::env::var("SPECLIB_PATH") else {
        println!("\nSPECLIB_PATH not set — ran only the synthetic dc=1 check.");
        return;
    };
    let n: usize = std::env::var("N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    let file = std::fs::File::open(&path).expect("open SPECLIB_PATH");
    let (entries, stats, at_eof) =
        parse_speclib_reader(std::io::BufReader::new(file)).expect("parse speclib");

    println!("\n== {} ==", path);
    println!("entries={}  landed_on_eof={}", entries.len(), at_eof);
    println!(
        "exclude-flagged(kept)={}  dropped: neutral-loss={} unknown-ion={} duplicate-label={}  embedded-decoys={}",
        stats.exclude_flagged,
        stats.loss_dropped,
        stats.unknown_ion_dropped,
        stats.dedup_dropped,
        stats.decoys_dropped
    );

    for (i, (eg, extra)) in entries.iter().take(n).enumerate() {
        println!(
            "\nentry[{i}] {}  mz={:.5} z={} rt(iRT)={:.4} im(iIM)={:.5}  frags={}",
            extra.modified_peptide,
            eg.precursor_mz(),
            eg.precursor_charge(),
            eg.rt_seconds(),
            eg.mobility_ook0(),
            eg.fragment_count(),
        );
        for (label, mz) in eg.iter_fragments() {
            let inten = extra
                .relative_intensities
                .iter()
                .find(|(ion, _)| ion == label)
                .map(|(_, v)| *v)
                .unwrap_or(f32::NAN);
            println!("    {label}  mz={mz:.5}  i={inten:.4}");
        }
    }

    if !at_eof {
        println!("\nWARNING: parse did not land on EOF — structural desync suspected.");
    }
}

// ---- synthetic v-3 buffer with one dc=1 entry ------------------------------

fn push_i32(b: &mut Vec<u8>, v: i32) {
    b.extend_from_slice(&v.to_le_bytes());
}
fn push_f32(b: &mut Vec<u8>, v: f32) {
    b.extend_from_slice(&v.to_le_bytes());
}
fn push_f64(b: &mut Vec<u8>, v: f64) {
    b.extend_from_slice(&v.to_le_bytes());
}
fn push_str(b: &mut Vec<u8>, s: &str) {
    push_i32(b, s.len() as i32);
    b.extend_from_slice(s.as_bytes());
}
/// A 12-byte Product: mz, height, charge, type, index, loss.
fn push_frag(b: &mut Vec<u8>, mz: f32, height: f32, charge: u8, typ: u8, index: u8, loss: u8) {
    push_f32(b, mz);
    push_f32(b, height);
    b.push(charge);
    b.push(typ);
    b.push(index);
    b.push(loss);
}
/// A v<=-2 / v<=-3 Peptide record. `nfrag` fragments of a placeholder shape.
fn push_peptide(b: &mut Vec<u8>, charge: i32, length: i32, mz: f32, nfrag: usize) {
    push_i32(b, 0); // index
    push_i32(b, charge);
    push_i32(b, length);
    push_f32(b, mz); // mz
    push_f32(b, 12.5); // iRT
    push_f32(b, 0.0); // sRT
    push_f32(b, 0.001); // lib_qvalue (v<=-2)
    push_f32(b, 0.95); // iIM       (v<=-2)
    push_f32(b, 0.0); // sIM        (v<=-2)
    push_i32(b, nfrag as i32);
    for i in 0..nfrag {
        // y(length-1), y(length-2), ... : type=2 (y), index = i+1
        push_frag(
            b,
            500.0 + i as f32,
            1.0 - 0.1 * i as f32,
            1,
            2,
            (i + 1) as u8,
            0,
        );
    }
}

fn check_embedded_decoy_branch() {
    let mut b: Vec<u8> = Vec::new();

    // header (v-3)
    push_i32(&mut b, -3); // version
    push_i32(&mut b, 1); // gen_decoys
    push_i32(&mut b, 0); // gen_charges
    push_i32(&mut b, 0); // infer_proteotypicity
    push_str(&mut b, "synthetic"); // name
    push_str(&mut b, ""); // fasta_names

    push_i32(&mut b, 0); // section 1: 0 isoforms
    push_i32(&mut b, 0); // section 2: 0 PGs
    push_i32(&mut b, 0); // section 3: 0 precursor strings
    push_i32(&mut b, 0); // section 4: 0 names
    push_i32(&mut b, 0); // section 5: 0 genes
    push_f64(&mut b, 0.0); // section 6: iRT min
    push_f64(&mut b, 100.0); // iRT max

    // section 7: 1 entry, dc = 1
    push_i32(&mut b, 1);
    push_peptide(&mut b, 2, 6, 700.5, 2); // target peptide, 2 frags
    push_i32(&mut b, 1); // dc != 0
    push_peptide(&mut b, 2, 6, 700.5, 1); // embedded DECOY peptide (must be read to stay synced)
    push_i32(&mut b, 0); // entry_flags
    push_i32(&mut b, 1); // proteotypic
    push_i32(&mut b, -1); // pid_index
    push_str(&mut b, "PEPTIK2"); // name (charge suffix "2" -> modified "PEPTIK")
    push_f32(&mut b, 0.01); // pg_qvalue (v<=-3)
    push_f32(&mut b, 1.0); // ptm_qvalue
    push_f32(&mut b, 1.0); // site_conf

    // section 8 (v<=-1): elution_groups
    push_i32(&mut b, 1);
    push_i32(&mut b, 0);

    let (entries, stats, at_eof) = parse_speclib_reader(Cursor::new(b)).expect("parse synthetic");

    assert!(
        at_eof,
        "synthetic dc=1 buffer must land on EOF (stream stayed synced)"
    );
    assert_eq!(entries.len(), 1, "only the target must be emitted");
    assert_eq!(
        stats.decoys_dropped, 1,
        "the embedded decoy must be read-then-discarded"
    );
    let (eg, extra) = &entries[0];
    assert_eq!(extra.modified_peptide, "PEPTIK");
    assert_eq!(
        eg.fragment_count(),
        2,
        "target frags kept, decoy frags not leaked"
    );

    println!(
        "synthetic dc=1 check OK: 1 target emitted, decoy discarded, stream synced to EOF (target frags={})",
        eg.fragment_count()
    );
}
