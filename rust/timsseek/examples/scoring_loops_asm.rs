//! Standalone playground for chunked/SIMD variants of the scoring hot
//! loops. Throwaway — delete once the winning shape is folded into the
//! real `apex_finding.rs` paths.
//!
//! Run:          cargo run -r -p timsseek --example scoring_loops_asm
//! Dump ASM:     cargo rustc -r -p timsseek --example scoring_loops_asm -- --emit=asm
//!   -> ASM lands under target/release/examples/scoring_loops_asm-*.s
//!
//! Layout: two micro-benches, each comparing a scalar form of a real
//! hot loop against a chunked form. Each variant runs over ~200
//! independent synthetic `Vec<f32>` sets so the timer averages over
//! realistic data variation. Sizes (8 fragments × 64 cycles) match
//! production: TOP_N_FRAGMENTS=8, RT window ≈30-60s / cycle ≈0.7s.

use std::hint::black_box;
use std::time::Instant;

const N_CASES: usize = 200;
const N_CYCLES: usize = 64;
const N_ITERS: usize = 10_000;

// ---------------------------------------------------------------------------
// Bench 1 — compute_pass_1 inner loop: accumulate 4 parallel f32 vecs
// with a gating `if intensity > 0.0`.
// ---------------------------------------------------------------------------

/// Mirrors the current production body at apex_finding.rs:545-557.
/// Branching `if intensity > 0.0 { ... }` on every element — blocks autovec.
#[inline(never)]
pub fn pass1_scalar(
    chrom: &[f32],
    sqrt_exp: f32,
    ms2_dot_prod: &mut [f32],
    ms2_norm_sq_obs: &mut [f32],
    sqrt_sum: &mut [f32],
    raw_sum: &mut [f32],
) {
    for (i, &intensity) in chrom.iter().enumerate() {
        if intensity > 0.0 {
            ms2_dot_prod[i] += intensity * sqrt_exp;
            ms2_norm_sq_obs[i] += intensity * intensity;
            sqrt_sum[i] += intensity.sqrt();
        }
        raw_sum[i] += intensity.max(0.0);
    }
}

/// Chunked variant generic over chunk width. Compute predicated
/// arithmetic via `v.max(0.0)` instead of a branch. All four
/// accumulators written in lockstep.
#[inline(always)]
fn pass1_chunked_n<const N: usize>(
    chrom: &[f32],
    sqrt_exp: f32,
    ms2_dot_prod: &mut [f32],
    ms2_norm_sq_obs: &mut [f32],
    sqrt_sum: &mut [f32],
    raw_sum: &mut [f32],
) {
    let (c_ch, c_tail) = chrom.as_chunks::<N>();
    let (d_ch, d_tail) = ms2_dot_prod.as_chunks_mut::<N>();
    let (n_ch, n_tail) = ms2_norm_sq_obs.as_chunks_mut::<N>();
    let (s_ch, s_tail) = sqrt_sum.as_chunks_mut::<N>();
    let (r_ch, r_tail) = raw_sum.as_chunks_mut::<N>();
    for ((((c, d), n), s), r) in c_ch.iter().zip(d_ch).zip(n_ch).zip(s_ch).zip(r_ch) {
        for k in 0..N {
            let v = c[k].max(0.0);
            d[k] += v * sqrt_exp;
            n[k] += v * v;
            s[k] += v.sqrt();
            r[k] += v;
        }
    }
    for ((((c, d), n), s), r) in c_tail
        .iter()
        .zip(d_tail)
        .zip(n_tail)
        .zip(s_tail)
        .zip(r_tail)
    {
        let v = c.max(0.0);
        *d += v * sqrt_exp;
        *n += v * v;
        *s += v.sqrt();
        *r += v;
    }
}

#[inline(never)]
pub fn pass1_chunked4(
    chrom: &[f32],
    sqrt_exp: f32,
    ms2_dot_prod: &mut [f32],
    ms2_norm_sq_obs: &mut [f32],
    sqrt_sum: &mut [f32],
    raw_sum: &mut [f32],
) {
    pass1_chunked_n::<4>(
        chrom,
        sqrt_exp,
        ms2_dot_prod,
        ms2_norm_sq_obs,
        sqrt_sum,
        raw_sum,
    );
}

#[inline(never)]
pub fn pass1_chunked8(
    chrom: &[f32],
    sqrt_exp: f32,
    ms2_dot_prod: &mut [f32],
    ms2_norm_sq_obs: &mut [f32],
    sqrt_sum: &mut [f32],
    raw_sum: &mut [f32],
) {
    pass1_chunked_n::<8>(
        chrom,
        sqrt_exp,
        ms2_dot_prod,
        ms2_norm_sq_obs,
        sqrt_sum,
        raw_sum,
    );
}

#[inline(never)]
pub fn pass1_chunked16(
    chrom: &[f32],
    sqrt_exp: f32,
    ms2_dot_prod: &mut [f32],
    ms2_norm_sq_obs: &mut [f32],
    sqrt_sum: &mut [f32],
    raw_sum: &mut [f32],
) {
    pass1_chunked_n::<16>(
        chrom,
        sqrt_exp,
        ms2_dot_prod,
        ms2_norm_sq_obs,
        sqrt_sum,
        raw_sum,
    );
}

#[inline(never)]
pub fn pass1_chunked32(
    chrom: &[f32],
    sqrt_exp: f32,
    ms2_dot_prod: &mut [f32],
    ms2_norm_sq_obs: &mut [f32],
    sqrt_sum: &mut [f32],
    raw_sum: &mut [f32],
) {
    pass1_chunked_n::<32>(
        chrom,
        sqrt_exp,
        ms2_dot_prod,
        ms2_norm_sq_obs,
        sqrt_sum,
        raw_sum,
    );
}

// ---------------------------------------------------------------------------
// Bench 2 — compute_main_score_trace Pass 1 elementwise multiply.
// ---------------------------------------------------------------------------

/// Scalar form matching the pre-C8 body.
#[inline(never)]
pub fn cms_pass1_scalar(a: &[f32], b: &[f32], out: &mut [f32]) {
    for i in 0..out.len() {
        out[i] = a[i] * b[i];
    }
}

/// Chunked form.
#[inline(never)]
pub fn cms_pass1_chunked(a: &[f32], b: &[f32], out: &mut [f32]) {
    let (a_ch, a_tail) = a.as_chunks::<8>();
    let (b_ch, b_tail) = b.as_chunks::<8>();
    let (o_ch, o_tail) = out.as_chunks_mut::<8>();
    for ((a, b), o) in a_ch.iter().zip(b_ch).zip(o_ch) {
        for k in 0..8 {
            o[k] = a[k] * b[k];
        }
    }
    for ((a, b), o) in a_tail.iter().zip(b_tail).zip(o_tail) {
        *o = *a * *b;
    }
}

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

fn xorshift32(state: &mut u32) -> u32 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    x
}

fn rand_f32(state: &mut u32) -> f32 {
    (xorshift32(state) >> 9) as f32 / (1 << 23) as f32
}

struct Case {
    chrom: Vec<f32>,
    dot_prod: Vec<f32>,
    norm_sq: Vec<f32>,
    sqrt_sum: Vec<f32>,
    raw_sum: Vec<f32>,
    sqrt_exp: f32,
}

fn build_cases(n: usize, len: usize, seed: u32) -> Vec<Case> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            let mut chrom = vec![0.0f32; len];
            // ~60% of cycles carry signal, 40% zero, to exercise the gated branch.
            for v in chrom.iter_mut() {
                let r = rand_f32(&mut state);
                *v = if r < 0.4 { 0.0 } else { r * 1e5 };
            }
            Case {
                chrom,
                dot_prod: vec![0.0f32; len],
                norm_sq: vec![0.0f32; len],
                sqrt_sum: vec![0.0f32; len],
                raw_sum: vec![0.0f32; len],
                sqrt_exp: rand_f32(&mut state),
            }
        })
        .collect()
}

/// Run scalar + each chunked variant on the same data and panic on any
/// divergence. Uses chunks that don't divide the length to exercise the
/// remainder tail path.
fn check_parity() {
    let odd_len = 53;
    let mut state = 0xC0FFEEu32;
    let chrom: Vec<f32> = (0..odd_len)
        .map(|_| {
            let r = rand_f32(&mut state);
            if r < 0.4 { 0.0 } else { r * 1e5 }
        })
        .collect();
    let sqrt_exp = rand_f32(&mut state);

    let mut ref_d = vec![0.0f32; odd_len];
    let mut ref_n = vec![0.0f32; odd_len];
    let mut ref_s = vec![0.0f32; odd_len];
    let mut ref_r = vec![0.0f32; odd_len];
    pass1_scalar(
        &chrom, sqrt_exp, &mut ref_d, &mut ref_n, &mut ref_s, &mut ref_r,
    );

    for (name, f) in [
        (
            "pass1_chunked4",
            pass1_chunked4 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
        (
            "pass1_chunked8",
            pass1_chunked8 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
        (
            "pass1_chunked16",
            pass1_chunked16 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
        (
            "pass1_chunked32",
            pass1_chunked32 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
    ] {
        let mut d = vec![0.0f32; odd_len];
        let mut n = vec![0.0f32; odd_len];
        let mut s = vec![0.0f32; odd_len];
        let mut r = vec![0.0f32; odd_len];
        f(&chrom, sqrt_exp, &mut d, &mut n, &mut s, &mut r);
        for i in 0..odd_len {
            let tol = 1e-4;
            assert!(
                (d[i] - ref_d[i]).abs() <= tol * ref_d[i].abs().max(1.0),
                "{name}: dot_prod[{i}] {} vs ref {}",
                d[i],
                ref_d[i],
            );
            assert!(
                (n[i] - ref_n[i]).abs() <= tol * ref_n[i].abs().max(1.0),
                "{name}: norm_sq[{i}] {} vs ref {}",
                n[i],
                ref_n[i],
            );
            assert!(
                (s[i] - ref_s[i]).abs() <= tol * ref_s[i].abs().max(1.0),
                "{name}: sqrt_sum[{i}] {} vs ref {}",
                s[i],
                ref_s[i],
            );
            assert!(
                (r[i] - ref_r[i]).abs() <= tol * ref_r[i].abs().max(1.0),
                "{name}: raw_sum[{i}] {} vs ref {}",
                r[i],
                ref_r[i],
            );
        }
    }
    eprintln!("parity OK: scalar == chunked{{4,8,16,32}} at len={odd_len}");
}

fn main() {
    check_parity();

    // --- compute_pass_1 inner ---
    eprintln!(
        "=== pass1 inner: {} cases × {} cycles × {} iters ===",
        N_CASES, N_CYCLES, N_ITERS
    );
    let mut cases = build_cases(N_CASES, N_CYCLES, 0xdead_beef);
    {
        // scalar
        for c in cases.iter_mut() {
            c.dot_prod.iter_mut().for_each(|v| *v = 0.0);
            c.norm_sq.iter_mut().for_each(|v| *v = 0.0);
            c.sqrt_sum.iter_mut().for_each(|v| *v = 0.0);
            c.raw_sum.iter_mut().for_each(|v| *v = 0.0);
        }
        let start = Instant::now();
        for _ in 0..N_ITERS {
            for c in cases.iter_mut() {
                pass1_scalar(
                    black_box(&c.chrom),
                    black_box(c.sqrt_exp),
                    &mut c.dot_prod,
                    &mut c.norm_sq,
                    &mut c.sqrt_sum,
                    &mut c.raw_sum,
                );
                black_box(&c.dot_prod);
            }
        }
        let e = start.elapsed();
        let total = N_ITERS * N_CASES * N_CYCLES;
        eprintln!(
            "pass1_scalar         {:?}  ({:.2} ns/elem)",
            e,
            e.as_nanos() as f64 / total as f64
        );
    }
    for (label, f) in &[
        (
            "pass1_chunked4",
            pass1_chunked4 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
        (
            "pass1_chunked8",
            pass1_chunked8 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
        (
            "pass1_chunked16",
            pass1_chunked16 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
        (
            "pass1_chunked32",
            pass1_chunked32 as fn(&[f32], f32, &mut [f32], &mut [f32], &mut [f32], &mut [f32]),
        ),
    ] {
        for c in cases.iter_mut() {
            c.dot_prod.iter_mut().for_each(|v| *v = 0.0);
            c.norm_sq.iter_mut().for_each(|v| *v = 0.0);
            c.sqrt_sum.iter_mut().for_each(|v| *v = 0.0);
            c.raw_sum.iter_mut().for_each(|v| *v = 0.0);
        }
        let start = Instant::now();
        for _ in 0..N_ITERS {
            for c in cases.iter_mut() {
                f(
                    black_box(&c.chrom),
                    black_box(c.sqrt_exp),
                    &mut c.dot_prod,
                    &mut c.norm_sq,
                    &mut c.sqrt_sum,
                    &mut c.raw_sum,
                );
                black_box(&c.dot_prod);
            }
        }
        let e = start.elapsed();
        let total = N_ITERS * N_CASES * N_CYCLES;
        eprintln!(
            "{label:<20} {:?}  ({:.2} ns/elem)",
            e,
            e.as_nanos() as f64 / total as f64
        );
    }

    // --- compute_main_score_trace Pass 1 (elementwise mul) ---
    eprintln!(
        "\n=== cms_pass1: {} cases × {} cycles × {} iters ===",
        N_CASES, N_CYCLES, N_ITERS
    );
    let mut state = 0x1234_5678u32;
    let cms_cases: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = (0..N_CASES)
        .map(|_| {
            let a: Vec<f32> = (0..N_CYCLES).map(|_| rand_f32(&mut state)).collect();
            let b: Vec<f32> = (0..N_CYCLES).map(|_| rand_f32(&mut state)).collect();
            let o: Vec<f32> = vec![0.0; N_CYCLES];
            (a, b, o)
        })
        .collect();
    {
        let mut cc = cms_cases.clone();
        let start = Instant::now();
        for _ in 0..N_ITERS {
            for (a, b, o) in cc.iter_mut() {
                cms_pass1_scalar(black_box(a), black_box(b), o);
                black_box(&o);
            }
        }
        let e = start.elapsed();
        let total = N_ITERS * N_CASES * N_CYCLES;
        eprintln!(
            "cms_pass1_scalar     {:?}  ({:.2} ns/elem)",
            e,
            e.as_nanos() as f64 / total as f64
        );
    }
    {
        let mut cc = cms_cases.clone();
        let start = Instant::now();
        for _ in 0..N_ITERS {
            for (a, b, o) in cc.iter_mut() {
                cms_pass1_chunked(black_box(a), black_box(b), o);
                black_box(&o);
            }
        }
        let e = start.elapsed();
        let total = N_ITERS * N_CASES * N_CYCLES;
        eprintln!(
            "cms_pass1_chunked    {:?}  ({:.2} ns/elem)",
            e,
            e.as_nanos() as f64 / total as f64
        );
    }
}
