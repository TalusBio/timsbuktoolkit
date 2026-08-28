//! Inline-capacity bench for `ExpectedIntensities`.
//!
//! Answers two questions:
//!
//! 1. What did switching the scoring hot path from `self.expected = try_from_pairs(..)`
//!    to `refill_from_pairs` buy?
//! 2. Is `TinyVec`'s inline storage earning its keep, and at what capacity?
//!
//! The fragment-count distribution is taken from the DIA-NN HeLa `.speclib`
//! fixture (n=3192 rows): min 6, median 12, p95 12, max 32, mean 11.5, with
//! 1.1% of rows above the current inline capacity of 13. The precursor
//! envelope is always exactly 3 entries, because `isotope_dist_or_averagine`
//! returns `[f32; 3]` and every `IsotopeStrategy::FromComposition` site
//! hardcodes `n_isotopes: 3`.
//!
//! Run:
//!   cargo run -p timsseek --release --example expected_intensities_bench
//!   cargo run -p timsseek --release --example expected_intensities_bench -- 20000

use std::hint::black_box;
use std::time::Instant;

use micromzpaf::IonAnnot;
use timsquery::tinyvec::{
    Array,
    TinyVec,
};

/// Fragment counts matching the fixture's measured percentiles.
fn fragment_counts(n: usize) -> Vec<usize> {
    (0..n)
        .map(|i| match i % 1000 {
            0..=10 => 14 + (i % 19), // 1.1% spill past 13, up to 32
            11..=110 => 6 + (i % 6), // left tail down to the observed min of 6
            _ => 12,                 // the mode, and p25 through p95
        })
        .collect()
}

fn make_items(counts: &[usize]) -> Vec<Vec<(IonAnnot, f32)>> {
    counts
        .iter()
        .map(|&c| {
            (1..=c)
                .map(|i| {
                    let ord = u8::try_from(i).unwrap();
                    (
                        IonAnnot::try_new('y', Some(ord), 1, 0).unwrap(),
                        i as f32 * 0.1,
                    )
                })
                .collect()
        })
        .collect()
}

const PRECS: [(i8, f32); 3] = [(0, 1.0), (1, 0.5), (2, 0.2)];

/// The uniqueness check `try_from_pairs`/`refill_from_pairs` runs: O(n^2) over
/// a linear scan, which is the same shape `linear_get` uses for lookup.
fn push_checked<V: Extend<(IonAnnot, f32)> + AsRef<[(IonAnnot, f32)]>>(
    dst: &mut V,
    src: &[(IonAnnot, f32)],
) {
    for &(k, v) in src {
        assert!(
            !dst.as_ref().iter().any(|(kk, _)| *kk == k),
            "bench inputs are unique"
        );
        dst.extend(std::iter::once((k, v)));
    }
}

/// Refill into a reused buffer: clear, then push with the uniqueness check.
fn bench_refill<A>(items: &[Vec<(IonAnnot, f32)>], reps: usize) -> f64
where
    A: Array<Item = (IonAnnot, f32)>,
{
    let mut frags: TinyVec<A> = TinyVec::new();
    let mut precs: TinyVec<[(i8, f32); 3]> = TinyVec::new();
    let t = Instant::now();
    for _ in 0..reps {
        for item in items {
            frags.clear();
            push_checked(&mut frags, item);
            precs.clear();
            precs.extend(PRECS);
            black_box((&frags, &precs));
        }
    }
    ns_per_item(t, items.len() * reps)
}

/// Build a fresh value per item and move it, which is what the hot path did
/// before `refill_from_pairs`.
fn bench_assign<A>(items: &[Vec<(IonAnnot, f32)>], reps: usize) -> f64
where
    A: Array<Item = (IonAnnot, f32)>,
{
    let t = Instant::now();
    for _ in 0..reps {
        for item in items {
            let mut frags: TinyVec<A> = TinyVec::new();
            push_checked(&mut frags, item);
            let mut precs: TinyVec<[(i8, f32); 3]> = TinyVec::new();
            precs.extend(PRECS);
            black_box((frags, precs));
        }
    }
    ns_per_item(t, items.len() * reps)
}

/// Capacity 0: a plain `Vec`, reused the same way `refill` reuses a `TinyVec`.
fn bench_refill_vec(items: &[Vec<(IonAnnot, f32)>], reps: usize) -> f64 {
    let mut frags: Vec<(IonAnnot, f32)> = Vec::new();
    let mut precs: Vec<(i8, f32)> = Vec::new();
    let t = Instant::now();
    for _ in 0..reps {
        for item in items {
            frags.clear();
            push_checked(&mut frags, item);
            precs.clear();
            precs.extend(PRECS);
            black_box((&frags, &precs));
        }
    }
    ns_per_item(t, items.len() * reps)
}

/// A fresh `Vec` per item: no reuse and no inline storage, the worst case.
fn bench_assign_vec(items: &[Vec<(IonAnnot, f32)>], reps: usize) -> f64 {
    let t = Instant::now();
    for _ in 0..reps {
        for item in items {
            let mut frags: Vec<(IonAnnot, f32)> = Vec::new();
            push_checked(&mut frags, item);
            let mut precs: Vec<(i8, f32)> = Vec::new();
            precs.extend(PRECS);
            black_box((frags, precs));
        }
    }
    ns_per_item(t, items.len() * reps)
}

/// A named variant plus the closure that times it.
type Variant<'a> = (&'static str, Box<dyn FnMut() -> f64 + 'a>);

fn ns_per_item(t: Instant, n: usize) -> f64 {
    t.elapsed().as_nanos() as f64 / n as f64
}

/// Min of `trials`, interleaved across variants so machine drift hits every
/// variant equally. Min rather than mean: the fast path is the real cost and
/// everything above it is interference.
fn best_of(trials: usize, mut f: impl FnMut() -> f64) -> f64 {
    (0..trials).map(|_| f()).fold(f64::INFINITY, f64::min)
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(10_000);
    let reps = 10;
    let trials = 9;

    let counts = fragment_counts(n);
    let items = make_items(&counts);
    let spilled = counts.iter().filter(|c| **c > 13).count();
    println!(
        "n={} items x {} reps, best of {} trials, {:.1}% above inline capacity 13\n",
        n,
        reps,
        trials,
        100.0 * spilled as f64 / n as f64
    );

    // Warm up: the first pass pays page faults the timed runs should not see.
    black_box(bench_refill::<[(IonAnnot, f32); 13]>(&items, 2));

    let mut variants: Vec<Variant<'_>> = vec![
        (
            "refill, TinyVec inline 13",
            Box::new(|| bench_refill::<[(IonAnnot, f32); 13]>(&items, reps)),
        ),
        (
            "refill, TinyVec inline 8",
            Box::new(|| bench_refill::<[(IonAnnot, f32); 8]>(&items, reps)),
        ),
        (
            "refill, TinyVec inline 16",
            Box::new(|| bench_refill::<[(IonAnnot, f32); 16]>(&items, reps)),
        ),
        (
            "refill, TinyVec inline 32",
            Box::new(|| bench_refill::<[(IonAnnot, f32); 32]>(&items, reps)),
        ),
        (
            "refill, Vec (current)",
            Box::new(|| bench_refill_vec(&items, reps)),
        ),
        (
            "assign, TinyVec inline 13 (old hot path)",
            Box::new(|| bench_assign::<[(IonAnnot, f32); 13]>(&items, reps)),
        ),
        (
            "assign, Vec (inline 0)",
            Box::new(|| bench_assign_vec(&items, reps)),
        ),
    ];

    let results: Vec<(&str, f64)> = variants
        .iter_mut()
        .map(|(label, f)| (*label, best_of(trials, f)))
        .collect();

    println!("=== ns per item ===");
    println!("  {:<42} {:>9}", "variant", "ns/item");
    for (label, ns) in &results {
        println!("  {label:<42} {ns:>9.1}");
    }

    println!("\n=== sizes (bytes) ===");
    for (label, sz) in [
        (
            "ExpectedIntensities<IonAnnot>",
            size_of::<timsseek::ExpectedIntensities<IonAnnot>>(),
        ),
        (
            "fragment TinyVec<13>",
            size_of::<TinyVec<[(IonAnnot, f32); 13]>>(),
        ),
        (
            "precursor TinyVec<13> (old)",
            size_of::<TinyVec<[(i8, f32); 13]>>(),
        ),
        (
            "precursor TinyVec<3> (measured max)",
            size_of::<TinyVec<[(i8, f32); 3]>>(),
        ),
    ] {
        println!("  {label:<42} {sz:>9}");
    }
}
