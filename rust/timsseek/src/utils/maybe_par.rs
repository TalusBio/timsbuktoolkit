//! Compile-time toggle between rayon-parallel and plain-serial data-parallel
//! shapes.
//!
//! Crate-wide concurrency utility: every `#[cfg(feature = "rayon")]` /
//! `#[cfg(not(...))]` pair that would otherwise be duplicated at a call site
//! lives here instead. The closures are identical across feature flags, so
//! call sites need only one copy of the logic, and the serial halves — which
//! only CI's `--no-default-features` leg ever compiles — cannot bit-rot
//! independently per call site.
//!
//! Three shapes, one per access pattern:
//!
//! * [`fold_reduce`] — per-ITEM fold over a slice with per-worker accumulator
//!   init. Reduction order depends on the worker count, so the accumulator
//!   must be order-insensitive (or the caller must not care).
//! * [`chunked_fold_reduce`] — per-CHUNK accumulate over a flat slice, then a
//!   sequential merge of the partials **in chunk order**. Because the chunk
//!   length is a caller-supplied constant, both the number of partials and
//!   the merge order are independent of the rayon thread pool: the result is
//!   bitwise-identical across runs and thread counts even for float sums.
//! * [`scatter_write`] — parallel scatter-write of `out[i] = f(i)`. No
//!   reduction at all, so trivially deterministic.

#[cfg(feature = "rayon")]
use rayon::prelude::*;

/// Execute a per-item fold over `items` with deterministic per-worker init.
///
/// * `init` runs once per rayon worker (parallel) or exactly once (serial).
/// * `fold` receives `(chunk_idx, &item)` — ignore the index with `|(_, item)|`
///   if the caller does not need it. Keeping it uniform lets both scoring
///   batch fns share one helper.
/// * `reduce` merges accumulators from different workers. In serial mode the
///   closure is never invoked but is still type-checked.
/// * `init()` must produce an identity element for `reduce`: for all `x`,
///   `reduce(init(), x) == x` and `reduce(x, init()) == x`. Rayon calls
///   `init` to seed BOTH the fold identity and the reduce identity, so a
///   non-identity init can silently miscount on empty / tiny inputs. The
///   serial branch ignores this, but writing closures that honor it keeps
///   parallel semantics correct.
///
/// The number of accumulators — and therefore the order in which `reduce`
/// combines them — follows the rayon thread pool. Use
/// [`chunked_fold_reduce`] when that order has to be reproducible.
pub(crate) fn fold_reduce<Item, Acc, Init, Fold, Reduce>(
    items: &[Item],
    init: Init,
    fold: Fold,
    reduce: Reduce,
) -> Acc
where
    // `Item: Sync` is required by rayon's par_iter; the serial branch does not
    // need it, but keeping the bound unified means callers don't have to think
    // about feature-flag-dependent trait bounds.
    Item: Sync,
    Acc: Send,
    Init: Fn() -> Acc + Send + Sync,
    Fold: Fn(Acc, (usize, &Item)) -> Acc + Send + Sync,
    Reduce: Fn(Acc, Acc) -> Acc + Send + Sync,
{
    #[cfg(feature = "rayon")]
    {
        items
            .par_iter()
            .enumerate()
            .fold(&init, fold)
            .reduce(&init, reduce)
    }
    #[cfg(not(feature = "rayon"))]
    {
        let _ = &reduce; // suppress unused-param lint in serial builds
        items
            .iter()
            .enumerate()
            .fold(init(), |acc, pair| fold(acc, pair))
    }
}

/// Split `items` into fixed `chunk_len` blocks, accumulate each block
/// independently, then merge the partials sequentially in chunk order.
///
/// * `zero` builds one fresh accumulator per CHUNK (not per worker), and must
///   be the identity for `merge`.
/// * `fold` receives `(&mut acc, chunk_idx, block)`. `chunk_idx` is the
///   position of `block` in `items.chunks(chunk_len)`, so callers that need a
///   global element index can recover it as `chunk_idx * chunk_len + local`.
///   The final block is short whenever `chunk_len` does not divide
///   `items.len()`.
/// * `merge` folds a partial into the running total, always left-to-right over
///   ascending `chunk_idx`, starting from `zero()`.
///
/// **Determinism**: the partition depends only on `items.len()` and
/// `chunk_len`, and the merge walks partials in index order, so the exact
/// sequence of float additions is the same on 1 thread and on 64. This is the
/// shape to use for reductions that must be bitwise-reproducible; picking
/// `chunk_len` per thread count would forfeit that.
///
/// Panics if `chunk_len == 0`.
pub(crate) fn chunked_fold_reduce<Item, Acc, Zero, Fold, Merge>(
    items: &[Item],
    chunk_len: usize,
    zero: Zero,
    fold: Fold,
    merge: Merge,
) -> Acc
where
    Item: Sync,
    Acc: Send,
    Zero: Fn() -> Acc + Send + Sync,
    Fold: Fn(&mut Acc, usize, &[Item]) + Send + Sync,
    Merge: Fn(&mut Acc, Acc),
{
    assert!(chunk_len > 0, "chunk_len must be non-zero");
    let per_chunk = |(ci, block): (usize, &[Item])| {
        let mut acc = zero();
        fold(&mut acc, ci, block);
        acc
    };

    #[cfg(feature = "rayon")]
    let partials: Vec<Acc> = items
        .par_chunks(chunk_len)
        .enumerate()
        .map(per_chunk)
        .collect();
    #[cfg(not(feature = "rayon"))]
    let partials: Vec<Acc> = items.chunks(chunk_len).enumerate().map(per_chunk).collect();

    let mut out = zero();
    for p in partials {
        merge(&mut out, p);
    }
    out
}

/// Fill `out` in parallel with `out[i] = f(i)`.
///
/// Each index is written exactly once by exactly one worker, so there is no
/// reduction and no ordering concern; the result matches the serial loop
/// element for element.
pub(crate) fn scatter_write<T, F>(out: &mut [T], f: F)
where
    T: Send,
    F: Fn(usize) -> T + Send + Sync,
{
    #[cfg(feature = "rayon")]
    {
        out.par_iter_mut().enumerate().for_each(|(i, o)| *o = f(i));
    }
    #[cfg(not(feature = "rayon"))]
    {
        for (i, o) in out.iter_mut().enumerate() {
            *o = f(i);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sum_matches_plain_fold() {
        let items: Vec<u32> = (0..1_000).collect();
        let got = fold_reduce(
            &items,
            || 0u64,
            |acc, (idx, v)| acc + *v as u64 + idx as u64,
            |a, b| a + b,
        );
        let expected: u64 = items
            .iter()
            .enumerate()
            .map(|(idx, v)| *v as u64 + idx as u64)
            .sum();
        assert_eq!(got, expected);
    }

    #[test]
    fn empty_slice_returns_init() {
        // init must be an identity element for reduce so that both serial
        // (init called once) and rayon (init may be called >1 for identity)
        // agree on the result.
        let items: Vec<u32> = Vec::new();
        let got = fold_reduce(&items, || 0u64, |acc, (_, v)| acc + *v as u64, |a, b| a + b);
        assert_eq!(got, 0);
    }

    #[test]
    fn index_is_passed_through() {
        let items: Vec<char> = vec!['a', 'b', 'c', 'd'];
        let mut got = fold_reduce(
            &items,
            Vec::<usize>::new,
            |mut acc, (idx, _)| {
                acc.push(idx);
                acc
            },
            |mut a, b| {
                a.extend(b);
                a
            },
        );
        got.sort_unstable();
        assert_eq!(got, vec![0, 1, 2, 3]);
    }

    /// The partition is a pure function of `len` and `chunk_len`, and the merge
    /// walks it in order: the chunk index handed to `fold` must reconstruct the
    /// exact `chunks()` split, trailing short block included.
    #[test]
    fn chunked_blocks_are_index_addressable() {
        let items: Vec<u32> = (0..10).collect();
        let seen = chunked_fold_reduce(
            &items,
            4,
            Vec::<(usize, Vec<u32>)>::new,
            |acc, ci, block| acc.push((ci, block.to_vec())),
            |a, b| a.extend(b),
        );
        assert_eq!(
            seen,
            vec![
                (0, vec![0, 1, 2, 3]),
                (1, vec![4, 5, 6, 7]),
                (2, vec![8, 9]),
            ]
        );
    }

    /// Float addition is not associative, so a reduction whose partial count or
    /// merge order tracked the thread pool would drift run to run. Fixing
    /// `chunk_len` pins both: the result must equal the hand-written
    /// chunk-then-merge sum bit for bit.
    #[test]
    fn chunked_sum_is_bitwise_stable() {
        let items: Vec<f64> = (0..5_000).map(|i| 1.0 / (i as f64 + 1.0)).collect();
        let got = chunked_fold_reduce(
            &items,
            64,
            || 0.0f64,
            |acc, _, block| {
                for v in block {
                    *acc += v;
                }
            },
            |a, b| *a += b,
        );
        let expected = items.chunks(64).fold(0.0f64, |acc, block| {
            acc + block.iter().fold(0.0f64, |a, v| a + v)
        });
        assert_eq!(got.to_bits(), expected.to_bits());
    }

    #[test]
    fn scatter_write_fills_every_index() {
        let mut out = vec![0usize; 257];
        scatter_write(&mut out, |i| i * 3);
        assert!(out.iter().enumerate().all(|(i, v)| *v == i * 3));
    }
}
