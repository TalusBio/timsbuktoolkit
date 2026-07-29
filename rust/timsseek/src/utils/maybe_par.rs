//! Rayon-or-serial data-parallel shapes, so the `#[cfg(feature = "rayon")]`
//! pairs live here instead of at every call site.
//!
//! The shapes, and what each one promises about ordering:
//!
//! * [`fold_reduce`] — per-item fold with per-worker accumulators. The
//!   accumulator count, and so the reduce order, follows the thread pool: NOT
//!   reproducible for a non-associative `reduce` (e.g. float addition).
//! * [`chunked_fold_reduce`] — per-chunk accumulate over a thread-count-
//!   independent partition, merged in ascending chunk order: bitwise
//!   reproducible.
//! * [`sort_by`] — in-place STABLE sort: ties keep input order, on any thread
//!   count. Reproducible, which matters when a later pass walks the order
//!   positionally (see `assign_qval`).

#[cfg(feature = "rayon")]
use rayon::prelude::*;

/// Per-item fold with per-worker accumulator init.
///
/// `init()` must be an identity for `reduce`: rayon uses it to seed both the
/// fold and the reduce, so a non-identity init miscounts on small inputs.
///
/// The accumulator count — and therefore the order `reduce` combines them —
/// follows the thread pool. Use [`chunked_fold_reduce`] when that order must
/// be reproducible.
pub(crate) fn fold_reduce<Item, Acc, Init, Fold, Reduce>(
    items: &[Item],
    init: Init,
    fold: Fold,
    reduce: Reduce,
) -> Acc
where
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

/// Per-chunk accumulate over fixed `chunk_len` blocks, then a sequential merge
/// of the partials in ascending chunk order.
///
/// `zero()` runs once per chunk and must be an identity for `merge`. `fold`
/// gets `(&mut acc, chunk_idx, block)` where `chunk_idx` is the position in
/// `items.chunks(chunk_len)`, so a global element index is
/// `chunk_idx * chunk_len + local`; the last block is short when `chunk_len`
/// does not divide `items.len()`.
///
/// Determinism: the partition depends only on `items.len()` and `chunk_len`,
/// and the merge walks partials in index order, so the float addition sequence
/// is identical on 1 thread and on 64. Unlike [`fold_reduce`], this is
/// bitwise-reproducible.
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

/// Sort `items` in place by `cmp`, in parallel when rayon is on.
///
/// Stable, so ties keep input order on any thread count — `assign_qval` walks
/// the sorted order positionally, so that reproducibility matters. Needs
/// scratch proportional to `items.len() * size_of::<T>()`.
pub(crate) fn sort_by<T, F>(items: &mut [T], cmp: F)
where
    T: Send,
    F: Fn(&T, &T) -> std::cmp::Ordering + Send + Sync,
{
    #[cfg(feature = "rayon")]
    {
        items.par_sort_by(cmp);
    }
    #[cfg(not(feature = "rayon"))]
    {
        items.sort_by(cmp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_descending_and_is_stable() {
        // 5 keys over 1000 elements, so every key is a long tie run.
        let mut items: Vec<(i32, usize)> = (0..1_000).map(|i| ((i % 5) as i32, i)).collect();
        sort_by(&mut items, |a, b| b.0.cmp(&a.0));
        assert!(items.windows(2).all(|w| w[0].0 >= w[1].0), "not descending");
        for run in items.chunk_by(|a, b| a.0 == b.0) {
            assert!(run.windows(2).all(|w| w[0].1 < w[1].1), "tie run reordered");
        }
    }

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

    /// Callers recover global element indices as `chunk_idx * chunk_len + local`
    /// (see `lda.rs` label lookup), so `chunk_idx` must track `chunks()` exactly.
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

    /// Float addition is not associative: a thread-count-dependent merge order
    /// would drift run to run.
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
}
