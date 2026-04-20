//! Compile-time toggle between rayon-parallel and plain-serial fold-reduce.
//!
//! `fold_reduce` is the only exit point. Both variants:
//! - call `init()` once per rayon worker (parallel) or exactly once (serial),
//! - feed the returned `Acc` through `fold(acc, (idx, &item))` for every item,
//! - merge per-worker accumulators via `reduce(a, b)` (parallel) or skip it (serial, single acc).
//!
//! The closures are identical across feature flags, so scoring call sites need
//! only one copy of the logic.

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
pub(super) fn fold_reduce<Item, Acc, Init, Fold, Reduce>(
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
            .fold(&init, |acc, pair| fold(acc, pair))
            .reduce(&init, |a, b| reduce(a, b))
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
}
