//! Rayon-or-serial data-parallel shapes, so the `#[cfg(feature = "rayon")]`
//! pairs live here instead of at every call site.
//!
//! The shapes, and what each one promises about ordering:
//!
//! * [`BlockWorkers`] -- scoped jobs over a thread-count-dependent partition,
//!   followed by an ordered merge.
//! * [`chunked_fold_reduce`] -- per-chunk accumulate over a thread-count-
//!   independent partition, merged in ascending chunk order: bitwise
//!   reproducible.
//! * [`sort_by`] -- in-place STABLE sort: ties keep input order, on any thread
//!   count. Reproducible, which matters when a later pass walks the order
//!   positionally (see `assign_qval`).

#[cfg(feature = "rayon")]
use rayon::prelude::*;

use std::sync::Mutex;

/// Blocks per pool thread and batch. Blocks are equal in item count, not in
/// cost, and the batch ends when its slowest block does, so finer blocks keep
/// that tail short on mixed performance/efficiency cores (4 per thread left
/// Phase 1 ~10% slower than the recursive split on an M2).
const CHUNKS_PER_THREAD: usize = 16;

/// Reuses one scratch slot per pool thread across batches of scoped jobs.
/// Jobs run to completion and take their slot once per block.
///
/// Each block is a flat FIFO job, not a node of a recursive split, so a
/// worker's stack depth does not grow with the thread count. Recursive
/// `fold`/`reduce` nested join frames per steal and overflowed 2 MiB stacks
/// on 48 workers.
pub(crate) struct BlockWorkers<State> {
    slots: Vec<Mutex<Option<State>>>,
}

impl<State: Send> BlockWorkers<State> {
    pub(crate) fn new() -> Self {
        Self {
            slots: Vec::with_capacity(1),
        }
    }

    fn chunk_len(&mut self, item_count: usize) -> usize {
        #[cfg(feature = "rayon")]
        let thread_count = rayon::current_num_threads();
        #[cfg(not(feature = "rayon"))]
        let thread_count = 1;
        if thread_count > self.slots.len() {
            self.slots.resize_with(thread_count, || Mutex::new(None));
        }
        if thread_count == 1 {
            item_count.max(1)
        } else {
            item_count.div_ceil(thread_count * CHUNKS_PER_THREAD).max(1)
        }
    }

    fn run<Block: Send>(
        &self,
        blocks: impl IntoIterator<Item = Block>,
        init: impl Fn() -> State + Sync,
        work: impl Fn(&mut State, Block) + Sync,
    ) {
        let run_block = |index: usize, block| {
            // A block job never waits on Rayon, so it never runs nested in
            // another job on the same thread.
            let mut slot = self.slots[index]
                .try_lock()
                .expect("one block job per thread at a time");
            work(slot.get_or_insert_with(&init), block);
        };
        #[cfg(feature = "rayon")]
        if rayon::current_num_threads() > 1 {
            rayon::in_place_scope_fifo(|scope| {
                let run_block = &run_block;
                for block in blocks {
                    scope.spawn_fifo(move |_| {
                        let index = rayon::current_thread_index().expect("scoped Rayon job");
                        run_block(index, block);
                    });
                }
            });
            return;
        }
        for block in blocks {
            run_block(0, block);
        }
    }

    /// Fold blocks independently, then merge partials in input order.
    pub(crate) fn fold_reduce<Item: Sync, Acc: Send>(
        &mut self,
        items: &[Item],
        init: impl Fn() -> State + Sync,
        fold: impl Fn(&mut State, &[Item]) -> Acc + Sync,
        merge: impl Fn(&mut Acc, Acc),
    ) -> Acc {
        let chunk_len = self.chunk_len(items.len());
        let chunk_count = items.len().div_ceil(chunk_len);
        if chunk_count <= 1 {
            let slot = self.slots[0].get_mut().expect("scoring scratch poisoned");
            return fold(slot.get_or_insert_with(init), items);
        }
        let mut partials = Vec::with_capacity(chunk_count);
        partials.resize_with(chunk_count, || None);
        self.run(
            items.chunks(chunk_len).zip(&mut partials),
            init,
            |state, (chunk, partial)| *partial = Some(fold(state, chunk)),
        );
        let mut partials = partials
            .into_iter()
            .map(|partial| partial.expect("scope completed every chunk"));
        let mut result = partials.next().expect("nonempty input has a chunk");
        for partial in partials {
            merge(&mut result, partial);
        }
        result
    }

    pub(crate) fn map_into<Item: Sync, Output: Send>(
        &mut self,
        items: &[Item],
        output: &mut [Output],
        init: impl Fn() -> State + Sync,
        map: impl Fn(&mut State, &Item) -> Output + Sync,
    ) {
        assert_eq!(items.len(), output.len());
        let chunk_len = self.chunk_len(items.len());
        self.run(
            items.chunks(chunk_len).zip(output.chunks_mut(chunk_len)),
            init,
            |state, (items, output)| {
                for (item, destination) in items.iter().zip(output) {
                    *destination = map(state, item);
                }
            },
        );
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
/// is identical on 1 thread and on 64. Unlike [`BlockWorkers`], this is
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
/// Stable, so ties keep input order on any thread count -- `assign_qval` walks
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
    fn empty_slice_is_folded_once() {
        let got = BlockWorkers::new().fold_reduce(
            &[] as &[u32],
            || (),
            |(), chunk| chunk.len(),
            |_, _| unreachable!(),
        );
        assert_eq!(got, 0);
    }

    #[test]
    fn chunks_merge_in_input_order() {
        let items: Vec<u32> = (0..1_003).collect();
        let got = BlockWorkers::new().fold_reduce(
            &items,
            || (),
            |(), chunk| chunk.to_vec(),
            |result, partial| {
                #[cfg(feature = "rayon")]
                assert_eq!(
                    rayon::current_thread_index(),
                    None,
                    "merge runs on the caller"
                );
                result.extend(partial);
            },
        );
        assert_eq!(got, items);
    }

    /// Regression: the recursive `par_iter().fold().reduce()` this replaced
    /// overflows these 64 KiB stacks on 48 workers, because every split level
    /// keeps partial accumulators on the stack.
    #[cfg(feature = "rayon")]
    #[test]
    fn scoring_fold_uses_bounded_jobs_on_small_stacks() {
        const THREAD_COUNT: usize = 48;
        const ITEM_COUNT: usize = 10_003;
        // Partials move by value, like the scoring heaps and timings.
        const PAYLOAD_WORDS: usize = 512;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(THREAD_COUNT)
            .stack_size(64 * 1024)
            .build()
            .expect("build a small-stack test pool");
        let items: Vec<usize> = (0..ITEM_COUNT).collect();
        pool.install(|| {
            use std::sync::atomic::{
                AtomicUsize,
                Ordering,
            };
            let constructions = AtomicUsize::new(0);
            let mut workers = BlockWorkers::new();
            for _ in 0..3 {
                let (counts, sum, jobs) = workers.fold_reduce(
                    &items,
                    || {
                        constructions.fetch_add(1, Ordering::Relaxed);
                    },
                    |(), chunk| {
                        let mut counts = [0usize; PAYLOAD_WORDS];
                        counts[chunk.len() % PAYLOAD_WORDS] = 1;
                        (counts, chunk.iter().sum::<usize>(), 1)
                    },
                    |result, partial| {
                        for (total, count) in result.0.iter_mut().zip(partial.0) {
                            *total += count;
                        }
                        result.1 += partial.1;
                        result.2 += partial.2;
                    },
                );
                assert_eq!(sum, ITEM_COUNT * (ITEM_COUNT - 1) / 2);
                assert_eq!(counts.iter().sum::<usize>(), jobs);
                assert!(
                    (THREAD_COUNT..=CHUNKS_PER_THREAD * THREAD_COUNT).contains(&jobs),
                    "{jobs} jobs"
                );
            }
            assert!((1..=THREAD_COUNT).contains(&constructions.into_inner()));
        });
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
