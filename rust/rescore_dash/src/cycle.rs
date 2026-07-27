//! Wrapping movement through a fixed list of variants.
//!
//! Four enums here are keyboard-cycled and stored by position in their own
//! `ALL` array. Each used to spell the same `position` + modular step by hand.

/// Position of `cur` in `all`, or `0` if it is somehow absent.
///
/// Absent is unreachable for the `ALL` arrays this is used with, and falling
/// back to `0` keeps callers total: a panic here would take down a redraw.
pub(crate) fn index_of<T: Copy + PartialEq>(all: &[T], cur: T) -> usize {
    all.iter().position(|x| *x == cur).unwrap_or(0)
}

/// `cur` moved `delta` positions through `all`, wrapping in either direction.
pub(crate) fn step<T: Copy + PartialEq>(all: &[T], cur: T, delta: isize) -> T {
    let i = index_of(all, cur) as isize;
    all[(i + delta).rem_euclid(all.len() as isize) as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_wraps_both_ways_and_round_trips() {
        let all = ['a', 'b', 'c'];
        assert_eq!(step(&all, 'a', 1), 'b');
        assert_eq!(step(&all, 'c', 1), 'a', "forward wraps");
        assert_eq!(step(&all, 'a', -1), 'c', "backward wraps");
        for &c in &all {
            assert_eq!(step(&all, step(&all, c, 1), -1), c);
        }
    }

    #[test]
    fn an_absent_item_restarts_from_the_front() {
        assert_eq!(index_of(&['a', 'b'], 'z'), 0);
        assert_eq!(step(&['a', 'b'], 'z', 1), 'b');
    }
}
