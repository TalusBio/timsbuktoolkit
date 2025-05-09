use std::cmp::Ordering;

/// Expand the scan offset slice to mobilities.
///
/// The scan offsets is in essence a run-length
/// encoded vector of scan numbers that can be converter to the 1/k0
/// values.
///
/// Essentially ... the slice [0,4,5,5], would expand to
/// [0,0,0,0,1]; 0 to 4 have index 0, 4 to 5 have index 1, 5 to 5 would
/// have index 2 but its empty!
///
/// Then this index can be converted using the Scan2ImConverter.convert
///
/// ... This should problably be implemented and exposed in timsrust.
fn expand_mobility_iter(scan_offsets: &'_ [usize]) -> impl Iterator<Item = u16> + '_ {
    let ims_iter = scan_offsets
        .windows(2)
        .enumerate()
        .filter_map(|(i, w)| {
            let num = w[1] - w[0];
            if num == 0 {
                return None;
            }
            let lo = w[0];
            let hi = w[1];

            Some((i as u16, lo, hi))
        })
        .flat_map(|(im, lo, hi)| (lo..hi).map(move |_| im));
    ims_iter
}

/// Explodes the compressed representation of a vector to its
/// original representation.
///
/// # Examples
/// ```
/// use timsquery::utils::compress_explode::explode_vec;
///
/// let input = vec![0, 0, 5, 5, 5, 7];
/// let out = explode_vec(&input);
/// assert_eq!(out, vec![1, 1, 1, 1, 1, 4, 4]);
/// ```
///
/// This function is the inverse of `compress_vec`.
pub fn explode_vec(input: &[usize]) -> Vec<u16> {
    let last_val = match input.last() {
        Some(last) => *last,
        None => return Vec::new(),
    };

    let mut out = Vec::with_capacity(last_val);
    out.extend(expand_mobility_iter(input));
    out
}

/// Compresses a monotonically increasing vector of usize.
///
/// Returns a vector where each index represents the starting position of that value
/// in the original slice. Output length is max(input) + 2.
///
/// # Panics
///
/// Panics if the input is not monotonically increasing.
///
/// # Examples
///
/// ```
/// use timsquery::utils::compress_explode::compress_vec;
///
/// let input = vec![1, 1, 1, 1, 1, 4, 4];
/// assert_eq!(compress_vec(&input), vec![0, 0, 5, 5, 5, 7]);
///
/// let input = vec![0, 0, 2, 2, 5, 5, 5];
/// assert_eq!(compress_vec(&input), vec![0, 2, 2, 4, 4, 4, 7]);
///
/// assert_eq!(compress_vec(&[]), vec![]);
/// ```
///
/// This function is the inverse of `explode_vec`.
pub fn compress_vec(input: &[u16]) -> Vec<usize> {
    if input.is_empty() {
        return vec![];
    }

    assert!(
        input.windows(2).all(|w| w[0] <= w[1]),
        "Input slice must be monotonically increasing, got {:?}",
        input
    );

    let max_value = (*input.last().unwrap()) as usize;
    // Here the output is actually max + 2 to account for the 0 value.
    let mut compressed = vec![0; max_value + 2];

    for value in 0..=max_value {
        // Straight up stolen from here: https://stackoverflow.com/a/75790348/4295016
        let res: usize = input
            .binary_search_by(|element| match element.cmp(&(value as u16)) {
                Ordering::Equal => Ordering::Less,
                ord => ord,
            })
            .unwrap_err();

        // Here we add 1 to account for the 0 value.
        compressed[value + 1] = res;
    }

    compressed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec_explode_small() {
        let data = vec![0, 0, 5, 5, 5, 7];
        let out = explode_vec(&data);
        assert_eq!(out, vec![1, 1, 1, 1, 1, 4, 4]);
    }

    #[test]
    fn test_vec_explode_empty() {
        let data = vec![];
        let out = explode_vec(&data);
        assert_eq!(out.len(), 0);
    }

    #[test]
    fn test_back_and_forth() {
        let gt_compression = vec![0, 0, 5, 5, 5, 7];
        let gt_exploded = vec![1, 1, 1, 1, 1, 4, 4];

        let real_compressed = compress_vec(&gt_exploded);
        assert_eq!(real_compressed, gt_compression);

        let real_exploded = explode_vec(&gt_compression);
        assert_eq!(real_exploded, gt_exploded);

        let bf = compress_vec(&real_exploded);
        assert_eq!(bf, gt_compression);

        let fb = explode_vec(&bf);
        assert_eq!(fb, real_exploded);
    }

    #[test]
    fn test_vec_compress_small() {
        let data = vec![1, 1, 1, 1, 1, 4, 4];
        let out = compress_vec(&data);
        assert_eq!(out, vec![0, 0, 5, 5, 5, 7]);
    }

    #[test]
    fn test_vec_compress_empty() {
        let data = vec![];
        let out = compress_vec(&data);
        assert_eq!(out.len(), 0);
    }

    #[test]
    fn test_vec_compress_large() {
        let mut data = vec![0];
        data.extend(vec![1; 5]);
        data.extend(vec![4; 19]);
        let out = compress_vec(&data);
        assert_eq!(out, vec![0, 1, 6, 6, 6, 19 + 6]);
    }
}
