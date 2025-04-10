use std::cmp::Ordering;

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
pub fn explode_vec(input: &[usize]) -> Vec<usize> {
    let last_val = match input.last() {
        Some(last) => *last,
        None => return Vec::new(),
    };

    let mut out = Vec::with_capacity(last_val);
    for i in 0..(input.len() - 1) {
        let start = input[i];
        let end = input[i + 1];
        for _ in start..end {
            out.push(i);
        }
    }
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
pub fn compress_vec(input: &[usize]) -> Vec<usize> {
    if input.is_empty() {
        return vec![];
    }

    assert!(
        input.windows(2).all(|w| w[0] <= w[1]),
        "Input slice must be monotonically increasing, got {:?}",
        input
    );

    let max_value = *input.last().unwrap();
    // Here the output is actually max + 2 to account for the 0 value.
    let mut compressed = vec![0; max_value + 2];

    for value in 0..=max_value {
        // Straight up stolen from here: https://stackoverflow.com/a/75790348/4295016
        let res: usize = input
            .binary_search_by(|element| match element.cmp(&value) {
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
