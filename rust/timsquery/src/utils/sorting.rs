/// Returns the top n elements of a slice.
///
/// The indices of the elements are also returned.
///
/// # Example
/// ```
/// use timsquery::utils::sorting::top_n;
///
/// let v = vec![1, 2, 10, 4, 5, 9, 7, 8, 2, 1];
/// let (top, indices) = top_n(&v, 3);
/// assert_eq!(top, vec![10, 9, 8]);
/// assert_eq!(indices, vec![2, 5, 7]);
/// ```
///
/// Reference: <https://users.rust-lang.org/t/solved-best-way-to-find-largest-three-values-in-unsorted-slice/34754/12>
pub fn top_n<T: PartialOrd + Copy>(slice: &[T], n: usize) -> (Vec<T>, Vec<usize>) {
    let mut result = Vec::with_capacity(n + 1);
    let mut result_indices = Vec::with_capacity(n + 1);
    for (ind, v) in slice.iter().copied().enumerate() {
        result.push(v);
        result_indices.push(ind);

        let mut i = result.len() - 1;
        while i > 0 {
            if result[i] <= result[i - 1] {
                break;
            }

            // swap
            let t = result[i];
            let tind = result_indices[i];
            result[i] = result[i - 1];
            result[i - 1] = t;
            result_indices[i] = result_indices[i - 1];
            result_indices[i - 1] = tind;

            i -= 1;
        }

        if result.len() > n {
            result.pop();
            result_indices.pop();
        }
    }
    (result, result_indices)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_top_n() {
        let v = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let (top, indices) = top_n(&v, 3);
        assert_eq!(top, vec![10, 9, 8]);
        assert_eq!(indices, vec![9, 8, 7]);
    }
}
