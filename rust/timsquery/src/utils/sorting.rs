/// Macro that sorts an arbitrary number of vecs by a the values
/// first one.
///
/// NOTE: This macro creates a new ordered vec for each one.
/// In theory its possible to have this happen in-place (see commit history)
/// but it seems ineficient for the most part when I benchmarked it.
///
/// TODO: Make a variant that sort in parallel.
/// TODO: Add a parameter to specify ascending or descending.
///
/// # Example
/// ```
/// use timsquery::sort_vecs_by_first;
///
/// let va = vec![9, 8, 7];
/// let vb = vec![1, 2, 3];
/// let vc = vec!['a', 'b', 'c'];
/// let out = sort_vecs_by_first!(&va, &vb, &vc);
///
/// assert_eq!(out.0, vec![7, 8, 9]);
/// assert_eq!(out.1, vec![3, 2, 1]);
/// assert_eq!(out.2, vec!['c', 'b', 'a']);
/// ```
#[macro_export]
macro_rules! sort_vecs_by_first {
    ($first:expr $(,$rest:expr)*) => {{
        let first_vec = $first;
        let len = first_vec.len();

        // Create and sort indices
        let mut indices: Vec<_> = (0..len).collect();
        indices.sort_unstable_by_key(|&i| &first_vec[i]);

        // Reorder first vector
        let sorted_first: Vec<_> = indices.iter().map(|&i| first_vec[i]).collect();

        // Reorder all other vectors
        (sorted_first, $( {
            let other_vec = $rest;
            assert_eq!(other_vec.len(), len, "All vectors must have the same length");
            indices.iter().map(|&i| other_vec[i]).collect::<Vec<_>>()
        }, )*)
    }};
}

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
/// https://users.rust-lang.org/t/solved-best-way-to-find-largest-three-values-in-unsorted-slice/34754/12
pub fn top_n(slice: &[u32], n: usize) -> (Vec<u32>, Vec<usize>) {
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
    fn test_sort_two_vecs() {
        let v1 = vec![3, 1, 4, 1, 5];
        let v2 = vec!['a', 'b', 'c', 'd', 'e'];

        let (sorted_v1, sorted_v2) = sort_vecs_by_first!(v1, v2);

        assert_eq!(sorted_v1, vec![1, 1, 3, 4, 5]);
        assert_eq!(sorted_v2, vec!['b', 'd', 'a', 'c', 'e']);
    }

    #[test]
    fn test_sort_three_vecs() {
        let v1 = vec![3, 1, 4];
        let v2 = vec!['x', 'y', 'z'];
        let v3 = vec![true, false, true];

        let (sorted_v1, sorted_v2, sorted_v3) = sort_vecs_by_first!(v1, v2, v3);

        assert_eq!(sorted_v1, vec![1, 3, 4]);
        assert_eq!(sorted_v2, vec!['y', 'x', 'z']);
        assert_eq!(sorted_v3, vec![false, true, true]);
    }

    #[test]
    fn test_sort_four_vecs() {
        let v1 = vec![3, 1, 4];
        let v2 = vec!['x', 'y', 'z'];
        let v3 = vec![true, false, true];
        let v4 = vec![1.0, 2.0, 3.0];

        let (sorted_v1, sorted_v2, sorted_v3, sorted_v4) = sort_vecs_by_first!(v1, v2, v3, v4);

        assert_eq!(sorted_v1, vec![1, 3, 4]);
        assert_eq!(sorted_v2, vec!['y', 'x', 'z']);
        assert_eq!(sorted_v3, vec![false, true, true]);
        assert_eq!(sorted_v4, vec![2.0, 1.0, 3.0]);
    }

    #[test]
    fn test_top_n() {
        let v = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let (top, indices) = top_n(&v, 3);
        assert_eq!(top, vec![10, 9, 8]);
        assert_eq!(indices, vec![9, 8, 7]);
    }
}
