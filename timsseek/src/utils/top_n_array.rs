/// Preserves in sorted order the top N elements.
///
/// Notes:
/// This requires the default value to be comparable.
/// Values that return None when comparing to the default value are ignored.
///
/// # Example
///
/// ```
/// use timsseek::utils::top_n_array::TopNArray;
///
/// let mut top_3 = TopNArray::<3, f64>::new();
/// top_3.push(1.0);
/// top_3.push(2.0);
/// top_3.push(3.0);
/// top_3.push(3.0);
/// top_3.push(12.0);
///
/// assert_eq!(top_3.get_values(), [12.0, 3.0, 3.0]);
/// assert_eq!(top_3.max_val(), 12.0);
/// assert_eq!(top_3.min_val(), 3.0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct TopNArray<const N: usize, T: Clone + Copy + PartialOrd + Default> {
    array: [T; N],
    len: usize,
}

impl<const N: usize, T: Clone + Copy + PartialOrd + Default> Default for TopNArray<N, T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize, T: Clone + Copy + PartialOrd + Default> TopNArray<N, T> {
    pub fn new() -> Self {
        Self {
            array: [T::default(); N],
            len: 0,
        }
    }

    pub fn max_val(&self) -> T {
        self.array[0]
    }

    pub fn min_val(&self) -> T {
        self.array[self.len - 1]
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn capacity(&self) -> usize {
        N
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn get_values(&self) -> [T; N] {
        self.array
    }

    pub fn push(&mut self, value: T) {
        if self.len < N {
            match self.array[self.len].partial_cmp(&value) {
                Some(_) => {}
                None => {
                    return;
                }
            }
            self.array[self.len] = value;
            self.sort_val(self.len);
            self.len += 1;
        } else {
            // Replace last value if it is lower
            // Consider NaN as always being lower
            let last = self.array[self.len - 1];

            match value.partial_cmp(&last) {
                Some(std::cmp::Ordering::Greater) => {
                    self.array[self.len - 1] = value;
                }
                Some(_) => {
                    return;
                }
                None => {
                    return;
                }
            }

            self.sort_val(self.len - 1);
        }
    }

    fn sort_val(&mut self, idx: usize) {
        // Keep swapping until we get to the right place
        // Loop from idx to 0
        for i in (1..=idx).rev() {
            let curr = self.array[i];
            let prev = self.array[i - 1];

            match curr.partial_cmp(&prev) {
                Some(std::cmp::Ordering::Greater) => {
                    self.array[i] = prev;
                    self.array[i - 1] = curr;
                }
                None => {
                    // Assume incomparables need to be swapped
                    self.array[i] = prev;
                    self.array[i - 1] = curr;
                }
                Some(_) => {
                    continue;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_push_incomplete() {
        let mut top_3 = TopNArray::<3, f64>::new();
        top_3.push(1.0);
        top_3.push(2.0);

        assert_eq!(top_3.array, [2.0, 1.0, 0.0]);

        let minval = top_3.min_val();
        assert_eq!(minval, 1.0);

        let maxval = top_3.max_val();
        assert_eq!(maxval, 2.0);
    }

    #[test]
    fn test_push() {
        let mut top_3 = TopNArray::<3, f64>::new();
        top_3.push(1.0);
        top_3.push(2.0);
        top_3.push(3.0);
        top_3.push(4.0);
        top_3.push(5.0);
        assert_eq!(top_3.array, [5.0, 4.0, 3.0]);

        let minval = top_3.min_val();
        assert_eq!(minval, 3.0);

        let maxval = top_3.max_val();
        assert_eq!(maxval, 5.0);
    }

    #[test]
    fn test_push_nan() {
        let mut top_3 = TopNArray::<3, f64>::new();
        top_3.push(f64::NAN);
        top_3.push(1.0);
        top_3.push(2.0);
        top_3.push(4.0);
        top_3.push(5.0);
        assert_eq!(top_3.array, [5.0, 4.0, 2.0]);

        let minval = top_3.min_val();
        assert_eq!(minval, 2.0);

        let maxval = top_3.max_val();
        assert_eq!(maxval, 5.0);
    }
}
