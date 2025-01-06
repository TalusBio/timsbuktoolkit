use crate::errors::{
    DataProcessingError,
    Result,
    TimsSeekError,
};
use std::cmp::PartialOrd;

/// Aligns a list of values to a list of reference positions.
///
/// In the simplest terms, provided two vecs that are the same length
/// one with values and one with reference positions, this function
/// returns a list of values with the same length as the reference positions
/// where the values are aligned to the reference positions.
///
/// This function preserves only the largest value at each reference position.
///
/// # Example
///
/// ```
/// use timsseek::utils::aligning::snap_to_reference;
///
/// let values = vec![1, 2, 3];
/// let value_positions = vec![1.0, 2.0, 3.0];
/// let reference_positions = vec![0.0, 1.0, 2.0, 3.0, 4.0];
/// let out = snap_to_reference(&values, &value_positions, &reference_positions).unwrap();
///
/// assert_eq!(out, vec![0, 1, 2, 3, 0]);
pub fn snap_to_reference<R: PartialOrd, V: Copy + Default + PartialOrd>(
    values: &[V],
    value_positions: &[R],
    reference_positions: &[R],
) -> Result<Vec<V>> {
    if values.len() != value_positions.len() {
        return Err(TimsSeekError::DataProcessingError(
            DataProcessingError::ExpectedSlicesSameLength {
                expected: value_positions.len(),
                other: values.len(),
                context: "snap_to_reference".to_string(),
            },
        ));
    }
    let mut out = vec![V::default(); reference_positions.len()];

    let mut ref_i = 0;
    let limit = reference_positions.len() - 1;
    for (vv, vp) in values.iter().zip(value_positions.iter()) {
        while ref_i < limit && reference_positions[ref_i] < *vp {
            ref_i += 1;
        }
        // Replace only if its larger
        let cmp = out[ref_i].partial_cmp(vv);
        match cmp {
            Some(std::cmp::Ordering::Less) => {
                out[ref_i] = *vv;
            }
            _ => {}
        }
    }
    Ok(out)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_snap_to_reference_same() {
        let values = vec![1, 2, 3, 4, 5];
        let value_positions = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let reference_positions = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        assert_eq!(
            snap_to_reference(&values, &value_positions, &reference_positions).unwrap(),
            vec![1, 2, 3, 4, 5],
        );
    }

    #[test]
    fn test_snap_to_reference_reflonger() {
        let values = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let value_positions = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let reference_positions = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(
            snap_to_reference(&values, &value_positions, &reference_positions).unwrap(),
            vec![1, 2, 3, 4, 5, 10],
        );
    }

    #[test]
    fn test_snap_to_reference_ref_shorter() {
        let values = vec![1, 2, 3];
        let value_positions = vec![3.0, 4.0, 5.0];
        let reference_positions = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        assert_eq!(
            snap_to_reference(&values, &value_positions, &reference_positions).unwrap(),
            vec![0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0],
        );
    }
}
