use serde::Serialize;
use std::fmt::Debug;
use std::hash::Hash;

/// Marker trait for types that can be used as keys in aggregator collections.
///
/// This trait is automatically implemented for all types that satisfy the required bounds.
/// Keys are used to identify and group data (e.g., precursor IDs, query identifiers).
///
/// # Required Bounds
///
/// - `Clone`: Keys must be cloneable for internal data structures
/// - `Eq`: Keys must be comparable for lookups and deduplication
/// - `Serialize`: Keys must be serializable for persistence and debugging
/// - `Hash`: Keys must be hashable for efficient map/set operations
/// - `Send + Sync`: Keys must be thread-safe for parallel query execution
/// - `Debug`: Keys must be debuggable for error messages and logging
///
/// # Examples
///
/// Common key types include:
/// ```
/// // Unsigned integers (precursor indices)
/// let _key: usize = 42;
///
/// // String identifiers
/// let _key: String = "precursor_001".to_string();
///
/// // Tuples for compound keys
/// let _key: (usize, u8) = (42, 2); // (precursor_id, isotope_offset)
/// ```
pub trait KeyLike: Clone + Eq + Serialize + Hash + Send + Sync + Debug + Default {}

/// Marker trait for types that can be used as values in aggregator results.
///
/// This trait is automatically implemented for all types that satisfy the required bounds.
/// Values represent aggregated data (e.g., intensities, counts, statistics).
///
/// # Required Bounds
///
/// - `Copy`: Values must be trivially copyable (typically numeric types)
/// - `Clone`: Implied by `Copy`, but explicit for clarity
/// - `Serialize`: Values must be serializable for output and analysis
/// - `Send + Sync`: Values must be thread-safe for parallel aggregation
/// - `Debug`: Values must be debuggable for verification
///
/// Note: The `Copy` bound restricts values to small, stack-allocated types like
/// numeric primitives (`f32`, `f64`, `u32`, etc.) and small aggregates.
///
/// # Examples
///
/// Common value types include:
/// ```
/// // Floating-point intensities
/// let _value: f32 = 1234.5;
///
/// // Integer counts
/// let _value: u32 = 42;
///
/// // Half-precision for memory efficiency
/// use half::f16;
/// let _value: f16 = f16::from_f32(123.4);
/// ```
pub trait ValueLike: Copy + Clone + Serialize + Send + Sync + Debug {}

impl<T: Clone + Eq + Serialize + Hash + Send + Sync + Debug + Default> KeyLike for T {}
impl<T: Copy + Clone + Serialize + Send + Sync + Debug> ValueLike for T {}
