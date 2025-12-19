use OptionallyRestricted::{
    Restricted,
    Unrestricted,
};
use core::f32;
use half::f16;
use serde::{
    Deserialize,
    Serialize,
};
use timscentroid::utils::{
    OptionallyRestricted,
    TupleRange,
};

/// Tolerance settings for the search.
///
/// This is meant to encapsulate the different tolerances needed
/// for all the dimensions in a search.
///
/// Example:
/// ```
/// use timsquery::Tolerance;
///
/// let tolerance = Tolerance::default();
/// ```
///
/// Since every dimenion has a different way of defining its tolerance
/// every dimension has its own type (highly encourage you to thech the options
/// for each dimension).
///
/// Convention:
/// In contrast with how some software defines tolerance, here we define the ranges
/// in terms of positive values. For instance, here a tolerance of (1,1) on a value
/// of 10 means a range of (9,11) while in some software the same range would be defined
/// as (-1,1).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Tolerance {
    pub ms: MzTolerance,
    #[serde(default)]
    pub rt: RtTolerance,
    pub mobility: MobilityTolerance,
    pub quad: QuadTolerance,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum MzTolerance {
    #[serde(rename = "da")]
    Absolute((f64, f64)),
    #[serde(rename = "ppm")]
    Ppm((f64, f64)),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
pub enum RtTolerance {
    #[serde(rename = "minutes")]
    Minutes((f32, f32)),
    #[serde(rename = "percent")]
    Pct((f32, f32)),
    #[default]
    Unrestricted,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum MobilityTolerance {
    #[serde(rename = "absolute")]
    Absolute((f32, f32)),
    #[serde(rename = "percent")]
    Pct((f32, f32)),
    Unrestricted,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum QuadTolerance {
    #[serde(rename = "absolute")]
    Absolute((f32, f32)),
}

impl Default for Tolerance {
    fn default() -> Self {
        Tolerance {
            ms: MzTolerance::Ppm((20.0, 20.0)),
            rt: RtTolerance::Minutes((5.0, 5.0)),
            mobility: MobilityTolerance::Pct((3.0, 3.0)),
            quad: QuadTolerance::Absolute((0.1, 0.1)),
        }
    }
}

impl Tolerance {
    // ============================================================================
    // M/Z Tolerance Methods
    // ============================================================================

    /// Calculate m/z tolerance range (primary method, returns f64).
    ///
    /// This is the canonical m/z range method. All other m/z methods delegate to this.
    ///
    /// # Arguments
    ///
    /// * `mz` - Target m/z value in daltons
    ///
    /// # Returns
    ///
    /// A `TupleRange<f64>` representing `[mz - tolerance_low, mz + tolerance_high]`.
    /// The range is always restricted (m/z must have bounds).
    ///
    /// # Tolerance Types
    ///
    /// - `MzTolerance::Absolute((low, high))`: Fixed dalton range
    /// - `MzTolerance::Ppm((low, high))`: Parts-per-million range
    ///
    /// # Example
    ///
    /// ```
    /// use timsquery::Tolerance;
    ///
    /// let tol = Tolerance::default();  // 20 ppm default
    /// let range = tol.mz_range(500.0);
    ///
    /// // For 500 Da at 20 ppm: ±0.01 Da
    /// assert!((range.start() - 499.99).abs() < 0.001);
    /// assert!((range.end() - 500.01).abs() < 0.001);
    ///
    /// // Range is always bounded for m/z
    /// assert!(range.start() < range.end());
    /// ```
    pub fn mz_range(&self, mz: f64) -> TupleRange<f64> {
        match self.ms {
            MzTolerance::Absolute((low, high)) => (mz - low, mz + high).try_into().expect(
                "mz tolerance should never result in an invalid range, since low and high are positive",
            ),
            MzTolerance::Ppm((low, high)) => {
                let low = mz * low / 1e6;
                let high = mz * high / 1e6;
                (mz - low, mz + high).try_into().expect(
                    "mz tolerance should never result in an invalid range, since low and high are positive",
                )
            }
        }
    }

    /// Calculate m/z tolerance range (convenience method, returns f32).
    ///
    /// Same as [`mz_range`](Self::mz_range) but accepts and returns `f32`.
    pub fn mz_range_f32(&self, mz: f32) -> TupleRange<f32> {
        let tmp = self.mz_range(mz as f64);
        (tmp.start() as f32, tmp.end() as f32).try_into().unwrap()
    }

    // ============================================================================
    // Retention Time Tolerance Methods
    // ============================================================================
    //
    // Note: RT methods return OptionallyRestricted because RT can be Unrestricted
    // (match all retention times), unlike m/z which must always have bounds.

    /// Calculate RT tolerance range in seconds as f16 (convenience method).
    ///
    /// Accepts RT in **seconds**, returns range in **seconds** as half-precision floats.
    ///
    /// Delegates to [`rt_range_minutes`](Self::rt_range_minutes) internally.
    pub fn rt_range_seconds_f16(&self, rt_seconds: f32) -> OptionallyRestricted<TupleRange<f16>> {
        let minutes = rt_seconds / 60.0;
        let tmp = self.rt_range_minutes(minutes);
        let sixty = f16::from_f32(60.0);
        // TODO: make this significantly more efficient

        match tmp {
            Restricted(x) => Restricted(
                (
                    f16::from_f32(x.start()) * sixty,
                    f16::from_f32(x.end()) * sixty,
                )
                    .try_into()
                    .unwrap(),
            ),
            Unrestricted => Unrestricted,
        }
    }

    /// Calculate RT tolerance range in minutes (primary method, returns f32).
    ///
    /// This is the canonical RT range method. All other RT methods delegate to this.
    ///
    /// # Arguments
    ///
    /// * `rt_minutes` - Target retention time in **minutes**
    ///
    /// # Returns
    ///
    /// - `Restricted(range)` - RT tolerance window in minutes
    /// - `Unrestricted` - No RT filtering (match all retention times)
    ///
    /// # Tolerance Types
    ///
    /// - `RtTolerance::Minutes((low, high))`: Fixed minute range
    /// - `RtTolerance::Pct((low, high))`: Percentage of RT value
    /// - `RtTolerance::Unrestricted`: No RT bounds
    pub fn rt_range_minutes(&self, rt_minutes: f32) -> OptionallyRestricted<TupleRange<f32>> {
        match self.rt {
            RtTolerance::Minutes((low, high)) => {
                Restricted((rt_minutes - low, rt_minutes + high).try_into().unwrap())
            }
            RtTolerance::Pct((low, high)) => {
                let low = rt_minutes * low / 100.0;
                let high = rt_minutes * high / 100.0;
                Restricted((rt_minutes - low, rt_minutes + high).try_into().unwrap())
            }
            RtTolerance::Unrestricted => Unrestricted,
        }
    }

    /// Calculate RT tolerance range in milliseconds as u32 (convenience method).
    ///
    /// Accepts RT in **seconds**, returns range in **milliseconds** as unsigned integers.
    /// Useful when working with frame/cycle indices that are stored as millisecond timestamps.
    ///
    /// **Unit conversions**: seconds → minutes → tolerance calc → seconds → milliseconds
    pub fn rt_range_as_milis(&self, rt_seconds: f32) -> OptionallyRestricted<TupleRange<u32>> {
        let minutes = rt_seconds / 60.0;
        let tmp = self.rt_range_minutes(minutes);
        match tmp {
            Restricted(x) => {
                let start_seconds = x.start() * 60.0;
                let end_seconds = x.end() * 60.0;
                Restricted(
                    (
                        (start_seconds * 1000.0) as u32,
                        (end_seconds * 1000.0) as u32,
                    )
                        .try_into()
                        .unwrap(),
                )
            }
            Unrestricted => Unrestricted,
        }
    }

    // ============================================================================
    // Ion Mobility Tolerance Methods
    // ============================================================================
    //
    // Note: Mobility methods return OptionallyRestricted because mobility can be
    // Unrestricted (match all mobilities).

    /// Calculate ion mobility tolerance range (primary method, returns f32).
    ///
    /// This is the canonical mobility range method.
    ///
    /// # Arguments
    ///
    /// * `mobility` - Target ion mobility value (1/K0 in V*s/cm^2)
    ///
    /// # Returns
    ///
    /// - `Restricted(range)` - Mobility tolerance window
    /// - `Unrestricted` - No mobility filtering
    ///
    /// # Tolerance Types
    ///
    /// - `MobilityTolerance::Absolute((low, high))`: Fixed mobility range
    /// - `MobilityTolerance::Pct((low, high))`: Percentage of mobility value
    /// - `MobilityTolerance::Unrestricted`: No mobility bounds
    pub fn mobility_range(&self, mobility: f32) -> OptionallyRestricted<TupleRange<f32>> {
        match self.mobility {
            MobilityTolerance::Absolute((low, high)) => {
                Restricted((mobility - low, mobility + high).try_into().unwrap())
            }
            MobilityTolerance::Pct((low, high)) => {
                let low = mobility * (low / 100.0);
                let high = mobility * (high / 100.0);
                Restricted((mobility - low, mobility + high).try_into().unwrap())
            }
            MobilityTolerance::Unrestricted => Unrestricted,
        }
    }

    /// Calculate ion mobility tolerance range (convenience method, returns f16).
    ///
    /// Same as [`mobility_range`](Self::mobility_range) but returns half-precision floats.
    /// Useful for memory-efficient storage of mobility ranges.
    pub fn mobility_range_f16(&self, mobility: f32) -> OptionallyRestricted<TupleRange<f16>> {
        let tmp = self.mobility_range(mobility);
        match tmp {
            Restricted(x) => Restricted(
                (f16::from_f32(x.start()), f16::from_f32(x.end()))
                    .try_into()
                    .unwrap(),
            ),
            Unrestricted => Unrestricted,
        }
    }

    // ============================================================================
    // Quadrupole Isolation Tolerance Methods
    // ============================================================================
    //
    // Quadrupole methods expand precursor isolation windows by the quad tolerance.
    // Always returns restricted ranges (quad isolation must have bounds).

    /// Calculate quadrupole isolation range (convenience method, returns f32).
    ///
    /// Accepts precursor m/z range as `(f32, f32)`, delegates to [`quad_range`](Self::quad_range).
    pub fn quad_range_f32(&self, precursor_mz_range: (f32, f32)) -> TupleRange<f32> {
        let tmp = self.quad_range((precursor_mz_range.0 as f64, precursor_mz_range.1 as f64));
        (tmp.start() as f32, tmp.end() as f32).try_into().unwrap()
    }

    /// Calculate quadrupole isolation range (primary method, returns f64).
    ///
    /// Expands the precursor m/z isolation window by the quad tolerance.
    ///
    /// # Arguments
    ///
    /// * `precursor_mz_range` - Tuple of `(min_mz, max_mz)` for precursor isolation window
    ///
    /// # Returns
    ///
    /// Expanded quadrupole range: `[(min - tol_low), (max + tol_high)]`
    ///
    /// # Tolerance Types
    ///
    /// - `QuadTolerance::Absolute((low, high))`: Fixed dalton expansion
    ///
    /// # Usage
    ///
    /// In dia PASEF, the quadrupole isolates a precursor m/z window.
    /// This method expands that window by the tolerance to account for
    /// quad isolation inaccuracy.
    pub fn quad_range(&self, precursor_mz_range: (f64, f64)) -> TupleRange<f64> {
        match self.quad {
            QuadTolerance::Absolute((low, high)) => {
                let mz_low = precursor_mz_range.0.min(precursor_mz_range.1) - (low as f64);
                let mz_high = precursor_mz_range.1.max(precursor_mz_range.0) + (high as f64);
                assert!(mz_low <= mz_high);
                assert!(
                    mz_low > 0.0,
                    "Precursor mz is 0 or less, inputs: self: {:?}, precursor_mz_range: {:?}",
                    self,
                    precursor_mz_range,
                );
                (mz_low, mz_high).try_into().unwrap()
            }
        }
    }

    // ============================================================================
    // Indexed Domain Conversion Methods
    // ============================================================================
    //
    // These methods convert tolerance ranges from physical units (Da, 1/K0)
    // to instrument-specific integer indices (TOF, scan number).

    // ============================================================================
    // Builder Methods
    // ============================================================================

    /// Create a new `Tolerance` with modified RT tolerance (builder pattern).
    ///
    /// # Example
    ///
    /// ```
    /// use timsquery::Tolerance;
    /// use timsquery::models::tolerance::RtTolerance;
    /// use timscentroid::utils::OptionallyRestricted;
    ///
    /// let tol = Tolerance::default()
    ///     .with_rt_tolerance(RtTolerance::Unrestricted);
    ///
    /// // Verify RT is now unrestricted
    /// let rt_range = tol.rt_range_minutes(10.0);
    /// assert!(matches!(rt_range, OptionallyRestricted::Unrestricted));
    /// ```
    pub fn with_rt_tolerance(self, rt: RtTolerance) -> Self {
        Self { rt, ..self }
    }

    /// Create a new `Tolerance` with modified mobility tolerance (builder pattern).
    ///
    /// # Example
    ///
    /// ```
    /// use timsquery::Tolerance;
    /// use timsquery::models::tolerance::MobilityTolerance;
    /// use timscentroid::utils::OptionallyRestricted;
    ///
    /// let tol = Tolerance::default()
    ///     .with_mobility_tolerance(MobilityTolerance::Absolute((0.05, 0.05)));
    ///
    /// // Verify mobility tolerance is applied
    /// let mob_range = tol.mobility_range(1.0);
    /// match mob_range {
    ///     OptionallyRestricted::Restricted(range) => {
    ///         assert_eq!(range.start(), 0.95);
    ///         assert_eq!(range.end(), 1.05);
    ///     }
    ///     OptionallyRestricted::Unrestricted => panic!("Expected restricted range"),
    /// }
    /// ```
    pub fn with_mobility_tolerance(self, tol: MobilityTolerance) -> Self {
        Self {
            mobility: tol,
            ..self
        }
    }
}
