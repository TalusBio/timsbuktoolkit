use std::fmt::Display;
use std::hash::Hash;

use timsrust::QuadrupoleSettings;

use crate::OptionallyRestricted;
use crate::utils::tolerance_ranges::IncludedRange;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SingleQuadrupoleSettingIndex {
    pub major_index: usize,
    pub sub_index: usize,
}

impl Display for SingleQuadrupoleSettingIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

// TOOD: There is a lot of instances where I use raw IncludedRange<u16> ...
// I think I can just alias that to a type ... 'IndexedMobilityScanRange'?

#[derive(Debug, Clone, Copy)]
pub struct SingleQuadrupoleSettingRanges {
    pub scan_start: u16,
    pub scan_end: u16,
    // TODO: Should I add here the start-end mobility??
    pub isolation_mz: f64,
    pub isolation_width: f64,
    pub isolation_high: f64,
    pub isolation_low: f64,
    pub collision_energy: f64,
}

impl Display for SingleQuadrupoleSettingRanges {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SingleQuadrupoleSettingRanges {{ scan_start: {}, scan_end: {}, isolation_mz: {}, isolation_width: {}, isolation_high: {}, isolation_low: {}, collision_energy: {} }}",
            self.scan_start,
            self.scan_end,
            self.isolation_mz,
            self.isolation_width,
            self.isolation_high,
            self.isolation_low,
            self.collision_energy
        )
    }
}

/// A single quadrupole setting.
///
/// NOTE: Only the index is used for the hashing of this struct.
/// So make sure that int he buildingn process there is no generation
/// of settings with the same index.
#[derive(Debug, Clone, Copy)]
pub struct SingleQuadrupoleSetting {
    pub index: SingleQuadrupoleSettingIndex,
    pub ranges: SingleQuadrupoleSettingRanges,
}

impl Display for SingleQuadrupoleSetting {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SingleQuadrupoleSetting {{\n index: {},\n ranges: {}\n }}",
            self.index, self.ranges
        )
    }
}

impl Hash for SingleQuadrupoleSetting {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.index.hash(state);
    }
}

impl PartialEq for SingleQuadrupoleSetting {
    fn eq(&self, other: &Self) -> bool {
        let out = self.index == other.index;
        if out && cfg!(debug_assertions) {
            debug_assert!(
                self.ranges.scan_start == other.ranges.scan_start,
                "{:?} != {:?} in scan_start",
                self.ranges.scan_start,
                other.ranges.scan_start
            );
            debug_assert!(self.ranges.scan_end == other.ranges.scan_end);
            debug_assert!(self.ranges.isolation_mz == other.ranges.isolation_mz);
            debug_assert!(self.ranges.isolation_width == other.ranges.isolation_width);
            debug_assert!(self.ranges.isolation_high == other.ranges.isolation_high);
            debug_assert!(self.ranges.isolation_low == other.ranges.isolation_low);
            debug_assert!(self.ranges.collision_energy == other.ranges.collision_energy);
        };
        out
    }
}

impl Eq for SingleQuadrupoleSetting {}

#[derive(Debug, Clone, Hash, PartialEq)]
pub enum ExpandedFrameQuadSettings {
    Unfragmented,
    Fragmented(Vec<SingleQuadrupoleSetting>),
}

pub fn expand_quad_settings(quad_settings: &QuadrupoleSettings) -> ExpandedFrameQuadSettings {
    if quad_settings.scan_ends.is_empty() {
        return ExpandedFrameQuadSettings::Unfragmented;
    }

    let mut out = Vec::with_capacity(quad_settings.scan_ends.len());
    for i in 0..quad_settings.scan_ends.len() {
        let isolation_width = quad_settings.isolation_width[i];
        let half_width = isolation_width / 2.0;
        let isolation_center = quad_settings.isolation_mz[i];
        let collision_energy = quad_settings.collision_energy[i];
        let scan_start = quad_settings.scan_starts[i]
            .try_into()
            .expect("Scan number to fit u16");
        let scan_end = quad_settings.scan_ends[i]
            .try_into()
            .expect("Scan number to fit u16");

        let index = SingleQuadrupoleSettingIndex {
            major_index: quad_settings.index,
            sub_index: i,
        };

        let ranges = SingleQuadrupoleSettingRanges {
            scan_start,
            scan_end,
            isolation_mz: quad_settings.isolation_mz[i],
            isolation_width,
            isolation_high: isolation_center + half_width,
            isolation_low: isolation_center - half_width,
            collision_energy,
        };

        out.push(SingleQuadrupoleSetting { index, ranges });
    }
    ExpandedFrameQuadSettings::Fragmented(out)
}

pub fn get_matching_quad_settings(
    flat_quad_settings: &[SingleQuadrupoleSetting],
    precursor_mz_range: IncludedRange<f64>,
    scan_range: OptionallyRestricted<IncludedRange<u16>>,
) -> impl Iterator<Item = SingleQuadrupoleSettingIndex> + '_ {
    flat_quad_settings
        .iter()
        .filter(move |qs| {
            qs.matches_iso_window(precursor_mz_range) && qs.matches_scan_range(scan_range)
        })
        .map(|e| e.index)
}

impl SingleQuadrupoleSetting {
    fn matches_iso_window(&self, precursor_mz_range: IncludedRange<f64>) -> bool {
        (self.ranges.isolation_low <= precursor_mz_range.end())
            && (self.ranges.isolation_high >= precursor_mz_range.start())
    }

    fn matches_scan_range(&self, scan_range: OptionallyRestricted<IncludedRange<u16>>) -> bool {
        match scan_range {
            OptionallyRestricted::Restricted(tmp_range) => {
                let min_scan = tmp_range.start();
                let max_scan = tmp_range.end();
                assert!(self.ranges.scan_start <= self.ranges.scan_end);

                // Above quad
                // Quad                   [----------]
                // Query                               [------]
                let above_quad = self.ranges.scan_end < min_scan;

                // Below quad
                // Quad                  [------]
                // Query       [------]
                let below_quad = self.ranges.scan_start > max_scan;

                if above_quad || below_quad {
                    // This quad is completely outside the scan range
                    false
                } else {
                    true
                }
            }
            OptionallyRestricted::Unrestricted => true,
        }
    }

    pub fn matches_quad_settings(
        &self,
        precursor_mz_range: IncludedRange<f64>,
        scan_range: OptionallyRestricted<IncludedRange<u16>>,
    ) -> bool {
        self.matches_iso_window(precursor_mz_range) && self.matches_scan_range(scan_range)
    }
}
