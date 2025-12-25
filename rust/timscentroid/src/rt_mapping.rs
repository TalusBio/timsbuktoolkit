//! The main purpose of this module is to provide some ... domain-knowledge
//! to map indices to retention times. (In essence make more idiomatic the conversion
//! between RT in milliseconds and index in the timsTOF data ... and to distinguish
//! between MS1 and MS2 RTs).
//!
//! I am aware that this seems like a single use bad abstraction but
//! on the longer term It will help with more complex RT mapping needs
//! (eg: schemas where MS1 and MS2 RTs are not aligned 1:1).

use serde::{
    Deserialize,
    Serialize,
};
use std::fmt::Debug;
use std::ops::{
    Range,
    RangeInclusive,
};

/// Represents an index into the MS1 cycles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MS1CycleIndex {
    pub index: u32,
}

/// Represents an index into the MS2 windows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct WindowCycleIndex {
    /// TODO: add a field for the window number?
    /// In theory an u16 should be enough here.
    pub index: u32,
}

pub trait RTIndex:
    Copy + PartialEq + Eq + PartialOrd + Ord + std::fmt::Debug + Send + Sync + Debug
{
    fn index(&self) -> usize;
    fn new(index: u32) -> Self;
    fn as_u32(&self) -> u32;
}
impl RTIndex for MS1CycleIndex {
    fn index(&self) -> usize {
        self.index as usize
    }

    fn new(index: u32) -> Self {
        MS1CycleIndex { index }
    }

    fn as_u32(&self) -> u32 {
        self.index
    }
}
impl RTIndex for WindowCycleIndex {
    fn index(&self) -> usize {
        self.index as usize
    }

    fn new(index: u32) -> Self {
        WindowCycleIndex { index }
    }

    fn as_u32(&self) -> u32 {
        self.index
    }
}

#[derive(Clone)]
pub struct CycleToRTMapping<T: RTIndex> {
    // Q: Do I want to newtype the u32 to make clear its meant
    // for RT in milliseconds?
    rt_milis: Vec<u32>,
    phantom: std::marker::PhantomData<T>,
}

fn glimpse_slc_u32(slice: &[u32]) -> String {
    if slice.len() <= 5 {
        format!("{:?}", slice)
    } else {
        format!(
            "[{}, {}, ..., {}, {}] (len={})",
            slice[0],
            slice[1],
            slice[slice.len() - 2],
            slice[slice.len() - 1],
            slice.len()
        )
    }
}

impl<T: RTIndex> Debug for CycleToRTMapping<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CycleToRTMapping<RTIndex>")
            .field("rt_milis", &glimpse_slc_u32(&self.rt_milis))
            .finish()
    }
}

impl<T: RTIndex> CycleToRTMapping<T> {
    pub fn new(rt_milis: Vec<u32>) -> Self {
        // TODO: check that rt_milis is sorted ascendingly?
        Self {
            rt_milis,
            phantom: std::marker::PhantomData,
        }
    }

    pub fn ms_to_closest_index(&self, rt_milis: u32) -> T {
        // Binary search for the closest RT in milliseconds
        let pos = self.rt_milis.partition_point(|x| *x <= rt_milis);
        if pos == 0 {
            T::new(0)
        } else if pos >= self.rt_milis.len() {
            T::new((self.rt_milis.len() - 1) as u32)
        } else {
            let before = self.rt_milis[pos - 1];
            let after = self.rt_milis[pos];
            if rt_milis - before <= after - rt_milis {
                T::new((pos - 1) as u32)
            } else {
                T::new(pos as u32)
            }
        }
    }

    pub fn range_milis(&self) -> (u32, u32) {
        (
            *self.rt_milis.first().unwrap(),
            *self.rt_milis.last().unwrap(),
        )
    }

    pub fn len(&self) -> usize {
        self.rt_milis.len()
    }

    pub fn unpack(self) -> Vec<u32> {
        self.rt_milis
    }

    pub fn rt_milis_for_index(&self, index: &T) -> Result<u32, RTMappingError> {
        self.rt_milis
            .get(index.index())
            .copied()
            .ok_or(RTMappingError::IndexOutOfBounds)
    }

    pub fn get_slice(&self, range: Range<T>) -> Result<&[u32], RTMappingError> {
        let start = range.start.index();
        let end = range.end.index();
        if start >= self.len() || end > self.len() || start > end {
            return Err(RTMappingError::IndexOutOfBounds);
        }
        Ok(&self.rt_milis[start..end])
    }

    pub fn get_inclusive_slice(&self, range: RangeInclusive<T>) -> Result<&[u32], RTMappingError> {
        let start = range.start().index();
        let end = range.end().index();
        if start >= self.len() || end >= self.len() || start > end {
            return Err(RTMappingError::IndexOutOfBounds);
        }
        Ok(&self.rt_milis[start..=end])
    }
}

#[derive(Debug)]
pub enum RTMappingError {
    IndexOutOfBounds,
}

impl Serialize for MS1CycleIndex {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_u32(self.index)
    }
}

impl<'de> Deserialize<'de> for MS1CycleIndex {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let index = u32::deserialize(deserializer)?;
        Ok(MS1CycleIndex { index })
    }
}

impl Serialize for WindowCycleIndex {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_u32(self.index)
    }
}

impl<'de> Deserialize<'de> for WindowCycleIndex {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let index = u32::deserialize(deserializer)?;
        Ok(WindowCycleIndex { index })
    }
}

impl<T: RTIndex> Serialize for CycleToRTMapping<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.rt_milis.serialize(serializer)
    }
}

impl<'de, T: RTIndex> Deserialize<'de> for CycleToRTMapping<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let rt_milis = Vec::<u32>::deserialize(deserializer)?;
        Ok(CycleToRTMapping::new(rt_milis))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ms1_cycle_index_serialization() {
        let index = MS1CycleIndex { index: 42 };
        let serialized = serde_json::to_string(&index).unwrap();
        assert_eq!(serialized, "42");
        let deserialized: MS1CycleIndex = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized.index, 42);
    }

    #[test]
    fn test_window_cycle_index_serialization() {
        let index = WindowCycleIndex { index: 99 };
        let serialized = serde_json::to_string(&index).unwrap();
        assert_eq!(serialized, "99");
        let deserialized: WindowCycleIndex = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized.index, 99);
    }

    #[test]
    fn test_container_serde() {
        let mapping: CycleToRTMapping<MS1CycleIndex> =
            CycleToRTMapping::new(vec![1000, 2000, 3000]);
        let serialized = serde_json::to_string(&mapping).unwrap();
        assert_eq!(serialized, r#"[1000,2000,3000]"#);
        let deserialized: CycleToRTMapping<MS1CycleIndex> =
            serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized.rt_milis, vec![1000, 2000, 3000]);
    }
}
