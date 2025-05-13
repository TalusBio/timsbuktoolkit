use crate::models::frames::ResolvedPeakInQuad;
use serde::Serialize;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::AddAssign;

pub trait KeyLike: Clone + Eq + Serialize + Hash + Send + Sync + Debug {}
pub trait ValueLike: Copy + Clone + Serialize + Send + Sync + Debug {}
pub trait PeakAddable: AddAssign<ResolvedPeakInQuad> + ValueLike + Default {}

impl<T: Clone + Eq + Serialize + Hash + Send + Sync + Debug> KeyLike for T {}
impl<T: Copy + Clone + Serialize + Send + Sync + Debug> ValueLike for T {}
impl<T: AddAssign<ResolvedPeakInQuad> + ValueLike + Default> PeakAddable for T {}
