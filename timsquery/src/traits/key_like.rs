use serde::Serialize;
use std::fmt::Debug;
use std::hash::Hash;

pub trait KeyLike: Clone + Eq + Serialize + Hash + Send + Sync + Debug {}

impl<T: Clone + Eq + Serialize + Hash + Send + Sync + Debug> KeyLike for T {}
