/// A trait that defines how to aggregate items.
///
/// The `Item` type is the type of the item that is being aggregated.
/// The `Output` type is the type of the output of the aggregation.
/// The `Context` type is the type of the context that is being aggregated.
///
/// The `add` method takes an item of type `Item` OR a type that
/// imlements `Into<Item>`.
///
/// The `set_context` method sets the context that is being aggregated.
/// Note that if the aggregator does not support contexts, this method
/// should never be called. (thus the default implementation is empty).
///
/// The `supports_context` method returns true if the aggregator supports
/// contexts.
///
/// The `finalize` method returns the output of the aggregation.
pub trait Aggregator: Send + Sync {
    type Item: Send + Sync + Clone;
    type Output: Send + Sync;
    type Context: Send + Sync + std::fmt::Debug;

    fn add(&mut self, item: impl Into<Self::Item>);
    fn finalize(self) -> Self::Output;
    fn supports_context(&self) -> bool;
    fn set_context(&mut self, context: Self::Context) {
        panic!(
            "Misconfigured aggregator, context not supported, got {:?}",
            context
        );
    }
    fn get_context(&self) -> Self::Context {
        panic!("Misconfigured aggregator, context not supported");
    }
}

/// Trait purely for the purpose of semantic meaning.
///
/// Every thread safe type automatically implements this trait.
/// For the context `NoContext`.
pub trait ProvidesContext {
    type Context: Send + Sync;
    fn provides_context(&self) -> bool {
        false
    }
}

/// Dummy type to denote that no context is provided.
///
/// Felt more semantically correct to make a type rather
/// than passing `()` to the context.
#[derive(Debug, Clone, Copy)]
pub struct NoContext {}
