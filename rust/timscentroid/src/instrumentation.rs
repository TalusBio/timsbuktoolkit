//! Instrumentation and observability for storage operations
//!
//! This module provides wrappers around ObjectStore that track metrics:
//! - Number of GET/PUT/HEAD operations
//! - Total bytes transferred
//! - Time spent in each operation
//! - Detailed operation logs
//!
//! # Note on Timing Metrics
//!
//! Individual operation times (get_time_us, etc.) represent cumulative time across all operations.
//! When operations execute concurrently, these cumulative times will exceed wall-clock time.
//! For wall-clock measurements, use external timers. Operation counts and bytes transferred
//! remain accurate regardless of concurrency.

use async_trait::async_trait;
use bytes::Bytes;
use futures::stream::BoxStream;
use object_store::path::Path;
use object_store::{
    GetOptions,
    GetResult,
    ListResult,
    MultipartUpload,
    ObjectMeta,
    ObjectStore,
    PutMultipartOptions,
    PutOptions,
    PutPayload,
    PutResult,
    Result,
};
use std::fmt::Display;
use std::ops::Range;
use std::sync::Arc;
use std::sync::atomic::{
    AtomicU64,
    AtomicUsize,
    Ordering,
};

/// Metrics collected during storage operations
#[derive(Debug, Clone)]
pub struct StorageMetrics {
    // Operation counts
    pub get_count: Arc<AtomicUsize>,
    pub put_count: Arc<AtomicUsize>,
    pub head_count: Arc<AtomicUsize>,
    pub delete_count: Arc<AtomicUsize>,
    pub list_count: Arc<AtomicUsize>,

    // Bytes transferred
    pub bytes_read: Arc<AtomicU64>,
    pub bytes_written: Arc<AtomicU64>,

    // Time spent (in microseconds)
    pub get_time_us: Arc<AtomicU64>,
    pub put_time_us: Arc<AtomicU64>,
    pub head_time_us: Arc<AtomicU64>,
}

impl Default for StorageMetrics {
    fn default() -> Self {
        Self::new()
    }
}

impl StorageMetrics {
    pub fn new() -> Self {
        Self {
            get_count: Arc::new(AtomicUsize::new(0)),
            put_count: Arc::new(AtomicUsize::new(0)),
            head_count: Arc::new(AtomicUsize::new(0)),
            delete_count: Arc::new(AtomicUsize::new(0)),
            list_count: Arc::new(AtomicUsize::new(0)),
            bytes_read: Arc::new(AtomicU64::new(0)),
            bytes_written: Arc::new(AtomicU64::new(0)),
            get_time_us: Arc::new(AtomicU64::new(0)),
            put_time_us: Arc::new(AtomicU64::new(0)),
            head_time_us: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Reset all metrics to zero
    pub fn reset(&self) {
        self.get_count.store(0, Ordering::SeqCst);
        self.put_count.store(0, Ordering::SeqCst);
        self.head_count.store(0, Ordering::SeqCst);
        self.delete_count.store(0, Ordering::SeqCst);
        self.list_count.store(0, Ordering::SeqCst);
        self.bytes_read.store(0, Ordering::SeqCst);
        self.bytes_written.store(0, Ordering::SeqCst);
        self.get_time_us.store(0, Ordering::SeqCst);
        self.put_time_us.store(0, Ordering::SeqCst);
        self.head_time_us.store(0, Ordering::SeqCst);
    }

    /// Get snapshot of current metrics
    pub fn snapshot(&self) -> MetricsSnapshot {
        MetricsSnapshot {
            get_count: self.get_count.load(Ordering::SeqCst),
            put_count: self.put_count.load(Ordering::SeqCst),
            head_count: self.head_count.load(Ordering::SeqCst),
            delete_count: self.delete_count.load(Ordering::SeqCst),
            list_count: self.list_count.load(Ordering::SeqCst),
            bytes_read: self.bytes_read.load(Ordering::SeqCst),
            bytes_written: self.bytes_written.load(Ordering::SeqCst),
            get_time_us: self.get_time_us.load(Ordering::SeqCst),
            put_time_us: self.put_time_us.load(Ordering::SeqCst),
            head_time_us: self.head_time_us.load(Ordering::SeqCst),
        }
    }
}

/// Immutable snapshot of metrics
#[derive(Debug, Clone, Copy)]
pub struct MetricsSnapshot {
    pub get_count: usize,
    pub put_count: usize,
    pub head_count: usize,
    pub delete_count: usize,
    pub list_count: usize,
    pub bytes_read: u64,
    pub bytes_written: u64,
    pub get_time_us: u64,
    pub put_time_us: u64,
    pub head_time_us: u64,
}

impl MetricsSnapshot {
    pub fn print_report(&self, label: &str) {
        println!("\n=== Storage Metrics: {} ===", label);
        println!("Operations:");
        println!(
            "  GET:    {:>6} calls, {:>10.2} MB, {:>8.2} ms total, {:>6.2} ms/call",
            self.get_count,
            self.bytes_read as f64 / 1_024_000.0,
            self.get_time_us as f64 / 1_000.0,
            if self.get_count > 0 {
                self.get_time_us as f64 / 1_000.0 / self.get_count as f64
            } else {
                0.0
            }
        );
        println!(
            "  HEAD:   {:>6} calls, {:>8.2} ms total, {:>6.2} ms/call",
            self.head_count,
            self.head_time_us as f64 / 1_000.0,
            if self.head_count > 0 {
                self.head_time_us as f64 / 1_000.0 / self.head_count as f64
            } else {
                0.0
            }
        );
        println!(
            "  PUT:    {:>6} calls, {:>10.2} MB, {:>8.2} ms total",
            self.put_count,
            self.bytes_written as f64 / 1_024_000.0,
            self.put_time_us as f64 / 1_000.0
        );
        println!("  DELETE: {:>6} calls", self.delete_count);
        println!("  LIST:   {:>6} calls", self.list_count);
        println!("\nTotals:");
        println!(
            "  Total operations: {}",
            self.get_count + self.put_count + self.head_count + self.delete_count + self.list_count
        );
        println!(
            "  Total bytes read: {:.2} MB",
            self.bytes_read as f64 / 1_024_000.0
        );
        println!(
            "  Total bytes written: {:.2} MB",
            self.bytes_written as f64 / 1_024_000.0
        );
        println!(
            "  Total I/O time: {:.2} ms",
            (self.get_time_us + self.put_time_us + self.head_time_us) as f64 / 1_000.0
        );
    }
}

/// ObjectStore wrapper that tracks metrics
pub struct InstrumentedStore {
    pub(crate) inner: Arc<dyn ObjectStore>,
    metrics: Arc<StorageMetrics>,
    label: String,
    /// Optional fake latency to simulate network delays (e.g., for testing S3 performance)
    fake_latency: Option<std::time::Duration>,
}

impl InstrumentedStore {
    pub fn new(
        inner: Arc<dyn ObjectStore>,
        metrics: Arc<StorageMetrics>,
        label: impl Into<String>,
    ) -> Self {
        Self {
            inner,
            metrics,
            label: label.into(),
            fake_latency: None,
        }
    }

    /// Set fake latency to simulate network delays (useful for testing cloud storage performance)
    pub fn with_fake_latency(mut self, latency: std::time::Duration) -> Self {
        self.fake_latency = Some(latency);
        self
    }

    pub fn metrics(&self) -> Arc<StorageMetrics> {
        self.metrics.clone()
    }

    pub fn print_metrics(&self) {
        self.metrics.snapshot().print_report(&self.label);
    }

    /// Simulate network latency if configured
    async fn simulate_latency(&self) {
        if let Some(latency) = self.fake_latency {
            tokio::time::sleep(latency).await;
        }
    }
}

impl std::fmt::Debug for InstrumentedStore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "InstrumentedStore({})", self.label)
    }
}

impl Display for InstrumentedStore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "InstrumentedStore({})", self.label)
    }
}

#[async_trait]
impl ObjectStore for InstrumentedStore {
    async fn put(&self, location: &Path, payload: PutPayload) -> Result<PutResult> {
        let start = std::time::Instant::now();
        let bytes_len = payload.content_length();

        let result = self.inner.put(location, payload).await;

        let elapsed = start.elapsed();
        self.metrics.put_count.fetch_add(1, Ordering::SeqCst);
        self.metrics
            .bytes_written
            .fetch_add(bytes_len as u64, Ordering::SeqCst);
        self.metrics
            .put_time_us
            .fetch_add(elapsed.as_micros() as u64, Ordering::SeqCst);

        result
    }

    async fn put_opts(
        &self,
        location: &Path,
        payload: PutPayload,
        opts: PutOptions,
    ) -> Result<PutResult> {
        let start = std::time::Instant::now();
        let bytes_len = payload.content_length();

        let result = self.inner.put_opts(location, payload, opts).await;

        let elapsed = start.elapsed();
        self.metrics.put_count.fetch_add(1, Ordering::SeqCst);
        self.metrics
            .bytes_written
            .fetch_add(bytes_len as u64, Ordering::SeqCst);
        self.metrics
            .put_time_us
            .fetch_add(elapsed.as_micros() as u64, Ordering::SeqCst);

        result
    }

    async fn get(&self, location: &Path) -> Result<GetResult> {
        let start = std::time::Instant::now();

        // Simulate latency first (included in timing)
        self.simulate_latency().await;

        let result = self.inner.get(location).await;
        let elapsed = start.elapsed();

        self.metrics.get_count.fetch_add(1, Ordering::SeqCst);
        self.metrics
            .get_time_us
            .fetch_add(elapsed.as_micros() as u64, Ordering::SeqCst);

        // Track bytes read when the result is consumed
        if let Ok(get_result) = &result {
            let bytes_len = get_result.meta.size;
            self.metrics
                .bytes_read
                .fetch_add(bytes_len, Ordering::SeqCst);
        }

        result
    }

    async fn get_opts(&self, location: &Path, options: GetOptions) -> Result<GetResult> {
        let start = std::time::Instant::now();

        // Simulate latency first (included in timing)
        self.simulate_latency().await;

        let result = self.inner.get_opts(location, options).await;
        let elapsed = start.elapsed();

        self.metrics.get_count.fetch_add(1, Ordering::SeqCst);
        self.metrics
            .get_time_us
            .fetch_add(elapsed.as_micros() as u64, Ordering::SeqCst);

        if let Ok(get_result) = &result {
            let bytes_len = get_result.meta.size;
            self.metrics
                .bytes_read
                .fetch_add(bytes_len, Ordering::SeqCst);
        }

        result
    }

    async fn get_range(&self, location: &Path, range: Range<u64>) -> Result<Bytes> {
        let start = std::time::Instant::now();

        // Simulate latency first (included in timing)
        self.simulate_latency().await;

        let result = self.inner.get_range(location, range.clone()).await;
        let elapsed = start.elapsed();

        self.metrics.get_count.fetch_add(1, Ordering::SeqCst);
        self.metrics
            .get_time_us
            .fetch_add(elapsed.as_micros() as u64, Ordering::SeqCst);

        if let Ok(bytes) = &result {
            self.metrics
                .bytes_read
                .fetch_add(bytes.len() as u64, Ordering::SeqCst);
        }

        result
    }

    async fn head(&self, location: &Path) -> Result<ObjectMeta> {
        let start = std::time::Instant::now();

        // Simulate latency first (included in timing)
        self.simulate_latency().await;

        let result = self.inner.head(location).await;
        let elapsed = start.elapsed();

        self.metrics.head_count.fetch_add(1, Ordering::SeqCst);
        self.metrics
            .head_time_us
            .fetch_add(elapsed.as_micros() as u64, Ordering::SeqCst);

        result
    }

    async fn delete(&self, location: &Path) -> Result<()> {
        self.metrics.delete_count.fetch_add(1, Ordering::SeqCst);
        self.inner.delete(location).await
    }

    fn list(&self, prefix: Option<&Path>) -> BoxStream<'static, Result<ObjectMeta>> {
        self.metrics.list_count.fetch_add(1, Ordering::SeqCst);
        self.inner.list(prefix)
    }

    fn list_with_offset(
        &self,
        prefix: Option<&Path>,
        offset: &Path,
    ) -> BoxStream<'static, Result<ObjectMeta>> {
        self.metrics.list_count.fetch_add(1, Ordering::SeqCst);
        self.inner.list_with_offset(prefix, offset)
    }

    async fn list_with_delimiter(&self, prefix: Option<&Path>) -> Result<ListResult> {
        self.metrics.list_count.fetch_add(1, Ordering::SeqCst);
        self.inner.list_with_delimiter(prefix).await
    }

    async fn copy(&self, from: &Path, to: &Path) -> Result<()> {
        self.inner.copy(from, to).await
    }

    async fn copy_if_not_exists(&self, from: &Path, to: &Path) -> Result<()> {
        self.inner.copy_if_not_exists(from, to).await
    }

    async fn rename_if_not_exists(&self, from: &Path, to: &Path) -> Result<()> {
        self.inner.rename_if_not_exists(from, to).await
    }

    async fn put_multipart(&self, location: &Path) -> Result<Box<dyn MultipartUpload>> {
        self.metrics.put_count.fetch_add(1, Ordering::SeqCst);
        self.inner.put_multipart(location).await
    }

    async fn put_multipart_opts(
        &self,
        location: &Path,
        opts: PutMultipartOptions,
    ) -> Result<Box<dyn MultipartUpload>> {
        self.metrics.put_count.fetch_add(1, Ordering::SeqCst);
        self.inner.put_multipart_opts(location, opts).await
    }
}
