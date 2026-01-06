# timscentroid

Efficient indexing and lazy loading for timsTOF mass spectrometry data.

## Quick Start

### Basic Usage (Local Files)

```rust
use timscentroid::{IndexedTimstofPeaks, CentroidingConfig};
use timscentroid::lazy::LazyIndexedTimstofPeaks;
use timscentroid::utils::{TupleRange, OptionallyRestricted::*};

// 1. Index your timsTOF data
let index = IndexedTimstofPeaks::from_timstof_file(
    &file,
    CentroidingConfig::default()
);

// 2. Save to disk
index.save_to_directory("./indexed_peaks")?;

// 3. Load lazily (only loads metadata, ~50ms)
let lazy_index = LazyIndexedTimstofPeaks::load_from_directory("./indexed_peaks")?;

// 4. Query peaks (loads relevant data on-demand)
let mz_range = TupleRange::try_new(400.0, 500.0)?;
for peak in lazy_index.query_peaks_ms1(mz_range, Unrestricted, Unrestricted) {
    println!("m/z: {}, intensity: {}", peak.mz, peak.intensity);
}
```

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
timscentroid = "0.21"
```

For cloud storage support, enable the appropriate features:

```toml
[dependencies]
timscentroid = { version = "0.21", features = ["aws"] }  # For S3
# or
timscentroid = { version = "0.21", features = ["gcp"] }  # For GCS
# or
timscentroid = { version = "0.21", features = ["azure"] }  # For Azure
```

## Cloud Storage

### Saving to Cloud

```rust
// Save to S3
index.save_to_url("s3://my-bucket/indexed_peaks/")?;

// Save to Google Cloud Storage
index.save_to_url("gs://my-bucket/indexed_peaks/")?;

// Save to Azure Blob Storage
index.save_to_url("az://my-container/indexed_peaks/")?;
```

### Loading from Cloud

```rust
// Load from S3
let lazy_index = LazyIndexedTimstofPeaks::load_from_url(
    "s3://my-bucket/indexed_peaks/"
)?;

// Query works the same as local
let peaks = lazy_index.query_peaks_ms1(mz_range, Unrestricted, Unrestricted);
```

### Authentication

Cloud storage uses default credential chains:

**AWS S3:**
```bash
# Option 1: AWS credentials file
cat ~/.aws/credentials
#[default]
#aws_access_key_id = YOUR_KEY
#aws_secret_access_key = YOUR_SECRET

# Option 2: Environment variables
export AWS_ACCESS_KEY_ID=YOUR_KEY
export AWS_SECRET_ACCESS_KEY=YOUR_SECRET

# Option 3: IAM role (when running on EC2)
```

**Google Cloud Storage:**
```bash
# Option 1: Service account key
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account-key.json

# Option 2: gcloud auth
gcloud auth application-default login
```

**Azure Blob Storage:**
```bash
# Option 1: Connection string
export AZURE_STORAGE_CONNECTION_STRING="DefaultEndpointsProtocol=https;..."

# Option 2: Account key
export AZURE_STORAGE_ACCOUNT=myaccount
export AZURE_STORAGE_KEY=mykey
```

## Async vs Sync API

The library provides both async and blocking APIs. Use async when you're already in an async context:

### Async (Recommended in async contexts)

```rust
use timscentroid::lazy::LazyIndexedTimstofPeaks;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load asynchronously (doesn't block executor thread)
    let lazy_index = LazyIndexedTimstofPeaks::load_from_url_async(
        "s3://my-bucket/indexed_peaks/"
    ).await?;

    // Save asynchronously
    index.save_to_url_async("s3://my-bucket/output/").await?;

    Ok(())
}
```

### Blocking (Simpler for scripts)

```rust
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Blocks current thread (spawns internal runtime)
    let lazy_index = LazyIndexedTimstofPeaks::load_from_url(
        "s3://my-bucket/indexed_peaks/"
    )?;

    // This also blocks
    index.save_to_url("s3://my-bucket/output/")?;

    Ok(())
}
```

## Configuration

### Serialization Settings

Control parquet file characteristics:

```rust
use timscentroid::serialization::SerializationConfig;
use parquet::basic::{Compression, ZstdLevel};

let config = SerializationConfig {
    row_group_size: 100_000,  // Peaks per row group
    compression: Compression::ZSTD(ZstdLevel::try_new(3)?),
};

index.save_to_directory_with_config("./indexed_peaks", config)?;
```

**Row group size recommendations:**
- Small datasets (< 1M peaks): 10k-50k per row group
- Medium datasets (1M-10M peaks): 50k-200k per row group
- Large datasets (> 10M peaks): 200k-1M per row group

Smaller row groups = more granular queries but more overhead.

### Centroiding Settings

```rust
use timscentroid::CentroidingConfig;

let config = CentroidingConfig {
    min_intensity: 100.0,
    min_peaks_per_group: 3,
    ..Default::default()
};

let index = IndexedTimstofPeaks::from_timstof_file(&file, config);
```

## Working with S3-Compatible Services

Digital Ocean Spaces, MinIO, Wasabi, and other S3-compatible services work via environment variables:

```bash
# Set the custom endpoint
export AWS_ENDPOINT_URL="https://nyc3.digitaloceanspaces.com"
export AWS_ACCESS_KEY_ID="your-spaces-key"
export AWS_SECRET_ACCESS_KEY="your-spaces-secret"
export AWS_REGION="nyc3"
```

```rust
// Use s3:// URL - will connect to custom endpoint
let index = LazyIndexedTimstofPeaks::load_from_url("s3://my-space/data/")?;
```

## Performance Notes

**Lazy loading initialization:**
- Local: ~20-50ms (reads metadata.json only)
- Cloud: ~100-150ms (includes network round-trip)

**First query (cold):**
- Local: 20-50ms per row group
- Cloud: ~200ms for concurrent row group fetching

**Concurrent row group fetching:**
The library automatically fetches multiple parquet row groups in parallel when querying cloud storage, providing ~5x speedup compared to sequential fetching.

## Examples

### Query with Filters

```rust
use timscentroid::utils::{TupleRange, OptionallyRestricted::*};

let mz_range = TupleRange::try_new(400.0, 500.0)?;
let rt_range = TupleRange::try_new(100, 200)?;  // Cycle indices
let im_range = TupleRange::try_new(
    half::f16::from_f32(0.8),
    half::f16::from_f32(1.2)
)?;

let peaks: Vec<_> = lazy_index
    .query_peaks_ms1(
        mz_range,
        Restricted(rt_range),
        Restricted(im_range)
    )
    .collect();

println!("Found {} peaks", peaks.len());
```

### Multiple Queries (Reuse Index)

```rust
// Load once
let lazy_index = LazyIndexedTimstofPeaks::load_from_directory("./indexed_peaks")?;

// Query multiple times (efficient)
for window_start in (400..=1000).step_by(100) {
    let mz_range = TupleRange::try_new(window_start as f32, (window_start + 100) as f32)?;
    let peaks: Vec<_> = lazy_index
        .query_peaks_ms1(mz_range, Unrestricted, Unrestricted)
        .collect();
    println!("m/z {}-{}: {} peaks", window_start, window_start + 100, peaks.len());
}
```

## Troubleshooting

### "Unsupported URL scheme" Error

Enable the appropriate feature flag in `Cargo.toml`:

```toml
timscentroid = { version = "0.21", features = ["aws"] }
```

### Authentication Errors

Verify credentials are configured:

```bash
# AWS
aws s3 ls s3://my-bucket/

# GCP
gcloud auth application-default login

# Azure
az storage blob list --account-name myaccount --container-name mycontainer
```

### Slow Cloud Performance

- Ensure you're in the same region as your data
- Check row group size configuration
- Monitor cloud provider logs for rate limiting
- Consider using a cloud VM near your data

## Architecture

**Storage abstraction:** Uses the `object_store` crate for unified local/cloud access

**Async runtime:** Single global Tokio current-thread runtime (~1.5MB), created lazily on first use

**API design:** Async methods are primary, blocking wrappers provided for convenience

**Zero-cost abstraction:** No runtime overhead for storage operations, feature flags only control URL parsing
