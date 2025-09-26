use bincode;
use serde::Serialize;
use std::fs::File;
use std::path::Path;
use timsquery::{CentroidingConfig, IndexedTimstofPeaks, TimsTofPath};
use tracing::{error, info};
use zstd::stream::read::Decoder;
use zstd::stream::write::Encoder;

// Save with compression
fn save_compressed<T: Serialize>(data: &T, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let mut encoder = Encoder::new(file, 3)?;
    bincode::serialize_into(&mut encoder, data)?;
    encoder.finish()?;
    Ok(())
}

// Load with decompression
fn load_compressed<T: serde::de::DeserializeOwned>(
    path: impl AsRef<Path>,
) -> Result<T, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let decoder = Decoder::new(file)?;
    let data = bincode::deserialize_from(decoder)?;
    Ok(data)
}

fn maybe_cache_load_index(index_cache_loc: impl AsRef<Path>) -> Option<IndexedTimstofPeaks> {
    info!(
        "Attempting to load index from cache at {:?}",
        index_cache_loc.as_ref()
    );
    match load_compressed(index_cache_loc.as_ref()) {
        Ok(idx) => {
            info!("Loaded index from cache at {:?}", index_cache_loc.as_ref());
            Some(idx)
        }
        Err(e) => {
            error!(
                "Failed to load index from cache at {:?}: {:?}",
                index_cache_loc.as_ref(),
                e
            );
            None
        }
    }
}

fn uncached_load_index(
    timstofpath: &TimsTofPath,
    cache_loc: &Option<std::path::PathBuf>,
) -> IndexedTimstofPeaks {
    let centroiding_config = CentroidingConfig {
        max_peaks: 50_000,
        mz_ppm_tol: 10.0,
        im_pct_tol: 5.0,
        early_stop_iterations: 200,
    };
    info!("Using centroiding config: {:#?}", centroiding_config);
    info!("Starting centroiging + load of the raw data (might take a min)");
    let (index, build_stats) =
        IndexedTimstofPeaks::from_timstof_file(&timstofpath, centroiding_config);
    info!("Index built with stats: {}", build_stats);

    // Save to cache
    if let Some(ref idx_path) = cache_loc {
        info!("Saving index to cache at {:?}", idx_path);
        if let Err(e) = save_compressed(&index, idx_path.to_str().unwrap()) {
            error!("Failed to save index to cache: {:?}", e);
        } else {
            info!("Saved index to cache");
        }
    }
    index
}

pub fn load_index_caching(
    file_location: impl AsRef<Path>,
) -> Result<IndexedTimstofPeaks, crate::errors::TimsSeekError> {
    let st = std::time::Instant::now();
    let timstofpath = TimsTofPath::new(file_location.as_ref())?;
    // TODO: make timstofpath expose the path so we can use that as a prefix ...
    // Since that would be the "cannonical" path
    let index_location = file_location
        .as_ref()
        .to_path_buf()
        .with_extension("idx.zst");

    let out = if let Some(idx) = maybe_cache_load_index(&index_location) {
        Ok(idx)
    } else {
        Ok(uncached_load_index(&timstofpath, &Some(index_location)))
    };
    let et = st.elapsed();
    info!("Loading index took: {:#?}", et);
    out
}
