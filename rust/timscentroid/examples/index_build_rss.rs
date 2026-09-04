//! Build a full in-memory index with the default config and print the peak
//! counts. Run under `/usr/bin/time -l` to get peak memory footprint.
//!
//! Usage: cargo run --release --example index_build_rss -- <file.d>

use timscentroid::{
    IndexedTimstofPeaks,
    IndexingCentroidingConfig,
};
use timsrust::TimsTofPath;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: index_build_rss <file.d>");
    let file = TimsTofPath::new(&path).unwrap();
    let cfg = IndexingCentroidingConfig::default();
    println!("file: {path}");

    let st = std::time::Instant::now();
    let (index, stats) = IndexedTimstofPeaks::from_timstof_file(&file, cfg);
    println!("\nbuild time: {:.1?}", st.elapsed());
    println!("{stats}");
    index.print_glimpse();
}
