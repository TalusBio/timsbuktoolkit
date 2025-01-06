
1. Redefine the elution group to separate the expected intensities and the actual query.
2. Split scoring + cosine sim into timsseek (out of timsquery)
3. Separate scoring into mz-major vs rt-major.
4. Make the elution group less generic. (single use abstraction for a query)

mz major (all elements of the same mz are contiguous in mem) scores:
- mz errors
- mobility errors
- 

rt major (all elements of the same rt are contiguous in mem) scores:
- cosine similarity
- npeaks
- summed intensity

pub struct ScoresAtTime {
    /// Gen 0
    pub retention_time_miliseconds: u32,
    pub transition_intensities: Vec<u64>,

    /// Gen 1
    // RT major
    pub lazyerscore: f64,
    pub lazyerscore_vs_baseline: f64,
    pub npeaks: u8,
    pub average_mobility: f64,
    pub summed_intensity: u64,
    pub cosine_similarity: f64,

    // mz major
    pub mz_errors: Vec<f64>,
    pub mobility_errors: Vec<f64>,

    /// Gen 2
    pub norm_lazyerscore_vs_baseline: f64,
}