use rand::{
    Rng,
    SeedableRng,
};
use rand_chacha::ChaCha8Rng;
use serde::Serialize;
use std::env;
use std::fs::File;
use std::path::{
    Path,
    PathBuf,
};
use std::sync::Arc;
use std::time::{
    Duration,
    Instant,
};
use timsquery::models::aggregators::PointIntensityAggregator;
use timsquery::models::indices::expanded_raw_index::ExpandedRawFrameIndex;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::models::tolerance::{
    MobilityTolerance,
    MzTolerance,
    QuadTolerance,
    RtTolerance,
};
use timsquery::traits::QueriableData;
use timsquery::{
    ElutionGroup,
    Tolerance,
};
use tracing::subscriber::set_global_default;
use tracing_bunyan_formatter::{
    BunyanFormattingLayer,
    JsonStorageLayer,
};
use tracing_chrome::ChromeLayerBuilder;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

const NUM_ELUTION_GROUPS: usize = 1000;
const NUM_ITERATIONS: usize = 1;

#[derive(Debug, Clone)]
struct EnvConfig {
    tims_data_file: PathBuf,
    file_stem: String,
    skip_highmem: bool,
    skip_slow: bool,
    skip_build: bool,
    skip_query: bool,
}

impl EnvConfig {
    fn new() -> Self {
        let mut tims_data_file = String::new();
        let mut skip_highmem = false;
        // let mut skip_slow = false;
        let mut skip_slow = true;
        let mut skip_build = false;
        let mut skip_query = false;

        for (key, value) in env::vars() {
            match key.as_str() {
                "TIMS_DATA_FILE" => tims_data_file = value,
                "SKIP_HIGHMEM" => skip_highmem = value == "1",
                // "SKIP_SLOW" => skip_slow = value == "1",
                // I hate this double negative but I want to make slow benches
                // opt-in rather than opt-out.
                "RUN_SLOW" => skip_slow = value != "1",
                "SKIP_BUILD" => skip_build = value == "1",
                "SKIP_QUERY" => skip_query = value == "1",
                _ => {}
            }
        }

        if tims_data_file.is_empty() {
            panic!("TIMS_DATA_FILE environment variable not set");
        }

        // Check that the file exists and is readable

        let tims_data_file = Path::new(&tims_data_file);
        if !tims_data_file.exists() {
            panic!("TIMS_DATA_FILE={} does not exist", tims_data_file.display());
        }
        let file_stem = tims_data_file.file_stem().unwrap().to_str().unwrap();

        let out = Self {
            tims_data_file: tims_data_file.to_path_buf(),
            file_stem: file_stem.to_string(),
            skip_highmem,
            skip_slow,
            skip_build,
            skip_query,
        };

        println!("Env config: {:#?}", out);

        out
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BenchmarkTag {
    HighMem,
    Slow,
    Build,
    Query,
}

#[derive(Debug, Serialize)]
struct BenchmarkIteration {
    iteration: usize,
    duration_seconds: f64,
}

#[derive(Debug, Serialize)]
struct BenchmarkResult {
    name: String,
    context: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    iterations: Option<Vec<BenchmarkIteration>>,
    mean_duration_seconds: f64,
    mean_duration_human_readable: String,
    setup_time_seconds: f64,
    setup_time_human_readable: String,
    note: Option<String>,
}

#[derive(Debug, Serialize)]
struct BenchmarkReport {
    settings: BenchmarkSettings,
    results: Vec<BenchmarkResult>,
    metadata: BenchmarkMetadata,
    full_benchmark_time_seconds: f64,
}

#[derive(Debug, Serialize)]
struct BenchmarkSettings {
    num_elution_groups: usize,
    num_iterations: usize,
}

#[derive(Debug, Serialize)]
struct BenchmarkMetadata {
    basename: String,
}

fn duration_to_seconds(duration: Duration) -> f64 {
    duration.as_secs() as f64 + duration.subsec_nanos() as f64 * 1e-9
}

fn build_elution_groups() -> Vec<Arc<ElutionGroup<u64>>> {
    const NUM_FRAGMENTS: usize = 10;
    const MAX_RT: f32 = 22.0 * 60.0;
    const MAX_MOBILITY: f32 = 1.5;
    const MIN_MOBILITY: f32 = 0.5;
    const MAX_MZ: f64 = 1000.0;
    const MIN_MZ: f64 = 300.0;

    let mut out_egs = Vec::with_capacity(NUM_ELUTION_GROUPS);
    let mut rng = ChaCha8Rng::seed_from_u64(43u64);

    for i in 1..NUM_ELUTION_GROUPS {
        let rt = rng.gen::<f32>() * MAX_RT;
        let mobility = rng.gen::<f32>() * (MAX_MOBILITY - MIN_MOBILITY) + MIN_MOBILITY;
        let mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;

        let mut fragments = Vec::with_capacity(NUM_FRAGMENTS);
        for ii in 0..NUM_FRAGMENTS {
            let fragment_mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;
            fragments.push((ii as u64, fragment_mz));
        }

        out_egs.push(Arc::new(ElutionGroup {
            id: i as u64,
            rt_seconds: rt,
            mobility,
            precursors: Arc::from(vec![(0, mz)]),
            fragments: Arc::from(fragments),
        }));
    }
    out_egs
}

impl EnvConfig {
    fn with_benchmark<T>(
        &self,
        name: &str,
        context: &str,
        setup_fn: impl Fn() -> T,
        f: impl Fn(&T, usize) -> Option<String>,
        tags: &[BenchmarkTag],
    ) -> Option<BenchmarkResult> {
        if self.skip_slow && tags.contains(&BenchmarkTag::Slow) {
            println!("Skipping slow benchmark for {} {}", name, context);
            return None;
        }

        if self.skip_highmem && tags.contains(&BenchmarkTag::HighMem) {
            println!("Skipping highmem benchmark {} {}", name, context);
            return None;
        }

        if self.skip_build && tags.contains(&BenchmarkTag::Build) {
            println!("Skipping build benchmark {} {}", name, context);
            return None;
        }

        if self.skip_query && tags.contains(&BenchmarkTag::Query) {
            println!("Skipping query benchmark {} {}", name, context);
            return None;
        }

        println!("{} {} Tagged with {:?}", name, context, tags);
        let sst = Instant::now();
        println!("Setting up {} {}", name, context);
        let setup = setup_fn();
        let setup_elap = sst.elapsed();
        println!("Setup took {:#?}", setup_elap);

        let mut iterations = Vec::with_capacity(NUM_ITERATIONS);
        let mut durations = Vec::with_capacity(NUM_ITERATIONS);
        let mut extras = None;
        let glob_start = Instant::now();
        let mut i = 0;
        while (glob_start.elapsed().as_millis() < 1000) || (i < NUM_ITERATIONS) {
            if i < 10 || i % 100 == 0 {
                println!("{name} iteration {i}");
            }
            let start = Instant::now();
            if let Some(out) = f(&setup, i) {
                if i < 3 {
                    println!("{}", out);
                }
                extras = Some(out);
            }

            let duration = start.elapsed();

            iterations.push(BenchmarkIteration {
                iteration: i + 1,
                duration_seconds: duration_to_seconds(duration),
            });
            durations.push(duration);
            i += 1;
        }

        let mean = iterations.iter().map(|i| i.duration_seconds).sum::<f64>() / (i as f64);
        let avg_duration: Duration = durations.iter().sum::<Duration>() / i as u32;
        let mean_duration_human_readable = format!("{:?}", avg_duration);

        println!("Mean duration: {:?}", mean_duration_human_readable);

        Some(BenchmarkResult {
            name: name.to_string(),
            context: context.to_string(),
            iterations: Some(iterations),
            mean_duration_seconds: mean,
            mean_duration_human_readable,
            setup_time_seconds: setup_elap.as_secs_f64(),
            setup_time_human_readable: format!("{:?}", setup_elap),
            note: extras,
        })
    }
}

fn run_encoding_benchmark(raw_file_path: &Path, env_config: EnvConfig) -> Vec<BenchmarkResult> {
    let raw_file_path = raw_file_path.to_str().unwrap();
    let mut out = vec![];
    let erfic = env_config.with_benchmark(
        "ExpandedRawFileIndexCentroided",
        "Encoding",
        || {},
        |_, _| {
            ExpandedRawFrameIndex::from_path_centroided(raw_file_path).unwrap();
            None
        },
        &[BenchmarkTag::Build],
    );
    out.push(erfic);

    let erfi = env_config.with_benchmark(
        "ExpandedRawFileIndex",
        "Encoding",
        || {},
        |_, _| {
            ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
            None
        },
        &[BenchmarkTag::HighMem, BenchmarkTag::Build],
    );
    out.push(erfi);
    let tqi = env_config.with_benchmark(
        "TransposedQuadIndex",
        "Encoding",
        || {},
        |_, _| {
            QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap();
            None
        },
        &[BenchmarkTag::HighMem, BenchmarkTag::Build],
    );
    out.push(tqi);
    let tqic = env_config.with_benchmark(
        "TransposedQuadIndexCentroided",
        "Encoding",
        || {},
        |_, _| {
            QuadSplittedTransposedIndex::from_path_centroided(raw_file_path).unwrap();
            None
        },
        &[BenchmarkTag::Build],
    );
    out.push(tqic);

    out.into_iter().flatten().collect()
}

fn run_batch_access_benchmark(raw_file_path: &Path, env_config: EnvConfig) -> Vec<BenchmarkResult> {
    let raw_file_path = raw_file_path.to_str().unwrap();
    let mut out = vec![];
    let query_groups = build_elution_groups();
    let tolerance_with_rt = Tolerance {
        ms: MzTolerance::Ppm((20.0, 20.0)),
        rt: RtTolerance::Minutes((5.0, 5.0)),
        mobility: MobilityTolerance::Pct((3.0, 3.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    };
    let tolerance_with_nort = Tolerance {
        ms: MzTolerance::Ppm((20.0, 20.0)),
        rt: RtTolerance::None,
        mobility: MobilityTolerance::Pct((3.0, 3.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    };
    let tolerances = [
        (tolerance_with_rt, "narrow_rt"),
        (tolerance_with_nort, "full_rt"),
    ];

    // TODO: Refactor this ... there is a lot of code duplication but over different types ... int
    // ion thoty this is a good place for a macro.

    for (tolerance, tol_name) in tolerances.clone() {
        let erfi = env_config.with_benchmark(
            "ExpandedRawFileIndex",
            &format!("BatchAccess_{}", tol_name),
            || ExpandedRawFrameIndex::from_path(raw_file_path).unwrap(),
            |index, _i| {
                let mut qa: Vec<_> = query_groups
                    .iter()
                    .map(|x| PointIntensityAggregator::new_with_elution_group(x.clone()))
                    .collect();
                index.par_add_query_multi(&mut qa, &tolerance);
                let tot: f64 = qa.into_iter().map(|x| x.intensity).sum();
                let out = format!(
                    "ExpandedRawFileIndex::query_multi_group aggregated {} ",
                    tot,
                );
                Some(out)
            },
            &[BenchmarkTag::HighMem, BenchmarkTag::Query],
        );
        out.push(erfi);
    }

    for (tolerance, tol_name) in tolerances.clone() {
        let erfic = env_config.with_benchmark(
            "ExpandedRawFileIndexCentroided",
            &format!("BatchAccess_{}", tol_name),
            || ExpandedRawFrameIndex::from_path_centroided(raw_file_path).unwrap(),
            |index, _i| {
                let mut tmp: Vec<_> = query_groups
                    .iter()
                    .map(|x| PointIntensityAggregator::new_with_elution_group(x.clone()))
                    .collect();
                index.par_add_query_multi(&mut tmp, &tolerance);
                let tot: f64 = tmp.into_iter().map(|x| x.intensity).sum();
                let out = format!(
                    "ExpandedRawFileIndexCentroided::query_multi_group aggregated {} ",
                    tot,
                );
                Some(out)
            },
            &[BenchmarkTag::Query],
        );
        out.push(erfic);
    }

    for (tolerance, tol_name) in tolerances.clone() {
        let tqi = env_config.with_benchmark(
            "TransposedQuadIndex",
            &format!("BatchAccess_{}", tol_name),
            || QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap(),
            |index, _i| {
                let mut tmp: Vec<_> = query_groups
                    .iter()
                    .map(|x| PointIntensityAggregator::new_with_elution_group(x.clone()))
                    .collect();
                index.par_add_query_multi(&mut tmp, &tolerance);
                let tot: f64 = tmp.into_iter().map(|x| x.intensity).sum();
                let out = format!("TransposedQuadIndex::query_multi_group aggregated {} ", tot,);
                Some(out)
            },
            &[BenchmarkTag::HighMem, BenchmarkTag::Query],
        );
        out.push(tqi);
    }

    for (tolerance, tol_name) in tolerances.clone() {
        let tqi = env_config.with_benchmark(
            "TransposedQuadIndexCentroided",
            &format!("BatchAccess_{}", tol_name),
            || QuadSplittedTransposedIndex::from_path_centroided(raw_file_path).unwrap(),
            |index, _i| {
                let mut queriable_aggregators: Vec<_> = query_groups
                    .iter()
                    .map(|x| PointIntensityAggregator::new_with_elution_group(x.clone()))
                    .collect();
                index.par_add_query_multi(&mut queriable_aggregators, &tolerance);
                let tot: f64 = queriable_aggregators.iter().map(|x| x.intensity).sum();
                let out = format!("TransposedQuadIndex::query_multi_group aggregated {} ", tot,);
                Some(out)
            },
            &[BenchmarkTag::Query],
        );
        out.push(tqi);
    }
    out.into_iter().flatten().collect()
}

fn write_results(
    mut report: BenchmarkReport,
    basename: &str,
    parent: &Path,
) -> std::io::Result<(PathBuf, String)> {
    let filepath = parent.join(format!("benchmark_results_{}.json", basename));
    let file = File::create(&filepath)?;
    serde_json::to_writer_pretty(file, &report)?;
    let _out = serde_json::to_string_pretty(&report)?;

    for result in report.results.iter_mut() {
        result.iterations = None;
    }
    let slim_filepath = parent.join(format!("slim_benchmark_results_{}.json", basename));
    let slim_file = File::create(&slim_filepath)?;
    serde_json::to_writer_pretty(slim_file, &report)?;
    let slim_out = serde_json::to_string_pretty(&report)?;

    Ok((filepath, slim_out))
}

fn main() {
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    let formatting_layer = BunyanFormattingLayer::new("timsquery".into(), std::io::stdout);
    let (chrome_layer, _guard) = ChromeLayerBuilder::new().build();

    let subscriber = Registry::default()
        .with(env_filter)
        .with(chrome_layer)
        .with(JsonStorageLayer)
        .with(formatting_layer);

    set_global_default(subscriber).expect("Setting default subscriber failed");
    let st = Instant::now();

    let env_config = EnvConfig::new();
    let raw_file_path = env_config.tims_data_file.clone();
    let basename = env_config.file_stem.clone();

    let file_parent = std::path::Path::new(&raw_file_path)
        .parent()
        .expect("Expected to find a parent directory");
    let mut all_results: Vec<BenchmarkResult> = vec![];

    println!(
        "Run Settings: Elution groups={}, Iterations={}, basename={}",
        NUM_ELUTION_GROUPS, NUM_ITERATIONS, basename
    );

    println!("Running encoding benchmarks...");
    let encoding_results = run_encoding_benchmark(&raw_file_path, env_config.clone());
    all_results.extend(encoding_results);

    println!("Running batch access benchmarks...");
    let batch_results = run_batch_access_benchmark(&raw_file_path, env_config.clone());
    all_results.extend(batch_results);
    let report = BenchmarkReport {
        settings: BenchmarkSettings {
            num_elution_groups: NUM_ELUTION_GROUPS,
            num_iterations: NUM_ITERATIONS,
        },
        results: all_results,
        metadata: BenchmarkMetadata {
            basename: basename.to_string(),
        },
        full_benchmark_time_seconds: st.elapsed().as_secs_f64(),
    };

    match write_results(report, &basename, file_parent) {
        Ok((fp, res)) => println!("Results written to {} \n{}", fp.to_string_lossy(), res),
        Err(e) => eprintln!("Error writing results: {}", e),
    }
}
