use std::env;
use std::time::Instant;
use timsquery::models::frames::expanded_frame::{
    FrameProcessingConfig,
    par_expand_and_arrange_frames,
    par_expand_and_centroid_frames,
    warn_and_skip_badframes,
};
use timsrust::MSLevel;
use timsrust::readers::{
    FrameReader,
    MetadataReader,
};
use tracing::subscriber::set_global_default;
use tracing_bunyan_formatter::{
    BunyanFormattingLayer,
    JsonStorageLayer,
};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::Registry;

fn main() {
    // "TIMS_DATA_FILE"
    // Read from env
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    let formatting_layer = BunyanFormattingLayer::new("timsquery".into(), std::io::stdout);
    let subscriber = Registry::default()
        .with(env_filter)
        .with(JsonStorageLayer)
        .with(formatting_layer);

    set_global_default(subscriber).expect("Setting default subscriber failed");

    let tims_data_file = env::var("TIMS_DATA_FILE").expect("TIMS_DATA_FILE not set");
    let sql_path = std::path::Path::new(&tims_data_file).join("analysis.tdf");
    let meta_converters = MetadataReader::new(&sql_path).unwrap();
    let centroiding_config = FrameProcessingConfig::default_centroided();
    let centroiding_config = centroiding_config
        .with_converters(meta_converters.im_converter, meta_converters.mz_converter);

    let frame_reader = FrameReader::new(&tims_data_file).unwrap();
    let ms1_iter = frame_reader.parallel_filter(|x| x.ms_level == MSLevel::MS1);
    let ms1_iter = warn_and_skip_badframes(ms1_iter);
    let start = Instant::now();
    println!("Starting expansion");
    let _expanded_ms1_frames = match centroiding_config {
        FrameProcessingConfig::Centroided {
            settings,
            ims_converter,
            mz_converter,
        } => par_expand_and_centroid_frames(
            ms1_iter,
            settings.ims_tol_pct,
            settings.mz_tol_ppm,
            settings.window_width,
            &ims_converter.unwrap(),
            &mz_converter.unwrap(),
        ),
        FrameProcessingConfig::NotCentroided => par_expand_and_arrange_frames(ms1_iter),
    };
    let elapsed = start.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}
