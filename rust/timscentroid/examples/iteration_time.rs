use timscentroid::CentroidingConfig;
use timsrust::converters::{
    ConvertableDomain,
    Tof2MzConverter2,
};
use timsrust::{
    Frame,
    TimsTofPath,
};

use timscentroid::centroiding::PeakCentroider;

const MAX_PEAKS: usize = 20_000;
const DIA_TEST: &str =
    "/Users/sebastianpaez/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_on_DIA.d/";

fn time_itercentroid<T1: ConvertableDomain, T2: ConvertableDomain>(
    frame: &Frame,
    centroider: &mut PeakCentroider<T1, T2>,
    mz_converter: &Tof2MzConverter2,
    ims_convertsion_fun: &impl Fn(f64) -> f64,
) -> usize {
    let inp_start = frame.iter_corrected_peaks().map(|p| {
        let mz = mz_converter.convert(p.tof_index as f64);
        let ims = ims_convertsion_fun(p.scan_index as f64);
        (mz as f32, ims as f32, p.corrected_intensity as f32)
    });
    let start = std::time::Instant::now();
    let (reason, peaks_iter) = centroider.centroid_frame(frame);
    let duration = start.elapsed();
    println!(
        "Iter centroiding took: {:?}, stop reason: {:?}",
        duration, reason
    );
    let out_peaks = peaks_iter
        .map(|p| {
            let mz = mz_converter.convert(p.tof_index as f64);
            let ims = ims_convertsion_fun(p.scan_index as f64);
            (mz as f32, ims as f32, p.corrected_intensity as f32)
        })
        .collect::<Vec<_>>();

    // Write to disk ...
    println!(
        "Writing {} peaks to disk from {:#?}",
        out_peaks.len(),
        frame.meta
    );
    serde_json::to_writer(
        std::fs::File::create("centroided_peaks.json").unwrap(),
        &out_peaks,
    )
    .map_err(|e| {
        eprintln!("Error writing centroided peaks: {}", e);
        e
    })
    .ok();
    serde_json::to_writer(
        std::fs::File::create("input_peaks.json").unwrap(),
        &inp_start.collect::<Vec<_>>(),
    )
    .map_err(|e| {
        eprintln!("Error writing input peaks: {}", e);
        e
    })
    .ok();
    let count = out_peaks.len();
    println!("Done");
    count
}

fn main() {
    let file = TimsTofPath::new(DIA_TEST).unwrap();
    let frame_reader = file.load_frame_reader().unwrap();
    let metadata = file.load_metadata().unwrap();
    let (mz_calib, tims_calib) = metadata.get_calibration().unwrap();

    let ims_convertsion_fun = tims_calib.get_by_id(1).unwrap().get_conversion_function();

    let mut centroider = PeakCentroider::with_capacity(
        50_000,
        CentroidingConfig {
            max_peaks: MAX_PEAKS,
            im_pct_tol: 3.0,
            mz_ppm_tol: 5.0,
            early_stop_iterations: 200,
        },
        metadata.mz_converter,
        metadata.im_converter,
    );

    // This is an ms2 in my file
    let sample_frame = frame_reader.get_by_internal_index(2000).unwrap();
    // In my file this is an MS1
    // let sample_frame = frame_reader.get_by_internal_index(2007).unwrap();
    let mz_converter = Tof2MzConverter2::try_from_calibration(
        mz_calib
            .get_by_id(sample_frame.meta.calibration.calibration_id)
            .unwrap(),
        sample_frame.meta.calibration.t1,
        sample_frame.meta.calibration.t2,
    )
    .unwrap(); // TODO: make this an error instead of an option...

    println!("Frame has {} peaks", sample_frame.peaks.len());
    let count = time_itercentroid(
        &sample_frame,
        &mut centroider,
        &mz_converter,
        &ims_convertsion_fun,
    );
    println!("CentroidedFrame has {} peaks", count);
}
