use timscentroid::CentroidingConfig;
use timscentroid::centroiding::{
    ConvertableDomain,
    PeakCentroider,
};
use timsrust::TimsTofPath;
use timsrust::core::Frame;
use timsrust_calibration::{
    CalibratedScan2ImConverter,
    CalibratedTof2MzConverter,
    RunCalibration,
};

const MAX_PEAKS: usize = 20_000;
const DIA_TEST: &str =
    "/Users/sebastianpaez/data/decompressed_timstof/250225_Desnaux_200ng_Hela_ICC_on_DIA.d/";

fn converted_input_peaks(
    frame: &Frame,
    mz_converter: &CalibratedTof2MzConverter,
    im_converter: &CalibratedScan2ImConverter,
) -> Vec<(f32, f32, f32)> {
    let factor = frame.info().intensity_correction_factor();
    let mut out = Vec::with_capacity(frame.len());
    for scan in 0..frame.ions().scan_count() {
        for (tof, intensity) in frame.ions().read_scan(scan) {
            out.push((
                ConvertableDomain::convert(mz_converter, u32::from(tof) as f64) as f32,
                ConvertableDomain::convert(im_converter, scan as f64) as f32,
                (factor * u32::from(intensity) as f64) as f32,
            ));
        }
    }
    out
}

fn time_itercentroid(
    frame: &Frame,
    centroider: &mut PeakCentroider<CalibratedTof2MzConverter, CalibratedScan2ImConverter>,
    mz_converter: &CalibratedTof2MzConverter,
    im_converter: &CalibratedScan2ImConverter,
) -> usize {
    let inp_start = converted_input_peaks(frame, mz_converter, im_converter);
    let start = std::time::Instant::now();
    let (reason, peaks_iter) = centroider.centroid_frame(frame);
    let duration = start.elapsed();
    println!("Iter centroiding took: {duration:?}, stop reason: {reason:?}");
    let out_peaks = peaks_iter
        .map(|p| {
            (
                ConvertableDomain::convert(mz_converter, p.tof_index as f64) as f32,
                ConvertableDomain::convert(im_converter, p.scan_index as f64) as f32,
                p.corrected_intensity as f32,
            )
        })
        .collect::<Vec<_>>();

    println!(
        "Writing {} peaks to disk from {:#?}",
        out_peaks.len(),
        frame.info()
    );
    serde_json::to_writer(
        std::fs::File::create("centroided_peaks.json").unwrap(),
        &out_peaks,
    )
    .unwrap();
    serde_json::to_writer(
        std::fs::File::create("input_peaks.json").unwrap(),
        &inp_start,
    )
    .unwrap();
    out_peaks.len()
}

fn main() {
    let file = TimsTofPath::new(DIA_TEST).unwrap();
    let frame_reader = file.frame_reader().unwrap();
    let calibration = RunCalibration::from_path(&file).unwrap();
    let mz_median = calibration.mz_converter_median().unwrap();
    let im_median = calibration.im_converter_median().unwrap();

    let mut centroider = PeakCentroider::with_capacity(
        50_000,
        CentroidingConfig {
            max_peaks: MAX_PEAKS,
            im_pct_tol: 3.0,
            mz_tol: timscentroid::centroiding::MzTolerance::Bins(0),
            early_stop_iterations: 200,
            window_cap: None,
        },
        mz_median,
        im_median,
    );

    let mut frame_indices: Vec<_> = frame_reader.iter_indices().collect();
    frame_indices.sort_unstable();
    let sample_frame = frame_reader.get_frame(frame_indices[2000]).unwrap();
    let mz_converter = calibration
        .mz_converter(sample_frame.info().index())
        .unwrap();
    let im_converter = calibration
        .im_converter(sample_frame.info().index())
        .unwrap();

    println!("Frame has {} peaks", sample_frame.len());
    let count = time_itercentroid(&sample_frame, &mut centroider, &mz_converter, &im_converter);
    println!("CentroidedFrame has {count} peaks");
}
