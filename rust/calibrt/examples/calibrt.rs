use calibrt::{
    CalibrationState,
    LibraryRT,
};
use rand::Rng;

fn setup_tracing() {
    let _ = tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .with_thread_ids(true)
        .with_thread_names(true)
        .try_init();
}

fn main() {
    println!("Running Calib-RT Example");
    setup_tracing();

    // 1. Generate some sample data
    let mut rng = rand::thread_rng();
    let mut real_x_to_y = |x| -> f64 { x + 10. + rng.gen_range(-5.0..5.0) };
    let mut points = Vec::new();
    // Generate some points along a rough line y = x + 10
    for i in 0..500 {
        let x = i as f64;
        let y = real_x_to_y(x);
        points.push((x, y));
    }
    // Add some random noise points
    let mut rng = rand::thread_rng();
    for _ in 0..50 {
        points.push((rng.gen_range(0.0..100.0), rng.gen_range(0.0..150.0)));
    }

    println!("Generated {} initial data points.", points.len());

    // 2. Define calibration parameters
    let grid_size = 100;

    // 3. Run the Calib-RT algorithm (ranges are derived from the points)
    let mut state = CalibrationState::deferred(grid_size, 30).unwrap();
    match state.refit(grid_size, points.iter().copied(), &mut ()) {
        Ok(ranges) => {
            println!("Calibration successful over {ranges:?}!");
            let calibration_curve = state.curve().expect("a successful fit has a curve");
            println!("CalibrationCurve: {:#?}", calibration_curve);
            // 4. Use the calibration curve
            let test_x_vals: [f64; 6] = [0.0, 25.5, 50.0, 99.0, -10.0, 110.0];
            println!("\n--- Predictions ---");
            for &x in &test_x_vals {
                let real_expect = real_x_to_y(x);
                match calibration_curve.predict(LibraryRT(x)) {
                    Ok(predicted_y) => println!(
                        "- For library RT {:.2}, predicted measured RT is {:.2}; expect {}",
                        x, predicted_y.0, real_expect
                    ),
                    Err(e) => eprintln!(
                        "- For library RT {:.2}, prediction failed: {:?} expected: {}",
                        x, e, real_expect
                    ),
                }
            }
        }
        Err(e) => {
            eprintln!("Calibration failed: {:?}", e);
        }
    }

    // Example with no points
    println!("\n--- Testing error case (no points) ---");
    match state.refit(grid_size, std::iter::empty(), &mut ()) {
        Ok(_) => println!("This should have failed!"),
        Err(e) => eprintln!("Correctly failed with error: {:?}", e),
    }
}
