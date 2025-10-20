use calibrt::{
    Point,
    calibrate,
};
use rand::Rng;
use tracing_subscriber;

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
        points.push(Point { x, y, weight: 1.0 });
    }
    // Add some random noise points
    let mut rng = rand::thread_rng();
    for _ in 0..50 {
        points.push(Point {
            x: rng.gen_range(0.0..100.0),
            y: rng.gen_range(0.0..150.0),
            weight: 1.0,
        });
    }

    println!("Generated {} initial data points.", points.len());

    // Pre-calculate the ranges, as required by the `calibrate` function signature
    let min_x = points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let max_x = points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let max_y = points.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

    // 2. Define calibration parameters
    let grid_size = 100;

    // 3. Run the Calib-RT algorithm
    match calibrate(&points, (min_x, max_x), (min_y, max_y), grid_size) {
        Ok(calibration_curve) => {
            println!("Calibration successful!");
            println!("CalibrationCurve: {:#?}", calibration_curve);
            // 4. Use the calibration curve
            let test_x_vals: [f64; 6] = [0.0, 25.5, 50.0, 99.0, -10.0, 110.0];
            println!("\n--- Predictions ---");
            for &x in &test_x_vals {
                let real_expect = real_x_to_y(x);
                match calibration_curve.predict(x) {
                    Ok(predicted_y) => println!(
                        "- For library RT {:.2}, predicted measured RT is {:.2}; expect {}",
                        x, predicted_y, real_expect
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
    let empty_points: Vec<Point> = vec![];
    match calibrate(&empty_points, (0.0, 100.0), (0.0, 100.0), grid_size) {
        Ok(_) => println!("This should have failed!"),
        Err(e) => eprintln!("Correctly failed with error: {:?}", e),
    }

    // Example with zero range
    println!("\n--- Testing error case (zero range) ---");
    match calibrate(&points, (50.0, 50.0), (0.0, 100.0), grid_size) {
        Ok(_) => println!("This should have failed!"),
        Err(e) => eprintln!("Correctly failed with error: {:?}", e),
    }
}
