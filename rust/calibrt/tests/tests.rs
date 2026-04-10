use calibrt::{
    Point,
    calibrate,
    calibrate_with_ranges,
};

#[test]
fn test_calibrate_with_linear_data() {
    // Test: Linear relationship y = x + 10 with slight noise
    let points: Vec<Point> = (0..100)
        .map(|i| Point {
            library: i as f64,
            observed: i as f64 + 10.0,
            weight: 1.0,
        })
        .collect();

    let result = calibrate(&points, 50);
    assert!(result.is_ok());

    let curve = result.unwrap();
    let predicted = curve.predict(50.0).unwrap();
    assert!((predicted - 60.0).abs() < 5.0); // Allow small error
}

#[test]
fn test_calibrate_empty_points() {
    // Test: Empty input should return error
    let empty: Vec<Point> = vec![];
    let result = calibrate(&empty, 50);
    assert!(result.is_err());
}

#[test]
fn test_calibrate_zero_x_range() {
    // Test: Zero x range should return error
    let points = vec![Point {
        library: 50.0,
        observed: 60.0,
        weight: 1.0,
    }];
    let result = calibrate_with_ranges(&points, (50.0, 50.0), (0.0, 100.0), 50, 30);
    assert!(result.is_err());
}

#[test]
fn test_calibrate_zero_y_range() {
    // Test: Zero y range should return error
    let points = vec![Point {
        library: 50.0,
        observed: 60.0,
        weight: 1.0,
    }];
    let result = calibrate_with_ranges(&points, (0.0, 100.0), (60.0, 60.0), 50, 30);
    assert!(result.is_err());
}

#[test]
fn test_predict_within_range() {
    // Test: Prediction within calibration range
    let points: Vec<Point> = (0..50)
        .map(|i| Point {
            library: i as f64,
            observed: i as f64 * 2.0,
            weight: 1.0,
        })
        .collect();

    let curve = calibrate(&points, 30).unwrap();
    let result = curve.predict(25.0);
    assert!(result.is_ok());
}

#[test]
fn test_predict_outside_range() {
    // Test: Prediction outside calibration range
    let points: Vec<Point> = (0..50)
        .map(|i| Point {
            library: i as f64,
            observed: i as f64,
            weight: 1.0,
        })
        .collect();

    let curve = calibrate(&points, 30).unwrap();
    let result = curve.predict(100.0);
    assert!(result.is_err());
}

#[test]
fn test_calibrate_single_point() {
    // Test: Single point should succeed or fail gracefully
    let points = vec![Point {
        library: 5.0,
        observed: 10.0,
        weight: 1.0,
    }];
    let result = calibrate(&points, 10);
    // Verify it doesn't panic
    let _ = result;
}

#[test]
fn test_calibrate_weighted_points() {
    // Test: Different weights affect calibration
    let mut points = vec![
        Point {
            library: 10.0,
            observed: 20.0,
            weight: 10.0,
        },
        Point {
            library: 10.0,
            observed: 30.0,
            weight: 1.0,
        },
    ];
    for i in 0..20 {
        points.push(Point {
            library: i as f64,
            observed: i as f64 * 2.0,
            weight: 1.0,
        });
    }

    let result = calibrate(&points, 20);
    assert!(result.is_ok());
}
