use calibrt::{
    LibraryRT,
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
    let predicted = curve.predict(LibraryRT(50.0)).unwrap();
    assert!((predicted.0 - 60.0).abs() < 5.0); // Allow small error
}

#[test]
fn test_calibrate_empty_points() {
    // Test: Empty input should return error
    let empty: Vec<Point> = vec![];
    let result = calibrate(&empty, 50);
    assert!(result.is_err());
}

/// A degenerate range on *either* axis is an error. Both, because the two are
/// separate checks in `grid::spans` — dropping the y one leaves the x case
/// green.
#[test]
fn test_calibrate_zero_range_on_either_axis() {
    let points = vec![Point {
        library: 50.0,
        observed: 60.0,
        weight: 1.0,
    }];
    for (x_range, y_range) in [((50.0, 50.0), (0.0, 100.0)), ((0.0, 100.0), (60.0, 60.0))] {
        let result = calibrate_with_ranges(&points, x_range, y_range, 50, 30);
        assert!(
            result.is_err(),
            "x {x_range:?} y {y_range:?} must be rejected"
        );
    }
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
    let result = curve.predict(LibraryRT(100.0));
    assert!(result.is_err());
}
