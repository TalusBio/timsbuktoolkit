use calibrt::{
    CalibRtError,
    CalibrationState,
    LibraryRT,
};

/// A degenerate range on *either* axis is an error. Both, because the two are
/// separate checks in `grid::spans` — dropping the y one leaves the x case
/// green.
#[test]
fn a_zero_range_on_either_axis_is_rejected() {
    for (x_range, y_range) in [((50.0, 50.0), (0.0, 100.0)), ((0.0, 100.0), (60.0, 60.0))] {
        let result = CalibrationState::new(50, x_range, y_range, 30);
        assert!(
            matches!(result, Err(CalibRtError::ZeroRange)),
            "x {x_range:?} y {y_range:?} must be rejected"
        );
    }
}

#[test]
fn predicting_outside_the_calibrated_range_is_out_of_bounds() {
    let mut state = CalibrationState::deferred(30, 30).unwrap();
    state
        .refit(30, (0..50).map(|i| (i as f64, i as f64)))
        .unwrap();

    let curve = state.curve().expect("the diagonal must fit");
    assert!((curve.predict(LibraryRT(25.0)).unwrap().0 - 25.0).abs() < 5.0);
    assert!(matches!(
        curve.predict(LibraryRT(100.0)),
        Err(CalibRtError::OutOfBounds(_))
    ));
}
