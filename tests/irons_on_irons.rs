use mag_calib::Irons;
use nalgebra::{Matrix3, Vector3};

#[test]
fn test_irons_on_irons_composition() {
    // Synthetic test point
    let point = [1.0, 2.0, 3.0];

    // First Irons (soft1 and hard1)
    let soft1 = Matrix3::new(1.1, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 1.05);
    let hard1 = Vector3::new(0.2, -0.3, 0.1);
    let irons1 = Irons {
        soft: soft1,
        hard: hard1,
    };

    // Second Irons, to be applied in already calibrated space
    let soft2 = Matrix3::new(0.95, 0.0, 0.0, 0.0, 1.1, 0.0, 0.0, 0.0, 0.98);
    let hard2 = Vector3::new(-0.05, 0.1, -0.2);
    let irons2 = Irons {
        soft: soft2,
        hard: hard2,
    };

    // Compose them
    let composed = irons2.clone().irons_on_irons(&irons1);

    // Apply transformations sequentially
    let once = irons1.calibrate_points(&[point]);
    let twice = irons2.calibrate_points(&once);

    // Apply composed transformation
    let direct = composed.calibrate_points(&[point]);

    let v_once = Vector3::new(once[0][0], once[0][1], once[0][2]);
    let diff = Vector3::new(
        twice[0][0] - direct[0][0],
        twice[0][1] - direct[0][1],
        twice[0][2] - direct[0][2],
    );
    assert!(
        diff.norm() / v_once.norm() < 0.01,
        "difference is too large: {}%",
        diff.norm() / v_once.norm() * 100.0
    );
}
