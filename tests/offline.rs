use csv::Reader;
use mag_calib::{
    calibration_mse, EllipsoidSpecificFitter, SphericalDownSampler, StatisticsAccumulator,
};
use std::collections::HashMap;
use std::fs;

mod plots;

const EARTH_MAGNETIC_FIELD_STRENGTH: f64 = 500.0;

#[test]
fn offline_calib() {
    let csv_filepath = "raw_data/imu_gnss_data.csv";
    let mut reader = Reader::from_path(csv_filepath).expect("Failed to open CSV file");
    let headers = reader.headers().expect("Failed to read headers");

    let mut header_map: HashMap<&str, usize> = HashMap::new();
    for (i, header) in headers.iter().enumerate() {
        header_map.insert(header.trim(), i);
    }
    let idx_x = *header_map.get("magx").expect("Missing 'magx' column");
    let idx_y = *header_map.get("magy").expect("Missing 'magy' column");
    let idx_z = *header_map.get("magz").expect("Missing 'magz' column");

    let mut raw_points = Vec::new();
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        let x: f64 = record[idx_x].parse().expect("Failed to parse magx");
        let y: f64 = record[idx_y].parse().expect("Failed to parse magy");
        let z: f64 = record[idx_z].parse().expect("Failed to parse magz");
        raw_points.push([x, y, z]);
    }

    let mut raw_pnt_fitter = EllipsoidSpecificFitter::new();
    let ret = raw_pnt_fitter.add_points(&raw_points).fit();
    assert!(ret.is_ok(), "Ellipsoid fitting failed");
    let raw_method_coeffs = ret.unwrap();
    let raw_method_irons = raw_method_coeffs
        .to_irons()
        .expect("Failed to convert to Irons")
        .scale_soft(EARTH_MAGNETIC_FIELD_STRENGTH);

    //-------------------compare to normalized points-------------------
    let mut stats = StatisticsAccumulator::new();
    stats.add_points(&raw_points);
    let normalized_points = stats.normalize_by_mean_and_ranges(&raw_points);
    let mut normalized_pnt_fitter = EllipsoidSpecificFitter::new();
    let ret = normalized_pnt_fitter.add_points(&normalized_points).fit();
    assert!(ret.is_ok(), "Ellipsoid fitting failed");
    let normalized_method_coeffs = ret.unwrap();
    let normalized_method_irons = normalized_method_coeffs
        .to_irons()
        .expect("Failed to convert to Irons")
        .revert_normalization(stats.mean(), stats.axis_ranges())
        .scale_soft(EARTH_MAGNETIC_FIELD_STRENGTH);

    assert!(
        (raw_method_irons.soft - normalized_method_irons.soft).norm()
            / raw_method_irons.soft.norm()
            < 0.05
    );
    assert!(
        (raw_method_irons.hard - normalized_method_irons.hard).norm()
            / raw_method_irons.hard.norm()
            < 0.05
    );

    let raw_method_calibrated_points = raw_method_irons.calibrate_points(&raw_points);
    let raw_method_mse =
        calibration_mse(&raw_method_calibrated_points, EARTH_MAGNETIC_FIELD_STRENGTH);
    let normalized_method_calibrated_points = normalized_method_irons.calibrate_points(&raw_points);
    let normalized_method_mse = calibration_mse(
        &normalized_method_calibrated_points,
        EARTH_MAGNETIC_FIELD_STRENGTH,
    );

    assert!((raw_method_mse - normalized_method_mse).abs() / raw_method_mse < 0.3);
    println!(
        "raw_method_mse: {}, normalized_method_mse: {}",
        raw_method_mse, normalized_method_mse
    );

    //create a directory for test outputs
    fs::create_dir_all("test_outputs/offline").expect("Failed to create output directory");

    let plot_html = plots::scatter3d(
        vec![
            &raw_method_calibrated_points,
            &normalized_method_calibrated_points,
        ],
        vec!["raw_method".to_string(), "normalized_method".to_string()],
    );
    fs::write("test_outputs/offline/calibrated_points.html", plot_html)
        .expect("Failed to write calibrated_points HTML");

    let mut downsampler = SphericalDownSampler::new(20, 10);
    downsampler
        .add_points(&raw_method_calibrated_points, [0.0; 3].into(), None)
        .unwrap();
    println!(
        "fill percentage of calibrated points: {}",
        downsampler.fill_percentage(0)
    );

    let heatmap_html = plots::heatmap_buckets(&downsampler, None);
    fs::write("test_outputs/offline/heatmap.html", heatmap_html)
        .expect("Failed to write heatmap HTML");

    let scatter3d_buckets_html =
        plots::scatter3d_buckets_with_graticules(&downsampler, EARTH_MAGNETIC_FIELD_STRENGTH);
    fs::write(
        "test_outputs/offline/scatter3d_buckets_with_graticules.html",
        scatter3d_buckets_html,
    )
    .expect("Failed to write scatter3d buckets HTML");
}
