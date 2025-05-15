use csv::Reader;
use mag_calib::online_calibrator::{
    CalibrationReportedProgress, OnlineMagCalibrator, OnlineMagCalibratorParameters,
};

mod plots;

use simple_logger::SimpleLogger;
use std::collections::HashMap;
use std::fs;

const BUFFER_SIZE: usize = 10;
const EARTH_MAGNETIC_FIELD_STRENGTH: f64 = 500.0;

#[test]
fn online_calib() {
    SimpleLogger::new()
        .with_level(log::LevelFilter::Debug)
        .init()
        .expect("Failed to initialize logger");

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

    let mut all_rawpoints = Vec::new();
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        let x: f64 = record[idx_x].parse().expect("Failed to parse magx");
        let y: f64 = record[idx_y].parse().expect("Failed to parse magy");
        let z: f64 = record[idx_z].parse().expect("Failed to parse magz");
        all_rawpoints.push([x, y, z]);
    }

    //ready outputs:
    fs::create_dir_all("test_outputs/online").expect("Failed to create output directory");

    let params = OnlineMagCalibratorParameters::default_from_earth_magnetic_field_strength(
        EARTH_MAGNETIC_FIELD_STRENGTH,
    );
    let mut online_calibrator = OnlineMagCalibrator::new(params);
    let mut calibration_progress = CalibrationReportedProgress::NotStarted;

    let mut inserted_points = 0;
    while !matches!(calibration_progress, CalibrationReportedProgress::Finished)
        && inserted_points + BUFFER_SIZE < all_rawpoints.len()
    {
        let new_rawpoints_to_insert =
            &all_rawpoints[inserted_points..inserted_points + BUFFER_SIZE];
        inserted_points += BUFFER_SIZE;
        calibration_progress = online_calibrator.add_points(new_rawpoints_to_insert);
    }

    match calibration_progress {
        CalibrationReportedProgress::Finished => {
            println!("Calibration finished successfully");
        }
        CalibrationReportedProgress::InProgress => {
            println!("Did not finish calibration");
        }
        _ => {}
    }
    println!(
        "inserted points: {}, all_points: {}",
        inserted_points,
        all_rawpoints.len()
    );
    println!(
        "cover percentage: {}",
        online_calibrator.get_cover_percentage()
    );
    println!(
        "calibration mse: {}",
        online_calibrator.calculate_mse().unwrap()
    );

    if let Ok(calibrated_points_slice) = online_calibrator.get_calibrated_points() {
        let plot_html = plots::scatter3d(
            vec![calibrated_points_slice],
            vec!["calibrated_points".to_string()],
        );
        fs::write("test_outputs/online/calibrated_points.html", plot_html)
            .expect("Failed to write calibrated_points HTML");
    }
}
