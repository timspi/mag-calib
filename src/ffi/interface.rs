//this file generates the interface
use super::super::*;

#[cfg(feature = "java_wrapper")]
use crate::ffi::jni_c_header::*;

//logger
foreign_class!(
    class Logger {
        fn init_logger_with_level(level: &str) {
            let level = match level.to_lowercase().as_str() {
                "error" => log::Level::Error,
                "warn" => log::Level::Warn,
                "info" => log::Level::Info,
                "debug" => log::Level::Debug,
                "trace" => log::Level::Trace,
                _ => {
                    panic!("
                    Logger level not recognized, please use one of the following:
                    Error
                    Warn
                    Info
                    Debug
                    Trace
                    ");
                }
            };
            simple_logger::init_with_level(level).unwrap();
        }
    }
);
use online_calibrator::CalibrationReportedProgress;
foreign_enum!(
    enum CalibrationReportedProgress {
        NotStarted = CalibrationReportedProgress::NotStarted,
        InProgress = CalibrationReportedProgress::InProgress,
        Finished = CalibrationReportedProgress::Finished,
    }
);

use online_calibrator::OnlineMagCalibratorParameters;
foreign_class!(
    class OnlineMagCalibratorParameters {
        self_type OnlineMagCalibratorParameters;
        constructor OnlineMagCalibratorParameters::new(
            min_points_for_first_fit: usize,
            min_covariance_for_first_fit: f64,
            min_per_axis_range_for_first_fit: f64,
            earth_magnetic_field_strength: f64,
            acceptable_cover_percentage: f64,
            downsampler_lon_div: usize,
            downsampler_height_div: usize,
            stall_points_on_fit_fail: usize,
            points_for_refit_on_sphere_covering: usize,
        ) -> OnlineMagCalibratorParameters {
            OnlineMagCalibratorParameters{
                min_points_for_first_fit,
                min_covariance_for_first_fit,
                min_per_axis_range_for_first_fit,
                earth_magnetic_field_strength,
                acceptable_cover_percentage,
                downsampler_lon_div,
                downsampler_height_div,
                stall_points_on_fit_fail,
                points_for_refit_on_sphere_covering,
        }
    }

    fn OnlineMagCalibratorParameters::default_from_earth_magnetic_field_strength(
        earth_magnetic_field_strength: f64,
    ) -> OnlineMagCalibratorParameters;
});

pub struct IronsFFI {
    soft_xx: f64,
    soft_yy: f64,
    soft_zz: f64,
    soft_xy: f64,
    soft_xz: f64,
    soft_yz: f64,
    hard_x: f64,
    hard_y: f64,
    hard_z: f64,
}

impl From<Irons> for IronsFFI {
    fn from(iron: Irons) -> Self {
        Self {
            soft_xx: iron.soft[(0, 0)],
            soft_yy: iron.soft[(1, 1)],
            soft_zz: iron.soft[(2, 2)],
            soft_xy: iron.soft[(0, 1)],
            soft_xz: iron.soft[(0, 2)],
            soft_yz: iron.soft[(1, 2)],
            hard_x: iron.hard[0],
            hard_y: iron.hard[1],
            hard_z: iron.hard[2],
        }
    }
}

impl IronsFFI {
    pub fn get_soft_xx(&self) -> f64 {
        self.soft_xx
    }
    pub fn get_soft_yy(&self) -> f64 {
        self.soft_yy
    }
    pub fn get_soft_zz(&self) -> f64 {
        self.soft_zz
    }
    pub fn get_soft_xy(&self) -> f64 {
        self.soft_xy
    }
    pub fn get_soft_xz(&self) -> f64 {
        self.soft_xz
    }
    pub fn get_soft_yz(&self) -> f64 {
        self.soft_yz
    }
    pub fn get_hard_x(&self) -> f64 {
        self.hard_x
    }
    pub fn get_hard_y(&self) -> f64 {
        self.hard_y
    }
    pub fn get_hard_z(&self) -> f64 {
        self.hard_z
    }
}

foreign_class!(
    class IronsFFI {
        self_type IronsFFI;
        constructor IronsFFI::new(
            soft_xx: f64,
            soft_yy: f64,
            soft_zz: f64,
            soft_xy: f64,
            soft_xz: f64,
            soft_yz: f64,
            hard_x: f64,
            hard_y: f64,
            hard_z: f64,
        ) -> IronsFFI {
            IronsFFI {
                soft_xx,
                soft_yy,
                soft_zz,
                soft_xy,
                soft_xz,
                soft_yz,
                hard_x,
                hard_y,
                hard_z
            }
        }
        fn IronsFFI::get_soft_xx(&self) -> f64;
        fn IronsFFI::get_soft_yy(&self) -> f64;
        fn IronsFFI::get_soft_zz(&self) -> f64;
        fn IronsFFI::get_soft_xy(&self) -> f64;
        fn IronsFFI::get_soft_xz(&self) -> f64;
        fn IronsFFI::get_soft_yz(&self) -> f64;
        fn IronsFFI::get_hard_x(&self) -> f64;
        fn IronsFFI::get_hard_y(&self) -> f64;
        fn IronsFFI::get_hard_z(&self) -> f64;
    }
);

pub struct TenPointContainer {
    points: [[f64; 3];10],
    counter: usize,
}

impl TenPointContainer {
    pub fn new() -> Self {
        Self {
            points: [[0.0; 3]; 10],
            counter: 0,
        }
    }

    pub fn add_point(&mut self, x: f64, y: f64, z: f64) -> usize  {
        if self.counter < 10 {
            self.points[self.counter] = [x, y, z];
            self.counter += 1;
            return self.counter
        } else {
            10
        }
    }

    pub fn clear(&mut self) {
        self.counter = 0;
    }
}

foreign_class!(
    class TenPointContainer {
    self_type TenPointContainer;
    constructor TenPointContainer::new() -> TenPointContainer;
    fn TenPointContainer::add_point(&mut self, x: f64, y: f64, z: f64) -> usize;
    fn TenPointContainer::clear(&mut self);
});

use online_calibrator::OnlineMagCalibrator;

impl OnlineMagCalibrator {
    pub fn take_points_from_container(
        &mut self,
        ten_point_container: &mut TenPointContainer,
    ) -> CalibrationReportedProgress {
        let calib_progress = self.add_points(&ten_point_container.points[..ten_point_container.counter]);
        ten_point_container.clear();
        calib_progress
    }
    pub fn get_irons_ffi(&self) -> Result<IronsFFI, &'static str> {
        let irons = self.get_irons()?;
        Ok(IronsFFI::from(irons))
    }
}

foreign_class!(
    class OnlineMagCalibrator {
        self_type OnlineMagCalibrator;
        constructor OnlineMagCalibrator::new(
            builder_params: OnlineMagCalibratorParameters,
        ) -> OnlineMagCalibrator {OnlineMagCalibrator::new(builder_params)}

        fn OnlineMagCalibrator::take_points_from_container(&mut self, ten_point_container : &mut TenPointContainer) -> CalibrationReportedProgress;
        fn OnlineMagCalibrator::get_irons_ffi(&self) -> Result<IronsFFI, &'static str>; alias get_irons;
        fn OnlineMagCalibrator::get_cover_percentage(&self) -> f64;
        fn OnlineMagCalibrator::calculate_mse(&self) -> Result<f64, &'static str>;
    }
);



