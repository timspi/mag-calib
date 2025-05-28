use crate::{
    calibration_mse, EllipsoidSpecificFitter, Irons, SphericalDownSampler, StatisticsAccumulator,
};
use log;

#[derive(Debug)]
enum CalibrationState {
    CantFit,
    ProduceFitForSphereCovering,
    SphereCovering {
        irons: Irons,
        points_to_refit: usize,
    },
    ProducingCandidateFit,
    ProducedSuccessfulFit {
        irons: Irons,
    },
}

pub enum CalibrationReportedProgress {
    NotStarted,
    InProgress,
    Finished,
}

pub struct OnlineMagCalibratorParameters {
    pub min_points_for_first_fit: usize,
    pub min_covariance_for_first_fit: f64,
    pub min_per_axis_range_for_first_fit: f64,
    pub earth_magnetic_field_strength: f64,
    pub acceptable_cover_percentage: f64,
    pub downsampler_lon_div: usize,
    pub downsampler_height_div: usize,
    pub stall_points_on_fit_fail: usize,
    pub points_for_refit_on_sphere_covering: usize,
}

impl OnlineMagCalibratorParameters {
    pub fn default_from_earth_magnetic_field_strength(earth_magnetic_field_strength: f64) -> Self {
        Self {
            min_points_for_first_fit: 100,
            min_covariance_for_first_fit: 100.0,
            min_per_axis_range_for_first_fit: earth_magnetic_field_strength / 2.0,
            earth_magnetic_field_strength: earth_magnetic_field_strength,
            acceptable_cover_percentage: 0.6,
            downsampler_lon_div: 10,
            downsampler_height_div: 10,
            stall_points_on_fit_fail: 50,
            points_for_refit_on_sphere_covering: 100,
        }
    }
}

pub struct OnlineMagCalibrator {
    calibration_state: CalibrationState,
    stats: StatisticsAccumulator,
    all_collected_rawpoints: Vec<[f64; 3]>,
    calibrated_points: Vec<[f64; 3]>, //store here to allow easy access for user if he wants 'em
    downsampler: SphericalDownSampler,
    fitter: EllipsoidSpecificFitter,
    params: OnlineMagCalibratorParameters,
    last_cover_percentage: f64,
    points_to_stall: usize,
}

impl OnlineMagCalibrator {
    pub fn new(params: OnlineMagCalibratorParameters) -> Self {
        assert!(
            params.acceptable_cover_percentage > 0.0 && params.acceptable_cover_percentage < 1.0,
            "acceptable_cover_percentage must be between 0.0 and 1.0"
        );
        Self {
            calibration_state: CalibrationState::CantFit,
            stats: StatisticsAccumulator::new(),
            all_collected_rawpoints: Vec::with_capacity(3000), //generally a big enough buffer
            calibrated_points: Vec::with_capacity(3000),
            downsampler: SphericalDownSampler::new(
                params.downsampler_lon_div,
                params.downsampler_height_div,
            ),
            fitter: EllipsoidSpecificFitter::new(),
            params: params.into(),
            last_cover_percentage: 0.0,
            points_to_stall: 0,
        }
    }

    pub fn add_points(&mut self, new_rawpoints_slice: &[[f64; 3]]) -> CalibrationReportedProgress {
        self.stats.add_points(new_rawpoints_slice);
        self.all_collected_rawpoints
            .extend_from_slice(new_rawpoints_slice);

        if self.points_to_stall > 0 {
            self.points_to_stall = self.points_to_stall.saturating_sub(new_rawpoints_slice.len());
            if self.points_to_stall != 0 {
                return CalibrationReportedProgress::InProgress;
            }
        }

        if matches!(self.calibration_state, CalibrationState::CantFit) {
            if self.all_collected_rawpoints.len() > self.params.min_points_for_first_fit
                && self
                    .stats
                    .axis_ranges()
                    .iter()
                    .all(|x| *x > self.params.min_per_axis_range_for_first_fit)
                && self.stats.cov_condition_number() < self.params.min_covariance_for_first_fit
            {
                log::info!(
                    "mag calibrator can produce first fit for sphere covering after {}",
                    self.all_collected_rawpoints.len(),
                );
                self.calibration_state = CalibrationState::ProduceFitForSphereCovering;
            } else {
                return CalibrationReportedProgress::InProgress;
            }
        }

        if matches!(
            self.calibration_state,
            CalibrationState::ProduceFitForSphereCovering
        ) {
            log::debug!(
                "producing fit for sphere covering with {} points as last sphere cover was {} / {}",
                self.all_collected_rawpoints.len(),
                self.last_cover_percentage,
                self.params.acceptable_cover_percentage,
            );
            let normalized_points = self
                .stats
                .normalize_by_mean_and_ranges(&self.all_collected_rawpoints);

            let irons = match self
                .fitter
                .reset()
                .add_points(&normalized_points)
                .fit()
                .and_then(|x| x.to_irons())
                .map(|x| x.revert_normalization(self.stats.mean(), self.stats.axis_ranges()))
            {
                Ok(ret) => ret,
                Err(e) => {
                    log::warn!(
                        "fitting failed when trying to produce fit for sphere covering with {} points: {}",
                        self.all_collected_rawpoints.len(),
                        e,
                    );
                    self.points_to_stall = self.params.stall_points_on_fit_fail;
                    return CalibrationReportedProgress::InProgress;
                }
            };

            self.downsampler.reset();
            self.calibrated_points.clear();
            irons.calibrate_points_into(&self.all_collected_rawpoints, &mut self.calibrated_points);
            self.downsampler
                .add_points(
                    &self.calibrated_points[..self.all_collected_rawpoints.len()],
                    [0.0; 3].into(),
                    None,
                )
                .expect("cant fail here, s_in_bin is None");
            self.last_cover_percentage = self.downsampler.fill_fraction();

            self.calibration_state = CalibrationState::SphereCovering {
                irons: irons,
                points_to_refit: self.params.points_for_refit_on_sphere_covering,
            };
            //we return here instead of going directly to sphere covering because all rawpoints were used
            return CalibrationReportedProgress::InProgress;
        }

        if let CalibrationState::SphereCovering {
            ref irons,
            ref mut points_to_refit,
        } = self.calibration_state
        {
            let new_calibrated_points = irons.calibrate_points(new_rawpoints_slice);
            self.downsampler
                .add_points(&new_calibrated_points, [0.0; 3].into(), None)
                .expect("cant fail here, s_in_bin is None");

            self.last_cover_percentage = self.downsampler.fill_fraction();

            if self.last_cover_percentage > self.params.acceptable_cover_percentage {
                self.calibration_state = CalibrationState::ProducingCandidateFit;
                //we dont return here because we dont need to update downsampler with new points, continue to produce candidate fit
            } else if *points_to_refit == 0 {
                self.calibration_state = CalibrationState::ProduceFitForSphereCovering;
                return CalibrationReportedProgress::InProgress;
            } else {
                *points_to_refit = points_to_refit.saturating_sub(new_rawpoints_slice.len());
                return CalibrationReportedProgress::InProgress;
            }
        }

        if matches!(
            self.calibration_state,
            CalibrationState::ProducingCandidateFit
                | CalibrationState::ProducedSuccessfulFit { .. }
        ) {
            log::info!(
                "producing candidate rough irons fit with {} points as last sphere cover was {}",
                self.all_collected_rawpoints.len(),
                self.last_cover_percentage,
            );
            //user may want to add points none the less
            let normalized_points = self
                .stats
                .normalize_by_mean_and_ranges(&self.all_collected_rawpoints);

            let rough_irons = match self
                .fitter
                .reset()
                .add_points(&normalized_points)
                .fit()
                .and_then(|x| x.to_irons())
                .map(|x| {
                    x.revert_normalization(self.stats.mean(), self.stats.axis_ranges())
                        .scale_soft(self.params.earth_magnetic_field_strength)
                }) {
                Ok(ret) => ret,
                Err(e) => {
                    //we dont have irons to go back to sphere covering... so lets just wait for more points
                    log::warn!(
                        "fitting failed when trying to produce candidate fit: {}, with {} points,
                         waiting for {} more points and trying again",
                        e,
                        self.all_collected_rawpoints.len(),
                        self.points_to_stall,
                    );
                    self.points_to_stall = self.params.stall_points_on_fit_fail;
                    return CalibrationReportedProgress::InProgress;
                }
            };

            //add all used points to downsampler
            self.downsampler.reset();
            self.calibrated_points.clear();
            rough_irons
                .calibrate_points_into(&self.all_collected_rawpoints, &mut self.calibrated_points);
            self.downsampler
                .add_points(
                    &self.calibrated_points[..self.all_collected_rawpoints.len()],
                    [0.0; 3].into(),
                    Some(&self.all_collected_rawpoints),
                )
                .expect("s_in_bin and s_out_bin should be of the same size");
            self.last_cover_percentage = self.downsampler.fill_fraction();
            if self.last_cover_percentage < self.params.acceptable_cover_percentage {
                log::warn!("mag calibrator tried to produce a rough iron fit, but the sphere cover was reduced from {} to {},
                 as such, more points are required, heading back to sphere covering",
                    self.last_cover_percentage,
                    self.downsampler.fill_fraction(),
                );
                self.calibration_state = CalibrationState::SphereCovering {
                    irons: rough_irons,
                    points_to_refit: self.params.points_for_refit_on_sphere_covering,
                };
                return CalibrationReportedProgress::InProgress;
            }

            let target_points = self.downsampler.sample_buckets(None);
            let downsampled_irons = match self
                .fitter
                .reset()
                .add_points(&target_points)
                .fit()
                .and_then(|x| x.to_irons())
                .map(|x| {x.scale_soft(self.params.earth_magnetic_field_strength)
                }) {
                Ok(ret) => ret,
                Err(e) => {
                    log::warn!(
                        "fitting failed when trying to produce downsampled irons when all collected points were {}: {}
                        heading back to sphere covering",
                        self.all_collected_rawpoints.len(),
                        e,
                    );
                    self.calibration_state = CalibrationState::SphereCovering {
                        irons: rough_irons,
                        points_to_refit: self.params.points_for_refit_on_sphere_covering,
                    };
                    return CalibrationReportedProgress::InProgress;
                }
            };

            log::info!(
                "produced downsampled irons with {} points",
                target_points.len(),
            );

            // calibrate all collected points with the downsampled irons
            self.calibrated_points.clear();
            downsampled_irons
                .calibrate_points_into(&self.all_collected_rawpoints, &mut self.calibrated_points);
            // --------------------------------------------------------------------------------------------------

            self.calibration_state = CalibrationState::ProducedSuccessfulFit {
                irons: downsampled_irons,
            };
            return CalibrationReportedProgress::Finished;
        }

        panic!("Should never get here!");
    }

    pub fn get_irons(&self) -> Result<Irons, &'static str> {
        match self.calibration_state {
            CalibrationState::ProducedSuccessfulFit { ref irons } => Ok(irons.clone()),
            _ => Err("Calibration not successful"),
        }
    }

    pub fn get_calibrated_points(&self) -> Result<&[[f64; 3]], &'static str> {
        match self.calibration_state {
            CalibrationState::ProducedSuccessfulFit { .. } => {
                Ok(&self.calibrated_points[..self.all_collected_rawpoints.len()])
            }
            _ => Err("Calibration not successful"),
        }
    }

    pub fn calculate_mse(&self) -> Result<f64, &'static str> {
        match self.calibration_state {
            CalibrationState::ProducedSuccessfulFit { .. } => Ok(calibration_mse(
                &self.calibrated_points[..self.all_collected_rawpoints.len()],
                self.params.earth_magnetic_field_strength,
            )),
            _ => Err("Calibration not successful"),
        }
    }

    pub fn get_cover_percentage(&self) -> f64 {
        self.last_cover_percentage
    }
}
