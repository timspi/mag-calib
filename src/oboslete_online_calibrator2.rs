use crate::{EllipsoidSpecificFitter, Irons, SphericalDownSampler, StatisticsAccumulator, calibration_mse};
use crate::online_calibrator::{OnlineMagCalibratorParameters, CalibrationReportedProgress};
use log;

//THIS FAILED MISERABLY
//The idea was that we can downsample before fits, and use irons_on_irons to improve an existing fit
//but I am not sure that irons_on_irons works as intended even though it was tested
//also, after downsampling, may times the points are not enough to fit

#[derive(Debug)]
enum CalibrationState {
    CantFit,
    ProduceFirstFitForSphereCovering,
    ProduceSubsequentFitForSphereCovering {
        irons: Irons,
    },
    SphereCovering {
        irons: Irons,
        points_to_refit: usize,
    },
    ProduceFinalFit {
        irons : Irons,
    },
}

pub struct OnlineMagCalibrator2 {
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

impl OnlineMagCalibrator2 {
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
            self.points_to_stall -= new_rawpoints_slice.len();
            if self.points_to_stall <= 0 {
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
                self.calibration_state = CalibrationState::ProduceFirstFitForSphereCovering;
                //we dont return here, we continue to next stage
            } else {
                return CalibrationReportedProgress::InProgress;
            }
        }

        if matches!(
            self.calibration_state,
            CalibrationState::ProduceFirstFitForSphereCovering
        ) {
            log::debug!(
                "producing first fit for sphere covering with {} points",
                self.all_collected_rawpoints.len(),
            );
            assert!(self.downsampler.fill_percentage(0) == 0.0); //must be empty at this point!

            let normalized_points = self
                .stats
                .normalize_by_mean_and_ranges(&self.all_collected_rawpoints);
            self.fitter.add_points(&normalized_points);

            let irons = match self.fitter.fit().and_then(|x| x.to_irons()) {
                Ok(ret) => ret.revert_normalization(self.stats.mean(),self.stats.axis_ranges()),
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

            self.calibrated_points.clear();
            irons.calibrate_points_into(&self.all_collected_rawpoints, &mut self.calibrated_points);
            self.downsampler
                .add_points(&self.calibrated_points[..self.all_collected_rawpoints.len()], [0.0; 3].into(), None).unwrap();
            self.last_cover_percentage = self.downsampler.fill_percentage(0);

            self.calibration_state = CalibrationState::SphereCovering {
                irons: irons,
                points_to_refit: self.params.points_for_refit_on_sphere_covering,
            };
            //we return here instead of going directly to sphere covering because all rawpoints were used
            return CalibrationReportedProgress::InProgress;
        }

        //we add points to downsampler after using current irons and check if we have enough points to fit
        if let CalibrationState::SphereCovering {
            ref irons,
            ref mut points_to_refit,
        } = self.calibration_state
        {
            let new_calibrated_points = irons.calibrate_points(new_rawpoints_slice);
            //downsampler is provided with calibrated points!
            self.downsampler
                .add_points(&new_calibrated_points, [0.0; 3].into(), None).unwrap();

            self.last_cover_percentage = self.downsampler.fill_percentage(0);

            if self.last_cover_percentage > self.params.acceptable_cover_percentage {
                self.calibration_state = CalibrationState::ProduceFinalFit {
                    irons: irons.clone(),
                };
                //we dont return, we continue to next stage
            } else if *points_to_refit == 0 {
                self.calibration_state = CalibrationState::ProduceSubsequentFitForSphereCovering {
                    irons: irons.clone(),
                };
                //we dont return, we continue to next stage
            } else {
                *points_to_refit = points_to_refit.saturating_sub(new_rawpoints_slice.len());
                return CalibrationReportedProgress::InProgress;
            }
        }

        if let CalibrationState::ProduceSubsequentFitForSphereCovering { ref irons } =
            self.calibration_state
        {
            log::debug!(
                "producing subsequent fit for sphere covering with {} points as last sphere cover was {}",
                self.all_collected_rawpoints.len(),
                self.last_cover_percentage,
            );
            
            //fit is on downsampled points, which are created from calibrated points, so we need to use irons_on_irons
            let target_points = self.downsampler.sample_buckets(None);
            let new_irons = match self.fitter.reset().add_points(&target_points).fit().and_then(|x| x.to_irons()) {
                Ok(ret) => ret.irons_on_irons(irons),
                Err(e) => {
                    log::warn!(
                        "fitting failed when trying to produce subsequent fit for sphere covering with {} points: {}",
                        self.all_collected_rawpoints.len(),
                        e,
                    );
                    self.points_to_stall = self.params.stall_points_on_fit_fail;
                    self.calibration_state = CalibrationState::SphereCovering {
                        irons: irons.clone(),
                        points_to_refit: self.params.points_for_refit_on_sphere_covering,
                    };
                    return CalibrationReportedProgress::InProgress;
                }
            };

            self.calibrated_points.clear();
            new_irons.calibrate_points_into(&self.all_collected_rawpoints, &mut self.calibrated_points);
            //downsampler is provided with calibrated points!
            self.downsampler.reset();
            self.downsampler
                .add_points(&self.calibrated_points[..self.all_collected_rawpoints.len()],
                            [0.0; 3].into(),
                            None).unwrap();
            println!("MSE: {}", calibration_mse(&self.calibrated_points, 1.0));
            
            self.calibration_state = CalibrationState::SphereCovering {
                irons: new_irons,
                points_to_refit: self.params.points_for_refit_on_sphere_covering,
            };
            return CalibrationReportedProgress::InProgress;
        }

        if let CalibrationState::ProduceFinalFit { ref irons } =
            self.calibration_state
        {
            log::debug!(
                "producing candidate fit with {} points as last sphere cover was {}",
                self.all_collected_rawpoints.len(),
                self.last_cover_percentage,
            );
           
            let target_points = self.downsampler.sample_buckets(None);
            let new_irons = match self.fitter.reset().add_points(&target_points).fit().and_then(|x| x.to_irons()) {
                Ok(ret) => ret.irons_on_irons(irons),
                Err(e) => {
                    log::warn!(
                        "fitting failed when trying to produce a final fit for sphere covering with {} points: {}",
                        self.all_collected_rawpoints.len(),
                        e,
                    );
                    self.points_to_stall = self.params.stall_points_on_fit_fail;
                    self.calibration_state = CalibrationState::SphereCovering {
                        irons: irons.clone(),
                        points_to_refit: self.params.points_for_refit_on_sphere_covering,
                    };
                    return CalibrationReportedProgress::InProgress;
                }
            };

            self.calibrated_points.clear();
            new_irons.calibrate_points_into(&self.all_collected_rawpoints, &mut self.calibrated_points);
            //downsampler is provided with calibrated points!
            self.downsampler.reset();
            self.downsampler
                .add_points(&self.calibrated_points[..self.all_collected_rawpoints.len()],
                            [0.0; 3].into(),
                            None).unwrap();

            let current_cover_precentage = self.downsampler.fill_percentage(0);
            if current_cover_precentage > self.params.acceptable_cover_percentage {
                log::info!("
                    mag calibrator found succesful fit with {} points and last sphere cover was {}",
                    self.all_collected_rawpoints.len(),
                    self.last_cover_percentage,
                );
                self.calibration_state = CalibrationState::ProduceFinalFit { irons: new_irons };
                return CalibrationReportedProgress::Finished;
            } else {
                log::warn!("
                    mag calibrator found candidate fit, but it downgraded sphere cover from {} to {}, continuing with sphere covering",
                    self.last_cover_percentage,
                    current_cover_precentage,
                );
                self.last_cover_percentage = current_cover_precentage;
                self.points_to_stall = self.params.stall_points_on_fit_fail;
                self.calibration_state = CalibrationState::SphereCovering {
                    irons: new_irons,
                    points_to_refit: self.params.points_for_refit_on_sphere_covering,
                };
                return CalibrationReportedProgress::InProgress;
            }
        }

        panic!("Should never get here!");
    }

    //we scale to irons to the earth magnetic field strength only on output
    pub fn get_irons(&self) -> Result<Irons, &'static str> {
        match self.calibration_state {
            CalibrationState::ProduceFinalFit {ref irons } => Ok(irons.clone()
                .scale_soft(self.params.earth_magnetic_field_strength)),
            _ => Err("Calibration not successful"),
        }
    }

    pub fn get_calibrated_points(&self) -> Result<&Vec<[f64; 3]>, &'static str> {
        panic!("this is wrong.. should scale calibrated points? they are using irons that were not scaled to earth magnetic field strength");
        match self.calibration_state {
            CalibrationState::ProduceFinalFit {..} => Ok(&self.calibrated_points),
            _ => Err("Calibration not successful"),
        }
    }

    pub fn get_cover_percentage(&self) -> f64 {
        self.last_cover_percentage
    }
}
