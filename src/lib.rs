use nalgebra::SMatrix;
use nalgebra::{DMatrix, Matrix3, Matrix6, SymmetricEigen, Vector3};
type Matrix10<T> = SMatrix<T, 10, 10>;
use std::collections::BTreeMap;
use std::fmt::Display;

mod ffi;
pub mod online_calibrator;
use std::f64::consts::PI;

#[derive(Debug)]
pub enum MagCalibErrors {
    MatrixSolveError,
    MatrixSolveErrorOnC1,
    ShouldNeverHappen,
    NotEnoughPoints,
    FailedFitting,
    NegativeEigenvalues,
    MatrixInverseError,
    WrongEtaParameterInput,
    FittedNonEllipsoidQuadratic,
    SizeMismatch,
    CoudlntGetMaxEigenValue,
}

impl Display for MagCalibErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MagCalibErrors::MatrixSolveError => write!(f, "Matrix solve error"),
            MagCalibErrors::MatrixSolveErrorOnC1 => write!(f, "Matrix solve error on C1"),
            MagCalibErrors::ShouldNeverHappen => write!(f, "This should never happen"),
            MagCalibErrors::NotEnoughPoints => write!(f, "Not enough points for fitting"),
            MagCalibErrors::FailedFitting => write!(f, "Failed to fit the ellipsoid"),
            MagCalibErrors::NegativeEigenvalues => write!(f, "Negative eigenvalues found"),
            MagCalibErrors::MatrixInverseError => write!(f, "Matrix inverse error"),
            MagCalibErrors::WrongEtaParameterInput => {
                write!(f, "Wrong eta parameter input, can be between 1.0 and INF. 
                                                                    1.0 is recommended as it enfroces ellipsoid specificity")
            }
            MagCalibErrors::FittedNonEllipsoidQuadratic => {
                write!(f, "Fitted a non-ellipsoid quadratic")
            }
            MagCalibErrors::SizeMismatch => write!(f, "Size mismatch"),
            MagCalibErrors::CoudlntGetMaxEigenValue => {
                write!(f, "Couldn't get max eigenvalue")
            }
        }
    }
}
#[derive(Debug, Clone)]
pub struct Irons {
    pub soft: Matrix3<f64>,
    pub hard: Vector3<f64>,
}
impl Irons {
    pub fn calibrate_points_into(&self, s: &[[f64; 3]], inout: &mut Vec<[f64; 3]>) {
        if inout.len() < s.len() {
            inout.resize(s.len()*2, [0.0, 0.0, 0.0]); //dont want to resize it every time
        }
        for i in 0..s.len() {
            let x = s[i][0];
            let y = s[i][1];
            let z = s[i][2];

            let col = Vector3::new(x, y, z);
            let calibrated = self.soft * (col - self.hard);
            inout[i] = [calibrated.x, calibrated.y, calibrated.z];
        }
    }

    pub fn calibrate_points(&self, s: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut out = Vec::with_capacity(s.len());
        for i in 0..s.len() {
            out.push(self.calibrate_point(&s[i]));
        }
        out
    }

    pub fn calibrate_point(&self, s: &[f64; 3]) -> [f64; 3] {
        let col = Vector3::from_column_slice(s);
        let calibrated = self.soft * (col - self.hard);
        [calibrated.x, calibrated.y, calibrated.z]
    }

    pub fn scale_soft(mut self, scale: f64) -> Self {
        self.soft.scale_mut(scale);
        self
    }

    //irons takes an ellipsoid points and transforms them to a sphere (hopefully)
    //we can take these sphere point and fit a new sphere to them, creating a new set of irons
    //we can then compose the two transformations so that we can go from the original ellipsoid to the new sphere

    //y = A * (x - b)
    //y' = A' * (y - b')
    //y' = A' * (A * (x - b) - b')
    //y' = A' * (A * x - A * b - A * A^-1 * b')
    //y' = A' A (x - (A * b + A^-1 * b'))
    pub fn irons_on_irons(mut self, old_irons: &Irons) -> Self {
        self.soft = self.soft * old_irons.soft;
        self.hard = old_irons.soft * old_irons.hard
            + self
                .soft
                .try_inverse()
                .expect("inversing soft should not fail")
                * self.hard;
        self
    }

        /*
    Reverts a linear model defined in normalized coordinates back to operate on original input.

    Original model (normalized space):
        y' = A * (y - b)
    where:
        y = (x - mean) / scale      // element-wise normalization

    Goal:
        Express y' directly in terms of unnormalized x:
            y' = A_new * (x - b_new)

    Derivation:
        y = (x - mean) / scale
        write S = diag(scale) and iS = diag(1/scale)
        y = iS * (x - mean)

        y' = A * (y - b)
           = A * (iS * (x - mean) - b)
           = A *(iS *x - iS *mean -  iS * S * b)
           = A * iS (x - mean - S * b)

    Therefore:
        A_new = A * iS
        b_new = mean + S * b

    we do this element wise just because we can
    */
    pub fn revert_normalization(mut self, mean: Vector3<f64>, scale: Vector3<f64>) -> Self {
        let recip_scale = scale.map(|x| 1.0 / x);
        for i in 0..3 {
            self.soft.column_mut(i).scale_mut(recip_scale[i]);
        }
        self.hard = mean + scale.component_mul(&self.hard);
        self
    }
}

#[derive(Debug)]
pub struct EllipsoidAlgebraicCoeffs {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub h: f64,
    pub g: f64,
    pub f: f64,
    pub p: f64,
    pub q: f64,
    pub r: f64,
    pub d: f64,
}

// ax^2 + by^2 + cy^2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0
impl EllipsoidAlgebraicCoeffs {
    /*
    Its still a bit of a mystery to me how we extract soft and hard. We followed the derivation from
    https://teslabs.com/articles/magnetometer-calibration/#mjx-eqn-eq_h_const
    There is another way to do it with 4x4 Q matrix and a 4x4 T that transforms Q into the origin.. this is more efficient
    */
    #[allow(non_snake_case)]
    pub fn to_irons(&self) -> Result<Irons, MagCalibErrors> {
        let m = Matrix3::new(
            self.a, self.h, self.g, self.h, self.b, self.f, self.g, self.f, self.c,
        );
        let n = Vector3::new(self.p, self.q, self.r);

        let inv_m = m.try_inverse().ok_or(MagCalibErrors::MatrixInverseError)?;
        //Note that m matrix follows Li's convention. often it will be parameterised differently and then 0.5 will apear here
        let hard = -inv_m * n;

        let eig = m.symmetric_eigen();
        if eig.eigenvalues.iter().any(|&x| x <= 0.0) {
            return Err(MagCalibErrors::NegativeEigenvalues);
        }

        let sqrt_d = Matrix3::from_diagonal(&eig.eigenvalues.map(|x| x.sqrt()));
        let v = eig.eigenvectors;
        let unscaled_transform = v * sqrt_d * v.transpose();

        let denom = ((n.transpose() * inv_m * n)[0] - self.d).sqrt();
        let soft = unscaled_transform / denom;

        Ok(Irons { soft, hard })
    }
}

pub fn calibration_mse(calibrated_points: &[[f64; 3]], mag_field_strength: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..calibrated_points.len() {
        let x = calibrated_points[i][0];
        let y = calibrated_points[i][1];
        let z = calibrated_points[i][2];

        sum += ((x * x + y * y + z * z).sqrt() - mag_field_strength).powi(2);
    }
    (sum / calibrated_points.len() as f64).sqrt()
}

/*
Based on 2004 paper by Qingde Li and Gareth Griffith named "Least squares ellipsoid specific fitting"
note: if points.len() > 100 it is recommended to normalize and then use revert_point_normalization of 'Irons' struct
*/
pub struct EllipsoidSpecificFitter {
    pub s: Matrix10<f64>,
    points_added: usize,
}

impl EllipsoidSpecificFitter {
    pub fn new() -> Self {
        Self {
            s: Matrix10::<f64>::zeros(),
            points_added: 0,
        }
    }

    pub fn reset(&mut self) -> &mut Self {
        self.s.fill(0.0);
        self.points_added = 0;
        self
    }

    pub fn add_points(&mut self, s: &[[f64; 3]]) -> &Self {
        let n = s.len();

        let mut d = DMatrix::<f64>::zeros(10, n);
        for j in 0..n {
            let x = s[j][0];
            let y = s[j][1];
            let z = s[j][2];

            d[(0, j)] = x * x;
            d[(1, j)] = y * y;
            d[(2, j)] = z * z;
            d[(3, j)] = 2.0 * y * z;
            d[(4, j)] = 2.0 * x * z;
            d[(5, j)] = 2.0 * x * y;
            d[(6, j)] = 2.0 * x;
            d[(7, j)] = 2.0 * y;
            d[(8, j)] = 2.0 * z;
            d[(9, j)] = 1.0;
        }

        //S  = d * d^T
        for i in 0..10 {
            for j in 0..10 {
                if i >= 6 && j < 4 {
                    continue;
                }; //no point filling s21 as s is symmetric
                for k in 0..n {
                    self.s[(i, j)] += d[(i, k)] * d[(j, k)]; //d(j,k) is the transpose thing
                }
            }
        }

        self.points_added += n;
        self
    }

    pub fn fit(&self) -> Result<EllipsoidAlgebraicCoeffs, MagCalibErrors> {
        if self.points_added <= 10 {
            return Err(MagCalibErrors::NotEnoughPoints);
        }

        let s11 = self.s.fixed_view::<6, 6>(0, 0);
        let s12 = self.s.fixed_view::<6, 4>(0, 6);
        let s22 = self.s.fixed_view::<4, 4>(6, 6);

        let k = 4.0;
        let mut c1 = Matrix6::<f64>::zeros();
        c1[(0, 0)] = -1.0;
        c1[(1, 1)] = -1.0;
        c1[(2, 2)] = -1.0;
        c1[(0, 1)] = k / 2.0 - 1.0;
        c1[(0, 2)] = k / 2.0 - 1.0;
        c1[(1, 0)] = k / 2.0 - 1.0;
        c1[(1, 2)] = k / 2.0 - 1.0;
        c1[(2, 0)] = k / 2.0 - 1.0;
        c1[(2, 1)] = k / 2.0 - 1.0;
        c1[(3, 3)] = -k;
        c1[(4, 4)] = -k;
        c1[(5, 5)] = -k;

        // Solve E = C1⁻¹ * (S11 - S12 * S22⁻¹ * S21)
        let s22_inv_times_s21 = s22
            .qr()
            .solve(&s12.transpose())
            .ok_or(MagCalibErrors::MatrixSolveError)?;
        let middle = s11 - s12 * s22_inv_times_s21;
        let e = c1
            .qr()
            .solve(&middle)
            .ok_or(MagCalibErrors::MatrixSolveErrorOnC1)?;

        // Solve eigenproblem
        let eig = SymmetricEigen::new(e);
        let max_idx = eig
            .eigenvalues
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(i, _)| i)
            .ok_or(MagCalibErrors::CoudlntGetMaxEigenValue)?;

        let max_val = eig.eigenvalues[max_idx];

        if max_val < -1e-6 {
            //in case of a perfect fit, the max eigenvalue can be slightly negative. Insane eh?
            return Err(MagCalibErrors::NegativeEigenvalues);
        }

        //v = [a, b, c, f, g, h, p, q, r, d].T
        let v1 = eig.eigenvectors.column(max_idx).clone_owned();
        let v2 = -s22_inv_times_s21 * &v1;

        let a = v1[0];
        let b = v1[1];
        let c = v1[2];
        let f = v1[3];
        let g = v1[4];
        let h = v1[5];
        let p = v2[0];
        let q = v2[1];
        let r = v2[2];
        let d = v2[3];

        #[allow(non_snake_case)]
        let I = a + b + c;
        #[allow(non_snake_case)]
        let J = a * b + b * c + a * c - h * h - g * g - f * f;
        if (4.0 * J - I * I) < 0.0 {
            return Err(MagCalibErrors::FittedNonEllipsoidQuadratic);
        }

        if a < 0.0 {
            return Ok(EllipsoidAlgebraicCoeffs {
                a: -a,
                b: -b,
                c: -c,
                h: -h,
                g: -g,
                f: -f,
                p: -p,
                q: -q,
                r: -r,
                d: -d,
            });
        } else {
            return Ok(EllipsoidAlgebraicCoeffs {
                a,
                b,
                c,
                h,
                g,
                f,
                p,
                q,
                r,
                d,
            });
        }
    }
}

pub struct SphericalCell {
    pub avg_radius: f64,
    pub radii: Vec<f64>,
    pub points: Vec<[f64; 3]>,
}

impl SphericalCell {
    pub fn new() -> Self {
        Self {
            avg_radius: 0.0,
            radii: Vec::new(),
            points: Vec::new(),
        }
    }
}

pub struct RadiiFilter {
    pub min: f64,
    pub max: f64,
}
//spherical downsampler cant be updated incrementally, because points may change bins
pub struct SphericalDownSampler {
    pub lon_div: usize,
    pub height_div: usize,
    pub buckets: BTreeMap<(usize, usize), SphericalCell>,
}

impl SphericalDownSampler {
    //if mean is not given, it will be calculated from the points
    pub fn new(lon_div: usize, height_div: usize) -> Self {
        Self {
            lon_div,
            height_div,
            buckets: BTreeMap::new(),
        }
    }

    pub fn bucket_is_filled(&self, lon_idx: usize, height_idx: usize) -> bool {
        self.buckets.contains_key(&(lon_idx, height_idx))
    }

    pub fn reset(&mut self) {
        self.buckets.clear();
    }

    //if s_in_bin is None, s_bin_by is used to fill the buckets
    pub fn add_points(
        &mut self,
        s_bin_by: &[[f64; 3]],
        center_ref: [f64; 3],
        s_in_bin: Option<&[[f64; 3]]>,
    ) -> Result<(), MagCalibErrors> {
        if let Some(s) = s_in_bin {
            if s.len() != s_bin_by.len() {
                return Err(MagCalibErrors::SizeMismatch);
            }
        }

        for (i, &[x, y, z]) in s_bin_by.iter().enumerate() {
            // Normalize
            let dx = x - center_ref[0];
            let dy = y - center_ref[1];
            let dz = z - center_ref[2];
            let norm = (dx * dx + dy * dy + dz * dz).sqrt();
            if norm == 0.0 {
                continue;
            } // skip degenerate points

            let ux = dx / norm;
            let uy = dy / norm;
            let uz = dz / norm;

            // Compute spherical coordinates
            let theta = (uy).atan2(ux); // azimuth

            let lon_idx = ((theta + PI) / (2.0 * PI) * self.lon_div as f64).floor() as usize;
            let height_idx = ((uz + 1.0) / 2.0 * self.height_div as f64).floor() as usize;

            let key = (
                lon_idx.min(self.lon_div - 1), //just for safety, we min.
                height_idx.min(self.height_div - 1),
            );

            let cell = self.buckets.entry(key).or_insert_with(SphericalCell::new);
            //wolfson's method for incremental mean
            let insert_idx: usize;
            if cell.points.len() == 0 {
                cell.avg_radius = norm;
                insert_idx = 0;
            } else {
                let delta = norm - cell.avg_radius;
                cell.avg_radius += delta / (cell.points.len() as f64 + 1.0);

                insert_idx = cell
                    .radii
                    .binary_search_by(|&x| {
                        x.partial_cmp(&norm)
                            .expect("Failed to compare radii, this should never happen")
                    })
                    .unwrap_or_else(|x| x); //if not found, insert at x
            }
            cell.radii.insert(insert_idx, norm);
            if let Some(s) = s_in_bin {
                cell.points.insert(insert_idx, s[i]);
            } else {
                cell.points.insert(insert_idx, s_bin_by[i]);
            }
        }
        Ok(())
    }

    pub fn sample_buckets(&self, radii_filter: Option<RadiiFilter>) -> Vec<[f64; 3]> {
        let mut target_points = Vec::new();
        for i in 0..self.lon_div {
            for j in 0..self.height_div {
                let key = (i, j);
                if let Some(cell) = self.buckets.get(&key) {
                    if cell.points.len() != 0 {
                        //was created but point was never pushed in
                        let target_radii: f64;
                        let target_point: [f64; 3];
                        if cell.points.len() > 3 {
                            //take most median point by radii
                            let mut median_idx = cell.points.len() / 2;
                            if cell.points.len() % 2 == 0 {
                                median_idx -= 1;
                            }
                            target_radii = cell.radii[median_idx];
                            target_point = cell.points[median_idx];
                        } else {
                            target_radii = cell.radii[0];
                            target_point = cell.points[0];
                        }

                        if let Some(filter) = &radii_filter {
                            if target_radii < filter.min || target_radii > filter.max {
                                continue;
                            }
                        }

                        target_points.push(target_point);
                    }
                }
            }
        }
        target_points
    }

    pub fn fill_fraction(&mut self) -> f64 {
        return self.buckets.len() as f64 / (self.lon_div * self.height_div) as f64;
    }
}

pub struct StatisticsAccumulator {
    mean: Vector3<f64>,
    cov_accumulator: Matrix3<f64>,
    count: usize,
    axis_min: Vector3<f64>,
    axis_max: Vector3<f64>,
}

impl StatisticsAccumulator {
    pub fn new() -> Self {
        Self {
            mean: Vector3::zeros(),
            cov_accumulator: Matrix3::zeros(),
            count: 0,
            axis_min: Vector3::new(f64::MAX, f64::MAX, f64::MAX),
            axis_max: Vector3::new(f64::MIN, f64::MIN, f64::MIN),
        }
    }

    pub fn add_points(&mut self, points: &[[f64; 3]]) {
        for p in points {
            self.add_point(Vector3::new(p[0], p[1], p[2]));
        }
    }

    pub fn add_point(&mut self, p: Vector3<f64>) {
        //using Welford's method for incremental mean and covariance
        self.count += 1;
        let delta = p - self.mean;
        self.mean += delta / self.count as f64;
        let delta2 = p - self.mean;

        if self.count > 1 {
            //notice that its delta*delta2!
            self.cov_accumulator += delta * delta2.transpose();
        }

        // Update axis min and max
        self.axis_min.x = self.axis_min.x.min(p.x);
        self.axis_min.y = self.axis_min.y.min(p.y);
        self.axis_min.z = self.axis_min.z.min(p.z);
        self.axis_max.x = self.axis_max.x.max(p.x);
        self.axis_max.y = self.axis_max.y.max(p.y);
        self.axis_max.z = self.axis_max.z.max(p.z);
    }
    pub fn mean(&self) -> Vector3<f64> {
        self.mean
    }

    pub fn cov(&self) -> Matrix3<f64> {
        if self.count < 2 {
            return Matrix3::zeros();
        }
        self.cov_accumulator / (self.count - 1) as f64
    }

    pub fn cov_condition_number(&self) -> f64 {
        if self.count < 2 {
            return f64::MAX;
        }

        let eig = self.cov().symmetric_eigen();
        let mut max = f64::MIN;
        let mut min = f64::MAX;
        for i in 0..3 {
            if eig.eigenvalues[i] > max {
                max = eig.eigenvalues[i];
            }
            if eig.eigenvalues[i] < min {
                min = eig.eigenvalues[i];
                if min <= 0.0 {
                    return f64::MAX; //this is not a valid covariance matrix, all eigenvalues should be positive
                }
            }
        }
        max / min
    }

    pub fn axis_ranges(&self) -> Vector3<f64> {
        self.axis_max - self.axis_min
    }

    pub fn normalize_by_mean_and_ranges(&self, raw_points: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let recip_ranges = self.axis_ranges().map(|x| 1.0 / x);

        let normalized_points = raw_points
            .iter()
            .map(|p| {
                let v3 = nalgebra::Vector3::new(p[0], p[1], p[2]);
                let normalized_v3 = (v3 - self.mean).component_mul(&recip_ranges);
                [normalized_v3.x, normalized_v3.y, normalized_v3.z]
            })
            .collect::<Vec<_>>();
        normalized_points
    }
}
