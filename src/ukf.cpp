#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double EPS = 0.0001; // Just a small number

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.4;//3;//1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.35;//0.7;//0.57;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  double weight = 0.5 / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = weight;
  }

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double ro = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      double ro_dot = meas_package.raw_measurements_(2);

      double px = ro * cos(theta);
      double py = ro * sin(theta);
      double vx = ro_dot * cos(theta);
      double vy = ro_dot * sin(theta);
      double v = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, theta, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      x_ << px, py, 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    cout << "initialized x_: " << x_ << endl;
    return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar updates
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    // Laser updates
    UpdateLidar(meas_package);
  }

  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  //set remaining sigma points
  MatrixXd s = sqrt(lambda_ + n_aug_) * L;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + s.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - s.col(i);
  }

  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(
    const MatrixXd& Xsig_aug, double delta_t, MatrixXd* Xsig_out) {

  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  double delta_t_sqr = delta_t * delta_t;
  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_aug = Xsig_aug.col(i);
    VectorXd x = VectorXd(5);
    double p_x = x_aug(0);
    double p_y = x_aug(1);
    double v = x_aug(2);
    double yaw = x_aug(3);
    double yawd = x_aug(4);
    double ni_a = x_aug(5);
    double ni_yawdd = x_aug(6);

    // Add noise
    x(0) = p_x + 0.5 * delta_t_sqr * cos(yaw) * ni_a;
    x(1) = p_y + 0.5 * delta_t_sqr * sin(yaw) * ni_a;
    x(2) = v + delta_t * ni_a;
    x(3) = yaw + 0.5 * delta_t_sqr * ni_yawdd;
    x(4) = yawd + delta_t * ni_yawdd;

    if (fabs(yawd) > EPS) {
        double v_yawd = v / yawd;
        x(0) += v_yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        x(1) += v_yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
        x(2) += 0;
        x(3) += yawd * delta_t;
        x(4) += 0;
    } else {
        x(0) += v * cos(yaw) * delta_t;
        x(1) += v * sin(yaw) * delta_t;
        x(2) += 0;
        x(3) += yawd * delta_t;
        x(4) += 0;
    }

    Xsig_pred.col(i) = x;
  }

  *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(
    const MatrixXd& Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    NormalizeAngle(&x_diff(3));

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  *x_out = x;
  *P_out = P;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate Sigma Points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  // input: x_, P_, output: Xsig_aug
  AugmentedSigmaPoints(&Xsig_aug);

  // Predict Sigma Points
  // input: Xsig_aug, output: Xsig_pred
  MatrixXd Xsig_pred;
  SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred);

  // Predict mean and covariance
  VectorXd x_pred = VectorXd(n_x_);
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  // input: Xsig_pred, output: x_, P_
  PredictMeanAndCovariance(Xsig_pred, &x_pred, &P_pred);

  x_ = x_pred;
  cout << "x_pred: x_: " << x_ << endl;
  P_ = P_pred;
  Xsig_pred_ = Xsig_pred;
}

void UKF::PredictMeanMeasurementAndCovarianceMetrix(
    MeasurementPackage::SensorType sensor_type,
    const MatrixXd& Zsig,
    VectorXd* z_out,
    MatrixXd* S_out) {
  int n_z = -1;
  if (sensor_type == MeasurementPackage::RADAR) {
    n_z = 3;
  } else if (sensor_type == MeasurementPackage::LASER) {
    n_z = 2;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  //calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(&z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }

  if (sensor_type == MeasurementPackage::RADAR) {
    S = S + R_radar_;
  } else if (sensor_type == MeasurementPackage::LASER) {
    S = S + R_lidar_;
  }

  *z_out = z_pred;
  *S_out = S;
}


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i=0; i<2 * n_aug_ + 1; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw  = Xsig_pred_(3, i);
    double yawd  = Xsig_pred_(4, i);

    double d = sqrt(px * px + py * py);

    if (abs(d) < EPS) {
      d = EPS;
    }

    if (abs(px) < EPS) {
      px = EPS;
    }

    Zsig(0, i) = d;
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / d;
  }
  *Zsig_out = Zsig;

  PredictMeanMeasurementAndCovarianceMetrix(MeasurementPackage::RADAR, Zsig, z_out, S_out);
}


void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i=0; i<2 * n_aug_ + 1; ++i) {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }
  *Zsig_out = Zsig;

  PredictMeanMeasurementAndCovarianceMetrix(MeasurementPackage::LASER, Zsig, z_out, S_out);
}

void UKF::UpdateUKF(
    MeasurementPackage::SensorType sensor_type,
    const VectorXd& z_measured,
    const MatrixXd& Zsig,
    const VectorXd& z_pred,
    const MatrixXd& S) {

  int n_z = -1;
  if (sensor_type == MeasurementPackage::RADAR) {
    n_z = 3;
  } else if (sensor_type == MeasurementPackage::LASER) {
    n_z = 2;
  }

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(&x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (sensor_type == MeasurementPackage::RADAR) {
      NormalizeAngle(&z_diff(1));
    }
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K * (z_measured - z_pred);
  P_ = P_ - K * S * K.transpose();
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z_pred;
  MatrixXd S;
  MatrixXd Zsig;
  // input: Xsig_pred_, output: Zsig
  // input: Zsig, output: z_pred
  // input: Zsig, z_pred, output: S
  PredictLidarMeasurement(&z_pred, &S, &Zsig);

  double px = meas_package.raw_measurements_(0);
  double py = meas_package.raw_measurements_(1);

  //create vector for incoming radar measurement
  VectorXd z_measured = VectorXd(2);
  z_measured << px, py;

  // update x_ and P_
  UpdateUKF(MeasurementPackage::LASER, z_measured, Zsig, z_pred, S);
  cout << "UKF::UpdateLidar: x_: " << x_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  VectorXd z_pred;
  MatrixXd S;
  MatrixXd Zsig;
  // input: Xsig_pred_, output: Zsig
  // input: Zsig, output: z_pred
  // input: Zsig, z_pred, output: S
  PredictRadarMeasurement(&z_pred, &S, &Zsig);

  double rho = meas_package.raw_measurements_(0);
  double theta = meas_package.raw_measurements_(1);
  double rho_dot = meas_package.raw_measurements_(2);

  //create vector for incoming radar measurement
  VectorXd z_measured = VectorXd(3);
  z_measured << rho, theta, rho_dot;

  // update x_ and P_
  UpdateUKF(MeasurementPackage::RADAR, z_measured, Zsig, z_pred, S);
  cout << "UKF::UpdateRadar: x_: " << x_ << endl;
}

void UKF::NormalizeAngle(double* angle) {
  while (*angle >  M_PI) {
    *angle -= 2. * M_PI;
  }
  while (*angle < -M_PI) {
    *angle += 2. * M_PI;
  }
}
