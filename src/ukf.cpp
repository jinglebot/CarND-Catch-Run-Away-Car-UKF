#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // iteration
  iteration = 0;

  // initialization
  is_initialized_ = false;

  // time step
  time_us_ = 0.0;

  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = n_x_ + 2;

  // sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // predicted sigma state dimension
  n_sig_ = 2 * n_aug_ + 1;

  // set weights
  weights_ = VectorXd(n_sig_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // noise covariance matrices
  // laser
  R_laser = MatrixXd(2, 2);
  R_laser <<  std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  // radar
  R_radar = MatrixXd(3, 3);
  R_radar <<  std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

void UKF::Normalize_angle(double a) {
  double a1 = a;
  while (a > M_PI) a -= 2. * M_PI;
  while (a < -M_PI) a += 2. * M_PI;
}

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
    // state vector
    x_.fill(0.0);
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0;

    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // initial state
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];

      // initial covariance matrix
      P_ << std_laspx_, 0.0, 0.0, 0.0, 0.0,
            0.0, std_laspy_, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0;
    } 
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // initial state
      double rho = meas_package.raw_measurements_[0];
      double theta = meas_package.raw_measurements_[1];
      
      x_[0] = rho * cos(theta);
      x_[1] = rho * sin(theta);

      // initial covariance matrix
      P_ << std_radr_, 0.0, 0.0, 0.0, 0.0,
            0.0, std_radphi_, 0.0, 0.0, 0.0,
            0.0, 0.0, std_radrd_, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0;
    }
    x_[2] = 1;
    x_[3] = 1;
    x_[4] = 0.1;

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }
  
  // prediction
  double delta_t = fabs(meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);
  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) UpdateLidar(meas_package);
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) UpdateRadar(meas_package);          
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

  // augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;

  // augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

  // root square of P_aug_
  MatrixXd L_ = P_aug_.llt().matrixL();

  // sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // sigma points matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_sig_);
  Xsig_aug_.fill(0.0);
  Xsig_aug_.col(0) = x_aug_;
  for (size_t i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L_.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L_.col(i);
  }

  // sigma points prediction
  for (size_t i = 0; i < n_sig_; i++) {
    double p_x        = Xsig_aug_(0, i);
    double p_y        = Xsig_aug_(1, i);
    double v          = Xsig_aug_(2, i);
    double yaw        = Xsig_aug_(3, i);
    double yawd       = Xsig_aug_(4, i);
    double nu_a       = Xsig_aug_(5, i);
    double nu_yawdd   = Xsig_aug_(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // predicted sigma points matrix
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (size_t i = 0; i < n_sig_; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // state mean prediction
  x_.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // state covariance prediction
  P_.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    Normalize_angle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
            
  // cout << "delta_t: " << delta_t << endl;
  // cout << "x_P: " << x_ << endl;
  // cout << "P_: " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // lidar measurement dimension
  int n_z_ = 2;

  // sigma points matrix in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);
  for (size_t i = 0; i < n_sig_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    VectorXd z = VectorXd(n_z_);
    z << px, py;
    Zsig.col(i) = z;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    VectorXd z_diff = VectorXd(n_z_);
    z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_laser;

  // cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // kalman gain
  MatrixXd K = Tc * S.inverse();

  // updated state mean
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  x_ = x_ + K * z_diff;

  // updated state covariance matrix
  P_ = P_ - K * S * K.transpose();

  // cout << "x_L: " << x_ << endl;
  // cout << "P_: " << P_ << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // radar measurement dimension
  int n_z_ = 3;

  // sigma points matrix in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);
  for (size_t i = 0; i < n_sig_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double psi = Xsig_pred_(3, i);

    double rho = sqrt(px * px + py * py);
    if (rho < 0.0001) rho = 0.01;
    double theta = atan2(py, px);
    double rho_dot = (px * cos(psi) * v + py * sin(psi) * v) / rho;

    VectorXd z = VectorXd(n_z_);
    z << rho, theta, rho_dot;
    Zsig.col(i) = z;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix 
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    Normalize_angle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_radar;

  // cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    Normalize_angle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Normalize_angle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // kalman gain
  MatrixXd K = Tc * S.inverse();

  // update state mean
  VectorXd z = VectorXd(n_z_);
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  Normalize_angle(z_diff(1));
  x_ = x_ + K * z_diff;

  // update state covariance matrix
  P_ = P_ - K * S * K.transpose();

  // cout << "x_R: " << x_ << endl;
  // cout << "P_: " << P_ << endl;
}


