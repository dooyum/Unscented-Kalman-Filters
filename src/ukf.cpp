#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  n_x_ = 5;

  n_aug_ = 7;

  alpha_ = .001;

  beta_ = 2;

  lambda_ = alpha_ * alpha_ * n_aug_ - n_aug_;

  weight_initial_ = .5 / (n_aug_ + lambda_);

  // initialize state vector
  x_ = VectorXd::Zero(n_x_);

  // initialize covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  // initialize weights
  weights_m_ = VectorXd::Zero(2 * n_aug_ + 1);
  weights_c_ = VectorXd::Zero(2 * n_aug_ + 1);

  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
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
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_) {
    P_ << 2,0,0,0,0,
          0,4,0,0,0,
          0,0,1,0,0,
          0,0,0,0.5,0,
          0,0,0,0,0.5;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      x_ = tools.PolarToCartesian(meas_package.raw_measurements_, n_x_);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

      if (fabs(x_(0)) < 0.001 and fabs(x_(1)) < 0.001){
        x_(0) = 0;
        x_(1) = 0;
      }
    }

    // set weights
    double const weight_m_0 = lambda_ / (lambda_ + n_aug_);
    double const weight_c_0 = lambda_ / (lambda_ + n_aug_) + (1 - alpha_ * alpha_ + beta_);
    weights_m_.fill(weight_initial_);
    weights_c_.fill(weight_initial_);
    weights_m_(0) = weight_m_0;
    weights_c_(0) = weight_c_0;
    
    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  float const delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1e6;
  previous_timestamp_ = meas_package.timestamp_;
  bool const is_radar = meas_package.sensor_type_ == MeasurementPackage::RADAR;
  Prediction(delta_t, is_radar);
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * @param {boolean} is_radar whether this measurement is from radar
 */
void UKF::Prediction(double delta_t, bool is_radar) {
  /**
  Estimates the object's location.
  Predict sigma points, the state, and the state covariance matrix.
  */

  // First generate and augment sigma points
  MatrixXd Xsig_aug = AugmentSigmaPoints();
  
  //predict sigma points
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double const p_x = Xsig_aug(0,i);
    double const p_y = Xsig_aug(1,i);
    double const v = Xsig_aug(2,i);
    double const yaw = Xsig_aug(3,i);
    double const yawd = Xsig_aug(4,i);
    double const nu_a = Xsig_aug(5,i);
    double const nu_yawdd = Xsig_aug(6,i);
    
    //predicted state values
    double px_p = 0;
    double py_p = 0;
    
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;
    
    //add noise
    px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;
    
    // write predicted sigma points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  PredictMeanAndCovariance(is_radar);
}

MatrixXd UKF::AugmentSigmaPoints() {
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  //calculate square root of P_aug
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  return Xsig_aug;
}

void UKF::PredictMeanAndCovariance(bool is_radar) {
  // Predicted state mean
  VectorXd x = VectorXd::Zero(n_x_);
  x = Xsig_pred_ * weights_m_; // vectorised sum

  //predicted state covariance matrix
  MatrixXd P = MatrixXd::Zero(n_x_, n_x_);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
//    if (is_radar) x_diff = tools.NormalizeRadians(x_diff, 3);
    P = P + weights_c_(i) * x_diff * x_diff.transpose();
  }
  x_ = x;
  P_ = P;
}


void UKF::UpdateState(VectorXd z, bool is_radar) {
  int const n_z = is_radar ? 3 : 2;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //residual
    VectorXd z_diff = Zsig_pred_.col(i) - z_pred_;
    if (is_radar) z_diff = tools.NormalizeRadians(z_diff, 1);

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    if (is_radar) x_diff = tools.NormalizeRadians(x_diff, 3);

    Tc = Tc + weights_c_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_;
  if (is_radar) z_diff = tools.NormalizeRadians(z_diff, 1);

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_ * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   Using lidar data to update the belief about the object's
   position, modify the state vector, x_, and covariance, P_.
   */
  PredictMeasurement(false);
  VectorXd z = VectorXd::Zero(2);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  UpdateState(z, false);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Using radar data to update the belief about the object's
  position, modify the state vector, x_, and covariance, P_.
   */
  PredictMeasurement(true);
  VectorXd z = VectorXd::Zero(3);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  UpdateState(z, true);
}

/**
 * @param {bool} is_radar is true if the measurement to be predicted is radar.
 */
void UKF::PredictMeasurement(bool is_radar) {
  //set measurement dimension, radar can measure rho, phi, and rho_dot
  int const n_z = is_radar ? 3 : 2;

  //create matrix for sigma points in predicted measurement space
  Zsig_pred_ = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    double const p_x = Xsig_pred_(0,i);
    double const p_y = Xsig_pred_(1,i);
    double const v  = Xsig_pred_(2,i);
    double const yaw = Xsig_pred_(3,i);
    double const v1 = cos(yaw) * v;
    double const v2 = sin(yaw) * v;

    if (is_radar){
      Zsig_pred_(0,i) = sqrt(p_x * p_x + p_y * p_y); //rho
      Zsig_pred_(1,i) = atan2(p_y, p_x); //phi
      Zsig_pred_(2,i) = (p_x * v1 + p_y * v2 ) / sqrt(p_x * p_x + p_y * p_y); //rho_dot
    } else {
      Zsig_pred_(0,i) = p_x;
      Zsig_pred_(1,i) = p_y;
    }
  }

  //mean predicted measurement
  z_pred_ = VectorXd::Zero(n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_ = z_pred_ + weights_m_(i) * Zsig_pred_.col(i);
  }

  //innovation covariance matrix S
  S_ = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig_pred_.col(i) - z_pred_;
    if (is_radar) z_diff = tools.NormalizeRadians(z_diff, 1);
    S_ = S_ + weights_c_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  if (is_radar) {
    R << std_radr_ * std_radr_, 0, 0,
         0, std_radphi_ * std_radphi_, 0,
         0, 0, std_radrd_ * std_radrd_;
  } else {
    R <<  std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;
  }
  S_ = S_ + R;
}
