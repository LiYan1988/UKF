#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <vector>
#include <fstream>

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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  is_initialized_ = false;

  x_ << 0, 0, 0, 0, 0;

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  n_x_ = 5;

  n_aug_ = 7;

  time_us_ = 0;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i<2 * n_aug_ + 1; i++) weights_(i) = 0.5 / (lambda_ + n_aug_);

  double lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  vector <double> nis_laser_;
  vector <double> nis_radar_;

  ofstream outfile1("NIS_Laser.csv");
  outfile1.close();

  ofstream outfile2("NIS_Radar.csv");
  outfile2.close();
//  nis_laser_.clear();
//  nis_radar_.clear();
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

  if(!is_initialized_){
    cout<<"Initialize UKF."<<endl;
    if (meas_package.sensor_type_==MeasurementPackage::LASER) {
      // laser initialization
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      x_(0) = px;
      x_(1) = py;
    }
    else if (meas_package.sensor_type_==MeasurementPackage::RADAR) {
       // radar initialization
       float rho = meas_package.raw_measurements_[0];
       float phi = meas_package.raw_measurements_[1];
       x_(0) = rho * cos(phi);
       x_(1) = rho * sin(phi);
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  if ((meas_package.sensor_type_==MeasurementPackage::LASER)&use_laser_) {
    Prediction(dt);
    UpdateLidar(meas_package);
  }
  else if ((meas_package.sensor_type_==MeasurementPackage::RADAR)&use_radar_) {
    Prediction(dt);
    UpdateRadar(meas_package);
  }

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

  // construct matrices
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;

  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();

  // augmented sigma points before prediction
  for(int i=0; i<n_aug_; i++){
    Xsig_aug.col(i+1)         = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_)  = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

//  cout<<"Xsig_aug = \n" << Xsig_aug<<endl;

  // prediction
  for (int i=0; i<2*n_aug_+1; i++) {
    // extract values
    double px     = Xsig_aug(0, i);
    double py     = Xsig_aug(1, i);
    double v      = Xsig_aug(2, i);
    double yaw    = Xsig_aug(3, i);
    double yawd   = Xsig_aug(4, i);
    double vd     = Xsig_aug(5, i);
    double yawdd  = Xsig_aug(6, i);

    double dpx    = delta_t * delta_t * cos(yaw) * vd / 2;
    double dpy    = delta_t * delta_t * sin(yaw) * vd / 2;
    double dv     = delta_t * vd;
    double dyaw   = delta_t * delta_t * yawdd / 2 + yawd * delta_t;
    double dyawd  = delta_t * yawdd;

    if (yawd==0) {
      dpx += v * cos(yaw) * delta_t;
      dpy += v * sin(yaw) * delta_t;
    }
    else {
      dpx += (sin(yaw + yawd * delta_t) - sin(yaw)) * v / yawd;
      dpy += (-cos(yaw + yawd * delta_t) + cos(yaw)) * v / yawd;
    }
    Xsig_pred_.col(i) << px + dpx, py + dpy, v + dv, yaw + dyaw, yawd + dyawd;
  }

//  cout<<"dt = " << delta_t<<endl;
//  cout<<"Xsig_pred_ = \n" <<Xsig_pred_<<endl;

  // calculate mean and covariance
  x_ = VectorXd::Zero(n_x_);
  for (int i=0; i<2 * n_aug_ + 1; i++) x_ += weights_(i) * Xsig_pred_.col(i);
//  cout<< "x = "<< x_<<endl;

  P_ = MatrixXd::Zero(n_x_, n_x_);
  for (int i=0; i<2 * n_aug_ + 1; i++) {
    VectorXd v = Xsig_pred_.col(i) - x_;
    P_ += weights_(i) * v * v.transpose();
  }
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
  // measured position
  double px = meas_package.raw_measurements_[0];
  double py = meas_package.raw_measurements_[1];
  VectorXd z = VectorXd(2);
  z << px, py;

  // predicted measurement mean
  VectorXd z_pred_meas_mean = VectorXd::Zero(2);
  MatrixXd z_pred = MatrixXd::Zero(2, 2 * n_aug_ + 1);
  for (int i=0; i<2 * n_aug_ + 1; i++) {
    z_pred.col(i) << Xsig_pred_(0, i), Xsig_pred_(1, i);
    z_pred_meas_mean += weights_(i) * z_pred.col(i);
  }

  MatrixXd S = MatrixXd::Zero(2, 2);
  S << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  for (int i=0; i<2 * n_aug_ + 1; i++)
    S += weights_(i) * (z_pred.col(i) - z_pred_meas_mean) * (z_pred.col(i) - z_pred_meas_mean).transpose();

  MatrixXd T = MatrixXd::Zero(n_x_, 2);
  for (int i=0; i<2 * n_aug_ + 1; i++)
    T += weights_(i) * (Xsig_pred_.col(i) - x_) * (z_pred.col(i) - z_pred_meas_mean).transpose();

  MatrixXd K = T * S.inverse();
  x_ += K * (z - z_pred_meas_mean);
  P_ -= K * S * K.transpose();

  nis_laser_.push_back((z - z_pred_meas_mean).transpose() * S.inverse() * (z - z_pred_meas_mean));

  ofstream outfile_laser;
  outfile_laser.open("NIS_Laser.csv", std::ios_base::app);
  outfile_laser << nis_laser_.back() << ", "<<endl;
  outfile_laser.close();

//  cout<<"K = \n" <<K<<endl;
//  cout<<"T = \n" <<T<<endl;
//  cout<<"S = \n" <<S<<endl;
//  cout<<"z = \n" <<z<<endl;
//  cout<<"z_pred = \n"<<z_pred<<endl;
//  cout<<"z_pred_meas_mean = \n" <<z_pred_meas_mean<<endl;
//  cout<<"x = \n" <<x_<<endl;
//  cout<<"P = \n" <<P_<<endl;
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

  // measurement data
  VectorXd z = VectorXd::Zero(3);
  z<<meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  // predicted measurement mean
  VectorXd z_pred_meas_mean = VectorXd::Zero(3);
  MatrixXd z_pred = MatrixXd::Zero(3, 2 * n_aug_ + 1);
  for (int i=0; i<2 * n_aug_ + 1; i++) {
    double px   = Xsig_pred_(0, i);
    double py   = Xsig_pred_(1, i);
    double v    = Xsig_pred_(2, i);
    double yaw  = Xsig_pred_(3, i);
    double yawd = Xsig_pred_(4, i);

    double rho  = sqrt(px * px + py * py);
    double phi  = atan2(py, px);
    double rhod = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;

    z_pred.col(i) << rho, phi, rhod;
    z_pred_meas_mean += weights_(i) * z_pred.col(i);
  }

  MatrixXd S = MatrixXd::Zero(3, 3);
  S(0, 0) = std_radr_ * std_radr_;
  S(1, 1) = std_radphi_ * std_radphi_;
  S(2, 2) = std_radrd_ * std_radrd_;
  for (int i=0; i<2 * n_aug_ + 1; i++)
    S += weights_(i) * (z_pred.col(i) - z_pred_meas_mean) * (z_pred.col(i) - z_pred_meas_mean).transpose();

  MatrixXd T = MatrixXd::Zero(n_x_, 3);
  for (int i=0; i<2 * n_aug_ + 1; i++)
    T += weights_(i) * (Xsig_pred_.col(i) - x_) * (z_pred.col(i) - z_pred_meas_mean).transpose();

  MatrixXd K = T * S.inverse();
  x_ += K * (z - z_pred_meas_mean);
  P_ -= K * S * K.transpose();
  nis_radar_.push_back((z - z_pred_meas_mean).transpose() * S.inverse() * (z - z_pred_meas_mean));

  ofstream outfile_radar;
  outfile_radar.open("NIS_Radar.csv", std::ios_base::app);
  outfile_radar << nis_radar_.back() << ", "<<endl;
  outfile_radar.close();
}
