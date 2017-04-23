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
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  n_x_=5;
  n_aug_=7;
  lambda_=3-n_aug_;
  // initial state vector
  x_ = VectorXd(n_x_);
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  //P_.setIdentity();
  //P_=P_*92;
  P_<<0.85,0,0,0,0,
      0,0.85,0,0,0,
      0,0,0.85,0,0,
      0,0,0,0.85,0,
      0,0,0,0,0.85;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;//30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.2;// 30;

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

  Xsig_pred_= MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {

  if (!is_initialized_) {
    // first measurement
    double px=0;
    double py=0;
    double phi=0;
    previous_timestamp_=measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      phi = measurement_pack.raw_measurements_[1];
      px =  measurement_pack.raw_measurements_[0]*cos(phi);
      py =  measurement_pack.raw_measurements_[0]*sin(phi);//py
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
    }
    if(px==0){
        px=0.0001;
    }
    if(py==0){
        py=0.0001;
    }
    x_<<px,py,0,0,0;
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  Prediction(dt);
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_==true) {
     UpdateRadar(measurement_pack.raw_measurements_);
  }
  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_==true) {
    UpdateLidar(measurement_pack.raw_measurements_);
  }
  previous_timestamp_ = measurement_pack.timestamp_;
}

void UKF::Prediction(double delta_t) {

  VectorXd x_aug = VectorXd(7);// Augmented state matrix;
  MatrixXd P_aug = MatrixXd(7, 7);//Augmented state covariance
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);//Augmented Sigma Points
  Xsig_pred_.fill(0.0);
  weights_ = VectorXd(2*n_aug_+1);
  GenerateAugmentedSigmaPoints(P_aug,x_aug,Xsig_aug);
  PredictSigmaPoints(delta_t,Xsig_aug);
  PredictMeanAndCovariance();

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF:: GenerateAugmentedSigmaPoints(MatrixXd &P_aug, VectorXd & x_aug, MatrixXd &Xsig_aug){
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  MatrixXd L = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_ ) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
}

void UKF::PredictSigmaPoints(double delta_t,MatrixXd &Xsig_aug){
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    //predicted state values
    double px_p, py_p;
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

}

void UKF::PredictMeanAndCovariance(){
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
  // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=NormalizeAngle(x_diff(3));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(VectorXd z) {
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
   //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = Tc * S.inverse();
  //residual
  VectorXd z_diff = z - z_pred;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  //Calculate Normalized Innovation for laser
  NIS_laser_=z_diff.transpose()*S.inverse()*z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(VectorXd z){
  int n_z=3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
   //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1)=NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1)=NormalizeAngle(z_diff(1));
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=NormalizeAngle(x_diff(3));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  //residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  z_diff(1)=NormalizeAngle(z_diff(1));
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  //Calculate Normalized Innovation Squared (NIS)
  NIS_radar_=z_diff.transpose()*S.inverse()*z_diff;

}

double UKF::NormalizeAngle(double angle){
  if(angle > M_PI){
    angle = angle - ceil(angle/(2*M_PI))*2*M_PI;
  }
  if(angle < -M_PI){
    angle = ceil(abs(angle)/(2*M_PI))*2*M_PI - abs(angle);
  }
  return angle;
}
