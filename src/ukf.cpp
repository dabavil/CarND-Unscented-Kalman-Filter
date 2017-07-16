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
  
  //Set to false when the instance is created
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // ***********For prediction state
  // State dim
  n_x_ = 5;

  // Augmented state dim
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Lambda and spread vector calc
  lambda_ = 3 - n_aug_;
  spread_vector_ = sqrt(lambda_ + n_aug_);

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0);

  //augmented state
  x_aug_ = VectorXd(n_aug_);

  // augmented state covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // augmented sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //predicted sig points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


  //matrix for lidar update
  R << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  

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

  //Initialize first

  if(!is_initialized_)
  {

    cout<<"Init started...\n";
    //check which time of measurement we have

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    cout<<"Covar matrix P initialized...\n";

    //Initialize time difference
    double delta_t = double(meas_package.timestamp_ - time_us_) / 1000000;
    //Update to most recent time stamp
      time_us_ = meas_package.timestamp_;


    //Initialize weights
    weights_(0) = lambda_/(lambda_ + n_aug_);
    for (int i=1; i<2*n_aug_+1; i++)
      {  //2n+1 weights
      weights_(i) = 0.5/(n_aug_+lambda_);
      }
    cout<<"weights_ initialized...\n";
    cout<<weights_<<endl;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) 
    {
      // TODO Initialize using radar coordinates

      cout<<"Initializing based on radar measurement...\n";

      float rho = meas_package.raw_measurements_[0]; // Range 
      float phi = meas_package.raw_measurements_[1]; // Bearing 
      float rho_dot = meas_package.raw_measurements_(2);

      float x_in = rho * cos(phi);
      float y_in = rho * sin(phi);
      x_ << x_in, y_in, 0, 0, 0;


    } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {
      // TODO Initialize using lidar coordinates

      cout<<"Initializing based on lidar measurement...\n";

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    

    if (fabs(x_(0)) < RSV and fabs(x_(1)) < RSV) 
    {
      cout<<"Adjusting for 0 x and y...\n";
      x_(0) = RSV;
      x_(1) = RSV;
    }

    
    is_initialized_ = true;
    cout<<"INITIALIZATION COMPLETE."<<endl;

     
  } //End of initialization
  else
  {

    delta_t = double(meas_package.timestamp_ - time_us_) / 1000000;

    time_us_ = meas_package.timestamp_;


  Prediction(delta_t);


  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) 
  {
    
    //Measurement update RADAR
    UpdateRadar(meas_package);
    
    
  } else 
    {  
    //Measurement update LIDAR
    UpdateLidar(meas_package);
    }

  } 

cout<<"End of process measurement cycle."<<endl;

} // End of process measurement function


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

// First, calculate the Sigma points
  //*************************************************************
//state and covariance must already be initialized
  
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0.0;
  x_aug_(6) = 0.0;
  
  //create augmented covariance matrix
  
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd L_ = P_aug_.llt().matrixL();


  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  
  for(int i = 0; i < n_aug_; i++)
  {
    Xsig_aug_.col(i+1) = x_aug_ + spread_vector_ * L_.col(i);
    Xsig_aug_.col(i+8) = x_aug_ - spread_vector_ * L_.col(i);
  }

  //now sigma points are stored in the Xsig_aug matrix
  cout<<"Xsig_aug_:  \n"<<Xsig_aug_<<"\n";
  //exit(1);

  
  //This section projects the sigma points through the process model
  //*************************************************************
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    const double p_x = Xsig_aug_(0,i);
    const double p_y = Xsig_aug_(1,i);
    const double v = Xsig_aug_(2,i);
    const double yaw = Xsig_aug_(3,i);
    const double yawd = Xsig_aug_(4,i);
    const double nu_a = Xsig_aug_(5,i);
    const double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values for x and y position
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) < RSV) 
    {
      px_p = p_x + v * cos(yaw) * delta_t;
      py_p = p_y + v * sin(yaw) * delta_t;
    }
    else 
    {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    //remaining predicted state values
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;


    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;


    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  


  cout<<"delta_t: "<<delta_t<<endl;
  cout<<"Xsig_pred_:  \n"<<Xsig_pred_<<"\n";
  

  // now we have a complete set of predicted sigma points stored in the Xsig_pred_ matrix


  //Next up - calc the mean and normal distr. from these predicted sig points
  //*******************************************


  //predicted state mean
  x_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)  //iterate over sigma points
    {  
    x_ += weights_(i) * Xsig_pred_.col(i);
    }

  //predicted state covariance matrix
  P_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
    { 
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }
cout<<"PREDICTION step done. This is the result: "<<endl;
cout<<"x_ : "<<x_<<endl;
cout<<"P_: "<<P_<<endl;
//exit(1);
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

  cout<<"Lidar update started."<<endl;

  
  
  
  //set the measurements
  VectorXd z = meas_package.raw_measurements_;
  
  //initialize Zsig
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  
  //initialize S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  
  //initialize measurement prediction vector
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //transform sig points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over simga points

    //sigma point predictions in process space
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);

    //sigma point predictions in measurement space
    Zsig(0,i) = px;                   
    Zsig(1,i) = py;                                 
  }

  //mean predicted measurement
  z_pred = Zsig * weights_; // Matrix Calc (2X15) * (15X1) = (2X1)

  //measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix

  
  S = S + R;
  
  //Now Update the State x_ and state covariance P_
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
  //NIS Lidar Update
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff; 





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

  cout<<"Radar update started."<<endl;

//set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    // extract values for better readibility
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

  cout<<"Projection into measurement model.";

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++)
    {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;


////********************************************
  //Now update the state vector and covariance P

VectorXd z = VectorXd(n_z);

float rho = meas_package.raw_measurements_[0]; // Range 
float phi = meas_package.raw_measurements_[1]; // Bearing 
float rho_dot = meas_package.raw_measurements_[2]; 

z << rho, phi, rho_dot; 



MatrixXd Tc = MatrixXd(n_x_, n_z);


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();


  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;



}


