#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_; //done

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_; //done

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_; //done

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_; //done

  ///* state covariance matrix
  MatrixXd P_; //done

  //augmented state
  VectorXd x_aug_;

  ///* augmented state covariance matrix
  MatrixXd P_aug_; //done

  ///* time when the state is true, in us
  long long time_us_; //done

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_ = 0.5; //done

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_ = 0.2; //done

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_ = 0.15; //done

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_ = 0.15; //done

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_ = 0.3; //done

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_ = 0.03; //done

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_ = 0.3; //done

  const double RSV = 0.001;

  //set measurement dimension, laser can measure px and py
  const int n_z_lid_ = 2;

  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z_rad_ = 3;

  //noise matrix for lidar update
  MatrixXd R_lid_ = MatrixXd(n_z_lid_,n_z_lid_);

  //noise matrix for radar update
  MatrixXd R_rad_ = MatrixXd(n_z_rad_,n_z_rad_);




  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_; //done

  ///* Augmented state dimension
  int n_aug_; //done

  ///* Sigma point spreading parameter
  double lambda_; //done
  double spread_vector_; //done

  float delta_t;
  float previous_timestamp_;

  // Matrix to store augmented sigma points
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_; 

  double NIS_radar_;
  double NIS_laser_;



  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
