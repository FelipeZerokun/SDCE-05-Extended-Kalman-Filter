#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  // first, initialize the measurement matrix H for Laser *Laser Measurements Part 1
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // I am going to initialize the Jacobian matrix
  Hj_ << 1, 1, 0, 0,
		 1, 1, 0, 0,
		 1, 1, 1, 1;

    // Now, the state covariance matrix *Laser Measurements Part 3
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  // the state transition Matrix F without dT 
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;  

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /**
   * Initialization
   */
  if (!is_initialized_)
  {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      
      // I think the best way to get the Px and Py is with a right triangle 
      // As seen in the lesson, a right triangle is formed with Px, Py and Rho
      // Px is the opposite side, Py is the adjacent side and Rho is the Hypotenuse
      // so sin(phi) = Py/Rho and cos(phi) = Px/Rho
      // There is not information in the polar coordinates to get the Vx and Vy values, so I set them to 0
      float Px_ = rho * cos(phi);
      float Py_ = rho * sin(phi);
      float Vx_ = 0;
      float Vy_ = 0;
      ekf_.x_ << Px_, Py_ , Vx_ , Vy_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0; 
    }
	// Initializing the time stamp so when the dt calculation happens, there is a previous timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // set the acceleration noise components.
  int noise_ax = 9;
  int noise_ay = 9;  
  
  // Here I will create the time difference dt. Divide by 1E6 to have it in seconds.
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;  
  
  //Modify the F matrix so that the time is included
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // Pre compute terms for the Q matrix
  const float dt_2 = dt * dt;
  const float dt_3 = (dt_2 * dt)/2;
  const float dt_4 = (dt_3 * dt)/4;
  
  //set the process covariance matrix Q with the noise and dt terms
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4*noise_ax, 0, dt_3*noise_ax, 0,
              0, dt_4*noise_ay, 0, dt_3*noise_ay,
              dt_3*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  
  else
  {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
