#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  VectorXd x_ = x_in;
  MatrixXd P_ = P_in;
  MatrixXd F_ = F_in;
  MatrixXd H_ = H_in;
  MatrixXd R_ = R_in;
  MatrixXd Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_; // First, predict the new State Vector
  MatrixXd Ft_ = F_.transpose(); // Create the transpose of the Matrix F
  P_ = F_ * P_ * Ft_ + Q_; 
}

void KalmanFilter::Update(const VectorXd &z)
{
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_predict_ = H_ * x_;
  VectorXd y_ = z - z_predict_;
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S_ = H_*P_*Ht_ + R_;
  MatrixXd S_inv_ = S_.inverse();
  MatrixXd PHt_ = P_*Ht_;
  MatrixXd K_ = PHt_*S_inv_;

  // New Estimate Equations
  x_ = x_ + (K_*y_);  // The actual new State Vector
  MatrixXd I = MatrixXd::Identity(x_.rows(), x_.rows());
  MatrixXd KH = K_*H_;
  P_ = (I - KH)*P_;  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  float rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));  //(px^2 + py^2)^1/2
  float phi = atan2(x_(1), x_(0));                  // arc tan (py/px)
  float rho_dot = 0;

  // For rho_dot, a division by zero would imply that rho is zero.
  // If rho is Zero, I will skip the update step
  if (fabs(rho) < 0.0001)
  {
    return;
  }
  else
  {
    rho_dot = (x_(0) * x_(2) + x_(1) * x_(3))/rho;  //(px*vx + py+vy/)/rho
  }
  
 // Now, I will define the predicted state
  VectorXd z_predict_(3);
  z_predict_ << rho, phi, rho_dot;
  // cout << z_predict_;
  // now, lets see the error between real and predicted sensor readings
  VectorXd y_ = z - z_predict_;

  // Normalizing the angles
  // In C++, atan2() returns values between -pi and pi. 
  // When calculating phi in y = z - h(x) for radar measurements
  // the resulting angle phi in the y vector should be adjusted so that it is between -pi and pi.
  
  // HINT: when working in radians, you can add 2π2\pi2π or subtract 2π2\pi2π
  // until the angle is within the desired range.
  
  // First, lets check angles below PI
  float pi = 3.14159;
  if (y_(1) < -pi)
  {
    while(y_(1) < -pi)
    {
      y_(1) = y_(1) + 2 * pi;
    }
  }
  // same for the angle above PI

  if( y_(1) > pi)
  {
    while(y_(1) > pi)
    {
      y_(1) = y_(1) - 2 * pi;
    }
  }
  
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht_ + R_;
  MatrixXd S_inv_ = S_.inverse();
  MatrixXd PHt = P_ * Ht_;
  MatrixXd K_ = PHt * S_inv_;

  //new estimate
  x_ = x_ + (K_ * y_);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  MatrixXd KH = K_*H_;
  P_ = (I - KH)*P_;  
  
}
