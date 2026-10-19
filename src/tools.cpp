#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  /**
   * TODO: Calculate the RMSE here.
   */

  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // First, I'll check that the estimations are not 0
  //and the size is equal to the ground_truth vector size
  if (estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // Now, I will accumulate the residuals
  for (int i=0; i < estimations.size(); ++i)
  {
    // calculate the residual in each iteration
    VectorXd residual = estimations[i] - ground_truth[i];
    // coefficient-wise multiplication to make the residual a matrix
    residual = residual.array()*residual.array();
    // add up the current iteration residual to the rmse variable
    rmse += residual;
  }

  // After adding the residuals, I will calculate the mean
  rmse = rmse/estimations.size();

  // and finally the square root
  rmse = rmse.array().sqrt();

  // Now, I will return the RMSE
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  // First I need the x_state
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute repetitive terms in the equations
  float term1 = px*px + py*py; // (px^2 + py^2)
  float term2 = sqrt(term1);// (px^2 + py^2)^1/2
  float term3 = term1*term2; // (px^2 + py^2)^3/2

  // check division by zero
  if (fabs(term2) < 0.0001)
  {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/term2), (py/term2), 0, 0,
        (-py/term1), (px/term1), 0, 0,
        (py*vx*py - py*vy*px)/term3, (px*vy*px - px*vx*py)/term3, (px/term2), (py/term2);
  return Hj;
}
