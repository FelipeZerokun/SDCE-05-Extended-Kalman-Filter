#ifndef TOOLS_H_
#define TOOLS_H_

#include "Eigen/Dense"

#include <vector>

class Tools {
public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * Computes component-wise RMSE for paired states [px, py, vx, vy].
   * Returns position errors in metres and velocity errors in metres/second.
   * Inputs are not modified; the caller must align corresponding observations.
   *
   * Throws std::invalid_argument if inputs are empty, sequence lengths differ,
   * any state has other than four components, or any input is non-finite.
   * Throws std::overflow_error if arithmetic produces a non-finite value.
   *
   * Uses direct squared-error accumulation. Extreme finite inputs may overflow
   * even when the mathematical RMSE is representable. Tiny squared errors may
   * underflow to zero.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd &x_state);
};

#endif // TOOLS_H_
