#include "kalman_filter.h"

#include <Eigen/Cholesky>
#include <cmath>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in,
                        MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  const VectorXd innovation = z - H_ * x_;
  ApplyInnovation(innovation);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  if (z.size() != 3 || !z.allFinite() || z(0) < 0.0) {
    throw std::invalid_argument(
        "Radar measurement requires three finite values and nonnegative range");
  }

  if (x_.size() != 4 || !x_.allFinite()) {
    throw std::invalid_argument("Radar update requires four finite state components");
  }

  const double px = x_(0);
  const double py = x_(1);
  const double range = std::hypot(px, py);

  if (!std::isfinite(range)) {
    throw std::overflow_error("Predicted radar range is not finite");
  }

  // Must match the Jacobian's minimum-range policy.
  constexpr double minimum_range = 1e-6;
  if (range <= minimum_range) {
    throw std::domain_error("Radar update is undefined near zero range");
  }

  const double bearing = std::atan2(py, px);
  const double radial_velocity = (px / range) * x_(2) + (py / range) * x_(3);

  const Eigen::Vector3d predicted_measurement(range, bearing, radial_velocity);

  VectorXd innovation = z - predicted_measurement;

  const double two_pi = 2.0 * std::acos(-1.0);
  innovation(1) = std::remainder(innovation(1), two_pi);

  ApplyInnovation(innovation);
}

void KalmanFilter::ApplyInnovation(const VectorXd &innovation) {
  if (!innovation.allFinite()) {
    throw std::runtime_error("Measurement innovation is not finite");
  }

  const MatrixXd covariance_times_h_transpose = P_ * H_.transpose();
  const MatrixXd innovation_covariance = H_ * covariance_times_h_transpose + R_;

  const Eigen::LDLT<MatrixXd> decomposition(innovation_covariance);

  if (decomposition.info() != Eigen::Success || !decomposition.isPositive() ||
      (decomposition.vectorD().array() <= 0.0).any()) {
    throw std::runtime_error("Innovation covariance must be positive definite");
  }

  const MatrixXd gain = decomposition.solve(covariance_times_h_transpose.transpose()).transpose();

  const VectorXd updated_state = x_ + gain * innovation;
  const MatrixXd identity = MatrixXd::Identity(x_.size(), x_.size());
  const MatrixXd correction = identity - gain * H_;

  const MatrixXd updated_covariance =
      correction * P_ * correction.transpose() + gain * R_ * gain.transpose();

  if (!updated_state.allFinite() || !updated_covariance.allFinite()) {
    throw std::runtime_error("Measurement update produced non-finite values");
  }

  x_ = updated_state;
  P_ = updated_covariance;
}
