#include "FusionEKF.h"

#include "Eigen/Dense"
#include "tools.h"

#include <cmath>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */

  const auto &measurement = measurement_pack.raw_measurements_;
  const auto sensor = measurement_pack.sensor_type_;

  if (sensor != MeasurementPackage::LASER && sensor != MeasurementPackage::RADAR) {
    throw std::invalid_argument("Unsupported sensor type");
  }

  const int expected_size = (sensor == MeasurementPackage::LASER) ? 2 : 3;

  if (measurement.size() != expected_size || !measurement.allFinite()) {
    throw std::invalid_argument(
        "Measurement must have the expected dimension and contain finite values, but got: " +
        std::to_string(measurement.size()));
  }

  if (measurement_pack.timestamp_ < 0) {
    throw std::invalid_argument("Measurement timestamp must be nonnegative");
  }

  if (sensor == MeasurementPackage::RADAR && measurement(0) < 0.0) {
    throw std::invalid_argument("Radar measurement range must be non-negative, but got: " +
                                std::to_string(measurement(0)));
  }
  if (!is_initialized_) {
    ekf_.x_ = VectorXd::Zero(4);

    if (sensor == MeasurementPackage::LASER) {
      ekf_.x_(0) = measurement(0);
      ekf_.x_(1) = measurement(1);
    } else if (sensor == MeasurementPackage::RADAR) {
      const double rho = measurement(0);
      const double phi = measurement(1);

      ekf_.x_(0) = rho * std::cos(phi);
      ekf_.x_(1) = rho * std::sin(phi);
    }

    ekf_.P_ = MatrixXd::Zero(4, 4);
    ekf_.P_.diagonal() << 1.0, 1.0, 1000.0, 1000.0;

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  if (measurement_pack.timestamp_ < previous_timestamp_) {
    throw std::invalid_argument("Measurement timestamp must not move backward");
  }

  const double dt =
      static_cast<double>(measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  const double dt2 = dt * dt;
  const double dt3 = dt2 * dt;
  const double dt4 = dt2 * dt2;

  constexpr double noise_ax = 9.0;
  constexpr double noise_ay = 9.0;

  // Work on a candidate so a failed update preserves the accepted estimate.
  KalmanFilter candidate = ekf_;

  candidate.F_ = MatrixXd::Identity(4, 4);
  candidate.F_(0, 2) = dt;
  candidate.F_(1, 3) = dt;

  candidate.Q_ = MatrixXd::Zero(4, 4);
  candidate.Q_(0, 0) = dt4 * noise_ax / 4.0;
  candidate.Q_(0, 2) = dt3 * noise_ax / 2.0;
  candidate.Q_(2, 0) = candidate.Q_(0, 2);
  candidate.Q_(2, 2) = dt2 * noise_ax;

  candidate.Q_(1, 1) = dt4 * noise_ay / 4.0;
  candidate.Q_(1, 3) = dt3 * noise_ay / 2.0;
  candidate.Q_(3, 1) = candidate.Q_(1, 3);
  candidate.Q_(3, 3) = dt2 * noise_ay;

  candidate.Predict();

  if (!candidate.x_.allFinite() || !candidate.P_.allFinite()) {
    throw std::runtime_error("Prediction produced non-finite values");
  }

  if (sensor == MeasurementPackage::RADAR) {
    candidate.H_ = tools.CalculateJacobian(candidate.x_);
    candidate.R_ = R_radar_;
    candidate.UpdateEKF(measurement);
  } else {
    candidate.H_ = H_laser_;
    candidate.R_ = R_laser_;
    candidate.Update(measurement);
  }

  ekf_ = candidate;
  previous_timestamp_ = measurement_pack.timestamp_;
}
