#include "kalman_filter.h"
#include "tools.h"

#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("Prediction advances position and propagates uncertainty", "[kalman_filter]") {
  KalmanFilter filter;

  filter.x_ = Eigen::Vector4d(1.0, 2.0, 3.0, -4.0);
  filter.P_ = Eigen::Matrix4d::Identity();

  // A time step of two seconds.
  filter.F_ = Eigen::Matrix4d::Identity();
  filter.F_(0, 2) = 2.0;
  filter.F_(1, 3) = 2.0;

  filter.Q_ = 0.5 * Eigen::Matrix4d::Identity();

  filter.Predict();

  const Eigen::Vector4d expected_state(7.0, -6.0, 3.0, -4.0);

  Eigen::Matrix4d expected_covariance;
  expected_covariance << 5.5, 0.0, 2.0, 0.0, 0.0, 5.5, 0.0, 2.0, 2.0, 0.0, 1.5, 0.0, 0.0, 2.0, 0.0,
      1.5;

  REQUIRE(filter.x_.isApprox(expected_state, 1e-12));
  REQUIRE(filter.P_.isApprox(expected_covariance, 1e-12));
}

TEST_CASE("Identity prediction preserves state and covariance", "[kalman_filter]") {
  KalmanFilter filter;

  const Eigen::Vector4d initial_state(-2.0, 5.0, 1.5, -0.5);

  Eigen::Matrix4d initial_covariance;
  initial_covariance << 4.0, 0.0, 1.0, 0.0, 0.0, 9.0, 0.0, -2.0, 1.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0,
      3.0;

  filter.x_ = initial_state;
  filter.P_ = initial_covariance;
  filter.F_ = Eigen::Matrix4d::Identity();
  filter.Q_ = Eigen::Matrix4d::Zero();

  filter.Predict();

  REQUIRE(filter.x_.isApprox(initial_state, 1e-12));
  REQUIRE(filter.P_.isApprox(initial_covariance, 1e-12));
}

TEST_CASE("Lidar update corrects position and correlated velocity", "[kalman_filter]") {
  KalmanFilter filter;
  filter.x_ = Eigen::Vector4d(0.0, 0.0, 1.0, -1.0);

  filter.P_ = Eigen::Matrix4d::Identity();
  filter.P_(0, 2) = 0.5;
  filter.P_(2, 0) = 0.5;
  filter.P_(1, 3) = 0.5;
  filter.P_(3, 1) = 0.5;

  filter.H_ = Eigen::MatrixXd::Zero(2, 4);
  filter.H_(0, 0) = 1.0;
  filter.H_(1, 1) = 1.0;
  filter.R_ = Eigen::Matrix2d::Identity();

  filter.Update(Eigen::Vector2d(2.0, -4.0));

  const Eigen::Vector4d expected_state(1.0, -2.0, 1.5, -2.0);

  Eigen::Matrix4d expected_covariance;
  expected_covariance << 0.5, 0.0, 0.25, 0.0, 0.0, 0.5, 0.0, 0.25, 0.25, 0.0, 0.875, 0.0, 0.0, 0.25,
      0.0, 0.875;

  REQUIRE(filter.x_.isApprox(expected_state, 1e-12));
  REQUIRE(filter.P_.isApprox(expected_covariance, 1e-12));
}

TEST_CASE("Matching lidar measurement reduces uncertainty without moving state",
          "[kalman_filter]") {
  KalmanFilter filter;
  const Eigen::Vector4d initial_state(2.0, -3.0, 1.0, -1.0);

  filter.x_ = initial_state;
  filter.P_ = Eigen::Matrix4d::Identity();
  filter.H_ = Eigen::MatrixXd::Zero(2, 4);
  filter.H_(0, 0) = 1.0;
  filter.H_(1, 1) = 1.0;
  filter.R_ = Eigen::Matrix2d::Identity();

  filter.Update(Eigen::Vector2d(2.0, -3.0));

  Eigen::Matrix4d expected_covariance = Eigen::Matrix4d::Identity();
  expected_covariance(0, 0) = 0.5;
  expected_covariance(1, 1) = 0.5;

  REQUIRE(filter.x_.isApprox(initial_state, 1e-12));
  REQUIRE(filter.P_.isApprox(expected_covariance, 1e-12));
}

TEST_CASE("Matching radar measurement preserves state and reduces uncertainty", "[kalman_filter]") {
  KalmanFilter filter;
  Tools tools;

  const Eigen::Vector4d initial_state(2.0, 0.0, 0.0, 0.0);

  filter.x_ = initial_state;
  filter.P_ = Eigen::Matrix4d::Identity();
  filter.H_ = tools.CalculateJacobian(filter.x_);
  filter.R_ = Eigen::Matrix3d::Identity();

  // Range 2, bearing 0, radial velocity 0.
  filter.UpdateEKF(Eigen::Vector3d(2.0, 0.0, 0.0));

  Eigen::Matrix4d expected_covariance = Eigen::Matrix4d::Zero();
  expected_covariance.diagonal() << 0.5, 0.8, 0.5, 1.0;

  REQUIRE(filter.x_.isApprox(initial_state, 1e-12));
  REQUIRE(filter.P_.isApprox(expected_covariance, 1e-12));
}

TEST_CASE("Radar update normalizes bearing innovation across the angle boundary",
          "[kalman_filter]") {
  Tools tools;
  const double pi = std::acos(-1.0);
  constexpr double offset = 0.01;

  for (const double direction : {-1.0, 1.0}) {
    CAPTURE(direction);

    KalmanFilter filter;
    const double predicted_bearing = direction * (pi - offset);

    const Eigen::Vector4d initial_state(std::cos(predicted_bearing), std::sin(predicted_bearing),
                                        0.0, 0.0);

    filter.x_ = initial_state;
    filter.P_ = Eigen::Matrix4d::Identity();
    filter.H_ = tools.CalculateJacobian(filter.x_);
    filter.R_ = Eigen::Matrix3d::Identity();

    const double measured_bearing = -direction * (pi - offset);
    filter.UpdateEKF(Eigen::Vector3d(1.0, measured_bearing, 0.0));

    // At unit range, the bearing gain is half the unit tangent.
    // The normalized bearing innovation is direction * 2 * offset.
    Eigen::Vector4d expected_state = initial_state;
    expected_state(0) -= direction * offset * initial_state(1);
    expected_state(1) += direction * offset * initial_state(0);

    REQUIRE(filter.x_.isApprox(expected_state, 1e-12));
  }
}