#include "FusionEKF.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
MeasurementPackage MakeLidar(double px, double py, long long timestamp) {
  MeasurementPackage measurement;
  measurement.sensor_type_ = MeasurementPackage::LASER;
  measurement.timestamp_ = timestamp;
  measurement.raw_measurements_ = Eigen::Vector2d(px, py);
  return measurement;
}
} // namespace

TEST_CASE("First lidar measurement initializes state and covariance", "[fusion_ekf]") {
  FusionEKF fusion;

  MeasurementPackage measurement;
  measurement.sensor_type_ = MeasurementPackage::LASER;
  measurement.timestamp_ = 1000000;
  measurement.raw_measurements_ = Eigen::Vector2d(1.25, -2.5);

  fusion.ProcessMeasurement(measurement);

  const auto &state = fusion.ekf_.x_;
  REQUIRE(state.size() == 4);
  REQUIRE(state(0) == 1.25);
  REQUIRE(state(1) == -2.5);
  REQUIRE(state(2) == 0.0);
  REQUIRE(state(3) == 0.0);

  const auto &covariance = fusion.ekf_.P_;
  REQUIRE(covariance.rows() == 4);
  REQUIRE(covariance.cols() == 4);

  Eigen::Matrix4d expected = Eigen::Matrix4d::Zero();
  expected.diagonal() << 1.0, 1.0, 1000.0, 1000.0;
  REQUIRE(covariance.isApprox(expected, 1e-12));
}

TEST_CASE("First radar measurement initializes Cartesian position", "[fusion_ekf]") {
  FusionEKF fusion;

  MeasurementPackage measurement;
  measurement.sensor_type_ = MeasurementPackage::RADAR;
  measurement.timestamp_ = 1000000;

  // Range 2, bearing pi/2, radial velocity 3.
  const double half_pi = std::acos(-1.0) / 2.0;
  measurement.raw_measurements_ = Eigen::Vector3d(2.0, half_pi, 3.0);

  fusion.ProcessMeasurement(measurement);

  const auto &state = fusion.ekf_.x_;
  REQUIRE(state.size() == 4);

  REQUIRE_THAT(state(0), Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(state(1), Catch::Matchers::WithinAbs(2.0, 1e-12));
  REQUIRE(state(2) == 0.0);
  REQUIRE(state(3) == 0.0);

  const auto &covariance = fusion.ekf_.P_;
  REQUIRE(covariance.rows() == 4);
  REQUIRE(covariance.cols() == 4);

  Eigen::Matrix4d expected = Eigen::Matrix4d::Zero();
  expected.diagonal() << 1.0, 1.0, 1000.0, 1000.0;
  REQUIRE(covariance.isApprox(expected, 1e-12));
}

TEST_CASE("Invalid initialization allows retry with a valid measurement", "[fusion_ekf]") {
  FusionEKF fusion;

  MeasurementPackage invalid;
  invalid.timestamp_ = 1000000;

  SECTION("Lidar measurement has incorrect dimensions") {
    invalid.sensor_type_ = MeasurementPackage::LASER;
    invalid.raw_measurements_ = Eigen::Vector3d(1.0, 2.0, 3.0);
  }

  SECTION("Radar range is negative") {
    invalid.sensor_type_ = MeasurementPackage::RADAR;
    invalid.raw_measurements_ = Eigen::Vector3d(-1.0, 0.0, 0.0);
  }

  REQUIRE_THROWS_AS(fusion.ProcessMeasurement(invalid), std::invalid_argument);

  MeasurementPackage valid;
  valid.sensor_type_ = MeasurementPackage::LASER;
  valid.timestamp_ = 2000000;
  valid.raw_measurements_ = Eigen::Vector2d(4.0, -3.0);

  REQUIRE_NOTHROW(fusion.ProcessMeasurement(valid));

  REQUIRE(fusion.ekf_.x_.size() == 4);
  REQUIRE(fusion.ekf_.x_(0) == 4.0);
  REQUIRE(fusion.ekf_.x_(1) == -3.0);
  REQUIRE(fusion.ekf_.x_(2) == 0.0);
  REQUIRE(fusion.ekf_.x_(3) == 0.0);
}

TEST_CASE("Two lidar measurements perform prediction and correction", "[fusion_ekf]") {
  FusionEKF fusion;
  fusion.ProcessMeasurement(MakeLidar(1.0, 2.0, 0));
  fusion.ProcessMeasurement(MakeLidar(2.0, 0.0, 1000000));

  // After one second, each position/velocity covariance block is:
  // [1003.25, 1004.5]
  // [1004.5,  1009.0]
  const double s = 1003.25 + 0.0225;
  const double position_gain = 1003.25 / s;
  const double velocity_gain = 1004.5 / s;

  const Eigen::Vector4d expected_state(1.0 + position_gain, 2.0 - 2.0 * position_gain,
                                       velocity_gain, -2.0 * velocity_gain);

  Eigen::Matrix4d expected_covariance = Eigen::Matrix4d::Zero();
  expected_covariance(0, 0) = 1003.25 * 0.0225 / s;
  expected_covariance(1, 1) = expected_covariance(0, 0);
  expected_covariance(0, 2) = 1004.5 * 0.0225 / s;
  expected_covariance(2, 0) = expected_covariance(0, 2);
  expected_covariance(1, 3) = expected_covariance(0, 2);
  expected_covariance(3, 1) = expected_covariance(0, 2);
  expected_covariance(2, 2) = 1009.0 - 1004.5 * 1004.5 / s;
  expected_covariance(3, 3) = expected_covariance(2, 2);

  REQUIRE(fusion.ekf_.x_.isApprox(expected_state, 1e-10));
  REQUIRE(fusion.ekf_.P_.isApprox(expected_covariance, 1e-10));
}

TEST_CASE("Equal timestamps allow correction without motion", "[fusion_ekf]") {
  FusionEKF fusion;
  fusion.ProcessMeasurement(MakeLidar(1.0, 2.0, 1000000));
  fusion.ProcessMeasurement(MakeLidar(2.0, 0.0, 1000000));

  const double gain = 1.0 / 1.0225;
  const Eigen::Vector4d expected_state(1.0 + gain, 2.0 - 2.0 * gain, 0.0, 0.0);

  Eigen::Matrix4d expected_covariance = Eigen::Matrix4d::Zero();
  expected_covariance.diagonal() << 0.0225 / 1.0225, 0.0225 / 1.0225, 1000.0, 1000.0;

  REQUIRE(fusion.ekf_.x_.isApprox(expected_state, 1e-12));
  REQUIRE(fusion.ekf_.P_.isApprox(expected_covariance, 1e-12));
}

TEST_CASE("Rejected measurements preserve the accepted estimate and timestamp", "[fusion_ekf]") {
  FusionEKF fusion;
  FusionEKF control;

  const auto first = MakeLidar(1.0, 2.0, 1000000);
  fusion.ProcessMeasurement(first);
  control.ProcessMeasurement(first);

  SECTION("Earlier timestamp") {
    REQUIRE_THROWS_AS(fusion.ProcessMeasurement(MakeLidar(3.0, 4.0, 500000)),
                      std::invalid_argument);
  }

  SECTION("Incorrect lidar dimensions") {
    auto invalid = MakeLidar(3.0, 4.0, 1500000);
    invalid.raw_measurements_ = Eigen::Vector3d(3.0, 4.0, 5.0);
    REQUIRE_THROWS_AS(fusion.ProcessMeasurement(invalid), std::invalid_argument);
  }

  SECTION("Negative radar range") {
    MeasurementPackage radar;
    radar.sensor_type_ = MeasurementPackage::RADAR;
    radar.timestamp_ = 1500000;
    radar.raw_measurements_ = Eigen::Vector3d(-2.0, 0.5, 1.0);

    REQUIRE_THROWS_AS(fusion.ProcessMeasurement(radar), std::invalid_argument);
  }

  SECTION("Correction overflows after prediction") {
    const double largest = std::numeric_limits<double>::max();

    fusion.ekf_.x_(0) = -largest;
    control.ekf_.x_(0) = -largest;

    REQUIRE_THROWS_AS(fusion.ProcessMeasurement(MakeLidar(largest, 2.0, 1500000)),
                      std::runtime_error);
  }

  // Rejection must preserve every state component exactly.
  REQUIRE((fusion.ekf_.x_.array() == control.ekf_.x_.array()).all());
  REQUIRE(fusion.ekf_.P_.isApprox(control.ekf_.P_, 1e-12));

  // Restore ordinary values in both filters for the recovery comparison.
  // The internal timestamps remain untouched.
  fusion.ekf_.x_ = Eigen::Vector4d(1.0, 2.0, 0.0, 0.0);
  control.ekf_.x_ = Eigen::Vector4d(1.0, 2.0, 0.0, 0.0);

  const auto next = MakeLidar(2.0, 0.0, 2000000);
  fusion.ProcessMeasurement(next);
  control.ProcessMeasurement(next);

  REQUIRE(fusion.ekf_.x_.isApprox(control.ekf_.x_, 1e-12));
  REQUIRE(fusion.ekf_.P_.isApprox(control.ekf_.P_, 1e-12));
}

TEST_CASE("Radar correction uses radar measurement noise", "[fusion_ekf]") {
  FusionEKF fusion;
  fusion.ProcessMeasurement(MakeLidar(2.0, 0.0, 1000000));

  MeasurementPackage radar;
  radar.sensor_type_ = MeasurementPackage::RADAR;
  radar.timestamp_ = 1000000;
  radar.raw_measurements_ = Eigen::Vector3d(2.0, 0.0, 0.0);

  fusion.ProcessMeasurement(radar);

  const Eigen::Vector4d expected_state(2.0, 0.0, 0.0, 0.0);

  Eigen::Matrix4d expected_covariance = Eigen::Matrix4d::Zero();
  expected_covariance.diagonal() << 0.09 / 1.09, 0.0009 / 0.2509, 1000.0 * 0.09 / 1000.09, 1000.0;

  REQUIRE(fusion.ekf_.x_.isApprox(expected_state, 1e-12));
  REQUIRE(fusion.ekf_.P_.isApprox(expected_covariance, 1e-12));
}

TEST_CASE("Radar update uses the predicted moving state", "[fusion_ekf]") {
  FusionEKF fusion;
  fusion.ProcessMeasurement(MakeLidar(2.0, 0.0, 1000000));

  // Supply a known velocity to isolate prediction and radar correction.
  fusion.ekf_.x_(3) = 1.0;

  // After one second: [px, py, vx, vy] = [2, 1, 0, 1].
  const double range = std::sqrt(5.0);

  MeasurementPackage radar;
  radar.sensor_type_ = MeasurementPackage::RADAR;
  radar.timestamp_ = 2000000;
  radar.raw_measurements_ = Eigen::Vector3d(range, std::atan2(1.0, 2.0), 1.0 / range);

  fusion.ProcessMeasurement(radar);

  const Eigen::Vector4d expected_state(2.0, 1.0, 0.0, 1.0);

  Eigen::Matrix<double, 3, 4> expected_jacobian;
  expected_jacobian << 2.0 / range, 1.0 / range, 0.0, 0.0, -0.2, 0.4, 0.0, 0.0,
      -2.0 / (5.0 * range), 4.0 / (5.0 * range), 2.0 / range, 1.0 / range;

  REQUIRE(fusion.ekf_.x_.isApprox(expected_state, 1e-12));
  REQUIRE(fusion.ekf_.H_.isApprox(expected_jacobian, 1e-12));
  REQUIRE(fusion.ekf_.P_.allFinite());
  REQUIRE(fusion.ekf_.P_.isApprox(fusion.ekf_.P_.transpose(), 1e-12));
}

TEST_CASE("Singular radar correction preserves the filter and allows recovery", "[fusion_ekf]") {
  FusionEKF fusion;
  FusionEKF control;

  const auto first = MakeLidar(0.0, 0.0, 1000000);
  fusion.ProcessMeasurement(first);
  control.ProcessMeasurement(first);

  MeasurementPackage radar;
  radar.sensor_type_ = MeasurementPackage::RADAR;
  radar.timestamp_ = 1500000;
  radar.raw_measurements_ = Eigen::Vector3d(1.0, 0.0, 0.0);

  REQUIRE_THROWS_AS(fusion.ProcessMeasurement(radar), std::domain_error);

  REQUIRE((fusion.ekf_.x_.array() == control.ekf_.x_.array()).all());
  REQUIRE((fusion.ekf_.P_.array() == control.ekf_.P_.array()).all());

  const auto next = MakeLidar(1.0, 0.5, 2000000);
  fusion.ProcessMeasurement(next);
  control.ProcessMeasurement(next);

  REQUIRE(fusion.ekf_.x_.isApprox(control.ekf_.x_, 1e-12));
  REQUIRE(fusion.ekf_.P_.isApprox(control.ekf_.P_, 1e-12));
}

TEST_CASE("Radar update corrects range and radial velocity", "[kalman_filter]") {
  KalmanFilter filter;
  Tools tools;

  filter.x_ = Eigen::Vector4d(2.0, 0.0, 0.0, 0.0);
  filter.P_ = Eigen::Matrix4d::Identity();
  filter.H_ = tools.CalculateJacobian(filter.x_);
  filter.R_ = Eigen::Matrix3d::Identity();

  // Predicted measurement: [2, 0, 0].
  // Innovation: [2, 0, 2].
  filter.UpdateEKF(Eigen::Vector3d(4.0, 0.0, 2.0));

  const Eigen::Vector4d expected_state(3.0, 0.0, 1.0, 0.0);

  Eigen::Matrix4d expected_covariance = Eigen::Matrix4d::Zero();
  expected_covariance.diagonal() << 0.5, 0.8, 0.5, 1.0;

  REQUIRE(filter.x_.isApprox(expected_state, 1e-12));
  REQUIRE(filter.P_.isApprox(expected_covariance, 1e-12));
}