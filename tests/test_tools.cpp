#include "tools.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

TEST_CASE("RMSE accumulates opposite-signed errors", "[rmse]") {
  Eigen::VectorXd first(4);
  first << 1.0, 2.0, 3.0, 4.0;

  Eigen::VectorXd second(4);
  second << -1.0, -2.0, -3.0, -4.0;

  const std::vector<Eigen::VectorXd> estimates{first, second};
  const std::vector<Eigen::VectorXd> ground_truth{Eigen::VectorXd::Zero(4),
                                                  Eigen::VectorXd::Zero(4)};

  Tools tools;
  const Eigen::VectorXd result = tools.CalculateRMSE(estimates, ground_truth);

  REQUIRE(result.size() == 4);

  // Absolute tolerance for these small, hand-calculated expected values.
  constexpr double tolerance = 1e-12;

  REQUIRE_THAT(result(0), Catch::Matchers::WithinAbs(1.0, tolerance));
  REQUIRE_THAT(result(1), Catch::Matchers::WithinAbs(2.0, tolerance));
  REQUIRE_THAT(result(2), Catch::Matchers::WithinAbs(3.0, tolerance));
  REQUIRE_THAT(result(3), Catch::Matchers::WithinAbs(4.0, tolerance));
}

TEST_CASE("RMSE is zero for identical sequences", "[rmse]") {
  Eigen::VectorXd first(4);
  first << 1.0, -2.0, 3.0, -4.0;

  Eigen::VectorXd second(4);
  second << -5.0, 6.0, -7.0, 8.0;

  const std::vector<Eigen::VectorXd> estimates{first, second};
  const std::vector<Eigen::VectorXd> ground_truth = estimates;

  Tools tools;
  const Eigen::VectorXd result = tools.CalculateRMSE(estimates, ground_truth);

  REQUIRE(result.size() == 4);

  for (Eigen::Index component = 0; component < result.size(); ++component) {
    REQUIRE(result(component) == 0.0);
  }
}

TEST_CASE("RMSE for one observation equals absolute error", "[rmse]") {
  Eigen::VectorXd estimate(4);
  estimate << 2.0, -3.0, 5.0, -7.0;

  Eigen::VectorXd truth(4);
  truth << 5.0, -5.0, 1.0, -2.0;

  const std::vector<Eigen::VectorXd> estimates{estimate};
  const std::vector<Eigen::VectorXd> ground_truth{truth};

  Tools tools;
  const Eigen::VectorXd result = tools.CalculateRMSE(estimates, ground_truth);

  REQUIRE(result.size() == 4);

  constexpr double tolerance = 1e-12;

  REQUIRE_THAT(result(0), Catch::Matchers::WithinAbs(3.0, tolerance));
  REQUIRE_THAT(result(1), Catch::Matchers::WithinAbs(2.0, tolerance));
  REQUIRE_THAT(result(2), Catch::Matchers::WithinAbs(4.0, tolerance));
  REQUIRE_THAT(result(3), Catch::Matchers::WithinAbs(5.0, tolerance));
}

TEST_CASE("RMSE matches hand-calculated multiple observations", "[rmse]") {
  Eigen::VectorXd first_estimate(4);
  first_estimate << 2.0, 4.0, 6.0, 8.0;

  Eigen::VectorXd second_estimate(4);
  second_estimate << 5.0, 3.0, 13.0, 1.0;

  Eigen::VectorXd first_truth(4);
  first_truth << 1.0, 2.0, 3.0, 4.0;

  Eigen::VectorXd second_truth(4);
  second_truth << -2.0, 1.0, 4.0, -3.0;

  const std::vector<Eigen::VectorXd> estimates{first_estimate, second_estimate};
  const std::vector<Eigen::VectorXd> ground_truth{first_truth, second_truth};

  Tools tools;
  const Eigen::VectorXd result = tools.CalculateRMSE(estimates, ground_truth);

  REQUIRE(result.size() == 4);

  constexpr double tolerance = 1e-12;

  REQUIRE_THAT(result(0), Catch::Matchers::WithinAbs(5.0, tolerance));
  REQUIRE_THAT(result(1), Catch::Matchers::WithinAbs(2.0, tolerance));
  REQUIRE_THAT(result(2), Catch::Matchers::WithinAbs(std::sqrt(45.0), tolerance));
  REQUIRE_THAT(result(3), Catch::Matchers::WithinAbs(4.0, tolerance));
}

TEST_CASE("RMSE rejects empty input sequences", "[rmse]") {
  const std::vector<Eigen::VectorXd> empty;
  const std::vector<Eigen::VectorXd> valid{Eigen::VectorXd::Zero(4)};

  Tools tools;

  REQUIRE_THROWS_AS(tools.CalculateRMSE(empty, empty), std::invalid_argument);
  REQUIRE_THROWS_AS(tools.CalculateRMSE(empty, valid), std::invalid_argument);
  REQUIRE_THROWS_AS(tools.CalculateRMSE(valid, empty), std::invalid_argument);
}

TEST_CASE("RMSE rejects mismatched sequence lengths", "[rmse]") {
  const std::vector<Eigen::VectorXd> one_state{Eigen::VectorXd::Zero(4)};
  const std::vector<Eigen::VectorXd> two_states{Eigen::VectorXd::Zero(4), Eigen::VectorXd::Zero(4)};

  Tools tools;

  REQUIRE_THROWS_AS(tools.CalculateRMSE(one_state, two_states), std::invalid_argument);
  REQUIRE_THROWS_AS(tools.CalculateRMSE(two_states, one_state), std::invalid_argument);
}

TEST_CASE("RMSE rejects states with incorrect dimensions", "[rmse]") {
  const std::vector<Eigen::VectorXd> valid{Eigen::VectorXd::Zero(4), Eigen::VectorXd::Zero(4)};

  Tools tools;

  for (const int dimension : {0, 3, 5}) {
    CAPTURE(dimension);

    const std::vector<Eigen::VectorXd> invalid{Eigen::VectorXd::Zero(4),
                                               Eigen::VectorXd::Zero(dimension)};

    REQUIRE_THROWS_AS(tools.CalculateRMSE(invalid, valid), std::invalid_argument);
    REQUIRE_THROWS_AS(tools.CalculateRMSE(valid, invalid), std::invalid_argument);
    REQUIRE_THROWS_AS(tools.CalculateRMSE(invalid, invalid), std::invalid_argument);
  }
}

TEST_CASE("RMSE rejects non-finite input components", "[rmse]") {
  const std::vector<Eigen::VectorXd> valid{Eigen::VectorXd::Zero(4), Eigen::VectorXd::Zero(4)};

  const double invalid_values[] = {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::infinity(),
                                   -std::numeric_limits<double>::infinity()};

  Tools tools;

  for (const double value : invalid_values) {
    for (Eigen::Index component = 0; component < 4; ++component) {
      CAPTURE(value, component);

      std::vector<Eigen::VectorXd> invalid = valid;
      invalid[1](component) = value;

      REQUIRE_THROWS_AS(tools.CalculateRMSE(invalid, valid), std::invalid_argument);
      REQUIRE_THROWS_AS(tools.CalculateRMSE(valid, invalid), std::invalid_argument);
    }
  }
}

TEST_CASE("RMSE rejects arithmetic overflow", "[rmse]") {
  Tools tools;

  SECTION("Subtracting finite inputs overflows") {
    const double largest = std::numeric_limits<double>::max();

    Eigen::VectorXd estimate = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd truth = Eigen::VectorXd::Zero(4);
    estimate(0) = largest;
    truth(0) = -largest;

    const std::vector<Eigen::VectorXd> estimates{estimate};
    const std::vector<Eigen::VectorXd> ground_truth{truth};

    REQUIRE_THROWS_AS(tools.CalculateRMSE(estimates, ground_truth), std::overflow_error);
  }

  SECTION("Squaring a finite residual overflows") {
    Eigen::VectorXd estimate = Eigen::VectorXd::Zero(4);
    estimate(0) = 1e200;

    const std::vector<Eigen::VectorXd> estimates{estimate};
    const std::vector<Eigen::VectorXd> ground_truth{Eigen::VectorXd::Zero(4)};

    REQUIRE_THROWS_AS(tools.CalculateRMSE(estimates, ground_truth), std::overflow_error);
  }

  SECTION("Accumulating finite squared errors overflows") {
    Eigen::VectorXd estimate = Eigen::VectorXd::Zero(4);
    estimate(0) = 1e154;

    const std::vector<Eigen::VectorXd> estimates{estimate, estimate};
    const std::vector<Eigen::VectorXd> ground_truth{Eigen::VectorXd::Zero(4),
                                                    Eigen::VectorXd::Zero(4)};

    REQUIRE_THROWS_AS(tools.CalculateRMSE(estimates, ground_truth), std::overflow_error);
  }
}