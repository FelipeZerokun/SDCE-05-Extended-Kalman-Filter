#include "dataset_parser.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string JoinFields(const std::vector<std::string> &fields) {
  std::string row;
  for (const auto &field : fields) {
    row += field + " ";
  }
  return row;
}
} // namespace

TEST_CASE("Dataset parser validates timestamps for both sensors", "[dataset_parser]") {
  for (const std::string prefix : {"L 1 2 ", "R 1 2 3 "}) {
    CAPTURE(prefix);
    const std::string suffix = " 4 5 6 7 8 9";
    for (const auto timestamp : {0LL, std::numeric_limits<long long>::max()}) {
      CAPTURE(timestamp);
      REQUIRE(
          ParseDatasetLine(prefix + std::to_string(timestamp) + suffix).measurement.timestamp_ ==
          timestamp);
    }
    for (const std::string timestamp :
         {"-1", "1.5", "1e6", "123abc", "invalid", "9223372036854775808"}) {
      CAPTURE(timestamp);
      REQUIRE_THROWS_AS(ParseDatasetLine(prefix + timestamp + suffix), std::invalid_argument);
    }
  }
}

TEST_CASE("Dataset parser validates every floating-point field", "[dataset_parser]") {
  for (const std::string sensor : {"L", "R"}) {
    const bool radar = sensor == "R";
    const std::vector<std::string> valid =
        radar ? std::vector<std::string>{"R", "1", "2", "3", "0", "4", "5", "6", "7", "8", "9"}
              : std::vector<std::string>{"L", "1", "2", "0", "4", "5", "6", "7", "8", "9"};
    for (std::size_t index = 1; index < valid.size(); ++index) {
      if (index == (radar ? 4U : 3U)) {
        continue;
      }
      for (const std::string token : {"invalid", "1.2abc", "nan", "inf", "-inf", "1e9999"}) {
        CAPTURE(sensor, index, token);
        auto fields = valid;
        fields[index] = token;
        REQUIRE_THROWS_AS(ParseDatasetLine(JoinFields(fields)), std::invalid_argument);
      }
    }
  }
}

TEST_CASE("Dataset parser reads a valid lidar row", "[dataset_parser]") {
  const DatasetRecord record =
      ParseDatasetLine("L 1.25 -2.5 1477010443000000 1.0 -2.0 3.5 -4.5 0.2 -0.1");

  REQUIRE(record.measurement.sensor_type_ == MeasurementPackage::LASER);
  REQUIRE(record.measurement.timestamp_ == 1477010443000000LL);
  REQUIRE(record.measurement.raw_measurements_.size() == 2);

  constexpr double tolerance = 1e-12;

  REQUIRE_THAT(record.measurement.raw_measurements_(0),
               Catch::Matchers::WithinAbs(1.25, tolerance));
  REQUIRE_THAT(record.measurement.raw_measurements_(1),
               Catch::Matchers::WithinAbs(-2.5, tolerance));

  REQUIRE_THAT(record.ground_truth(0), Catch::Matchers::WithinAbs(1.0, tolerance));
  REQUIRE_THAT(record.ground_truth(1), Catch::Matchers::WithinAbs(-2.0, tolerance));
  REQUIRE_THAT(record.ground_truth(2), Catch::Matchers::WithinAbs(3.5, tolerance));
  REQUIRE_THAT(record.ground_truth(3), Catch::Matchers::WithinAbs(-4.5, tolerance));
}

TEST_CASE("Dataset parser reads a valid radar row", "[dataset_parser]") {
  const DatasetRecord record =
      ParseDatasetLine("R 5.25 -0.5 -2.75 1477010443050000 1.0 -2.0 3.5 -4.5 0.2 -0.1");

  REQUIRE(record.measurement.sensor_type_ == MeasurementPackage::RADAR);
  REQUIRE(record.measurement.timestamp_ == 1477010443050000LL);
  REQUIRE(record.measurement.raw_measurements_.size() == 3);

  constexpr double tolerance = 1e-12;

  REQUIRE_THAT(record.measurement.raw_measurements_(0),
               Catch::Matchers::WithinAbs(5.25, tolerance));
  REQUIRE_THAT(record.measurement.raw_measurements_(1),
               Catch::Matchers::WithinAbs(-0.5, tolerance));
  REQUIRE_THAT(record.measurement.raw_measurements_(2),
               Catch::Matchers::WithinAbs(-2.75, tolerance));

  REQUIRE_THAT(record.ground_truth(0), Catch::Matchers::WithinAbs(1.0, tolerance));
  REQUIRE_THAT(record.ground_truth(1), Catch::Matchers::WithinAbs(-2.0, tolerance));
  REQUIRE_THAT(record.ground_truth(2), Catch::Matchers::WithinAbs(3.5, tolerance));
  REQUIRE_THAT(record.ground_truth(3), Catch::Matchers::WithinAbs(-4.5, tolerance));
}

TEST_CASE("Dataset parser rejects empty rows", "[dataset_parser]") {
  REQUIRE_THROWS_AS(ParseDatasetLine(""), std::invalid_argument);
  REQUIRE_THROWS_AS(ParseDatasetLine(" \t\r\n"), std::invalid_argument);
}

TEST_CASE("Dataset parser rejects unsupported sensor types", "[dataset_parser]") {
  REQUIRE_THROWS_AS(ParseDatasetLine("X 1.25 -2.5 1477010443000000 1.0 -2.0 3.5 -4.5 0.2 -0.1"),
                    std::invalid_argument);
}

TEST_CASE("Dataset parser rejects incorrect field counts", "[dataset_parser]") {
  SECTION("Lidar row is missing yaw rate") {
    REQUIRE_THROWS_AS(ParseDatasetLine("L 1.25 -2.5 1477010443000000 1.0 -2.0 3.5 -4.5 0.2"),
                      std::invalid_argument);
  }

  SECTION("Lidar row has an extra field") {
    REQUIRE_THROWS_AS(
        ParseDatasetLine("L 1.25 -2.5 1477010443000000 1.0 -2.0 3.5 -4.5 0.2 -0.1 99"),
        std::invalid_argument);
  }

  SECTION("Radar row is missing yaw rate") {
    REQUIRE_THROWS_AS(ParseDatasetLine("R 5.25 -0.5 -2.75 1477010443050000 1.0 -2.0 3.5 -4.5 0.2"),
                      std::invalid_argument);
  }

  SECTION("Radar row has an extra field") {
    REQUIRE_THROWS_AS(
        ParseDatasetLine("R 5.25 -0.5 -2.75 1477010443050000 1.0 -2.0 3.5 -4.5 0.2 -0.1 99"),
        std::invalid_argument);
  }
}

TEST_CASE("Dataset parser rejects invalid timestamps", "[dataset_parser]") {
  const auto make_row = [](const std::string &timestamp) {
    return "L 1.25 -2.5 " + timestamp + " 1.0 -2.0 3.5 -4.5 0.2 -0.1";
  };

  SECTION("Negative timestamp") {
    REQUIRE_THROWS_AS(ParseDatasetLine(make_row("-1")), std::invalid_argument);
  }

  SECTION("Fractional timestamp") {
    REQUIRE_THROWS_AS(ParseDatasetLine(make_row("123.5")), std::invalid_argument);
  }

  SECTION("Scientific notation") {
    REQUIRE_THROWS_AS(ParseDatasetLine(make_row("1e6")), std::invalid_argument);
  }

  SECTION("Trailing characters") {
    REQUIRE_THROWS_AS(ParseDatasetLine(make_row("123abc")), std::invalid_argument);
  }

  SECTION("Integer overflow") {
    REQUIRE_THROWS_AS(ParseDatasetLine(make_row("999999999999999999999999999999")),
                      std::invalid_argument);
  }
}

TEST_CASE("Dataset parser validates radar range", "[dataset_parser]") {
  SECTION("Zero range is accepted") {
    const auto record =
        ParseDatasetLine("R 0 -0.5 -2.75 1477010443050000 1.0 -2.0 3.5 -4.5 0.2 -0.1");

    REQUIRE(record.measurement.raw_measurements_(0) == 0.0);
  }

  SECTION("Negative range is rejected") {
    REQUIRE_THROWS_AS(
        ParseDatasetLine("R -1 -0.5 -2.75 1477010443050000 1.0 -2.0 3.5 -4.5 0.2 -0.1"),
        std::invalid_argument);
  }
}

TEST_CASE("Dataset parser rejects invalid floating-point fields", "[dataset_parser]") {
  SECTION("Measurement contains trailing characters") {
    REQUIRE_THROWS_AS(
        ParseDatasetLine("L 1.25abc -2.5 1477010443000000 1.0 -2.0 3.5 -4.5 0.2 -0.1"),
        std::invalid_argument);
  }

  SECTION("Radar measurement exceeds floating-point range") {
    REQUIRE_THROWS_AS(
        ParseDatasetLine("R 5.25 -0.5 1e9999 1477010443050000 1.0 -2.0 3.5 -4.5 0.2 -0.1"),
        std::invalid_argument);
  }

  SECTION("Ground truth contains NaN") {
    REQUIRE_THROWS_AS(ParseDatasetLine("L 1.25 -2.5 1477010443000000 nan -2.0 3.5 -4.5 0.2 -0.1"),
                      std::invalid_argument);
  }

  SECTION("Yaw contains infinity") {
    REQUIRE_THROWS_AS(ParseDatasetLine("L 1.25 -2.5 1477010443000000 1.0 -2.0 3.5 -4.5 inf -0.1"),
                      std::invalid_argument);
  }

  SECTION("Yaw rate is not numeric") {
    REQUIRE_THROWS_AS(
        ParseDatasetLine("R 5.25 -0.5 -2.75 1477010443050000 1.0 -2.0 3.5 -4.5 0.2 invalid"),
        std::invalid_argument);
  }
}

TEST_CASE("Dataset parser accepts whitespace and scientific notation", "[dataset_parser]") {
  const auto record =
      ParseDatasetLine(" \tL\t1.25e1  -2.5e-1\t0  1e0 -2e0 3.5e0 -4.5e0 2e-1 -1e-1 \r\n");

  REQUIRE(record.measurement.sensor_type_ == MeasurementPackage::LASER);
  REQUIRE(record.measurement.timestamp_ == 0);
  REQUIRE(record.measurement.raw_measurements_.size() == 2);

  constexpr double tolerance = 1e-12;
  REQUIRE_THAT(record.measurement.raw_measurements_(0),
               Catch::Matchers::WithinAbs(12.5, tolerance));
  REQUIRE_THAT(record.measurement.raw_measurements_(1),
               Catch::Matchers::WithinAbs(-0.25, tolerance));

  REQUIRE_THAT(record.ground_truth(0), Catch::Matchers::WithinAbs(1.0, tolerance));
  REQUIRE_THAT(record.ground_truth(1), Catch::Matchers::WithinAbs(-2.0, tolerance));
  REQUIRE_THAT(record.ground_truth(2), Catch::Matchers::WithinAbs(3.5, tolerance));
  REQUIRE_THAT(record.ground_truth(3), Catch::Matchers::WithinAbs(-4.5, tolerance));
}
