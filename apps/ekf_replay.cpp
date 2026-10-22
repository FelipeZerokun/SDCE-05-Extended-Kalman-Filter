#include "FusionEKF.h"
#include "dataset_parser.h"
#include "tools.h"

#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 3 || (argc == 3 && std::string(argv[2]) != "--check-rmse")) {
    std::cerr << "Usage: ekf_replay <dataset.txt> [--check-rmse]\n";
    return 1;
  }

  std::ifstream input(argv[1]);
  if (!input) {
    std::cerr << "Cannot open dataset: " << argv[1] << '\n';
    return 1;
  }

  FusionEKF fusion;
  Tools tools;
  std::vector<Eigen::VectorXd> estimates;
  std::vector<Eigen::VectorXd> ground_truth;
  std::string line;
  std::size_t line_number = 0;
  std::size_t lidar_count = 0;
  std::size_t radar_count = 0;

  while (std::getline(input, line)) {
    ++line_number;
    try {
      const auto record = ParseDatasetLine(line);
      fusion.ProcessMeasurement(record.measurement);
      estimates.push_back(fusion.ekf_.x_);
      ground_truth.push_back(record.ground_truth);
      if (record.measurement.sensor_type_ == MeasurementPackage::LASER) {
        ++lidar_count;
      } else {
        ++radar_count;
      }
    } catch (const std::exception &error) {
      std::cerr << argv[1] << ':' << line_number << ": " << error.what() << '\n';
      return 1;
    }
  }

  if (input.bad() || !input.eof()) {
    std::cerr << "Error reading dataset: " << argv[1] << '\n';
    return 1;
  }
  if (estimates.empty()) {
    std::cerr << "Dataset contains no measurements\n";
    return 1;
  }

  try {
    const Eigen::VectorXd rmse = tools.CalculateRMSE(estimates, ground_truth);
    std::cout << "Measurements: " << estimates.size() << " (lidar: " << lidar_count
              << ", radar: " << radar_count << ")\n"
              << "RMSE [px py vx vy]: " << std::fixed << std::setprecision(6) << rmse.transpose()
              << '\n';

    if (argc == 3) {
      const Eigen::Vector4d limits(0.11, 0.11, 0.52, 0.52);
      const bool passed = (rmse.array() <= limits.array()).all();
      std::cout << "Benchmark [0.11 0.11 0.52 0.52]: " << (passed ? "PASS" : "FAIL") << '\n';
      return passed ? 0 : 2;
    }
  } catch (const std::exception &error) {
    std::cerr << "RMSE evaluation failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
