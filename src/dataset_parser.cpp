#include "dataset_parser.h"

#include <cmath>
#include <cstddef>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double ParseFiniteDouble(const std::string &token) {
  std::istringstream input(token);
  input.imbue(std::locale::classic());

  double value = 0.0;

  if (!(input >> value) || !input.eof() || !std::isfinite(value)) {
    throw std::invalid_argument("Dataset field must be a finite number: " + token);
  }

  return value;
}

long long ParseTimestamp(const std::string &token) {
  std::istringstream input(token);
  input.imbue(std::locale::classic());

  long long timestamp = 0;

  if (!(input >> timestamp) || !input.eof() || timestamp < 0) {
    throw std::invalid_argument("Dataset timestamp must be a nonnegative integer: " + token);
  }

  return timestamp;
}

} // namespace

DatasetRecord ParseDatasetLine(const std::string &line) {
  std::istringstream input(line);
  input.imbue(std::locale::classic());

  std::vector<std::string> fields;
  std::string field;

  while (input >> field) {
    fields.push_back(field);
  }

  if (fields.empty()) {
    throw std::invalid_argument("Dataset row must not be empty");
  }

  const bool is_lidar = fields[0] == "L";
  const bool is_radar = fields[0] == "R";

  if (!is_lidar && !is_radar) {
    throw std::invalid_argument("Dataset row has an unsupported sensor type");
  }

  const std::size_t expected_fields = is_lidar ? 10 : 11;

  if (fields.size() != expected_fields) {
    throw std::invalid_argument("Dataset row has an incorrect number of fields");
  }

  DatasetRecord record;

  if (is_lidar) {
    record.measurement.sensor_type_ = MeasurementPackage::LASER;
  } else {
    record.measurement.sensor_type_ = MeasurementPackage::RADAR;
  }

  const std::size_t timestamp_index = is_lidar ? 3 : 4;
  record.measurement.timestamp_ = ParseTimestamp(fields[timestamp_index]);

  if (is_lidar) {
    record.measurement.raw_measurements_.resize(2);
    record.measurement.raw_measurements_(0) = ParseFiniteDouble(fields[1]);
    record.measurement.raw_measurements_(1) = ParseFiniteDouble(fields[2]);
  } else {
    record.measurement.raw_measurements_.resize(3);
    record.measurement.raw_measurements_(0) = ParseFiniteDouble(fields[1]);
    record.measurement.raw_measurements_(1) = ParseFiniteDouble(fields[2]);
    record.measurement.raw_measurements_(2) = ParseFiniteDouble(fields[3]);

    if (record.measurement.raw_measurements_(0) < 0.0) {
      throw std::invalid_argument("Dataset radar range must be nonnegative");
    }
  }

  const std::size_t ground_truth_index = timestamp_index + 1;
  record.ground_truth(0) = ParseFiniteDouble(fields[ground_truth_index]);
  record.ground_truth(1) = ParseFiniteDouble(fields[ground_truth_index + 1]);
  record.ground_truth(2) = ParseFiniteDouble(fields[ground_truth_index + 2]);
  record.ground_truth(3) = ParseFiniteDouble(fields[ground_truth_index + 3]);

  // Validate yaw and yaw rate even though we do not store them.
  ParseFiniteDouble(fields[ground_truth_index + 4]);
  ParseFiniteDouble(fields[ground_truth_index + 5]);

  return record;
}
