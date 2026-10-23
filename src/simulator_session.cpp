#include "simulator_session.h"

#include "dataset_parser.h"
#include "tools.h"

#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {
using nlohmann::json;
constexpr const char *manual = "42[\"manual\",{}]";

DatasetRecord ParseTelemetry(const std::string &row) {
  // Simulator builds may omit the dataset's unused yaw and yaw-rate fields.
  // Keep the offline parser strict; adapt only the simulator's short format.
  std::istringstream input(row);
  std::vector<std::string> fields;
  std::string field;
  while (input >> field) {
    fields.push_back(field);
  }
  if (!fields.empty() &&
      ((fields[0] == "L" && fields.size() == 8) || (fields[0] == "R" && fields.size() == 9))) {
    return ParseDatasetLine(row + " 0 0");
  }
  return ParseDatasetLine(row);
}
} // namespace

std::optional<std::string> SimulatorSession::HandleMessage(const std::string &message) {
  // Engine.IO text ping and Socket.IO default-namespace connection messages.
  if (!message.empty() && message[0] == '2') {
    return "3" + message.substr(1);
  }
  if (message == "40") {
    return "40";
  }
  if (message.compare(0, 2, "42") != 0) {
    return std::nullopt;
  }

  const auto event = json::parse(message.substr(2));
  if (event.is_null()) {
    return manual;
  }
  if (!event.is_array() || event.size() != 2 || !event[0].is_string()) {
    throw std::invalid_argument("Expected a Socket.IO [event, payload] array");
  }
  const auto name = event[0].get<std::string>();
  if (name == "reset") {
    *this = SimulatorSession{};
    return manual;
  }
  if (name != "telemetry") {
    return std::nullopt;
  }
  if (event[1].is_null()) {
    return manual;
  }
  if (!event[1].is_object() || !event[1].contains("sensor_measurement") ||
      !event[1]["sensor_measurement"].is_string()) {
    throw std::invalid_argument("Telemetry requires a sensor_measurement string");
  }
  const auto record = ParseTelemetry(event[1]["sensor_measurement"].get<std::string>());

  // Commit filter and evaluation history together only after successful processing.
  auto candidate = (last_timestamp_ && record.measurement.timestamp_ < *last_timestamp_)
                       ? SimulatorSession{}
                       : *this;
  candidate.filter_.ProcessMeasurement(record.measurement);
  candidate.estimates_.push_back(candidate.filter_.ekf_.x_);
  candidate.ground_truth_.push_back(record.ground_truth);
  Tools tools;
  const auto rmse = tools.CalculateRMSE(candidate.estimates_, candidate.ground_truth_);
  const auto &state = candidate.filter_.ekf_.x_;
  const json payload = {{"estimate_x", state(0)}, {"estimate_y", state(1)}, {"rmse_x", rmse(0)},
                        {"rmse_y", rmse(1)},      {"rmse_vx", rmse(2)},     {"rmse_vy", rmse(3)}};
  const std::string response = "42" + json::array({"estimate_marker", payload}).dump();
  candidate.last_timestamp_ = record.measurement.timestamp_;
  *this = std::move(candidate);
  return response;
}
