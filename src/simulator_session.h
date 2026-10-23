#ifndef SIMULATOR_SESSION_H_
#define SIMULATOR_SESSION_H_

#include "FusionEKF.h"

#include <optional>
#include <string>
#include <vector>

// One simulator connection/run. Invalid telemetry throws without changing state.
// A backward timestamp starts a new run (the simulator's dataset restart).
class SimulatorSession {
public:
  std::optional<std::string> HandleMessage(const std::string &message);

private:
  FusionEKF filter_;
  std::optional<long long> last_timestamp_;
  std::vector<Eigen::VectorXd> estimates_;
  std::vector<Eigen::VectorXd> ground_truth_;
};

#endif
