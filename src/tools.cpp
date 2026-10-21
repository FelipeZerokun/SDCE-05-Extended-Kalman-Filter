#include "tools.h"

#include <cstddef>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  if (estimations.empty() || ground_truth.empty()) {
    throw std::invalid_argument("RMSE requires nonempty input sequences");
  }

  if (estimations.size() != ground_truth.size()) {
    throw std::invalid_argument("RMSE requires equal estimation and ground-truth sequence lengths");
  }

  for (std::size_t i = 0; i < estimations.size(); ++i) {
    if (estimations[i].size() != 4 || ground_truth[i].size() != 4) {
      throw std::invalid_argument("RMSE requires four components in every state");
    }

    if (!estimations[i].allFinite() || !ground_truth[i].allFinite()) {
      throw std::invalid_argument("RMSE requires finite input components");
    }
  }

  VectorXd squared_error_sum = VectorXd::Zero(4);

  for (std::size_t i = 0; i < estimations.size(); ++i) {
    const VectorXd residual = estimations[i] - ground_truth[i];

    if (!residual.allFinite()) {
      throw std::overflow_error("RMSE residual overflow: non-finite result encountered");
    }

    const VectorXd squared_error = residual.array().square().matrix();

    if (!squared_error.allFinite()) {
      throw std::overflow_error("RMSE squared error overflow: non-finite result encountered");
    }

    squared_error_sum += squared_error;

    if (!squared_error_sum.allFinite()) {
      throw std::overflow_error("RMSE accumulation overflow: non-finite result encountered");
    }
  }

  const VectorXd mean_squared_error = squared_error_sum / static_cast<double>(estimations.size());

  const VectorXd rmse = mean_squared_error.array().sqrt().matrix();

  if (!rmse.allFinite()) {
    throw std::overflow_error("RMSE final calculation overflow: non-finite result encountered");
  }

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  throw std::logic_error("CalculateJacobian is not implemented yet");
}
