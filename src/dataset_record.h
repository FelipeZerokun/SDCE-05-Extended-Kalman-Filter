#ifndef DATASET_RECORD_H_
#define DATASET_RECORD_H_

#include "measurement_package.h"

#include <Eigen/Core>

// One dataset observation and its corresponding evaluation reference.
// Only measurement is passed to the estimator.
struct DatasetRecord {
  MeasurementPackage measurement;
  Eigen::Vector4d ground_truth;
};

#endif // DATASET_RECORD_H_