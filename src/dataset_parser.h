#ifndef DATASET_PARSER_H_
#define DATASET_PARSER_H_

#include "dataset_record.h"

#include <string>

// Parses one complete dataset row:
//
// L px py timestamp gt_px gt_py gt_vx gt_vy gt_yaw gt_yawrate
// R rho phi rhodot timestamp gt_px gt_py gt_vx gt_vy gt_yaw gt_yawrate
//
// Timestamp: nonnegative integer microseconds fitting in long long.
// All floating-point fields must be finite; radar range must be nonnegative.
// Yaw and yaw rate are validated but not stored.
// Leading/trailing whitespace and whitespace-separated fields are accepted.
//
// Throws std::invalid_argument for empty or malformed rows, unsupported
// sensors, missing/extra fields, or values violating these requirements.
//
// Does not check timestamp ordering across records.
DatasetRecord ParseDatasetLine(const std::string &line);

#endif // DATASET_PARSER_H_