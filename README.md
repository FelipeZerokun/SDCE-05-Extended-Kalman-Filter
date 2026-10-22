# Extended Kalman Filter — C++

A C++17 implementation of Udacity's lidar/radar sensor-fusion challenge, based on
the original starter structure, with CMake, Catch2, vcpkg, clang-format,
clang-tidy, and cross-platform CI.

**The core EKF and offline evaluation are implemented.** The supplied dataset
passes the historical RMSE benchmark on Windows. Simulator integration is not
part of the modern build.

## Results

Windows Debug evaluation on 2026-10-23, using all 500 observations
(250 lidar, 250 radar), including the initial estimate:

| Component | RMSE | Benchmark maximum | Result |
| --- | ---: | ---: | --- |
| px (m) | 0.097226 | 0.11 | Pass |
| py (m) | 0.085376 | 0.11 | Pass |
| vx (m/s) | 0.450855 | 0.52 | Pass |
| vy (m/s) | 0.439588 | 0.52 | Pass |

These are historical Dataset 1 benchmark values, not a formal submission
approval. The original rubric page could not be retrieved during evaluation.
The executable compares unrounded results.

The 40 core/parser tests passed on Windows and Ubuntu; Ubuntu clang-tidy
returned 0 with dependency diagnostics suppressed. The runner adds one CTest
integration test for mixed-sensor replay, malformed rows, empty input, and
benchmark failure. All 41 tests pass on Windows. The new runner's Ubuntu
verification and updated GitHub CI run are pending.

## Build and run

Requirements: C++17, CMake 3.22+, and vcpkg. Ubuntu presets use Ninja. Windows
presets require Visual Studio Build Tools 2026 and a compatible CMake version.
Verified development toolchains: MSVC 19.51/CMake 4.4.3 on Windows and
GCC 11/CMake 3.22.1 on Ubuntu 22.04.

The manifest supplies Eigen and Catch2 at vcpkg baseline
469a0504e4a6ebf6173d1c0f51a7fa34c382fd43.

Windows PowerShell:

~~~powershell
$env:VCPKG_ROOT = "$env:USERPROFILE/tools/vcpkg"
cmake --preset windows-debug
cmake --build --preset windows-debug
ctest --preset windows-debug
.\build\windows\Debug\ekf_replay.exe data/obj_pose-laser-radar-synthetic-input.txt --check-rmse
~~~

Ubuntu:

~~~bash
export VCPKG_ROOT="$HOME/tools/vcpkg"
cmake --preset ubuntu-debug
cmake --build --preset ubuntu-debug
ctest --preset ubuntu-debug
./build/ubuntu-debug/ekf_replay data/obj_pose-laser-radar-synthetic-input.txt --check-rmse
~~~

Existing windows-local and ubuntu-local presets remain supported. Their ignored
CMakeUserPresets.json inherits shared presets and supplies local VCPKG_ROOT
values. Local build/test presets must point to the local configure preset.
VS Code settings are also ignored; use ms-vscode.cmake-tools as the C/C++
configuration provider.

If vcpkg cannot find the baseline, fetch history before configuring:

~~~bash
git -C "$VCPKG_ROOT" fetch origin
~~~

Do not change generators in an existing build directory. For the course's Make
workflow, use a separate directory (not independently verified for the runner):

~~~bash
cmake -S . -B build/ubuntu-make -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
cmake --build build/ubuntu-make --parallel 2
ctest --test-dir build/ubuntu-make --output-on-failure
~~~

## Offline runner

~~~text
ekf_replay <dataset.txt> [--check-rmse]
~~~

Rows are processed in order. Only measurements enter FusionEKF; ground truth
is used afterward for RMSE. Every observation contributes one estimate,
including initialization. The runner prints counts and RMSE in [px, py, vx, vy]
order and does not write output files.

- Success returns 0.
- With --check-rmse, values exceeding [0.11, 0.11, 0.52, 0.52] return 2.
  This benchmark is intended for the reference dataset, not arbitrary trajectories.
- Usage, file, parsing, filtering, or evaluation errors return 1.
- Invalid rows or rejected measurements stop replay with a file/line diagnostic.
  Blank rows are rejected; nothing is silently skipped.

## Dataset and provenance

Place obj_pose-laser-radar-synthetic-input.txt from the
[Udacity starter data directory](https://github.com/udacity/CarND-Extended-Kalman-Filter-Project/tree/master/data)
in local data/. The verified copy matches the local
udacity_baseline/CarND-Extended-Kalman-Filter-Project/data/ copy byte-for-byte.
SHA-256:

~~~text
9879be6c2530833197f1a649cd0118deaf4de86cf90f9b07cda5a01493bdf2e5
~~~

Data/reference folders are ignored and not supplied by cloning this repository.
Automated tests use synthetic rows and need no downloaded dataset.

Accepted rows:

~~~text
L px py timestamp gt_px gt_py gt_vx gt_vy gt_yaw gt_yawrate
R rho phi rhodot timestamp gt_px gt_py gt_vx gt_vy gt_yaw gt_yawrate
~~~

Lidar requires 10 fields and radar 11, including the sensor tag. Abbreviated
8/9-field rows are unsupported. Timestamps are nonnegative integer microseconds
fitting in long long. Numeric fields must parse completely as finite doubles;
scientific notation and whitespace separators are accepted using the classic locale.
Position/range use metres, velocities metres/second, and angles radians.
Radar range must be nonnegative; negative bearing/radial velocity are valid.
Yaw and yaw rate are validated but not stored. DatasetRecord keeps ground truth
separate from MeasurementPackage.

## Estimator behavior

The state is [px, py, vx, vy], with a constant-velocity motion model.

- Lidar initializes Cartesian position directly; radar converts range/bearing.
  Initial velocity is zero and covariance is diag(1, 1, 1000, 1000).
- Elapsed time determines F and Q. Acceleration noise variances remain
  noise_ax = noise_ay = 9.
- Measurement covariance is diag(0.0225, 0.0225) for lidar and
  diag(0.09, 0.0009, 0.09) for radar.
- Radar evaluates the Jacobian at the predicted state and normalizes bearing
  innovation to [-pi, pi]. Corrections use an LDLT solve and Joseph covariance form.
- Equal timestamps allow correction without motion. Earlier timestamps are
  rejected; zero is valid.
- Radar correction/Jacobian evaluation rejects predicted ranges at or below
  1e-6 metres. Initialization at the origin is allowed.
- Rejection preserves the accepted estimate and timestamp.

The original public KalmanFilter matrix interface is retained. Direct callers
must supply consistent dimensions and valid covariance matrices.
FusionEKF::ProcessMeasurement is the tested high-level entry point.

RMSE rejects empty/mismatched sequences, incorrect state dimensions, and
non-finite values. Non-finite arithmetic throws. Direct squared-error accumulation
can overflow for extreme values even when a scaled algorithm could compute RMSE.

## Structure

| File or directory | Responsibility |
| --- | --- |
| src/FusionEKF.cpp | Initialization, timestamps, prediction/update orchestration |
| src/kalman_filter.cpp | Prediction and shared correction mathematics |
| src/tools.cpp | RMSE and radar Jacobian |
| src/dataset_parser.cpp, src/dataset_record.h | Input validation and evaluation reference |
| apps/ekf_replay.cpp | Offline runner |
| tests/ | Catch2 tests and CMake runner integration check |
| src/main.cpp, src/json.hpp | Original simulator entry point/dependency; not built |
| legacy/, udacity_baseline/, simulator/ | Ignored local reference material |
| data/, build/, outputs/ | Ignored local/generated material |

The starter's bundled Eigen is not used. The build produces ekf_core, ekf_replay,
and (with BUILD_TESTING enabled) ekf_tests.

## Quality checks

Formatting uses clang-format 23.1.3, LLVM style, two-space indentation, and a
100-column target. Install with:

~~~text
uv tool install "clang-format==23.1.3"
~~~

PowerShell:

~~~powershell
& "$env:USERPROFILE/.local/bin/clang-format.exe" --style=file --dry-run --Werror src/tools.cpp src/tools.h src/dataset_parser.cpp src/dataset_parser.h src/dataset_record.h src/FusionEKF.cpp src/kalman_filter.cpp src/kalman_filter.h apps/ekf_replay.cpp tests/test_tools.cpp tests/test_dataset_parser.cpp tests/test_FusionEKF.cpp tests/test_kalman_filter.cpp
~~~

On Ubuntu use ~/.local/bin/clang-format with the same arguments. Replace
--dry-run --Werror with -i to format maintained files. Do not format bundled
JSON or reference copies.

Compiler warnings are /W4 for MSVC and -Wall -Wextra -Wpedantic for GCC,
applied to project targets. Ubuntu static analysis uses clang-tidy 14:

~~~bash
clang-tidy -p build/ubuntu-debug --warnings-as-errors='*' \
  src/tools.cpp src/dataset_parser.cpp src/FusionEKF.cpp src/kalman_filter.cpp \
  apps/ekf_replay.cpp tests/test_tools.cpp tests/test_dataset_parser.cpp \
  tests/test_FusionEKF.cpp tests/test_kalman_filter.cpp
~~~

On the Ubuntu development machine, Clang 14 selects GCC 12 standard-library
headers while the build uses GCC 11. Installing libstdc++-12-dev resolved
the earlier missing cmath diagnostic. Suppressed dependency warnings do not
mean dependencies are warning-free.

GitHub Actions checks formatting and builds/tests on Ubuntu 22.04 and Windows
Server 2025 with Visual Studio 2026. Ubuntu additionally runs static analysis.
vcpkg is pinned; hosted images can change. The runner is included in the updated
workflow, whose remote result has not been verified. Commits use Conventional Commits.

## Completion and scope

The core challenge passes the reference dataset's historical accuracy benchmark
in offline evaluation. Remaining verification is the new runner on Ubuntu and
updated GitHub CI.

The original simulator entry point depends on the old uWebSockets API and is
not configured by the modern build. No simulator end-to-end result is claimed.
If submission requires simulator operation, that adapter needs a separate step.

Reset APIs, plots, sensor-mode experiments, configurable tuning, and expanded
frameworks are deferred. They are not prerequisites for the offline result.

| Milestone | Status |
| --- | --- |
| 1. Starter baseline | Complete |
| 2. RMSE and tests | Complete |
| 3. Reproducible workflow | Established on Windows and Ubuntu |
| 4. Dataset parser | Complete |
| 5. Core EKF | Implemented; 40 component tests passed on both platforms |
| 6. Offline evaluation | Implemented; Windows benchmark passed; Ubuntu/CI pending |

## Attribution

Based on the
[Udacity starter](https://github.com/udacity/CarND-Extended-Kalman-Filter-Project)
and completed coursework in local legacy/.
The [original simulator](https://github.com/udacity/self-driving-car-sim/releases)
is maintained separately. Applicable upstream notices must be preserved.
No new license or permission to redistribute course reference material is asserted.
