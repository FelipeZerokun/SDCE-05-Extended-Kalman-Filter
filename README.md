# Extended Kalman Filter — C++

A modern engineering rebuild of the Extended Kalman Filter project from
Udacity's Self-Driving Car Engineer Nanodegree.

The project will estimate a moving object's two-dimensional position and
velocity by combining asynchronous lidar and radar measurements. Its focus is
a reusable C++ estimator, reproducible builds, automated verification, and
quantitative evaluation, with an optional connection to the Udacity simulator.

**Current stage: project definition.** The original implementation was completed
successfully as coursework. The modern implementation, build configuration,
tests, and CI have not been created yet. Features and acceptance criteria below
describe planned work, not verified capabilities.

**Development platforms:** Ubuntu is the primary development environment.
The core library, offline application, and automated tests are intended to
build and run on both Ubuntu and Windows. Platform support will be claimed only
after verification. Simulator compatibility is tracked separately.

## Goal

Turn the original coursework into a maintainable sensor-fusion project that
another developer can build, test, understand, and evaluate.

The modernization will retain the original problem and reference data while
improving architecture, numerical reliability, dependency management, testing,
and documentation. Repeating the starter-code exercise is not the objective.

## Planned capabilities

- A reusable C++ library for initialization, prediction, lidar updates, radar
  updates, and explicit reset.
- An offline command-line application that replays timestamped measurements
  without requiring Unity or a network connection.
- Validated parsing with clear diagnostics for malformed or unsupported input.
- Exported estimates, component-wise RMSE, and reproducible run metadata.
- Lidar-only, radar-only, and fused evaluation modes.
- Trajectory and error plots for comparing estimates with ground truth.
- A separate adapter for the original simulator's communication protocol.
- Automated mathematical, regression, and integration checks.

## Estimation model

The state is `[px, py, vx, vy]`: Cartesian position in metres and velocity
in metres per second. The baseline uses a constant-velocity motion model with
process noise representing unmodelled acceleration.

| Sensor | Observation | Update |
| --- | --- | --- |
| Lidar | Cartesian position `[px, py]` | Linear Kalman filter |
| Radar | Range, bearing, and radial velocity `[rho, phi, rhodot]` | Extended Kalman filter |

The first accepted measurement initializes the state, covariance, and timestamp.
Subsequent measurements determine the elapsed time, update the transition and
process covariance matrices, trigger prediction, and apply the sensor update.

Radar position is converted from polar to Cartesian coordinates for
initialization. A single radial-velocity measurement does not uniquely determine
both Cartesian velocity components. Initial velocity and covariance will be
documented design choices.

The radar innovation uses the nonlinear measurement prediction `z - h(x)`.
Its bearing residual is normalized across the angular wrap boundary. The
Jacobian is evaluated at the predicted state for the radar update.

The baseline retains the supplied noise settings:

| Setting | Value |
| --- | --- |
| Lidar measurement covariance | `diag(0.0225, 0.0225)` |
| Radar measurement covariance | `diag(0.09, 0.0009, 0.09)` |
| Acceleration noise variances | `noise_ax = 9`, `noise_ay = 9` |

Parameter tuning is not required for the baseline. Any later experiment must
record its changes and compare against the unchanged baseline.

## Data contract

The initial reference dataset is
`obj_pose-laser-radar-synthetic-input.txt`. Each row represents one measurement:

```text
L px_measured py_measured timestamp px_gt py_gt vx_gt vy_gt yaw_gt yawrate_gt
R rho_measured phi_measured rhodot_measured timestamp px_gt py_gt vx_gt vy_gt yaw_gt yawrate_gt
```

Only sensor observations and timestamps enter the estimator. Ground truth is
kept separate and used exclusively for evaluation. Yaw and yaw rate are not
part of this project's estimated state.

Original timestamps will be stored as integers and elapsed time converted to
seconds. Timestamp units, coordinate conventions, dataset provenance, and
actual available row formats will be verified against the reference files
before implementing the parser.

Input validation will cover unknown sensor types, malformed rows, non-finite
values, invalid ranges, and timestamp ordering. Policies for duplicate
timestamps and near-zero radar range must be explicit and tested.

## Architecture

| Component | Responsibility |
| --- | --- |
| Estimation library | State, covariance, motion and measurement models, filtering lifecycle |
| Input adapters | Convert dataset rows or simulator messages into validated measurements |
| Offline application | Configure and replay a run, coordinate evaluation, write artifacts |
| Evaluation | Align estimates with ground truth and calculate metrics |
| Simulator adapter | Handle connections, protocol messages, and run resets |
| Analysis tools | Produce plots and comparisons from exported results |

The estimation library must not depend on Unity, network transport, plotting,
or ground truth. Offline replay and simulator integration will use the same
estimation interfaces.

The simulator is a client; the C++ adapter is its server. Incoming messages carry
`sensor_measurement`. Responses contain `estimate_x`, `estimate_y`, `rmse_x`,
`rmse_y`, `rmse_vx`, and `rmse_vy`. Compatibility includes the actual message
framing, not only these field names.

A new run must reset estimator state, timestamps, and evaluation history
together. Reconnection and simulator restart behavior will be defined and
verified during integration.

## Development tools and platform strategy

| Area | Planned approach |
| --- | --- |
| Build | Target-based CMake, with documented presets and out-of-source builds |
| C++ dependencies | Versioned, reproducible dependency resolution; manager to be selected |
| C++ tests | A dedicated test framework integrated with CTest; framework to be selected |
| Formatting and analysis | clang-format, clang-tidy, and compiler warnings |
| Runtime checks | Sanitizer-enabled checks where the selected toolchain supports them |
| Continuous integration | GitHub Actions builds and tests on Ubuntu and Windows |
| Optional Python analysis | uv for dependencies; pytest, Ruff, and mypy for maintained Python code |

The C++ standard, minimum tool versions, dependency versions, and supported
compiler matrix will be selected during environment setup. Python is supporting
tooling, not a dependency of the core estimator.

Ubuntu is the first development target. Windows must have a documented native
build and test workflow; WSL2 alone will not count as native Windows validation.
Paths and application interfaces must not assume a particular shell or machine.
CMake must retain a Make-compatible build path for the original course workflow,
while allowing appropriate generators on Windows.

Docker is optional. It will be introduced only if it solves a concrete
reproducibility or deployment need. Building and testing the estimator must not
require a graphical desktop or simulator installation.

## Verification and evaluation

Verification will include:

- Small analytical examples for prediction and measurement updates.
- Jacobian checks against numerical derivatives away from singularities.
- Angle-wrap, near-zero-range, initialization, and reset behavior.
- Invalid input and timestamp-policy tests.
- Finite state values and covariance symmetry and positive semidefiniteness
  within documented numerical tolerances.
- Dataset replay regression tests with fixed configuration and explicit
  tolerances, including comparisons across supported platforms.
- Offline application integration tests and simulator protocol checks.

Numerical design will prefer solving linear systems to explicitly inverting
matrices and will document the covariance-update strategy.

Evaluation will report separate RMSE values for `px`, `py`, `vx`, and `vy`,
along with sample counts and the evaluation interval. Sensor-mode comparisons
must define a common scoring timeline and initialization policy; comparing
metrics over unrelated subsets would be misleading.

Each recorded experiment should identify the source revision, dataset
fingerprint, configuration, toolchain, sensor mode, and relevant dependency
versions. Results should include estimates, metrics, and trajectory/error plots.
Existing run artifacts should be protected from accidental overwrite.

The legacy implementation is a comparison reference, not an unquestioned
correctness oracle. The original rubric's thresholds may be retained as
historical benchmarks once verified. Modern numerical acceptance thresholds
will be documented before regression tests are finalized.

**Results:** no modern implementation has been evaluated yet. Numerical
agreement across platforms will use tolerances rather than a promise of
bit-for-bit identical floating-point output.

## Milestones and acceptance criteria

| Milestone | Completion criteria |
| --- | --- |
| 1. Repository and toolchain | Repository ownership and branch layout settled; dependencies and tool versions selected; minimal build and test workflow works on Ubuntu and Windows |
| 2. Data and interfaces | Dataset provenance and units verified; measurement/state contracts and validation policies documented; parser tests pass |
| 3. Estimator | Prediction and both updates implemented; analytical, edge-case, covariance, and lifecycle checks pass |
| 4. Offline evaluation | Dataset replay exports reproducible artifacts; numerical regression criteria are defined and met; sensor-mode comparisons are documented |
| 5. Simulator adapter | Measurements and responses integrate with the same estimator; restart behavior verified; working simulator platform and version recorded |
| 6. Release documentation | Clean-checkout instructions verified on both platforms; CI passes; results, design decisions, attribution, and limitations documented |

The simulator adapter is optional for users of the library and offline
application, but remains a planned project deliverable. Simulator failures must
not block progress on the estimator or its automated validation.

## Current project structure

The existing reference folders are:

```text
legacy/             Original completed implementation
udacity_baseline/   Udacity starter templates and reference data
simulator/          Local simulator distributions
README.md           Modern project definition
```

This is the current reference layout, not the final source layout. The modern
library, applications, tests, configuration, and analysis directories will be
introduced as their responsibilities are implemented.

The outer repository and modernization branch arrangement still need to be
verified and finalized. Reference folders must not accidentally become embedded
Git repositories in the modern project's history. Their inclusion or exclusion
policy will be decided explicitly.

## Setup and usage

There are no modern build or execution commands yet. Verified Ubuntu and Windows
instructions will be added as the toolchain and applications become available.
Legacy commands describe the historical implementation only.

The Windows simulator currently closes during startup on the development PC.
Testing its Linux build on the Ubuntu laptop is planned; compatibility has not
yet been confirmed. Simulator binaries are local resources and are not intended
for inclusion in the modern repository.

## Scope and working approach

- The project owner writes implementation code and manages project files and
  folders. The assistant guides design, reviews code, and maintains authorized
  documentation.
- Work proceeds in small, reviewable steps with verification appropriate to each
  change.
- The original implementation and Udacity baseline remain reference material.
  Modernization does not require copying their architecture or bundled build
  artifacts.
- Infrastructure must serve correctness, maintainability, or reproducibility.
- This is an educational tracking system, not a safety-certified automotive
  component. Multi-object tracking, perception from raw sensor signals, and
  vehicle control are outside the initial scope.

## Attribution and licensing

- [Original implementation](https://github.com/FelipeZerokun/SDCE-05-Extended-Kalman-Filter)
- [Udacity starter project](https://github.com/udacity/CarND-Extended-Kalman-Filter-Project)
- [Udacity simulator releases](https://github.com/udacity/self-driving-car-sim/releases)

Applicable upstream notices must be preserved when reusing code or data.
The modern project's license and redistribution permissions for reference
materials will be reviewed before publication; no new license is asserted here.
