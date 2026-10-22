# Extended Kalman Filter — C++

An incremental modernization of the Extended Kalman Filter project from
Udacity's Self-Driving Car Engineer Nanodegree.

The project will estimate a moving object's two-dimensional position and
velocity by combining asynchronous lidar and radar measurements. Its focus is
a reusable C++ estimator, reproducible builds, automated verification, and
quantitative evaluation, with an optional connection to the Udacity simulator.

**Current stage: Milestone 2 complete — RMSE tested on Windows and Ubuntu.**
The project began with nine unchanged Udacity starter files. RMSE is now
implemented in `src/tools.cpp`, with its contract documented in `src/tools.h`.
The starter structure is retained; the legacy solution remains a reference.

The `ekf_core` library and `ekf_tests` executable configure and build in Debug
mode on Windows and Ubuntu. All nine RMSE test cases pass through CTest on both
platforms; the recorded Ubuntu run completed in 0.09 seconds.
Full filtering, offline replay, simulator integration, and CI remain planned.

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
| C++ dependencies | vcpkg manifest mode with a recorded registry baseline; Eigen |
| C++ tests | Catch2 v3-compatible API with CTest; ekf_tests target implemented |
| Formatting and analysis | clang-format, clang-tidy, and compiler warnings |
| Runtime checks | Sanitizer-enabled checks where the selected toolchain supports them |
| Continuous integration | GitHub Actions builds and tests on Ubuntu and Windows |
| Optional Python analysis | uv for dependencies; pytest, Ruff, and mypy for maintained Python code |

C++17 is required, with compiler extensions disabled; CMake 3.22 is the minimum.
The Ubuntu development environment uses Ubuntu 22.04.5 LTS, GCC/G++ 11.4.0, CMake 3.22.1,
clang-tidy/LLVM 14, and uv 0.13.3. Project formatting uses clang-format 23.1.3 via uv. GCC 12.3.0 is also installed
but is not selected. Python tooling will be introduced when analysis requires
it; no Python project is needed for this milestone.

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

**Results:** all nine RMSE test cases pass in Debug mode on Windows and Ubuntu.
The recorded Ubuntu run completed in 0.09 seconds. No full-estimator dataset
evaluation has been performed. Numerical
agreement across platforms will use tolerances rather than a promise of
bit-for-bit identical floating-point output.

## Milestones and acceptance criteria

| Milestone | Scope | Status |
| --- | --- | --- |
| 0. Project definition | Scope, architecture, boundaries, README | Complete |
| 1. Modernization baseline | Copy the Udacity starter sources unchanged into src; verify the starting point | Complete: all nine copied files verified byte-for-byte against the local starter on Windows |
| 2. First tested component | Review references; migrate and test RMSE; establish library and test targets | Complete: documented RMSE contract, library/test targets, nine passing tests on Windows and Ubuntu |
| 3. Reproducible development workflow | Shared presets, formatting, static analysis, Ubuntu/Windows CI | Shared formatting configuration and clang-format 23.1.3 verified on both platforms; configure/build/test presets verified on both platforms; MSVC and GCC warning builds validated; static analysis and CI pending |
| 4. Measurements and dataset input | Types, parser, validation, provenance | Planned |
| 5. Core estimation | Initialization, prediction, lidar update, radar update | Planned |
| 6. Offline tracking application | Replay, configuration, artifacts, regression criteria | Planned |
| 7. Evaluation and analysis | Sensor-mode comparisons, plots, error analysis | Planned |
| 8. Simulator integration and release | Adapter, restart behavior, documentation and final verification | Planned |

Ubuntu environment setup is a completed prerequisite.
The modernization starts from the Udacity starter, with the completed legacy
implementation retained as a reference for review and selective reuse. The
bundled Eigen directory was not copied; Eigen is supplied through vcpkg.
The original source baseline was verified before edits. RMSE builds and tests
pass on Windows and Ubuntu.

Milestone 2 requires a documented RMSE contract and failure behavior, a reusable
library and test executable that compile and link, and passing CTest-discovered
checks for analytical examples and invalid inputs. These tests must run without
the simulator or full estimator.

The simulator adapter is optional for users of the library and offline
application, but remains a planned project deliverable. Simulator failures must
not block progress on the estimator or its automated validation.

## Current project structure

The modernization branch, `modern/cpp`, contains:

```text
include/ekf/        Empty modern public-header directory
src/               Starter structure with RMSE implementation and documented contract
apps/              Empty modern application directory
tests/             test_tools.cpp: nine RMSE test cases
docs/              Local reference documentation
CMakeLists.txt     C++17, ekf_core library, optional ekf_tests and CTest discovery
vcpkg.json         Eigen and Catch2 dependencies with a registry baseline
legacy/            Original completed implementation (ignored)
udacity_baseline/  Udacity starter templates and reference data (ignored)
simulator/         Local simulator distributions (ignored)
```

The copied starter files retain their original responsibilities:

| Files in `src/` | Starting responsibility |
| --- | --- |
| `main.cpp` | Simulator transport and parsing; invokes filtering and RMSE |
| `FusionEKF.cpp`, `FusionEKF.h` | Filter initialization and prediction/update orchestration |
| `kalman_filter.cpp`, `kalman_filter.h` | State and covariance storage; prediction and update interfaces/stubs |
| `tools.cpp`, `tools.h` | Implemented RMSE; Jacobian placeholder throws std::logic_error |
| `measurement_package.h` | Sensor measurement representation |
| `json.hpp` | Bundled JSON dependency used by the starter simulator entry point |

The starter's bundled Eigen directory and old build artifacts were not copied.
Eigen is provided through the existing vcpkg/CMake configuration.
`main.cpp` still depends on the old uWebSockets API; that dependency is not
configured in the modern build. The current source snapshot is therefore a
starting point for development, not a runnable application. RMSE tests run independently of simulator dependencies.

The outer modern repository excludes the three local reference directories,
`data/`, `build/`, `outputs/`, `CMakeUserPresets.json`, and `.vscode/`. Empty directories
are not tracked by Git. Shared settings are tracked in `CMakePresets.json`; machine-specific overrides belong in ignored `CMakeUserPresets.json`.
Local reference directories are maintained separately on Windows and Ubuntu;
they are not supplied by cloning the modern repository.

## Setup and usage

vcpkg is installed outside the Ubuntu repository at `$HOME/tools/vcpkg`.
The manifest records baseline `469a0504e4a6ebf6173d1c0f51a7fa34c382fd43`.
The existing manifest uses `"version": "0.1.0"`.

Install dependencies and configure the Ubuntu Debug build:

```sh
cmake -S . -B build/ubuntu-debug \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_TOOLCHAIN_FILE="$HOME/tools/vcpkg/scripts/buildsystems/vcpkg.cmake"
```

VS Code CMake Tools can use the same build directory and vcpkg toolchain.

Build and run the tests on Ubuntu:

```sh
cmake --build build/ubuntu-debug --parallel 2
ctest --test-dir build/ubuntu-debug --output-on-failure
```

All nine Ubuntu tests pass, with a recorded total runtime of 0.09 seconds. This validates
the current RMSE targets only; full estimator and application support remain planned.

The following native Windows PowerShell workflow was validated with CMake
4.4.3, Visual Studio Build Tools 2026, MSVC 19.51.36256.0, and Windows SDK
10.0.26100.0. vcpkg resolved Catch2 3.16.0 and Eigen 5.0.1.
Adjust the toolchain path if vcpkg is installed elsewhere.

```powershell
cmake -S . -B build/windows -G "Visual Studio 18 2026" -A x64 -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE/tools/vcpkg/scripts/buildsystems/vcpkg.cmake" -DBUILD_TESTING=ON
cmake --build build/windows --config Debug
ctest --test-dir build/windows -C Debug --output-on-failure
```

If the local vcpkg checkout lacks the pinned baseline commit, fetch its remote
history with `git -C "$env:USERPROFILE/tools/vcpkg" fetch origin`, then retry
configuration without changing the manifest baseline.

### Implemented RMSE contract and test coverage

`Tools::CalculateRMSE` accepts constant references to paired sequences of
four-component dynamic Eigen vectors in `[px, py, vx, vy]` order. The caller
aligns observations. It returns component-wise RMSE in metres and metres/second.
Empty inputs, unequal sequence lengths, incorrect state dimensions, and
non-finite components throw `std::invalid_argument`. Non-finite arithmetic
intermediates or results throw `std::overflow_error`.

Direct squared-error accumulation can overflow even when the mathematical RMSE
is representable; tiny squares can underflow to zero. Fixed-size vectors were
considered, but the starter interface is retained with explicit dimension checks.

Nine test cases cover identical sequences, one observation, opposite-signed
errors, a hand-calculated multi-observation example, empty sequences, unequal
lengths, incorrect state dimensions, NaN and both infinities in either input,
and subtraction/squaring/accumulation overflow from finite inputs. Analytical
comparisons use an absolute tolerance of `1e-12`, except exact zero for identical
inputs. This validates RMSE only, not the unfinished filtering implementation.

The Linux simulator launches successfully on Ubuntu.
The Windows simulator crashes during startup; further troubleshooting is outside
this milestone. Simulator integration will be verified on Ubuntu separately from
native Windows core, offline application, and test support.

### CMake presets

`CMakePresets.json` provides shared `windows-debug` and `ubuntu-debug`
configure, build, and test presets. Ubuntu uses Ninja (required on PATH);
Windows uses the Visual Studio generator. Format version 3 supports CMake 3.22;
the Windows preset additionally requires a CMake version supporting the
Visual Studio 18 2026 generator and that toolchain installed.

Create an ignored `CMakeUserPresets.json` in the project root on each machine.
For Windows:

```json
{
    "version": 3,
    "configurePresets": [
        {
            "name": "windows-local",
            "inherits": "windows-debug",
            "environment": {
                "VCPKG_ROOT": "C:/Users/epile/tools/vcpkg"
            }
        }
    ],
    "buildPresets": [
        {
            "name": "windows-local",
            "inherits": "windows-debug",
            "configurePreset": "windows-local"
        }
    ],
    "testPresets": [
        {
            "name": "windows-local",
            "inherits": "windows-debug",
            "configurePreset": "windows-local"
        }
    ]
}
```

Adjust `VCPKG_ROOT` to the local installation. On Ubuntu, replace every
`windows-local` with `ubuntu-local`, every `windows-debug` with `ubuntu-debug`,
and set `VCPKG_ROOT` to `$penv{HOME}/tools/vcpkg`. Keep this local file out of Git.
Local build and test presets inherit the local configure environment, including
vcpkg's path, so automatic CMake regeneration uses the same toolchain settings.

Run from the project root, continuing only after each command succeeds:

```powershell
cmake --preset windows-local
cmake --build --preset windows-local
ctest --preset windows-local
```

This Windows Debug preset workflow passes all nine RMSE tests. The equivalent
Ubuntu Debug preset workflow also configures successfully and passes all nine tests:

```bash
cmake --preset ubuntu-local
cmake --build --preset ubuntu-local
ctest --preset ubuntu-local
```

Presets reuse `build/windows` and `build/ubuntu-debug`. The Ubuntu directory
uses Ninja; a build reporting `no work to do` is already up to date. An existing
build directory must use the same generator as its preset. Test presets display
failure output and fail when no tests are discovered. Shared presets can also
be used directly with `VCPKG_ROOT` supplied by the calling environment.

## Compiler warnings

Warning options apply privately to `ekf_core` and, when `BUILD_TESTING` is
enabled, `ekf_tests`: MSVC uses `/W4`; GCC uses `-Wall -Wextra -Wpedantic`.
Dependency targets are not modified, and warning options are not propagated
to consumers. Warnings are not currently treated as errors.

Windows (MSVC) and Ubuntu (GCC) Debug builds complete without warnings, and
all nine RMSE tests pass on both platforms with these warning options enabled.
The unimplemented Jacobian definition deliberately omits its parameter name
until the implementation needs it; its interface and throwing behavior remain
unchanged.

## Source formatting

The shared `.clang-format` uses LLVM style, two-space indentation, a 100-column
wrapping target, spaces instead of tabs, and regrouped includes. Formatting was
applied to `src/tools.cpp`, `src/tools.h`, and `tests/test_tools.cpp`. All nine
Windows Debug tests passed after rebuilding the formatted sources. Legacy, starter-reference, and bundled third-party files were not
formatted.

The project formatter is **clang-format 23.1.3** on Windows and Ubuntu.
Install the exact version as an isolated tool using uv (no Python project is
required). The community-maintained `clang-format` PyPI distribution packages
LLVM's formatter:

```text
uv tool install "clang-format==23.1.3"
uv tool dir --bin
```

Windows PowerShell, from the project root:

```powershell
$toolBin = uv tool dir --bin
& "$toolBin/clang-format.exe" --version
& "$toolBin/clang-format.exe" --style=file --dry-run --Werror src/tools.cpp src/tools.h tests/test_tools.cpp
```

Ubuntu, from the project root:

```bash
tool_bin="$(uv tool dir --bin)"
"$tool_bin/clang-format" --version
"$tool_bin/clang-format" --style=file --dry-run --Werror src/tools.cpp src/tools.h tests/test_tools.cpp
```

Both platforms pass the formatting check with version 23.1.3. Replace
`--dry-run --Werror` with `-i` to apply formatting to the listed files, then
review the diff. Keep the file list scoped to maintained project sources;
do not format bundled dependencies or local reference copies.

In each machine's ignored `.vscode/settings.json`, set `C_Cpp.formatting` to
`clangFormat` and `C_Cpp.clang_format_path` to the absolute executable path
inside the directory printed by `uv tool dir --bin` (`clang-format.exe` on
Windows, `clang-format` on Ubuntu). This selects the pinned formatter instead
of the extension-bundled or system version. Local paths may differ.

CI enforcement remains planned. `.gitattributes` normalizes C++ source files
to LF line endings.

## Scope and working approach

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
