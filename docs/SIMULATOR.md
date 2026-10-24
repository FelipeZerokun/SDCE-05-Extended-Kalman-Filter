# Udacity Term 2 simulator

The optional `ekf_simulator` executable connects the existing EKF to the
Udacity **Term 2, Project 1/2 EKF and UKF** simulator. Both programs must run
on the Ubuntu machine. Windows can build and test the server without running Unity.

## Build on Ubuntu

From the project root, with the simulator integration checked out:

```bash
cmake --preset ubuntu-local -DBUILD_SIMULATOR=ON
cmake --build --preset ubuntu-local
ctest --preset ubuntu-local
```

The first configure installs Boost.Beast and nlohmann-json through the optional
vcpkg `simulator` feature. It may take several minutes. If local user presets
are unavailable, set `VCPKG_ROOT="$HOME/tools/vcpkg"` and substitute the shared
`ubuntu-debug` preset. No system-wide uWebSockets installation is needed.

## Run the actual simulator

1. In a terminal, start the server:

   ```bash
   ./build/ubuntu-debug/ekf_simulator
   ```

2. Wait for `Listening on ws://127.0.0.1:4567`.
3. Open your existing Ubuntu Term 2 simulator executable. Choose the
   **Project 1/2 EKF and UKF** scene if presented, then **Dataset 1 / Start**.
   The exact labels depend on the simulator release.
4. Confirm the terminal prints `Simulator connected; starting a fresh run.`
5. Let the entire dataset finish. The estimate marker should track the target,
   and the displayed RMSE should finish near:

   ```text
   px: 0.097226  py: 0.085376  vx: 0.450855  vy: 0.439588
   ```

   These offline reference values apply to the supplied 500-row Dataset 1.
   RMSE appears in the simulator window; the console reports connections and
   errors. Actual simulator values can differ from offline replay; the verified
   results are recorded below. The benchmark is [0.11, 0.11, 0.52, 0.52].
6. Restart Dataset 1 and check that the run begins fresh. Disconnect/reconnect
   also starts a new EKF and RMSE history. Stop the server with Ctrl+C.

Ubuntu Unity validation is complete, including Dataset 1 and a reset/rerun.
The displayed results from both runs are recorded below.

## Protocol and lifecycle

The adapter preserves the original starter's protocol:

- Receives `42["telemetry",{"sensor_measurement":"L ..."}]` (or a radar row).
- Sends `42["estimate_marker",{...}]` with `estimate_x`, `estimate_y`,
  `rmse_x`, `rmse_y`, `rmse_vx`, and `rmse_vy`.
- Accepts full 10/11-field dataset rows and the simulator's 8/9-field rows
  without unused yaw/yaw-rate fields. The offline parser's contract is unchanged.
- Evaluates all accepted measurements, including initialization.
- Creates a fresh filter and evaluation history for every connection.
- A timestamp moving backward is interpreted as a dataset restart. This is an
  adapter-specific policy: the offline runner still rejects out-of-order data.
- An explicit `42["reset",{}]` event also resets the session.
- Null telemetry receives `42["manual",{}]`. Malformed or invalid telemetry is
  logged and receives the same response without updating the accepted run.
- Handles default-namespace `40` and Engine.IO text ping `2` / pong `3` messages.

This is a small adapter for the starter's direct WebSocket protocol, not a
general Socket.IO server (HTTP polling and arbitrary namespaces are unsupported).
It serves one local simulator connection at a time, with a 64 KiB message limit.
The listener binds only to IPv4 loopback. An alternative port can be passed as
`ekf_simulator 4568`, but the stock simulator expects 4567.

## Troubleshooting

- **Connection refused:** start the server first and run it on the same Ubuntu
  machine as the simulator. A server on Windows cannot serve Ubuntu's localhost.
- **Address already in use:** stop the previous server on port 4567. Do not run
  the legacy and modern servers simultaneously.
- **Rejected simulator message:** inspect the terminal diagnostic. If a release
  uses a different payload, record that error before changing the parser.
- **Restart does not reset the display:** reconnect the simulator, or restart
  the server and reconnect. Same-timestamp measurements alone do not signal restart.
- **Unknown preset:** keep the ignored Ubuntu user presets from the offline setup,
  or use the shared preset with VCPKG_ROOT set as above.

## Development and verification

Windows build commands:

```powershell
cmake --preset windows-local -DBUILD_SIMULATOR=ON
cmake --build --preset windows-local
ctest --preset windows-local
```

The simulator adds protocol tests and a real loopback WebSocket exchange test;
these do not require the Unity executable. CI enables the simulator option on
Windows and Ubuntu. The offline build remains available by configuring with
`-DBUILD_SIMULATOR=OFF` (the default for a new build directory).

All 46 CTest cases pass with the simulator enabled on Windows and Ubuntu CI.
Formatting and Ubuntu static analysis also pass. An
independent .NET WebSocket client also sent all 500 reference dataset rows to
the actual server executable and received the same passing RMSE as offline replay:
`[0.097226, 0.085376, 0.450855, 0.439588]`.

The actual Ubuntu Unity simulator connected successfully and completed Dataset 1
twice, with a reset between runs. The simulator displayed:

| Component | First run | After reset and rerun | Benchmark maximum |
| --- | ---: | ---: | ---: |
| X (m) | 0.0973 | 0.0974 | 0.11 |
| Y (m) | 0.0855 | 0.0855 | 0.11 |
| VX (m/s) | 0.4513 | 0.4517 | 0.52 |
| VY (m/s) | 0.4399 | 0.4404 | 0.52 |

Both runs pass every threshold. These manual observations establish successful
simulator operation before and after reset, not exact numerical repeatability.
The small differences from offline replay and between runs have not been
diagnosed; they are not attributed solely to display rounding.

The server uses Boost.Beast instead of the obsolete uWebSockets dependency.
`src/main.cpp` remains as the original protocol reference; `apps/ekf_simulator.cpp`
is the maintained simulator entry point. No EKF equations are changed.
