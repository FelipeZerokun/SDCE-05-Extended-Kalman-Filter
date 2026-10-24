# EKF visualization

The Python renderer creates a video, looping GIF, and static figure from actual
C++ filter output. It starts a temporary `ekf_simulator`, replays the local Udacity
dataset over WebSocket, and closes that server after collecting estimates.
The EKF is not reimplemented in Python. Unity is not required for rendering.

## Generate the media

Use the simulator integration branch until it is merged into main. Build
`ekf_simulator` using the [simulator instructions](SIMULATOR.md). With Python
3.12 or newer, create a separate environment from the repository root:

```bash
python -m venv build/visualization-env
```

Windows:

```powershell
build/visualization-env/Scripts/python.exe -m pip install -r scripts/requirements-visualization.txt
build/visualization-env/Scripts/python.exe scripts/visualize_ekf.py data/obj_pose-laser-radar-synthetic-input.txt --server build/windows/Debug/ekf_simulator.exe
```

Ubuntu:

```bash
build/visualization-env/bin/python -m pip install -r scripts/requirements-visualization.txt
build/visualization-env/bin/python scripts/visualize_ekf.py data/obj_pose-laser-radar-synthetic-input.txt --server build/ubuntu-debug/ekf_simulator
```

The default local port is 14568; use `--port` to choose another unused port.
Use `--output` to change the default `outputs/visualization` destination.
Python dependencies are separate from the CMake/vcpkg build. imageio-ffmpeg
supplies the FFmpeg executable used to encode the MP4.

## Outputs and interpretation

- `ekf_sensor_fusion.mp4`: 1280 × 960 H.264 video, 12 frames/second.
- `ekf_sensor_fusion.gif`: smaller looping preview, 5 frames/second.
- `ekf_summary.png`: 2304 × 1728 final trajectory and RMSE figure.
- `ekf_telemetry.csv`: measurements, reference positions, C++ estimated positions,
  timestamps, and cumulative RMSE used in the rendering.

Playback follows dataset timestamps at approximately real time, with a two-second
hold on the completed result. Every observation contributes to the filter and
metrics; animation frames sample the evolving history at the chosen frame rate.

Blue dots represent lidar positions. Orange crosses show radar range/bearing
converted to Cartesian coordinates, not a separate radar-only filter. The white
dashed line is ground truth, and the green line is the fused EKF estimate.
X/Y RMSE uses metres; VX/VY RMSE uses metres per second. All four metrics come
from the C++ server. Position RMSE is independently checked against the exported
estimates and ground truth before rendering.

This visualizes the synthetic reference dataset, not a Unity recording or a
real-world driving trial. Ground truth is used only for evaluation. Generated
outputs and the source dataset remain local in ignored directories. See the
[README](../README.md) for dataset provenance.

## Suggested post caption

Two sensors, one estimate: a C++17 Extended Kalman Filter fusing lidar and radar
measurements to track position and velocity. This animation shows the actual
filter output on Udacity's synthetic dataset, alongside noisy measurements and
ground truth. The project includes automated tests, Windows/Ubuntu CI, offline
evaluation, and integration with the Udacity simulator. Python visualizes the
results.
