"""Render real C++ EKF telemetry as a video, GIF, and documentation figure."""

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, PillowWriter
import numpy as np
import websocket


def collect(dataset, server, port):
    rows = [line.strip() for line in dataset.read_text().splitlines() if line.strip()]
    records = []
    for row in rows:
        fields = row.split()
        radar = fields[0] == "R"
        if fields[0] not in ("L", "R") or len(fields) != (11 if radar else 10):
            raise ValueError("Expected full Udacity L/R dataset rows")
        values = np.array(fields[1:], dtype=float)
        measurement = (values[0] * np.cos(values[1]), values[0] * np.sin(values[1])) if radar else values[:2]
        offset = 4 if radar else 3
        records.append(dict(sensor=fields[0], timestamp=int(fields[offset]),
                            mx=float(measurement[0]), my=float(measurement[1]),
                            gx=float(fields[offset + 1]), gy=float(fields[offset + 2])))
    if not records:
        raise ValueError("Dataset is empty")
    if any(b["timestamp"] < a["timestamp"] for a, b in zip(records, records[1:])):
        raise ValueError("Visualization requires chronological measurements")
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    # Silence per-connection logs; never kill an existing server on the chosen port.
    with subprocess.Popen([str(server.resolve()), str(port)], stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL, creationflags=flags) as process:
        connection = None
        try:
            for _ in range(50):
                if process.poll() is not None:
                    raise RuntimeError("Server exited; check the executable and choose an unused port")
                try:
                    connection = websocket.create_connection(
                        f"ws://127.0.0.1:{port}/socket.io/?EIO=3&transport=websocket",
                        timeout=5, http_no_proxy=["127.0.0.1"])
                    break
                except OSError:
                    time.sleep(0.1)
            if connection is None:
                raise RuntimeError("Server did not become ready")
            for row, record in zip(rows, records):
                connection.send("42" + json.dumps(["telemetry", {"sensor_measurement": row}]))
                reply = connection.recv()
                event, payload = json.loads(reply[2:])
                if event != "estimate_marker":
                    raise RuntimeError(f"Unexpected server response: {reply}")
                record.update(ex=payload["estimate_x"], ey=payload["estimate_y"],
                              rmse_x=payload["rmse_x"], rmse_y=payload["rmse_y"],
                              rmse_vx=payload["rmse_vx"], rmse_vy=payload["rmse_vy"])
        finally:
            if connection is not None:
                connection.close()
            process.terminate()
            process.wait(timeout=5)
    return records


def render(records, output):
    data = {key: np.array([r[key] for r in records]) for key in records[0]}
    t = (data["timestamp"] - data["timestamp"][0]) / 1e6
    if t[-1] <= 0:
        raise ValueError("Animation needs a positive time span")
    if not all(np.isfinite(v).all() for k, v in data.items() if k != "sensor"):
        raise ValueError("Non-finite telemetry")
    position_rmse = np.sqrt(np.mean((np.c_[data["ex"], data["ey"]] -
                                    np.c_[data["gx"], data["gy"]]) ** 2, axis=0))
    np.testing.assert_allclose(position_rmse, [data["rmse_x"][-1], data["rmse_y"][-1]], atol=1e-12)
    bg, fg, muted = "#0b1220", "#edf3fa", "#a5b5cb"
    cyan, amber, green = "#44c7f4", "#ffb454", "#63e6ac"
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": fg,
                         "axes.labelcolor": muted, "xtick.color": muted,
                         "ytick.color": muted, "font.size": 11})
    fig = plt.figure(figsize=(12.8, 9.6), facecolor=bg)
    fig.text(.07, .944, "TWO SENSORS. ONE ESTIMATE.", fontsize=25, weight="bold")
    fig.text(.07, .907, "Extended Kalman Filter  /  Lidar + radar sensor fusion", color=muted, fontsize=14)
    ax = fig.add_axes((.08, .26, .58, .58), facecolor=bg)
    ax.set_aspect("equal", adjustable="box")
    for axis, keys in (("x", ("mx", "gx", "ex")), ("y", ("my", "gy", "ey"))):
        values = np.concatenate([data[k] for k in keys])
        lo, hi = values.min(), values.max()
        getattr(ax, f"set_{axis}lim")(lo - .7, hi + .7)
    ax.set(xlabel="X position (m)", ylabel="Y position (m)")
    ax.grid(alpha=.12, color=muted)
    for spine in ax.spines.values():
        spine.set_color("#29364a")
    truth, = ax.plot([], [], color=fg, lw=2, ls="--", alpha=.85, label="Ground truth")
    estimate, = ax.plot([], [], color=green, lw=2.2, label="EKF estimate")
    lidar = ax.scatter([], [], s=23, color=cyan, alpha=.6, label="Lidar", zorder=2)
    radar = ax.scatter([], [], s=25, color=amber, marker="x", alpha=.65, label="Radar → Cartesian", zorder=2)
    head, = ax.plot([], [], "o", color=green, ms=8, markeredgecolor=bg, zorder=5)
    fig.legend(handles=[truth, estimate, lidar, radar], loc="lower left",
               bbox_to_anchor=(.07, .17), ncol=2, frameon=False, labelcolor=fg, fontsize=11)
    fig.text(.72, .81, "CUMULATIVE RMSE", fontsize=12, color=muted, weight="bold")
    metric_text = fig.text(.72, .76, "", va="top", fontsize=18, linespacing=1.65, fontfamily="DejaVu Sans Mono")
    fig.text(.72, .55, "Benchmark limits\nX / Y       0.11 m\nVX / VY  0.52 m/s", va="top", color=muted, linespacing=1.6)
    status = fig.text(.72, .41, "", color=green, weight="bold", fontsize=12)
    counter = fig.text(.72, .35, "", va="top", fontsize=12, linespacing=1.6)
    fig.text(.07, .117, "Actual C++ EKF output • Udacity synthetic Dataset 1", fontsize=12)
    fig.text(.07, .083, "Radar positions: x = ρ cos φ, y = ρ sin φ. Ground truth is used only for evaluation.",
             fontsize=10, color=muted)
    fig.text(.07, .046, "FelipeZerokun  /  Extended Kalman Filter", fontsize=10, color=muted)
    timeline = fig.add_axes((.08, .15, .83, .005), facecolor="#29364a")
    timeline.set(xlim=(0, t[-1]), ylim=(0, 1))
    timeline.axis("off")
    progress, = timeline.plot([], [], color=green, lw=5)

    def draw(index):
        stop = index + 1
        truth.set_data(data["gx"][:stop], data["gy"][:stop])
        estimate.set_data(data["ex"][:stop], data["ey"][:stop])
        for scatter, sensor in ((lidar, "L"), (radar, "R")):
            mask = data["sensor"][:stop] == sensor
            scatter.set_offsets(np.c_[data["mx"][:stop][mask], data["my"][:stop][mask]])
        head.set_data([data["ex"][index]], [data["ey"][index]])
        metric_text.set_text("\n".join(f"{label:<2}  {data[key][index]:.4f}" for label, key in
                                     (("X", "rmse_x"), ("Y", "rmse_y"), ("VX", "rmse_vx"), ("VY", "rmse_vy"))))
        passed = all(data[key][index] <= limit for key, limit in
                     zip(("rmse_x", "rmse_y", "rmse_vx", "rmse_vy"), (.11, .11, .52, .52)))
        status.set_text(("BENCHMARK PASS" if passed else "BENCHMARK NOT MET") if stop == len(records) else "ESTIMATING…")
        status.set_color(green if passed or stop < len(records) else amber)
        counter.set_text(f"{t[index]:05.2f} s / {t[-1]:.2f} s\n{stop:03d} / {len(records)} observations")
        progress.set_data([0, t[index]], [.5, .5])

    draw(len(records) - 1)
    fig.savefig(output / "ekf_summary.png", dpi=180, facecolor=bg)
    matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    for filename, fps, dpi, writer in (
        ("ekf_sensor_fusion.mp4", 12, 100, FFMpegWriter(fps=12, codec="libx264", bitrate=2400, extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"])),
        ("ekf_sensor_fusion.gif", 5, 65, PillowWriter(fps=5)),
    ):
        sample_times = np.r_[np.arange(0, t[-1], 1 / fps), np.repeat(t[-1], 2 * fps)]
        indices = np.maximum(0, np.searchsorted(t, sample_times, side="right") - 1)
        with writer.saving(fig, str(output / filename), dpi=dpi):
            for index in indices:
                draw(int(index))
                writer.grab_frame(facecolor=bg)
        print(f"Saved {output / filename}", flush=True)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--server", type=Path, required=True, help="Built ekf_simulator executable")
    parser.add_argument("--port", type=int, default=14568, help="Unused local port for the temporary server")
    parser.add_argument("--output", type=Path, default=Path("outputs/visualization"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = collect(args.dataset, args.server, args.port)
    with (args.output / "ekf_telemetry.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"Collected {len(records)} actual C++ estimates; rendering…", flush=True)
    render(records, args.output)


if __name__ == "__main__":
    main()
