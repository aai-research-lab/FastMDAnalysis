#!/usr/bin/env python3
"""Plot scaling metrics produced by benchmark_scaling.py."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext

LIBRARY_COLORS = {
    "fastmdanalysis": "#4472C4",  # blue
    "mdtraj": "#ED7D31",          # orange
    "mdanalysis": "#A5A5A5",      # gray
}

TICK_FONT_SIZE = 18
LABEL_FONT_SIZE = 21
TITLE_FONT_SIZE = 24

METRIC_PLOTS = {
    "runtime": {
        "ylabel": "Runtime (seconds)",
        "title": "Runtime Scaling",
    },
    "memory": {
        "ylabel": "Peak Memory (MB)",
        "title": "Memory Scaling",
    },
}


def load_records(path: Path) -> Dict[str, List[dict]]:
    data = json.loads(path.read_text())
    records = data.get("records", [])
    grouped: Dict[str, List[dict]] = {}
    for record in records:
        library_key = record.get("library_key")
        grouped.setdefault(library_key, []).append(record)
    return grouped


def _sort_points(points: List[dict]) -> List[dict]:
    return sorted(points, key=lambda rec: rec["frames"])


def plot_single_series(library_key: str, metric: str, points: List[dict], output_dir: Path) -> Path:
    library_label = points[0]["library"] if points else library_key
    metadata = METRIC_PLOTS[metric]
    frames = [rec["frames"] for rec in points]
    means = [rec["mean"] for rec in points]
    errors = [rec.get("stderr", 0.0) for rec in points]

    color = LIBRARY_COLORS.get(library_key, "#555555")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(frames, means, marker="o", linewidth=2.5, color=color)
    if any(err > 0 for err in errors):
        lower = [m - e for m, e in zip(means, errors)]
        upper = [m + e for m, e in zip(means, errors)]
        ax.fill_between(frames, lower, upper, color=color, alpha=0.2, label="± stderr")

    ax.set_xlabel("Frames")
    ax.set_ylabel(metadata["ylabel"])
    ax.set_title(f"{library_label} — {metadata['title']}")
    ax.grid(alpha=0.3, linestyle="--")
    ax.set_xticks(frames)
    ax.tick_params(axis="both", labelsize=11)

    if any(err > 0 for err in errors):
        ax.legend()

    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{library_key}_{metric}_scaling.png"
    path = output_dir / filename
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def _all_metric_values(grouped: Dict[str, List[dict]], metric: str) -> List[float]:
    values: List[float] = []
    for records in grouped.values():
        for rec in records:
            if rec.get("metric") == metric:
                values.append(float(rec.get("mean", 0.0)))
    return values


def _format_log_axis(ax, metric: str, values: List[float]) -> None:
    if not values:
        return
    min_val = max(min(values), 1e-3)
    max_val = max(values)
    lower_power = min(-1, math.floor(math.log10(min_val))) if min_val < 1 else math.floor(math.log10(min_val))
    upper_power = max(1, math.ceil(math.log10(max_val)))
    desired = {10.0, 100.0, 1000.0}
    ticks = {10 ** p for p in range(lower_power, upper_power + 1)}
    ticks.update(desired)
    ticks = sorted(ticks)
    ax.set_yticks([t for t in ticks if t > 0])
    formatter = LogFormatterMathtext()
    ax.yaxis.set_major_formatter(formatter)
    lower = min(min_val, 10.0)
    upper = max(max_val, 1000.0)
    ax.set_ylim(lower, upper)
    if metric == "memory":
        ax.set_ylim(10.0, max(upper, 1000.0))


def _apply_slide_fonts(ax) -> None:
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.xaxis.label.set_fontsize(LABEL_FONT_SIZE)
    ax.yaxis.label.set_fontsize(LABEL_FONT_SIZE)
    ax.title.set_fontsize(TITLE_FONT_SIZE)


def plot_combined(metric: str, grouped: Dict[str, List[dict]], output_dir: Path, y_scale: str, scale_label: str) -> Path:
    metadata = METRIC_PLOTS[metric]
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    metric_values = []

    for library_key, points in grouped.items():
        metric_points = _sort_points([p for p in points if p["metric"] == metric])
        if not metric_points:
            continue
        frames = [rec["frames"] for rec in metric_points]
        means = [rec["mean"] for rec in metric_points]
        metric_values.extend(means)
        errors = [rec.get("stderr", 0.0) for rec in metric_points]
        label = metric_points[0]["library"]
        color = LIBRARY_COLORS.get(library_key, "#555555")

        ax.plot(frames, means, marker="o", linewidth=2.2, color=color, label=label)
        if any(err > 0 for err in errors):
            lower = [m - e for m, e in zip(means, errors)]
            upper = [m + e for m, e in zip(means, errors)]
            ax.fill_between(frames, lower, upper, color=color, alpha=0.15)

    ax.set_xlabel("Frames")
    ax.set_ylabel(metadata["ylabel"])
    suffix = " (Log Scale)" if y_scale == "log" else " (Linear Scale)"
    ax.set_title(f"All Libraries — {metadata['title']}{suffix}")
    ax.grid(alpha=0.3, linestyle="--")
    if grouped:
        frame_values = sorted({rec["frames"] for records in grouped.values() for rec in records if rec.get("metric") == metric})
        if frame_values:
            ax.set_xticks(frame_values)
    if y_scale:
        ax.set_yscale(y_scale)
        if y_scale == "log":
            _format_log_axis(ax, metric, metric_values or _all_metric_values(grouped, metric))
    _apply_slide_fonts(ax)
    ax.legend(fontsize=TICK_FONT_SIZE - 2)

    output_dir.mkdir(parents=True, exist_ok=True)
    suffix_token = f"_{scale_label}" if scale_label else ""
    filename = f"combined_{metric}_scaling{suffix_token}.png"
    path = output_dir / filename
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def render_plots(
    json_path: Path,
    output_dir: Path,
    modes: Sequence[str],
    combined_scales: Sequence[str],
) -> List[Path]:
    grouped = load_records(json_path)
    generated: List[Path] = []
    if "per-library" in modes:
        for library_key, points in grouped.items():
            for metric in METRIC_PLOTS:
                metric_points = _sort_points([p for p in points if p["metric"] == metric])
                if not metric_points:
                    continue
                generated.append(plot_single_series(library_key, metric, metric_points, output_dir))
    if "combined" in modes:
        for metric in METRIC_PLOTS:
            for scale in combined_scales:
                generated.append(plot_combined(metric, grouped, output_dir, scale, scale))
    return generated


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot scaling figures from benchmark_scaling output.")
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("assets/benchmarks/scaling_metrics.json"),
        help="Path to scaling JSON file (default: assets/benchmarks/scaling_metrics.json).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("assets/benchmarks"),
        help="Directory to store generated PNG files (default: assets/benchmarks).",
    )
    parser.add_argument(
        "--mode",
        choices=["per-library", "combined", "both"],
        default="per-library",
        help="Which plots to generate (default: per-library).",
    )
    parser.add_argument(
        "--combined-y-scale",
        choices=["linear", "log"],
        nargs="+",
        default=["log", "linear"],
        help="Y-axis scales for combined plots (default: log linear).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "both":
        modes = ["per-library", "combined"]
    else:
        modes = [args.mode]
    paths = render_plots(args.data, args.output_dir, modes, args.combined_y_scale)
    if not paths:
        print("No plots generated — verify the data file contains metrics.")
    else:
        print("Generated plots:")
        for path in paths:
            print(f"  • {path}")


if __name__ == "__main__":
    main()
