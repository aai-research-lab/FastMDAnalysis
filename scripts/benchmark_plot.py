#!/usr/bin/env python3
"""Scaling plots of metrics produced by benchmark_scaling.py.

This version:
  - distinguishes libraries using BOTH color (matching the LOC / workflow plots)
    and marker shape / line style,
  - keeps the existing marker shapes and linestyles,
  - uses ONLY linear y-scales (no log plots),
  - labels the x-axis as "Number of Frames",
  - saves PNGs with simple names like combined_runtime.png, fastmdanalysis_runtime.png.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Same colors as benchmark_full_workflow / LOC plots
TOOL_COLORS: Dict[str, str] = {
    "fastmdanalysis": "#4472C4",
    "mdtraj": "#ED7D31",
    "mdanalysis": "#A5A5A5",
}

# Styles: marker + linestyle + color
LIBRARY_STYLES: Dict[str, Dict[str, str]] = {
    "fastmdanalysis": {
        "marker": "o",
        "linestyle": "--",
        "color": TOOL_COLORS["fastmdanalysis"],
    },
    "mdtraj": {
        "marker": "s",
        "linestyle": "-.",
        "color": TOOL_COLORS["mdtraj"],
    },
    "mdanalysis": {
        "marker": "^",
        "linestyle": ":",
        "color": TOOL_COLORS["mdanalysis"],
    },
}

# Global font sizes for all plots
TICK_FONT_SIZE = 18
LABEL_FONT_SIZE = 21
TITLE_FONT_SIZE = 24  # kept for completeness, but we do not draw titles


METRIC_PLOTS = {
    "runtime": {
        "ylabel": "Runtime (seconds)",
    },
    "memory": {
        "ylabel": "Peak Memory (MB)",
    },
}


def load_records(path: Path) -> Dict[str, List[dict]]:
    """Load and group records by library_key from a scaling_metrics.json file."""
    data = json.loads(path.read_text())
    records = data.get("records", [])
    grouped: Dict[str, List[dict]] = {}
    for record in records:
        library_key = record.get("library_key")
        grouped.setdefault(library_key, []).append(record)
    return grouped


def _sort_points(points: List[dict]) -> List[dict]:
    """Sort metric records by frame count."""
    return sorted(points, key=lambda rec: rec["frames"])


def _style_for(library_key: str) -> Dict[str, str]:
    """Return marker + linestyle + color for a given library.

    Falls back to circle + dashed + black if an unknown key appears.
    """
    return LIBRARY_STYLES.get(
        library_key,
        {
            "marker": "o",
            "linestyle": "--",
            "color": "#000000",
        },
    )


def _apply_slide_fonts(ax) -> None:
    """Apply slide-friendly font sizes to an axis."""
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.xaxis.label.set_fontsize(LABEL_FONT_SIZE)
    ax.yaxis.label.set_fontsize(LABEL_FONT_SIZE)
    # No title — we keep plots clean for embedding.


def _place_legend_top(ax, ncol: int | None = None) -> None:
    """Place a horizontal legend above the plot frame, centered, not overly wide."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    if ncol is None:
        ncol = len(labels)
    ax.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=ncol,
        frameon=False,
        fontsize=TICK_FONT_SIZE - 2,
        handlelength=2.5,
        borderaxespad=0.2,
    )


def plot_single_series(
    library_key: str,
    metric: str,
    points: List[dict],
    output_dir: Path,
) -> Path:
    """Plot per-library scaling for a single metric (runtime or memory)."""
    library_label = points[0]["library"] if points else library_key
    metadata = METRIC_PLOTS[metric]
    frames = [rec["frames"] for rec in points]
    means = [rec["mean"] for rec in points]
    errors = [rec.get("stderr", 0.0) for rec in points]

    style = _style_for(library_key)
    marker = style["marker"]
    linestyle = style["linestyle"]
    color = style["color"]

    fig, ax = plt.subplots(figsize=(8, 5))

    if any(err > 0 for err in errors):
        ax.errorbar(
            frames,
            means,
            yerr=errors,
            fmt=marker,
            linestyle=linestyle,
            color=color,
            linewidth=2.5,
            markersize=7,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=2.0,
            capsize=4,
            label=library_label,
        )
    else:
        ax.plot(
            frames,
            means,
            marker=marker,
            linestyle=linestyle,
            color=color,
            linewidth=2.5,
            markersize=7,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=2.0,
            label=library_label,
        )

    ax.set_xlabel("Number of Frames")
    ax.set_ylabel(metadata["ylabel"])
    ax.grid(alpha=0.3, linestyle="--")
    ax.set_xticks(frames)

    _apply_slide_fonts(ax)
    _place_legend_top(ax, ncol=1)

    output_dir.mkdir(parents=True, exist_ok=True)
    # e.g. fastmdanalysis_runtime.png, mdtraj_memory.png
    filename = f"{library_key}_{metric}.png"
    path = output_dir / filename
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_combined(
    metric: str,
    grouped: Dict[str, List[dict]],
    output_dir: Path,
) -> Path:
    """Plot all libraries on the same axes for a single metric (linear y-scale)."""
    metadata = METRIC_PLOTS[metric]
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    for library_key, points in grouped.items():
        metric_points = _sort_points([p for p in points if p["metric"] == metric])
        if not metric_points:
            continue

        frames = [rec["frames"] for rec in metric_points]
        means = [rec["mean"] for rec in metric_points]
        errors = [rec.get("stderr", 0.0) for rec in metric_points]

        style = _style_for(library_key)
        marker = style["marker"]
        linestyle = style["linestyle"]
        color = style["color"]
        label = metric_points[0]["library"]

        if any(err > 0 for err in errors):
            ax.errorbar(
                frames,
                means,
                yerr=errors,
                fmt=marker,
                linestyle=linestyle,
                color=color,
                linewidth=2.2,
                markersize=7,
                markerfacecolor="none",
                markeredgecolor=color,
                markeredgewidth=2.0,
                capsize=3,
                label=label,
            )
        else:
            ax.plot(
                frames,
                means,
                marker=marker,
                linestyle=linestyle,
                color=color,
                linewidth=2.2,
                markersize=7,
                markerfacecolor="none",
                markeredgecolor=color,
                markeredgewidth=2.0,
                label=label,
            )

    ax.set_xlabel("Number of Frames")
    ax.set_ylabel(metadata["ylabel"])
    ax.grid(alpha=0.3, linestyle="--")

    if grouped:
        frame_values = sorted(
            {
                rec["frames"]
                for records in grouped.values()
                for rec in records
                if rec.get("metric") == metric
            }
        )
        if frame_values:
            ax.set_xticks(frame_values)

    _apply_slide_fonts(ax)
    _place_legend_top(ax, ncol=None)

    output_dir.mkdir(parents=True, exist_ok=True)
    # e.g. combined_runtime.png, combined_memory.png
    filename = f"combined_{metric}.png"
    path = output_dir / filename
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def render_plots(
    json_path: Path,
    output_dir: Path,
    modes: Sequence[str],
) -> List[Path]:
    """Render per-library and/or combined plots from a JSON metrics file."""
    grouped = load_records(json_path)
    generated: List[Path] = []
    if "per-library" in modes:
        for library_key, points in grouped.items():
            for metric in METRIC_PLOTS:
                metric_points = _sort_points([p for p in points if p["metric"] == metric])
                if not metric_points:
                    continue
                generated.append(
                    plot_single_series(library_key, metric, metric_points, output_dir)
                )
    if "combined" in modes:
        for metric in METRIC_PLOTS:
            generated.append(
                plot_combined(metric, grouped, output_dir)
            )
    return generated


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scaling plots from benchmark_scaling output."
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("benchmark_output/scaling_metrics.json"),
        help="Path to scaling JSON file (default: benchmark_output/scaling_metrics.json).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmark_output"),
        help="Directory to store generated PNG files (default: benchmark_output).",
    )
    parser.add_argument(
        "--mode",
        choices=["per-library", "combined", "both"],
        default="both",
        help="Which plots to generate (default: both).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "both":
        modes = ["per-library", "combined"]
    else:
        modes = [args.mode]

    paths = render_plots(args.data, args.output_dir, modes)
    if not paths:
        print("No plots generated — verify the data file contains metrics.")
    else:
        print("Generated scaling plots:")
        for path in paths:
            print(f"  • {path}")


if __name__ == "__main__":
    main()
