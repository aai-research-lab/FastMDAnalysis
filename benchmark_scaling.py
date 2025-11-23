#!/usr/bin/env python3
"""Scaling benchmark for FastMDAnalysis, MDTraj, and MDAnalysis.

This script reuses the full-workflow benchmark implementations and simply varies
how many frames are processed so we can study scaling behaviour. Results are
persisted in JSON and CSV so downstream tooling (e.g., plot_scaling.py) can draw
line charts without rerunning expensive workloads.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as mp
import queue
import sys
import tempfile
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

import benchmark_full_workflow as base_benchmark

BENCHMARK_TOOLS = [
    "fastmdanalysis",
    "mdtraj",
    "mdanalysis",
]
DEFAULT_FRAME_COUNTS = [500, 1000, 2000, 4000, 5000]
DEFAULT_ITERATIONS = 3
RUNTIME_UNIT = "seconds"
MEMORY_UNIT = "MB"
WARMUP_TIMEOUT_SECONDS = 3600


def _available_tools(requested_tools: Sequence[str]) -> List[str]:
    if not requested_tools or requested_tools == ["all"]:
        requested_tools = BENCHMARK_TOOLS

    tools = []
    for tool in requested_tools:
        tk = tool.lower()
        if tk not in BENCHMARK_TOOLS:
            raise ValueError(f"Unknown tool '{tool}'. Valid options: {BENCHMARK_TOOLS}.")
        if tk == "mdanalysis" and not base_benchmark.HAS_MDANALYSIS:
            print("[warn] MDAnalysis not available. Skipping MDAnalysis runs.")
            continue
        tools.append(tk)
    if not tools:
        raise RuntimeError("No tools remain after filtering availability.")
    return tools


def _set_frame_slice(frame_count: int) -> Tuple[int, int, int]:
    previous = base_benchmark.FRAME_SLICE
    base_benchmark.FRAME_SLICE = (0, frame_count, 1)
    return previous


def _run_single_iteration(tool_key: str, frame_count: int) -> Tuple[float, float]:
    previous_slice = _set_frame_slice(frame_count)
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            if tool_key == "fastmdanalysis":
                runtime, memory = base_benchmark.benchmark_fastmda_single(tmpdir)
            elif tool_key == "mdtraj":
                runtime, memory = base_benchmark.benchmark_mdtraj_single(tmpdir)
            elif tool_key == "mdanalysis":
                runtime, memory = base_benchmark.benchmark_mdanalysis_single(tmpdir)
            else:
                raise ValueError(f"Unsupported tool '{tool_key}'.")
    finally:
        base_benchmark.FRAME_SLICE = previous_slice

    if runtime is None or memory is None:
        raise RuntimeError(f"{tool_key} run did not produce metrics. Check installation.")
    return runtime, memory


def _worker_run(tool_key: str, frame_count: int, warmup: bool, result_queue: "mp.Queue[dict]") -> None:
    try:
        if warmup:
            _run_single_iteration(tool_key, frame_count)
        runtime, memory = _run_single_iteration(tool_key, frame_count)
        result_queue.put({"ok": True, "runtime": runtime, "memory": memory})
    except Exception:
        result_queue.put({"ok": False, "traceback": traceback.format_exc()})


def _run_single_isolated(tool_key: str, frame_count: int, warmup: bool) -> Tuple[float, float]:
    ctx = mp.get_context("spawn")
    result_queue: mp.Queue = ctx.Queue()
    process = ctx.Process(target=_worker_run, args=(tool_key, frame_count, warmup, result_queue))
    process.start()
    result = None
    try:
        result = result_queue.get(timeout=WARMUP_TIMEOUT_SECONDS)
    except queue.Empty:
        process.terminate()
        process.join()
        raise RuntimeError(
            f"Timed out waiting for {tool_key} scaling run (frames={frame_count})."
        )
    process.join()
    if process.exitcode != 0:
        trace = ""
        if isinstance(result, dict):
            trace = result.get("traceback", "")
        raise RuntimeError(
            f"{tool_key} scaling process exited with code {process.exitcode}.\n{trace}"
        )
    if not isinstance(result, dict) or not result.get("ok"):
        trace = result.get("traceback", "Unknown error.") if isinstance(result, dict) else "Unknown error."
        raise RuntimeError(f"{tool_key} scaling run failed.\n{trace}")
    return float(result["runtime"]), float(result["memory"])


def _run_single(tool_key: str, frame_count: int, *, warmup: bool, isolate: bool) -> Tuple[float, float]:
    if isolate:
        return _run_single_isolated(tool_key, frame_count, warmup)
    if warmup:
        _run_single_iteration(tool_key, frame_count)
    return _run_single_iteration(tool_key, frame_count)


@dataclass
class MetricSummary:
    mean: float
    std: float
    stderr: float

    @classmethod
    def from_samples(cls, samples: Sequence[float]) -> "MetricSummary":
        if not samples:
            raise ValueError("No samples provided.")
        mean = float(np.mean(samples))
        if len(samples) == 1:
            std = 0.0
            stderr = 0.0
        else:
            std = float(np.std(samples, ddof=1))
            stderr = std / math.sqrt(len(samples))
        return cls(mean=mean, std=std, stderr=stderr)


def run_scaling_benchmark(
    frame_counts: Sequence[int],
    tools: Sequence[str],
    iterations: int,
    warmup: bool,
    isolate: bool,
) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []

    for tool in tools:
        display_name = {
            "fastmdanalysis": "FastMDAnalysis",
            "mdtraj": "MDTraj",
            "mdanalysis": "MDAnalysis",
        }[tool]
        print(f"\n=== {display_name} scaling benchmark ===")
        for frames in frame_counts:
            runtimes: List[float] = []
            memories: List[float] = []
            for iteration in range(iterations):
                include_warmup = warmup and iteration == 0
                label = f"[{display_name}] Frames={frames} Iteration {iteration + 1}/{iterations}"
                if include_warmup:
                    label += " (includes warm-up pass)"
                print(label)
                runtime, memory = _run_single(tool, frames, warmup=include_warmup, isolate=isolate)
                print(f"    runtime={runtime:.3f}s memory={memory:.1f} MB")
                runtimes.append(runtime)
                memories.append(memory)
            runtime_summary = MetricSummary.from_samples(runtimes)
            memory_summary = MetricSummary.from_samples(memories)
            records.append({
                "library": display_name,
                "library_key": tool,
                "metric": "runtime",
                "frames": frames,
                "mean": runtime_summary.mean,
                "stderr": runtime_summary.stderr,
                "std": runtime_summary.std,
                "unit": RUNTIME_UNIT,
                "iterations": iterations,
            })
            records.append({
                "library": display_name,
                "library_key": tool,
                "metric": "memory",
                "frames": frames,
                "mean": memory_summary.mean,
                "stderr": memory_summary.stderr,
                "std": memory_summary.std,
                "unit": MEMORY_UNIT,
                "iterations": iterations,
            })
    return records


def write_json(records: Sequence[Dict[str, object]], path: Path) -> None:
    ordered_frames: List[int] = []
    seen_frames = set()
    for rec in records:
        frames = rec["frames"]
        if frames not in seen_frames:
            ordered_frames.append(frames)
            seen_frames.add(frames)
    payload = {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "frame_counts": ordered_frames,
            "iterations": records[0]["iterations"] if records else 0,
            "dataset": "Ubiquitin",
            "notes": "Full workflow benchmark scaling study (single-iteration by default).",
        },
        "records": list(records),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"Saved JSON metrics -> {path}")


def write_csv(records: Sequence[Dict[str, object]], path: Path) -> None:
    if not records:
        raise ValueError("No records to persist.")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "library",
        "library_key",
        "metric",
        "frames",
        "mean",
        "stderr",
        "std",
        "unit",
        "iterations",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)
    print(f"Saved tabular metrics -> {path}")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate scaling metrics across frame counts.")
    parser.add_argument(
        "--frames",
        type=int,
        nargs="+",
        default=DEFAULT_FRAME_COUNTS,
        help="Frame counts to benchmark (default: 500 1000 2000 4000 5000).",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=DEFAULT_ITERATIONS,
        help="Number of iterations per point (default: 3).",
    )
    parser.add_argument(
        "--tools",
        nargs="+",
        default=["all"],
        help="Subset of tools to run (fastmdanalysis mdtraj mdanalysis).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("assets/benchmarks"),
        help="Directory to store metrics (default: assets/benchmarks).",
    )
    parser.add_argument(
        "--json-name",
        default="scaling_metrics.json",
        help="Filename for JSON output (default: scaling_metrics.json).",
    )
    parser.add_argument(
        "--csv-name",
        default="scaling_metrics.csv",
        help="Filename for CSV output (default: scaling_metrics.csv).",
    )
    parser.add_argument(
        "--warmup",
        dest="warmup",
        action="store_true",
        help="(Deprecated) Explicitly enable warm-up runs (default: enabled).",
    )
    parser.add_argument(
        "--no-warmup",
        dest="warmup",
        action="store_false",
        help="Disable warm-up passes before measurements.",
    )
    parser.add_argument(
        "--reuse-process",
        action="store_true",
        help=(
            "Reuse the current Python process for all runs (legacy behavior)."
            " Default spawns a fresh subprocess per run to keep memory baselines consistent."
        ),
    )
    parser.set_defaults(warmup=True)
    return parser.parse_args(argv)


def _ordered_positive_counts(values: Sequence[int]) -> List[int]:
    seen = set()
    ordered: List[int] = []
    for value in values:
        if value <= 0 or value in seen:
            continue
        ordered.append(value)
        seen.add(value)
    return ordered


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv or sys.argv[1:])
    frame_counts = _ordered_positive_counts(args.frames)
    if not frame_counts:
        raise ValueError("Provide at least one positive frame count.")
    tools = _available_tools(args.tools)
    iterations = max(1, args.iterations)

    print("=== Scaling Benchmark Configuration ===")
    print(f"Tools: {', '.join(tools)}")
    print(f"Frame counts: {frame_counts}")
    print(f"Iterations per point: {iterations}")

    isolate = not args.reuse_process
    records = run_scaling_benchmark(frame_counts, tools, iterations, args.warmup, isolate)

    output_dir: Path = args.output_dir
    write_json(records, output_dir / args.json_name)
    write_csv(records, output_dir / args.csv_name)

    print("\nBenchmarking complete. Use scripts/plot_scaling.py to generate plots.")


if __name__ == "__main__":
    main()
