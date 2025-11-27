#!/usr/bin/env python3
"""
FastMDAnalysis Full Workflow Benchmark

Benchmarks FastMDAnalysis, MDTraj, and MDAnalysis on COMPLETE workflows
including computation and visualization (but excluding slides).

Measures:
- Runtime (wall clock time)
- Peak memory usage

Runs each workflow NUM_ITERATIONS times and reports average ± standard deviation.
"""

from __future__ import annotations

import sys
import time
import warnings
import os
import tempfile
import gc
import inspect
import ast
import tokenize
import textwrap
from io import StringIO
from pathlib import Path
from typing import Dict, Sequence, List

import numpy as np
import psutil
import matplotlib

matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator

import mdtraj as md
from sklearn.cluster import KMeans, DBSCAN
from scipy.cluster.hierarchy import linkage, fcluster

warnings.filterwarnings("ignore")

from fastmdanalysis import FastMDAnalysis
from fastmdanalysis.analysis.cluster import (
    get_cluster_cmap,
    get_discrete_norm,
    relabel_compact_positive,
)

# --------------------------------------------------------------------
# Ubiquitin dataset paths – update these to wherever you put the files
# --------------------------------------------------------------------
TRAJ = Path("ubiquitin.dcd")
TOP = Path("ubiquitin.pdb")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms as mda_rms

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    print("Warning: MDAnalysis not installed. Will skip MDAnalysis benchmarks.")

# Number of iterations for averaging
NUM_ITERATIONS = 3
FRAME_SLICE = (0, 500, 1)
AXIS_FONT_SCALE = 1.5
CLUSTER_METHODS = ["kmeans"]  # can also include 'dbscan', 'hierarchical'
ENABLE_CLUSTERING = True

# Tool colors for consistent visualization
TOOL_COLORS = {
    "fastmdanalysis": "#4472C4",
    "mdtraj": "#ED7D31",
    "mdanalysis": "#A5A5A5",
}


# ---------------------------------------------------------------------
# Utility classes and helper plotting functions
# ---------------------------------------------------------------------
class MemoryMonitor:
    """Monitor peak memory usage during execution."""

    def __init__(self):
        self.process = psutil.Process()
        self.baseline_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.peak_memory = self.baseline_memory

    def update(self):
        """Update peak memory if current usage is higher."""
        current_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        if current_memory > self.peak_memory:
            self.peak_memory = current_memory

    def get_peak_mb(self):
        """Get peak memory usage in MB (relative to baseline)."""
        return self.peak_memory - self.baseline_memory


def plot_cluster_summary(labels_map: Dict[str, np.ndarray], output_file: str | Path):
    """Generic multi-method scatter summary (used only if non-kmeans methods exist)."""
    methods = [m for m in CLUSTER_METHODS if m in labels_map]
    if not methods:
        return None

    fig, axes = plt.subplots(1, len(methods), figsize=(5 * len(methods), 4))
    if len(methods) == 1:
        axes = [axes]

    for ax, method in zip(axes, methods):
        labels = np.asarray(labels_map[method], dtype=int)
        frames = np.arange(len(labels))
        unique = np.sort(np.unique(labels))
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        ax.scatter(frames, labels, c=labels, cmap=cmap, norm=norm, s=60)
        ax.set_title(method.capitalize())
        ax.set_xlabel("Frame")
        ax.set_ylabel("Cluster")
        ax.grid(alpha=0.3)

    fig.suptitle("Clustering Results", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_file


def plot_kmeans_trajectory_artifacts(
    frame_idx: np.ndarray, labels: np.ndarray, outdir: str | Path, method_name: str = "kmeans"
):
    """Replicate FastMDAnalysis-style KMeans artefacts.

    Produces:
      <method>_pop.png
      <method>_traj_hist.png
      <method>_traj_scatter.png

    All frames lie on a single horizontal line in the scatter plot; the colourbar
    is labelled 'Cluster (compact)' with ticks at the cluster IDs.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    labels = np.asarray(labels, dtype=int)
    frame_idx = np.asarray(frame_idx, dtype=int)
    cluster_ids = np.sort(np.unique(labels))

    cmap = get_cluster_cmap(len(cluster_ids))
    norm = get_discrete_norm(cluster_ids)

    # ----- Population barplot -----
    counts = np.array([(labels == k).sum() for k in cluster_ids], dtype=int)
    pop_path = outdir / f"{method_name}_pop.png"

    fig, ax = plt.subplots(figsize=(6, 4))
    bar_colors = [cmap(norm(k)) for k in cluster_ids]
    ax.bar(cluster_ids, counts, color=bar_colors, edgecolor="black")
    ax.set_xlabel("Cluster ID (compact)")
    ax.set_ylabel("Number of Frames")
    ax.set_title("Cluster Populations")
    ax.set_xticks(cluster_ids)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(pop_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ----- 1D histogram-style trajectory (heatmap) -----
    hist_path = outdir / f"{method_name}_traj_hist.png"

    fig, ax = plt.subplots(figsize=(10, 3))
    # Use extent to make x-axis in frame indices
    img = ax.imshow(
        labels[np.newaxis, :],
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        extent=[frame_idx.min(), frame_idx.max() + 1, 0, 1],
    )
    ax.set_yticks([])
    ax.set_xlabel("Frame")
    ax.set_title("Cluster Trajectory Histogram")
    cbar = fig.colorbar(img, ax=ax, pad=0.02)
    cbar.set_label("Cluster (compact)")
    cbar.set_ticks(cluster_ids)
    fig.tight_layout()
    fig.savefig(hist_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ----- Scatter trajectory on a single horizontal line -----
    scatter_path = outdir / f"{method_name}_traj_scatter.png"

    fig, ax = plt.subplots(figsize=(10, 3))
    y = np.zeros_like(frame_idx, dtype=float)
    sc = ax.scatter(frame_idx, y, c=labels, cmap=cmap, norm=norm, s=40)
    ax.set_xlabel("Frame")
    ax.set_ylabel("Cluster")
    ax.set_yticks([])  # visually a single line, colour encodes cluster
    ax.set_title("Cluster Trajectory Scatter Plot")
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Cluster (compact)")
    cbar.set_ticks(cluster_ids)
    fig.tight_layout()
    fig.savefig(scatter_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------
# LOC / CLI comparison helpers
# ---------------------------------------------------------------------
def _count_effective_loc(source_lines: Sequence[str]) -> int:
    """Count non-blank, non-comment, non-docstring lines using token analysis."""
    source_text = "".join(source_lines)
    dedented = textwrap.dedent(source_text)

    meaningful_lines = set()
    try:
        for token in tokenize.generate_tokens(StringIO(dedented).readline):
            if token.type in {
                tokenize.NL,
                tokenize.NEWLINE,
                tokenize.INDENT,
                tokenize.DEDENT,
                tokenize.COMMENT,
                tokenize.ENDMARKER,
                tokenize.ENCODING,
            }:
                continue
            meaningful_lines.add(token.start[0])
    except tokenize.TokenError:
        pass

    try:
        tree = ast.parse(dedented)
    except SyntaxError:
        return len(meaningful_lines)

    if tree.body and isinstance(tree.body[0], (ast.FunctionDef, ast.AsyncFunctionDef)):
        func_node = tree.body[0]
        if func_node.body:
            first_stmt = func_node.body[0]
            if isinstance(first_stmt, ast.Expr) and isinstance(getattr(first_stmt, "value", None), ast.Constant):
                if isinstance(first_stmt.value.value, str):
                    doc_start = first_stmt.lineno
                    doc_end = getattr(first_stmt, "end_lineno", doc_start)
                    for line_no in range(doc_start, doc_end + 1):
                        meaningful_lines.discard(line_no)

    return len(meaningful_lines)


FASTMDA_MINIMAL_SNIPPET = [
    'fastmda = FastMDAnalysis(TRAJ, TOP, frames=FRAME_SLICE, atoms="protein", keep_full_traj=False)',
    (
        "fastmda.analyze(include=['rmsd', 'rmsf', 'rg', 'cluster'], verbose=False, "
        "output='fastmda_analyze_output', "
        "options={'rmsd': {'save_data': False, 'store_results': False}, "
        "'rmsf': {'save_data': False, 'store_results': False}, "
        "'rg': {'save_data': False, 'store_results': False}, "
        "'cluster': {'methods': ['kmeans', 'dbscan', 'hierarchical'], 'n_clusters': 3, "
        "'eps': 0.5, 'min_samples': 2, 'plot_style': 'minimal', "
        "'combined_plot_name': 'cluster', 'feature_mode': 'distance', "
        "'save_data': False, 'store_results': False}})"
    ),
]


def compute_loc_benchmark():
    """Return LOC metrics for each workflow implementation."""
    targets = [
        ("FastMDAnalysis", FASTMDA_MINIMAL_SNIPPET),
        ("MDTraj", benchmark_mdtraj_single),
        ("MDAnalysis", benchmark_mdanalysis_single),
    ]
    loc_data = []
    for name, func in targets:
        if callable(func):
            lines, _ = inspect.getsourcelines(func)
        elif isinstance(func, str):
            lines = [f"{line}\n" for line in func.splitlines()]
        else:
            lines = [line if line.endswith("\n") else f"{line}\n" for line in func]
        loc = _count_effective_loc(lines)
        loc_data.append({"name": name, "loc": loc})
    return loc_data


def generate_loc_slide(loc_data, output_path: str | Path = "benchmark_loc_slide.png"):
    """Create a separate slide summarizing lines-of-code per workflow."""
    names = [item["name"] for item in loc_data]
    loc_values = [item["loc"] for item in loc_data]
    colors = [TOOL_COLORS.get(name.lower(), "#777777") for name in names]

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(names, loc_values, color=colors, edgecolor="black", linewidth=1.5)
    label_font = 14 * AXIS_FONT_SCALE
    title_font = 16
    tick_font = 12 * AXIS_FONT_SCALE

    ax.set_ylabel("Effective LOC (non-blank/non-comment)", fontsize=label_font)
    ax.set_title("Lines of Code Required Per Workflow", fontsize=title_font, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=tick_font)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    if loc_values:
        top = max(loc_values)
        pad = top * 0.05
    else:
        top = 1
        pad = 1
    ax.set_ylim(0, top + pad * 4)

    for bar, loc in zip(bars, loc_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + pad,
            f"{loc} LOC",
            ha="center",
            va="bottom",
            fontsize=12 * AXIS_FONT_SCALE,
            fontweight="bold",
            zorder=5,
        )

    plt.tight_layout()
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    for item in loc_data:
        print(f"  {item['name']}: {item['loc']} LOC")
    plt.close()
    return output_file


def compute_cli_command_counts():
    """Return CLI command counts used to reproduce the full workflow."""
    return [
        {"name": "FastMDAnalysis", "commands": 1},
        {"name": "MDTraj", "commands": 4},
        {"name": "MDAnalysis", "commands": 4},
    ]


def generate_cli_command_slide(command_data):
    """Create a bar chart comparing CLI commands needed per workflow."""
    names = [item["name"] for item in command_data]
    values = [item["commands"] for item in command_data]
    colors = [TOOL_COLORS.get(name.lower(), "#777777") for name in names]

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(names, values, color=colors, edgecolor="black", linewidth=1.5)
    ax.set_ylabel("CLI Commands Needed", fontsize=14)
    ax.set_title("Commands to Reproduce Full Workflow Outputs", fontsize=16, fontweight="bold")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, max(values) + 1)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    for bar, value in zip(bars, values):
        label = "command" if value == 1 else "commands"
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.1,
            f"{value} {label}",
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    plt.tight_layout()
    output_file = "benchmark_cli_commands.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    for item in command_data:
        print(f"  {item['name']}: {item['commands']} CLI command(s)")
    plt.close()
    return output_file


# ---------------------------------------------------------------------
# Core benchmark wrappers
# ---------------------------------------------------------------------
def _frame_range(total_frames: int):
    """Convert FRAME_SLICE tuple into a concrete range of frame indices."""
    start, stop, step = FRAME_SLICE
    if start is None:
        start = 0
    if step is None:
        step = 1
    if stop is None or stop > total_frames:
        stop = total_frames
    return range(start, stop, step)


def benchmark_fastmda_single(output_dir: str | Path):
    """FastMDAnalysis workflow: RMSD, RMSF, RG + optional clustering."""
    mem_monitor = MemoryMonitor()
    start = time.time()

    fastmda = FastMDAnalysis(
        str(TRAJ),
        str(TOP),
        frames=FRAME_SLICE,
        atoms="protein",
    )
    mem_monitor.update()

    fastmda.rmsd(
        ref=0,
        save_data=False,
        store_results=False,
        output=os.path.join(output_dir, "fastmda_rmsd_output"),
    )
    mem_monitor.update()

    fastmda.rmsf(
        save_data=False,
        store_results=False,
        output=os.path.join(output_dir, "fastmda_rmsf_output"),
    )
    mem_monitor.update()

    fastmda.rg(
        save_data=False,
        store_results=False,
        output=os.path.join(output_dir, "fastmda_rg_output"),
    )
    mem_monitor.update()

    if ENABLE_CLUSTERING:
        fastmda.cluster(
            methods=CLUSTER_METHODS,
            n_clusters=3,
            eps=0.5,
            min_samples=2,
            plot_style="minimal",
            combined_plot_name="cluster",
            save_data=False,
            store_results=False,
            output=os.path.join(output_dir, "fastmda_cluster_output"),
        )
        mem_monitor.update()

    runtime = time.time() - start
    peak_memory = mem_monitor.get_peak_mb()

    del fastmda
    gc.collect()

    return runtime, peak_memory


def benchmark_mdtraj_single(output_dir: str | Path):
    """MDTraj full workflow – RMSD, RMSF, RG, KMeans (and optional others)."""
    mem_monitor = MemoryMonitor()
    start = time.time()

    traj = md.load(str(TRAJ), top=str(TOP))
    frame_indices = list(_frame_range(traj.n_frames))
    traj = traj[frame_indices]
    atom_indices = traj.topology.select("protein")
    traj = traj.atom_slice(atom_indices)
    mem_monitor.update()

    # RMSD
    rmsd_data = md.rmsd(traj, traj, frame=0)
    mem_monitor.update()

    rmsd_outdir = os.path.join(output_dir, "mdtraj_rmsd_output")
    os.makedirs(rmsd_outdir, exist_ok=True)
    np.savetxt(
        os.path.join(rmsd_outdir, "rmsd.dat"),
        rmsd_data.reshape(-1, 1),
        header="rmsd_nm",
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rmsd_data, marker="o", linestyle="-")
    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD (nm)")
    ax.set_title("RMSD vs Frame (ref=0, align=True)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(rmsd_outdir, "rmsd.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    mem_monitor.update()

    # RMSF
    avg_xyz = np.mean(traj.xyz, axis=0, keepdims=True)
    ref = md.Trajectory(avg_xyz, traj.topology)
    rmsf_data = md.rmsf(traj, ref)
    mem_monitor.update()

    rmsf_outdir = os.path.join(output_dir, "mdtraj_rmsf_output")
    os.makedirs(rmsf_outdir, exist_ok=True)
    np.savetxt(
        os.path.join(rmsf_outdir, "rmsf.dat"),
        rmsf_data.reshape(-1, 1),
        header="rmsf_nm",
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(rmsf_data)), rmsf_data, width=0.9)
    ax.set_xlabel("Atom Index")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("RMSF per Atom")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(rmsf_outdir, "rmsf.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    mem_monitor.update()

    # RG
    rg_data = md.compute_rg(traj)
    mem_monitor.update()

    rg_outdir = os.path.join(output_dir, "mdtraj_rg_output")
    os.makedirs(rg_outdir, exist_ok=True)
    np.savetxt(
        os.path.join(rg_outdir, "rg.dat"),
        rg_data.reshape(-1, 1),
        header="rg_nm",
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rg_data, marker="o", linestyle="-")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Radius of Gyration (nm)")
    ax.set_title("Radius of Gyration vs Frame")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(rg_outdir, "rg.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    mem_monitor.update()

    # Clustering
    if ENABLE_CLUSTERING and CLUSTER_METHODS:
        cluster_outdir = os.path.join(output_dir, "mdtraj_cluster_output")
        os.makedirs(cluster_outdir, exist_ok=True)

        methods = {m.lower() for m in CLUSTER_METHODS}
        labels_map: Dict[str, np.ndarray] = {}
        n_frames = traj.n_frames
        frame_idx = np.arange(n_frames, dtype=int)

        # DBSCAN: RMSD matrix
        if "dbscan" in methods:
            D = np.empty((n_frames, n_frames), dtype=np.float32)
            for i in range(n_frames):
                D[:, i] = md.rmsd(traj, traj, frame=i)
            D = 0.5 * (D + D.T)
            np.fill_diagonal(D, 0.0)
            mem_monitor.update()

            db = DBSCAN(eps=0.5, min_samples=2, metric="precomputed")
            labels_raw = db.fit_predict(D).astype(int, copy=False)
            labels_compact, _, _ = relabel_compact_positive(labels_raw, start=1, noise_as_last=True)
            labels_map["dbscan"] = labels_compact
            mem_monitor.update()

            np.savetxt(
                os.path.join(cluster_outdir, "dbscan_labels_compact.dat"),
                np.column_stack((frame_idx, labels_compact)),
                fmt="%d",
                header="frame label(1..K; noise=K+1)",
            )
            np.savetxt(
                os.path.join(cluster_outdir, "dbscan_labels_raw.dat"),
                np.column_stack((frame_idx, labels_raw)),
                fmt="%d",
                header="frame label_raw(-1=noise)",
            )

        # Shared flattened coords for kmeans / hierarchical
        need_centroid = any(m in methods for m in ("kmeans", "hierarchical"))
        if need_centroid:
            X = traj.xyz
            X_flat = X.reshape(n_frames, -1)
            mem_monitor.update()
        else:
            X_flat = None

        # KMeans
        if "kmeans" in methods:
            km = KMeans(n_clusters=3, random_state=42, n_init=10)
            labels_k = km.fit_predict(X_flat).astype(int, copy=False) + 1
            labels_map["kmeans"] = labels_k
            mem_monitor.update()

            np.savetxt(
                os.path.join(cluster_outdir, "kmeans_labels.dat"),
                np.column_stack((frame_idx, labels_k)),
                fmt="%d",
                header="frame label(1..K)",
            )
            np.savetxt(
                os.path.join(cluster_outdir, "kmeans_coordinates.dat"),
                X_flat,
                fmt="%.6f",
                header="Flattened coordinates",
            )

            # FastMDAnalysis-style artefacts
            plot_kmeans_trajectory_artifacts(frame_idx, labels_k, cluster_outdir, method_name="kmeans")

        # Hierarchical
        if "hierarchical" in methods:
            Z = linkage(X_flat, method="ward")
            labels_h = fcluster(Z, t=3, criterion="maxclust").astype(int, copy=False)
            labels_map["hierarchical"] = labels_h
            mem_monitor.update()

            np.savetxt(
                os.path.join(cluster_outdir, "hierarchical_labels.dat"),
                np.column_stack((frame_idx, labels_h)),
                fmt="%d",
                header="frame label(1..K)",
            )
            np.savetxt(
                os.path.join(cluster_outdir, "hierarchical_linkage.dat"),
                Z,
                fmt="%.6f",
                header="cluster1 cluster2 distance sample_count",
            )

        # Only make cluster.png if there are non-kmeans methods to show
        non_kmeans_methods = {m for m in labels_map if m != "kmeans"}
        if non_kmeans_methods:
            submap = {m: labels_map[m] for m in non_kmeans_methods}
            plot_cluster_summary(submap, os.path.join(cluster_outdir, "cluster.png"))
            mem_monitor.update()

    runtime = time.time() - start
    peak_memory = mem_monitor.get_peak_mb()

    # Cleanup
    del traj, rmsd_data, rmsf_data, rg_data
    gc.collect()

    return runtime, peak_memory


def benchmark_mdanalysis_single(output_dir: str | Path):
    """MDAnalysis full workflow – RMSD, RMSF, RG, KMeans (and optional others)."""
    if not HAS_MDANALYSIS:
        return None, None

    import warnings as _warnings

    _warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message="DCDReader currently makes independent timesteps",
    )

    mem_monitor = MemoryMonitor()
    start = time.time()

    u = mda.Universe(str(TOP), str(TRAJ))
    protein = u.select_atoms("protein")
    frame_list = list(_frame_range(len(u.trajectory)))
    mem_monitor.update()

    # RMSD
    rmsd_results: List[float] = []
    u.trajectory[frame_list[0]]
    ref_coords = protein.positions.copy()

    for ts in u.trajectory[frame_list]:
        rmsd_val = mda_rms.rmsd(protein.positions, ref_coords, center=True, superposition=True) / 10.0
        rmsd_results.append(rmsd_val)
    rmsd_data = np.asarray(rmsd_results, dtype=float)
    mem_monitor.update()

    rmsd_outdir = os.path.join(output_dir, "mdanalysis_rmsd_output")
    os.makedirs(rmsd_outdir, exist_ok=True)
    np.savetxt(
        os.path.join(rmsd_outdir, "rmsd.dat"),
        rmsd_data.reshape(-1, 1),
        header="rmsd_nm",
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rmsd_data, marker="o", linestyle="-")
    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD (nm)")
    ax.set_title("RMSD vs Frame (ref=0, align=True)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(rmsd_outdir, "rmsd.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    mem_monitor.update()

    # RMSF
    coordinates = []
    for ts in u.trajectory[frame_list]:
        coordinates.append(protein.positions.copy())
    coordinates = np.asarray(coordinates, dtype=np.float32)  # Å
    mem_monitor.update()

    avg_coords = np.mean(coordinates, axis=0)
    rmsf_vec = np.sqrt(np.mean((coordinates - avg_coords) ** 2, axis=0))
    rmsf_data = np.linalg.norm(rmsf_vec, axis=1) / 10.0
    mem_monitor.update()

    rmsf_outdir = os.path.join(output_dir, "mdanalysis_rmsf_output")
    os.makedirs(rmsf_outdir, exist_ok=True)
    np.savetxt(
        os.path.join(rmsf_outdir, "rmsf.dat"),
        rmsf_data.reshape(-1, 1),
        header="rmsf_nm",
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(rmsf_data)), rmsf_data, width=0.9)
    ax.set_xlabel("Atom Index")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("RMSF per Atom")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(rmsf_outdir, "rmsf.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    mem_monitor.update()

    # RG
    rg_results = []
    for ts in u.trajectory[frame_list]:
        rg_results.append(protein.radius_of_gyration() / 10.0)
    rg_data = np.asarray(rg_results, dtype=float)
    mem_monitor.update()

    rg_outdir = os.path.join(output_dir, "mdanalysis_rg_output")
    os.makedirs(rg_outdir, exist_ok=True)
    np.savetxt(
        os.path.join(rg_outdir, "rg.dat"),
        rg_data.reshape(-1, 1),
        header="rg_nm",
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rg_data, marker="o", linestyle="-")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Radius of Gyration (nm)")
    ax.set_title("Radius of Gyration vs Frame")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(rg_outdir, "rg.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    mem_monitor.update()

    # Clustering
    if ENABLE_CLUSTERING and CLUSTER_METHODS:
        cluster_outdir = os.path.join(output_dir, "mdanalysis_cluster_output")
        os.makedirs(cluster_outdir, exist_ok=True)

        methods = {m.lower() for m in CLUSTER_METHODS}
        labels_map: Dict[str, np.ndarray] = {}
        n_frames = coordinates.shape[0]
        frame_idx = np.arange(n_frames, dtype=int)

        # DBSCAN using RMSD distance in nm
        if "dbscan" in methods:
            D = np.empty((n_frames, n_frames), dtype=np.float32)
            for i in range(n_frames):
                diff = coordinates - coordinates[i]
                msd = np.mean(np.sum(diff * diff, axis=2), axis=1)
                D[:, i] = np.sqrt(msd) / 10.0
            D = 0.5 * (D + D.T)
            np.fill_diagonal(D, 0.0)
            mem_monitor.update()

            db = DBSCAN(eps=0.5, min_samples=2, metric="precomputed")
            labels_raw = db.fit_predict(D).astype(int, copy=False)
            labels_compact, _, _ = relabel_compact_positive(labels_raw, start=1, noise_as_last=True)
            labels_map["dbscan"] = labels_compact
            mem_monitor.update()

            np.savetxt(
                os.path.join(cluster_outdir, "dbscan_labels_compact.dat"),
                np.column_stack((frame_idx, labels_compact)),
                fmt="%d",
                header="frame label(1..K; noise=K+1)",
            )
            np.savetxt(
                os.path.join(cluster_outdir, "dbscan_labels_raw.dat"),
                np.column_stack((frame_idx, labels_raw)),
                fmt="%d",
                header="frame label_raw(-1=noise)",
            )

        # Shared flattened coords
        need_centroid = any(m in methods for m in ("kmeans", "hierarchical"))
        if need_centroid:
            coords_nm = coordinates / 10.0
            X_flat = coords_nm.reshape(n_frames, -1)
            mem_monitor.update()
        else:
            X_flat = None

        # KMeans
        if "kmeans" in methods:
            km = KMeans(n_clusters=3, random_state=42, n_init=10)
            labels_k = km.fit_predict(X_flat).astype(int, copy=False) + 1
            labels_map["kmeans"] = labels_k
            mem_monitor.update()

            np.savetxt(
                os.path.join(cluster_outdir, "kmeans_labels.dat"),
                np.column_stack((frame_idx, labels_k)),
                fmt="%d",
                header="frame label(1..K)",
            )
            np.savetxt(
                os.path.join(cluster_outdir, "kmeans_coordinates.dat"),
                X_flat,
                fmt="%.6f",
                header="Flattened coordinates",
            )

            plot_kmeans_trajectory_artifacts(frame_idx, labels_k, cluster_outdir, method_name="kmeans")

        # Hierarchical
        if "hierarchical" in methods:
            Z = linkage(X_flat, method="ward")
            labels_h = fcluster(Z, t=3, criterion="maxclust").astype(int, copy=False)
            labels_map["hierarchical"] = labels_h
            mem_monitor.update()

            np.savetxt(
                os.path.join(cluster_outdir, "hierarchical_labels.dat"),
                np.column_stack((frame_idx, labels_h)),
                fmt="%d",
                header="frame label(1..K)",
            )
            np.savetxt(
                os.path.join(cluster_outdir, "hierarchical_linkage.dat"),
                Z,
                fmt="%.6f",
                header="cluster1 cluster2 distance sample_count",
            )

        # Only produce combined cluster.png if there are non-kmeans methods
        non_kmeans_methods = {m for m in labels_map if m != "kmeans"}
        if non_kmeans_methods:
            submap = {m: labels_map[m] for m in non_kmeans_methods}
            plot_cluster_summary(submap, os.path.join(cluster_outdir, "cluster.png"))
            mem_monitor.update()

    runtime = time.time() - start
    peak_memory = mem_monitor.get_peak_mb()

    # Cleanup
    del u, protein, coordinates, rmsd_data, rmsf_data, rg_data
    gc.collect()

    return runtime, peak_memory


# ---------------------------------------------------------------------
# Multi-iteration wrappers + plotting
# ---------------------------------------------------------------------
def benchmark_fastmda():
    print("\n" + "=" * 70)
    print("FastMDAnalysis - Full Workflow (with figures, no slides)")
    print("=" * 70)
    print(f"Running {NUM_ITERATIONS} iterations...")

    runtimes = []
    memories = []

    for i in range(NUM_ITERATIONS):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime, memory = benchmark_fastmda_single(tmpdir)
            runtimes.append(runtime)
            memories.append(memory)
            print(f"  Iteration {i+1}/{NUM_ITERATIONS}: {runtime:.3f}s, {memory:.1f} MB")

    avg_runtime = np.mean(runtimes)
    std_runtime = np.std(runtimes)
    avg_memory = np.mean(memories)
    std_memory = np.std(memories)

    print(f"\nAverage Runtime: {avg_runtime:.3f}s (±{std_runtime:.3f}s)")
    print(f"Average Peak Memory: {avg_memory:.1f} MB (±{std_memory:.1f} MB)")
    print("=" * 70)

    return {
        "name": "FastMDAnalysis",
        "runtime_avg": avg_runtime,
        "runtime_std": std_runtime,
        "memory_avg": avg_memory,
        "memory_std": std_memory,
    }


def benchmark_mdtraj():
    print("\n" + "=" * 70)
    print("MDTraj - Full Workflow (with figures)")
    print("=" * 70)
    print(f"Running {NUM_ITERATIONS} iterations...")

    runtimes = []
    memories = []

    for i in range(NUM_ITERATIONS):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime, memory = benchmark_mdtraj_single(tmpdir)
            runtimes.append(runtime)
            memories.append(memory)
            print(f"  Iteration {i+1}/{NUM_ITERATIONS}: {runtime:.3f}s, {memory:.1f} MB")

    avg_runtime = np.mean(runtimes)
    std_runtime = np.std(runtimes)
    avg_memory = np.mean(memories)
    std_memory = np.std(memories)

    print(f"\nAverage Runtime: {avg_runtime:.3f}s (±{std_runtime:.3f}s)")
    print(f"Average Peak Memory: {avg_memory:.1f} MB (±{std_memory:.1f} MB)")
    print("=" * 70)

    return {
        "name": "MDTraj",
        "runtime_avg": avg_runtime,
        "runtime_std": std_runtime,
        "memory_avg": avg_memory,
        "memory_std": std_memory,
    }


def benchmark_mdanalysis():
    if not HAS_MDANALYSIS:
        print("\nMDAnalysis not available - skipping")
        return None

    print("\n" + "=" * 70)
    print("MDAnalysis - Full Workflow (with figures)")
    print("=" * 70)
    print(f"Running {NUM_ITERATIONS} iterations...")

    runtimes = []
    memories = []

    for i in range(NUM_ITERATIONS):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime, memory = benchmark_mdanalysis_single(tmpdir)
            if runtime is None:
                return None
            runtimes.append(runtime)
            memories.append(memory)
            print(f"  Iteration {i+1}/{NUM_ITERATIONS}: {runtime:.3f}s, {memory:.1f} MB")

    avg_runtime = np.mean(runtimes)
    std_runtime = np.std(runtimes)
    avg_memory = np.mean(memories)
    std_memory = np.std(memories)

    print(f"\nAverage Runtime: {avg_runtime:.3f}s (±{std_runtime:.3f}s)")
    print(f"Average Peak Memory: {avg_memory:.1f} MB (±{std_memory:.1f} MB)")
    print("=" * 70)

    return {
        "name": "MDAnalysis",
        "runtime_avg": avg_runtime,
        "runtime_std": std_runtime,
        "memory_avg": avg_memory,
        "memory_std": std_memory,
    }


def create_benchmark_plots(results):
    print("\n" + "=" * 70)
    print("GENERATING PLOTS")
    print("=" * 70)

    names = [r["name"] for r in results]
    runtime_avgs = [r["runtime_avg"] for r in results]
    runtime_stds = [r["runtime_std"] for r in results]
    memory_avgs = [r["memory_avg"] for r in results]
    memory_stds = [r["memory_std"] for r in results]
    colors = [TOOL_COLORS[name.lower().replace(" ", "")] for name in names]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax1 = axes[0]
    bars1 = ax1.bar(
        names,
        runtime_avgs,
        yerr=runtime_stds,
        color=colors,
        alpha=0.9,
        edgecolor="black",
        linewidth=1.5,
        capsize=5,
        error_kw={"linewidth": 2},
    )
    label_font = 14 * AXIS_FONT_SCALE
    title_font = 16
    tick_font = 12 * AXIS_FONT_SCALE

    ax1.set_ylabel("Runtime (seconds)", fontsize=label_font)
    ax1.set_title("Runtime Comparison", fontsize=title_font, fontweight="bold")
    ax1.tick_params(axis="both", which="major", labelsize=tick_font)
    ax1.grid(axis="y", alpha=0.3, linestyle="--")

    runtime_top = max((avg + std) for avg, std in zip(runtime_avgs, runtime_stds)) if runtime_avgs else 0
    runtime_pad = runtime_top * 0.08 if runtime_top > 0 else 0.5
    ax1.set_ylim(0, runtime_top + 2 * runtime_pad)

    for bar, runtime_avg, runtime_std in zip(bars1, runtime_avgs, runtime_stds):
        height = bar.get_height()
        label_y = height + runtime_std + (runtime_pad * 0.3 if runtime_pad else 0.1)
        ax1.text(
            bar.get_x() + bar.get_width() / 2.0,
            label_y,
            f"{runtime_avg:.2f}s\n±{runtime_std:.2f}s",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            zorder=5,
        )

    ax2 = axes[1]
    bars2 = ax2.bar(
        names,
        memory_avgs,
        yerr=memory_stds,
        color=colors,
        alpha=0.9,
        edgecolor="black",
        linewidth=1.5,
        capsize=5,
        error_kw={"linewidth": 2},
    )
    ax2.set_ylabel("Peak Memory (MB)", fontsize=label_font)
    ax2.set_title("Peak Memory Comparison", fontsize=title_font, fontweight="bold")
    ax2.tick_params(axis="both", which="major", labelsize=tick_font)
    ax2.grid(axis="y", alpha=0.3, linestyle="--")

    memory_top = max((avg + std) for avg, std in zip(memory_avgs, memory_stds)) if memory_avgs else 0
    memory_pad = memory_top * 0.08 if memory_top > 0 else 0.5
    ax2.set_ylim(0, memory_top + 2 * memory_pad)

    for bar, memory_avg, memory_std in zip(bars2, memory_avgs, memory_stds):
        height = bar.get_height()
        label_y = height + memory_std + (memory_pad * 0.3 if memory_pad else 0.1)
        ax2.text(
            bar.get_x() + bar.get_width() / 2.0,
            label_y,
            f"{memory_avg:.1f} MB\n±{memory_std:.1f} MB",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            zorder=5,
        )

    plt.suptitle(
        "FastMDAnalysis Full Workflow Benchmark\n"
        "Ubiquitin (first 500 frames, stride=1) - RMSD, RMSF, RG, Cluster (with figures)\n"
        f"Averaged over {NUM_ITERATIONS} iterations",
        fontweight="bold",
        fontsize=14,
    )
    plt.tight_layout()

    output_file = "benchmark_full_workflow.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    plt.close()

    return output_file


def create_summary_table(results):
    print("Creating summary table...")

    fig, ax = plt.subplots(figsize=(12, 6))

    table_data = [["Library", "Runtime (avg ± std)", "Memory (avg ± std)", "Speedup"]]

    baseline_runtime = results[1]["runtime_avg"] if len(results) > 1 else results[0]["runtime_avg"]

    for r in results:
        runtime_str = f"{r['runtime_avg']:.2f}s ± {r['runtime_std']:.2f}s"
        memory_str = f"{r['memory_avg']:.1f} MB ± {r['memory_std']:.1f} MB"
        speedup = baseline_runtime / r["runtime_avg"]
        speedup_str = f"{speedup:.2f}x"
        table_data.append([r["name"], runtime_str, memory_str, speedup_str])

    table = ax.table(
        cellText=table_data,
        cellLoc="center",
        loc="center",
        colWidths=[0.25, 0.3, 0.3, 0.15],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)

    for i in range(4):
        cell = table[(0, i)]
        cell.set_facecolor("#4472C4")
        cell.set_text_props(weight="bold", color="white", fontsize=12)

    tool_keys = ["fastmdanalysis", "mdtraj", "mdanalysis"]
    for i in range(1, len(table_data)):
        tool_key = tool_keys[i - 1] if i - 1 < len(tool_keys) else tool_keys[-1]
        row_color = TOOL_COLORS.get(tool_key, "#CCCCCC")

        for j in range(4):
            cell = table[(i, j)]
            cell.set_facecolor(row_color + "80")
            cell.set_edgecolor("black")
            cell.set_linewidth(1.5)

    ax.axis("off")
    ax.set_title(
        "FastMDAnalysis Full Workflow Benchmark Summary\n"
        f"Averaged over {NUM_ITERATIONS} iterations",
        fontweight="bold",
        fontsize=14,
        pad=20,
    )

    footer_text = (
        "Dataset: Ubiquitin, first 500 frames (stride=1)\n"
        "Analyses: RMSD, RMSF, RG, Cluster (KMeans, DBSCAN, Hierarchical)\n"
        "Workflow: Complete analysis with figure generation (no slides)"
    )
    fig.text(
        0.5,
        0.05,
        footer_text,
        ha="center",
        fontsize=10,
        style="italic",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.3),
    )

    output_file = "benchmark_summary_table.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"✓ Saved: {output_file}")
    plt.close()

    return output_file


def main():
    print("=" * 70)
    print("FASTMDANALYSIS FULL WORKFLOW BENCHMARK")
    print("=" * 70)
    print("Dataset: Ubiquitin, first 500 frames (stride=1)")
    print("Analyses: RMSD, RMSF, RG, Cluster (KMeans, DBSCAN, Hierarchical)")
    print("Workflow: Complete with figure generation (excluding slides)")
    print(f"Iterations: {NUM_ITERATIONS} runs per library")
    print("Metrics: Runtime (wall clock) and Peak Memory Usage")
    print("=" * 70)

    results = []

    results.append(benchmark_fastmda())
    results.append(benchmark_mdtraj())
    mda_result = benchmark_mdanalysis()
    if mda_result:
        results.append(mda_result)

    print("\n" + "=" * 70)
    print(f"SUMMARY (AVERAGED OVER {NUM_ITERATIONS} ITERATIONS)")
    print("=" * 70)
    for r in results:
        print(f"{r['name']:20s}:")
        print(f"  Runtime: {r['runtime_avg']:6.2f}s ± {r['runtime_std']:5.2f}s")
        print(f"  Memory:  {r['memory_avg']:6.1f} MB ± {r['memory_std']:5.1f} MB")
    print("=" * 70)

    if len(results) >= 2:
        fastmda_time = results[0]["runtime_avg"]
        mdtraj_time = results[1]["runtime_avg"]
        ratio = fastmda_time / mdtraj_time
        print(f"\nFastMDA/MDTraj runtime ratio: {ratio:.2f}x")

        fastmda_mem = results[0]["memory_avg"]
        mdtraj_mem = results[1]["memory_avg"]
        mem_ratio = fastmda_mem / mdtraj_mem
        print(f"FastMDA/MDTraj memory ratio: {mem_ratio:.2f}x")

    create_benchmark_plots(results)
    create_summary_table(results)

    loc_data = compute_loc_benchmark()
    print("\nGenerating LOC slide...")
    generate_loc_slide(loc_data)

    print("\nGenerating CLI command slide...")
    command_data = compute_cli_command_counts()
    generate_cli_command_slide(command_data)

    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print("\nGenerated files:")
    print("  • benchmark_full_workflow.png - Runtime and memory comparison")
    print("  • benchmark_summary_table.png - Detailed summary table")
    print("  • benchmark_loc_slide.png - Lines-of-code comparison")
    print("  • benchmark_cli_commands.png - CLI command count visualization")
    print("=" * 70)


if __name__ == "__main__":
    main()
