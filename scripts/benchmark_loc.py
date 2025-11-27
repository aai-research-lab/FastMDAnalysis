#!/usr/bin/env python3
"""
Lines-of-code benchmark for the full workflow implementations.

Counts *statements* (AST-level) rather than physical lines, and:
  - Ignores MemoryMonitor for all tools (no helper, no instrumentation).
  - For MDTraj and MDAnalysis, drops any code involving DBSCAN or hierarchical
    clustering (the corresponding if-blocks).
  - For MDTraj and MDAnalysis, only counts helpers:
        _frame_range
        plot_kmeans_trajectory_artifacts
  - For FastMDAnalysis, counts only the workflow function.

The goal is to measure how many distinct statements a user must implement to
reproduce the FastMDAnalysis k-means workflow (RMSD, RMSF, RG, cluster) and its
artifacts.
"""

from __future__ import annotations

import ast
import inspect
import textwrap
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import benchmark_full_workflow as base


# ---------------------------------------------------------------------------
# Statement counting utilities
# ---------------------------------------------------------------------------

# Python 3.9-safe list of statement node types
_STMT_TYPES = (
    ast.FunctionDef,
    ast.AsyncFunctionDef,
    ast.For,
    ast.AsyncFor,
    ast.While,
    ast.If,
    ast.With,
    ast.AsyncWith,
    ast.Try,
    ast.Expr,
    ast.Assign,
    ast.AnnAssign,
    ast.AugAssign,
    ast.Return,
    ast.Raise,
    ast.Assert,
    ast.Import,
    ast.ImportFrom,
    ast.Global,
    ast.Nonlocal,
    ast.Pass,
    ast.Break,
    ast.Continue,
)


def _extract_function_source(fn: Callable) -> Tuple[str, ast.FunctionDef]:
    """Return dedented source for a function and its AST node."""
    source_lines, _ = inspect.getsourcelines(fn)
    dedented = textwrap.dedent("".join(source_lines))
    module = ast.parse(dedented)
    func_node = next(
        n for n in module.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
    )
    return dedented, func_node


def _stmt_source(block_lines: List[str], node: ast.stmt) -> str:
    start = node.lineno - 1
    end = getattr(node, "end_lineno", node.lineno) - 1
    return "\n".join(block_lines[start : end + 1])


def _test_contains_string(node: ast.AST, target: str) -> bool:
    """Return True if any Constant string in the subtree contains `target`."""
    for child in ast.walk(node):
        if isinstance(child, ast.Constant) and isinstance(child.value, str):
            if target in child.value:
                return True
    return False


def _should_skip_instrumentation(src: str) -> bool:
    """
    Skip instrumentation-related statements that are not logically required
    for the scientific workflow (MemoryMonitor, timing, explicit gc, etc.).
    """
    tokens = (
        "MemoryMonitor",
        "mem_monitor",
        "time.time(",
        "gc.collect",
        "peak_memory",
        "del ",
    )
    if any(tok in src for tok in tokens):
        return True
    # Skip the final "return runtime, peak_memory" / "return runtime, memory"
    if "return runtime, peak_memory" in src or "return runtime, memory" in src:
        return True
    return False


def _count_statements_in_function(
    fn: Callable,
    *,
    ignore_instrumentation: bool,
    ignore_dbscan: bool,
    ignore_hierarchical: bool,
) -> int:
    """
    Count AST statements inside a function, with options to ignore specific
    branches / instrumentation.

    - Multi-line calls (e.g., fastmda.cluster(...)) count as ONE statement.
    - If blocks for DBSCAN / hierarchical are dropped entirely when requested.
    - Instrumentation statements (MemoryMonitor, timing, gc, etc.) can be dropped.
    """
    dedented, func_node = _extract_function_source(fn)
    lines = dedented.splitlines()

    def count_stmt(node: ast.stmt) -> int:
        # For the outer function definition, do not count the def line itself;
        # only the body.
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            total = 0
            for child in node.body:
                if isinstance(child, _STMT_TYPES):
                    total += count_stmt(child)
            return total

        # For If statements, optionally drop DBSCAN / hierarchical branches.
        if isinstance(node, ast.If):
            if ignore_dbscan and _test_contains_string(node.test, "dbscan"):
                return 0
            if ignore_hierarchical and _test_contains_string(node.test, "hierarchical"):
                return 0

        # Instrumentation skipping (only applied in workflow functions).
        src = _stmt_source(lines, node)
        if ignore_instrumentation and _should_skip_instrumentation(src):
            return 0

        # Count this statement itself.
        total = 1

        # Recurse into compound statement bodies.
        if isinstance(
            node,
            (
                ast.If,
                ast.For,
                ast.AsyncFor,
                ast.While,
                ast.With,
                ast.AsyncWith,
                ast.Try,
            ),
        ):
            for child in getattr(node, "body", []):
                if isinstance(child, _STMT_TYPES):
                    total += count_stmt(child)
            for child in getattr(node, "orelse", []):
                if isinstance(child, _STMT_TYPES):
                    total += count_stmt(child)
            if isinstance(node, ast.Try):
                for handler in node.handlers:
                    for child in getattr(handler, "body", []):
                        if isinstance(child, _STMT_TYPES):
                            total += count_stmt(child)
                for child in getattr(node, "finalbody", []):
                    if isinstance(child, _STMT_TYPES):
                        total += count_stmt(child)

        return total

    return count_stmt(func_node)


# ---------------------------------------------------------------------------
# Library-specific configuration
# ---------------------------------------------------------------------------

def _fastmda_loc() -> Dict[str, object]:
    """LOC breakdown for FastMDAnalysis workflow."""
    workflow_fn = base.benchmark_fastmda_single

    workflow_loc = _count_statements_in_function(
        workflow_fn,
        ignore_instrumentation=True,
        ignore_dbscan=False,
        ignore_hierarchical=False,
    )

    helpers: List[Tuple[str, int]] = []  # MemoryMonitor explicitly ignored

    total = workflow_loc + sum(loc for _, loc in helpers)

    return {
        "name": "FastMDAnalysis",
        "workflow_loc": workflow_loc,
        "helpers": helpers,
        "total": total,
    }


def _mdtraj_loc() -> Dict[str, object]:
    """LOC breakdown for MDTraj workflow."""
    workflow_fn = base.benchmark_mdtraj_single

    workflow_loc = _count_statements_in_function(
        workflow_fn,
        ignore_instrumentation=True,
        ignore_dbscan=True,
        ignore_hierarchical=True,
    )

    helpers: List[Tuple[str, int]] = []

    # _frame_range helper
    helpers.append(
        (
            "_frame_range",
            _count_statements_in_function(
                base._frame_range,
                ignore_instrumentation=False,
                ignore_dbscan=False,
                ignore_hierarchical=False,
            ),
        )
    )

    # plot_kmeans_trajectory_artifacts helper
    helpers.append(
        (
            "plot_kmeans_trajectory_artifacts",
            _count_statements_in_function(
                base.plot_kmeans_trajectory_artifacts,
                ignore_instrumentation=False,
                ignore_dbscan=False,
                ignore_hierarchical=False,
            ),
        )
    )

    # MemoryMonitor + plot_cluster_summary intentionally NOT included.

    total = workflow_loc + sum(loc for _, loc in helpers)

    return {
        "name": "MDTraj",
        "workflow_loc": workflow_loc,
        "helpers": helpers,
        "total": total,
    }


def _mdanalysis_loc() -> Dict[str, object]:
    """LOC breakdown for MDAnalysis workflow."""
    if not getattr(base, "HAS_MDANALYSIS", False):
        return {
            "name": "MDAnalysis",
            "workflow_loc": 0,
            "helpers": [],
            "total": 0,
            "available": False,
        }

    workflow_fn = base.benchmark_mdanalysis_single

    workflow_loc = _count_statements_in_function(
        workflow_fn,
        ignore_instrumentation=True,
        ignore_dbscan=True,
        ignore_hierarchical=True,
    )

    helpers: List[Tuple[str, int]] = []

    helpers.append(
        (
            "_frame_range",
            _count_statements_in_function(
                base._frame_range,
                ignore_instrumentation=False,
                ignore_dbscan=False,
                ignore_hierarchical=False,
            ),
        )
    )

    helpers.append(
        (
            "plot_kmeans_trajectory_artifacts",
            _count_statements_in_function(
                base.plot_kmeans_trajectory_artifacts,
                ignore_instrumentation=False,
                ignore_dbscan=False,
                ignore_hierarchical=False,
            ),
        )
    )

    total = workflow_loc + sum(loc for _, loc in helpers)

    return {
        "name": "MDAnalysis",
        "workflow_loc": workflow_loc,
        "helpers": helpers,
        "total": total,
        "available": True,
    }


# ---------------------------------------------------------------------------
# Reporting + plotting
# ---------------------------------------------------------------------------

def _print_report(results: List[Dict[str, object]]) -> None:
    print("=" * 70)
    print("LINES-OF-CODE BENCHMARK (WORKFLOW + HELPERS)")
    print("=" * 70)
    print()

    for entry in results:
        name = entry["name"]
        workflow_loc = entry["workflow_loc"]
        helpers: List[Tuple[str, int]] = entry["helpers"]
        total = entry["total"]

        print(name)
        print("-" * len(name))
        print(f"  Workflow function: {workflow_loc} LOC")

        if helpers:
            print("  Helpers:")
            for helper_name, loc in helpers:
                print(f"    {helper_name}: {loc} LOC")
        else:
            print("  Helpers: (none)")

        print(f"  TOTAL: {total} LOC")
        print()


def _make_bar_plot(results: List[Dict[str, object]]) -> Path:
    names = [r["name"] for r in results]
    totals = [r["total"] for r in results]

    colors = [base.TOOL_COLORS.get(name.lower(), "#777777") for name in names]

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(names, totals, color=colors, edgecolor="black", linewidth=1.5)

    axis_font_scale = getattr(base, "AXIS_FONT_SCALE", 1.5)
    label_font = 14 * axis_font_scale
    title_font = 12
    tick_font = 12 * axis_font_scale

    ax.set_ylabel("Effective Lines of Code", fontsize=label_font)
    ax.set_title(
        "Statements Required Per Workflow (RMSD, RMSF, RG, KMeans)",
        fontsize=title_font,
        fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    ax.tick_params(axis="x", labelsize=tick_font)
    ax.tick_params(axis="y", labelsize=tick_font)

    if totals:
        top = max(totals)
        pad = max(1, int(top * 0.08))
        ax.set_ylim(0, top + 3 * pad)

    for bar, total in zip(bars, totals):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + (0.3 if total < 10 else pad),
            f"{total} LOC",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            zorder=5,
        )

    plt.tight_layout()
    outdir = Path("benchmark_output")
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "benchmark_loc.png"
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"\n✓ Saved LOC plot -> {outpath}")
    return outpath


def main() -> None:
    results: List[Dict[str, object]] = []

    fastmda = _fastmda_loc()
    results.append(fastmda)

    mdtraj = _mdtraj_loc()
    results.append(mdtraj)

    mda = _mdanalysis_loc()
    if mda.get("available", True):
        results.append(mda)

    _print_report(results)
    _make_bar_plot(results)


if __name__ == "__main__":
    main()
