"""Energy CSV diagnostics and simulation-health plots."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from fastmdxplora.analysis.plotting import new_figure, save_figure


def _normalize_column_name(name: object) -> str:
    text = str(name).strip()
    text = text.lstrip("#").strip()
    text = text.strip('"').strip("'").strip()
    return " ".join(text.split()).lower()


def _find_column(df: pd.DataFrame, *needles: str) -> str | None:
    normalized = {_normalize_column_name(col): str(col) for col in df.columns}
    for norm, original in normalized.items():
        if all(needle.lower() in norm for needle in needles):
            return original
    return None


def _column_map(df: pd.DataFrame) -> dict[str, str | None]:
    return {
        "step": _find_column(df, "step"),
        "time": _find_column(df, "time"),
        "potential_energy": _find_column(df, "potential", "energy"),
        "kinetic_energy": _find_column(df, "kinetic", "energy"),
        "total_energy": _find_column(df, "total", "energy"),
        "temperature": _find_column(df, "temperature"),
        "volume": _find_column(df, "volume"),
        "density": _find_column(df, "density"),
        "speed": _find_column(df, "speed"),
        "progress": _find_column(df, "progress"),
    }


def read_energy_csv(path: str | Path) -> pd.DataFrame:
    """Read an OpenMM/FastMDXplora energy.csv file."""
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"energy.csv not found: {csv_path}")

    try:
        df = pd.read_csv(csv_path)
    except EmptyDataError as exc:
        raise ValueError(f"energy.csv contains no columns: {csv_path}") from exc

    df.columns = [
        str(col).lstrip("#").strip().strip('"').strip("'").strip()
        for col in df.columns
    ]

    for col in df.columns:
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.notna().any():
            df[col] = converted

    if df.empty:
        raise ValueError(f"energy.csv contains no rows: {csv_path}")

    return df


def _x_axis(df: pd.DataFrame, columns: dict[str, str | None]) -> tuple[np.ndarray, str]:
    time_col = columns.get("time")
    step_col = columns.get("step")

    if time_col and time_col in df:
        time = pd.to_numeric(df[time_col], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(time).any():
            max_time = np.nanmax(time)
            if max_time >= 1000:
                return time / 1000.0, "Time (ns)"
            return time, "Time (ps)"

    if step_col and step_col in df:
        step = pd.to_numeric(df[step_col], errors="coerce").to_numpy(dtype=float)
        return step, "Step"

    return np.arange(len(df), dtype=float), "Frame"


def _plot_available_series(
    df: pd.DataFrame,
    *,
    series: list[tuple[str, str | None]],
    title: str,
    ylabel: str,
    output_path: Path,
) -> Path | None:
    columns = _column_map(df)
    x, xlabel = _x_axis(df, columns)

    available = [(label, col) for label, col in series if col and col in df.columns]
    if not available:
        return None

    fig, ax = new_figure(figsize=(7.4, 4.6), title=title)
    marker = "o" if len(x) <= 25 else None

    for label, col in available:
        y = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        ax.plot(x, y, marker=marker, linewidth=1.8, markersize=4, label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    ax.legend()

    if len(x) == 1:
        ax.set_xlim(float(x[0]) - 1.0, float(x[0]) + 1.0)

    return save_figure(fig, output_path)


def plot_energy_trace(
    energy_csv: str | Path,
    output_path: str | Path | None = None,
) -> Path | None:
    """Plot potential, kinetic, and total energy from energy.csv."""
    csv_path = Path(energy_csv)
    df = read_energy_csv(csv_path)
    columns = _column_map(df)

    output = Path(output_path) if output_path else csv_path.with_name("energy_trace.png")
    return _plot_available_series(
        df,
        series=[
            ("Potential energy", columns["potential_energy"]),
            ("Kinetic energy", columns["kinetic_energy"]),
            ("Total energy", columns["total_energy"]),
        ],
        title="Simulation Energy Trace",
        ylabel="Energy (kJ/mol)",
        output_path=output,
    )


def plot_simulation_health(
    energy_csv: str | Path,
    output_path: str | Path | None = None,
) -> Path | None:
    """Plot temperature, density, and volume from energy.csv."""
    csv_path = Path(energy_csv)
    df = read_energy_csv(csv_path)
    columns = _column_map(df)

    output = Path(output_path) if output_path else csv_path.with_name("simulation_health.png")
    return _plot_available_series(
        df,
        series=[
            ("Temperature", columns["temperature"]),
            ("Density", columns["density"]),
            ("Volume", columns["volume"]),
        ],
        title="Simulation Health Trace",
        ylabel="Reported value",
        output_path=output,
    )


def summarize_energy_csv(
    energy_csv: str | Path,
    output_path: str | Path | None = None,
) -> pd.DataFrame:
    """Write summary statistics for numeric columns in energy.csv."""
    csv_path = Path(energy_csv)
    df = read_energy_csv(csv_path)

    rows: list[dict[str, Any]] = []
    for col in df.columns:
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        if values.empty:
            continue

        rows.append(
            {
                "metric": col,
                "count": int(values.count()),
                "first": float(values.iloc[0]),
                "last": float(values.iloc[-1]),
                "drift": float(values.iloc[-1] - values.iloc[0]),
                "mean": float(values.mean()),
                "std": float(values.std(ddof=0)),
                "min": float(values.min()),
                "max": float(values.max()),
            }
        )

    summary = pd.DataFrame(rows)
    output = Path(output_path) if output_path else csv_path.with_name("energy_summary.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, index=False)
    return summary


def ensure_energy_report_assets(
    energy_csv: str | Path,
    output_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Generate energy diagnostic assets if energy.csv exists."""
    csv_path = Path(energy_csv)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return {}

    out_dir = Path(output_dir) if output_dir else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    assets: dict[str, Path] = {}

    energy_trace = plot_energy_trace(csv_path, out_dir / "energy_trace.png")
    if energy_trace is not None:
        assets["energy_trace"] = energy_trace

    health_trace = plot_simulation_health(csv_path, out_dir / "simulation_health.png")
    if health_trace is not None:
        assets["simulation_health"] = health_trace

    summarize_energy_csv(csv_path, out_dir / "energy_summary.csv")
    assets["energy_summary"] = out_dir / "energy_summary.csv"

    return assets
