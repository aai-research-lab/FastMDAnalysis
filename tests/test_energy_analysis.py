from __future__ import annotations

from pathlib import Path

import pandas as pd

from fastmdxplora.analysis.energy import (
    ensure_energy_report_assets,
    plot_energy_trace,
    plot_simulation_health,
    read_energy_csv,
    summarize_energy_csv,
)


def _write_energy_csv(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                '"#Step","Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"',
                "500,0.25,-1000.0,250.0,-750.0,299.0,10.0,0.95",
                "1000,0.50,-1002.0,252.0,-750.0,300.0,10.1,0.96",
                "1500,0.75,-1003.0,253.0,-750.0,301.0,10.2,0.97",
            ]
        ),
        encoding="utf-8",
    )


def test_read_energy_csv_cleans_openmm_header(tmp_path: Path) -> None:
    csv_path = tmp_path / "energy.csv"
    _write_energy_csv(csv_path)

    df = read_energy_csv(csv_path)

    assert not df.empty
    assert "Step" in df.columns
    assert "Potential Energy (kJ/mole)" in df.columns


def test_energy_plots_and_summary_are_written(tmp_path: Path) -> None:
    csv_path = tmp_path / "energy.csv"
    _write_energy_csv(csv_path)

    energy_plot = plot_energy_trace(csv_path)
    health_plot = plot_simulation_health(csv_path)
    summary = summarize_energy_csv(csv_path)

    assert energy_plot is not None
    assert health_plot is not None
    assert energy_plot.exists()
    assert health_plot.exists()
    assert (tmp_path / "energy_summary.csv").exists()
    assert isinstance(summary, pd.DataFrame)
    assert "metric" in summary.columns


def test_ensure_energy_report_assets(tmp_path: Path) -> None:
    csv_path = tmp_path / "energy.csv"
    _write_energy_csv(csv_path)

    assets = ensure_energy_report_assets(csv_path)

    assert assets["energy_trace"].exists()
    assert assets["simulation_health"].exists()
    assert assets["energy_summary"].exists()
