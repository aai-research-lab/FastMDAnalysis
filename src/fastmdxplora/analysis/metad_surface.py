"""The free energy surface from a metadynamics run.

Written to `metadynamics_surface.json` and stopped there: no figure, no entry
in the analysis manifest, no mention in the report. Ten analyses of the
trajectory each produced a curve and a plot, and the one result the run
existed to produce did not -- the same gap the umbrella study had.

Provisional surfaces are drawn too, and said to be provisional. A run whose
bias has not settled still describes the landscape it has filled so far, and
the refusal beside it explains what is missing; withholding the picture as
well leaves a reader with a sentence where they could have had both.

The trajectory of a metadynamics run is not a Boltzmann ensemble -- that is
the point of the method -- so this reads what the simulation phase computed
from the hills rather than analysing the frames.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis


class MetadynamicsSurface(Analysis):
    """The free energy along a metadynamics run's collective variable."""

    name = "metad_surface"
    description = "Metadynamics free energy surface"
    #: The coordinate was chosen when the bias was planned; a selection here
    #: would describe a different question.
    honours_selection = False
    default_selection = None
    time_series = False
    #: Only where a metadynamics run produced a record.
    requires_metadynamics = True

    def _record_path(self) -> Path | None:
        """Where the simulation phase left the surface it computed."""
        # `self.output_dir` already has this analysis's own name appended,
        # so the run's directory is two levels up: <run>/analysis/<name>.
        output = Path(self.output_dir)
        for candidate in (
            output.parent.parent / "simulation" / "metadynamics_surface.json",
            output.parent.parent / "metadynamics_surface.json",
            output.parent / "simulation" / "metadynamics_surface.json",
            output.parent / "metadynamics_surface.json",
        ):
            if candidate.is_file():
                return candidate
        return None

    def compute(self, traj: Any) -> dict[str, Any]:
        path = self._record_path()
        if path is None:
            raise FileNotFoundError(
                "No metadynamics_surface.json beside this run, so there is "
                "no surface to draw. This analysis reports what the "
                "simulation phase computed from the hills; it does not "
                "recompute it."
            )
        record = json.loads(path.read_text(encoding="utf-8"))
        grid = np.asarray(record.get("grid") or [], dtype=float)
        energy = np.asarray(
            [np.nan if value is None else float(value)
             for value in (record.get("free_energy_kjmol") or [])],
            dtype=float)

        self._refused = record.get("refused")
        self._provisional = bool(record.get("provisional"))
        self._evidence = record.get("evidence") or {}
        if not len(grid) or not len(energy):
            raise ValueError(
                "The metadynamics record holds no surface: "
                + (str(self._refused) if self._refused
                   else "the coordinate did not move, so there is nothing "
                        "along it to draw.")
            )
        return {"coordinate": grid, "free_energy_kjmol": energy}

    def plot(self, result: dict[str, Any], ax: plt.Axes) -> None:
        coordinate = result["coordinate"]
        energy = result["free_energy_kjmol"]
        ax.plot(coordinate, energy, linewidth=1.6)

        sampled = np.isfinite(energy)
        if sampled.any():
            lowest = int(np.nanargmin(np.where(sampled, energy, np.inf)))
            ax.axvline(coordinate[lowest], color="#888888", linestyle=":",
                       linewidth=1.0,
                       label=f"minimum at {coordinate[lowest]:.3g}")

        # Where the drift was judged, because the surface above it is read
        # off the shape rather than trusted point by point: those regions are
        # visited rarely and their estimate moves by several kJ/mol however
        # long the run.
        ceiling = self._evidence.get("drift_ceiling_kjmol")
        if ceiling and sampled.any():
            floor = float(np.nanmin(energy))
            ax.axhline(floor + float(ceiling), color="#cccccc",
                       linestyle="--", linewidth=0.8,
                       label=f"judged below {ceiling:g} kJ/mol")
        if sampled.any():
            ax.legend(loc="best")

        if self._provisional:
            ax.set_title(f"{self.figure_title()} (provisional)")

    def default_ylabel(self) -> str:
        return "Free energy (kJ/mol)"

    def default_xlabel(self) -> str:
        return "Collective variable"

    def save_data(self, result: dict[str, Any], path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        header = "# coordinate free_energy_kjmol"
        if self._provisional:
            header += "\n# provisional: the bias had not settled when this was cut"
        lines = [header]
        for x, y in zip(result["coordinate"], result["free_energy_kjmol"]):
            lines.append(f"{x:.6f} {'nan' if np.isnan(y) else f'{y:.6f}'}")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path


register_analysis(MetadynamicsSurface.name, MetadynamicsSurface)
