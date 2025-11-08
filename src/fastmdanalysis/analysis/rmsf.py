# FastMDAnalysis/src/fastmdanalysis/analysis/rmsf.py
"""
RMSF Analysis Module

Calculates the Root-Mean-Square Fluctuation (RMSF) for each atom in an MD
trajectory. If an atom selection is provided, only those atoms are analyzed;
otherwise, all atoms are used. The analysis computes the fluctuations relative
to the average structure, saves the computed data, and automatically generates
a bar plot.

Plotting note:
- The x-axis shows ONLY atom indices (numeric), no residue/atom codes.
- Tick labels are auto-thinned to stay readable.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence, Union
from pathlib import Path
import logging

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .base import BaseAnalysis, AnalysisError

logger = logging.getLogger(__name__)


def _auto_tick_step(n: int, max_ticks: int) -> int:
    """Compute a thinning step so that ~max_ticks or fewer labels are drawn."""
    if n <= 0:
        return 1
    if n <= max_ticks:
        return 1
    return int(np.ceil(n / float(max_ticks)))


class RMSFAnalysis(BaseAnalysis):
    """
    Per-atom RMSF analysis with a readable, numeric-only x-axis (atom index).
    """

    def __init__(self, trajectory: md.Trajectory, atoms: Optional[str] = None, **kwargs):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            MDTraj trajectory to analyze.
        atoms : str or None
            MDTraj selection string (e.g., "protein and name CA"). If None, all atoms are used.
        kwargs : dict
            Passed through to BaseAnalysis (e.g., output).
        """
        super().__init__(trajectory, **kwargs)
        self.atoms: Optional[str] = atoms

        # Populated during run()
        self.data: Optional[np.ndarray] = None              # shape (N, 1)
        self.results: Dict[str, np.ndarray] = {}

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute RMSF (nm) for each selected atom relative to the average structure.

        Returns
        -------
        dict
            {"rmsf": (N, 1) array of per-atom RMSF values in nm}
        """
        try:
            # Atom selection (global indices)
            if self.atoms:
                sel = self.traj.topology.select(self.atoms)
                if sel is None or len(sel) == 0:
                    raise AnalysisError(f"No atoms selected using selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(sel)
            else:
                subtraj = self.traj

            # Average structure as reference
            avg_xyz = np.mean(subtraj.xyz, axis=0, keepdims=True)
            ref = md.Trajectory(avg_xyz, subtraj.topology)

            # Per-atom RMSF (nm) relative to average structure
            rmsf_values = md.rmsf(subtraj, ref)  # shape (N,)
            result_array = np.asarray(rmsf_values, dtype=float).reshape(-1, 1)

            if self.save_data:
                self._save_data(result_array, "rmsf")

            self.plot(result_array)

            if self.store_results:
                self.data = result_array
                self.results = {"rmsf": result_array}
                return self.results

            # Release references when not storing results
            self.data = None
            self.results = {}
            return {"rmsf": result_array}
        except AnalysisError:
            raise
        except Exception as e:
            raise AnalysisError(f"RMSF analysis failed: {e}")

    def plot(
        self,
        data: Optional[Union[Sequence[float], np.ndarray]] = None,
        *,
        max_ticks: int = 30,        # hard cap on number of labeled x ticks
        tick_step: Optional[int] = None,  # show every Nth tick (overrides max_ticks)
        rotate: int = 45,           # tick label rotation
        figsize=(12, 6),
        title: str = "RMSF per Atom",
        xlabel: str = "Atom Index",
        ylabel: str = "RMSF (nm)",
        color: Optional[str] = None,
    ) -> Path:
        """
        Generate a bar plot of RMSF with a numeric-only x-axis.

        Parameters
        ----------
        data
            RMSF values to plot. If None, uses computed data from run().
        max_ticks
            Target maximum number of x tick labels (used when tick_step is None).
        tick_step
            Force showing every Nth tick. If provided, overrides max_ticks heuristic.
        rotate
            Rotation angle for x-tick labels (degrees).
        figsize, title, xlabel, ylabel, color
            Usual matplotlib/IO controls.

        Returns
        -------
        Path
            File path of the saved plot image.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RMSF data available to plot. Run the analysis first.")

        y = np.asarray(data, dtype=float).flatten()
        n = int(y.size)
        x = np.arange(n)

        # Numeric-only atom indices for x-axis (1-based for readability)
        labels_all = [str(i + 1) for i in x]

        # Determine ticks
        step = tick_step if tick_step is not None else _auto_tick_step(n, max_ticks)
        ticks = x[::step]
        ticklabels = [labels_all[i] for i in ticks]

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        bar_kwargs = {"width": 0.9}
        if color is not None:
            bar_kwargs["color"] = color
        ax.bar(x, y, **bar_kwargs)

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels, rotation=rotate, ha="right")

        ax.grid(axis="y", alpha=0.3)
        ax.grid(axis="x", alpha=0.12)

        fig.tight_layout()
        outpath = self._save_plot(fig, "rmsf")
        plt.close(fig)
        return Path(outpath)

