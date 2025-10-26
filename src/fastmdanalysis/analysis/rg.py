# FastMDAnalysis/src/fastmdanalysis/analysis/rg.py
"""
Radius of Gyration (RG) Analysis Module

Calculates the radius of gyration for each frame of an MD trajectory.
Optionally accepts an MDTraj atom selection string so RG can be computed on a subset
of atoms. If no selection is provided, the calculation uses all atoms.

Outputs
-------
- rg.dat  : (N, 1) array of RG values per frame (nm)
- rg.png  : line plot of RG vs frame

Notes
-----
- Units are nanometers (nm), consistent with MDTraj.
"""

from __future__ import annotations

from typing import Optional, Dict
import logging
import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

logger = logging.getLogger(__name__)


class RGAnalysis(BaseAnalysis):
    def __init__(self, trajectory, atoms: Optional[str] = None, **kwargs):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            MDTraj atom selection string (e.g., "protein and name CA").
            If None, all atoms are used.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        super().__init__(trajectory, **kwargs)
        self.atoms = atoms
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, np.ndarray] = {}

    def _subset_traj(self):
        """Return trajectory sliced by atom selection if provided."""
        if self.atoms:
            sel = self.traj.topology.select(self.atoms)
            if sel is None or len(sel) == 0:
                raise AnalysisError(f"No atoms selected using: '{self.atoms}'")
            return self.traj.atom_slice(sel)
        return self.traj

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute the radius of gyration (nm) for each frame.

        Returns
        -------
        dict
            {"rg": (N, 1) array of RG values per frame in nm}
        """
        try:
            subtraj = self._subset_traj()
            logger.info(
                "RG: starting (atoms=%s, n_frames=%d, n_atoms=%d)",
                self.atoms if self.atoms else "ALL",
                subtraj.n_frames,
                subtraj.n_atoms,
            )

            rg_values = md.compute_rg(subtraj)  # shape (N,), units nm
            self.data = np.asarray(rg_values, dtype=float).reshape(-1, 1)
            self.results = {"rg": self.data}

            # Save data and default plot
            self._save_data(self.data, "rg", header="rg_nm", fmt="%.6f")
            self.plot()

            logger.info("RG: done.")
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            raise AnalysisError(f"Radius of gyration analysis failed: {e}")

    def plot(self, data: Optional[np.ndarray] = None, **kwargs):
        """
        Generate a plot of radius of gyration versus frame number.

        Parameters
        ----------
        data : array-like, optional
            RG data to plot; if None, uses data from `run()`.
        kwargs : dict
            Matplotlib options:
              - title (str): default "Radius of Gyration vs Frame"
              - xlabel (str): default "Frame"
              - ylabel (str): default "Radius of Gyration (nm)"
              - color (str): line/marker color
              - linestyle (str): default "-"
              - marker (str): default "o"

        Returns
        -------
        pathlib.Path
            Path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RG data available to plot. Run the analysis first.")

        y = np.asarray(data, dtype=float).reshape(-1)
        x = np.arange(y.size, dtype=int)

        title = kwargs.get("title", "Radius of Gyration vs Frame")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "Radius of Gyration (nm)")
        color = kwargs.get("color", None)
        linestyle = kwargs.get("linestyle", "-")
        marker = kwargs.get("marker", "o")

        fig, ax = plt.subplots(figsize=(10, 6))
        line_kwargs = {"linestyle": linestyle, "marker": marker}
        if color is not None:
            line_kwargs["color"] = color

        ax.plot(x, y, **line_kwargs)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        fig.tight_layout()

        out = self._save_plot(fig, "rg")
        plt.close(fig)
        return out