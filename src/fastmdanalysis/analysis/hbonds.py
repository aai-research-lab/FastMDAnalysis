# FastMDAnalysis/src/fastmdanalysis/analysis/hbonds.py
"""
Hydrogen Bonds Analysis Module

Detects hydrogen bonds in an MD trajectory using the Baker–Hubbard algorithm.

Behavior
--------
- Optional atom selection via MDTraj DSL (e.g., "protein").
- Counts H-bonds **per frame** by running Baker–Hubbard on each frame.
- Saves:
    * hbonds_counts.dat  : (frame, n_hbonds)
    * hbonds.png         : line plot of H-bonds vs frame
- Returns in `results`:
    * "hbonds_counts": (T, 1) array (n_hbonds per frame)
    * "hbonds_per_frame": list of lists of (donor, hydrogen, acceptor) indices per frame
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Optional
import numpy as np
import mdtraj as md
import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from .base import BaseAnalysis, AnalysisError

logger = logging.getLogger(__name__)


class HBondsAnalysis(BaseAnalysis):
    def __init__(self, trajectory, atoms: Optional[str] = None, **kwargs):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            Trajectory to analyze.
        atoms : str or None
            MDTraj atom selection string to subset the trajectory. If None, use all atoms.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        super().__init__(trajectory, **kwargs)
        self.atoms = atoms
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, object] = {}

    def _subset_traj(self):
        """Return a trajectory possibly sliced by atom selection."""
        if self.atoms:
            sel = self.traj.topology.select(self.atoms)
            if sel is None or len(sel) == 0:
                raise AnalysisError(f"No atoms selected using: '{self.atoms}'")
            return self.traj.atom_slice(sel)
        return self.traj

    def run(self) -> Dict[str, object]:
        """
        Compute hydrogen bonds per frame using Baker–Hubbard.

        Returns
        -------
        dict
            {
              "hbonds_counts": (T, 1) array of per-frame counts,
              "hbonds_per_frame": list[ list[tuple(int,int,int)] ]
            }
        """
        try:
            subtraj = self._subset_traj()

            # Ensure standard bonds exist (required by some topologies for H-bond detection).
            try:
                subtraj.topology.create_standard_bonds()
            except Exception:
                # Not all topologies need this; ignore if unsupported.
                pass

            T = subtraj.n_frames
            counts = np.zeros(T, dtype=int)
            hbonds_per_frame: List[List[Tuple[int, int, int]]] = []

            # Robust approach: evaluate per frame (correct and fast for typical test-sized data).
            for i in range(T):
                hb = md.baker_hubbard(subtraj[i], periodic=False)
                hb_list = [(int(d), int(h), int(a)) for (d, h, a) in hb]
                hbonds_per_frame.append(hb_list)
                counts[i] = len(hb_list)

            # Store results
            self.data = counts.reshape(-1, 1)
            self.results = {
                "hbonds_counts": self.data,
                "hbonds_per_frame": hbonds_per_frame,
            }

            # Save data and plot
            # Two-column table: frame index, n_hbonds
            frames = np.arange(T, dtype=int).reshape(-1, 1)
            self._save_data(
                np.hstack([frames, self.data]),
                "hbonds_counts",
                header="frame n_hbonds",
                fmt="%d",
            )

            self.plot()  # ensure figure is produced by default
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            raise AnalysisError(f"Hydrogen bonds analysis failed: {e}")

    def plot(self, data: Optional[np.ndarray] = None, **kwargs):
        """
        Generate a plot of hydrogen bonds vs frame.

        Parameters
        ----------
        data : (T, 1) array-like, optional
            If None, uses data from `run()`.
        kwargs : dict
            Matplotlib options:
              - title (str): default "Hydrogen Bonds per Frame"
              - xlabel (str): default "Frame"
              - ylabel (str): default "Number of H-Bonds"
              - color (str): line/marker color
              - linestyle (str): default "-"
              - marker (str): default "o"

        Returns
        -------
        Path
            File path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No hydrogen bonds data to plot. Run the analysis first.")

        y = np.asarray(data).reshape(-1)
        x = np.arange(y.size)

        title = kwargs.get("title", "Hydrogen Bonds per Frame")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "Number of H-Bonds")
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
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        fig.tight_layout()
        outpath = self._save_plot(fig, "hbonds")
        plt.close(fig)
        return outpath
