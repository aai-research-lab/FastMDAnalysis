# FastMDAnalysis/src/fastmdanalysis/analysis/rmsd.py
"""
RMSD Analysis Module

Calculates the Root-Mean-Square Deviation (RMSD) of an MD trajectory relative to a reference frame.
Optionally accepts an MDTraj atom selection and an `align` switch:
  - align=True  (default): classical RMSD with optimal superposition (Kabsch) via mdtraj.rmsd
  - align=False: no-fit RMSD (raw coordinate differences, no superposition)

Outputs
-------
- rmsd.dat : (N, 1) array of RMSD values per frame (nm)
- rmsd.png : line plot of RMSD vs frame
"""

from __future__ import annotations

from typing import Optional, Dict
import logging
import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

logger = logging.getLogger(__name__)


def _rmsd_no_fit(traj: md.Trajectory, ref: md.Trajectory, atom_indices=None) -> np.ndarray:
    """
    Compute RMSD per frame without optimal superposition (no-fit).

    Parameters
    ----------
    traj : md.Trajectory
        Trajectory with shape (T, A, 3)
    ref : md.Trajectory
        Single-frame reference with shape (1, A, 3)
    atom_indices : array-like or None
        Optional atom indices to select before computing.

    Returns
    -------
    np.ndarray shape (T,)
        RMSD in nm for each frame.
    """
    X = traj.xyz
    R = ref.xyz
    if atom_indices is not None:
        X = X[:, atom_indices, :]
        R = R[:, atom_indices, :]

    # no-fit RMSD = sqrt(mean(||x_i - y_i||^2)) over atoms and xyz
    diff = X - R  # (T, n, 3)
    # mean over atom and spatial dimensions; keep T
    msd = np.mean(np.sum(diff * diff, axis=2), axis=1)
    return np.sqrt(msd).astype(np.float64, copy=False)


class RMSDAnalysis(BaseAnalysis):
    def __init__(
        self,
        trajectory,
        reference_frame: int = 0,
        atoms: Optional[str] = None,
        align: bool = True,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        reference_frame : int
            Reference frame index (default 0). Negative indices allowed.
        atoms : str or None
            MDTraj atom selection string (e.g., "protein and name CA"). If None, all atoms are used.
        align : bool
            If True, compute classical RMSD with optimal superposition (mdtraj.rmsd).
            If False, compute no-fit RMSD (raw differences).
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        super().__init__(trajectory, **kwargs)
        self.reference_frame = 0 if reference_frame is None else int(reference_frame)
        self.atoms = atoms
        self.align = bool(align)
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, np.ndarray] = {}

    def _select_atoms(self) -> Optional[np.ndarray]:
        """Return atom indices for selection, or None for all atoms."""
        if self.atoms:
            sel = self.traj.topology.select(self.atoms)
            if sel is None or len(sel) == 0:
                raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
            return sel
        return None

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute RMSD for each frame relative to the reference frame.

        Returns
        -------
        dict
            {"rmsd": (N, 1) array of RMSD values in nm}
        """
        try:
            # Reference frame (mdtraj supports negative indices)
            try:
                ref = self.traj[self.reference_frame]
            except Exception as e:
                raise AnalysisError(f"Invalid reference frame index: {self.reference_frame}") from e

            atom_indices = self._select_atoms()

            logger.info(
                "RMSD: starting (ref=%d, atoms=%s, align=%s, n_frames=%d, n_atoms=%d)",
                self.reference_frame,
                self.atoms if self.atoms else "ALL",
                self.align,
                self.traj.n_frames,
                self.traj.n_atoms if atom_indices is None else len(atom_indices),
            )

            if self.align:
                # md.rmsd performs optimal superposition internally
                rmsd_values = md.rmsd(self.traj, ref, atom_indices=atom_indices)
            else:
                # No-fit RMSD
                rmsd_values = _rmsd_no_fit(self.traj, ref, atom_indices=atom_indices)

            result_array = np.asarray(rmsd_values, dtype=float).reshape(-1, 1)

            if self.save_data:
                self._save_data(result_array, "rmsd", header="rmsd_nm", fmt="%.6f")

            self.plot(result_array)
            logger.info("RMSD: done.")

            if self.store_results:
                self.data = result_array
                self.results = {"rmsd": result_array}
                return self.results

            # Release references when not storing results
            self.data = None
            self.results = {}
            return {"rmsd": result_array}

        except AnalysisError:
            raise
        except Exception as e:
            raise AnalysisError(f"RMSD analysis failed: {e}")

    def plot(self, data: Optional[np.ndarray] = None, **kwargs):
        """
        Generate a plot of RMSD versus frame number.

        Parameters
        ----------
        data : array-like, optional
            RMSD data to plot; if None, uses self.data.
        kwargs : dict
            Matplotlib options:
              - title (str): default "RMSD vs Frame (ref=<idx>, align=<bool>)"
              - xlabel (str): default "Frame"
              - ylabel (str): default "RMSD (nm)"
              - color (str): line color
              - linestyle (str): default "-"
              - marker (str): default "o"

        Returns
        -------
        Path
            Path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RMSD data available to plot. Please run the analysis first.")

        y = np.asarray(data, dtype=float).reshape(-1)
        x = np.arange(y.size, dtype=int)

        title = kwargs.get("title", f"RMSD vs Frame (ref={self.reference_frame}, align={self.align})")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "RMSD (nm)")
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

        out = self._save_plot(fig, "rmsd")
        plt.close(fig)
        return out

