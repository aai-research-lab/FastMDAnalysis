# FastMDAnalysis/src/fastmdanalysis/analysis/sasa.py
"""
SASA Analysis Module

Computes Solvent Accessible Surface Area (SASA) for an MD trajectory using
MDTraj's Shrake–Rupley algorithm.

Outputs
-------
Data tables (.dat):
  - total_sasa.dat              : (T, 1) total SASA per frame (nm^2)
  - residue_sasa.dat            : (T, R) per-residue SASA per frame (rows=frames, cols=residues; nm^2)
  - average_residue_sasa.dat    : (R, 1) mean SASA per residue across frames (nm^2)

Figures (.png):
  - total_sasa.png              : total SASA vs frame
  - residue_sasa.png            : heatmap (residue index × frame)
  - average_residue_sasa.png    : bar plot per residue
"""

from __future__ import annotations

from typing import Dict, Optional
import logging
import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

logger = logging.getLogger(__name__)


def _auto_tick_step(n: int, max_ticks: int) -> int:
    if n <= 0:
        return 1
    if n <= max_ticks:
        return 1
    return int(np.ceil(n / float(max_ticks)))


class SASAAnalysis(BaseAnalysis):
    def __init__(self, trajectory, probe_radius: float = 0.14, atoms: Optional[str] = None, **kwargs):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            Trajectory to analyze.
        probe_radius : float
            Probe radius in nm (default 0.14 nm).
        atoms : str or None
            MDTraj atom selection string. If None, uses all atoms.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        super().__init__(trajectory, **kwargs)
        self.probe_radius = float(probe_radius)
        self.atoms = atoms
        self.data: Optional[Dict[str, np.ndarray]] = None
        self.results: Dict[str, np.ndarray] = {}

    # --------------------------------------------------------------------- run

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute SASA datasets and generate default plots.

        Returns
        -------
        dict
            {
              "total_sasa": (T,),
              "residue_sasa": (T, R),
              "average_residue_sasa": (R,)
            }
        """
        try:
            # Subset trajectory by atom selection if provided
            if self.atoms:
                sel = self.traj.topology.select(self.atoms)
                if sel is None or len(sel) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(sel)
            else:
                subtraj = self.traj

            T = subtraj.n_frames
            logger.info(
                "SASA: starting (atoms=%s, n_frames=%d, n_atoms=%d, probe=%.3f nm)",
                self.atoms if self.atoms else "ALL",
                T, subtraj.n_atoms, self.probe_radius
            )

            # --- Compute per-residue SASA (robust to MDTraj versions)
            residue_sasa = None
            try:
                # Newer MDTraj versions support mode="residue"
                residue_sasa = md.shrake_rupley(subtraj, probe_radius=self.probe_radius, mode="residue")
                # shape (T, R)
            except TypeError:
                # Fallback: compute per-atom SASA then sum by residue
                atom_sasa = md.shrake_rupley(subtraj, probe_radius=self.probe_radius)  # (T, A)
                # Map atoms -> residue index (0..R-1 within subtraj topology)
                atom_res = np.array([a.residue.index for a in subtraj.topology.atoms], dtype=int)
                R = int(max(atom_res) + 1) if atom_res.size else 0
                residue_sasa = np.zeros((T, R), dtype=np.float32)
                for r in range(R):
                    residue_sasa[:, r] = atom_sasa[:, atom_res == r].sum(axis=1)

            if residue_sasa.ndim != 2:
                raise AnalysisError("Unexpected residue_sasa shape; expected 2D (T, R).")
            T2, R = residue_sasa.shape
            if T2 != T:
                raise AnalysisError("residue_sasa frame dimension mismatch.")

            total_sasa = residue_sasa.sum(axis=1)                  # (T,)
            average_residue_sasa = residue_sasa.mean(axis=0)       # (R,)

            self.data = {
                "total_sasa": total_sasa,
                "residue_sasa": residue_sasa,
                "average_residue_sasa": average_residue_sasa,
            }
            self.results = self.data

            # Save data
            self._save_data(total_sasa.reshape(-1, 1), "total_sasa", header="total_sasa_nm2", fmt="%.6f")
            # Rows=frames, Cols=residue index (1-based)
            self._save_data(
                residue_sasa,
                "residue_sasa",
                header="residue_sasa_nm2 (rows=frames, cols=residue index 1..R)",
                fmt="%.6f"
            )
            self._save_data(average_residue_sasa.reshape(-1, 1), "average_residue_sasa", header="avg_residue_sasa_nm2", fmt="%.6f")

            # Default plots
            self.plot()

            logger.info("SASA: done.")
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            raise AnalysisError(f"SASA analysis failed: {e}")

    # -------------------------------------------------------------------- plot

    def plot(self, data: Optional[Dict[str, np.ndarray]] = None, option: str = "all", **kwargs):
        """
        Replot SASA analysis outputs with customizable options.

        Parameters
        ----------
        data : dict or None
            If None, uses self.data.
        option : {'all','total','residue','average'}
            Which plots to generate.

        Returns
        -------
        dict | str
            Dict of paths for "all" else a single path.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No SASA data available. Run the analysis first.")

        plots = {}
        if option in ("all", "total"):
            plots["total"] = self._plot_total_sasa(data["total_sasa"], **kwargs)
        if option in ("all", "residue"):
            plots["residue"] = self._plot_residue_sasa(data["residue_sasa"], **kwargs)
        if option in ("all", "average"):
            plots["average"] = self._plot_average_residue_sasa(data["average_residue_sasa"], **kwargs)

        return plots if option == "all" else plots[option]

    def _plot_total_sasa(self, total_sasa: np.ndarray, **kwargs):
        """Total SASA vs frame."""
        x = np.arange(total_sasa.shape[0], dtype=int)
        title = kwargs.get("title_total", "Total SASA vs Frame")
        xlabel = kwargs.get("xlabel_total", "Frame")
        ylabel = kwargs.get("ylabel_total", "Total SASA (nm²)")
        color = kwargs.get("color_total", None)
        linestyle = kwargs.get("linestyle_total", "-")
        marker = kwargs.get("marker_total", "o")

        fig, ax = plt.subplots(figsize=(10, 6))
        line_kwargs = {"linestyle": linestyle, "marker": marker}
        if color is not None:
            line_kwargs["color"] = color

        ax.plot(x, total_sasa, **line_kwargs)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        out = self._save_plot(fig, "total_sasa")
        plt.close(fig)
        return out

    def _plot_residue_sasa(self, residue_sasa: np.ndarray, **kwargs):
        """Per-residue SASA heatmap (rows=residues, cols=frames)."""
        title = kwargs.get("title_residue", "Per-Residue SASA vs Frame")
        xlabel = kwargs.get("xlabel_residue", "Frame")
        ylabel = kwargs.get("ylabel_residue", "Residue Index")
        cmap = kwargs.get("cmap", "viridis")
        max_y_ticks = kwargs.get("max_y_ticks", 40)
        tick_step = kwargs.get("tick_step", None)  # overrides max_y_ticks if provided

        # Prepare data: (R, T), origin at lower so residue 1 is at bottom
        R = residue_sasa.shape[1]
        data = residue_sasa.T  # (R, T)

        fig, ax = plt.subplots(figsize=(12, 8))
        im = ax.imshow(data, aspect="auto", interpolation="none", cmap=cmap, origin="lower")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("SASA (nm²)")

        # Y ticks: 1-based residue indices with thinning
        step = tick_step if tick_step is not None else _auto_tick_step(R, max_y_ticks)
        ticks = np.arange(0, R, step, dtype=int)
        labels = (ticks + 1).astype(int).tolist()
        ax.set_yticks(ticks)
        ax.set_yticklabels([str(v) for v in labels])

        fig.tight_layout()
        out = self._save_plot(fig, "residue_sasa")
        plt.close(fig)
        return out

    def _plot_average_residue_sasa(self, average_sasa: np.ndarray, **kwargs):
        """Average per-residue SASA bar plot."""
        R = int(average_sasa.shape[0])
        x = np.arange(R, dtype=int)
        title = kwargs.get("title_avg", "Average per-Residue SASA")
        xlabel = kwargs.get("xlabel_avg", "Residue")
        ylabel = kwargs.get("ylabel_avg", "Average SASA (nm²)")
        color = kwargs.get("color_avg", None)
        max_x_ticks = kwargs.get("max_x_ticks", 40)
        tick_step = kwargs.get("tick_step_avg", None)

        fig, ax = plt.subplots(figsize=(12, 6))
        bar_kwargs = {}
        if color is not None:
            bar_kwargs["color"] = color
        ax.bar(x + 1, average_sasa.flatten(), **bar_kwargs)  # 1-based residue labels on axis

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.3)

        # X ticks thinning (1-based labels)
        step = tick_step if tick_step is not None else _auto_tick_step(R, max_x_ticks)
        ticks = np.arange(0, R, step, dtype=int)
        labels = (ticks + 1).astype(int).tolist()
        ax.set_xticks(ticks + 1)
        ax.set_xticklabels([str(v) for v in labels], rotation=45, ha="right")

        fig.tight_layout()
        out = self._save_plot(fig, "average_residue_sasa")
        plt.close(fig)
        return out
