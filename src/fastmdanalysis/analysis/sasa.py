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
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style, auto_ticks

logger = logging.getLogger(__name__)


class SASAAnalysis(BaseAnalysis):
    _ALIASES = {
        "probe_radius_nm": "probe_radius",
        "atom_indices": "atoms",
        "selection": "atoms",
    }
    
    def __init__(
        self, 
        trajectory, 
        probe_radius: float = 0.14, 
        atoms: Optional[str] = None,
        n_sphere_points: Optional[int] = None,
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            Trajectory to analyze.
        probe_radius : float
            Probe radius in nm (default 0.14 nm). Alias: probe_radius_nm
        atoms : str or None
            MDTraj atom selection string. If None, uses all atoms.
            Aliases: atom_indices, selection
        n_sphere_points : int, optional
            Number of sphere points for SASA calculation (default uses MDTraj default).
        strict : bool
            If True, raise errors for unknown options. If False, log warnings.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "probe_radius": probe_radius,
            "atoms": atoms,
            "n_sphere_points": n_sphere_points,
            "strict": strict,
        }
        analysis_opts.update(kwargs)
        
        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"probe_radius", "atoms", "n_sphere_points", "strict", "output"},
            context="sasa",
            warn=warn_unknown,
        )

        probe_radius = resolved.get("probe_radius", 0.14)
        atoms = resolved.get("atoms", None)
        n_sphere_points = resolved.get("n_sphere_points", None)
        base_kwargs = {k: v for k, v in resolved.items() 
                      if k not in ("probe_radius", "atoms", "n_sphere_points", "strict")}

        super().__init__(trajectory, **base_kwargs)
        self.probe_radius = float(probe_radius)
        self.atoms = atoms
        self.n_sphere_points = int(n_sphere_points) if n_sphere_points is not None else None
        self.strict = strict
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

            residue_sasa = None
            sasa_kwargs = {"probe_radius": self.probe_radius}
            if self.n_sphere_points is not None:
                sasa_kwargs["n_sphere_points"] = self.n_sphere_points
            
            try:
                # Newer MDTraj versions support mode="residue"
                residue_sasa = md.shrake_rupley(subtraj, mode="residue", **sasa_kwargs)
                # shape (T, R)
            except TypeError:
                # Fallback: compute per-atom SASA then sum by residue
                atom_sasa = md.shrake_rupley(subtraj, **sasa_kwargs)  # (T, A)
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
        apply_slide_style(
            ax,
            x_values=x,
            y_values=total_sasa,
            integer_x=True,
            x_max_ticks=10,
            y_max_ticks=6,
             zero_x=True,
             zero_y=True,
        )
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

        if tick_step is not None:
            step = max(1, int(tick_step))
            res_ticks = np.arange(0, R, step, dtype=int)
        else:
            auto = auto_ticks(np.arange(R, dtype=int), max_ticks=max_y_ticks, integer=True)
            res_ticks = auto.astype(int) if auto is not None and auto.size > 0 else np.arange(R, dtype=int)

        apply_slide_style(
            ax,
            x_values=np.arange(residue_sasa.shape[0], dtype=int),
            y_ticks=res_ticks,
            integer_x=True,
            integer_y=True,
            zero_x=True,
            zero_y=True,
        )
        tick_font = ax.get_yticklabels()[0].get_fontsize() if ax.get_yticklabels() else None
        ax.set_yticks(res_ticks)
        ax.set_yticklabels(
            [str(int(v) + 1) for v in res_ticks],
            rotation_mode="anchor",
            fontsize=tick_font,
        )

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

        if tick_step is not None:
            step = max(1, int(tick_step))
            ticks = np.arange(1, R + 1, step, dtype=int)
        else:
            auto = auto_ticks(np.arange(1, R + 1, dtype=int), max_ticks=max_x_ticks, integer=True)
            ticks = auto.astype(int) if auto is not None and auto.size > 0 else np.arange(1, R + 1, dtype=int)

        apply_slide_style(
            ax,
            x_ticks=ticks,
            y_values=average_sasa.flatten(),
            integer_x=True,
            zero_y=True,
            x_tick_rotation=45,
        )
        tick_font = ax.get_xticklabels()[0].get_fontsize() if ax.get_xticklabels() else None
        ax.set_xticks(ticks)
        ax.set_xticklabels(
            [str(int(v)) for v in ticks],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=tick_font,
        )

        fig.tight_layout()
        out = self._save_plot(fig, "average_residue_sasa")
        plt.close(fig)
        return out

