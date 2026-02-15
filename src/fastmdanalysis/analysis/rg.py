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
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style

logger = logging.getLogger(__name__)


class RGAnalysis(BaseAnalysis):
    _ALIASES = {
        "atom_indices": "atoms",
        "selection": "atoms",
    }
    
    def __init__(
        self, 
        trajectory, 
        atoms: Optional[str] = None, 
        mass_weighted: bool = False,
        compute_stat: bool = False,
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            MDTraj atom selection string (e.g., "protein and name CA").
            If None, all atoms are used.
            Aliases: atom_indices, selection
        mass_weighted : bool
            If True and MDTraj supports it, compute mass-weighted RG.
            Note: MDTraj's compute_rg is already mass-weighted by default.
        strict : bool
            If True, raise errors for unknown options. If False, log warnings.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "atoms": atoms,
            "mass_weighted": mass_weighted,
            "compute_stat": compute_stat,
            "strict": strict,
        }
        analysis_opts.update(kwargs)
        
        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"atoms", "mass_weighted", "compute_stat", "strict", "output"},
            context="rg",
            warn=warn_unknown,
        )

        atoms = resolved.get("atoms", None)
        mass_weighted = resolved.get("mass_weighted", False)
        compute_stat = resolved.get("compute_stat", False)
        base_kwargs = {
            k: v
            for k, v in resolved.items()
            if k not in ("atoms", "mass_weighted", "compute_stat", "strict")
        }

        super().__init__(trajectory, **base_kwargs)
        self.atoms = atoms
        self.mass_weighted = mass_weighted
        self.compute_stat = bool(compute_stat)
        self.strict = strict
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, np.ndarray] = {}

        logger.info("Initialized RG analysis: atoms=%s, mass_weighted=%s",
                   self.atoms if self.atoms else "ALL", self.mass_weighted)

    def _subset_traj(self):
        """Return trajectory sliced by atom selection if provided."""
        if self.atoms:
            logger.debug("Selecting atoms: %s", self.atoms)
            sel = self.traj.topology.select(self.atoms)
            if sel is None or len(sel) == 0:
                raise AnalysisError(f"No atoms selected using: '{self.atoms}'")
            subtraj = self.traj.atom_slice(sel)
            logger.debug("Atom selection yielded %d atoms", len(sel))
            return subtraj
        logger.debug("Using all %d atoms", self.traj.n_atoms)
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
                "Starting RG calculation: atoms=%s, n_frames=%d, n_atoms=%d",
                self.atoms if self.atoms else "ALL",
                subtraj.n_frames,
                subtraj.n_atoms,
            )

            logger.debug("Computing radius of gyration...")
            rg_values = md.compute_rg(subtraj)  # shape (N,), units nm
            self.data = np.asarray(rg_values, dtype=float).reshape(-1, 1)
            self.results = {"rg": self.data}

            if self.compute_stat:
                mean_val = float(np.nanmean(rg_values))
                std_val = float(np.nanstd(rg_values))
                self.results["rg_stats"] = {"mean": mean_val, "std": std_val}
                self._save_data(
                    np.array([[mean_val, std_val]], dtype=float),
                    "rg_stats",
                    header="mean_nm std_nm",
                    fmt="%.6f",
                )

            # Save data and default plot
            logger.info("Saving RG data...")
            self._save_data(self.data, "rg", header="rg_nm", fmt="%.6f")
            
            logger.info("Generating RG plot...")
            self.plot()

            rg_range = (self.data.min(), self.data.max())
            logger.info("RG analysis complete: range [%.3f, %.3f] nm", rg_range[0], rg_range[1])
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            logger.exception("Radius of gyration analysis failed")
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

        logger.debug("Generating RG plot")
        y = np.asarray(data, dtype=float).reshape(-1)
        x = np.arange(1, y.size + 1, dtype=int)

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
        if self.compute_stat:
            mean_val = float(np.nanmean(y))
            std_val = float(np.nanstd(y))
            ax.axhline(mean_val, color="black", linestyle="--", linewidth=1.2, label="mean")
            ax.fill_between(
                [x.min(), x.max()],
                mean_val - std_val,
                mean_val + std_val,
                color="gray",
                alpha=0.2,
                label="±1 std",
            )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        if self.compute_stat:
            ax.legend(loc="best")

        apply_slide_style(
            ax,
            x_values=x,
            y_values=y,
            x_max_ticks=8,
            y_max_ticks=6,
            zero_x=True,
            zero_y=True,
        )
        fig.tight_layout()

        out = self._save_plot(fig, "rg")
        plt.close(fig)
        logger.debug("RG plot saved to %s", out)
        return out

    def _save_plot(self, fig, name: str):
        """Save the figure as a PNG file in the output directory and log its path."""
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        logger.info("Plot saved to %s", plot_path)
        return plot_path
