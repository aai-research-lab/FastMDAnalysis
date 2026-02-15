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
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style

logger = logging.getLogger(__name__)


class RMSFAnalysis(BaseAnalysis):
    """
    Per-atom RMSF analysis with a readable, numeric-only x-axis (atom index).
    """

    _ALIASES = {
        "atom_indices": "atoms",
        "selection": "atoms",
        "reference": "reference_frame",
    }

    def __init__(
        self, 
        trajectory: md.Trajectory, 
        atoms: Optional[str] = None, 
        per_residue: bool = False,
        compute_stat: bool = False,
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            MDTraj trajectory to analyze.
        atoms : str or None
            MDTraj selection string (e.g., "protein and name CA"). If None, all atoms are used.
            Aliases: atom_indices, selection
        per_residue : bool
            If True, aggregate per-atom RMSF to per-residue (mean over atoms in each residue).
        strict : bool
            If True, raise errors for unknown options. If False, log warnings.
        kwargs : dict
            Passed through to BaseAnalysis (e.g., output).
        """
        logger.info("Initializing RMSF analysis")
        logger.debug("Input parameters: atoms=%s, per_residue=%s, strict=%s", 
                    atoms, per_residue, strict)
        
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "atoms": atoms,
            "per_residue": per_residue,
            "compute_stat": compute_stat,
            "strict": strict,
        }
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"atoms", "per_residue", "compute_stat", "strict", "output", "reference_frame"},
            context="rmsf",
            warn=warn_unknown,
        )

        atoms = resolved.get("atoms", None)
        per_residue = resolved.get("per_residue", False)
        compute_stat = resolved.get("compute_stat", False)
        base_kwargs = {k: v for k, v in resolved.items() 
                      if k not in ("atoms", "per_residue", "compute_stat", "strict", "reference_frame")}

        super().__init__(trajectory, **base_kwargs)
        self.atoms: Optional[str] = atoms
        self.per_residue: bool = bool(per_residue)
        self.compute_stat: bool = bool(compute_stat)
        self.strict = strict

        # Populated during run()
        self.data: Optional[np.ndarray] = None              # shape (N, 1)
        self.results: Dict[str, np.ndarray] = {}
        
        logger.info("RMSF analysis initialized with %d frames, %d atoms", 
                   trajectory.n_frames, trajectory.n_atoms)
        if atoms:
            logger.info("Atom selection: %s", atoms)
        if per_residue:
            logger.info("Per-residue aggregation enabled")

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute RMSF (nm) for each selected atom relative to the average structure.

        Returns
        -------
        dict
            {"rmsf": (N, 1) array of per-atom RMSF values in nm}
            If per_residue=True, also includes {"rmsf_per_residue": (R, 1) array}
        """
        logger.info("Starting RMSF analysis")
        try:
            # Atom selection (global indices)
            if self.atoms:
                logger.debug("Selecting atoms with: %s", self.atoms)
                sel = self.traj.topology.select(self.atoms)
                if sel is None or len(sel) == 0:
                    raise AnalysisError(f"No atoms selected using selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(sel)
                logger.info("Selected %d atoms from trajectory", subtraj.n_atoms)
            else:
                subtraj = self.traj
                logger.debug("Using all %d atoms", subtraj.n_atoms)

            # Average structure as reference
            logger.debug("Computing average structure from %d frames", subtraj.n_frames)
            avg_xyz = np.mean(subtraj.xyz, axis=0, keepdims=True)
            ref = md.Trajectory(avg_xyz, subtraj.topology)

            # Per-atom RMSF (nm) relative to average structure
            logger.info("Computing RMSF values")
            rmsf_values = md.rmsf(subtraj, ref)  # shape (N,)
            self.data = np.asarray(rmsf_values, dtype=float).reshape(-1, 1)
            self.results = {"rmsf": self.data}

            if self.compute_stat:
                mean_val = float(np.nanmean(rmsf_values))
                std_val = float(np.nanstd(rmsf_values))
                self.results["rmsf_stats"] = {"mean": mean_val, "std": std_val}
                self._save_data(
                    np.array([[mean_val, std_val]], dtype=float),
                    "rmsf_stats",
                    header="mean_nm std_nm",
                    fmt="%.6f",
                )
            
            logger.info("RMSF computation completed - mean: %.4f nm, std: %.4f nm, range: [%.4f, %.4f] nm",
                       np.mean(rmsf_values), np.std(rmsf_values), 
                       np.min(rmsf_values), np.max(rmsf_values))

            if self.per_residue:
                logger.debug("Computing per-residue RMSF aggregation")
                atom_to_residue = np.array([a.residue.index for a in subtraj.topology.atoms], dtype=int)
                n_residues = int(max(atom_to_residue) + 1) if atom_to_residue.size else 0
                per_residue_rmsf = np.zeros(n_residues, dtype=float)
                for r in range(n_residues):
                    mask = atom_to_residue == r
                    if np.any(mask):
                        per_residue_rmsf[r] = np.mean(rmsf_values[mask])

                self.results["rmsf_per_residue"] = per_residue_rmsf.reshape(-1, 1)
                self._save_data(
                    per_residue_rmsf.reshape(-1, 1),
                    "rmsf_per_residue",
                    header="rmsf_per_residue_nm",
                    fmt="%.6f",
                )
                logger.info("Per-residue RMSF computed for %d residues - mean: %.4f nm, range: [%.4f, %.4f] nm",
                           n_residues, np.mean(per_residue_rmsf), 
                           np.min(per_residue_rmsf), np.max(per_residue_rmsf))

            # Save data and a default plot
            logger.debug("Saving RMSF data")
            self._save_data(self.data, "rmsf")
            
            logger.info("Generating RMSF plot")
            plot_path = self.plot()
            logger.info("RMSF plot saved to: %s", plot_path)

            return self.results
            
        except AnalysisError:
            logger.error("RMSF analysis failed with AnalysisError")
            raise
        except Exception as e:
            logger.error("RMSF analysis failed with unexpected error: %s", str(e))
            raise AnalysisError(f"RMSF analysis failed: {e}")

    def plot(
        self,
        data: Optional[Union[Sequence[float], np.ndarray]] = None,
        *,
        max_ticks: int = 8,        # Consistent with RMSD
        tick_step: Optional[int] = None,  # show every Nth tick (overrides max_ticks)
        rotate: int = 0,            # tick label rotation (default horizontal)
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
            Target maximum number of x tick labels (consistent with RMSD).
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
        logger.debug("Generating RMSF plot with parameters: max_ticks=%d, tick_step=%s, rotate=%d",
                    max_ticks, tick_step, rotate)
        
        if data is None:
            data = self.data
        if data is None:
            logger.error("No RMSF data available for plotting")
            raise AnalysisError("No RMSF data available to plot. Run the analysis first.")

        y = np.asarray(data, dtype=float).flatten()
        n = int(y.size)
        x = np.arange(n)
        
        logger.debug("Plotting %d RMSF values", n)

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        bar_kwargs = {"width": 0.9}
        if color is not None:
            bar_kwargs["color"] = color
            logger.debug("Using custom color: %s", color)
            
        ax.bar(x, y, **bar_kwargs)

        if self.compute_stat:
            mean_val = float(np.nanmean(y))
            std_val = float(np.nanstd(y))
            ax.axhline(mean_val, color="black", linestyle="--", linewidth=1.2, label="mean")
            ax.fill_between(
                [x.min() if x.size else 0, x.max() if x.size else 1],
                mean_val - std_val,
                mean_val + std_val,
                color="gray",
                alpha=0.2,
                label="±1 std",
            )

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if self.compute_stat:
            ax.legend(loc="best")

        ax.grid(axis="y", alpha=0.3)
        ax.grid(axis="x", alpha=0.12)

        # Use tick_step if provided, otherwise let auto_ticks handle it
        if tick_step is not None:
            step = max(1, int(tick_step))
            tick_positions = np.arange(0, n, step, dtype=float)
            if tick_positions.size == 0 or tick_positions[-1] != n - 1:
                tick_positions = np.append(tick_positions, float(n - 1))
            logger.debug("Using manual tick step: %d, resulting in %d tick positions", 
                        step, len(tick_positions))
        else:
            tick_positions = None
            logger.debug("Using automatic tick positioning with max_ticks=%d", max_ticks)

        # Let apply_slide_style handle tick generation cleanly
        applied = apply_slide_style(
            ax,
            x_values=x,
            y_values=y,
            x_ticks=tick_positions,
            x_max_ticks=max_ticks,
            zero_x=True,
            zero_y=True,
            x_tick_rotation=rotate,
        )
        logger.debug("Applied slide style with %d x-tick labels", len(ax.get_xticklabels()))

        # Apply rotation to x-tick labels if needed
        if rotate != 0:
            for label in ax.get_xticklabels():
                label.set_rotation(rotate)
                label.set_horizontalalignment("right")
            logger.debug("Applied %d degree rotation to x-tick labels", rotate)

        fig.tight_layout()
        outpath = self._save_plot(fig, "rmsf")
        plt.close(fig)
        
        logger.debug("RMSF plot saved to: %s", outpath)
        return Path(outpath)

