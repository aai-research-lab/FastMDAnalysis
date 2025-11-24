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
from ..utils.plotting import apply_slide_style, auto_ticks

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
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "atoms": atoms,
            "per_residue": per_residue,
            "strict": strict,
        }
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"atoms", "per_residue", "strict", "output", "reference_frame"},
            context="rmsf",
            warn=warn_unknown,
        )

        atoms = resolved.get("atoms", None)
        per_residue = resolved.get("per_residue", False)
        base_kwargs = {k: v for k, v in resolved.items() 
                      if k not in ("atoms", "per_residue", "strict", "reference_frame")}

        super().__init__(trajectory, **base_kwargs)
        self.atoms: Optional[str] = atoms
        self.per_residue: bool = bool(per_residue)
        self.strict = strict

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
            If per_residue=True, also includes {"rmsf_per_residue": (R, 1) array}
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

            if self.per_residue:
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
                logger.info("RMSF: per-residue aggregation computed (%d residues)", n_residues)

            # Save data and a default plot
            self._save_data(self.data, "rmsf")
            self.plot()

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

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        bar_kwargs = {"width": 0.9}
        if color is not None:
            bar_kwargs["color"] = color
        ax.bar(x, y, **bar_kwargs)

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.grid(axis="y", alpha=0.3)
        ax.grid(axis="x", alpha=0.12)

        if tick_step is not None:
            step = max(1, int(tick_step))
            tick_positions = x[::step]
        else:
            tick_positions = auto_ticks(x, max_ticks=max_ticks, integer=True)
            if tick_positions is None or tick_positions.size == 0:
                tick_positions = x
        tick_positions = np.asarray(tick_positions, dtype=float)
        tick_positions = np.clip(tick_positions, 0, n - 1)
        ticklabels = [labels_all[int(pos)] for pos in tick_positions.astype(int)]

        ax.set_xticks(tick_positions)
        apply_slide_style(
            ax,
            x_ticks=tick_positions,
            y_values=y,
            integer_x=True,
            zero_x=True,
            zero_y=True,
            x_tick_rotation=rotate,
        )
        tick_font = ax.get_xticklabels()[0].get_fontsize() if ax.get_xticklabels() else None
        ax.set_xticklabels(
            ticklabels,
            rotation=rotate,
            ha="right",
            rotation_mode="anchor",
            fontsize=tick_font,
        )

        fig.tight_layout()
        outpath = self._save_plot(fig, "rmsf")
        plt.close(fig)
        return Path(outpath)

