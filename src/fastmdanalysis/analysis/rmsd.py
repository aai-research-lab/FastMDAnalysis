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
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style

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
    _ALIASES = {
        "ref": "reference_frame",
        "reference": "reference_frame",
        "atom_indices": "atoms",
        "selection": "atoms",
    }
    
    def __init__(
        self,
        trajectory,
        reference_frame: int = 0,
        atoms: Optional[str] = None,
        align: bool = True,
        compute_stat: bool = False,
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        reference_frame : int
            Reference frame index (default 0). Negative indices allowed.
            Aliases: ref, reference
        atoms : str or None
            MDTraj atom selection string (e.g., "protein and name CA"). If None, all atoms are used.
            Aliases: atom_indices, selection
        align : bool
            If True, compute classical RMSD with optimal superposition (mdtraj.rmsd).
            If False, compute no-fit RMSD (raw differences).
        strict : bool
            If True, raise errors for unknown options. If False, log warnings.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "reference_frame": reference_frame,
            "atoms": atoms,
            "align": align,
            "compute_stat": compute_stat,
            "strict": strict,
        }
        analysis_opts.update(kwargs)
        
        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"reference_frame", "atoms", "align", "compute_stat", "strict", "output"},
            context="rmsd",
            warn=warn_unknown,
        )
        
        reference_frame = resolved.get("reference_frame", 0)
        atoms = resolved.get("atoms", None)
        align = resolved.get("align", True)
        compute_stat = resolved.get("compute_stat", False)
        base_kwargs = {
            k: v
            for k, v in resolved.items()
            if k not in ("reference_frame", "atoms", "align", "compute_stat", "strict")
        }
        
        super().__init__(trajectory, **base_kwargs)
        self.reference_frame = 0 if reference_frame is None else int(reference_frame)
        self.atoms = atoms
        self.align = bool(align)
        self.compute_stat = bool(compute_stat)
        self.strict = strict
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, np.ndarray] = {}

        logger.info("Initialized RMSD analysis: reference_frame=%d, align=%s, atoms=%s",
                   self.reference_frame, self.align, self.atoms if self.atoms else "ALL")

    def _select_atoms(self) -> Optional[np.ndarray]:
        """Return atom indices for selection, or None for all atoms."""
        if self.atoms:
            logger.debug("Selecting atoms: %s", self.atoms)
            sel = self.traj.topology.select(self.atoms)
            if sel is None or len(sel) == 0:
                raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
            logger.debug("Atom selection yielded %d atoms", len(sel))
            return sel
        logger.debug("Using all %d atoms", self.traj.n_atoms)
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
                logger.debug("Reference frame %d loaded successfully", self.reference_frame)
            except Exception as e:
                raise AnalysisError(f"Invalid reference frame index: {self.reference_frame}") from e

            atom_indices = self._select_atoms()

            logger.info(
                "Starting RMSD calculation: ref=%d, atoms=%s, align=%s, n_frames=%d, n_atoms=%d",
                self.reference_frame,
                self.atoms if self.atoms else "ALL",
                self.align,
                self.traj.n_frames,
                self.traj.n_atoms if atom_indices is None else len(atom_indices),
            )

            if self.align:
                logger.debug("Computing aligned RMSD")
                # md.rmsd performs optimal superposition internally
                rmsd_values = md.rmsd(self.traj, ref, atom_indices=atom_indices)
            else:
                logger.debug("Computing no-fit RMSD")
                # No-fit RMSD
                rmsd_values = _rmsd_no_fit(self.traj, ref, atom_indices=atom_indices)

            self.data = np.asarray(rmsd_values, dtype=float).reshape(-1, 1)
            self.results = {"rmsd": self.data}

            if self.compute_stat:
                mean_val = float(np.nanmean(rmsd_values))
                std_val = float(np.nanstd(rmsd_values))
                self.results["rmsd_stats"] = {"mean": mean_val, "std": std_val}
                self._save_data(
                    np.array([[mean_val, std_val]], dtype=float),
                    "rmsd_stats",
                    header="mean_nm std_nm",
                    fmt="%.6f",
                )

            # Save data and default plot
            logger.info("Saving RMSD data...")
            self._save_data(self.data, "rmsd", header="rmsd_nm", fmt="%.6f")
            
            logger.info("Generating RMSD plot...")
            self.plot()

            rmsd_range = (self.data.min(), self.data.max())
            logger.info("RMSD analysis complete: range [%.3f, %.3f] nm", rmsd_range[0], rmsd_range[1])
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            logger.exception("RMSD analysis failed")
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

        logger.debug("Generating RMSD plot")
        y = np.asarray(data, dtype=float).reshape(-1)
        x = np.arange(1, y.size + 1, dtype=int)

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

        out = self._save_plot(fig, "rmsd")
        plt.close(fig)
        logger.debug("RMSD plot saved to %s", out)
        return out

    def _save_plot(self, fig, name: str):
        """Save the figure as a PNG file in the output directory and log its path."""
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        logger.info("Plot saved to %s", plot_path)
        return plot_path