# FastMDAnalysis/src/fastmdanalysis/analysis/qvalue.py

"""
Fraction of Native Contacts (Q-Value) Analysis Module

Calculates the fraction of native contacts (Q-value) according to the
Best-Hummer-Eaton definition:

    Q(X) = 1/|S| * sum_{(i,j) in S} 1 / (1 + exp[beta * (rij(X) - lambda * r0ij)])

where:
    - X is a conformation
    - r_ij(X) is the distance between atoms i and j in conformation X
    - r0_ij is the distance in the reference (native) state
    - S is the set of contact pairs (heavy atoms >3 residues apart, distance < cutoff)
    - beta and lambda are tuning constants

Reference: Best, Hummer, and Eaton, "Native contacts determine protein folding 
mechanisms in atomistic simulations" PNAS (2013) 10.1073/pnas.1311599110

Outputs
-------
- qvalue.dat           : (N,) array of Q values per frame
- qvalue_metadata.json : JSON file with native contacts count and parameters
- qvalue.png           : line plot of Q vs frame with metadata annotation

Notes
-----
- Default parameters: beta=50 nm^-1, lambda=1.8, native_cutoff=0.45 nm
- Native contacts are identified between heavy atoms >3 residues apart
- Only atoms within native_cutoff distance in the reference frame are considered
"""

from __future__ import annotations

from typing import Optional, Dict, Any
import logging
import json
from itertools import combinations
import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style

logger = logging.getLogger(__name__)


class QAnalysis(BaseAnalysis):
    """Fraction of Native Contacts (Q-Value) Analysis."""
    
    _ALIASES = {
        "ref": "reference_frame",
        "reference": "reference_frame",
        "atom_indices": "atom_selection",
        "selection": "atom_selection",
    }
    
    def __init__(
        self,
        trajectory,
        reference_frame: int = 0,
        beta_const: float = 50.0,
        lambda_const: float = 1.8,
        native_cutoff: float = 0.45,
        atom_selection: Optional[str] = None,
        compute_stat: bool = False,
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        reference_frame : int, optional
            Frame index to use as the native/reference state (default: 0).
            Aliases: ref, reference
        beta_const : float, optional
            Beta constant in nanometers^-1 (default: 50.0 nm^-1).
            Controls steepness of the sigmoid function.
        lambda_const : float, optional
            Lambda constant (default: 1.8).
            Scaling factor for native distances.
        native_cutoff : float, optional
            Cutoff distance in nanometers for identifying native contacts (default: 0.45 nm).
            Heavy atom pairs within this distance in the reference frame are considered.
        atom_selection : str, optional
            MDTraj atom selection string for analysis (default: None, uses heavy atoms).
            Aliases: atom_indices, selection
        strict : bool, optional
            If True, raise errors for unknown options. If False, log warnings.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "reference_frame": reference_frame,
            "beta_const": beta_const,
            "lambda_const": lambda_const,
            "native_cutoff": native_cutoff,
            "atom_selection": atom_selection,
            "compute_stat": compute_stat,
            "strict": strict,
        }
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {
                "reference_frame",
                "beta_const",
                "lambda_const",
                "native_cutoff",
                "atom_selection",
                "compute_stat",
                "strict",
                "output",
            },
            context="qvalue",
            warn=warn_unknown,
        )

        reference_frame = resolved.get("reference_frame", 0)
        beta_const = resolved.get("beta_const", 50.0)
        lambda_const = resolved.get("lambda_const", 1.8)
        native_cutoff = resolved.get("native_cutoff", 0.45)
        atom_selection = resolved.get("atom_selection", None)
        compute_stat = resolved.get("compute_stat", False)
        base_kwargs = {
            k: v
            for k, v in resolved.items()
            if k
            not in (
                "reference_frame",
                "beta_const",
                "lambda_const",
                "native_cutoff",
                "atom_selection",
                "compute_stat",
                "strict",
            )
        }

        super().__init__(trajectory, **base_kwargs)
        self.reference_frame = reference_frame
        self.beta_const = beta_const
        self.lambda_const = lambda_const
        self.native_cutoff = native_cutoff
        self.atom_selection = atom_selection
        self.compute_stat = bool(compute_stat)
        self.strict = strict
        
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, np.ndarray] = {}
        self.native_contacts_count: int = 0
        self.metadata: Dict[str, Any] = {}

        logger.info(
            "Initialized Q-value analysis: reference_frame=%d, beta=%.1f, lambda=%.1f, cutoff=%.2f nm",
            self.reference_frame, self.beta_const, self.lambda_const, self.native_cutoff
        )

    def _get_heavy_pairs(self) -> np.ndarray:
        """
        Get pairs of heavy atoms that are >3 residues apart.

        Returns
        -------
        np.ndarray
            Array of shape (N_pairs, 2) with atom indices.
        """
        heavy = self.traj.topology.select_atom_indices("heavy")
        heavy_pairs = np.array(
            [
                (i, j) for (i, j) in combinations(heavy, 2)
                if abs(
                    self.traj.topology.atom(i).residue.index
                    - self.traj.topology.atom(j).residue.index
                ) > 3
            ]
        )
        logger.debug("Found %d heavy atom pairs >3 residues apart", len(heavy_pairs))
        return heavy_pairs

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute the fraction of native contacts (Q-value) for each frame.

        Returns
        -------
        dict
            {"qvalue": (N,) array of Q values per frame}
        """
        try:
            if self.reference_frame < 0 or self.reference_frame >= self.traj.n_frames:
                raise AnalysisError(
                    f"Reference frame {self.reference_frame} out of range [0, {self.traj.n_frames})"
                )

            logger.info(
                "Starting Q-value calculation: reference_frame=%d, n_frames=%d",
                self.reference_frame, self.traj.n_frames
            )

            # Get heavy atom pairs >3 residues apart
            heavy_pairs = self._get_heavy_pairs()

            # Compute distances in the reference frame
            logger.debug("Computing distances in reference frame %d...", self.reference_frame)
            heavy_pairs_distances = md.compute_distances(
                self.traj[self.reference_frame], heavy_pairs
            )[0]

            # Identify native contacts (distance < cutoff)
            native_contacts = heavy_pairs[heavy_pairs_distances < self.native_cutoff]
            self.native_contacts_count = len(native_contacts)
            logger.info("Identified %d native contacts", self.native_contacts_count)

            if self.native_contacts_count == 0:
                raise AnalysisError(
                    f"No native contacts found with cutoff {self.native_cutoff} nm. "
                    "Try increasing the native_cutoff parameter."
                )

            # Compute distances for all frames
            logger.debug("Computing distances for all %d frames...", self.traj.n_frames)
            r = md.compute_distances(self.traj, native_contacts)  # shape (N_frames, N_contacts)
            r0 = md.compute_distances(self.traj[self.reference_frame], native_contacts)[0]  # shape (N_contacts,)

            # Compute Q-values using Best-Hummer-Eaton formula
            logger.debug("Computing Q-values...")
            q_values = np.mean(
                1.0 / (1.0 + np.exp(self.beta_const * (r - self.lambda_const * r0))),
                axis=1
            )

            self.data = np.asarray(q_values, dtype=float).reshape(-1, 1)
            self.results = {"qvalue": self.data}

            if self.compute_stat:
                mean_val = float(np.nanmean(q_values))
                std_val = float(np.nanstd(q_values))
                self.results["qvalue_stats"] = {"mean": mean_val, "std": std_val}
                self._save_data(
                    np.array([[mean_val, std_val]], dtype=float),
                    "qvalue_stats",
                    header="mean std",
                    fmt="%.6f",
                )

            # Save data
            logger.info("Saving Q-value data...")
            self._save_data(self.data, "qvalue", header="qvalue", fmt="%.6f")

            # Save metadata
            self._save_metadata()

            # Generate plot
            logger.info("Generating Q-value plot...")
            self.plot()

            q_range = (self.data.min(), self.data.max())
            logger.info(
                "Q-value analysis complete: range [%.4f, %.4f], native_contacts=%d",
                q_range[0], q_range[1], self.native_contacts_count
            )
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            logger.exception("Q-value analysis failed")
            raise AnalysisError(f"Q-value analysis failed: {e}")

    def _save_metadata(self):
        """Save analysis metadata to JSON file."""
        self.metadata = {
            "native_contacts_count": int(self.native_contacts_count),
            "reference_frame": self.reference_frame,
            "beta_const_nm_inv": self.beta_const,
            "lambda_const": self.lambda_const,
            "native_cutoff_nm": self.native_cutoff,
            "atom_selection": self.atom_selection if self.atom_selection else "heavy",
            "n_frames": int(self.traj.n_frames),
            "n_atoms": int(self.traj.n_atoms),
        }

        metadata_path = self.outdir / "qvalue_metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(self.metadata, f, indent=2)
        logger.info("Metadata saved to %s", metadata_path)

    def plot(self, data: Optional[np.ndarray] = None, **kwargs):
        """
        Generate a plot of Q-value versus frame number.

        Parameters
        ----------
        data : array-like, optional
            Q-value data to plot; if None, uses data from `run()`.
        kwargs : dict
            Matplotlib options:
              - title (str): default "Fraction of Native Contacts (Q-Value) vs Frame"
              - xlabel (str): default "Frame"
              - ylabel (str): default "Q-Value"
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
            raise AnalysisError("No Q-value data available to plot. Run the analysis first.")

        logger.debug("Generating Q-value plot")
        y = np.asarray(data, dtype=float).reshape(-1)
        x = np.arange(1, y.size + 1, dtype=int)

        title = kwargs.get("title", "Fraction of Native Contacts (Q-Value) vs Frame")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "Q-Value")
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

        # Add metadata annotation
        if self.metadata:
            info_text = f"Native Contacts: {self.native_contacts_count}\n"
            info_text += f"Reference: Frame {self.reference_frame}"
            ax.text(
                0.98, 0.02,
                info_text,
                transform=ax.transAxes,
                fontsize=9,
                verticalalignment="bottom",
                horizontalalignment="right",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
            )

        apply_slide_style(
            ax,
            x_values=x,
            y_values=y,
            x_max_ticks=8,
            y_max_ticks=6,
            zero_x=True,
            zero_y=False,
        )
        if self.compute_stat:
            ax.legend(loc="best")
        fig.tight_layout()

        out = self._save_plot(fig, "qvalue")
        plt.close(fig)
        logger.debug("Q-value plot saved to %s", out)
        return out
