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
from ..utils.plotting import apply_slide_style, match_colorbar_font

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
        compute_stat: bool = False,
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
            "compute_stat": compute_stat,
            "strict": strict,
        }
        analysis_opts.update(kwargs)
        
        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"probe_radius", "atoms", "n_sphere_points", "compute_stat", "strict", "output"},
            context="sasa",
            warn=warn_unknown,
        )

        probe_radius = resolved.get("probe_radius", 0.14)
        atoms = resolved.get("atoms", None)
        n_sphere_points = resolved.get("n_sphere_points", None)
        compute_stat = resolved.get("compute_stat", False)
        base_kwargs = {
            k: v
            for k, v in resolved.items()
            if k not in ("probe_radius", "atoms", "n_sphere_points", "compute_stat", "strict")
        }

        super().__init__(trajectory, **base_kwargs)
        self.probe_radius = float(probe_radius)
        self.atoms = atoms
        self.n_sphere_points = int(n_sphere_points) if n_sphere_points is not None else None
        self.compute_stat = bool(compute_stat)
        self.strict = strict
        self.data: Optional[Dict[str, np.ndarray]] = None
        self.results: Dict[str, np.ndarray] = {}

        logger.info("Initialized SASA analysis: probe_radius=%.3f nm, atoms=%s", 
                   self.probe_radius, self.atoms if self.atoms else "ALL")

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
                logger.info("Selecting atoms: %s", self.atoms)
                sel = self.traj.topology.select(self.atoms)
                if sel is None or len(sel) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(sel)
                logger.info("Atom selection yielded %d atoms", len(sel))
            else:
                subtraj = self.traj
                logger.info("Using all %d atoms", self.traj.n_atoms)

            T = subtraj.n_frames
            logger.info(
                "Starting SASA calculation: %d frames, %d atoms, probe=%.3f nm",
                T, subtraj.n_atoms, self.probe_radius
            )

            residue_sasa = None
            sasa_kwargs = {"probe_radius": self.probe_radius}
            if self.n_sphere_points is not None:
                sasa_kwargs["n_sphere_points"] = self.n_sphere_points
                logger.debug("Using %d sphere points", self.n_sphere_points)
            
            try:
                # Newer MDTraj versions support mode="residue"
                logger.debug("Attempting residue-level SASA calculation...")
                residue_sasa = md.shrake_rupley(subtraj, mode="residue", **sasa_kwargs)
                logger.info("Residue-level SASA calculation successful")
                # shape (T, R)
            except TypeError:
                # Fallback: compute per-atom SASA then sum by residue
                logger.info("Falling back to atom-level SASA calculation with residue aggregation")
                atom_sasa = md.shrake_rupley(subtraj, **sasa_kwargs)  # (T, A)
                # Map atoms -> residue index (0..R-1 within subtraj topology)
                atom_res = np.array([a.residue.index for a in subtraj.topology.atoms], dtype=int)
                R = int(max(atom_res) + 1) if atom_res.size else 0
                residue_sasa = np.zeros((T, R), dtype=np.float32)
                for r in range(R):
                    residue_sasa[:, r] = atom_sasa[:, atom_res == r].sum(axis=1)
                logger.info("Aggregated SASA to %d residues", R)

            if residue_sasa.ndim != 2:
                raise AnalysisError("Unexpected residue_sasa shape; expected 2D (T, R).")
            T2, R = residue_sasa.shape
            if T2 != T:
                raise AnalysisError("residue_sasa frame dimension mismatch.")

            logger.info("Computing total and average SASA...")
            total_sasa = residue_sasa.sum(axis=1)                  # (T,)
            average_residue_sasa = residue_sasa.mean(axis=0)       # (R,)

            self.data = {
                "total_sasa": total_sasa,
                "residue_sasa": residue_sasa,
                "average_residue_sasa": average_residue_sasa,
            }
            self.results = self.data

            if self.compute_stat:
                mean_val = float(np.nanmean(total_sasa))
                std_val = float(np.nanstd(total_sasa))
                self.results["total_sasa_stats"] = {"mean": mean_val, "std": std_val}
                self._save_data(
                    np.array([[mean_val, std_val]], dtype=float),
                    "total_sasa_stats",
                    header="mean_nm2 std_nm2",
                    fmt="%.6f",
                )

            # Save data
            logger.info("Saving SASA data files...")
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
            logger.info("Generating SASA plots...")
            self.plot()

            logger.info("SASA analysis complete: %d frames, %d residues, total SASA range [%.3f, %.3f] nm²", 
                       T, R, total_sasa.min(), total_sasa.max())
            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            logger.exception("SASA analysis failed")
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

        logger.info("Generating SASA plots for option: %s", option)
        plots = {}
        if option in ("all", "total"):
            logger.debug("Generating total SASA plot")
            plots["total"] = self._plot_total_sasa(data["total_sasa"], **kwargs)
        if option in ("all", "residue"):
            logger.debug("Generating residue SASA heatmap")
            plots["residue"] = self._plot_residue_sasa(data["residue_sasa"], **kwargs)
        if option in ("all", "average"):
            logger.debug("Generating average residue SASA plot")
            plots["average"] = self._plot_average_residue_sasa(data["average_residue_sasa"], **kwargs)

        logger.info("SASA plots generated: %s", list(plots.keys()))
        return plots if option == "all" else plots[option]

    def _plot_total_sasa(self, total_sasa: np.ndarray, **kwargs):
        """Total SASA vs frame."""
        logger.debug("Plotting total SASA for %d frames", len(total_sasa))
        x = np.arange(1, total_sasa.shape[0] + 1, dtype=int)
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
        if self.compute_stat:
            mean_val = float(np.nanmean(total_sasa))
            std_val = float(np.nanstd(total_sasa))
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
            y_values=total_sasa,
            x_max_ticks=10,
            y_max_ticks=6,
            zero_x=True,
            zero_y=True,
        )
        fig.tight_layout()
        out = self._save_plot(fig, "total_sasa")
        plt.close(fig)
        logger.debug("Total SASA plot saved to %s", out)
        return out

    def _plot_residue_sasa(self, residue_sasa: np.ndarray, **kwargs):
        """Per-residue SASA heatmap (rows=residues, cols=frames)."""
        R = residue_sasa.shape[1]
        T = residue_sasa.shape[0]
        logger.debug("Plotting residue SASA heatmap: %d residues × %d frames", R, T)
        title = kwargs.get("title_residue", "Per-Residue SASA vs Frame")
        xlabel = kwargs.get("xlabel_residue", "Frame")
        ylabel = kwargs.get("ylabel_residue", "Residue Index")
        cmap = kwargs.get("cmap", "viridis")
        max_y_ticks = kwargs.get("max_y_ticks", 40)
        tick_step = kwargs.get("tick_step", None)  # overrides max_y_ticks if provided

        # Prepare data: (R, T), origin at lower so residue 1 is at bottom
        data = residue_sasa.T  # (R, T)

        fig, ax = plt.subplots(figsize=(12, 8))
        im = ax.imshow(data, aspect="auto", interpolation="none", cmap=cmap, origin="lower")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("SASA (nm²)")

        # Apply consistent colorbar font sizes
        match_colorbar_font(cbar, ax)

        # Improved residue tick generation to prevent overlapping
        # Use adaptive tick spacing based on number of residues
        max_ticks = 20  # Maximum number of y-axis ticks to prevent overlap
        
        if tick_step is not None:
            # User provided tick step
            step = max(1, int(tick_step))
            res_ticks = np.arange(0, R, step, dtype=int)
            if res_ticks.size == 0 or res_ticks[-1] != R - 1:
                res_ticks = np.append(res_ticks, R - 1)
            logger.debug("Using user-provided tick step %d for residue axis", step)
        elif R <= max_ticks:
            # For small proteins, show all residues
            res_ticks = np.arange(0, R, dtype=int)
            logger.debug("Using all %d residue ticks (<=%d residues)", R, max_ticks)
        else:
            # For larger proteins, calculate step size to get approximately max_ticks
            step = max(1, int(np.ceil(R / max_ticks)))
            res_ticks = np.arange(0, R, step, dtype=int)
            
            # Ensure we include the last residue
            if res_ticks.size == 0 or res_ticks[-1] != R - 1:
                res_ticks = np.append(res_ticks, R - 1)
            
            # If we still have too many ticks, increase step size
            while len(res_ticks) > max_ticks + 5:  # Allow some flexibility
                step += 1
                res_ticks = np.arange(0, R, step, dtype=int)
                if res_ticks.size == 0 or res_ticks[-1] != R - 1:
                    res_ticks = np.append(res_ticks, R - 1)
            
            logger.debug("Using %d residue ticks with adaptive step %d", len(res_ticks), step)

        frame_values = np.arange(1, residue_sasa.shape[0] + 1, dtype=int)
        apply_slide_style(
            ax,
            x_values=frame_values,
            y_ticks=res_ticks,
            zero_x=True,
            zero_y=True,
        )

        # Set y-ticks with 1-based residue numbering
        ax.set_yticks(res_ticks)
        ax.set_yticklabels(
            [str(int(v) + 1) for v in res_ticks],
            rotation=0,
            rotation_mode="anchor",
        )

        fig.tight_layout()
        out = self._save_plot(fig, "residue_sasa")
        plt.close(fig)
        logger.debug("Residue SASA heatmap saved to %s", out)
        return out

    def _plot_average_residue_sasa(self, average_sasa: np.ndarray, **kwargs):
        """Average per-residue SASA bar plot."""
        R = int(average_sasa.shape[0])
        logger.debug("Plotting average residue SASA for %d residues", R)
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

        residue_positions = np.arange(1, R + 1, dtype=int)
        y_values = average_sasa.flatten()

        # Improved residue tick generation for bar plot x-axis
        max_ticks = 20  # Maximum number of x-axis ticks to prevent overlapping
        
        if tick_step is not None:
            # User provided tick step
            step = max(1, int(tick_step))
            tick_positions = np.arange(1, R + 1, step, dtype=int)
            if tick_positions.size == 0 or tick_positions[-1] != R:
                tick_positions = np.append(tick_positions, R)
            logger.debug("Using user-provided tick step %d for residue axis", step)
        elif R <= max_ticks:
            # For small proteins, show all residues
            tick_positions = np.arange(1, R + 1, dtype=int)
            logger.debug("Using all %d residue ticks (<=%d residues)", R, max_ticks)
        else:
            # For larger proteins, calculate step size to get approximately max_ticks
            step = max(1, int(np.ceil(R / max_ticks)))
            tick_positions = np.arange(1, R + 1, step, dtype=int)
            
            # Ensure we include the last residue
            if tick_positions.size == 0 or tick_positions[-1] != R:
                tick_positions = np.append(tick_positions, R)
            
            # If we still have too many ticks, increase step size
            while len(tick_positions) > max_ticks + 5:  # Allow some flexibility
                step += 1
                tick_positions = np.arange(1, R + 1, step, dtype=int)
                if tick_positions.size == 0 or tick_positions[-1] != R:
                    tick_positions = np.append(tick_positions, R)
            
            logger.debug("Using %d residue ticks with adaptive step %d", len(tick_positions), step)

        apply_slide_style(
            ax,
            x_values=residue_positions,
            y_values=y_values,
            x_ticks=tick_positions,
            zero_y=True,
        )

        # Set x-ticks with residue numbers
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(
            [str(int(v)) for v in tick_positions],
            rotation=0,
            ha="center",
            rotation_mode="anchor",
        )

        fig.tight_layout()
        out = self._save_plot(fig, "average_residue_sasa")
        plt.close(fig)
        logger.debug("Average residue SASA plot saved to %s", out)
        return out

    def _save_plot(self, fig, name: str):
        """Save the figure as a PNG file in the output directory and log its path."""
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        logger.info("Plot saved to %s", plot_path)
        return plot_path
