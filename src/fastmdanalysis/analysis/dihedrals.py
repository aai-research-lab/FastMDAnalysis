# FastMDAnalysis/src/fastmdanalysis/analysis/dihedrals.py

"""
Dihedral Angles Analysis Module

Computes backbone dihedral angles (phi, psi, omega) from MD trajectories,
averages them circularly per residue, and generates plots including
Ramachandran plots for phi/psi analysis.

Uses MDTraj's compute_phi, compute_psi, compute_omega functions.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Union, List
from pathlib import Path
import logging

import numpy as np
import mdtraj as md
from scipy.stats import circmean, circstd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .base import BaseAnalysis, AnalysisError
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style

logger = logging.getLogger(__name__)


class PhiAnalysis(BaseAnalysis):
    """
    Phi backbone dihedral angle analysis with per-residue averaging and plotting.
    """

    _ALIASES = {
        "residue": "residues",
        "residue_selection": "residues",
    }

    def __init__(
        self,
        trajectory: md.Trajectory,
        residues: Optional[Union[int, Sequence[int]]] = None,
        units: str = "degrees",
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            MDTraj trajectory to analyze.
        residues : int or sequence of int, optional
            Residue indices to analyze (0-based). If None, analyze all residues.
        units : str
            Units for output angles: 'degrees' or 'radians'.
        strict : bool
            If True, raise errors for unknown options.
        kwargs : dict
            Passed through to BaseAnalysis (e.g., output).
        """
        logger.info("Initializing Phi analysis")
        logger.debug("Input parameters: residues=%s, units=%s, strict=%s",
                    residues, units, strict)

        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "residues": residues,
            "units": units,
            "strict": strict,
        }
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"residues", "units", "strict", "output"},
            context="phi",
            warn=warn_unknown,
        )

        residues = resolved.get("residues", None)
        units = resolved.get("units", "degrees")
        base_kwargs = {k: v for k, v in resolved.items()
                      if k not in ("residues", "units", "strict")}

        super().__init__(trajectory, **base_kwargs)
        self.residues: Optional[Union[int, Sequence[int]]] = residues
        self.units: str = units
        self.strict = strict

        # Residue indices corresponding to the rows in self.data
        self.residue_indices: Optional[np.ndarray] = None

        # Populated during run()
        self.data: Optional[np.ndarray] = None
        self.results: Dict[str, np.ndarray] = {}

        logger.info("Phi analysis initialized")

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute phi angles and average per residue.

        Returns
        -------
        dict
            {"phi_avg": (N, 1) array of average phi angles per residue}
        """
        logger.info("Starting Phi analysis")
        try:
            indices, angles = md.compute_phi(self.traj)
            if angles.size == 0:
                raise AnalysisError("No phi angles found in trajectory (no protein?)")

            # Restrict computation to selected residues (0-based dihedral/residue index)
            full_res_indices = np.arange(angles.shape[1])
            if self.residues is not None:
                res_list = [self.residues] if isinstance(self.residues, int) else list(self.residues)
                angles = angles[:, res_list]
                self.residue_indices = np.asarray(res_list, dtype=int)
            else:
                self.residue_indices = full_res_indices

            # Circular mean and std per residue
            n_residues = angles.shape[1]
            avg_angles = np.zeros(n_residues)
            std_angles = np.zeros(n_residues)
            for i in range(n_residues):
                res_angles = angles[:, i]
                # Remove NaN values (missing angles)
                valid_angles = res_angles[~np.isnan(res_angles)]
                if len(valid_angles) > 0:
                    avg_angles[i] = circmean(valid_angles, high=np.pi, low=-np.pi)
                    std_angles[i] = circstd(valid_angles, high=np.pi, low=-np.pi)
                else:
                    avg_angles[i] = np.nan
                    std_angles[i] = np.nan

            # Convert units
            if self.units == "degrees":
                avg_angles = np.degrees(avg_angles)
                std_angles = np.degrees(std_angles)

            self.data = avg_angles.reshape(-1, 1)
            # If residues were provided, self.data is already filtered.
            self.results = {
                "phi_avg": self.data,
                "phi_avg_filtered": self.data,
                "phi_residues": self.residue_indices,
                "phi_std": std_angles.reshape(-1, 1),
            }

            # Save data
            self._save_data(self.data, "phi_avg", header=f"phi_avg_{self.units}")

            logger.info("Phi analysis completed - %d residues analyzed", n_residues)

            # Generate plot
            plot_path = self.plot()
            logger.info("Phi plot saved to: %s", plot_path)

            return self.results

        except AnalysisError:
            raise
        except Exception as e:
            logger.error("Phi analysis failed: %s", str(e))
            raise AnalysisError(f"Phi analysis failed: {e}")

    def plot(
        self,
        data: Optional[Union[Sequence[float], np.ndarray]] = None,
        *,
        std_data: Optional[Union[Sequence[float], np.ndarray]] = None,
        residues: Optional[Union[int, Sequence[int]]] = None,
        highlight_residues: Optional[Union[int, Sequence[int]]] = None,
        figsize=(12, 6),
        title: str = "Average Phi Angles per Residue",
        xlabel: str = "Residue Index",
        ylabel: Optional[str] = None,
        color: str = "blue",
        **kwargs
    ) -> str:
        """
        Generate bar plot of average phi angles per residue.

        Parameters
        ----------
        data : array-like, optional
            Data to plot. If None, uses computed data.
        residues : int or sequence, optional
            Residue indices to plot. If None, plot all.
        highlight_residues : int or sequence, optional
            Residues to highlight in different color.
        figsize, title, xlabel, ylabel, color : matplotlib options

        Returns
        -------
        str
            Path to saved plot.
        """
        if data is None:
            data = self.data
        if std_data is None:
            std_data = self.results.get("phi_std")
        if data is None:
            raise AnalysisError("No phi data available to plot.")

        y = np.asarray(data, dtype=float).flatten()
        yerr = None
        if std_data is not None:
            yerr = np.asarray(std_data, dtype=float).flatten()

        # X-axis should reflect residue indices of the computed data (not 0..N-1)
        if self.residue_indices is not None and len(self.residue_indices) == len(y):
            x = self.residue_indices.astype(int)
        else:
            x = np.arange(len(y))

        # Filter residues
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(x, residues)
            x = x[mask]
            y = y[mask]
            if yerr is not None:
                yerr = yerr[mask]

        plot_matrix = np.column_stack([x, y])
        header = f"residue_index phi_mean_{self.units}"
        if yerr is not None:
            plot_matrix = np.column_stack([plot_matrix, yerr])
            header = f"{header} phi_std_{self.units}"
        self._save_data(plot_matrix, "phi_avg_plot", header=header)

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            color=color,
            ecolor=kwargs.pop("ecolor", "gray"),
            capsize=kwargs.pop("capsize", 3),
            **kwargs,
        )

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        if ylabel is None:
            ylabel = f"Phi Angle ({self.units})"
        ax.set_ylabel(ylabel)

        limit = 180.0 if self.units == "degrees" else np.pi
        pad = 20.0 if self.units == "degrees" else (np.pi / 9.0)
        ax.set_ylim(-(limit + pad), limit + pad)

        # Highlight specific residues
        if highlight_residues is not None:
            if isinstance(highlight_residues, int):
                highlight_residues = [highlight_residues]
            for res in highlight_residues:
                if res in x:
                    idx = np.where(x == res)[0][0]
                    ax.scatter([x[idx]], [y[idx]], color="red", alpha=0.7, zorder=3)

        apply_slide_style(ax, x_values=x, y_values=y)
        fig.tight_layout()
        outpath = self._save_plot(fig, "phi_avg")
        plt.close(fig)

        return outpath


class PsiAnalysis(BaseAnalysis):
    """
    Psi backbone dihedral angle analysis.
    """

    # Similar structure to PhiAnalysis, but for psi angles
    _ALIASES = {
        "residue": "residues",
        "residue_selection": "residues",
    }

    def __init__(
        self,
        trajectory: md.Trajectory,
        residues: Optional[Union[int, Sequence[int]]] = None,
        units: str = "degrees",
        strict: bool = False,
        **kwargs
    ):
        logger.info("Initializing Psi analysis")
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {"residues": residues, "units": units, "strict": strict}
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved, {"residues", "units", "strict", "output"}, context="psi", warn=warn_unknown
        )

        residues = resolved.get("residues", None)
        units = resolved.get("units", "degrees")
        base_kwargs = {k: v for k, v in resolved.items()
                      if k not in ("residues", "units", "strict")}

        super().__init__(trajectory, **base_kwargs)
        self.residues = residues
        self.units = units
        self.strict = strict
        self.data = None
        self.results = {}

        self.residue_indices: Optional[np.ndarray] = None

    def run(self) -> Dict[str, np.ndarray]:
        logger.info("Starting Psi analysis")
        try:
            indices, angles = md.compute_psi(self.traj)
            if angles.size == 0:
                raise AnalysisError("No psi angles found in trajectory")

            full_res_indices = np.arange(angles.shape[1])
            if self.residues is not None:
                res_list = [self.residues] if isinstance(self.residues, int) else list(self.residues)
                angles = angles[:, res_list]
                self.residue_indices = np.asarray(res_list, dtype=int)
            else:
                self.residue_indices = full_res_indices

            n_residues = angles.shape[1]
            avg_angles = np.zeros(n_residues)
            std_angles = np.zeros(n_residues)
            for i in range(n_residues):
                res_angles = angles[:, i]
                valid_angles = res_angles[~np.isnan(res_angles)]
                if len(valid_angles) > 0:
                    avg_angles[i] = circmean(valid_angles, high=np.pi, low=-np.pi)
                    std_angles[i] = circstd(valid_angles, high=np.pi, low=-np.pi)
                else:
                    avg_angles[i] = np.nan
                    std_angles[i] = np.nan

            if self.units == "degrees":
                avg_angles = np.degrees(avg_angles)
                std_angles = np.degrees(std_angles)

            self.data = avg_angles.reshape(-1, 1)
            self.results = {
                "psi_avg": self.data,
                "psi_avg_filtered": self.data,
                "psi_residues": self.residue_indices,
                "psi_std": std_angles.reshape(-1, 1),
            }

            self._save_data(self.data, "psi_avg", header=f"psi_avg_{self.units}")
            plot_path = self.plot()
            logger.info("Psi plot saved to: %s", plot_path)

            return self.results

        except Exception as e:
            raise AnalysisError(f"Psi analysis failed: {e}")

    def plot(self, **kwargs) -> str:
        # Similar to PhiAnalysis.plot but with psi-specific defaults
        kwargs.setdefault("title", "Average Psi Angles per Residue")
        kwargs.setdefault("ylabel", f"Psi Angle ({self.units})")
        kwargs.setdefault("color", "green")
        
        if "data" not in kwargs:
            kwargs["data"] = self.data
        if "std_data" not in kwargs:
            kwargs["std_data"] = self.results.get("psi_std")
        if kwargs["data"] is None:
            raise AnalysisError("No psi data available to plot.")

        y = np.asarray(kwargs["data"], dtype=float).flatten()
        yerr = None
        if kwargs.get("std_data") is not None:
            yerr = np.asarray(kwargs["std_data"], dtype=float).flatten()
        if self.residue_indices is not None and len(self.residue_indices) == len(y):
            x = self.residue_indices.astype(int)
        else:
            x = np.arange(len(y))

        # Filter residues
        residues = kwargs.get("residues")
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(x, residues)
            x = x[mask]
            y = y[mask]
            if yerr is not None:
                yerr = yerr[mask]

        plot_matrix = np.column_stack([x, y])
        header = f"residue_index psi_mean_{self.units}"
        if yerr is not None:
            plot_matrix = np.column_stack([plot_matrix, yerr])
            header = f"{header} psi_std_{self.units}"
        self._save_data(plot_matrix, "psi_avg_plot", header=header)

        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (12, 6)))
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            color=kwargs.get("color", "green"),
            ecolor=kwargs.get("ecolor", "gray"),
            capsize=kwargs.get("capsize", 3),
            **{
                k: v
                for k, v in kwargs.items()
                if k
                not in [
                    "data",
                    "std_data",
                    "residues",
                    "figsize",
                    "title",
                    "xlabel",
                    "ylabel",
                    "color",
                    "highlight_residues",
                    "ecolor",
                    "capsize",
                ]
            },
        )

        ax.set_title(kwargs.get("title", "Average Psi Angles per Residue"))
        ax.set_xlabel(kwargs.get("xlabel", "Residue Index"))
        ax.set_ylabel(kwargs.get("ylabel", f"Psi Angle ({self.units})"))

        limit = 180.0 if self.units == "degrees" else np.pi
        pad = 20.0 if self.units == "degrees" else (np.pi / 9.0)
        ax.set_ylim(-(limit + pad), limit + pad)

        # Highlight specific residues
        highlight_residues = kwargs.get("highlight_residues")
        if highlight_residues is not None:
            if isinstance(highlight_residues, int):
                highlight_residues = [highlight_residues]
            for res in highlight_residues:
                if res in x:
                    idx = np.where(x == res)[0][0]
                    ax.scatter([x[idx]], [y[idx]], color="red", alpha=0.7, zorder=3)

        apply_slide_style(ax, x_values=x, y_values=y)
        fig.tight_layout()
        outpath = self._save_plot(fig, "psi_avg")
        plt.close(fig)

        return outpath


class OmegaAnalysis(BaseAnalysis):
    """
    Omega backbone dihedral angle analysis.
    """

    _ALIASES = {
        "residue": "residues",
        "residue_selection": "residues",
    }

    def __init__(
        self,
        trajectory: md.Trajectory,
        residues: Optional[Union[int, Sequence[int]]] = None,
        units: str = "degrees",
        strict: bool = False,
        **kwargs
    ):
        logger.info("Initializing Omega analysis")
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {"residues": residues, "units": units, "strict": strict}
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved, {"residues", "units", "strict", "output"}, context="omega", warn=warn_unknown
        )

        residues = resolved.get("residues", None)
        units = resolved.get("units", "degrees")
        base_kwargs = {k: v for k, v in resolved.items()
                      if k not in ("residues", "units", "strict")}

        super().__init__(trajectory, **base_kwargs)
        self.residues = residues
        self.units = units
        self.strict = strict
        self.data = None
        self.results = {}

        self.residue_indices: Optional[np.ndarray] = None

    def run(self) -> Dict[str, np.ndarray]:
        logger.info("Starting Omega analysis")
        try:
            indices, angles = md.compute_omega(self.traj)
            if angles.size == 0:
                raise AnalysisError("No omega angles found in trajectory")

            full_res_indices = np.arange(angles.shape[1])
            if self.residues is not None:
                res_list = [self.residues] if isinstance(self.residues, int) else list(self.residues)
                angles = angles[:, res_list]
                self.residue_indices = np.asarray(res_list, dtype=int)
            else:
                self.residue_indices = full_res_indices

            n_residues = angles.shape[1]
            avg_angles = np.zeros(n_residues)
            std_angles = np.zeros(n_residues)
            for i in range(n_residues):
                res_angles = angles[:, i]
                valid_angles = res_angles[~np.isnan(res_angles)]
                if len(valid_angles) > 0:
                    avg_angles[i] = circmean(valid_angles, high=np.pi, low=-np.pi)
                    std_angles[i] = circstd(valid_angles, high=np.pi, low=-np.pi)
                else:
                    avg_angles[i] = np.nan
                    std_angles[i] = np.nan

            if self.units == "degrees":
                avg_angles = np.degrees(avg_angles)
                std_angles = np.degrees(std_angles)

            self.data = avg_angles.reshape(-1, 1)
            self.results = {
                "omega_avg": self.data,
                "omega_avg_filtered": self.data,
                "omega_residues": self.residue_indices,
                "omega_std": std_angles.reshape(-1, 1),
            }

            self._save_data(self.data, "omega_avg", header=f"omega_avg_{self.units}")
            plot_path = self.plot()
            logger.info("Omega plot saved to: %s", plot_path)

            return self.results

        except Exception as e:
            raise AnalysisError(f"Omega analysis failed: {e}")

    def plot(self, **kwargs) -> str:
        kwargs.setdefault("title", "Average Omega Angles per Residue")
        kwargs.setdefault("ylabel", f"Omega Angle ({self.units})")
        kwargs.setdefault("color", "orange")
        
        if "data" not in kwargs:
            kwargs["data"] = self.data
        if "std_data" not in kwargs:
            kwargs["std_data"] = self.results.get("omega_std")
        if kwargs["data"] is None:
            raise AnalysisError("No omega data available to plot.")

        y = np.asarray(kwargs["data"], dtype=float).flatten()
        yerr = None
        if kwargs.get("std_data") is not None:
            yerr = np.asarray(kwargs["std_data"], dtype=float).flatten()
        if self.residue_indices is not None and len(self.residue_indices) == len(y):
            x = self.residue_indices.astype(int)
        else:
            x = np.arange(len(y))

        # Filter residues
        residues = kwargs.get("residues")
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(x, residues)
            x = x[mask]
            y = y[mask]
            if yerr is not None:
                yerr = yerr[mask]

        plot_matrix = np.column_stack([x, y])
        header = f"residue_index omega_mean_{self.units}"
        if yerr is not None:
            plot_matrix = np.column_stack([plot_matrix, yerr])
            header = f"{header} omega_std_{self.units}"
        self._save_data(plot_matrix, "omega_avg_plot", header=header)

        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (12, 6)))
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            color=kwargs.get("color", "orange"),
            ecolor=kwargs.get("ecolor", "gray"),
            capsize=kwargs.get("capsize", 3),
            **{
                k: v
                for k, v in kwargs.items()
                if k
                not in [
                    "data",
                    "std_data",
                    "residues",
                    "figsize",
                    "title",
                    "xlabel",
                    "ylabel",
                    "color",
                    "highlight_residues",
                    "ecolor",
                    "capsize",
                ]
            },
        )

        ax.set_title(kwargs.get("title", "Average Omega Angles per Residue"))
        ax.set_xlabel(kwargs.get("xlabel", "Residue Index"))
        ax.set_ylabel(kwargs.get("ylabel", f"Omega Angle ({self.units})"))

        limit = 180.0 if self.units == "degrees" else np.pi
        pad = 20.0 if self.units == "degrees" else (np.pi / 9.0)
        ax.set_ylim(-(limit + pad), limit + pad)

        # Highlight specific residues
        highlight_residues = kwargs.get("highlight_residues")
        if highlight_residues is not None:
            if isinstance(highlight_residues, int):
                highlight_residues = [highlight_residues]
            for res in highlight_residues:
                if res in x:
                    idx = np.where(x == res)[0][0]
                    ax.scatter([x[idx]], [y[idx]], color="red", alpha=0.7, zorder=3)

        apply_slide_style(ax, x_values=x, y_values=y)
        fig.tight_layout()
        outpath = self._save_plot(fig, "omega_avg")
        plt.close(fig)

        return outpath


class DihedralsAnalysis(BaseAnalysis):
    """
    Combined dihedral analysis for phi, psi, omega with Ramachandran plotting.
    """

    _ALIASES = {
        "residue": "residues",
        "residue_selection": "residues",
    }

    def __init__(
        self,
        trajectory: md.Trajectory,
        types: Sequence[str] = ["phi", "psi", "omega"],
        residues: Optional[Union[int, Sequence[int]]] = None,
        units: str = "degrees",
        strict: bool = False,
        **kwargs
    ):
        logger.info("Initializing Dihedrals analysis")
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {"types": types, "residues": residues, "units": units, "strict": strict}
        analysis_opts.update(kwargs)

        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved, {"types", "residues", "units", "strict", "output"}, context="dihedrals", warn=warn_unknown
        )

        types = resolved.get("types", ["phi", "psi", "omega"])
        residues = resolved.get("residues", None)
        units = resolved.get("units", "degrees")
        base_kwargs = {k: v for k, v in resolved.items()
                      if k not in ("types", "residues", "units", "strict")}

        super().__init__(trajectory, **base_kwargs)
        self.types = types
        self.residues = residues
        self.units = units
        self.strict = strict
        self.data = None
        self.results = {}

    def run(self) -> Dict[str, np.ndarray]:
        logger.info("Starting Dihedrals analysis")
        results = {}

        base_outdir = self.outdir

        for angle_type in self.types:
            if angle_type == "phi":
                analysis = PhiAnalysis(
                    self.traj,
                    residues=self.residues,
                    units=self.units,
                    output=str(base_outdir / "phi"),
                )
            elif angle_type == "psi":
                analysis = PsiAnalysis(
                    self.traj,
                    residues=self.residues,
                    units=self.units,
                    output=str(base_outdir / "psi"),
                )
            elif angle_type == "omega":
                analysis = OmegaAnalysis(
                    self.traj,
                    residues=self.residues,
                    units=self.units,
                    output=str(base_outdir / "omega"),
                )
            else:
                continue

            analysis.run()
            results.update(analysis.results)

        self.results = results
        self.data = results  # Store combined results

        # Generate Ramachandran if phi and psi are computed
        if "phi_avg" in results and "psi_avg" in results:
            plot_path = self.plot_ramachandran()
            logger.info("Ramachandran plot saved to: %s", plot_path)

        return results

    def plot_ramachandran(
        self,
        phi_data: Optional[np.ndarray] = None,
        psi_data: Optional[np.ndarray] = None,
        phi_std: Optional[np.ndarray] = None,
        psi_std: Optional[np.ndarray] = None,
        residues: Optional[Union[int, Sequence[int]]] = None,
        figsize=(8, 8),
        title: str = "Ramachandran Plot (Average Angles)",
        **kwargs
    ) -> str:
        """
        Generate Ramachandran plot of average psi vs phi angles.
        """
        if phi_data is None and "phi_avg" in self.results:
            phi_data = self.results["phi_avg"].flatten()
        if psi_data is None and "psi_avg" in self.results:
            psi_data = self.results["psi_avg"].flatten()
        if phi_std is None and "phi_std" in self.results:
            phi_std = self.results["phi_std"].flatten()
        if psi_std is None and "psi_std" in self.results:
            psi_std = self.results["psi_std"].flatten()

        if phi_data is None or psi_data is None:
            raise AnalysisError("Phi and psi data required for Ramachandran plot")

        x = phi_data
        y = psi_data

        # Use residue indices when the per-angle analyses were residue-filtered
        res_indices = self.results.get("phi_residues")
        if res_indices is None:
            res_indices = np.arange(len(x))
        else:
            res_indices = np.asarray(res_indices, dtype=int)

        if residues is None:
            residues = self.residues
        res_list = None
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            res_list = list(residues)
            mask = np.isin(res_indices, res_list)
            x = x[mask]
            y = y[mask]
            res_indices = res_indices[mask]
            if phi_std is not None:
                phi_std = phi_std[mask]
            if psi_std is not None:
                psi_std = psi_std[mask]
        avg_matrix = np.column_stack([res_indices, x, y])
        header = f"residue_index phi_mean_{self.units} psi_mean_{self.units}"
        if phi_std is not None and psi_std is not None:
            avg_matrix = np.column_stack([avg_matrix, phi_std, psi_std])
            header = f"{header} phi_std_{self.units} psi_std_{self.units}"
        self._save_data(avg_matrix, "ramachandran_avg", header=header)

        avg_matrix = np.column_stack([res_indices, x, y])
        header = f"residue_index phi_mean_{self.units} psi_mean_{self.units}"
        if phi_std is not None and psi_std is not None:
            avg_matrix = np.column_stack([avg_matrix, phi_std, psi_std])
            header = f"{header} phi_std_{self.units} psi_std_{self.units}"
        self._save_data(avg_matrix, "ramachandran_avg", header=header)

        fig, ax = plt.subplots(figsize=figsize)
        cmap = plt.get_cmap("viridis")
        norm = plt.Normalize(vmin=res_indices.min(), vmax=res_indices.max()) if len(res_indices) else None

        for i, res in enumerate(res_indices):
            color = cmap(norm(res)) if norm is not None else "blue"
            ax.errorbar(
                x[i],
                y[i],
                xerr=None if phi_std is None else phi_std[i],
                yerr=None if psi_std is None else psi_std[i],
                fmt="o",
                color=color,
                ecolor=color,
                capsize=3,
                **{k: v for k, v in kwargs.items() if k not in ["capsize"]},
            )

        ax.set_title(title)
        ax.set_xlabel(f"Phi ({self.units})")
        ax.set_ylabel(f"Psi ({self.units})")
        ax.grid(True, alpha=0.3)
        limit = 180.0 if self.units == "degrees" else np.pi
        pad = 20.0 if self.units == "degrees" else (np.pi / 9.0)
        ax.set_xlim(-(limit + pad), limit + pad)
        ax.set_ylim(-(limit + pad), limit + pad)
        if len(res_indices):
            mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
            cbar = plt.colorbar(mappable, ax=ax)
            cbar.set_label("Residue Index")

        fig.tight_layout()
        outpath = self._save_plot(fig, "ramachandran")
        plt.close(fig)

        # If residues are selected, also generate per-residue frame-level plots
        if res_list:
            try:
                _, phi_angles = md.compute_phi(self.traj)
                _, psi_angles = md.compute_psi(self.traj)
                if phi_angles.size and psi_angles.size:
                    phi_angles = phi_angles[:, res_list]
                    psi_angles = psi_angles[:, res_list]

                    if self.units == "degrees":
                        phi_angles = np.degrees(phi_angles)
                        psi_angles = np.degrees(psi_angles)

                    per_residue: Dict[int, Path] = {}
                    for idx, res in enumerate(res_list):
                        fig_res, ax_res = plt.subplots(figsize=figsize)
                        ax_res.scatter(
                            phi_angles[:, idx],
                            psi_angles[:, idx],
                            color=kwargs.get("color", "blue"),
                            alpha=kwargs.get("alpha", 0.7),
                        )
                        ax_res.set_title(f"Ramachandran Plot — Residue {res}")
                        ax_res.set_xlabel(f"Phi ({self.units})")
                        ax_res.set_ylabel(f"Psi ({self.units})")
                        ax_res.grid(True, alpha=0.3)
                        ax_res.set_xlim(-(limit + pad), limit + pad)
                        ax_res.set_ylim(-(limit + pad), limit + pad)
                        fig_res.tight_layout()

                        frame_matrix = np.column_stack([phi_angles[:, idx], psi_angles[:, idx]])
                        frame_header = f"phi_{self.units} psi_{self.units}"
                        self._save_data(
                            frame_matrix,
                            f"ramachandran_res{res}",
                            header=frame_header,
                        )
                        per_path = self._save_plot(
                            fig_res,
                            "ramachandran",
                            filename=f"ramachandran_res{res}",
                        )
                        plt.close(fig_res)
                        per_residue[res] = per_path

                    self.results["ramachandran_per_residue"] = per_residue
            except Exception as exc:
                logger.warning("Per-residue Ramachandran plots failed: %s", exc)

        return outpath

    def plot(self, **kwargs) -> str:
        # Default to Ramachandran if available
        if "phi_avg" in self.results and "psi_avg" in self.results:
            return self.plot_ramachandran(**kwargs)
        else:
            raise AnalysisError("No suitable data for plotting")