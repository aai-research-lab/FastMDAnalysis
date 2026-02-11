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
import logging

import numpy as np
import mdtraj as md
from scipy.stats import circmean

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

            # Circular mean per residue
            n_residues = angles.shape[1]
            avg_angles = np.zeros(n_residues)
            for i in range(n_residues):
                res_angles = angles[:, i]
                # Remove NaN values (missing angles)
                valid_angles = res_angles[~np.isnan(res_angles)]
                if len(valid_angles) > 0:
                    avg_angles[i] = circmean(valid_angles, high=np.pi, low=-np.pi)
                else:
                    avg_angles[i] = np.nan

            # Convert units
            if self.units == "degrees":
                avg_angles = np.degrees(avg_angles)

            self.data = avg_angles.reshape(-1, 1)
            self.results = {"phi_avg": self.data}

            # Filter by residues if specified
            if self.residues is not None:
                if isinstance(self.residues, int):
                    res_list = [self.residues]
                else:
                    res_list = list(self.residues)
                filtered_data = self.data[res_list]
                self.results["phi_avg_filtered"] = filtered_data.reshape(-1, 1)

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
        if data is None:
            raise AnalysisError("No phi data available to plot.")

        y = np.asarray(data, dtype=float).flatten()
        n = len(y)
        x = np.arange(n)

        # Filter residues
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(x, residues)
            x = x[mask]
            y = y[mask]

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(x, y, color=color, **kwargs)

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        if ylabel is None:
            ylabel = f"Phi Angle ({self.units})"
        ax.set_ylabel(ylabel)

        # Highlight specific residues
        if highlight_residues is not None:
            if isinstance(highlight_residues, int):
                highlight_residues = [highlight_residues]
            for res in highlight_residues:
                if res in x:
                    idx = np.where(x == res)[0][0]
                    ax.bar([x[idx]], [y[idx]], color="red", alpha=0.7)

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

    def run(self) -> Dict[str, np.ndarray]:
        logger.info("Starting Psi analysis")
        try:
            indices, angles = md.compute_psi(self.traj)
            if angles.size == 0:
                raise AnalysisError("No psi angles found in trajectory")

            n_residues = angles.shape[1]
            avg_angles = np.zeros(n_residues)
            for i in range(n_residues):
                res_angles = angles[:, i]
                valid_angles = res_angles[~np.isnan(res_angles)]
                if len(valid_angles) > 0:
                    avg_angles[i] = circmean(valid_angles, high=np.pi, low=-np.pi)
                else:
                    avg_angles[i] = np.nan

            if self.units == "degrees":
                avg_angles = np.degrees(avg_angles)

            self.data = avg_angles.reshape(-1, 1)
            self.results = {"psi_avg": self.data}

            if self.residues is not None:
                res_list = [self.residues] if isinstance(self.residues, int) else list(self.residues)
                filtered_data = self.data[res_list]
                self.results["psi_avg_filtered"] = filtered_data.reshape(-1, 1)

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
        if kwargs["data"] is None:
            raise AnalysisError("No psi data available to plot.")

        y = np.asarray(kwargs["data"], dtype=float).flatten()
        n = len(y)
        x = np.arange(n)

        # Filter residues
        residues = kwargs.get("residues")
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(x, residues)
            x = x[mask]
            y = y[mask]

        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (12, 6)))
        ax.bar(x, y, color=kwargs.get("color", "green"), **{k: v for k, v in kwargs.items() if k not in ["data", "residues", "figsize", "title", "xlabel", "ylabel", "color", "highlight_residues"]})

        ax.set_title(kwargs.get("title", "Average Psi Angles per Residue"))
        ax.set_xlabel(kwargs.get("xlabel", "Residue Index"))
        ax.set_ylabel(kwargs.get("ylabel", f"Psi Angle ({self.units})"))

        # Highlight specific residues
        highlight_residues = kwargs.get("highlight_residues")
        if highlight_residues is not None:
            if isinstance(highlight_residues, int):
                highlight_residues = [highlight_residues]
            for res in highlight_residues:
                if res in x:
                    idx = np.where(x == res)[0][0]
                    ax.bar([x[idx]], [y[idx]], color="red", alpha=0.7)

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

    def run(self) -> Dict[str, np.ndarray]:
        logger.info("Starting Omega analysis")
        try:
            indices, angles = md.compute_omega(self.traj)
            if angles.size == 0:
                raise AnalysisError("No omega angles found in trajectory")

            n_residues = angles.shape[1]
            avg_angles = np.zeros(n_residues)
            for i in range(n_residues):
                res_angles = angles[:, i]
                valid_angles = res_angles[~np.isnan(res_angles)]
                if len(valid_angles) > 0:
                    avg_angles[i] = circmean(valid_angles, high=np.pi, low=-np.pi)
                else:
                    avg_angles[i] = np.nan

            if self.units == "degrees":
                avg_angles = np.degrees(avg_angles)

            self.data = avg_angles.reshape(-1, 1)
            self.results = {"omega_avg": self.data}

            if self.residues is not None:
                res_list = [self.residues] if isinstance(self.residues, int) else list(self.residues)
                filtered_data = self.data[res_list]
                self.results["omega_avg_filtered"] = filtered_data.reshape(-1, 1)

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
        if kwargs["data"] is None:
            raise AnalysisError("No omega data available to plot.")

        y = np.asarray(kwargs["data"], dtype=float).flatten()
        n = len(y)
        x = np.arange(n)

        # Filter residues
        residues = kwargs.get("residues")
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(x, residues)
            x = x[mask]
            y = y[mask]

        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (12, 6)))
        ax.bar(x, y, color=kwargs.get("color", "orange"), **{k: v for k, v in kwargs.items() if k not in ["data", "residues", "figsize", "title", "xlabel", "ylabel", "color", "highlight_residues"]})

        ax.set_title(kwargs.get("title", "Average Omega Angles per Residue"))
        ax.set_xlabel(kwargs.get("xlabel", "Residue Index"))
        ax.set_ylabel(kwargs.get("ylabel", f"Omega Angle ({self.units})"))

        # Highlight specific residues
        highlight_residues = kwargs.get("highlight_residues")
        if highlight_residues is not None:
            if isinstance(highlight_residues, int):
                highlight_residues = [highlight_residues]
            for res in highlight_residues:
                if res in x:
                    idx = np.where(x == res)[0][0]
                    ax.bar([x[idx]], [y[idx]], color="red", alpha=0.7)

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

        if phi_data is None or psi_data is None:
            raise AnalysisError("Phi and psi data required for Ramachandran plot")

        x = phi_data
        y = psi_data
        n = len(x)
        res_indices = np.arange(n)

        # Filter residues
        if residues is not None:
            if isinstance(residues, int):
                residues = [residues]
            mask = np.isin(res_indices, residues)
            x = x[mask]
            y = y[mask]
            res_indices = res_indices[mask]

        fig, ax = plt.subplots(figsize=figsize)
        scatter = ax.scatter(x, y, c=res_indices, cmap="viridis", **kwargs)
        ax.set_title(title)
        ax.set_xlabel(f"Phi ({self.units})")
        ax.set_ylabel(f"Psi ({self.units})")
        ax.grid(True, alpha=0.3)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Residue Index")

        fig.tight_layout()
        outpath = self._save_plot(fig, "ramachandran")
        plt.close(fig)

        return outpath

    def plot(self, **kwargs) -> str:
        # Default to Ramachandran if available
        if "phi_avg" in self.results and "psi_avg" in self.results:
            return self.plot_ramachandran(**kwargs)
        else:
            raise AnalysisError("No suitable data for plotting")