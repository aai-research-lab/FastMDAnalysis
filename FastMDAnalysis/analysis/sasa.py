"""
SASA Analysis Module

Computes the Solvent Accessible Surface Area (SASA) for an MD trajectory.
This module uses MDTraj's Shrake–Rupley algorithm with a specified probe radius to compute SASA.
It calculates three types of SASA data:
  1. Total SASA vs. Frame: The sum of per-residue SASA values for each frame.
  2. Per-Residue SASA vs. Frame: A 2D array showing each residue’s SASA across frames (heatmap).
  3. Average per-Residue SASA: The average SASA value for each residue over all frames.

If an atom selection string is provided (via the constructor or overridden in method calls),
only those atoms are used for the analysis. Otherwise, the entire trajectory is used.

The computed data are saved to files and default plots are generated automatically.
Users may later replot individual outputs with customizable options.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib import Path
from .base import BaseAnalysis, AnalysisError

class SASAAnalysis(BaseAnalysis):
    def __init__(self, trajectory, probe_radius: float = 0.14, atoms: str = None, **kwargs):
        """
        Initialize SASA analysis.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        probe_radius : float, optional
            The probe radius (in nm) for the Shrake–Rupley algorithm (default: 0.14).
        atoms : str, optional
            An MDTraj atom selection string (e.g., "protein", "protein and name CA") specifying which atoms to use.
            If None, all atoms are used.
        kwargs : dict
            Additional keyword arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        self.probe_radius = probe_radius
        self.atoms = atoms
        self.data = None

    def run(self) -> dict:
        """
        Compute SASA using MDTraj's Shrake–Rupley algorithm.

        If an atom selection is provided, the trajectory is subset accordingly.
        Computes:
          - total_sasa: a 1D array with the sum of SASA for each frame.
          - residue_sasa: a 2D array (n_frames x n_residues) of per-residue SASA.
          - average_residue_sasa: a 1D array with the average SASA per residue over all frames.

        Saves all computed datasets to files and automatically generates default plots.

        Returns
        -------
        dict
            Dictionary with keys "total_sasa", "residue_sasa", and "average_residue_sasa".
        """
        try:
            # Subset trajectory if a specific atom selection is given.
            if self.atoms is not None:
                atom_indices = self.traj.topology.select(self.atoms)
                if atom_indices is None or len(atom_indices) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(atom_indices)
            else:
                subtraj = self.traj

            # Compute per-residue SASA (returns an array of shape (n_frames, n_residues)).
            residue_sasa = md.shrake_rupley(subtraj, probe_radius=self.probe_radius, mode='residue')
            # Total SASA per frame is the sum over residues.
            total_sasa = residue_sasa.sum(axis=1)
            # Average per-residue SASA is computed over frames.
            average_residue_sasa = residue_sasa.mean(axis=0)

            # Store the computed data.
            self.data = {
                "total_sasa": total_sasa,                       # shape: (n_frames,)
                "residue_sasa": residue_sasa,                   # shape: (n_frames, n_residues)
                "average_residue_sasa": average_residue_sasa    # shape: (n_residues,)
            }
            self.results = self.data

            # Save the computed data to files.
            self._save_data(total_sasa.reshape(-1, 1), "total_sasa")
            self._save_data(residue_sasa, "residue_sasa")
            self._save_data(average_residue_sasa.reshape(-1, 1), "average_residue_sasa")

            # Generate default plots for all SASA outputs.
            self.plot()
            return self.data
        except Exception as e:
            raise AnalysisError(f"SASA analysis failed: {e}")

    def plot(self, data=None, option="all", **kwargs):
        """
        Replot SASA analysis outputs with customizable options.

        Parameters
        ----------
        data : dict, optional
            A dictionary containing SASA data (with keys "total_sasa", "residue_sasa", "average_residue_sasa").
            If None, the internal self.data is used.
        option : str, optional
            Which plot to generate. Options:
                "total"   - Replot total SASA vs. frame.
                "residue" - Replot per-residue SASA vs. frame heatmap.
                "average" - Replot average per-residue SASA as a bar plot.
                "all"     - Generate all plots. (Default)
        kwargs : dict
            Customizable plot options. For each type you can define:
              - For total SASA plot: title_total, xlabel_total, ylabel_total.
              - For residue SASA heatmap: title_residue, xlabel_residue, ylabel_residue, cmap (colormap).
              - For average SASA plot: title_avg, xlabel_avg, ylabel_avg.

        Returns
        -------
        dict or str
            If option is "all", returns a dictionary with keys "total", "residue", "average"
            mapping to file paths for each plot. If one option is specified, returns the file path as a string.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No SASA data available. Please run the analysis first.")

        plots = {}

        if option in ["all", "total"]:
            plots["total"] = self._plot_total_sasa(data["total_sasa"], **kwargs)
        if option in ["all", "residue"]:
            plots["residue"] = self._plot_residue_sasa(data["residue_sasa"], **kwargs)
        if option in ["all", "average"]:
            plots["average"] = self._plot_average_residue_sasa(data["average_residue_sasa"], **kwargs)
        
        if option == "all":
            return plots
        else:
            # option is a string like "total", "residue", or "average"
            return plots[option]

    def _plot_total_sasa(self, total_sasa, **kwargs):
        """Generate a plot for Total SASA vs. Frame."""
        frames = np.arange(self.traj.n_frames)
        title = kwargs.get("title_total", "Total SASA vs. Frame")
        xlabel = kwargs.get("xlabel_total", "Frame")
        ylabel = kwargs.get("ylabel_total", "Total SASA (nm²)")
        
        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, total_sasa, marker="o", linestyle="-")
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, "total_sasa")
        plt.close(fig)
        return plot_path

    def _plot_residue_sasa(self, residue_sasa, **kwargs):
        """Generate a heatmap for Per-Residue SASA vs. Frame."""
        title = kwargs.get("title_residue", "Per-Residue SASA vs. Frame")
        xlabel = kwargs.get("xlabel_residue", "Frame")
        ylabel = kwargs.get("ylabel_residue", "Residue Index")
        cmap = kwargs.get("cmap", "viridis")
        
        fig = plt.figure(figsize=(12, 8))
        # Transpose so that rows become residues and columns become frames.
        im = plt.imshow(residue_sasa.T, aspect="auto", interpolation="none", cmap=cmap)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        cbar = plt.colorbar(im)
        cbar.set_label("SASA (nm²)")
        # Set y-axis tick labels to whole numbers starting from 1.
        n_residues = residue_sasa.shape[1]
        plt.yticks(ticks=np.arange(n_residues), labels=np.arange(1, n_residues + 1))
        plot_path = self._save_plot(fig, "residue_sasa")
        plt.close(fig)
        return plot_path

    def _plot_average_residue_sasa(self, average_sasa, **kwargs):
        """Generate a bar plot for Average per-Residue SASA."""
        n_residues = average_sasa.shape[0]
        residues = np.arange(1, n_residues + 1)
        title = kwargs.get("title_avg", "Average per-Residue SASA")
        xlabel = kwargs.get("xlabel_avg", "Residue")
        ylabel = kwargs.get("ylabel_avg", "Average SASA (nm²)")
        
        fig = plt.figure(figsize=(12, 6))
        plt.bar(residues, average_sasa.flatten())
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(residues)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, "average_residue_sasa")
        plt.close(fig)
        return plot_path

