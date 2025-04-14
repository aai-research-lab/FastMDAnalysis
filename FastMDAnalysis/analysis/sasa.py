"""
SASA Analysis Module

Computes the Solvent Accessible Surface Area (SASA) for an MD trajectory.
This module uses the MDTraj Shrake-Rupley algorithm to compute SASA with a given probe radius.

It computes and saves three types of SASA data:
  1. Total SASA vs. Frame: Sum of SASA over all residues for each frame.
  2. Residue SASA vs. Frame: A heatmap of per-residue SASA over time.
  3. Average per-Residue SASA: The time-average SASA of each residue.

Default plots are generated automatically. Users may replot data with custom options
via the plot() method.

Usage Example (API):

    from FastMDAnalysis import FastMDAnalysis

    fastmda = FastMDAnalysis()
    sasa_analysis = fastmda.sasa("protein_traj.dcd", top="protein.pdb", output="sasa_output")
    data = sasa_analysis.data
    # Optionally, replot with custom settings:
    plots = sasa_analysis.plot(option="all", title_total="My Total SASA", xlabel_total="Frame",
                               ylabel_total="Total SASA (nm^2)",
                               title_residue="Per-Residue SASA Heatmap",
                               title_avg="Average SASA per Residue")
    print("Plots generated:", plots)
"""

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib import Path

from .base import BaseAnalysis, AnalysisError

class SASAAnalysis(BaseAnalysis):
    def __init__(self, trajectory, probe_radius: float = 0.14, **kwargs):
        """
        Initialize SASA analysis.

        Args:
            trajectory (mdtraj.Trajectory): The MD trajectory.
            probe_radius (float): The probe radius (in nm) to use in the Shrake–Rupley algorithm.
                                    Default is 0.14 nm.
            kwargs: Additional keyword arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        self.probe_radius = probe_radius
        self.results = {}  # Will hold computed SASA data.
        self.data = None

    def run(self) -> dict:
        """
        Compute SASA using MDTraj's Shrake–Rupley algorithm.
        
        Computes:
          - Per-residue SASA for each frame (using mode='residue').
          - Total SASA per frame by summing per-residue SASA.
          - Average per-residue SASA over all frames.
        
        Saves results to files and stores them in self.results.

        Returns:
            dict: A dictionary with keys:
                "total_sasa":  (n_frames, 1) array of total SASA per frame.
                "residue_sasa": (n_frames, n_residues) array of per-residue SASA.
                "average_residue_sasa": (n_residues, 1) array of average SASA per residue.
        """
        try:
            # Compute per-residue SASA (shape: n_frames x n_residues)
            residue_sasa = md.shrake_rupley(self.traj, probe_radius=self.probe_radius, mode='residue')
            # Total SASA is the sum over residues per frame.
            total_sasa = residue_sasa.sum(axis=1)
            # Average per-residue SASA is computed over frames.
            average_residue_sasa = residue_sasa.mean(axis=0)
            
            self.results = {
                "total_sasa": total_sasa,                      # shape: (n_frames,)
                "residue_sasa": residue_sasa,                  # shape: (n_frames, n_residues)
                "average_residue_sasa": average_residue_sasa   # shape: (n_residues,)
            }
            self.data = self.results

            # Save computed SASA data.
            self._save_data(total_sasa.reshape(-1, 1), "total_sasa")
            self._save_data(residue_sasa, "residue_sasa")
            self._save_data(average_residue_sasa.reshape(-1, 1), "average_residue_sasa")
            
            # Generate default plots.
            self._plot_total_sasa(total_sasa)
            self._plot_residue_sasa(residue_sasa)
            self._plot_average_residue_sasa(average_residue_sasa)
            
            return self.results
        except Exception as e:
            raise AnalysisError(f"SASA analysis failed: {e}")

    def plot(self, option=None, **kwargs):
        """
        Re-plot SASA data with customizable options.

        Args:
            option (str, optional): One of:
                "total" - replot total SASA vs. frame.
                "residue" - replot the heatmap for per-residue SASA vs. frame.
                "average" - replot average per-residue SASA (bar plot).
                If not provided (or option="all"), all three plots are re-generated.
            kwargs: Customizable matplotlib keyword arguments.
                For total SASA plot, you may use:
                    title_total, xlabel_total, ylabel_total.
                For residue SASA heatmap, you may use:
                    title_residue, xlabel_residue, ylabel_residue.
                For average SASA plot, you may use:
                    title_avg, xlabel_avg, ylabel_avg.
                    
        Returns:
            dict or str: If option is None or "all", returns a dictionary with keys
                         "total", "residue", and "average" mapping to plot file paths.
                         Otherwise, returns the plot file path (str) for the specified option.
        """
        if self.data is None:
            raise AnalysisError("No SASA data available. Please run analysis first.")

        plots = {}
        if option is None or option == "all":
            plots["total"] = self._plot_total_sasa(self.data["total_sasa"], **kwargs)
            plots["residue"] = self._plot_residue_sasa(self.data["residue_sasa"], **kwargs)
            plots["average"] = self._plot_average_residue_sasa(self.data["average_residue_sasa"], **kwargs)
            return plots
        elif option.lower() == "total":
            return self._plot_total_sasa(self.data["total_sasa"], **kwargs)
        elif option.lower() == "residue":
            return self._plot_residue_sasa(self.data["residue_sasa"], **kwargs)
        elif option.lower() == "average":
            return self._plot_average_residue_sasa(self.data["average_residue_sasa"], **kwargs)
        else:
            raise AnalysisError("Unknown plot option. Choose 'total', 'residue', 'average', or 'all'.")

    def _plot_total_sasa(self, total_sasa, **kwargs):
        """
        Plot Total SASA vs. Frame.
        """
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
        """
        Plot Per-Residue SASA vs. Frame as a heatmap.
        The heatmap shows residue index (starting at 1) on the y-axis and frame number on the x-axis.
        """
        title = kwargs.get("title_residue", "Per-Residue SASA vs. Frame")
        xlabel = kwargs.get("xlabel_residue", "Frame")
        ylabel = kwargs.get("ylabel_residue", "Residue Index")
        
        fig = plt.figure(figsize=(12, 8))
        # Transpose so that y-axis represents residues.
        im = plt.imshow(residue_sasa.T, aspect="auto", interpolation="none", cmap=kwargs.get("cmap", "viridis"))
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        cbar = plt.colorbar(im)
        cbar.set_label("SASA (nm²)")
        # Set y-axis tick labels to be whole numbers starting from 1.
        n_residues = residue_sasa.shape[1]
        plt.yticks(ticks=np.arange(n_residues), labels=np.arange(1, n_residues + 1))
        plot_path = self._save_plot(fig, "residue_sasa")
        plt.close(fig)
        return plot_path

    def _plot_average_residue_sasa(self, average_sasa, **kwargs):
        """
        Plot Average per-Residue SASA as a bar plot.
        X-axis: Residue index (starting at 1); Y-axis: Average SASA (nm²)
        """
        n_residues = average_sasa.shape[0]
        residues = np.arange(1, n_residues + 1)
        title = kwargs.get("title_avg", "Average per-Residue SASA")
        xlabel = kwargs.get("xlabel_avg", "Residue Index")
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

