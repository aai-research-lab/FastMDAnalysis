"""
Hydrogen Bonds Analysis Module
Detects hydrogen bonds in the trajectory using the Baker-Hubbard algorithm.
Counts the number of H-bonds per frame, saves the data,
and automatically generates a default plot.
Users can replot with custom options via keyword arguments.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError


class HBondsAnalysis(BaseAnalysis):
    def __init__(self, trajectory, **kwargs):
        """
        Initialize Hydrogen Bonds analysis.

        Args:
            trajectory: Input MD trajectory.
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        self.data = None

    def run(self):
        """
        Detect hydrogen bonds using Baker-Hubbard,
        count the number of H-bonds per frame, save the data,
        and automatically generate a default plot.
        """
        try:
            hbonds = md.baker_hubbard(self.traj, periodic=False)
            counts = np.zeros(self.traj.n_frames, dtype=int)
            for h in hbonds:
                frame = h[0]
                counts[frame] += 1
            self.data = counts.reshape(-1, 1)
            self.results = {"hbonds_counts": self.data, "raw_hbonds": hbonds}
            self._save_data(self.data, "hbonds_counts")
            # Automatically generate and save the default plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"Hydrogen bonds analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Plot the number of hydrogen bonds per frame.

        If no data is provided, uses self.data.

        Optional keyword arguments include:
          - title (str): Plot title (default: "Hydrogen Bonds per Frame").
          - xlabel (str): X-axis label (default: "Frame").
          - ylabel (str): Y-axis label (default: "Number of H-Bonds").
          - color (str or list): Color(s) for the plot.
          - linestyle (str): Line style (default: solid line).
          - filename (str): Base filename for the plot (default: "hbonds").

        Returns:
            Path: The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No hydrogen bonds data available to plot. Please run analysis first.")

        title = kwargs.get('title', "Hydrogen Bonds per Frame")
        xlabel = kwargs.get('xlabel', "Frame")
        ylabel = kwargs.get('ylabel', "Number of H-Bonds")
        color = kwargs.get('color', None)
        linestyle = kwargs.get('linestyle', '-')  # Solid line default.
        filename = kwargs.get('filename', "hbonds")
        
        frames = np.arange(len(data))
        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, data.flatten(), marker="o", linestyle=linestyle, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

