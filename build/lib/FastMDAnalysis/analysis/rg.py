"""
Radius of Gyration Analysis Module
Calculates the radius of gyration for each frame over the trajectory.
By default, when the analysis is run via fastmda.rg(...), the computed data 
are automatically saved and a default plot is generated.
Users can optionally replot the data with custom matplotlib-style options.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError


class RGAnalysis(BaseAnalysis):
    def __init__(self, trajectory, **kwargs):
        """
        Initialize Radius of Gyration analysis.

        Args:
            trajectory: Input MD trajectory.
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        self.data = None

    def run(self):
        """
        Compute the radius of gyration (Rg) for each frame.

        The computed Rg data are saved and a default plot is automatically generated.
        """
        try:
            # Compute the radius of gyration using MDTraj function.
            rg_values = md.compute_rg(self.traj)
            self.data = rg_values.reshape(-1, 1)
            self.results = {"rg": self.data}
            # Save the Rg data to file.
            self._save_data(self.data, "rg")
            # Automatically generate and save the default Rg plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"Radius of gyration analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Plot Radius of Gyration vs Frame number using the computed data.

        If no data is provided, self.data is used.

        Optional keyword arguments allow customization similar to matplotlib:
            - title (str): Plot title (default: "Radius of Gyration vs Frame").
            - xlabel (str): Label for the x-axis (default: "Frame").
            - ylabel (str): Label for the y-axis (default: "Radius of Gyration (nm)").
            - color (str or list): Color for the plot elements (default: matplotlib's default).
            - linestyle (str): Line style (default: solid line, "-").
            - filename (str): Base filename for saving the plot (default: "rg").

        Returns:
            Path: The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RG data available to plot. Please run the analysis first.")
        
        # Define default plotting parameters.
        title = kwargs.get('title', "Radius of Gyration vs Frame")
        xlabel = kwargs.get('xlabel', "Frame")
        ylabel = kwargs.get('ylabel', "Radius of Gyration (nm)")
        color = kwargs.get('color', None)  # Use default matplotlib colors if not specified.
        linestyle = kwargs.get('linestyle', '-')  # Solid line by default.
        filename = kwargs.get('filename', "rg")
        
        # Create the plot.
        frames = np.arange(len(data))
        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, data, marker="o", linestyle=linestyle, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        
        # Save the plot using the BaseAnalysis method.
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

