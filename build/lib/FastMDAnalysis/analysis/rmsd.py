"""
RMSD Analysis Module
Calculates the Root Mean Square Deviation (RMSD) of each frame in a trajectory
relative to a reference frame. By default, the reference frame is the frame at index 0.
The analysis automatically saves the computed data and produces a default plot.
Users can also replot the data with customized matplotlib arguments.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError


class RMSDAnalysis(BaseAnalysis):
    def __init__(self, trajectory, ref_frame: int = 0, **kwargs):
        """
        Initialize RMSD analysis.

        Args:
            trajectory: Input MD trajectory.
            ref_frame (int): The index of the reference frame. Default is 0.
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        self.ref_frame = ref_frame
        self.data = None

    def run(self):
        """
        Compute RMSD for each frame relative to the reference frame.
        
        The RMSD is computed using MDTraj's rmsd function. The computed data are saved
        in a file and a default plot is generated automatically.
        """
        try:
            # Extract the reference structure as a single-frame trajectory.
            ref = self.traj[self.ref_frame]
            # Compute RMSD values relative to the reference structure.
            rmsd_values = md.rmsd(self.traj, ref)
            self.data = rmsd_values.reshape(-1, 1)
            self.results = {"rmsd": self.data}
            # Save the computed RMSD data.
            self._save_data(self.data, "rmsd")
            # Automatically generate and save the default plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"RMSD analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Plot RMSD vs Frame number for the computed RMSD data.

        If no data is provided, the default self.data is used.
        Customization options (via keyword arguments) include:
          - title (str): Plot title (default: "RMSD vs Frame (Ref: {ref_frame})").
          - xlabel (str): X-axis label (default: "Frame").
          - ylabel (str): Y-axis label (default: "RMSD (nm)").
          - color (str or list): Color(s) for the plot (default: matplotlib's default).
          - linestyle (str): Line style for the plot (default: solid line).
          - filename (str): Base filename for saving the plot (default: "rmsd").

        Returns:
            Path: The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RMSD data available to plot. Please run the analysis first.")

        # Define default customization parameters.
        title = kwargs.get('title', f"RMSD vs Frame (Ref: {self.ref_frame})")
        xlabel = kwargs.get('xlabel', "Frame")
        ylabel = kwargs.get('ylabel', "RMSD (nm)")
        color = kwargs.get('color', None)  # Use default matplotlib color if not specified.
        linestyle = kwargs.get('linestyle', '-')  # Solid line by default.
        filename = kwargs.get('filename', "rmsd")

        # Create the plot.
        frames = np.arange(len(data))
        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, data, marker="o", linestyle=linestyle, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)

        # Save the plot in the output directory using the BaseAnalysis save method.
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

