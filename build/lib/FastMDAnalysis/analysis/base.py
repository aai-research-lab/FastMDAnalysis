"""
Base analysis class and custom exception for FastMDAnalysis.
"""

import os
from pathlib import Path
import numpy as np


class AnalysisError(Exception):
    """Custom exception class for analysis errors."""
    pass


class BaseAnalysis:
    """
    BaseAnalysis class that provides common functionality for all analysis types.
    """
    def __init__(self, trajectory, output=None, **kwargs):
        """
        Initialize the analysis.

        Args:
            trajectory: An MDTraj Trajectory object.
            output (str): Optional output directory.
            kwargs: Additional keyword arguments.
        """
        self.traj = trajectory
        self.output = output or self.__class__.__name__.replace("Analysis", "").lower() + "_output"
        self.results = {}
        self.data = None
        # Create the output directory if it doesn't exist.
        self._create_output_dir()

    def _create_output_dir(self):
        """Creates the output directory if it does not exist."""
        self.outdir = Path(self.output)
        self.outdir.mkdir(parents=True, exist_ok=True)

    def _save_plot(self, fig, name: str) -> Path:
        """
        Save a plot to the output directory.

        Args:
            fig: A matplotlib figure object.
            name (str): Base name for the plot file.

        Returns:
            Path: The file path to the saved plot.
        """
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        return plot_path

    def _save_data(self, data, filename: str) -> Path:
        """
        Save numpy data array to a .dat file.

        Args:
            data: Data to be saved (numpy array or similar).
            filename (str): Base name for the data file.

        Returns:
            Path: The file path to the saved data.
        """
        data_path = self.outdir / f"{filename}.dat"
        if isinstance(data, np.ndarray):
            # Save with a header for clarity.
            header = " ".join([f"col{i}" for i in range(data.shape[1])]) if data.ndim == 2 else "data"
            np.savetxt(data_path, data, header=header)
        else:
            with open(data_path, "w") as f:
                f.write(str(data))
        return data_path

    def run(self):
        """
        Method to perform the analysis.
        To be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement the run() method.")

    def plot(self):
        """
        Method to produce plots.
        To be implemented by subclasses if plotting is available.
        """
        raise NotImplementedError("Subclasses must implement the plot() method.")

