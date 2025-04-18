"""
Radius of Gyration (RG) Analysis Module

Calculates the radius of gyration for each frame of an MD trajectory.
This analysis accepts an optional atom selection string so that the RG can be computed on a subset
of atoms. If no atom selection is provided, the calculation is done on the entire trajectory.
The computed RG values are saved to a file and a default plot of RG vs. frame is automatically generated.
Users can replot the data with customizable plotting options.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

class RGAnalysis(BaseAnalysis):
    def __init__(self, trajectory, atoms: str = None, **kwargs):
        """
        Initialize RG analysis.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            An MDTraj atom selection string specifying which atoms to use for the analysis.
            For example, "protein and name CA". If None, all atoms are used.
        kwargs : dict
            Additional keyword arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        self.atoms = atoms
        self.data = None

    def run(self) -> dict:
        """
        Compute the radius of gyration for each frame.

        If an atom selection is provided, the trajectory is subset accordingly using atom_slice.
        The computed RG values are then saved and a default plot is generated.

        Returns
        -------
        dict
            A dictionary with the key "rg" mapping to a column vector of RG values.
        """
        try:
            # If an atom selection is specified, subset the trajectory.
            if self.atoms is not None:
                atom_indices = self.traj.topology.select(self.atoms)
                if atom_indices is None or len(atom_indices) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(atom_indices)
            else:
                subtraj = self.traj

            # Compute the radius of gyration for each frame using the (possibly) subset trajectory.
            rg_values = md.compute_rg(subtraj)
            self.data = rg_values.reshape(-1, 1)
            self.results = {"rg": self.data}

            # Save the computed data.
            self._save_data(self.data, "rg")
            # Automatically generate and save the default plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"Radius of gyration analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Generate a plot of radius of gyration versus frame number.

        Parameters
        ----------
        data : array-like, optional
            The RG data to plot. If None, uses the data computed by run().
        kwargs : dict
            Customizable matplotlib-style keyword arguments, including:
                - title: Plot title (default: "Radius of Gyration vs Frame").
                - xlabel: x-axis label (default: "Frame").
                - ylabel: y-axis label (default: "Radius of Gyration (nm)").
                - color: Line or marker color.
                - linestyle: Line style (default: "-").

        Returns
        -------
        pathlib.Path
            The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RG data available to plot. Please run analysis first.")

        frames = np.arange(len(data))
        title = kwargs.get("title", "Radius of Gyration vs Frame")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "Radius of Gyration (nm)")
        color = kwargs.get("color", None)
        linestyle = kwargs.get("linestyle", "-")

        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, data, marker="o", linestyle=linestyle, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, "rg")
        plt.close(fig)
        return plot_path

