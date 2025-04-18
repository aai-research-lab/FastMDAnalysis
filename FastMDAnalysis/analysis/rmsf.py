"""
RMSF Analysis Module

Calculates the Root-Mean-Square Fluctuation (RMSF) for each atom in an MD trajectory.
If an atom selection is provided, only those atoms are analyzed; otherwise, all atoms are used.
The analysis computes the fluctuations relative to the average structure, saves the computed data,
and automatically generates a default bar plot of RMSF per atom.
Users may replot the data with customizable plotting options.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

class RMSFAnalysis(BaseAnalysis):
    def __init__(self, trajectory, atoms: str = None, **kwargs):
        """
        Initialize RMSF analysis.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            MDTraj atom selection string specifying which atoms to use.
            If None, all atoms in the trajectory are used.
        kwargs : dict
            Additional keyword arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        self.atoms = atoms
        self.data = None

    def run(self) -> dict:
        """
        Compute RMSF for each atom relative to the average structure.

        If an atom selection is provided, the analysis is performed on the subtrajectory defined by that selection.
        The method calculates the per-atom RMSF, saves the resulting data, and generates a plot.

        Returns
        -------
        dict
            A dictionary containing the computed RMSF data.
        """
        try:
            # If an atom selection is specified, subset the trajectory accordingly.
            if self.atoms is not None:
                atom_indices = self.traj.topology.select(self.atoms)
                if atom_indices is None or len(atom_indices) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(atom_indices)
            else:
                subtraj = self.traj

            # Compute the average structure (mean coordinates) for the selected atoms.
            avg_xyz = np.mean(subtraj.xyz, axis=0, keepdims=True)
            ref = md.Trajectory(avg_xyz, subtraj.topology)

            # Compute per-atom RMSF values relative to the average structure.
            rmsf_values = md.rmsf(subtraj, ref)

            # Reshape the result to a column vector.
            self.data = rmsf_values.reshape(-1, 1)
            self.results = {"rmsf": self.data}

            # Save the computed RMSF data to a file.
            self._save_data(self.data, "rmsf")

            # Automatically generate and save the default RMSF plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"RMSF analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Generate a bar plot of RMSF per atom.

        Parameters
        ----------
        data : array-like, optional
            The RMSF data to plot. If None, uses the data computed by run().
        kwargs : dict
            Customizable matplotlib keyword arguments to style the plot. For example:
                - title: Plot title (default: "RMSF per Atom").
                - xlabel: x-axis label (default: "Atom Index").
                - ylabel: y-axis label (default: "RMSF (nm)").
                - color: Bar color.

        Returns
        -------
        Path
            The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RMSF data available to plot. Please run the analysis first.")

        # Prepare the x-axis (atom index)
        x = np.arange(len(data))
        title = kwargs.get("title", "RMSF per Atom")
        xlabel = kwargs.get("xlabel", "Atom Index")
        ylabel = kwargs.get("ylabel", "RMSF (nm)")
        color = kwargs.get("color", None)

        fig = plt.figure(figsize=(10, 6))
        plt.bar(x, data.flatten(), color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, "rmsf")
        plt.close(fig)
        return plot_path

