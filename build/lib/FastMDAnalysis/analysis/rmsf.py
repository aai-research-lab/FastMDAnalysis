"""
RMSF Analysis Module
Calculates the Root Mean Square Fluctuation (RMSF) for each atom.
The user can choose the set of atoms to be used in the calculation via the "selection" parameter.
Options include:
  - "all" or "all atoms": all atoms in the trajectory.
  - "c-alpha" (or variants): protein C-alpha atoms.
  - "backbone": protein backbone atoms.
  - "heavy" (or "heavy atoms"): non-hydrogen atoms.
  - Any valid MDTraj selection string.
  
By default, when run (via fastmda.rmsf), the analysis automatically plots and saves the default plot.
Users may re-plot with customized options by calling analysis.plot(data, **kwargs).
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError


class RMSFAnalysis(BaseAnalysis):
    def __init__(self, trajectory, selection: str = "all", **kwargs):
        """
        Initialize RMSF analysis.

        Args:
            trajectory: Input MD trajectory.
            selection (str): Atom selection for the RMSF calculation. Options include "all", "c-alpha",
                             "backbone", "heavy", or any valid MDTraj selection string. Default is "all".
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        self.selection = selection
        self.data = None

    def _get_selection_indices(self):
        """Determine the atom indices based on the selection string."""
        sel = self.selection.lower()
        if sel in ["all", "all atoms"]:
            indices = list(range(self.traj.n_atoms))
        elif sel in ["c-alpha", "calpha", "c_alpha"]:
            indices = self.traj.topology.select("protein and name CA")
        elif sel == "backbone":
            indices = self.traj.topology.select("protein and backbone")
        elif sel in ["heavy", "heavy atoms"]:
            indices = self.traj.topology.select("not element H")
        else:
            # Assume it's a custom MDTraj selection string.
            indices = self.traj.topology.select(self.selection)
            
        if indices is None or len(indices) == 0:
            raise AnalysisError(f"No atoms selected with the selection: {self.selection}")
        return indices

    def run(self):
        """
        Compute the RMSF.

        The RMSF is calculated as the square root of the average squared deviation of atom positions
        relative to the mean structure computed for the selected atoms.
        
        The method automatically generates a default plot and saves it.
        """
        try:
            # Get selected atom indices.
            indices = self._get_selection_indices()
            # Create a sub-trajectory containing only the selected atoms.
            subtraj = self.traj.atom_slice(indices)
            
            # Compute the average coordinates across all frames for the selected atoms.
            avg_xyz = np.mean(subtraj.xyz, axis=0, keepdims=True)
            # Create a single-frame trajectory from the average coordinates.
            ref = md.Trajectory(avg_xyz, subtraj.topology)
            # Calculate RMSF relative to the reference (average structure).
            rmsf_values = md.rmsf(subtraj, ref)
            self.data = rmsf_values.reshape(-1, 1)
            self.results = {"rmsf": self.data, "selected_indices": indices}
            # Save the RMSF data.
            self._save_data(self.data, "rmsf")
            # Automatically generate and save the default plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"RMSF analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Plot RMSF per atom for the selected atoms.

        If no data is provided, self.data is used.

        Optional keyword arguments (kwargs) allow customization similar to matplotlib:
          - title (str): Title of the plot.
          - xlabel (str): Label for the x-axis.
          - ylabel (str): Label for the y-axis.
          - color (str or list): Color(s) for the plot elements.
          - linestyle (str): Line style (if applicable).
          - filename (str): Base filename for saving the plot (default: "rmsf").

        Returns:
            Path: The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RMSF data available to plot. Please run the analysis first.")

        # Set default plot parameters.
        title = kwargs.get('title', f"RMSF per Atom (Selection: {self.selection})")
        xlabel = kwargs.get('xlabel', "Atom Index (in selected subset)")
        ylabel = kwargs.get('ylabel', "RMSF (nm)")
        color = kwargs.get('color', None)  # Use matplotlib default if not provided.
        # Note: 'linestyle' is not directly used in bar plots.
        filename = kwargs.get('filename', "rmsf")

        # Create the plot.
        atom_indices = np.arange(len(data))
        fig = plt.figure(figsize=(10, 6))
        plt.bar(atom_indices, data.flatten(), width=0.8, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)

        # Save the figure using BaseAnalysis method.
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

