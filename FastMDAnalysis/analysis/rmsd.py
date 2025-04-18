"""
RMSD Analysis Module

Calculates the Root-Mean-Square Deviation (RMSD) of an MD trajectory relative to a reference frame.
This analysis supports an optional atom selection so that RMSD can be computed on a subset of atoms.
The module automatically saves the RMSD data and generates a default plot of RMSD vs. frame number.
Users can optionally replot the data with customized plotting options.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

class RMSDAnalysis(BaseAnalysis):
    def __init__(self, trajectory, ref_frame: int = 0, atoms: str = None, **kwargs):
        """
        Initialize RMSD analysis.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        ref_frame : int, optional
            The index of the reference frame to which the RMSD is computed (default: 0).
        atoms : str, optional
            MDTraj atom selection string specifying which atoms to use.
            If None, all atoms are used.
        kwargs : dict
            Additional keyword arguments passed to the BaseAnalysis class.
        """
        super().__init__(trajectory, **kwargs)
        self.ref_frame = ref_frame
        self.atoms = atoms  # Optional atom selection string.
        self.data = None

    def run(self) -> dict:
        """
        Compute RMSD for each frame relative to the reference frame.

        If an atom selection is provided, only those atoms are used for the calculation.

        Returns
        -------
        dict
            A dictionary containing the computed RMSD data.
        """
        try:
            # Extract the reference frame.
            ref = self.traj[self.ref_frame]
            # Determine the atom indices, if an atom selection is provided.
            if self.atoms is not None:
                atom_indices = self.traj.topology.select(self.atoms)
                if atom_indices is None or len(atom_indices) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
            else:
                atom_indices = None

            # Compute RMSD values relative to the reference.
            rmsd_values = md.rmsd(self.traj, ref, atom_indices=atom_indices)
            self.data = rmsd_values.reshape(-1, 1)
            self.results = {"rmsd": self.data}

            # Save the computed data to a file.
            self._save_data(self.data, "rmsd")
            # Automatically generate the default RMSD plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"RMSD analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Generate a plot of RMSD versus frame number.

        Parameters
        ----------
        data : array-like, optional
            RMSD data to plot. If None, uses the data computed by run().
        kwargs : dict
            Matplotlib-style keyword arguments to customize the plot, e.g.:
            - title: plot title (default: "RMSD vs Frame (Reference Frame: ref_frame)").
            - xlabel: x-axis label (default: "Frame").
            - ylabel: y-axis label (default: "RMSD (nm)").
            - color: line color.
            - linestyle: line style.

        Returns
        -------
        Path
            The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No RMSD data available to plot. Please run the analysis first.")

        frames = np.arange(len(data))
        title = kwargs.get("title", f"RMSD vs Frame (Reference Frame: {self.ref_frame})")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "RMSD (nm)")
        color = kwargs.get("color", None)
        linestyle = kwargs.get("linestyle", "-")

        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, data, marker="o", linestyle=linestyle, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, "rmsd")
        plt.close(fig)
        return plot_path

