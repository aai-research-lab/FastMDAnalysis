"""
Hydrogen Bonds Analysis Module

Detects hydrogen bonds in an MD trajectory using the Baker-Hubbard algorithm.
If an atom selection is provided (via the 'atoms' parameter), the trajectory is subset accordingly.
The analysis computes the number of hydrogen bonds for each frame, saves the resulting data,
and automatically generates a default plot of hydrogen bonds versus frame.
Users can later replot the data with customizable plotting options.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError

class HBondsAnalysis(BaseAnalysis):
    def __init__(self, trajectory, atoms: str = None, **kwargs):
        """
        Initialize Hydrogen Bonds analysis.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            An MDTraj atom selection string specifying which atoms to use.
            If provided, the trajectory will be subset using this selection.
            If None, all atoms in the trajectory are used.
        kwargs : dict
            Additional keyword arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        self.atoms = atoms
        self.data = None

    def run(self) -> dict:
        """
        Compute hydrogen bonds for each frame using the Baker-Hubbard algorithm.

        If an atom selection is provided, only the selected atoms are used for the calculation.
        The algorithm returns an array of detected hydrogen bonds (tuples),
        and we count the number of hydrogen bonds per frame.

        Returns
        -------
        dict
            A dictionary containing:
              "hbonds_counts": a column vector with the number of H-bonds per frame,
              "raw_hbonds": the raw hydrogen bond data as returned by md.baker_hubbard.
        """
        try:
            # If an atom selection is provided, subset the trajectory accordingly.
            if self.atoms is not None:
                atom_indices = self.traj.topology.select(self.atoms)
                if atom_indices is None or len(atom_indices) == 0:
                    raise AnalysisError(f"No atoms selected using the selection: '{self.atoms}'")
                subtraj = self.traj.atom_slice(atom_indices)
            else:
                subtraj = self.traj

            # Use the Baker-Hubbard algorithm to detect hydrogen bonds.
            hbonds = md.baker_hubbard(subtraj, periodic=False)
            counts = np.zeros(self.traj.n_frames, dtype=int)
            for bond in hbonds:
                # bond[0] is the frame index where the H-bond is detected.
                frame_index = bond[0]
                counts[frame_index] += 1

            self.data = counts.reshape(-1, 1)
            self.results = {"hbonds_counts": self.data, "raw_hbonds": hbonds}

            # Save the hydrogen bonds counts.
            self._save_data(self.data, "hbonds_counts")
            return self.results
        except Exception as e:
            raise AnalysisError(f"Hydrogen bonds analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Generate a plot of hydrogen bonds versus frame.

        Parameters
        ----------
        data : array-like, optional
            The hydrogen bond count data to plot. If None, uses the data computed by run().
        kwargs : dict
            Customizable matplotlib-style keyword arguments. For example:
                - title: Plot title (default: "Hydrogen Bonds per Frame").
                - xlabel: x-axis label (default: "Frame").
                - ylabel: y-axis label (default: "Number of H-Bonds").
                - color: Line or marker color.
                - linestyle: Line style (default: "-" for solid line).

        Returns
        -------
        Path
            The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No hydrogen bonds data available to plot. Please run analysis first.")

        frames = np.arange(len(data))
        title = kwargs.get("title", "Hydrogen Bonds per Frame")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "Number of H-Bonds")
        color = kwargs.get("color", None)
        linestyle = kwargs.get("linestyle", "-")

        fig = plt.figure(figsize=(10, 6))
        plt.plot(frames, data.flatten(), marker="o", linestyle=linestyle, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, "hbonds")
        plt.close(fig)
        return plot_path

