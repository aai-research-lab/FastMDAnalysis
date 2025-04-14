"""
SS Analysis Module
Computes secondary structure assignments for each frame using DSSP.
Saves the assignments to a text file, automatically generates a heatmap plot,
and writes an ss_README.md file that explains the secondary structure (SS) letter codes.
The heatmap uses a discrete colormap with very distinct colors,
and the residue index axis is labeled with whole numbers starting from 1.
Users can replot with custom options via keyword arguments.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib import Path
from .base import BaseAnalysis, AnalysisError

class SSAnalysis(BaseAnalysis):
    def __init__(self, trajectory, **kwargs):
        """
        Initialize SS analysis.

        Args:
            trajectory: Input MD trajectory.
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        self.data = None

    def _generate_readme(self):
        """
        Generates an ss_README.md file in the output directory explaining the SS letter codes.
        """
        content = (
            "# Secondary Structure (SS) Letter Codes\n\n"
            "This document explains the secondary structure codes used by DSSP and displayed in the \n"
            "FastMDAnalysis SS heatmap.\n\n"
            "| Code | Description                                      |\n"
            "|------|--------------------------------------------------|\n"
            "| H    | Alpha helix                                      |\n"
            "| B    | Isolated beta-bridge                             |\n"
            "| E    | Extended strand (beta sheet)                     |\n"
            "| G    | 3-10 helix                                       |\n"
            "| I    | Pi helix                                         |\n"
            "| T    | Turn                                             |\n"
            "| S    | Bend                                             |\n"
            "| C or (space) | Coil / Loop (no regular secondary structure) |\n"
        )
        readme_path = self.outdir / "ss_README.md"
        with open(readme_path, "w") as f:
            f.write(content)
        return readme_path

    def run(self) -> dict:
        """
        Compute secondary structure assignments using DSSP.
        Save the assignments to a text file, automatically generate a heatmap plot,
        and write the ss_README.md file.
        """
        try:
            dssp = md.compute_dssp(self.traj)
            self.data = dssp  # shape: (n_frames, n_residues)
            self.results = {"ss_data": self.data}
            # Save the secondary structure assignments.
            data_path = self.outdir / "ss.dat"
            with open(data_path, "w") as f:
                for frame_idx, ss in enumerate(dssp):
                    f.write(f"Frame {frame_idx}: {', '.join(ss)}\n")
            # Generate the README file.
            self._generate_readme()
            # Automatically generate and save the default heatmap plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"SS analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Plot a heatmap of SS assignments over frames. The heatmap uses a discrete colormap
        with distinct colors so that each secondary structure letter can be easily differentiated.
        The colorbar tick labels show the SS letter codes, and the y-axis (residue index)
        is labeled with whole numbers starting from 1.
        
        Optional keyword arguments include:
          - title (str): Plot title (default: "SS Heatmap").
          - xlabel (str): X-axis label (default: "Frame").
          - ylabel (str): Y-axis label (default: "Residue Index").
          - cmap (str): Colormap (default: a custom discrete colormap).
          - filename (str): Base filename for saving the plot (default: "ss").
        
        Returns:
            The file path (Path object) to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No SS data available to plot. Please run analysis first.")
        
        # Define a mapping from SS letters to numeric values.
        mapping = {"H": 1, "B": 2, "E": 3, "G": 4, "I": 5, "T": 6, "S": 7, "C": 0, " ": 0}
        # Convert the SS letters to a numeric array.
        numeric = np.array([[mapping.get(s, 0) for s in frame] for frame in data])
        
        # Create a discrete colormap with very distinct colors.
        # For 8 levels (0 to 7), we choose 8 distinct colors.
        from matplotlib.colors import ListedColormap, BoundaryNorm
        distinct_colors = ['#AAAAAA', '#FF0000', '#FFA500', '#0000FF', '#008000', '#FF00FF', '#FFFF00', '#00FFFF']
        cmap = ListedColormap(distinct_colors)
        boundaries = np.arange(-0.5, 8, 1)  # boundaries for 0 to 7
        norm = BoundaryNorm(boundaries, cmap.N)
        
        title = kwargs.get("title", "SS Heatmap")
        xlabel = kwargs.get("xlabel", "Frame")
        ylabel = kwargs.get("ylabel", "Residue Index")
        filename = kwargs.get("filename", "ss")
        
        fig = plt.figure(figsize=(12, 8))
        im = plt.imshow(numeric.T, aspect="auto", interpolation="none", cmap=cmap, norm=norm)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        cbar = plt.colorbar(im, ticks=np.arange(0, 8))
        # Set the colorbar tick labels to the corresponding SS letters.
        tick_labels = ["C", "H", "B", "E", "G", "I", "T", "S"]
        cbar.set_ticklabels(tick_labels)
        cbar.set_label("SS Code")
        # Ensure the y-axis shows residue indices as whole numbers, starting at 1.
        n_residues = numeric.shape[1]
        plt.yticks(ticks=np.arange(n_residues), labels=np.arange(1, n_residues + 1))
        
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

