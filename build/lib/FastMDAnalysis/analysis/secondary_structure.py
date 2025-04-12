"""
Secondary Structure Analysis Module
Computes secondary structure assignments for each frame using DSSP,
saves the assignments to a text file, automatically generates a heatmap plot,
and writes a README file (secondary_structure_README.md) explaining the letter codes.
Users can replot with custom options via keyword arguments.
"""

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from .base import BaseAnalysis, AnalysisError


class SecondaryStructureAnalysis(BaseAnalysis):
    def __init__(self, trajectory, **kwargs):
        """
        Initialize Secondary Structure analysis.

        Args:
            trajectory: Input MD trajectory.
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        self.data = None

    def _generate_readme(self):
        """
        Generates a secondary_structure_README.md file in the output directory
        that explains the secondary structure letter codes.
        """
        content = (
            "# Secondary Structure Letter Codes\n\n"
            "This document explains the secondary structure codes used by DSSP and displayed in the \n"
            "FastMDAnalysis secondary structure heatmap.\n\n"
            "| Code | Description                                      |\n"
            "|------|--------------------------------------------------|\n"
            "| H    | Alpha helix                                      |\n"
            "| B    | Isolated beta-bridge                             |\n"
            "| E    | Extended strand (beta sheet)                     |\n"
            "| G    | 3-10 helix                                       |\n"
            "| I    | Pi helix                                         |\n"
            "| T    | Turn                                             |\n"
            "| S    | Bend                                             |\n"
            "| C    | Coil / Loop (no regular secondary structure)     |\n\n"
            "The FastMDAnalysis secondary structure module maps these codes to numeric values as follows:\n\n"
            "- **C** or blank → 0 (Coil / Loop)\n"
            "- **H** → 1 (Alpha helix)\n"
            "- **B** → 2 (Isolated beta-bridge)\n"
            "- **E** → 3 (Extended strand)\n"
            "- **G** → 4 (3-10 helix)\n"
            "- **I** → 5 (Pi helix)\n"
            "- **T** → 6 (Turn)\n"
            "- **S** → 7 (Bend)\n"
        )
        readme_path = self.outdir / "secondary_structure_README.md"
        with open(readme_path, "w") as f:
            f.write(content)
        return readme_path

    def run(self):
        """
        Compute secondary structure assignments using DSSP.
        Save the assignments to a text file, automatically generate a heatmap plot,
        and write the secondary_structure_README.md file.
        """
        try:
            dssp = md.compute_dssp(self.traj)
            self.data = dssp  # shape: (n_frames, n_residues)
            self.results = {"secondary_structure": self.data}
            # Save secondary structure assignments in a text file.
            data_path = self.outdir / "secondary_structure.dat"
            with open(data_path, "w") as f:
                for frame_idx, ss in enumerate(dssp):
                    f.write(f"Frame {frame_idx}: {', '.join(ss)}\n")
            # Generate the README file.
            self._generate_readme()
            # Automatically generate and save the default heatmap plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"Secondary structure analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Plot a heatmap of secondary structure assignments over frames.
        The heatmap uses a discrete colormap with very distinct colors so that each
        secondary structure letter is clearly differentiable. The colorbar is labeled
        with the secondary structure letter codes.
        
        The y-axis (residue index) is set to whole numbers starting from 1.
        
        Optional keyword arguments include:
          - title (str): Plot title (default: "Secondary Structure Heatmap").
          - xlabel (str): X-axis label (default: "Frame").
          - ylabel (str): Y-axis label (default: "Residue Index").
          - filename (str): Base filename for saving the plot (default: "secondary_structure").
        
        Returns:
            Path: The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No secondary structure data to plot. Please run analysis first.")

        # Define the discrete mapping from secondary structure letters to integers.
        mapping = {"H": 1, "B": 2, "E": 3, "G": 4, "I": 5, "T": 6, "S": 7, "C": 0, " ": 0}
        # Convert the DSSP letter codes to a numeric array.
        numeric = np.array([[mapping.get(s, 0) for s in frame] for frame in data])
        
        title = kwargs.get('title', "Secondary Structure Heatmap")
        xlabel = kwargs.get('xlabel', "Frame")
        ylabel = kwargs.get('ylabel', "Residue Index")
        filename = kwargs.get('filename', "secondary_structure")
        
        # Create a discrete colormap with distinct, easily differentiated colors.
        colors = ['lightgray', 'red', 'orange', 'blue', 'green', 'magenta', 'yellow', 'cyan']
        cmap = ListedColormap(colors)
        boundaries = np.arange(-0.5, 8, 1)  # for 8 discrete levels (0 to 7)
        norm = BoundaryNorm(boundaries, cmap.N)
        
        fig = plt.figure(figsize=(12, 8))
        im = plt.imshow(numeric.T, aspect="auto", interpolation="none", cmap=cmap, norm=norm)
        cbar = plt.colorbar(im, ticks=np.arange(0, 8))
        # Fix colorbar tick labels to display the corresponding secondary structure letters.
        tick_labels = ["C", "H", "B", "E", "G", "I", "T", "S"]
        cbar.set_ticklabels(tick_labels)
        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        # Set the y-axis (residue index) ticks to be whole numbers starting from 1.
        n_residues = numeric.shape[1]
        plt.yticks(ticks=np.arange(n_residues), labels=np.arange(1, n_residues + 1))
        
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

