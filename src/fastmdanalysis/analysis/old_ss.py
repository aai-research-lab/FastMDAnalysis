"""
SS Analysis Module

Computes secondary structure assignments for each frame using DSSP.
Saves the SS assignments to a text file and automatically generates a heatmap plot.
The heatmap uses a discrete colormap with very distinct colors so that each SS letter is easily differentiated.
The colorbar tick labels display the SS letter codes.
The residue index axis is labeled with whole numbers starting at 1.
An ss_README.md file and companion PNG are also generated to explain the SS letter codes.

Usage:
    from fastmdanalysis import SSAnalysis
    analysis = SSAnalysis(trajectory, atoms="protein")
    analysis.run()         # Computes SS and generates default plots and README file.
    analysis.plot()        # Replot if needed with customization options.
"""
from __future__ import annotations

import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib import Path
from .base import BaseAnalysis, AnalysisError


SS_CODE_ROWS = [
    ("H", "Alpha helix"),
    ("B", "Isolated beta-bridge"),
    ("E", "Extended strand (beta sheet)"),
    ("G", "3-10 helix"),
    ("I", "Pi helix"),
    ("T", "Turn"),
    ("S", "Bend"),
    ("C or (space)", "Coil / Loop (no regular secondary structure)"),
]


class SSAnalysis(BaseAnalysis):
    def __init__(self, trajectory, atoms: str = None, **kwargs):
        """
        Initialize SS analysis.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            An MDTraj atom selection string to specify which atoms to consider.
            If None, all atoms are used.
        kwargs : dict
            Additional keyword arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        self.atoms = atoms
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

    def _generate_reference_image(self):
        """Create a PNG table summarizing the SS letter codes for quick reference."""
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.axis("off")
        fig.suptitle("Secondary Structure (SS) Letter Codes", fontsize=14, fontweight="bold", y=0.95)
        fig.text(
            0.5,
            0.88,
            "DSSP codes used in the FastMDAnalysis SS heatmap.",
            ha="center",
            va="center",
            fontsize=10,
        )
        table = ax.table(
            cellText=SS_CODE_ROWS,
            colLabels=["Code", "Description"],
            colWidths=[0.18, 0.82],
            loc="center",
            cellLoc="left",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 1.4)
        for cell in table.get_celld().values():
            cell.get_text().set_ha("left")
            cell.get_text().set_va("center")
            cell.get_text().set_wrap(True)
        image_path = self._save_plot(fig, "ss_letter_codes")
        plt.close(fig)
        return image_path

    def run(self) -> dict:
        """
        Compute SS assignments using DSSP.
        Saves the results to a text file, generates a discrete heatmap plot,
        writes an ss_README.md file explaining the SS letter codes,
        and saves a PNG of the SS code table for quick reference.

        Returns
        -------
        dict
            A dictionary containing the SS assignments with key "ss_data".
        """
        try:
            # Compute DSSP assignments for each frame.
            dssp = md.compute_dssp(self.traj)
            self.data = dssp  # shape: (n_frames, n_residues)
            self.results = {"ss_data": self.data}

            # Save the SS assignments to a file.
            data_path = self.outdir / "ss.dat"
            with open(data_path, "w") as f:
                for frame_idx, ss in enumerate(dssp):
                    f.write(f"Frame {frame_idx}: {', '.join(ss)}\n")

            # Generate the ss_README.md file.
            readme_path = self._generate_readme()
            codes_png_path = self._generate_reference_image()
            self.results.update({
                "ss_readme": readme_path,
                "ss_codes_plot": codes_png_path,
            })
            # Automatically generate the default heatmap plot.
            self.plot()
            return self.results
        except Exception as e:
            raise AnalysisError(f"SS analysis failed: {e}")

    def plot(self, data=None, **kwargs):
        """
        Generate a heatmap plot of SS assignments over frames.

        The heatmap displays a discrete colormap where each color corresponds to a specific SS code.
        The colorbar tick labels show the corresponding SS letters.
        The y-axis tick labels represent residue indices as whole numbers starting from 1.

        Parameters
        ----------
        data : array-like, optional
            SS data to plot. If None, uses self.data.
        kwargs : dict
            Customizable options:
              - title: Plot title (default "SS Heatmap").
              - xlabel: X-axis label (default "Frame").
              - ylabel: Y-axis label (default "Residue Index").
              - filename: Base filename for the plot (default "ss").
              - cmap: Optionally override the discrete colormap.

        Returns
        -------
        Path
            The file path to the saved plot.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No SS data available to plot. Please run analysis first.")

        # Map SS letters to numeric values.
        mapping = {"H": 1, "B": 2, "E": 3, "G": 4, "I": 5, "T": 6, "S": 7, "C": 0, " ": 0}
        numeric = np.array([[mapping.get(s, 0) for s in frame] for frame in data])
        
        # Define a discrete colormap with 8 distinct colors.
        from matplotlib.colors import ListedColormap, BoundaryNorm
        distinct_colors = ['#AAAAAA', '#FF0000', '#FFA500', '#0000FF', '#008000', '#FF00FF', '#FFFF00', '#00FFFF']
        cmap = kwargs.get("cmap", ListedColormap(distinct_colors))
        boundaries = np.arange(-0.5, 8, 1)
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
        # Define tick labels corresponding to the mapped values.
        tick_labels = ["C", "H", "B", "E", "G", "I", "T", "S"]
        cbar.set_ticklabels(tick_labels)
        cbar.set_label("SS Code")
        # Ensure y-axis tick labels are whole numbers starting from 1.
        n_residues = numeric.shape[1]
        plt.yticks(ticks=np.arange(n_residues), labels=np.arange(1, n_residues + 1))

        plt.tight_layout()

        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

