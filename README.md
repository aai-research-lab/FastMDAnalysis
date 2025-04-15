# FastMDAnalysis

FastMDAnalysis is a software for fast and automated analysis of molecular dynamics (MD) trajectories. It provides a unified Python API as well as a command‐line interface (CLI) for performing a variety of MD analyses. This package is designed to simplify your workflow by allowing you to load a trajectory once (with options for frame and atom selection) and then run multiple analyses without repeating input file details.

## Features

- **rmsd**: Calculate Root-Mean-Square Deviation relative to a reference frame.
- **rmsf**: Compute per-atom Root-Mean-Square Fluctuation.
- **rg**: Determine the Radius of Gyration for each frame.
- **hbonds**: Detect and count hydrogen bonds using the Baker-Hubbard algorithm.
- **cluster**: Perform clustering (DBSCAN and/or KMeans) on trajectory frames.
- **ss**: Compute Secondary Structure (SS) assignments using DSSP with a discrete heatmap.
- **sasa**: Compute Solvent Accessible Surface Area (SASA) in multiple ways:
  - Total SASA vs. frame.
  - Per-residue SASA vs. frame (heatmap).
  - Average per-residue SASA (bar plot).
- **dimred**: Perform dimensionality reduction (PCA, MDS, t-SNE) to project high-dimensional data into 2D.

## Installation

Navigate to the root directory of the package (the directory containing `setup.py`).

For a standard installation, run:
```bash
pip install .
```

For development (editable) mode, run:
```bash
pip install -e .
```

## Usage

### Python API

Instantiate a FastMDAnalysis object with your trajectory and topology file paths. Optionally, specify frame selection and atom selection. Frame selection is provided as a tuple (start, stop, stride). Negative indices (e.g., -1 for the last frame) are supported. If no options are provided, the entire trajectory and all atoms are used by default.

#### Example:

```python
from FastMDAnalysis import FastMDAnalysis

# Load the trajectory from file, selecting frames 0 to the last frame with a stride of 10,
# and use only "protein" atoms as the default selection.
fastmda = FastMDAnalysis("path/to/trajectory.dcd", "path/to/topology.pdb", frames=(0, -1, 10), atoms="protein")

# Run RMSD analysis (uses the default frames and atom selection unless overridden):
rmsd_analysis = fastmda.rmsd(ref=0)
print("RMSD Data:", rmsd_analysis.data)
rmsd_plot_file = rmsd_analysis.plot()
print("RMSD plot saved to:", rmsd_plot_file)

# Run Secondary Structure (SS) analysis:
ss_analysis = fastmda.ss()
print("SS Data:", ss_analysis.data)
ss_plot = ss_analysis.plot()
print("SS plot saved to:", ss_plot)

# Run SASA analysis with a probe radius of 0.14 nm:
sasa_analysis = fastmda.sasa(probe_radius=0.14)
print("SASA Data:", sasa_analysis.data)
sasa_plots = sasa_analysis.plot(option="all")
print("SASA plots saved to:", sasa_plots)

# Run Dimensionality Reduction using PCA and t-SNE:
dimred_analysis = fastmda.dimred(methods=["pca", "tsne"], atom_selection="protein and name CA")
print("DimRed Data:", dimred_analysis.data)
dimred_plots = dimred_analysis.plot()
print("DimRed plots saved to:", dimred_plots)
```

### Command-Line Interface (CLI)

After installation, you can run FastMDAnalysis from the command line using the fastmda command. Global options allow you to specify the trajectory, topology, frame selection, and atom selection.

#### Examples:

- **RMSD Analysis:**

```bash
fastmda rmsd -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rmsd_output --frames 0,-1,10 --atoms "protein" --ref 0
```

## License

FastMDAnalysis is licensed under the MIT License. See the LICENSE file for more details.


## Acknowledgements

FastMDAnalysis leverages MDTraj for trajectory analysis. It also relies on popular Python libraries such as NumPy, scikit-learn, and Matplotlib for data processing and visualization. Special thanks to the community for their continuous support and contributions.
