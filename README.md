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

#### Examples:


- **RMSD Analysis:**

```python
from FastMDAnalysis import FastMDAnalysis

fastmda = FastMDAnalysis("traj.dcd", "top.pdb")

# Run RMSD analysis (uses the default frames and atom selection unless overridden):
rmsd_analysis = fastmda.rmsd(ref=0)

# Optionally retrieve rmsd data 
print("RMSD Data:", rmsd_analysis.data)

# Optionally customize rmsd plot
rmsd_analysis.plot()

```

### Command-Line Interface (CLI)
After installation, you can run FastMDAnalysis from the command line using the fastmda command. Global options allow you to specify the trajectory, topology, frame selection, and atom selection.

#### Examples:

- **RMSD Analysis:**

```bash
fastmda rmsd -traj traj.dcd -top top.pdb 
```

## License

FastMDAnalysis is licensed under the MIT License. See the LICENSE file for more details.


## Acknowledgements

FastMDAnalysis leverages MDTraj for trajectory analysis. It also relies on popular Python libraries such as NumPy, scikit-learn, and Matplotlib for data processing and visualization. Special thanks to the community for their continuous support and contributions.
