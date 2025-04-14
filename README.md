# FastMDAnalysis

FastMDAnalysis is a Python package for fast and automated MD trajectory analysis using MDTraj. It provides a user-friendly Python API as well as a command-line interface (CLI) to perform standard analyses including:

- **rmsd**: Root-Mean-Square Deviation analysis
- **rmsf**: Root-Mean-Square Fluctuation analysis
- **rg**: Radius of Gyration analysis
- **hbonds**: Hydrogen Bonds analysis
- **cluster**: Cluster Analysis (supports methods such as DBSCAN and KMeans)
- **ss**: Secondary Structure analysis (using DSSP, with results displayed as an SS heatmap)
- **sasa**: Solvent Accessible Surface Area analysis (total SASA, per-residue SASA, and average per-residue SASA)
- **dimred**: Dimensionality Reduction analysis (using PCA, MDS, and t-SNE)

All analyses can be run via a simple Python API or through the command-line interface.

## Installation

To install FastMDAnalysis, navigate to the root directory of the project and run:

```bash
pip install .

```

## Development 

To install in editable (development) mode:

```bash
pip install -e .
