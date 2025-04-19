# FastMDAnalysis Package Structure

Below is the directory structure of the FastMDAnalysis package:

FastMDAnalysis/
├── FastMDAnalysis/
│   ├── __init__.py         # Package initialization. Uses load_trajectory from utils and exposes all analyses.
│   ├── cli.py              # Command-line interface supporting all subcommands (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, dimred)
│   ├── datasets.py         # New module providing paths for example datasets (e.g., ubiquitin, trp-cage)
│   ├── analysis/
│   │   ├── __init__.py     # Imports and exposes individual analysis modules.
│   │   ├── base.py         # Base analysis class and common exceptions.
│   │   ├── rmsd.py         # RMSD analysis implementation.
│   │   ├── rmsf.py         # RMSF analysis implementation.
│   │   ├── rg.py           # Radius of gyration analysis implementation.
│   │   ├── hbonds.py       # Hydrogen bonds analysis implementation.
│   │   ├── cluster.py      # Clustering analysis implementation (with DBSCAN, KMeans, hierarchical, dendrogram, etc.)
│   │   ├── ss.py           # Secondary structure (SS) analysis implementation.
│   │   ├── dimred.py       # Dimensionality reduction analysis implementation.
│   │   └── sasa.py         # Solvent accessible surface area (SASA) analysis implementation.
│   └── utils.py            # Utility functions (e.g., load_trajectory supports multiple inputs, create_dummy_trajectory, etc.)
├── tests/
│   └── tests.py            # Unit tests that use real datasets (e.g., ubiquitin) from the datasets module.
├── examples/
│   └── (example scripts showing API and CLI usage)
├── data/
│   ├── ubiquitin.dcd       # Example ubiquitin trajectory file.
│   ├── ubiquitin.pdb       # Example ubiquitin topology file.
│   ├── trp_cage.dcd        # Example trp-cage trajectory file.
│   └── trp_cage.pdb        # Example trp-cage topology file.
├── setup.py                # Package setup file.
├── requirements.txt        # List of Python dependencies.
└── README.md               # Package README with usage, API, CLI, testing, contributing, and license information.

