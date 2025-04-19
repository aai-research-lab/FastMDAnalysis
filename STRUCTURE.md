# FastMDAnalysis Package Structure

This document explains the organization of the FastMDAnalysis package. The directory tree below summarizes the layout of the code, tests, data, and supporting files.

FastMDAnalysis/ ├── FastMDAnalysis/ │ ├── init.py # Package initialization. Uses load_trajectory from utils and exposes all analyses. │ ├── cli.py # Command-line interface supporting all subcommands (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, dimred) │ ├── datasets.py # Module providing paths and metadata for example datasets (e.g., ubiquitin, trp-cage) │ ├── analysis/ │ │ ├── init.py # Imports and exposes individual analysis modules. │ │ ├── base.py # Base analysis class and common exceptions. │ │ ├── rmsd.py # RMSD analysis implementation. │ │ ├── rmsf.py # RMSF analysis implementation. │ │ ├── rg.py # Radius of gyration analysis implementation. │ │ ├── hbonds.py # Hydrogen bonds analysis implementation. │ │ ├── cluster.py # Clustering analysis implementation (DBSCAN, KMeans, hierarchical, dendrogram, etc.) │ │ ├── ss.py # Secondary structure (SS) analysis implementation. │ │ ├── dimred.py # Dimensionality reduction analysis implementation. │ │ └── sasa.py # Solvent accessible surface area (SASA) analysis implementation. │ └── utils.py # Utility functions (e.g., load_trajectory supports multiple inputs, create_dummy_trajectory, etc.) ├── tests/ │ └── tests.py # Unit tests that use real datasets (e.g., ubiquitin) from the datasets module. ├── examples/ │ └── (example scripts showing API and CLI usage) ├── data/ │ ├── ubiquitin.dcd # Example ubiquitin trajectory file. │ ├── ubiquitin.pdb # Example ubiquitin topology file. │ ├── trp_cage.dcd # Example trp-cage trajectory file. │ └── trp_cage.pdb # Example trp-cage topology file. ├── setup.py # Package setup file. ├── requirements.txt # List of Python dependencies. └── README.md # Package README with usage, API, CLI, testing, contributing, and license information.


## Detailed Description


```css
- **FastMDAnalysis/**  
  This is the core package directory:
  - **`__init__.py`** initializes the package, loads the trajectory via the extended `load_trajectory` function (which now supports multiple input types), and exposes the analysis classes (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, dimred).
  - **`cli.py`** implements the command-line interface, supporting subcommands to run different analyses while managing logging.
  - **`datasets.py`** provides a central location for dataset file paths and associated simulation metadata (such as time step, force field, integrator, temperature, pressure, and run-time tags).
  - **`analysis/`** contains modules for each type of analysis:
    - **`base.py`** holds the BaseAnalysis class and common exceptions.
    - **`rmsd.py`**, **`rmsf.py`**, **`rg.py`**, **`hbonds.py`**, **`cluster.py`**, **`ss.py`**, **`dimred.py`**, and **`sasa.py`** implement various analysis methods.
  - **`utils.py`** has helper functions such as `load_trajectory` (now extended to accept lists, comma-separated strings, or glob patterns) and (if needed) a function to create dummy trajectories.

- **tests/**  
  Contains unit tests (in `tests.py`) that verify the functionality of the package using real datasets (e.g., the ubiquitin dataset loaded from `datasets.py`).

- **examples/**  
  Provides example scripts that illustrate how to use the API and command-line interface.

- **data/**  
  Contains the actual dataset files used for testing and examples (such as ubiquitin and trp-cage trajectories and topologies). Ensure that these files are correctly placed so that the paths in `datasets.py` are valid.

- **setup.py & requirements.txt**  
  Standard files to allow installation and dependency management.

- **README.md**  
  The main documentation file with usage instructions, API descriptions, and contribution guidelines.

## How to Use This Structure

- **For Developers:**  
  This structure cleanly separates the core package code, tests, example usage, and data. It makes it easy to run unit tests and to extend functionality.

- **For Users:**  
  The README and datasets modules provide enough documentation and example data to help users get started with running analyses using either the CLI or the API.

You can save this file as `STRUCTURE.md` in the root directory of your repository and reference it in your README or documentation.

Enjoy building with FastMDAnalysis!
```
