# FastMDAnalysis Package Structure

This document explains the organization of the FastMDAnalysis package. The [**`directory tree`**](./tree) summarizes the layout of the code, tests, data, and supporting files.


## Detailed Description

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

