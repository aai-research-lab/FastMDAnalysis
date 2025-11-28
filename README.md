# FastMDAnalysis

FastMDAnalysis is a high‑level Python toolkit for **reproducible, end‑to‑end analysis of molecular dynamics (MD) trajectories**.  
It provides a unified API and command‑line interface that wrap and standardize core functionality from
MDTraj, MDAnalysis, scikit‑learn, and SciPy into a single, workflow‑oriented package.

FastMDAnalysis is designed for:

- **Rapid, scriptable analysis** of MD trajectories (RMSD, RMSF, Rg, SASA, H‑bonds, secondary structure, clustering, etc.).  
- **Standardized data structures** and figure‑ready outputs.  
- **Reproducible pipelines** that can be run from the command line or within Python.  

> **Reference dataset:** Many examples and validation tests use the TrpCage mini‑protein dataset bundled with the package.

---

## Installation

FastMDAnalysis is distributed on PyPI.

```bash
pip install fastmdanalysis
```

We recommend using a fresh virtual environment:

```bash
python -m venv fastmda_env
source fastmda_env/bin/activate      # Windows: fastmda_env\Scripts\activate

pip install --upgrade pip
pip install fastmdanalysis
```

Additional optional dependencies (for some analysis modules and validation) include:

- `mdtraj`
- `MDAnalysis`
- `scikit-learn`
- `scipy`
- `matplotlib`
- `seaborn` (for some plotting recipes)

Install them as needed:

```bash
pip install mdtraj MDAnalysis scikit-learn scipy matplotlib seaborn
```

---

## Quick Start (Python API)

The central entry point is the `FastMDAnalysis` class, which loads your trajectory and provides
convenience methods for common analyses.

```python
from fastmdanalysis import FastMDAnalysis

fmda = FastMDAnalysis(
    traj="traj.dcd",
    top="topology.pdb",
    frames=(0, -1, 10),      # (start, stop, stride); -1 means "last frame"
    atoms="protein"          # MDTraj/MDAnalysis-style atom selection
)

# Core structural analyses
rmsd = fmda.rmsd(ref=0)
rmsf = fmda.rmsf()
rg   = fmda.rg()

# Hydrogen bonds, secondary structure, SASA
hbonds = fmda.hbonds()
ss     = fmda.ss()
sasa   = fmda.sasa(probe_radius=0.14)

# Dimensionality reduction (PCA / MDS / t-SNE) and clustering
dimred   = fmda.dimred(methods=["pca", "mds", "tsne"])
clusters = fmda.cluster(methods=["kmeans", "dbscan", "hierarchical"])
```

Most analysis objects expose their results via a `.data` attribute (arrays or dictionaries)
and provide helper methods for plotting and exporting.

---

## Command‑Line Interface (CLI)

FastMDAnalysis also provides a CLI for running common pipelines without writing Python code.

After installation, the `fastmda` command should be available in your environment. Typical usage
follows the pattern:

```bash
fastmda run   --traj traj.dcd   --top top.pdb   --frames 0:-1:10   --atoms "protein"   --config analysis.yaml
```

Where `analysis.yaml` describes which modules to run (e.g. RMSD, RMSF, Rg, SASA, clustering)
and how to configure them. Example configuration files and workflows are provided in the
repository under `examples/` (or `docs/examples/`, depending on layout).

See:

```bash
fastmda --help
fastmda run --help
```

for up‑to‑date CLI options.

---

## Built‑in Datasets

FastMDAnalysis ships with small reference datasets for testing and examples. For instance,
the **TrpCage** mini‑protein dataset is exposed via:

```python
from fastmdanalysis.datasets import TrpCage

print(TrpCage.traj)  # path to the bundled trajectory (e.g., .dcd)
print(TrpCage.top)   # path to the bundled topology   (e.g., .pdb)
```

These paths can be passed directly into `FastMDAnalysis` or used with MDTraj/MDAnalysis.

---

## Validation Against the PyPI Release

This repository includes a standalone validation script, `validate_fastmda.py`, that
compares the **published** `fastmdanalysis` package on PyPI against reference
implementations (MDTraj, MDAnalysis, scikit‑learn, SciPy) on the TrpCage dataset.

The goal is to answer:

> “Does the version of `fastmdanalysis` installed from PyPI numerically agree with
> standard libraries on a well‑defined benchmark?”

### 1. Create a fresh environment and install from PyPI

From the repository root:

```bash
python -m venv .venv-pypi-validation
source .venv-pypi-validation/bin/activate  # Windows: .venv-pypi-validation\Scripts\activate

pip install --upgrade pip
pip install fastmdanalysis
pip install mdtraj MDAnalysis scikit-learn scipy
```

> **Important:** Do **not** install the local repo in this environment
> (no `pip install -e .`). This ensures that `import fastmdanalysis` uses the
> **PyPI** release, not the local source tree.

### 2. Run the validation script

From the repository root (adjust the path if `validate_fastmda.py` lives in `scripts/`):

```bash
python validate_fastmda.py
```

The default settings are:

- Frames: `0:-1:10` (start:stop:stride)
- Atoms: `protein`
- Output directory: `validation_output/`

You can override these, for example:

```bash
python validate_fastmda.py   --frames 0:-1:5   --atoms "protein"   --output-dir validation_output_fast_stride5
```

### 3. Outputs

The script produces:

- `validation_output/validation_report.json`  
  Detailed per‑metric results (RMSD, RMSF, Rg, H‑bonds, secondary structure,
  SASA, dimensionality reduction, clustering), including shapes, RMSE, and
  summary statistics.

- `validation_output/validation_summary.csv`  
  CSV summary (one row per metric) with:
  - `analysis_name`, `backend`, `metric`, `status` (`pass/warn/fail/error/info`)
  - `max_abs_diff`, `mean_abs_diff`, `rmse`, `mismatch_count`
  - `fastmda_*` and `ref_*` statistics (min, max, mean, std)
  - `fastmda_shape`, `ref_shape`

- `validation_results.csv` in the repository root  
  A mirrored copy of the CSV for quick inspection and CI hooks.

### 4. Interpreting the results

Each row in the CSV corresponds to one validation check, e.g.:

- `RMSD` vs MDTraj (`metric = rmsd`)
- `RMSF` vs MDTraj (`metric = rmsf`)
- `Radius of Gyration` vs MDTraj (`metric = rg`)
- `Hydrogen Bonds` vs MDTraj (`metric = hbond_per_frame`)
- `Secondary Structure` vs MDTraj (`metric = dssp`)
- `SASA` total / per‑residue / average per‑residue vs MDTraj
- Dimensionality reduction (PCA, MDS, t‑SNE) vs direct scikit‑learn
- Clustering (k‑means, DBSCAN, hierarchical) vs direct scikit‑learn/SciPy

A typical successful run (for a consistent PyPI release) should show:

- All core structural metrics (`RMSD`, `RMSF`, `Rg`, `SASA`) with `status = pass`
  and RMSE close to zero.
- Secondary structure and H‑bond metrics with very high agreement.
- Dimensionality reduction and clustering consistent with direct sklearn/SciPy
  calls under fixed random seeds.

Use this script whenever you:

- Publish a new version to PyPI.
- Change numerical kernels or analysis defaults.
- Want to confirm that the public release matches the validated behavior in this repo.

---

## Contributing

Contributions are welcome. Typical contribution workflow:

1. Fork the repository and create a feature branch:
   ```bash
   git checkout -b feature/my-improvement
   ```
2. Implement your changes with tests where appropriate.
3. Ensure the test suite passes and that `validate_fastmda.py` still reports
   consistent behavior (for public releases).
4. Open a pull request with a clear description and, if relevant, a short
   example of the new behavior.

Please follow the existing code style and docstring conventions where possible.

---

## Citing FastMDAnalysis

If you use FastMDAnalysis in published work, please cite the associated preprint:

> Aina, A. O., *et al.* **FastMDAnalysis: Reproducible, End‑to‑End Molecular Dynamics Analysis in Python.**  
> ChemRxiv (2025). doi:10.26434/chemrxiv-2025-x8xnq

(Replace with the final journal citation once available.)

---

## License

FastMDAnalysis is distributed under an open‑source license. See the `LICENSE` file in this
repository for details.

