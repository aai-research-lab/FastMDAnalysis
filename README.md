# FastMDAnalysis Validation Suite

This validation suite quantifies the numerical agreement between **FastMDAnalysis**
(as published on PyPI) and established reference libraries (MDTraj, scikit-learn,
SciPy) on the Trp-cage miniprotein benchmark used in the JCC submission.

The procedure below reproduces the **Accuracy Validation** results reported in the
manuscript and supporting information.

---

## Reproducing JCC Validation Results

### Step 1: Get the Validation Code

Create and activate a clean virtual environment, then check out the
`validation` branch of the lab repository.

Using `conda` (Recommended):

```bash
conda create -n fastmda_validation_env python=3.9
conda activate fastmda_validation_env
```

Clone the lab repository and switch to the validation branch:

```bash
git clone -b validation https://github.com/aai-research-lab/FastMDAnalysis.git
cd FastMDAnalysis
```

If you already have the repository cloned, you can instead run:

```bash
git fetch origin
git checkout validation
```

### Step 2: Install Dependencies

Install the **published** FastMDAnalysis package from PyPI and the minimal
dependencies required by the validation script.

```bash
pip install --upgrade pip
pip install fastmdanalysis mdtraj numpy scipy scikit-learn
```

Notes:

- Do **not** run `pip install -e .` in this environment. This ensures that
  `import fastmdanalysis` resolves to the **PyPI** installation, not the
  local source tree.
- MDAnalysis is **not** required for this validation workflow.

### Step 3: Run the Validation Script

From the `FastMDAnalysis` repository root (on the `validation` branch):

```bash
python validate_fastmda.py
```

By default, this reproduces the JCC validation configuration:

- Dataset: Trp-cage miniprotein (PDB ID: 1L2Y)
- Trajectory: bundled Trp-cage MD trajectory from the PyPI package
- Frames: `0:-1:10` (start at frame 0, go to last frame, stride = 10)  
  → 500 frames sampled from a 5000-frame trajectory
- Atom selection: `protein` (304 protein atoms, 20 residues)
- Output directory: `validation_output/`

You can optionally change the frame range, atom selection, or output directory,
for example:

```bash
python validate_fastmda.py   --frames 0:-1:5   --atoms "protein"   --output-dir validation_output_stride5
```

---

## Validation Output Files

After a successful run (using the default settings), the following files are
generated:

- `validation_output/validation_report.json`  
  Detailed, machine-readable report of all validation comparisons.

- `validation_output/validation_summary.csv`  
  Human-readable summary (one row per comparison) with key statistics.

These files correspond to the numerical validation results summarized in the
JCC manuscript (Accuracy Validation section and Table~1). 

### JSON Report (`validation_report.json`)

Each entry in the JSON report includes (fields may vary slightly by metric):

- `name` – analysis module (e.g., `RMSD`, `SASA (total)`, `Clustering (kmeans)`)
- `backend` – reference implementation (`mdtraj`, `sklearn`, `scipy`)
- `metric` – metric identifier (e.g., `rmsd`, `total_sasa`, `dimred_pca`)
- `status` – qualitative outcome (`pass`, `warn`, `fail`, `error`, `info`)
- `shape_match` – whether the result shapes are identical
- `max_abs_diff`, `mean_abs_diff`, `rmse`, `mismatch_count`
- `fastmda_stats`, `ref_stats` – min, max, mean, and std for each array
- `fastmda_shape`, `ref_shape`
- `detail` – explanatory message (e.g., “Excellent agreement (RMSE=0.00e+00)”).

### CSV Summary (`validation_summary.csv`)

The CSV file provides a compact view suitable for quick inspection,
spreadsheets, and automated regression checks. It includes columns:

- `analysis_name`, `backend`, `metric`, `status`, `shape_match`
- `max_abs_diff`, `mean_abs_diff`, `rmse`, `mismatch_count`
- `fastmda_min`, `fastmda_max`, `fastmda_mean`, `fastmda_std`
- `ref_min`, `ref_max`, `ref_mean`, `ref_std`
- `fastmda_shape`, `ref_shape`

Each row corresponds to a specific analysis/metric combination, e.g., RMSD vs
MDTraj or PCA vs scikit-learn.

---

## Validation Details

### Analyses Performed

The validation script tests the main analysis modules provided by
FastMDAnalysis against direct calls to the underlying reference libraries:

**Structural metrics (MDTraj backends)**

- RMSD – time series of backbone RMSD relative to a reference frame.
- RMSF – per-atom fluctuation relative to the average structure.
- Radius of gyration (Rg) – time series of mass-weighted Rg.
- SASA – total SASA per frame, per-residue SASA per frame, and average
  per-residue SASA over the trajectory.
- Hydrogen bonds – per-frame hydrogen-bond counts using the Baker–Hubbard
  geometric criteria.
- Secondary structure – DSSP assignments (simplified alphabet) per residue
  and frame.

**Statistical learning (scikit-learn and SciPy backends)**

- Dimensionality reduction:
  - PCA
  - Multidimensional scaling (MDS)
  - t-distributed stochastic neighbor embedding (t-SNE)
- Clustering:
  - K-means
  - DBSCAN (density-based clustering)
  - Hierarchical clustering (SciPy `linkage` with Ward method and `fcluster`).

For each module, FastMDAnalysis is configured to use the same trajectory,
selection, and hyperparameters as the corresponding direct MDTraj / scikit-learn /
SciPy call.

### Comparison Metrics and Criteria

For each analysis, FastMDAnalysis results are compared directly to the reference
results using:

- Root mean square error (RMSE)
- Maximum absolute difference
- Mean absolute difference
- Number of elements exceeding a fixed tolerance
- Shape consistency checks

The script assigns qualitative labels based on RMSE and shape agreement:

- **Excellent agreement** – typically RMSE < 1×10⁻⁴
- **Good agreement** – RMSE < 1×10⁻²
- `warn`, `fail`, or `error` – used when discrepancies are larger, shapes
  differ, or errors occur during computation.

For the Trp-cage benchmark used in the JCC submission, all core modules
(RMSD, RMSF, Rg, SASA, hydrogen bonds, secondary structure, dimensionality
reduction, and clustering) achieve **excellent agreement**, with differences
at or below numerical precision for most metrics.

### Reproducibility Notes

- **System:** Trp-cage miniprotein (PDB ID: 1L2Y)  
  100 ns trajectory, 5000 frames; validation uses 500 frames (stride = 10)
  and 304 protein atoms across 20 residues (“protein” selection).

- **Reference libraries:**  
  The JCC submission used versions consistent with:

  - `mdtraj` (e.g., 1.11.0)
  - `scikit-learn` (e.g., 1.7.2)
  - `scipy==1.13.1`

To closely match the published results, pin these versions (or the ones
reported in the manuscript) when installing dependencies.
