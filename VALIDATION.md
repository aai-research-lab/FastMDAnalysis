# FastMDAnalysis Validation

This document describes the validation script `validate_fastmda.py` which compares FastMDAnalysis calculations with reference implementations from MDTraj and MDAnalysis.

## Overview

The validation script performs the following analyses on the Trp-cage test dataset:

- **RMSD** (Root Mean Square Deviation)
- **RMSF** (Root Mean Square Fluctuation)
- **Radius of Gyration (Rg)**
- **Hydrogen Bonds**
- **Secondary Structure (DSSP)**
- **SASA** (Solvent Accessible Surface Area)
  - Total SASA
  - Per-residue SASA
  - Average per-residue SASA
- **Dimensionality Reduction** (PCA, MDS, t-SNE)
- **Clustering** (KMeans, DBSCAN, Hierarchical)

## Usage

### Basic Usage

```bash
python validate_fastmda.py
```

This runs validation with default settings:
- Frame selection: `0:-1:10` (every 10th frame, from start to end)
- Atom selection: `protein` (all protein atoms)
- Output directory: `validation_output`

### Custom Frame Selection

```bash
python validate_fastmda.py --frames 0:100:5
```

Frame selection format: `START:STOP:STRIDE`
- `START`: First frame (0-indexed)
- `STOP`: Last frame (exclusive, use `-1` for last frame)
- `STRIDE`: Step size

Examples:
- `0:100:1` - Frames 0-99, every frame
- `0:-1:10` - All frames, every 10th frame
- `10:50:2` - Frames 10-49, every 2nd frame

### Custom Atom Selection

```bash
python validate_fastmda.py --atoms "name CA"
```

Uses MDTraj selection syntax. Examples:
- `"protein"` - All protein atoms
- `"name CA"` - Only alpha carbons
- `"protein and name CA"` - Protein alpha carbons
- `"backbone"` - Backbone atoms

### Custom Output Directory

```bash
python validate_fastmda.py --output-dir my_validation
```

### Complete Example

```bash
python validate_fastmda.py --frames 0:200:5 --atoms "protein and name CA" --output-dir validation_CA
```

## Output Files

The validation script generates two output files:

### 1. validation_report.json

A detailed JSON report containing all validation results with full statistics and comparison metrics.

### 2. validation_summary.csv

A CSV table with the following columns:

| Column | Description |
|--------|-------------|
| `analysis_name` | Name of the analysis (e.g., "RMSD", "SASA (total)") |
| `backend` | Reference backend used (mdtraj/mdanalysis) |
| `metric` | Specific metric name (e.g., "rmsd", "total_sasa") |
| `status` | Validation status: pass/warn/fail/error/info |
| `shape_match` | Whether array shapes match (True/False) |
| `max_abs_diff` | Maximum absolute difference |
| `mean_abs_diff` | Mean absolute difference |
| `rmse` | Root Mean Square Error |
| `mismatch_count` | Number of elements with differences > 1e-6 |
| `detail` | Detailed description of the result |
| `fastmda_min` | Minimum value in FastMDAnalysis result |
| `fastmda_max` | Maximum value in FastMDAnalysis result |
| `fastmda_mean` | Mean value in FastMDAnalysis result |
| `fastmda_std` | Standard deviation in FastMDAnalysis result |
| `ref_min` | Minimum value in reference result |
| `ref_max` | Maximum value in reference result |
| `ref_mean` | Mean value in reference result |
| `ref_std` | Standard deviation in reference result |
| `fastmda_shape` | Shape of FastMDAnalysis array |
| `ref_shape` | Shape of reference array |

## Validation Criteria

### Pass Criteria

- **RMSE < 1e-4**: Excellent agreement
- **RMSE < 1e-2**: Good agreement
- **Shape match**: Arrays must have the same shape
- **Secondary Structure**: >95% match rate for string arrays

### Status Levels

- **pass**: Validation succeeded, results match within tolerance
- **warn**: Results differ more than expected but are reasonable
- **fail**: Significant differences detected
- **error**: Validation could not be completed due to errors
- **info**: Informational comparison (e.g., counts that can't be directly compared)

## Example Results

With default settings (500 frames, protein atoms), typical results are:

```
RMSD: Excellent agreement (RMSE=0.00e+00) ✓
RMSF: Excellent agreement (RMSE=0.00e+00) ✓
Radius of Gyration: Excellent agreement (RMSE~1e-09) ✓
Secondary Structure: 100% match with MDTraj DSSP ✓
SASA (total): Good agreement (RMSE~1e-04) ✓
SASA (per-residue): Excellent agreement (RMSE~3e-05) ✓
SASA (avg per-residue): Excellent agreement (RMSE~1e-06) ✓
Clustering (KMeans): Excellent agreement (RMSE=0.00e+00) ✓
Clustering (DBSCAN): Excellent agreement (RMSE=0.00e+00) ✓
Clustering (Hierarchical): Excellent agreement (RMSE=0.00e+00) ✓
```

## Clustering Validation Details

Clustering validation compares FastMDAnalysis results against direct sklearn calls:

- **KMeans**: Validates against sklearn.cluster.KMeans with the same parameters (random_state=42, n_init=10). Labels are compared after accounting for 0-based (sklearn) vs 1-based (FastMDAnalysis) indexing.

- **DBSCAN**: Validates against sklearn.cluster.DBSCAN with precomputed RMSD distance matrix. Raw labels (before FastMDAnalysis relabeling) are compared to ensure exact match.

- **Hierarchical**: Validates against scipy.cluster.hierarchy linkage with ward method and fcluster. Labels should match exactly since the algorithm is deterministic.

All three methods should achieve RMSE=0.00e+00, indicating perfect agreement with sklearn.

## Requirements

The validation script requires:

- FastMDAnalysis (with TrpCage dataset)
- MDTraj >= 1.9.7
- NumPy >= 1.21.0
- MDAnalysis (optional, for additional comparisons)

## Notes

- MDAnalysis comparisons are currently commented out in the RMSD validation for performance reasons
- Hydrogen bonds validation provides informational counts rather than direct comparison
- Dimensionality reduction results are validated for correct shape/format rather than exact values (due to randomness in algorithms)
- **Clustering results are validated against sklearn directly** - FastMDAnalysis and sklearn should produce identical labels with the same parameters and data
- Some analyses may require trajectory bonds to be defined; the script handles these cases gracefully

## Troubleshooting

### "No bonds found in topology"

Some analyses (like hydrogen bonds) require bond information. With certain atom selections (e.g., only CA atoms), bonds may not be present. This is expected behavior and reported as an error for that specific analysis.

### "perplexity must be less than n_samples"

t-SNE requires at least `perplexity + 1` samples. When testing with very few frames, t-SNE may fail. Use more frames or the validation will skip this method.

### Network timeouts during installation

If pip installation times out, use:
```bash
pip install --no-cache-dir --index-url https://pypi.python.org/simple/ mdtraj numpy ...
```

## Development

To modify the validation script:

1. Edit `validate_fastmda.py`
2. Test with a small frame subset first:
   ```bash
   python validate_fastmda.py --frames 0:20:5
   ```
3. Run full validation:
   ```bash
   python validate_fastmda.py
   ```

Each validation function follows the pattern:
```python
def validate_<analysis>(fastmda, traj) -> List[Dict[str, Any]]:
    # Run FastMDAnalysis
    fmda_result = fastmda.<analysis>()
    
    # Run reference implementation
    ref_result = mdtraj.<analysis>(traj)
    
    # Compare and return results
    comparison = compare_arrays(fmda_result, ref_result, name)
    return [comparison]
```
