# FastMDAnalysis Benchmark

This benchmark compares FastMDAnalysis, MDTraj, and MDAnalysis on **complete workflow performance** including computation and figure generation.

## Overview

The benchmark measures performance on the following:
- **Runtime**: Complete workflow time averaged over 5 iterations (including figure generation)
- **Peak Memory**: Maximum memory usage during execution
- **Standard Deviation**: Consistency of performance across multiple runs

## Dataset

- **TrpCage**: 500 frames selected using `frames=(0, -1, 10)` from 4999 total frames
- **Atom Selection**: Protein atoms

## Analyses Performed

All three libraries perform the same analyses:
- **RMSD** (Root Mean Square Deviation)
- **RMSF** (Root Mean Square Fluctuation)
- **RG** (Radius of Gyration)
- **Cluster** (KMeans, DBSCAN, Hierarchical)

## Usage

### Prerequisites

```bash
pip install mdtraj numpy scikit-learn scipy MDAnalysis matplotlib psutil
```

### Running the Benchmark

#### Full Workflow Benchmark (Recommended)
```bash
# Run the complete workflow benchmark including figure generation
python benchmark_full_workflow.py
```

This generates:
- `benchmark_full_workflow.png` - Runtime and memory comparison with error bars
- `benchmark_summary_table.png` - Detailed summary table with statistics

#### Legacy Benchmarks

For historical comparison, older benchmark scripts are also available:

```bash
# Pure computation only (no plotting)
python benchmark_performance.py

# Computation with visualizations
python benchmark_with_visualization.py
```

## Results

### Full Workflow (Averaged over 5 iterations)

| Library | Runtime (avg) | Std Dev | Memory (avg) | Std Dev |
|---------|---------------|---------|--------------|---------|
| **FastMDAnalysis** | ~6.4s | ±0.04s | ~101 MB | ±13 MB |
| **MDTraj** | ~2.1s | ±0.12s | ~15 MB | ±9 MB |
| **MDAnalysis** | ~2.3s | ±0.15s | ~14 MB | ±18 MB |

**FastMDA/MDTraj runtime ratio: ~3.0x** - FastMDA is ~3x slower but includes automatic figure generation and workflow organization.

## Key Insights

### 1. Complete Workflow Performance
- FastMDA provides a complete workflow with automatic figure generation
- Total runtime includes computation + figure generation + file I/O
- FastMDA is ~3x slower than raw MDTraj, but provides automated visualization

### 2. Memory Usage
- FastMDA uses ~7x more memory than MDTraj/MDAnalysis
- This is due to FastMDA's internal caching and convenience features
- Memory usage is still reasonable for typical MD analysis tasks

### 3. Trade-offs
- **FastMDA**: Best for interactive analysis, teaching, and rapid prototyping
  - Single-line API with automatic figures
  - ~3x slower but significantly more convenient
- **MDTraj**: Best for production pipelines requiring maximum speed
  - Minimal memory footprint
  - Manual figure generation required
- **MDAnalysis**: Best for complex trajectory manipulations
  - Similar performance to MDTraj
  - More flexible for custom analyses

## Important Notes

### What This Benchmark Measures
- **Complete workflow**: Computation + figure generation + file I/O
- **Realistic usage**: Reflects typical user workflows
- **Apples-to-apples**: All libraries generate identical outputs (4 figures each)
- **Multiple iterations**: Averaged over 5 runs with standard deviation

### Figures Generated (Per Library)
All three libraries generate exactly 4 figures for fair comparison:
1. **RMSD**: Line plot of RMSD vs Frame
2. **RMSF**: Bar plot of RMSF per Atom
3. **RG**: Line plot of Radius of Gyration vs Frame
4. **Cluster**: Combined plot showing KMeans, DBSCAN, and Hierarchical clustering results

## Interpretation

The benchmark demonstrates that:

1. **FastMDA's Design Philosophy**
   - Optimized for ease of use and convenience
   - Automatic figure generation and workflow organization
   - ~3x slowdown is acceptable for the convenience gained
   - ~6.4s total for complete workflow vs ~2.1s for manual approach

2. **Performance Breakdown**
   - Computation: FastMDA uses MDTraj backend (minimal overhead)
   - Figure generation: Adds ~2-3 seconds
   - File I/O and organization: Adds ~1-2 seconds
   - Total overhead: ~4 seconds for convenience features

3. **When to Use Each Library**
   - **FastMDAnalysis**: Interactive analysis, teaching, exploratory research
     - Single-line API with automatic visualization
     - Worth the 3x slowdown for convenience
   - **MDTraj**: Production pipelines, batch processing, performance-critical tasks
     - Minimal overhead, manual control
     - Best for experienced users
   - **MDAnalysis**: Complex trajectory manipulations, custom analyses
     - Most flexible, similar performance to MDTraj
     - Steeper learning curve

## Example Usage

All three libraries perform the same analyses with figure generation:

### FastMDAnalysis (Simplest)
```python
from fastmdanalysis import FastMDAnalysis
fastmda = FastMDAnalysis(traj_file, top_file, frames=(0, -1, 10), atoms="protein")
fastmda.rmsd(ref=0)  # Auto-generates rmsd.png
fastmda.rmsf()       # Auto-generates rmsf.png
fastmda.rg()         # Auto-generates rg.png
fastmda.cluster(methods=['kmeans', 'dbscan', 'hierarchical'])  # Auto-generates cluster plots
```

### MDTraj (Manual Control)
```python
import mdtraj as md
import numpy as np
from sklearn.cluster import KMeans, DBSCAN
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

# Load and prepare
traj = md.load('traj.dcd', top='top.pdb')
traj = traj[0::10]
atom_indices = traj.topology.select('protein')
traj = traj.atom_slice(atom_indices)

# Analyses
rmsd = md.rmsd(traj, traj, frame=0)
avg_xyz = np.mean(traj.xyz, axis=0, keepdims=True)
ref = md.Trajectory(avg_xyz, traj.topology)
rmsf = md.rmsf(traj, ref)
rg = md.compute_rg(traj)

# Clustering
rmsd_matrix = np.empty((traj.n_frames, traj.n_frames))
for i in range(traj.n_frames):
    rmsd_matrix[i] = md.rmsd(traj, traj, frame=i)
    
kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
kmeans_labels = kmeans.fit_predict(rmsd_matrix)

dbscan = DBSCAN(eps=0.5, min_samples=2, metric='precomputed')
dbscan_labels = dbscan.fit_predict(rmsd_matrix)

rmsd_matrix_sym = (rmsd_matrix + rmsd_matrix.T) / 2
np.fill_diagonal(rmsd_matrix_sym, 0)
condensed_dist = squareform(rmsd_matrix_sym)
linkage_matrix = linkage(condensed_dist, method='ward')
```

## Extending the Benchmark

To benchmark additional analyses or datasets:

1. Add new analysis functions to each benchmark function
2. Update the dataset by changing TrpCage to another dataset
3. Modify frame selection with different stride values

## License

This benchmark script is part of the FastMDAnalysis project and follows the same MIT license.
