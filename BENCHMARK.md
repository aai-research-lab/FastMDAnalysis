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
- `benchmark_loc_slide.png` - LOC comparison slide (effective non-blank, non-comment lines)
- `benchmark_cli_commands.png` - FastMDA CLI command count (analysis vs extra plotting)

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
| **FastMDAnalysis** | 2.47 s | ± 0.81 s | 82.6 MB | ± 3.9 MB |
| **MDTraj** | 2.33 s | ± 0.01 s | 26.1 MB | ± 8.6 MB |
| **MDAnalysis** | 7.65 s | ± 0.03 s | 19.3 MB | ± 2.4 MB |

**FastMDA/MDTraj runtime ratio: 1.06×** – FastMDA now runs within ~6 % of MDTraj while still handling figure generation. Memory remains ~3.2× higher, primarily due to distance-matrix workloads retained for clustering.

### Lines of Code Benchmark (Workflow Implementations)

To quantify ergonomics, we also measured the effective lines of code (LOC) needed to express each workflow within `benchmark_full_workflow.py`. LOC counts exclude blank lines, comments, and docstrings, relying on Python's token stream for accuracy.

| Library | Effective LOC |
|---------|---------------|
| **FastMDAnalysis** | 2 |
| **MDTraj** | 72 |
| **MDAnalysis** | 94 |

**FastMDA/MDTraj LOC ratio: 0.03×** – the FastMDA workflow needs just two executable lines to reproduce the full workflow, versus dozens for the MDTraj and MDAnalysis scripts. A dedicated visualization (`benchmark_loc_slide.png`) is generated alongside the runtime/memory plots.

#### FastMDAnalysis CLI Equivalent

The same workflow can be triggered from the bundled CLI with a single command; the benchmark script simply expands it into Python so we can capture runtime, memory, and LOC alongside the other libraries.

```bash
fastmda analyze -traj src/fastmdanalysis/data/trp_cage.dcd -top src/fastmdanalysis/data/trp_cage.pdb --frames 0,-1,10 --atoms protein --include rmsd rmsf rg cluster -o fastmda_cli_output
```

`benchmark_cli_commands.png` compares the number of CLI invocations each workflow needs to reproduce the full set of figures. FastMDA finishes in **1 command** because plots are emitted automatically. MDTraj and MDAnalysis lack an equivalent workflow CLI, so reaching parity requires four manual steps (one per analysis/plot), which is what we visualize for consistency.

#### How the 2 LOC Were Measured

`benchmark_full_workflow.py` feeds `_count_effective_loc` with the literal two-line snippet a user types to reproduce the workflow:

```python
fastmda = FastMDAnalysis(traj, top, frames=(0, -1, 10), atoms="protein", keep_full_traj=False)
fastmda.analyze(include=["rmsd", "rmsf", "rg", "cluster"], verbose=False, output="fastmda_analyze_output", options={...})
```

The helper removes blanks, comments, and docstrings before counting tokens, so FastMDAnalysis reports exactly **2** effective lines. MDTraj and MDAnalysis use the full benchmark implementations (manual loading, plotting, clustering), which is why their counts stay at 72 and 94 lines respectively.
The elided `options={...}` simply disables data dumps and matches the clustering parameters used throughout the benchmark.

## Key Insights

### 1. Complete Workflow Performance
- FastMDA provides a complete workflow with automatic figure generation
- Total runtime includes computation + figure generation + file I/O
- After aligning workloads (shared RMSD distance matrices, identical clustering), FastMDA now runs within ~6 % of MDTraj while still producing figures automatically

### 2. Memory Usage
- FastMDA uses ~3.2× more memory than MDTraj and ~4.3× more than MDAnalysis in this configuration
- The extra footprint stems from retaining the pairwise RMSD matrix for clustering (also constructed in MDTraj/MDAnalysis for parity) and intermediate analysis buffers
- Memory usage remains manageable for typical MD analysis tasks, but further reduction is an active optimization target

### 3. Trade-offs
- **FastMDA**: Best for interactive analysis, teaching, and rapid prototyping
  - Single-line API with automatic figures
  - Runtime now near MDTraj; still trades higher memory for convenience
- **MDTraj**: Best for production pipelines requiring maximum speed
  - Minimal memory footprint
  - Manual figure generation required
- **MDAnalysis**: Best for complex trajectory manipulations
  - When forced to compute the same RMSD distance matrix, runtime increases to ~7.5 s but feature parity is maintained

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
   - With distance-matrix reuse and disabled result caching, runtime overhead has effectively vanished relative to MDTraj while keeping the streamlined API

2. **Performance Breakdown**
   - Computation: FastMDA continues to rely on MDTraj primitives under the hood
   - Clustering: All libraries now build the same RMSD distance matrix, ensuring an apples-to-apples comparison
   - Figure generation/file I/O: Adds modest overhead but no longer dominates runtime

3. **When to Use Each Library**
   - **FastMDAnalysis**: Interactive analysis, teaching, exploratory research
     - One-stop workflow with figures; near-MDTraj runtime, higher memory
   - **MDTraj**: Production pipelines, batch processing, performance-critical tasks
     - Lowest memory footprint, manual control of outputs
   - **MDAnalysis**: Complex trajectory manipulations, custom analyses
     - Highly flexible; when matched to MDTraj’s workload, incurs higher runtime but identical outputs

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
