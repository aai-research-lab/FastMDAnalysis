# MD Analysis Performance Benchmark Comparison

This benchmark compares FastMDAnalysis, MDTraj, and MDAnalysis performance on molecular dynamics trajectory analysis tasks.

## Overview

The benchmark measures performance on the following metrics:
- **Runtime**: Total execution time including computation and plotting
- **Memory**: Peak memory usage during execution  
- **Lines of Code (LOC)**: Code complexity for equivalent functionality

## Dataset

- **TrpCage**: 500 frames selected using `--frames 0,-1,10` from 4999 total frames
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
pip install mdtraj numpy matplotlib scikit-learn scipy MDAnalysis python-pptx Pillow cairosvg PyYAML
```

### Running the Benchmark

```bash
# Set PYTHONPATH to include FastMDAnalysis src directory (if not installed)
export PYTHONPATH=/path/to/FastMDAnalysis/src:$PYTHONPATH

# Run the benchmark
python benchmark_performance.py
```

## Approaches Compared

### 1. FastMDAnalysis (CLI)
Single command (1 LOC):
```bash
fastmda analyze -traj <traj.dcd> -top <top.pdb> --frames 0,-1,10 --include cluster rmsd rg rmsf
```

### 2. MDTraj (Manual)
Manual analysis + plotting (~50 LOC):
- Load trajectory with frame selection
- Compute each analysis separately
- Create plots for each analysis
- Organize output files

### 3. MDAnalysis (Manual)
Manual analysis + plotting (~60 LOC):
- Load trajectory with Universe
- Iterate through frames for each analysis
- Compute results manually
- Create plots for each analysis
- Organize output files

## Expected Output

The benchmark produces:

1. **Console Output**: Real-time progress and final results
2. **benchmark_results.png**: Comparison visualization showing runtime, memory, and LOC
3. **benchmark_summary.txt**: Detailed results summary
4. **analyze_output/**: FastMDAnalysis output directory
5. **mdtraj_output/**: MDTraj analysis results and plots
6. **mdanalysis_output/**: MDAnalysis analysis results and plots

### Sample Results

```
MD Analysis Performance Benchmark Results
======================================================================

Dataset: TrpCage (500 frames with frames=0,-1,10)
Analyses: RMSD, RMSF, RG, Cluster

Results:
----------------------------------------------------------------------

FastMDAnalysis:
  Runtime (computation + plotting): ~16s
  Peak Memory: ~65 MB
  Lines of Code: 1

MDTraj:
  Runtime (computation + plotting): ~2.5s
  Peak Memory: ~20 MB
  Lines of Code: 50

MDAnalysis:
  Runtime (computation + plotting): ~3s
  Peak Memory: ~10 MB
  Lines of Code: 60

----------------------------------------------------------------------

Key Findings:
• FastMDAnalysis: Single CLI command (1 LOC) with automatic workflow
• MDTraj: Manual analysis and plotting (~50 LOC)
• MDAnalysis: Manual analysis and plotting (~60 LOC)
• All measurements include computation + plotting time
• FastMDAnalysis provides simplest API with comparable performance
```

## Interpretation

The benchmark demonstrates that:

### FastMDAnalysis
- **Simplest interface**: 1 line of code using CLI
- **Complete workflow**: Computation + plotting + organization in one command
- **Publication-quality outputs**: Automatic generation of analysis figures
- **Trade-off**: Additional runtime due to comprehensive output generation and organization

### MDTraj
- **Fast computation**: Lightweight and efficient
- **Manual work required**: ~50 lines of code for equivalent functionality
- **No automatic organization**: User must manage outputs

### MDAnalysis
- **Fast computation**: Efficient for basic analyses
- **Most complex**: ~60 lines of code required
- **Different paradigm**: Trajectory iteration model
- **No automatic organization**: User must manage outputs

## Key Insights

1. **Simplicity vs. Performance Trade-off**
   - FastMDAnalysis: Prioritizes ease of use (1 LOC) over raw speed
   - MDTraj/MDAnalysis: Faster but require significantly more code

2. **Time Investment**
   - FastMDAnalysis: Immediate results with publication-quality figures
   - MDTraj/MDAnalysis: Faster computation but requires time to write plotting code

3. **Use Cases**
   - **FastMDAnalysis**: Ideal for exploratory analysis and rapid prototyping
   - **MDTraj**: Best for custom workflows requiring fine-grained control
   - **MDAnalysis**: Best for complex trajectory manipulations

## Notes

- Benchmark uses `tracemalloc` for memory profiling
- Runtime includes both computation and figure generation for all libraries
- Results may vary based on system specifications
- All libraries produce equivalent analysis results

## Extending the Benchmark

To benchmark additional analyses or datasets:

1. Modify the analyses in each benchmark function
2. Update the frame selection with different parameters
3. Use `--atoms` to specify different atom selections
4. Add more datasets by updating the dataset paths

## License

This benchmark script is part of the FastMDAnalysis project and follows the same MIT license.
