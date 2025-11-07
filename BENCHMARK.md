# FastMDAnalysis Performance Benchmark

This benchmark evaluates FastMDAnalysis performance using the CLI on molecular dynamics trajectory analysis tasks.

## Overview

The benchmark measures FastMDAnalysis performance on the following metrics:
- **Runtime**: Total execution time including computation and plotting
- **Memory**: Peak memory usage during execution  
- **Lines of Code (LOC)**: Code complexity (1 LOC using CLI)

## Dataset

- **TrpCage**: 500 frames selected using `--frames 0,-1,10` from 4999 total frames
- **Atom Selection**: All atoms (default)

## Analyses Performed

The benchmark runs the following analyses:
- **RMSD** (Root Mean Square Deviation)
- **RMSF** (Root Mean Square Fluctuation)
- **RG** (Radius of Gyration)
- **Cluster** (KMeans, DBSCAN, Hierarchical)

## Usage

### Prerequisites

```bash
pip install mdtraj numpy matplotlib scikit-learn scipy python-pptx Pillow cairosvg PyYAML
```

### Running the Benchmark

```bash
# Set PYTHONPATH to include FastMDAnalysis src directory (if not installed)
export PYTHONPATH=/path/to/FastMDAnalysis/src:$PYTHONPATH

# Run the benchmark
python benchmark_performance.py
```

## CLI Command

The benchmark executes the following single-line command (1 LOC):

```bash
fastmda analyze -traj <traj.dcd> -top <top.pdb> --frames 0,-1,10 --include cluster rmsd rg rmsf
```

Where `<traj.dcd>` and `<top.pdb>` are the trajectory and topology files (e.g., from TrpCage dataset).

## Expected Output

The benchmark produces:

1. **Console Output**: Real-time progress and final results
2. **benchmark_results.png**: Visualization showing runtime, memory, and LOC
3. **benchmark_summary.txt**: Detailed results summary
4. **analyze_output/**: FastMDAnalysis output directory with analysis results and figures

### Sample Results

```
FastMDAnalysis Performance Benchmark Results
======================================================================

Dataset: TrpCage (500 frames with frames=0,-1,10)
Analyses: RMSD, RMSF, RG, Cluster
CLI Command: fastmda analyze -traj <traj.dcd> -top <top.pdb> --frames 0,-1,10 --include cluster rmsd rg rmsf

Results:
----------------------------------------------------------------------
Runtime (computation + plotting): ~30-40s
Peak Memory: ~200-300 MB
Lines of Code: 1
----------------------------------------------------------------------

Key Findings:
• Single-line CLI command provides complete analysis workflow
• Includes automatic computation, plotting, and output organization
• Total time includes both computation and figure generation
• Ideal for rapid exploratory analysis and publication-quality outputs
```

## Interpretation

The benchmark demonstrates that FastMDAnalysis provides:
- **Simplest possible interface**: 1 line of code using CLI
- **Complete workflow**: Computation + plotting + organization in one command
- **Publication-quality outputs**: Automatic generation of analysis figures
- **Efficient execution**: Reasonable runtime and memory for comprehensive analysis

The single CLI command replaces what would typically require:
- Loading trajectory and topology
- Configuring each analysis separately
- Running each analysis manually
- Creating plots for each analysis
- Organizing output files
- Typically 50-100+ lines of code in traditional approaches

## Notes

- Benchmark uses `tracemalloc` for memory profiling
- Runtime includes both computation and figure generation
- Results may vary based on system specifications
- The benchmark focuses on FastMDAnalysis only (no comparison with other libraries)
- Custom plotting script creates visualization of benchmark results

## Extending the Benchmark

To benchmark additional analyses or datasets:

1. Modify the `--include` flag in the CLI command to add/remove analyses
2. Update the frame selection with `--frames` parameter
3. Use `--atoms` to specify atom selection
4. Modify dataset in the script by changing `TrpCage` to other available datasets

Example with different parameters:
```bash
fastmda analyze -traj <traj.dcd> -top <top.pdb> --frames 0,1000,5 --atoms "protein" --include rmsd rmsf rg hbonds ss sasa
```

## License

This benchmark script is part of the FastMDAnalysis project and follows the same MIT license.
