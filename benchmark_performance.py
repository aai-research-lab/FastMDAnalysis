#!/usr/bin/env python3
"""
FastMDAnalysis Performance Benchmark

This script benchmarks FastMDAnalysis using the CLI with RMSD, RMSF, RG, and Cluster
analyses on the TrpCage dataset with 500 frames (frames 0,-1,10).

It measures:
- Total runtime (computation + plotting)
- Peak memory usage
- Lines of code (LOC = 1, using CLI command)

The benchmark runs the FastMDAnalysis CLI command and then creates custom benchmark plots.

Usage:
    python benchmark_performance.py
"""

import sys
import time
import warnings
from pathlib import Path
import tracemalloc
import shutil

import numpy as np
import matplotlib.pyplot as plt

# Filter out benign warnings
warnings.filterwarnings('ignore', message='Unlikely unit cell vectors detected')
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Import FastMDAnalysis components for dataset paths
try:
    sys.path.insert(0, str(Path(__file__).parent / 'src'))
    from fastmdanalysis.datasets import TrpCage
except ImportError as e:
    print(f"Error importing FastMDAnalysis: {e}", file=sys.stderr)
    print("Make sure FastMDAnalysis is installed or src is in PYTHONPATH", file=sys.stderr)
    sys.exit(1)


def format_memory(bytes_val):
    """Format memory in human-readable form."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} TB"


def format_time(seconds):
    """Format time in human-readable form."""
    if seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = seconds % 60
        return f"{mins}m {secs:.2f}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}h {mins}m {secs:.2f}s"


def run_fastmda_benchmark():
    """
    Run FastMDAnalysis benchmark using CLI approach.
    
    This represents a single line of code usage:
    fastmda analyze -traj <traj.dcd> -top <top.pdb> --frames 0,-1,10 --include cluster rmsd rg rmsf
    """
    # Define CLI command arguments
    cli_args = {
        'traj': TrpCage.traj,
        'top': TrpCage.top,
        'frames': '0,-1,10',
        'include': ['cluster', 'rmsd', 'rg', 'rmsf']
    }
    
    print("\n" + "="*70)
    print("FastMDAnalysis Performance Benchmark")
    print("="*70)
    print(f"Dataset: TrpCage")
    print(f"Frame selection: {cli_args['frames']} -> ~500 frames")
    print(f"Analyses: RMSD, RMSF, RG, Cluster")
    print(f"Command (1 LOC): fastmda analyze -traj {cli_args['traj']} -top {cli_args['top']} --frames {cli_args['frames']} --include {' '.join(cli_args['include'])}")
    print("="*70)
    
    # Clean up any previous output
    output_dir = Path("analyze_output")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    
    # Import CLI function
    from fastmdanalysis.cli.main import main as cli_main
    
    # Start tracking
    tracemalloc.start()
    start_time = time.time()
    
    print("\nRunning FastMDAnalysis CLI command...")
    print("-" * 70)
    
    # Simulate CLI invocation by setting sys.argv
    original_argv = sys.argv
    try:
        sys.argv = [
            'fastmda', 'analyze',
            '-traj', cli_args['traj'],
            '-top', cli_args['top'],
            '--frames', cli_args['frames'],
            '--include'
        ] + cli_args['include']
        
        # Run the CLI
        cli_main()
        
    except SystemExit as e:
        if e.code != 0:
            print(f"✗ FastMDAnalysis command failed with exit code {e.code}")
            return None
    except Exception as e:
        print(f"✗ FastMDAnalysis command failed: {e}")
        import traceback
        traceback.print_exc()
        return None
    finally:
        sys.argv = original_argv
    
    # Stop tracking
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    runtime = end_time - start_time
    
    print("-" * 70)
    print(f"✓ FastMDAnalysis completed successfully")
    print(f"  Runtime (computation + plotting): {format_time(runtime)}")
    print(f"  Peak Memory: {format_memory(peak)}")
    print(f"  Lines of Code: 1")
    print("="*70)
    
    return {
        'runtime': runtime,
        'memory_peak': peak,
        'loc': 1,
        'success': True
    }


def create_benchmark_plots(result):
    """
    Create custom benchmark visualization plots.
    """
    if result is None or not result['success']:
        return
    
    print("\nCreating benchmark visualization plots...")
    
    # Create a figure with benchmark results
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('FastMDAnalysis Performance Benchmark\n(TrpCage, 500 frames, RMSD + RMSF + RG + Cluster)', 
                 fontsize=14, fontweight='bold')
    
    # Runtime plot
    ax1 = axes[0]
    ax1.bar(['FastMDAnalysis'], [result['runtime']], color='#2E86AB', alpha=0.8)
    ax1.set_ylabel('Runtime (seconds)', fontsize=12)
    ax1.set_title('Total Runtime\n(Computation + Plotting)', fontsize=11, fontweight='bold')
    ax1.text(0, result['runtime'] + result['runtime']*0.05, 
             f"{format_time(result['runtime'])}", 
             ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax1.set_ylim(0, result['runtime'] * 1.2)
    ax1.grid(axis='y', alpha=0.3)
    
    # Memory plot
    ax2 = axes[1]
    memory_mb = result['memory_peak'] / (1024 * 1024)
    ax2.bar(['FastMDAnalysis'], [memory_mb], color='#A23B72', alpha=0.8)
    ax2.set_ylabel('Peak Memory (MB)', fontsize=12)
    ax2.set_title('Peak Memory Usage', fontsize=11, fontweight='bold')
    ax2.text(0, memory_mb + memory_mb*0.05, 
             f"{format_memory(result['memory_peak'])}", 
             ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax2.set_ylim(0, memory_mb * 1.2)
    ax2.grid(axis='y', alpha=0.3)
    
    # LOC plot
    ax3 = axes[2]
    ax3.bar(['FastMDAnalysis'], [result['loc']], color='#F18F01', alpha=0.8)
    ax3.set_ylabel('Lines of Code', fontsize=12)
    ax3.set_title('Code Complexity', fontsize=11, fontweight='bold')
    ax3.text(0, result['loc'] + 0.1, 
             f"{result['loc']} LOC", 
             ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax3.set_ylim(0, 2)
    ax3.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = 'benchmark_results.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Benchmark plot saved to: {output_file}")
    plt.close()
    
    # Create a summary text file
    summary_file = 'benchmark_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("FastMDAnalysis Performance Benchmark Results\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Dataset: TrpCage (500 frames with frames=0,-1,10)\n")
        f.write(f"Analyses: RMSD, RMSF, RG, Cluster\n")
        f.write(f"CLI Command: fastmda analyze -traj <traj.dcd> -top <top.pdb> --frames 0,-1,10 --include cluster rmsd rg rmsf\n\n")
        f.write("Results:\n")
        f.write("-" * 70 + "\n")
        f.write(f"Runtime (computation + plotting): {format_time(result['runtime'])}\n")
        f.write(f"Peak Memory: {format_memory(result['memory_peak'])}\n")
        f.write(f"Lines of Code: {result['loc']}\n")
        f.write("-" * 70 + "\n\n")
        f.write("Key Findings:\n")
        f.write("• Single-line CLI command provides complete analysis workflow\n")
        f.write("• Includes automatic computation, plotting, and output organization\n")
        f.write("• Total time includes both computation and figure generation\n")
        f.write("• Ideal for rapid exploratory analysis and publication-quality outputs\n")
    
    print(f"✓ Benchmark summary saved to: {summary_file}")


def main():
    """Main benchmark function."""
    # Run the benchmark
    result = run_fastmda_benchmark()
    
    if result is None:
        print("\n✗ Benchmark failed")
        sys.exit(1)
    
    # Create visualization plots
    create_benchmark_plots(result)
    
    print("\n" + "="*70)
    print("Benchmark Complete!")
    print("="*70)
    print(f"\nResults:")
    print(f"  Runtime: {format_time(result['runtime'])}")
    print(f"  Memory: {format_memory(result['memory_peak'])}")
    print(f"  LOC: {result['loc']}")
    print("\nOutput files:")
    print("  - benchmark_results.png (visualization)")
    print("  - benchmark_summary.txt (detailed results)")
    print("  - analyze_output/ (FastMDAnalysis output directory)")
    print("="*70)


if __name__ == '__main__':
    main()
