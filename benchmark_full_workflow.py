#!/usr/bin/env python3
"""
FastMDAnalysis Full Workflow Benchmark

Benchmarks FastMDAnalysis, MDTraj, and MDAnalysis on COMPLETE workflows
including computation and visualization (but excluding slides).

Measures:
- Runtime (wall clock time)
- Peak memory usage

Runs each workflow 5 times and reports average ± standard deviation.
"""

import sys
import time
import warnings
import os
import tempfile
import shutil
from pathlib import Path
import numpy as np
import psutil
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import mdtraj as md
from sklearn.cluster import KMeans, DBSCAN
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform, pdist

warnings.filterwarnings('ignore')

# Add src to path for FastMDA
sys.path.insert(0, str(Path(__file__).parent / 'src'))
from fastmdanalysis.datasets import TrpCage

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms as mda_rms
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    print("Warning: MDAnalysis not installed. Will skip MDAnalysis benchmarks.")

# Number of iterations for averaging
NUM_ITERATIONS = 5

# Tool colors for consistent visualization
TOOL_COLORS = {
    "fastmdanalysis": "#4472C4",
    "mdtraj": "#ED7D31",
    "mdanalysis": "#A5A5A5"
}


class MemoryMonitor:
    """Monitor peak memory usage during execution"""
    def __init__(self):
        self.process = psutil.Process()
        self.baseline_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.peak_memory = self.baseline_memory
        
    def update(self):
        """Update peak memory if current usage is higher"""
        current_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        if current_memory > self.peak_memory:
            self.peak_memory = current_memory
    
    def get_peak_mb(self):
        """Get peak memory usage in MB (relative to baseline)"""
        return self.peak_memory - self.baseline_memory


def benchmark_fastmda_single(output_dir):
    """FastMDA full workflow - single run with basic figures only"""
    from fastmdanalysis import FastMDAnalysis
    
    mem_monitor = MemoryMonitor()
    start = time.time()
    
    # Initialize FastMDA
    fastmda = FastMDAnalysis(
        TrpCage.traj, 
        TrpCage.top, 
        frames=(0, -1, 10),
        atoms="protein"
    )
    mem_monitor.update()
    
    # Run analyses - FastMDA automatically generates ONE plot per analysis
    # RMSD: generates rmsd.png (line plot)
    rmsd_result = fastmda.rmsd(ref=0)
    mem_monitor.update()
    
    # RMSF: generates rmsf.png (bar plot)
    rmsf_result = fastmda.rmsf()
    mem_monitor.update()
    
    # RG: generates rg.png (line plot)
    rg_result = fastmda.rg()
    mem_monitor.update()
    
    # Cluster: generates basic cluster plots (one per method)
    # We'll just run it without extra plots - FastMDA generates compact outputs
    cluster_result = fastmda.cluster(
        methods=['kmeans', 'dbscan', 'hierarchical'],
        n_clusters=3,
        eps=0.5,
        min_samples=2
    )
    mem_monitor.update()
    
    runtime = time.time() - start
    peak_memory = mem_monitor.get_peak_mb()
    
    return runtime, peak_memory


def benchmark_mdtraj_single(output_dir):
    """MDTraj full workflow - single run with basic figures only (matching FastMDA)"""
    mem_monitor = MemoryMonitor()
    start = time.time()
    
    # Load and prepare trajectory
    traj = md.load(TrpCage.traj, top=TrpCage.top)
    traj = traj[0::10]
    atom_indices = traj.topology.select('protein')
    traj = traj.atom_slice(atom_indices)
    mem_monitor.update()
    
    # RMSD computation and plot (1 figure)
    rmsd_data = md.rmsd(traj, traj, frame=0)
    mem_monitor.update()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rmsd_data, marker='o', linestyle='-')
    ax.set_xlabel('Frame')
    ax.set_ylabel('RMSD (nm)')
    ax.set_title('RMSD vs Frame (ref=0, align=True)')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsd.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    # RMSF computation and plot (1 figure)
    avg_xyz = np.mean(traj.xyz, axis=0, keepdims=True)
    ref = md.Trajectory(avg_xyz, traj.topology)
    rmsf_data = md.rmsf(traj, ref)
    mem_monitor.update()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(rmsf_data)), rmsf_data, width=0.9)
    ax.set_xlabel('Atom Index')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title('RMSF per Atom')
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsf.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    # RG computation and plot (1 figure)
    rg_data = md.compute_rg(traj)
    mem_monitor.update()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rg_data, marker='o', linestyle='-')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Radius of Gyration (nm)')
    ax.set_title('Radius of Gyration vs Frame')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rg.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    # Clustering computation (3 methods) and plot (1 combined figure)
    rmsd_matrix = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        rmsd_matrix[i] = md.rmsd(traj, traj, frame=i)
    mem_monitor.update()
    
    # KMeans
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(rmsd_matrix)
    
    # DBSCAN
    dbscan = DBSCAN(eps=0.5, min_samples=2, metric='precomputed')
    dbscan_labels = dbscan.fit_predict(rmsd_matrix)
    
    # Hierarchical
    rmsd_matrix_sym = (rmsd_matrix + rmsd_matrix.T) / 2
    np.fill_diagonal(rmsd_matrix_sym, 0)
    condensed_dist = squareform(rmsd_matrix_sym)
    linkage_matrix = linkage(condensed_dist, method='ward')
    hierarchical_labels = fcluster(linkage_matrix, 3, criterion='maxclust')
    mem_monitor.update()
    
    # Single clustering plot with 3 subplots (1 figure total)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    axes[0].scatter(range(len(kmeans_labels)), kmeans_labels, c=kmeans_labels, cmap='viridis', s=50)
    axes[0].set_title('KMeans')
    axes[0].set_xlabel('Frame')
    axes[0].set_ylabel('Cluster')
    axes[0].grid(alpha=0.3)
    
    axes[1].scatter(range(len(dbscan_labels)), dbscan_labels, c=dbscan_labels, cmap='viridis', s=50)
    axes[1].set_title('DBSCAN')
    axes[1].set_xlabel('Frame')
    axes[1].set_ylabel('Cluster')
    axes[1].grid(alpha=0.3)
    
    axes[2].scatter(range(len(hierarchical_labels)), hierarchical_labels, c=hierarchical_labels, cmap='viridis', s=50)
    axes[2].set_title('Hierarchical')
    axes[2].set_xlabel('Frame')
    axes[2].set_ylabel('Cluster')
    axes[2].grid(alpha=0.3)
    
    plt.suptitle('Clustering Results', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cluster.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    runtime = time.time() - start
    peak_memory = mem_monitor.get_peak_mb()
    
    return runtime, peak_memory


def benchmark_mdanalysis_single(output_dir):
    """MDAnalysis full workflow - single run with basic figures only (matching FastMDA)"""
    if not HAS_MDANALYSIS:
        return None, None
    
    mem_monitor = MemoryMonitor()
    start = time.time()
    
    # Load trajectory
    u = mda.Universe(TrpCage.top, TrpCage.traj)
    protein = u.select_atoms('protein')
    frame_list = list(range(0, len(u.trajectory), 10))
    mem_monitor.update()
    
    # RMSD computation and plot (1 figure)
    rmsd_results = []
    ref_coords = protein.positions.copy()
    for ts in u.trajectory[frame_list]:
        rmsd_val = mda_rms.rmsd(protein.positions, ref_coords, center=True) / 10.0
        rmsd_results.append(rmsd_val)
    rmsd_data = np.array(rmsd_results)
    mem_monitor.update()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rmsd_data, marker='o', linestyle='-')
    ax.set_xlabel('Frame')
    ax.set_ylabel('RMSD (nm)')
    ax.set_title('RMSD vs Frame (ref=0, align=True)')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsd.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    # RMSF computation and plot (1 figure)
    coordinates = []
    for ts in u.trajectory[frame_list]:
        coordinates.append(protein.positions.copy())
    coordinates = np.array(coordinates)
    avg_coords = np.mean(coordinates, axis=0)
    rmsf_data = np.sqrt(np.mean((coordinates - avg_coords) ** 2, axis=0))
    rmsf_data = np.linalg.norm(rmsf_data, axis=1) / 10.0
    mem_monitor.update()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(rmsf_data)), rmsf_data, width=0.9)
    ax.set_xlabel('Atom Index')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title('RMSF per Atom')
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rmsf.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    # RG computation and plot (1 figure)
    rg_results = []
    for ts in u.trajectory[frame_list]:
        rg_results.append(protein.radius_of_gyration() / 10.0)
    rg_data = np.array(rg_results)
    mem_monitor.update()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rg_data, marker='o', linestyle='-')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Radius of Gyration (nm)')
    ax.set_title('Radius of Gyration vs Frame')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'rg.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    # Clustering computation (3 methods) and plot (1 combined figure)
    coords_flat = coordinates.reshape(len(frame_list), -1)
    
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(coords_flat)
    
    distances = pdist(coords_flat)
    dist_matrix = squareform(distances)
    dbscan = DBSCAN(eps=50.0, min_samples=2, metric='precomputed')
    dbscan_labels = dbscan.fit_predict(dist_matrix)
    
    linkage_matrix = linkage(distances, method='ward')
    hierarchical_labels = fcluster(linkage_matrix, 3, criterion='maxclust')
    mem_monitor.update()
    
    # Single clustering plot with 3 subplots (1 figure total)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    axes[0].scatter(range(len(kmeans_labels)), kmeans_labels, c=kmeans_labels, cmap='viridis', s=50)
    axes[0].set_title('KMeans')
    axes[0].set_xlabel('Frame')
    axes[0].set_ylabel('Cluster')
    axes[0].grid(alpha=0.3)
    
    axes[1].scatter(range(len(dbscan_labels)), dbscan_labels, c=dbscan_labels, cmap='viridis', s=50)
    axes[1].set_title('DBSCAN')
    axes[1].set_xlabel('Frame')
    axes[1].set_ylabel('Cluster')
    axes[1].grid(alpha=0.3)
    
    axes[2].scatter(range(len(hierarchical_labels)), hierarchical_labels, c=hierarchical_labels, cmap='viridis', s=50)
    axes[2].set_title('Hierarchical')
    axes[2].set_xlabel('Frame')
    axes[2].set_ylabel('Cluster')
    axes[2].grid(alpha=0.3)
    
    plt.suptitle('Clustering Results', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cluster.png'), dpi=300, bbox_inches='tight')
    plt.close()
    mem_monitor.update()
    
    runtime = time.time() - start
    peak_memory = mem_monitor.get_peak_mb()
    
    return runtime, peak_memory


def benchmark_fastmda():
    """Run FastMDA workflow multiple times"""
    print("\n" + "="*70)
    print("FastMDAnalysis - Full Workflow (with figures, no slides)")
    print("="*70)
    print(f"Running {NUM_ITERATIONS} iterations...")
    
    runtimes = []
    memories = []
    
    for i in range(NUM_ITERATIONS):
        # Create temporary directory for output
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime, memory = benchmark_fastmda_single(tmpdir)
            runtimes.append(runtime)
            memories.append(memory)
            print(f"  Iteration {i+1}/{NUM_ITERATIONS}: {runtime:.3f}s, {memory:.1f} MB")
    
    avg_runtime = np.mean(runtimes)
    std_runtime = np.std(runtimes)
    avg_memory = np.mean(memories)
    std_memory = np.std(memories)
    
    print(f"\nAverage Runtime: {avg_runtime:.3f}s (±{std_runtime:.3f}s)")
    print(f"Average Peak Memory: {avg_memory:.1f} MB (±{std_memory:.1f} MB)")
    print("="*70)
    
    return {
        'name': 'FastMDAnalysis',
        'runtime_avg': avg_runtime,
        'runtime_std': std_runtime,
        'memory_avg': avg_memory,
        'memory_std': std_memory
    }


def benchmark_mdtraj():
    """Run MDTraj workflow multiple times"""
    print("\n" + "="*70)
    print("MDTraj - Full Workflow (with figures)")
    print("="*70)
    print(f"Running {NUM_ITERATIONS} iterations...")
    
    runtimes = []
    memories = []
    
    for i in range(NUM_ITERATIONS):
        # Create temporary directory for output
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime, memory = benchmark_mdtraj_single(tmpdir)
            runtimes.append(runtime)
            memories.append(memory)
            print(f"  Iteration {i+1}/{NUM_ITERATIONS}: {runtime:.3f}s, {memory:.1f} MB")
    
    avg_runtime = np.mean(runtimes)
    std_runtime = np.std(runtimes)
    avg_memory = np.mean(memories)
    std_memory = np.std(memories)
    
    print(f"\nAverage Runtime: {avg_runtime:.3f}s (±{std_runtime:.3f}s)")
    print(f"Average Peak Memory: {avg_memory:.1f} MB (±{std_memory:.1f} MB)")
    print("="*70)
    
    return {
        'name': 'MDTraj',
        'runtime_avg': avg_runtime,
        'runtime_std': std_runtime,
        'memory_avg': avg_memory,
        'memory_std': std_memory
    }


def benchmark_mdanalysis():
    """Run MDAnalysis workflow multiple times"""
    if not HAS_MDANALYSIS:
        print("\nMDAnalysis not available - skipping")
        return None
    
    print("\n" + "="*70)
    print("MDAnalysis - Full Workflow (with figures)")
    print("="*70)
    print(f"Running {NUM_ITERATIONS} iterations...")
    
    runtimes = []
    memories = []
    
    for i in range(NUM_ITERATIONS):
        # Create temporary directory for output
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime, memory = benchmark_mdanalysis_single(tmpdir)
            if runtime is None:
                return None
            runtimes.append(runtime)
            memories.append(memory)
            print(f"  Iteration {i+1}/{NUM_ITERATIONS}: {runtime:.3f}s, {memory:.1f} MB")
    
    avg_runtime = np.mean(runtimes)
    std_runtime = np.std(runtimes)
    avg_memory = np.mean(memories)
    std_memory = np.std(memories)
    
    print(f"\nAverage Runtime: {avg_runtime:.3f}s (±{std_runtime:.3f}s)")
    print(f"Average Peak Memory: {avg_memory:.1f} MB (±{std_memory:.1f} MB)")
    print("="*70)
    
    return {
        'name': 'MDAnalysis',
        'runtime_avg': avg_runtime,
        'runtime_std': std_runtime,
        'memory_avg': avg_memory,
        'memory_std': std_memory
    }


def create_benchmark_plots(results):
    """Create benchmark comparison plots"""
    print("\n" + "="*70)
    print("GENERATING PLOTS")
    print("="*70)
    
    # Extract data
    names = [r['name'] for r in results]
    runtime_avgs = [r['runtime_avg'] for r in results]
    runtime_stds = [r['runtime_std'] for r in results]
    memory_avgs = [r['memory_avg'] for r in results]
    memory_stds = [r['memory_std'] for r in results]
    
    # Get colors
    colors = [TOOL_COLORS[name.lower().replace(" ", "")] for name in names]
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Runtime plot
    ax1 = axes[0]
    bars1 = ax1.bar(names, runtime_avgs, yerr=runtime_stds, 
                     color=colors, alpha=0.9, edgecolor='black', 
                     linewidth=1.5, capsize=5, error_kw={'linewidth': 2})
    ax1.set_ylabel('Runtime (seconds)', fontsize=14)
    ax1.set_title('Runtime Comparison', fontsize=16, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels
    for bar, runtime_avg, runtime_std in zip(bars1, runtime_avgs, runtime_stds):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + runtime_std,
                f'{runtime_avg:.2f}s\n±{runtime_std:.2f}s',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Memory plot
    ax2 = axes[1]
    bars2 = ax2.bar(names, memory_avgs, yerr=memory_stds, 
                     color=colors, alpha=0.9, edgecolor='black', 
                     linewidth=1.5, capsize=5, error_kw={'linewidth': 2})
    ax2.set_ylabel('Peak Memory (MB)', fontsize=14)
    ax2.set_title('Peak Memory Comparison', fontsize=16, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels
    for bar, memory_avg, memory_std in zip(bars2, memory_avgs, memory_stds):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + memory_std,
                f'{memory_avg:.1f} MB\n±{memory_std:.1f} MB',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.suptitle(f'FastMDAnalysis Full Workflow Benchmark\n'
                 f'TrpCage (500 frames) - RMSD, RMSF, RG, Cluster (with figures)\n'
                 f'Averaged over {NUM_ITERATIONS} iterations',
                 fontweight='bold', fontsize=14)
    plt.tight_layout()
    
    output_file = 'benchmark_full_workflow.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()
    
    return output_file


def create_summary_table(results):
    """Create summary table visualization"""
    print("Creating summary table...")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Prepare table data
    table_data = [['Library', 'Runtime (avg ± std)', 'Memory (avg ± std)', 'Speedup']]
    
    baseline_runtime = results[1]['runtime_avg'] if len(results) > 1 else results[0]['runtime_avg']
    
    for r in results:
        runtime_str = f"{r['runtime_avg']:.2f}s ± {r['runtime_std']:.2f}s"
        memory_str = f"{r['memory_avg']:.1f} MB ± {r['memory_std']:.1f} MB"
        speedup = baseline_runtime / r['runtime_avg']
        speedup_str = f"{speedup:.2f}x"
        table_data.append([r['name'], runtime_str, memory_str, speedup_str])
    
    # Create table
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.25, 0.3, 0.3, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)
    
    # Style header row
    for i in range(4):
        cell = table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white', fontsize=12)
    
    # Style data rows
    tool_keys = ['fastmdanalysis', 'mdtraj', 'mdanalysis']
    for i in range(1, len(table_data)):
        tool_key = tool_keys[i-1] if i-1 < len(tool_keys) else tool_keys[-1]
        row_color = TOOL_COLORS.get(tool_key, '#CCCCCC')
        
        for j in range(4):
            cell = table[(i, j)]
            cell.set_facecolor(row_color + '80')  # Add transparency
            cell.set_edgecolor('black')
            cell.set_linewidth(1.5)
    
    ax.axis('off')
    ax.set_title(f'FastMDAnalysis Full Workflow Benchmark Summary\n'
                 f'Averaged over {NUM_ITERATIONS} iterations',
                 fontweight='bold', fontsize=14, pad=20)
    
    # Add footer
    footer_text = (
        f'Dataset: TrpCage, 500 frames (frames 0,-1,10)\n'
        f'Analyses: RMSD, RMSF, RG, Cluster (KMeans, DBSCAN, Hierarchical)\n'
        f'Workflow: Complete analysis with figure generation (no slides)'
    )
    fig.text(0.5, 0.05, footer_text, ha='center', fontsize=10, style='italic',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    output_file = 'benchmark_summary_table.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()
    
    return output_file


def main():
    print("="*70)
    print("FASTMDANALYSIS FULL WORKFLOW BENCHMARK")
    print("="*70)
    print("Dataset: TrpCage, 500 frames (frames 0,-1,10)")
    print("Analyses: RMSD, RMSF, RG, Cluster (KMeans, DBSCAN, Hierarchical)")
    print("Workflow: Complete with figure generation (excluding slides)")
    print(f"Iterations: {NUM_ITERATIONS} runs per library")
    print("Metrics: Runtime (wall clock) and Peak Memory Usage")
    print("="*70)
    
    results = []
    
    # Run benchmarks
    results.append(benchmark_fastmda())
    results.append(benchmark_mdtraj())
    mda_result = benchmark_mdanalysis()
    if mda_result:
        results.append(mda_result)
    
    # Print summary
    print("\n" + "="*70)
    print(f"SUMMARY (AVERAGED OVER {NUM_ITERATIONS} ITERATIONS)")
    print("="*70)
    for r in results:
        print(f"{r['name']:20s}:")
        print(f"  Runtime: {r['runtime_avg']:6.2f}s ± {r['runtime_std']:5.2f}s")
        print(f"  Memory:  {r['memory_avg']:6.1f} MB ± {r['memory_std']:5.1f} MB")
    print("="*70)
    
    # Calculate speedup
    if len(results) >= 2:
        fastmda_time = results[0]['runtime_avg']
        mdtraj_time = results[1]['runtime_avg']
        ratio = fastmda_time / mdtraj_time
        print(f"\nFastMDA/MDTraj runtime ratio: {ratio:.2f}x")
        
        fastmda_mem = results[0]['memory_avg']
        mdtraj_mem = results[1]['memory_avg']
        mem_ratio = fastmda_mem / mdtraj_mem
        print(f"FastMDA/MDTraj memory ratio: {mem_ratio:.2f}x")
    
    # Generate plots
    create_benchmark_plots(results)
    create_summary_table(results)
    
    print("\n" + "="*70)
    print("BENCHMARK COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  • benchmark_full_workflow.png - Runtime and memory comparison")
    print("  • benchmark_summary_table.png - Detailed summary table")
    print("="*70)


if __name__ == '__main__':
    main()
