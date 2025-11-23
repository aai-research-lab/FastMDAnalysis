#!/usr/bin/env python3
"""
Pure Computation Benchmark - No Plotting Overhead

Compares FastMDAnalysis, MDTraj, and MDAnalysis on pure computation
without any plotting, file I/O, or other overhead.
"""

import sys
import time
import warnings
from pathlib import Path
import numpy as np
import mdtraj as md

# Filter warnings
warnings.filterwarnings('ignore')

# Import FastMDAnalysis
sys.path.insert(0, str(Path(__file__).parent / 'src'))
from fastmdanalysis.datasets import Ubiquitin

# Try importing MDAnalysis
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms as mda_rms
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

def benchmark_fastmda_computation_only():
    """Pure computation benchmark for FastMDA"""
    from fastmdanalysis import FastMDAnalysis
    
    print("\n" + "="*70)
    print("FastMDAnalysis - Pure Computation (No Plotting)")
    print("="*70)
    
    start = time.time()
    
    # Initialize
    fastmda = FastMDAnalysis(
        Ubiquitin.traj,
        Ubiquitin.top,
        frames=(0, -1, 10),
        atoms="protein"
    )
    
    # Run analyses - computation only
    rmsd_result = fastmda.rmsd(ref=0)
    rmsf_result = fastmda.rmsf()
    rg_result = fastmda.rg()
    cluster_result = fastmda.cluster(
        methods=['kmeans', 'dbscan', 'hierarchical'],
        n_clusters=3,
        eps=0.5,
        min_samples=2
    )
    
    runtime = time.time() - start
    
    print(f"Runtime: {runtime:.3f}s")
    print(f"RMSD data shape: {rmsd_result.results['rmsd'].shape}")
    print(f"RMSF data shape: {rmsf_result.results['rmsf'].shape}")
    print(f"RG data shape: {rg_result.results['rg'].shape}")
    print(f"Cluster methods: {list(cluster_result.results.keys())}")
    print("="*70)
    
    return runtime

def benchmark_mdtraj_computation_only():
    """Pure computation benchmark for MDTraj"""
    print("\n" + "="*70)
    print("MDTraj - Pure Computation (No Plotting)")
    print("="*70)
    
    start = time.time()
    
    # Load and prepare trajectory
    traj = md.load(Ubiquitin.traj, top=Ubiquitin.top)
    traj = traj[0::10]  # frames 0,-1,10
    atom_indices = traj.topology.select('protein')
    traj = traj.atom_slice(atom_indices)
    
    # RMSD
    rmsd_data = md.rmsd(traj, traj, frame=0, atom_indices=None)
    
    # RMSF
    avg_xyz = np.mean(traj.xyz, axis=0, keepdims=True)
    ref = md.Trajectory(avg_xyz, traj.topology)
    rmsf_data = md.rmsf(traj, ref)
    
    # RG
    rg_data = md.compute_rg(traj)
    
    # Clustering - compute distance matrix
    rmsd_matrix = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        rmsd_matrix[i] = md.rmsd(traj, traj, frame=i)
    
    # KMeans
    from sklearn.cluster import KMeans, DBSCAN
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    
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
    
    runtime = time.time() - start
    
    print(f"Runtime: {runtime:.3f}s")
    print(f"RMSD data shape: {rmsd_data.shape}")
    print(f"RMSF data shape: {rmsf_data.shape}")
    print(f"RG data shape: {rg_data.shape}")
    print(f"Cluster results: KMeans, DBSCAN, Hierarchical")
    print("="*70)
    
    return runtime

def benchmark_mdanalysis_computation_only():
    """Pure computation benchmark for MDAnalysis"""
    if not HAS_MDANALYSIS:
        print("\nMDAnalysis not available - skipping")
        return None
        
    print("\n" + "="*70)
    print("MDAnalysis - Pure Computation (No Plotting)")
    print("="*70)
    
    start = time.time()
    
    # Load trajectory
    u = mda.Universe(Ubiquitin.top, Ubiquitin.traj)
    protein = u.select_atoms('protein')
    
    # Frame selection
    frame_list = list(range(0, len(u.trajectory), 10))
    
    # RMSD
    rmsd_results = []
    ref_coords = protein.positions.copy()
    for ts in u.trajectory[frame_list]:
        current_coords = protein.positions
        rmsd_val = mda_rms.rmsd(current_coords, ref_coords, center=True) / 10.0
        rmsd_results.append(rmsd_val)
    rmsd_data = np.array(rmsd_results)
    
    # RMSF
    coordinates = []
    for ts in u.trajectory[frame_list]:
        coordinates.append(protein.positions.copy())
    coordinates = np.array(coordinates)
    avg_coords = np.mean(coordinates, axis=0)
    rmsf_data = np.sqrt(np.mean((coordinates - avg_coords) ** 2, axis=0))
    rmsf_data = np.linalg.norm(rmsf_data, axis=1) / 10.0
    
    # RG
    rg_results = []
    for ts in u.trajectory[frame_list]:
        rg_val = protein.radius_of_gyration() / 10.0
        rg_results.append(rg_val)
    rg_data = np.array(rg_results)
    
    # Clustering
    from sklearn.cluster import KMeans, DBSCAN
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist, squareform
    
    coords_flat = coordinates.reshape(len(frame_list), -1)
    
    # KMeans
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(coords_flat)
    
    # DBSCAN
    distances = pdist(coords_flat)
    dist_matrix = squareform(distances)
    dbscan = DBSCAN(eps=50.0, min_samples=2, metric='precomputed')
    dbscan_labels = dbscan.fit_predict(dist_matrix)
    
    # Hierarchical
    linkage_matrix = linkage(distances, method='ward')
    
    runtime = time.time() - start
    
    print(f"Runtime: {runtime:.3f}s")
    print(f"RMSD data shape: {rmsd_data.shape}")
    print(f"RMSF data shape: {rmsf_data.shape}")
    print(f"RG data shape: {rg_data.shape}")
    print(f"Cluster results: KMeans, DBSCAN, Hierarchical")
    print("="*70)
    
    return runtime

def main():
    print("="*70)
    print("PURE COMPUTATION BENCHMARK - NO PLOTTING OVERHEAD")
    print("="*70)
    print("Dataset: Ubiquitin Q99, frames 0:-1:10")
    print("Analyses: RMSD, RMSF, RG, Cluster (KMeans, DBSCAN, Hierarchical)")
    print("="*70)
    
    results = {}
    
    # Run benchmarks
    results['FastMDAnalysis'] = benchmark_fastmda_computation_only()
    results['MDTraj'] = benchmark_mdtraj_computation_only()
    results['MDAnalysis'] = benchmark_mdanalysis_computation_only()
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY - Pure Computation Times")
    print("="*70)
    for name, runtime in results.items():
        if runtime is not None:
            print(f"{name:20s}: {runtime:6.3f}s")
    print("="*70)
    
    # Analysis
    if results['FastMDAnalysis'] and results['MDTraj']:
        ratio = results['FastMDAnalysis'] / results['MDTraj']
        print(f"\nFastMDA/MDTraj ratio: {ratio:.2f}x")
        print("\nNote: FastMDA uses MDTraj backend, so they should be very close.")
        print("Any difference is due to Python wrapper overhead and initialization.")

if __name__ == '__main__':
    main()
