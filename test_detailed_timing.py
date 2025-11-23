import sys, time
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
sys.path.insert(0, str(Path('.') / 'src'))
from fastmdanalysis import FastMDAnalysis
from fastmdanalysis.datasets import Ubiquitin
import warnings
warnings.filterwarnings('ignore')

total_start = time.time()

print("Initializing...")
start = time.time()
fastmda = FastMDAnalysis(Ubiquitin.traj, Ubiquitin.top, frames=(0, -1, 10), atoms="protein")
print(f"  Init: {time.time() - start:.2f}s")

print("RMSD analysis...")
start = time.time()
rmsd_result = fastmda.rmsd(ref=0)
print(f"  Computation: {time.time() - start:.2f}s")
start = time.time()
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(rmsd_result.results['rmsd'])
plt.savefig('/tmp/rmsd.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Plotting: {time.time() - start:.2f}s")

print("RMSF analysis...")
start = time.time()
rmsf_result = fastmda.rmsf()
print(f"  Computation: {time.time() - start:.2f}s")
start = time.time()
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(rmsf_result.results['rmsf'])
plt.savefig('/tmp/rmsf.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Plotting: {time.time() - start:.2f}s")

print("RG analysis...")
start = time.time()
rg_result = fastmda.rg()
print(f"  Computation: {time.time() - start:.2f}s")
start = time.time()
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(rg_result.results['rg'])
plt.savefig('/tmp/rg.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Plotting: {time.time() - start:.2f}s")

print("Cluster analysis...")
start = time.time()
cluster_result = fastmda.cluster(
    methods=['kmeans', 'dbscan', 'hierarchical'],
    n_clusters=3,
    eps=0.5,
    min_samples=2
)
print(f"  Computation: {time.time() - start:.2f}s")

start = time.time()
if 'kmeans' in cluster_result.results:
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = cluster_result.results['kmeans']['labels']
    ax.scatter(range(len(labels)), labels, c=labels, cmap='viridis', alpha=0.6)
    plt.savefig('/tmp/cluster_kmeans.png', dpi=150, bbox_inches='tight')
    plt.close()
if 'dbscan' in cluster_result.results:
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = cluster_result.results['dbscan']['labels']
    ax.scatter(range(len(labels)), labels, c=labels, cmap='viridis', alpha=0.6)
    plt.savefig('/tmp/cluster_dbscan.png', dpi=150, bbox_inches='tight')
    plt.close()
if 'hierarchical' in cluster_result.results:
    fig, ax = plt.subplots(figsize=(10, 6))
    linkage_matrix = cluster_result.results['hierarchical']['linkage']
    dendrogram(linkage_matrix, ax=ax, no_labels=True)
    plt.savefig('/tmp/cluster_hierarchical.png', dpi=150, bbox_inches='tight')
    plt.close()
print(f"  Plotting (all 3): {time.time() - start:.2f}s")

print(f"\nTotal time: {time.time() - total_start:.2f}s")
