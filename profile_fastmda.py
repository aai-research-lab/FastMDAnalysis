import sys, time
from pathlib import Path
sys.path.insert(0, str(Path('.') / 'src'))
from fastmdanalysis import FastMDAnalysis
from fastmdanalysis.datasets import Ubiquitin
import warnings
warnings.filterwarnings('ignore')

print("Profiling FastMDA individual analyses:\n")

# Initialize
start = time.time()
fastmda = FastMDAnalysis(Ubiquitin.traj, Ubiquitin.top, frames=(0, -1, 10), atoms="protein")
print(f"Initialization: {time.time() - start:.3f}s")

# RMSD
start = time.time()
rmsd_result = fastmda.rmsd(ref=0)
print(f"RMSD: {time.time() - start:.3f}s")

# RMSF  
start = time.time()
rmsf_result = fastmda.rmsf()
print(f"RMSF: {time.time() - start:.3f}s")

# RG
start = time.time()
rg_result = fastmda.rg()
print(f"RG: {time.time() - start:.3f}s")

# Cluster - all methods
start = time.time()
cluster_result = fastmda.cluster(
    methods=['kmeans', 'dbscan', 'hierarchical'],
    n_clusters=3,
    eps=0.5,
    min_samples=2
)
print(f"Cluster (all 3): {time.time() - start:.3f}s")

# Test individual cluster methods
print("\nIndividual cluster methods:")
start = time.time()
km = fastmda.cluster(methods=['kmeans'], n_clusters=3)
print(f"  KMeans only: {time.time() - start:.3f}s")

start = time.time()
db = fastmda.cluster(methods=['dbscan'], eps=0.5, min_samples=2)
print(f"  DBSCAN only: {time.time() - start:.3f}s")

start = time.time()
hier = fastmda.cluster(methods=['hierarchical'], n_clusters=3)
print(f"  Hierarchical only: {time.time() - start:.3f}s")
