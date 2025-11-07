import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path('.') / 'src'))
from fastmdanalysis import FastMDAnalysis
from fastmdanalysis.datasets import TrpCage
import warnings
warnings.filterwarnings('ignore')

# Test raw computation without plotting
start = time.time()
fastmda = FastMDAnalysis(TrpCage.traj, TrpCage.top, frames=(0, -1, 10), atoms="protein")

# Run analyses without plotting
rmsd_result = fastmda.rmsd(ref=0)
rmsf_result = fastmda.rmsf()
rg_result = fastmda.rg()
cluster_result = fastmda.cluster(methods='kmeans', n_clusters=3)

end = time.time()
print(f"FastMDA raw computation time: {end - start:.2f}s")
