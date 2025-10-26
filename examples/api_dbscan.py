from fastmdanalysis import FastMDAnalysis
from fastmdanalysis.datasets import TrpCage

fmda = FastMDAnalysis(TrpCage.traj, TrpCage.top, frames=(0, -1, 10), atoms="protein")
a = fmda.cluster(methods="dbscan", eps=0.17, min_samples=5)  # 0.2 nm ≈ 2 Å
print(a.results["dbscan"]["n_clusters"], a.results["dbscan"]["distance_percentiles_nm"])

