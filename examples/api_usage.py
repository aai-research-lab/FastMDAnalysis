"""
Example API usage for FastMDAnalysis.

This example demonstrates the use of FastMDAnalysis with the new API design. 
The FastMDAnalysis object is created once by providing the trajectory and topology file paths,
as well as optional frame selection (with support for negative indices) and atom selection.
Subsequent analysis methods (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, and dimred) use the preloaded
trajectory and default atom selection unless overridden.

Usage:
  - Update the file paths below to point to your trajectory and topology files.
  - Run this script to perform various MD analyses.
"""

from FastMDAnalysis import FastMDAnalysis

def main():
    # Specify the trajectory and topology file paths.
    traj_path = "path/to/trajectory.dcd"
    top_path  = "path/to/topology.pdb"
    
    # Instantiate FastMDAnalysis.
    # 'frames' is specified as (start, stop, stride). Negative indices are allowed;
    # using -1 for stop means the last frame is included.
    # 'atoms' specifies the default atom selection; here we use "protein".
    fastmda = FastMDAnalysis(traj_path, top_path, frames=(0, -1, 10), atoms="protein")
    
    # ---------------------------
    # Run RMSD Analysis
    # ---------------------------
    rmsd_analysis = fastmda.rmsd(ref=0)
    print("RMSD Data:", rmsd_analysis.data)
    rmsd_plot_file = rmsd_analysis.plot()
    print("RMSD plot saved to:", rmsd_plot_file)
    
    # ---------------------------
    # Run RMSF Analysis
    # ---------------------------
    rmsf_analysis = fastmda.rmsf()
    print("RMSF Data:", rmsf_analysis.data)
    rmsf_plot_file = rmsf_analysis.plot()
    print("RMSF plot saved to:", rmsf_plot_file)
    
    # ---------------------------
    # Run Radius of Gyration (RG) Analysis
    # ---------------------------
    rg_analysis = fastmda.rg()
    print("RG Data:", rg_analysis.data)
    rg_plot_file = rg_analysis.plot()
    print("RG plot saved to:", rg_plot_file)
    
    # ---------------------------
    # Run Hydrogen Bonds (HBonds) Analysis
    # ---------------------------
    hbonds_analysis = fastmda.hbonds()
    print("Hydrogen Bonds Data:", hbonds_analysis.data)
    hbonds_plot_file = hbonds_analysis.plot()
    print("Hydrogen Bonds plot saved to:", hbonds_plot_file)
    
    # ---------------------------
    # Run Clustering Analysis (using DBSCAN)
    # ---------------------------
    cluster_analysis = fastmda.cluster(methods="dbscan", eps=0.5, min_samples=5)
    print("Clustering Results (DBSCAN):", cluster_analysis.results)
    cluster_plot_file = cluster_analysis.plot(method="dbscan")
    print("Cluster population plot (DBSCAN) saved to:", cluster_plot_file)
    
    # ---------------------------
    # Run Secondary Structure (SS) Analysis
    # ---------------------------
    ss_analysis = fastmda.ss()
    print("Secondary Structure Data:", ss_analysis.data)
    ss_plot_file = ss_analysis.plot()
    print("SS plot saved to:", ss_plot_file)
    
    # ---------------------------
    # Run SASA Analysis
    # ---------------------------
    sasa_analysis = fastmda.sasa(probe_radius=0.14)
    print("SASA Data:", sasa_analysis.data)
    sasa_plots = sasa_analysis.plot(option="all")
    print("SASA plots saved to:", sasa_plots)
    
    # ---------------------------
    # Run Dimensionality Reduction Analysis (using PCA and t-SNE)
    # ---------------------------
    dimred_analysis = fastmda.dimred(methods=["pca", "tsne"], atom_selection="protein and name CA")
    print("Dimensionality Reduction Data:", dimred_analysis.data)
    dimred_plots = dimred_analysis.plot()
    print("Dimensionality Reduction plots saved to:", dimred_plots)

if __name__ == "__main__":
    main()

