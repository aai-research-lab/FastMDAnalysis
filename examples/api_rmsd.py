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
    traj_path = "protein_traj.dcd"
    top_path  = "protein.pdb"
    
    # Instantiate FastMDAnalysis.
    # 'frames' is specified as (start, stop, stride). Negative indices are allowed;
    # using -1 for stop means the last frame is included.
    # 'atoms' specifies the default atom selection; here we use "protein".
    fastmda = FastMDAnalysis(traj_path, top_path, frames=(0, -1, 10), atoms="protein")
    
    # ---------------------------
    # Run RMSD Analysis
    # ---------------------------
    rmsd_analysis = fastmda.rmsd(ref=0)
    

    # Optinally, retrieve data and/or replot/customized figures
    #print("RMSD Data:", rmsd_analysis.data)
    #rmsd_plot_file = rmsd_analysis.plot()
    #print("RMSD plot saved to:", rmsd_plot_file)
    

if __name__ == "__main__":
    main()

