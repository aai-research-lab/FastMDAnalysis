"""
Example API usage for FastMDAnalysis.
"""

from FastMDAnalysis import FastMDAnalysis

# Initialize the FastMDAnalysis object.
fastmda = FastMDAnalysis()

# Run RMSD analysis (update the file paths accordingly).
rmsd_analysis = fastmda.rmsd("path/to/trajectory.dcd", top="path/to/topology.pdb", output="rmsd_output", ref=0)
# Get RMSD data.
data = rmsd_analysis.data
print("RMSD Data:", data)
# Plot and save RMSD.
plot_file = rmsd_analysis.plot()
print("RMSD plot saved to:", plot_file)

