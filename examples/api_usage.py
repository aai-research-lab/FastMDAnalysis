"""
Example API usage for FastMDAnalysis.

This example demonstrates how to run an RMSD analysis using the FastMDAnalysis API.
The analysis calculates RMSD values relative to a reference frame and automatically
saves the results and a default plot. Accessing the analysis data and replotting is optional;
you may use these features to inspect or customize the output further.
"""

from FastMDAnalysis import FastMDAnalysis

# Initialize the FastMDAnalysis object.
fastmda = FastMDAnalysis()

# Run RMSD analysis. Update the file paths to your trajectory and topology files.
rmsd_analysis = fastmda.rmsd("path/to/trajectory.dcd", top="path/to/topology.pdb", 
                             output="rmsd_output", ref=0)

# (Optional) Retrieve the RMSD data.
#rmsd_data = rmsd_analysis.data
#print("RMSD Data:", rmsd_data)

# (Optional) Replot the RMSD data with default settings (or pass in custom matplotlib options).
#rmsd_plot_file = rmsd_analysis.plot()
#print("RMSD plot saved to:", rmsd_plot_file)

