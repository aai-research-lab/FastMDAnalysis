Analysis Modules
================

FastMDAnalysis implements each metric discussed in the manuscript as a dedicated
class that inherits from :class:`FastMDAnalysis.analysis.base.BaseAnalysis`.
Every module accepts an :class:`mdtraj.Trajectory`, honours optional atom/frame
selections, writes structured outputs, and provides one or more publication-ready
plots. The pages below summarise the scientific background, implementation
details, and usage patterns for each analysis.

.. toctree::
   :maxdepth: 1

   rmsd
   rmsf
   rg
   hbonds
   cluster
   sasa
   ss
   dihedrals
   qvalue
   dimred
