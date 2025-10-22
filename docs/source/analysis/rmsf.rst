.. currentmodule:: FastMDAnalysis.analysis.rmsf

RMSF Analysis
=============

The :class:`RMSFAnalysis` module quantifies site-specific flexibility by
calculating the root-mean-square fluctuation (RMSF) of each atom about the
trajectory-average structure. This maps to Section 3.2 of the manuscript and is
useful for identifying mobile loops, termini, or ligand-binding hotspots.

Theory
------

Let :math:`\mathbf{r}_i(t)` denote the Cartesian coordinates of atom
:math:`i` in frame :math:`t`, and let :math:`\bar{\mathbf{r}}_i` be the mean
coordinate of that atom across all frames. The RMSF is

.. math::

   \mathrm{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T} \left\| \mathbf{r}_i(t) - \bar{\mathbf{r}}_i \right\|^2 }

Implementation Highlights
-------------------------

- Uses :func:`mdtraj.rmsf` after constructing an average structure via
  :class:`mdtraj.Trajectory`, reproducing the manuscript pipeline.
- Optional atom selection (for example, ``"protein and name CA"``) restricts
  the calculation to a subset while retaining the original residue numbering in
  output files.
- Automatically writes the column vector of RMSF values and generates a bar
  chart with one bar per atom.

Usage
-----

**API**

.. code-block:: python

   from FastMDAnalysis import FastMDAnalysis

  fastmda = FastMDAnalysis("traj.dcd", "top.pdb")
  rmsf = fastmda.rmsf(atoms="protein and name CA")
  data = rmsf.run()["rmsf"]
  rmsf.plot(title="Cα RMSF")

**CLI**

.. code-block:: bash

  fastmda rmsf -traj traj.dcd -top top.pdb --atoms "protein and name CA" -o analysis/rmsf

Outputs
-------

- ``rmsf_output/rmsf.dat`` — per-atom RMSF values (nm).
- ``rmsf_output/rmsf.png`` — default bar chart with atom index on the x-axis.
- ``rmsf_output/rmsf.log`` — execution log when launched from the CLI.

Interpretation Tips
-------------------

- Plateaus near zero indicate rigid regions; spikes mark flexible segments that
  may correlate with functionally important motions.
- Aggregate RMSF values by residue by grouping three consecutive atoms or by
  replotting after slicing the trajectory to the backbone (``"protein and name CA"``).

Troubleshooting
---------------

- ``AnalysisError: No atoms selected`` signals that the MDTraj selection did
  not match the topology; double-check residue/atom naming and capitalization.
- RMSF magnitudes scale with the frame stride supplied to ``FastMDAnalysis``;
  confirm that equilibration frames are excluded if you see inflated baseline
  noise.
