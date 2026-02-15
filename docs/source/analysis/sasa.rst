.. currentmodule:: FastMDAnalysis.analysis.sasa

SASA Analysis
=============

The :class:`SASAAnalysis` module applies the Shrake–Rupley algorithm to compute
solvent accessible surface area (SASA) statistics discussed in Section 3.5 of
the manuscript. It generates three complementary datasets: total SASA per frame,
per-residue SASA tables, and residue-averaged SASA profiles.

Theory
------

Shrake–Rupley places :math:`N_p` probe points on each atom's van der Waals
surface and counts the fraction exposed to solvent. The SASA of residue
:math:`j` in frame :math:`t` is

.. math::

   A_j(t) = \frac{k}{N_p} \times 4 \pi r_j^2

where :math:`k` is the number of probe points not occluded by neighbouring
atoms and :math:`r_j` is the probe-adjusted radius. Summing over residues yields
the total SASA per frame. MDTraj implements this procedure in
:func:`mdtraj.shrake_rupley` using a probe radius of 0.14 nm by default.

Implementation Highlights
-------------------------

- Supports optional atom selections and custom probe radii via the constructor
  parameters ``atoms`` and ``probe_radius``.
- Stores three numpy arrays in ``self.data``: ``total_sasa`` (shape ``(n_frames,)``),
  ``residue_sasa`` (``n_frames x n_residues``), and ``average_residue_sasa``
  (``n_residues,``).
- :meth:`SASAAnalysis.plot` can regenerate any combination of the default plots
  by passing ``option="total"``, ``"residue"``, ``"average"``, or
  ``"all"`` (the default).

Usage
-----

**API**

.. code-block:: python

   from fastmdanalysis import FastMDAnalysis

   fastmda = FastMDAnalysis("traj.dcd", "top.pdb", atoms="protein")
   sasa = fastmda.sasa(probe_radius=0.14)
   datasets = sasa.run()
   sasa.plot(option="all", cmap="magma")

**CLI**

.. code-block:: bash

   fastmda sasa -traj traj.dcd -top top.pdb --atoms "protein" --probe_radius 0.14 -o analysis/sasa

Outputs
-------

- ``sasa_output/total_sasa.dat`` — total SASA per frame (nm²).
- ``sasa_output/total_sasa_stats.dat`` — mean and standard deviation when
  ``compute_stat=True``.
- ``sasa_output/residue_sasa.dat`` — matrix of per-residue SASA values.
- ``sasa_output/average_residue_sasa.dat`` — mean SASA per residue.
- ``sasa_output/total_sasa.png`` — total SASA trace.
- ``sasa_output/residue_sasa.png`` — residue-by-frame heatmap.
- ``sasa_output/average_residue_sasa.png`` — bar chart of average residue SASA.
- ``sasa_output/sasa.log`` — command log when executed via ``fastmda``.

Interpretation Tips
-------------------

- Decreasing total SASA typically signals folding or binding events, whereas
  sustained spikes may indicate partial unfolding or solvent exposure of
  hydrophobic cores.
- The per-residue heatmap is an excellent companion to :doc:`ss` and
  :doc:`rmsf` for detecting persistent solvent hot spots.

Troubleshooting
---------------

- ``AnalysisError: No atoms selected`` means the supplied selection did not
  match any topology atoms; confirm residue names and chain IDs.
- The Shrake–Rupley routine expects all heavy atoms; missing hydrogens are fine
  but truncated side chains can yield underestimates.
- SASA scales with the probe radius: use ``--probe_radius`` to match the probe
  used in your experimental reference before drawing quantitative comparisons.
