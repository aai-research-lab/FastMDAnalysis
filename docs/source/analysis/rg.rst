.. currentmodule:: FastMDAnalysis.analysis.rg

Radius of Gyration
==================

The :class:`RGAnalysis` module reports the radius of gyration (Rg) trajectory
described in Section 3.3 of the manuscript. Rg tracks the mass-weighted spatial
extent of the selected atoms, making it a convenient scalar descriptor of global
compaction or expansion events.

Theory
------

For :math:`N` atoms with coordinates :math:`\mathbf{r}_i(t)` and a center of
geometry :math:`\mathbf{r}_{\mathrm{COM}}(t)`, the instantaneous radius of
gyration is

.. math::

   R_g(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N \left\| \mathbf{r}_i(t) - \mathbf{r}_{\mathrm{COM}}(t) \right\|^2 }

The MDTraj implementation used here is equivalent to the gyration tensor trace
and operates in nanometres.

Implementation Highlights
-------------------------

- Delegates to :func:`mdtraj.compute_rg`, respecting any atom selection supplied
  either globally (``FastMDAnalysis(..., atoms=...)``) or locally via
  ``fastmda.rg(atoms=...)``.
- Results are stored as a column vector in ``self.data`` and persisted via
  :meth:`~FastMDAnalysis.analysis.base.BaseAnalysis._save_data`.
- :meth:`RGAnalysis.plot` renders an Rg-versus-frame trace with optional custom
  styling (colour, line style, labels).

Usage
-----

**API**

.. code-block:: python

   from fastmdanalysis import FastMDAnalysis

   fastmda = FastMDAnalysis("traj.dcd", "top.pdb", atoms="protein")
   rg = fastmda.rg()
   rg_results = rg.run()["rg"]
   rg.plot(color="#9467bd", title="Protein Radius of Gyration")

**CLI**

.. code-block:: bash

   fastmda rg -traj traj.dcd -top top.pdb --atoms "protein" -o analysis/rg

Outputs
-------

- ``rg_output/rg.dat`` — one value per frame (nm).
- ``rg_output/rg_stats.dat`` — mean and standard deviation when
  ``compute_stat=True``.
- ``rg_output/rg.png`` — default line plot written by :meth:`RGAnalysis.plot`.
- ``rg_output/rg.log`` — execution log captured by the CLI wrapper.

Interpretation Tips
-------------------

- Sudden drops point to folding or collapse events; plateaus imply conformational stability.
- Compare Rg trends with :doc:`sasa` to distinguish between compaction and solvent exposure changes.

Troubleshooting
---------------

- ``AnalysisError: No atoms selected`` indicates the atom selection did not
  match the topology; verify residue naming, chain IDs, and capitalization.
- For multi-chain systems ensure the topology contains all chains referenced by
  the selection string; MDTraj silently excludes missing chains which can bias
  Rg downward.
