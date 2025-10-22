.. currentmodule:: FastMDAnalysis.analysis.dimred

Dimensionality Reduction
========================

The :class:`DimRedAnalysis` module projects high-dimensional coordinates onto
two dimensions using PCA, classical MDS, or t-SNE as detailed in Section 3.8 of
the manuscript. Colour-coded scatter plots reveal slow collective motions and
metastable states.

Feature Matrix
--------------

- By default the analysis selects ``"protein and name CA"`` atoms and flattens
  the resulting coordinates into a ``(n_frames, n_atoms * 3)`` matrix.
- Supply ``atoms`` to override the default (for example, include ligand
  heavy atoms when probing binding pathways).

Available Methods
-----------------

- ``pca`` — Principal component analysis via :class:`sklearn.decomposition.PCA`.
- ``mds`` — Metric MDS (:class:`sklearn.manifold.MDS`) using Euclidean distances.
- ``tsne`` — t-distributed stochastic neighbour embedding.
- ``all`` — Convenience option that runs PCA, MDS, and t-SNE sequentially.

Implementation Highlights
-------------------------

- Accepts a string or list for ``methods``; mixed requests (for example,
  ``["pca", "tsne"]``) are supported.
- Stores each embedding as a two-column array saved to
  ``dimred_<method>.dat`` and caches plots under the same stem.
- :meth:`DimRedAnalysis.plot` colours points by frame index so recent frames
  appear warm-toned and early frames cool-toned, matching the manuscript
  figures.

Usage
-----

**API**

.. code-block:: python

   from FastMDAnalysis import FastMDAnalysis

   fastmda = FastMDAnalysis("traj.dcd", "top.pdb")
  dimred = fastmda.dimred(methods=["pca", "tsne"], atoms="protein and resname LIG")
   embeddings = dimred.run()
   dimred.plot(method="pca", cmap="viridis")

**CLI**

.. code-block:: bash

  fastmda dimred -traj traj.dcd -top top.pdb --methods pca tsne \
    --atoms "protein and resname LIG" -o analysis/dimred

Outputs
-------

- ``dimred_output/dimred_<method>.dat`` — two-component embedding for each requested method.
- ``dimred_output/dimred_<method>.png`` — scatter plot coloured by frame index.
- ``dimred_output/dimred.log`` — consolidated log file.

Interpretation Tips
-------------------

- PCA captures variance linearly; inspect explained variance ratios (available
  via ``embeddings["pca"]`` and the parent ``PCA`` estimator) to ensure two
  components suffice.
- t-SNE magnifies local neighbourhoods. Align the scatter with clustering
  results by colouring points using ``cluster.results[...]`` for richer context.

Troubleshooting
---------------

- ``AnalysisError: ... atom selection`` indicates the selection returned zero
  atoms; verify residue naming and chain IDs.
- t-SNE is stochastic; rerun with ``random_state`` overrides by subclassing or
  edit ``DimRedAnalysis`` if you require reproducibility across sessions.
- t-SNE and MDS scale poorly with frame count; down-sample long trajectories or
  pre-filter frames via clustering before running dimensionality reduction.
