CLI Usage
=========

The ``fastmda`` command-line interface mirrors the API methods described in the
FastMDAnalysis manuscript. It is ideal for batch pipelines, HPC clusters, or
users who prefer not to write Python scripts. Every plot produced through the
CLI runs through the same styling helpers as the API, so PNG outputs are
publication-ready without manual Matplotlib tweaks.

Global options
--------------

``--frames start,stop,stride``
	Apply trajectory subsampling. Negative indices follow Python slicing rules;
	``0,-1,10`` matches the examples in the paper by taking every 10th frame.

``--atoms <selection>``
	Default MDTraj selection string (e.g., ``"protein and name CA"``) that all
	analyses inherit unless a subcommand overrides it.

``--verbose``
	Elevates console logging from ``WARNING`` to ``INFO`` while keeping detailed
	``DEBUG`` output in the ``<analysis>_output/<command>.log`` file.

Subcommands
-----------

Each analysis discussed in the manuscript is available as a subcommand. Common
arguments:

``--trajectory, -traj``
	Path to a trajectory file, list, or glob.

``--topology, -top``
	Path to the corresponding topology file.

``--output, -o``
	Optional explicit output directory; defaults to ``<command>_output``.

``--slides``
	(Orchestrator only) write a PowerPoint deck combining every figure. Accepts a
	boolean flag (timestamped filename) or a path; deck slides embed the same
	publication-ready PNGs saved during each analysis.

Quick-reference table
~~~~~~~~~~~~~~~~~~~~~

================  =========================================  ===============================
Subcommand        Highlights                                   Extra options
================  =========================================  ===============================
``rmsd``          Reference frame alignment, RMSD vs. time     ``--reference-frame`` ``--atoms``
``rmsf``          Per-atom/per-residue flexibility             ``--atoms``
``rg``            Radius of gyration timeseries                ``--atoms``
``hbonds``        Baker–Hubbard hydrogen bonds                 ``--atoms``
``ss``            DSSP-based secondary structure               ``--atoms``
``sasa``          Shrake–Rupley solvent exposure               ``--probe_radius`` ``--atoms``
``cluster``       DBSCAN / KMeans / hierarchical clustering    ``--methods`` ``--n_clusters`` ``--atoms``
``dimred``        PCA / MDS / t-SNE embeddings                 ``--methods`` ``--atoms``
================  =========================================  ===============================

Analyze orchestrator
--------------------

``fastmda analyze`` runs multiple analyses in one pass, sharing the cached
trajectory and styling settings. Use ``--include``/``--exclude`` to control the
portfolio and ``--options`` to point at a YAML/JSON file with per-analysis
keywords (same schema as the Python ``options`` dict). The resulting PNGs and
optional deck preserve the publication-ready defaults.
See :doc:`usage/plotting` for helper details and option descriptions referenced
by each analysis module.

Use ``--compute-stat`` to compute global mean/standard deviation values for
time-series analyses (RMSD, RMSF, Rg, total SASA, Q-value). This adds
``*_stats.dat`` files and overlays mean/±1σ bands on the saved plots.

.. code-block:: bash

	fastmda analyze \
	    --traj data/trp_cage.dcd \
	    --top data/trp_cage.pdb \
	    --include rmsd rg sasa \
	    --options configs/slide_ready.yaml \
	    --slides reports/trp_cage_slides.pptx

Example ``options`` excerpt:

.. code-block:: yaml

	rmsf:
	  tick_step: 3
	  rotate: 45
	sasa:
	  tick_step_avg: 5
	  color_total: "#2c3e50"
	dimred:
	  title_pca: "PCA (Publication)"

Example session
---------------

1. RMSD with logging and custom output directory::

		fastmda rmsd \
			--trajectory data/trp_cage.dcd \
			--topology data/trp_cage.pdb \
			--frames 0,-1,10 \
			--atoms "protein" \
			--reference-frame 0 \
			--output results/trp_rmsd

	Inspect ``results/trp_rmsd/rmsd.dat`` and ``results/trp_rmsd/rmsd.png`` and
	review ``results/trp_rmsd/rmsd.log`` for provenance.

2. Clustering with multiple algorithms, mirroring the manuscript benchmarks::

		fastmda cluster \
			--trajectory data/trp_cage.dcd \
			--topology data/trp_cage.pdb \
			--methods dbscan hierarchical \
			--n_clusters 5 \
			--eps 0.3 \
			--min_samples 5 \
			--verbose

	Output files include population plots, trajectory histograms, scatter plots,
	distance matrices, and dendrograms for each requested algorithm.

3. Dimensionality reduction for rapid visualisation::

	fastmda dimred \
		--trajectory data/trp_cage.dcd \
		--topology data/trp_cage.pdb \
		--methods pca mds tsne \
		--atoms "protein and name CA"

	This yields ``dimred_output/pca_embedding.dat``, ``dimred_output/mds_embedding.dat``,
	and ``dimred_output/tsne_embedding.dat`` alongside matching ``.png`` scatter plots
	for each method (exact filenames depend on the module).

Batch execution tips
--------------------

* Combine commands with shell scripting to reproduce the workflow described in
  the paper (e.g., run RMSD, RMSF, SASA, and clustering sequentially).
* Use ``--output`` to separate experiments, especially when comparing different
  atom selections or frame windows.
* Capture non-zero exit codes in CI/CD pipelines to ensure analyses succeed
  before publishing plots.

For advanced configuration and the analysis APIs behind each subcommand, see
:doc:`analysis/index`, :doc:`usage/api`, and :doc:`usage/plotting`.
