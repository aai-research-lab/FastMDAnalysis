API Usage
=========

The Python API mirrors the architecture described in the FastMDAnalysis paper:
instantiate one controller object, reuse the cached trajectory, and compose
analyses that return self-contained result objects.

1. Import and initialisation
----------------------------

.. code-block:: python

	from fastmdanalysis import FastMDAnalysis
	from fastmdanalysis.datasets import trp_cage

	fastmda = FastMDAnalysis(
		 trp_cage.traj,
		 trp_cage.top,
		 frames=(0, -1, 5),       
		 atoms="protein and name CA"  
	)

``FastMDAnalysis`` immediately loads the trajectory using
``utils.load_trajectory``. Lists, comma-separated strings, and glob patterns are
accepted; negative indices in ``frames`` are resolved to absolute positions, as
highlighted in the manuscript.

2. Running core analyses
------------------------

Each method constructs a dedicated analysis object, invokes ``run()``, and
returns the populated instance. Example:

.. code-block:: python

	rmsd = fastmda.rmsd(reference_frame=0)
	rg = fastmda.rg()
	hb = fastmda.hbonds()

All three calls create ``<analysis>_output`` directories, save ``.dat`` files,
write PNG plots, and drop a log file named after the command (e.g.,
``rmsd_output/rmsd.log``).

3. Inspecting results programmatically
--------------------------------------

Result objects expose three primary hooks:

``analysis.data``
	Numpy array (often columnar) holding the primary metric.

``analysis.results``
	Dictionary containing secondary data or plot paths (e.g., cluster labels,
	dendrogram image).

``analysis.plot(...)``
	Re-generates plots with optional custom styling.

.. code-block:: python

	import numpy as np

	rmsd_nm = rmsd.data[:, 0]
	print(np.max(rmsd_nm))

	custom_path = rmsd.plot(color="firebrick", title="Backbone RMSD")
	print(custom_path)

4. Chaining analyses efficiently
--------------------------------

Because the trajectory is cached, sequential runs avoid reloads:

.. code-block:: python

	analyses = {
		 "rmsd": fastmda.rmsd(reference_frame=0),
		 "rmsf": fastmda.rmsf(atoms="protein"),
		 "sasa": fastmda.sasa(probe_radius=0.14),
	}
	summary = {name: float(a.data.mean()) for name, a in analyses.items()}
	print(summary)

5. Orchestrated runs with ``analyze``
------------------------------------

The ``FastMDAnalysis.analyze`` helper mirrors the CLI orchestrator. It accepts
``include``/``exclude`` lists, per-analysis ``options`` (same schema as the CLI
YAML/JSON file), and a ``slides`` flag or path to emit a PowerPoint deck. Use
``compute_stat=True`` to compute mean/standard deviation summaries for time-series
analyses and overlay the bands on plots.

.. code-block:: python

	result = fastmda.analyze(
		include=["rmsd", "rg"],
		exclude=None,
		options={
			"rmsd": {"ref": 0, "align": True},
			"sasa": {"tick_step_avg": 4, "color_total": "#2c3e50"},
		},
		slides="results.pptx",
		compute_stat=True,
	)

Each entry in ``result`` is an :class:`analysis.AnalysisResult` exposing ``ok``
and ``value`` fields. Figures generated through ``analyze`` inherit the same
publication-focused styling defaults (tick thinning, font balancing,
synchronized colorbars), so CLI and API workflows stay visually consistent.
See :doc:`usage/plotting` for the helper reference and option mapping.

6. Advanced: clustering and embeddings
--------------------------------------

The manuscript dedicates significant space to the clustering workflow. The API
call mirrors that description:

.. code-block:: python

	cluster = fastmda.cluster(
		 methods=["dbscan", "hierarchical"],
		 eps=0.3,
		 min_samples=5,
		 n_clusters=5
	)
	results = cluster.results
	dbscan_labels = results["dbscan"]["labels"]
	dendrogram = results["hierarchical"]["dendrogram_plot"]

Dimensionality reduction leans on scikit-learn as well:

.. code-block:: python

	dim = fastmda.dimred(methods=["pca", "tsne"], atoms="protein")
	embeddings = dict(dim.results)
	pca_coords = embeddings["pca"]["embedding"]

7. Logging, slides, and provenance
---------------------------------

``FastMDAnalysis`` configures Python’s ``logging`` module so each run emits a
timestamped ``.log`` file inside the chosen output directory. Inspect these logs
to see the parameter set, environment versions, and any warnings – an important
reproducibility aid emphasised in the paper.

When ``slides`` is truthy, the orchestrator writes ``fastmda_slides_<ts>.pptx``
alongside the usual analysis artifacts (or honors a custom path). Decks embed
the same publication-ready PNGs saved under ``<analysis>_output``.

Next steps
----------

* Dive into :doc:`analysis/index` for per-module formulas and visualisations.
* See :doc:`usage/cli` if you prefer scripted batch processing.
* Review :doc:`usage/plotting` for deeper control over ticks, fonts, and
	colorbars exposed by ``options`` dictionaries.
* Consult :doc:`contributing` when you are ready to add a new analysis class.
