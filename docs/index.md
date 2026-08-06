# FastMDXplora

> **F**ully **A**utomated **Sy**s**T**em for **M**olecular **D**ynamics e**Xplora**tion

Give it a structure — or just a PDB identifier — and it prepares the system,
runs the dynamics, measures the trajectory, and writes the study up. The parts
that usually need an expert are done for you, and refused where the structure
does not say enough to do them.

## The quickest way in

```bash
conda create -n fastmdxplora -c conda-forge fastmdxplora
conda activate fastmdxplora
fastmdx gui
```

A tab opens. Type `1L2Y`, press **Run**, and watch a protein fold. Everything
FastMDXplora does is on that page: designing a run, starting it, watching the
molecule move, reading the results.

If you would rather use a terminal:

```bash
fastmdx explore --system 1L2Y --output runs/trpcage
```

[Your first run](getting_started.md) takes either route from nothing to a
finished report in about ten minutes.

## Finding your way around

- **New here?** [Installation](installation.md), then
  [Your first run](getting_started.md).
- **Want to see what it can do?** [The FastMDXplora GUI](gui.md) and
  [The four phases](phases.md). Every step explains itself as it runs, so
  [Your first run](getting_started.md) is also the shortest way to learn what
  the steps are for.
- **Running something real?** [Beyond a box of water](simulations.md) for
  restraints, membranes, and enhanced sampling — umbrella, steered and
  metadynamics; [Production and GPUs](production.md)
  for long runs.
- **Reading results?** [Protein-ligand interactions](interactions.md).
- **Looking for a flag or a setting?** [CLI reference](cli_reference.md),
  [Configuration](configuration.md), or `fastmdx explore --help`, which is
  generated and therefore never out of date.

```{toctree}
:maxdepth: 2
:caption: Start here

installation
getting_started
gui
phases
```

```{toctree}
:maxdepth: 2
:caption: Running simulations

simulations
production
```

```{toctree}
:maxdepth: 2
:caption: Reading the results

interactions
```

```{toctree}
:maxdepth: 2
:caption: Driving it

cli_reference
configuration
usage_examples
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
```

```{toctree}
:maxdepth: 2
:caption: For maintainers

interactions_design
```
