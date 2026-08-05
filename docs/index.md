# FastMDXplora documentation

> **F**ully **A**utomated **Sy**s**T**em for **M**olecular **D**ynamics e**Xplora**tion

FastMDXplora explores a protein's behavior end to end from a single command.
Given a structure (or just a PDB ID) it prepares the system, runs OpenMM
molecular dynamics, analyzes the trajectory, and writes publication-ready
reports.

## Start here

New to FastMDXplora? Follow the [Beginner's guide](getting_started.md). It
walks through installation, backend checks, a short smoke test, a complete
CLI run, the live dashboard, the Python API, protein-ligand/PLUMED runs, and
the `v1` compatibility profile.

The shortest safe first run is:

```bash
fastmdx explore --system 1L2Y --output runs/first_smoke \
  --include setup simulation --simulate-preset gentle --simulate-platform CPU
```

Then use the same output directory with `--dashboard` for live monitoring or
`fastmdx gui --output ...` to reopen an existing run. For long
GPU jobs, read [Production and GPU runs](production.md) first.

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
region_highlights
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
interactions_design
pdb_smoke_campaign
```
