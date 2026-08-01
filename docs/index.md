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
:caption: Getting started

installation
getting_started
production
phases
cli_reference
usage_examples
configuration
region_highlights
gui
pdb_smoke_campaign
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
```

## Quick links

- [GitHub repository](https://github.com/aai-research-lab/FastMDXplora)
- [PyPI: fastmdxplora](https://pypi.org/project/fastmdxplora/)
- [PyPI: fastmdx (alias)](https://pypi.org/project/fastmdx/)
- [Foundational paper (JCC 2026)](https://doi.org/10.1002/jcc.70350)

## Citation

If you use FastMDXplora in your work, please cite:

> Aina, A.; Kwan, D. *FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics Trajectories.* J. Comput. Chem. **2026**, 47, e70350. DOI: [10.1002/jcc.70350](https://doi.org/10.1002/jcc.70350)
