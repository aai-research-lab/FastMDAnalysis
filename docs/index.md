# FastMDXplora documentation

> **F**ully **A**utomated **Sy**s**T**em for **M**olecular **D**ynamics e**Xplora**tion

FastMDXplora explores a protein's behavior end to end from a single command.
Given a structure (or just a PDB ID) it prepares the system, runs OpenMM
molecular dynamics, analyzes the trajectory, and writes publication-ready
reports.

## Start here

If you are new to FastMDXplora, follow these steps:

1. [Install the package and chemistry backends](installation.md).
2. Run `fastmdx info` to confirm which backends are available.
3. Run a conservative CPU smoke test:

   ```bash
   fastmdx explore --system 1L2Y --include setup simulation \
     --simulate-preset gentle --simulate-platform CPU
   ```

4. When that succeeds, run the complete workflow:

   ```bash
   fastmdx explore --system 1L2Y --dashboard
   ```

The dashboard opens at a local URL printed in the terminal. For a command
reference, see [CLI commands](cli_reference.md). For an explanation of every
output directory, see [The four phases](phases.md).

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
phases
cli_reference
usage_examples
configuration
region_highlights
live_dashboard
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
