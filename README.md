<div align="center">

# FastMDXplora

**Molecular dynamics from a PDB code to a finished study — in one command.**

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fjcc.70350-blue)](https://doi.org/10.1002/jcc.70350)
[![PyPI](https://img.shields.io/pypi/v/fastmdxplora?label=pypi)](https://pypi.org/project/fastmdxplora/)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/fastmdxplora?label=conda-forge&color=44A833)](https://anaconda.org/conda-forge/fastmdxplora)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://pypi.org/project/fastmdxplora/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Tests](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml/badge.svg)](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/aai-research-lab/FastMDXplora/branch/main/graph/badge.svg)](https://codecov.io/gh/aai-research-lab/FastMDXplora)
[![Docs](https://img.shields.io/readthedocs/fastmdxplora?label=docs)](https://fastmdxplora.readthedocs.io)
[![conda downloads](https://img.shields.io/conda/dn/conda-forge/fastmdxplora?label=conda%20downloads&color=44A833)](https://anaconda.org/conda-forge/fastmdxplora)
[![OpenMM](https://img.shields.io/badge/engine-OpenMM-orange)](https://openmm.org)

[**Documentation**](https://fastmdxplora.readthedocs.io) ·
[**Quick start**](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html) ·
[**GUI**](https://fastmdxplora.readthedocs.io/en/latest/gui.html) ·
[**Cite**](#citation)

</div>

---

```bash
fastmdx explore --system 181L
```

Four characters of input. FastMDXplora fetches T4 lysozyme, parameterises the
benzene bound in its cavity, runs the dynamics, measures the trajectory
fifteen ways, works out which residues hold the ligand in place, and writes
the whole thing up as a PDF you could put in front of a supervisor.

Or open the GUI and watch it happen:

```bash
fastmdx gui
```

## Install

```bash
conda install -c conda-forge fastmdxplora
```

Everything included: OpenMM, ligand parameterisation, enhanced sampling, PDF
reports. `pip install fastmdxplora` works too, but two dependencies are
conda-only — `fastmdx info` tells you which and how to get them.

## Why FastMDXplora

**Stops instead of guessing.** An ambiguous ligand charge, a protein pointed
the wrong way for a membrane, a metadynamics run that could never converge —
FastMDXplora halts and names what it could not decide. You never publish a
number it invented.

**Tells you how much the run supports.** Every report states how many
*independent* observations the trajectory holds, which is far fewer than the
frame count. Interaction occupancies do the same: 450 consecutive frames and
450 alternating frames are both "50%", and FastMDXplora shows you which one
you have.

**Writes the methods paragraph for you.** Assembled from what the run actually
recorded, against the reporting checklists journals apply. Nothing invented,
and anything the run did not record is named as missing.

**A GUI that leaves nothing out.** Design a run, start it, watch the molecule
move, read the results — with every setting the command line has, because both
are generated from one declaration.

## What it does

| | |
|---|---|
| **Setup** | A PDB code or a file becomes a simulation-ready system: repaired, protonated, solvated, parameterised. Ligands are identified and their chemistry resolved; membranes are built and the protein's orientation checked. |
| **Simulation** | OpenMM, with restraints released in stages, a barostat that knows a bilayer when it sees one, and metadynamics you set up by naming a coordinate. |
| **Analysis** | Fifteen measures, including eight types of protein-ligand interaction, each against a published criterion. |
| **Report** | Markdown, PDF, slides, a browser dashboard, and a bundle you can send to a collaborator. |

## Documentation

**Start here** — [Install](https://fastmdxplora.readthedocs.io/en/latest/installation.html) ·
[Your first run](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html) ·
[The GUI](https://fastmdxplora.readthedocs.io/en/latest/gui.html) ·
[The four phases](https://fastmdxplora.readthedocs.io/en/latest/phases.html)

**Going further** — [Restraints, membranes, metadynamics](https://fastmdxplora.readthedocs.io/en/latest/simulations.html) ·
[Production and GPUs](https://fastmdxplora.readthedocs.io/en/latest/production.html) ·
[Protein-ligand interactions](https://fastmdxplora.readthedocs.io/en/latest/interactions.html)

**Reference** — [CLI](https://fastmdxplora.readthedocs.io/en/latest/cli_reference.html) ·
[Configuration](https://fastmdxplora.readthedocs.io/en/latest/configuration.html) ·
[Examples](https://fastmdxplora.readthedocs.io/en/latest/usage_examples.html) ·
[Python API](https://fastmdxplora.readthedocs.io/en/latest/api.html)

## Citation

> Aina, A.; Kwan, D. *FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics Trajectories.* J. Comput. Chem. **2026**, 47, e70350. DOI: [10.1002/jcc.70350](https://doi.org/10.1002/jcc.70350)

```bibtex
@article{aina2026fastmd,
  author  = {Aina, Adekunle and Kwan, Derrick},
  title   = {FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics Trajectories},
  journal = {Journal of Computational Chemistry},
  volume  = {47},
  number  = {8},
  pages   = {e70350},
  year    = {2026},
  doi     = {10.1002/jcc.70350},
}
```

## Contributing

Contributions welcome — see [CONTRIBUTING.md](CONTRIBUTING.md). FastMDXplora
follows the [Contributor Covenant](CODE_OF_CONDUCT.md).

## License

MIT. See [LICENSE](LICENSE).

---

<div align="center">

Built in the [AAI Research Lab](https://aai-research-lab.github.io) at
California State University Dominguez Hills, on MDTraj, OpenMM, PDBFixer,
OpenFF, RDKit, NumPy, SciPy, scikit-learn and Matplotlib.

</div>
