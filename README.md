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

That prepares T4 lysozyme with its ligand, parameterises the benzene, runs the
dynamics, measures fifteen things about the trajectory, works out what holds
the ligand in the pocket, and writes the study up as a PDF — from four
characters of input.

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

## Why this one

**It refuses.** An ambiguous ligand charge, a protein pointed the wrong way
for a membrane, a metadynamics run that would never converge — FastMDXplora
stops and says what it could not decide, instead of returning a number that
looks fine.

**It says how much to trust the answer.** Every report carries a convergence
section stating how many *independent* observations the trajectory holds — far
fewer than the frame count, and usually few enough to matter. Interaction
occupancies carry the same: 450 consecutive frames and 450 alternating frames
are both "50%", and only one has an error bar.

**It writes the methods paragraph.** Against the reporting checklists journals
actually apply, from what the run recorded. Nothing invented; anything missing
named as missing.

**The GUI is the whole tool.** Design a run, start it, watch the molecule move,
read the results — every setting the command line has, because both are built
from one declaration.

## What it does

| | |
|---|---|
| **Setup** | PDB code or file → repaired, protonated, solvated, parameterised. Ligands identified and their chemistry resolved. Membranes built and the orientation checked. |
| **Simulation** | OpenMM. Restraints released in stages, membrane-aware barostat, metadynamics from a named coordinate rather than PLUMED syntax. |
| **Analysis** | Fifteen measures, including eight types of protein-ligand interaction against published criteria. |
| **Report** | Markdown, PDF, slides, a browser dashboard, and a shareable bundle. |

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
