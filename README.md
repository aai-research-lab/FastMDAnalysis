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
benzene bound in its cavity, runs the dynamics, analyses the trajectory
fifteen ways, works out which residues hold the ligand in place, and writes
the whole study up as a PDF.

Or open the GUI and watch it happen:

```bash
fastmdx gui
```

## Install

```bash
conda install -c conda-forge fastmdxplora
fastmdx info
```

`fastmdx info` lists every backend and how to get anything missing.

## What you can study

| | |
|---|---|
| **A protein on its own** | Fold, flexibility, secondary structure, native contacts, conformational clustering — from a PDB code. |
| **A protein with a ligand** | The ligand is found, its chemistry resolved, its protonation settled in the binding site. Eight interaction types against published criteria tell you what *holds* it, not just what it touches. |
| **A membrane protein** | Embedded in one of seven bilayers, with the orientation checked rather than assumed and pressure coupling that suits a lipid system. |
| **Free energy along a coordinate** | Metadynamics on a ligand distance, an RMSD, a torsion or a radius of gyration — with walls and funnels, which is what makes an absolute binding free energy recoverable. |
| **A trajectory from another engine** | Skip the simulation and analyse what you already have, from anything MDTraj reads. |
| **Many systems at once** | Mutants against wild type, a sweep across a setting, runs pinned one per GPU, and a comparison report across all of them. |

And where a structure does not say enough — an ambiguous ligand charge, a
protein pointed the wrong way for a membrane — FastMDXplora stops and names
what it could not decide, rather than returning a number that looks fine.

## How a run works

```
  setup  →  simulation  →  analysis  →  report
```

Four phases. Run all of them, or any one on its own — `fastmdx setup`,
`simulate`, `analyze`, `report`. Each records what it did, so a run can be
picked up, repeated or explained afterwards.

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
