<div align="center">

# FastMDXplora

**Molecular dynamics from a PDB code to a finished study — in one command.**

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fjcc.70350-blue)](https://doi.org/10.1002/jcc.70350)
[![PyPI](https://img.shields.io/pypi/v/fastmdxplora?label=pypi)](https://pypi.org/project/fastmdxplora/)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/fastmdxplora?label=conda-forge&color=44A833)](https://anaconda.org/conda-forge/fastmdxplora)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/fastmdxplora?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/fastmdxplora)
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

```
setup  →  simulation  →  analysis  →  report
```

Four characters of PDB ID as input. FastMDXplora fetches T4 lysozyme,
parameterises the benzene bound in its cavity, runs the dynamics, analyses the
trajectory, works out which residues hold the ligand in place, and writes the
whole study up as a PDF.

Run all four phases, or any one on its own — `fastmdx setup`, `simulate`,
`analyze`, `report`. Each records what it did, so a run can be picked up,
repeated or explained afterwards.

Or do the whole thing in the GUI:

```bash
fastmdx gui
```

The GUI is not a viewer bolted onto a command-line tool. Its form is generated
from the same schema the CLI validates against, so every one of the hundred
settings is reachable, every system FastMDXplora can study can be built there,
and nothing can be configured in one interface and not the other. What it hands
back is a config file — checked by the same validator before you ever see it —
and `fastmdx explore --config` runs those same bytes on a laptop or a cluster.
Build the study where it is convenient to think; run it where the compute is.

## Install

```bash
conda create -n fastmdxplora -c conda-forge fastmdxplora
conda activate fastmdxplora
```

`fastmdx info` lists every backend and how to get anything missing:

```bash
fastmdx info
```

## What you can study

| | |
|---|---|
| **A protein on its own** | Fold, flexibility, secondary structure, native contacts, conformational clustering — from a PDB code. |
| **A protein with a ligand** | The ligand is found, its chemistry resolved, its protonation settled in the binding site. Eight interaction types against published criteria tell you what *holds* it, not just what it touches. |
| **A membrane protein** | Embedded in one of seven bilayers, with the orientation checked rather than assumed and pressure coupling that suits a lipid system. |
| **Free energy along a coordinate** | Umbrella sampling, metadynamics and steered MD from a named collective variable — eight of them — without writing PLUMED input. Each says what its output is and is not: a surface if the bias converged, a pathway and the work along it, a potential of mean force if the windows overlap. |
| **A trajectory from another engine** | Skip the simulation and analyse what you already have — GROMACS `.xtc` and `.trr`, Amber `.nc`, NAMD and CHARMM `.dcd`, LAMMPS `.lammpstrj`, or `.pdb`, `.cif` and `.h5` that carry their own topology. |
| **Many systems at once** | Mutants against wild type, a sweep across a setting, runs pinned one per GPU, and a comparison report across all of them. |

**It refuses rather than guesses.** An ambiguous ligand charge, or a protein
pointed the wrong way into a membrane, stops the run and gets named. A
metadynamics run that crossed its barrier once when you asked for four gets no
free energy surface. Averages from a biased trajectory are corrected back to
equilibrium where the bias allows it, and labelled where it does not.

Every step says why it is happening, and cites the paper worth reading. What
comes out is something you can defend — or is marked clearly as something you
cannot.

## The config is the study

A FastMDXplora config is the whole description of a molecular dynamics study:
the system, how it is prepared, how it is simulated, what is measured, and how
it is written up. Capture that, and the four phases — setup, simulation,
analysis, report — run themselves.

```yaml
systems:
  - system: 181L
simulation:
  duration_ns: 100
```

That is a complete study. Everything unnamed takes a documented default, and
every run writes `resolved_config.yml` with defaults, file and command line
merged, so the exact study can be run again by anyone holding that one file.

**Three ways to build a config**, and each of them also runs all four phases:

| | |
|---|---|
| **The GUI** | `fastmdx gui`. A form generated from the schema, so every system and every setting is reachable. Worth using even for a command-line or Python workflow: build the study where the options are visible and explained, then take the file away. |
| **The CLI** | `fastmdx explore --config study.yml`, or `fastmdx init-config` for a commented template, or a flag for any setting. |
| **The Python API** | `FastMDXplora(config="study.yml").explore()`, or the same blocks passed as options. |

A fourth way is to write the YAML by hand, which is short and readable enough
that people do.

None of these is the primary interface and none is a subset of another. The
form, the flags and the API are generated from one declaration of what the
software can do, so a study you can express in one you can express in all of
them. A study designed on a laptop in the GUI runs unchanged on a cluster from
the command line, because what travels between them is the config.

## Documentation

**Start here** — [Install](https://fastmdxplora.readthedocs.io/en/latest/installation.html) ·
[Your first run](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html) ·
[The GUI](https://fastmdxplora.readthedocs.io/en/latest/gui.html) ·
[The four phases](https://fastmdxplora.readthedocs.io/en/latest/phases.html)

**Going further** — [Restraints, membranes, enhanced sampling](https://fastmdxplora.readthedocs.io/en/latest/simulations.html) ·
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
