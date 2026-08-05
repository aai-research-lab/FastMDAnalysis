# FastMDXplora

> **F**ully **A**utomated **Sy**s**T**em for **M**olecular **D**ynamics e**Xplora**tion

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fjcc.70350-blue)](https://doi.org/10.1002/jcc.70350)
[![PyPI version](https://img.shields.io/pypi/v/fastmdxplora)](https://pypi.org/project/fastmdxplora/)
[![Python versions](https://img.shields.io/pypi/pyversions/fastmdxplora)](https://pypi.org/project/fastmdxplora/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml/badge.svg)](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/aai-research-lab/FastMDXplora/branch/main/graph/badge.svg)](https://codecov.io/gh/aai-research-lab/FastMDXplora)

---

**FastMDXplora** explores a protein's behavior end to end from a single
command. Given a structure (or just a PDB ID) it performs molecular dynamics
exploration all the way through setup, simulation, analysis, and reporting,
then hands back publication-ready results.

## What it does

**One command, a PDB identifier, a finished study.** Setup, simulation,
analysis and a report — with the parts that usually need an expert done for
you, and refused where the structure does not say enough to do them.

- **A GUI that does everything.** Design a run, start it, watch the molecule
  move, read the results — with every setting the CLI has, because both are
  built from the same declaration.
- **Protein-ligand systems from an identifier.** The ligand is found, its
  chemistry retrieved, its protonation settled in the binding site.
- **Membrane proteins too.** Embedded in a bilayer, with the orientation
  checked rather than assumed and the right barostat chosen automatically.
- **What holds a ligand, not just what it touches.** Eight interaction types
  against published criteria, each occupancy carrying the observation behind
  it.
- **Enhanced sampling without PLUMED syntax.** Metadynamics from a named
  coordinate, with walls and funnels.
- **A report you can submit.** Methods written against the checklists journals
  apply, and a convergence assessment that says what the run cannot support.
- **One config, three ways in.** What you build in the GUI runs on a cluster;
  what runs on a cluster opens in the GUI.

## The four phases

```
  setup  →  simulation  →  analysis  →  report
```

**setup** builds a simulation-ready system from your structure — fixing,
solvating, and deciding what each non-standard residue is for.
**simulation** runs it. **analysis** measures it. **report** writes it up.

Each is documented in [The four phases](https://fastmdxplora.readthedocs.io/en/latest/phases.html),
which also says what every measure computes.

## Install

```bash
conda create -n fastmdxplora -c conda-forge fastmdxplora
conda activate fastmdxplora
```

conda-forge is the recommended route because it brings every backend with it.
OpenMM, PDBFixer, the OpenFF toolkit, RDKit, PROPKA and PLUMED are all
conda-forge packages, and two of them are not on PyPI at all — so a pip
install gives you part of the pipeline and tells you which part.

```bash
pip install fastmdxplora            # analysis and reporting
pip install "fastmdxplora[md]"      # adds setup and simulation
```

`fastmdx info` says which backends are present and how to get any that are
not. Extras, Windows and WSL2, and troubleshooting a partial install are in
the [installation guide](https://fastmdxplora.readthedocs.io/en/latest/installation.html).

## Quick start

```bash
fastmdx explore --system 1L2Y      # a full exploration, from a PDB ID
fastmdx gui                        # the FastMDXplora GUI: design, start, and watch one
fastmdx info                       # what is installed and which backends are present
```

From Python:

```python
import fastmdxplora as fastmdx

runs = fastmdx.FastMDXplora(system="1L2Y", output_dir="trpcage").explore()
print(runs[0].output_dir)
```

For anything beyond a quick run, capture the whole exploration in a YAML file.
`fastmdx init-config` writes a commented template, and the same file drives
both the CLI and the Python API:

```bash
fastmdx explore --config study.yml
```


## Outputs by phase

Each phase writes to its own subdirectory under the output root, with a
parameters manifest so every artifact is traceable to the options that
produced it.

| Phase | Key outputs |
|---|---|
| `setup` | `prepared.pdb`, `solvated.pdb`, `setup_parameters.json` |
| `simulation` | `production.dcd`, `topology.pdb`, `simulation_parameters.json` |
| `analysis` | `<analysis>/*.dat`, `<analysis>/*.png`, `analysis_manifest.json` |
| `report` | `report.md`, `report.pdf`, `dashboard.html`, `slides.pptx`, `project_bundle.zip` |

## Documentation

Everything is at
[fastmdxplora.readthedocs.io](https://fastmdxplora.readthedocs.io).

**Start here**

- [Installation](https://fastmdxplora.readthedocs.io/en/latest/installation.html) — every route, and what each extra needs
- [Your first run](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html) — start to finish, with a real structure
- [The FastMDXplora GUI](https://fastmdxplora.readthedocs.io/en/latest/gui.html) — design a run, watch it happen, read the results
- [The four phases](https://fastmdxplora.readthedocs.io/en/latest/phases.html) — what setup, simulation, analysis and report each do, and what every measure means

**Running simulations**

- [Beyond a box of water](https://fastmdxplora.readthedocs.io/en/latest/simulations.html) — restraints, membrane proteins, metadynamics with walls and funnels, and what the phase says when a run fails
- [Production and GPUs](https://fastmdxplora.readthedocs.io/en/latest/production.html) — long trajectories, platforms, scaling

**Reading the results**

- [Protein-ligand interactions](https://fastmdxplora.readthedocs.io/en/latest/interactions.html) — what holds a ligand, how confidently it is known, and why some measurements are refused

**Driving it**

- [CLI reference](https://fastmdxplora.readthedocs.io/en/latest/cli_reference.html) — every command and flag
- [Configuration](https://fastmdxplora.readthedocs.io/en/latest/configuration.html) — the YAML file, option by option
- [Worked examples](https://fastmdxplora.readthedocs.io/en/latest/usage_examples.html) — recipes for common studies

## Citation

If you use FastMDXplora in your work, please cite:

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

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md). FastMDXplora follows the [Contributor Covenant](CODE_OF_CONDUCT.md).

## License

MIT. See [LICENSE](LICENSE).

## Acknowledgements

FastMDXplora is developed in the [AAI Research Lab](https://aai-research-lab.github.io) at California State University Dominguez Hills. It builds on a deep ecosystem of open-source scientific Python: MDTraj, OpenMM, PDBFixer, NumPy, SciPy, scikit-learn, Matplotlib, python-pptx, and many others.
