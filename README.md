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

```
  setup  →  simulation  →  analysis  →  report
```

## Highlights

- Explore a protein's full dynamics with a single command, covering setup, simulation, analysis, and reporting
- Build a protein-ligand system from a PDB identifier alone: the ligand is identified, its chemistry retrieved, and its protonation settled in the binding site, with a refusal rather than a guess where the structure is ambiguous
- Put a membrane protein in a lipid bilayer from a PDB identifier in one command, with the orientation checked rather than assumed and the barostat chosen from the topology
- Hold parts of a system still while the solvent settles around them, released in stages and reported in the methods
- Bias a run along a coordinate you name rather than PLUMED input you write
- Find out what holds a ligand in place, not just how much of the protein it touches: hydrogen bonds, hydrophobic contacts, salt bridges, pi-stacking, pi-cation, halogen bonds, metal coordination and water bridges, each against a published criterion
- Read an occupancy that says how much watching it rests on: a contact present in 450 consecutive frames and one present in 450 alternating frames are both fifty per cent, and only the second has an error bar
- Reach beyond plain MD with built-in PLUMED enhanced sampling (metadynamics, umbrella sampling, steered MD)
- Design, start, watch, and review an exploration from a browser, with a 3D viewer and live telemetry
- Analyse a trajectory you already have, from any engine, without running the simulation again
- Build a config in the browser and run it there, or take it to a cluster — the same file either way
- Scale from a quick single-protein exploration to large parallel campaigns, driven the same way from the CLI or the Python API

## Phases of FastMDXplora

| Phase | What it does |
|---|---|
| **setup** | Cleans up your structure and builds a simulation-ready system: fixes missing atoms, adds hydrogens, solvates, and adds ions. Decides what each non-standard residue means — a bound ligand is parameterized, a cryoprotectant discarded, a coordinated metal kept — and stops where the structure does not say. |
| **simulation** | Runs the molecular dynamics (energy minimization, equilibration, and production), with optional enhanced sampling. |
| **analysis** | Computes the standard structural and dynamic metrics (and protein-ligand metrics when a ligand is present), with figures ready to use. |
| **report** | Packages everything into a written report (Markdown and PDF), a slide deck, a browser dashboard, and a self-contained bundle you can share. The report includes a methods paragraph written against the reporting checklists journals apply, and a convergence assessment that says what the run can and cannot support. |

## Install

From conda-forge, which brings every backend with it:

```bash
conda create -n fastmdxplora -c conda-forge fastmdxplora
conda activate fastmdxplora
```

This is the recommended route. OpenMM, PDBFixer, OpenFF, RDKit, and PROPKA are
conda-forge packages, and the ligand path needs all of them; installing them
any other way is more work for the same result.

From PyPI, if you already manage those yourself or only need part of the
pipeline:

```bash
pip install fastmdxplora            # analysis and reporting
pip install "fastmdxplora[md]"      # adds setup and simulation
pip install "fastmdxplora[ligand]"  # adds most of protein-ligand preparation
pip install "fastmdxplora[pdf]"     # adds the PDF report
```

Two of those are partial from PyPI, which is why conda-forge is the
recommendation rather than a preference.

The **ligand** path also needs the OpenFF toolkit, which has no PyPI
distribution at all — it is a conda-forge package, so no pip command reaches
it:

```bash
conda install -c conda-forge openff-toolkit openmmforcefields
```

The **pdf** extra needs Pango and Cairo present as system libraries. The
conda-forge package brings its own; from PyPI they have to be installed
separately. Where they are missing the run says so and writes the other
formats rather than failing.

Run `fastmdx info` afterwards to see which backends were found.

Working on FastMDXplora itself, Windows and WSL2, and troubleshooting a partial
install are covered in the
[installation guide](https://fastmdxplora.readthedocs.io/en/latest/installation.html).

## Quick start

```bash
fastmdx explore --system 1L2Y      # a full exploration, from a PDB ID
fastmdx gui                        # design, start, and watch one in a browser
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

Parameter sweeps, multi-system campaigns, cross-run comparison, parallel
execution, and the full flag list are covered in the
[usage examples](https://fastmdxplora.readthedocs.io/en/latest/usage_examples.html)
and the
[CLI reference](https://fastmdxplora.readthedocs.io/en/latest/cli_reference.html).

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

Full documentation is at
[fastmdxplora.readthedocs.io](https://fastmdxplora.readthedocs.io):

- [Beginner's guide](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html): first run, start to finish
- [Protein-ligand interactions](https://fastmdxplora.readthedocs.io/en/latest/interactions.html): what holds a ligand, and how confidently it is known
- [Simulations beyond a box of water](https://fastmdxplora.readthedocs.io/en/latest/simulations.html): restraints, membranes, metadynamics, and what to do when a run fails
- [Installation](https://fastmdxplora.readthedocs.io/en/latest/installation.html): every route, platform notes, troubleshooting
- [GUI](https://fastmdxplora.readthedocs.io/en/latest/gui.html): the browser interface, and watching a cluster run
- [Configuration](https://fastmdxplora.readthedocs.io/en/latest/configuration.html): the YAML file, option by option
- [CLI reference](https://fastmdxplora.readthedocs.io/en/latest/cli_reference.html): every command and flag
- [Production runs](https://fastmdxplora.readthedocs.io/en/latest/production.html): GPUs, long trajectories, scaling

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
