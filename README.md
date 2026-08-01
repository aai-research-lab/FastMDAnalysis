# FastMDXplora

> **F**ully **A**utomated **Sy**s**T**em for **M**olecular **D**ynamics e**Xplora**tion

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fjcc.70350-blue)](https://doi.org/10.1002/jcc.70350)
[![PyPI version](https://img.shields.io/pypi/v/fastmdxplora)](https://pypi.org/project/fastmdxplora/)
[![Python versions](https://img.shields.io/pypi/pyversions/fastmdxplora)](https://pypi.org/project/fastmdxplora/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml/badge.svg)](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/aai-research-lab/FastMDXplora/branch/main/graph/badge.svg)](https://codecov.io/gh/aai-research-lab/FastMDXplora)

---

**FastMDXplora** explores a protein's behavior end to end from a single command. Given a structure (or just a PDB ID) it performs molecular dynamics exploration all the way through setup, simulation, analysis, and reporting, then hands back publication-ready results. Start with the [Beginner's guide](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html) for the complete installation, CLI, dashboard, Python, GPU, and analysis workflow. The [CLI reference](https://fastmdxplora.readthedocs.io/en/latest/cli_reference.html) lists the current flags, including the `v1` compatibility profile:

```
  setup  →  simulation  →  analysis  →  report
```

## Highlights

- Explore a protein's full dynamics with a single command, covering setup, simulation, analysis, and reporting
- Probe protein-ligand binding automatically with analyses for pose stability, contacts, and protein-ligand hydrogen bonds
- Reach beyond plain MD with built-in PLUMED enhanced sampling (metadynamics, umbrella sampling, steered MD) and a full analysis suite that turns trajectories into slide-ready, publication-quality figures
- Highlight user-defined residue regions on RMSF report figures for loops, helices, active-site neighborhoods, or other known intervals
- Scale from a quick single-protein exploration to large-scale parallel campaigns, driven the same way from the CLI or the Python API

## Phases of FastMDXplora

| Phase | What it does |
|---|---|
| **setup** | Cleans up your structure and builds a simulation-ready system: fixes missing atoms, adds hydrogens, solvates, and adds ions. |
| **simulation** | Runs the molecular dynamics (energy minimization, equilibration, and production), with optional enhanced sampling. |
| **analysis** | Computes the standard structural and dynamic metrics (and protein-ligand metrics when a ligand is present), with figures ready to use. |
| **report** | Packages everything into a slide deck, a written report, and a self-contained bundle you can share. |

## Installation

FastMDXplora has two dependency footprints. The **analysis** and **report**
phases are pure pip. The **setup** and **simulation** phases additionally need
OpenMM and PDBFixer, which are distributed through conda-forge. Pick the route
that matches what you need.

Analysis and reporting run on Linux, macOS, and Windows. For the full
OpenMM/PDBFixer workflow use Linux, macOS, or WSL2; on Windows, native
PowerShell is fine for analysis and reporting.

### Analysis and reporting only (pip)

```bash
pip install fastmdxplora
fastmdx explore --system 1L2Y --include analysis report
```

FastMDXplora is published under two names that resolve to the same package:
`fastmdxplora` (canonical) and `fastmdx` (an alias). Either works. Both install
a real `fastmdx` console script on `PATH`, declared by `[project.scripts]`, so
it behaves the same on every platform.

### All four phases (pip plus the chemistry stack)

Create a conda environment, add OpenMM and PDBFixer from conda-forge, then pip
install FastMDXplora into it:

```bash
conda create -n fastmdxplora "python>=3.9,<3.14"
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer openmmforcefields
pip install fastmdxplora
```

Then run your first study:

```bash
fastmdx explore --system 1L2Y
```

`1L2Y` is a small trp-cage structure that exercises every phase quickly.
Replace it with any 4-character PDB ID or a path to a local `.pdb` / `.cif`.

If you already have a conda environment you want to use, only the
`conda install` and `pip install` lines are needed.

`mamba` is a drop-in, faster replacement for the conda solver and is worth
using here, since solving the OpenMM stack is exactly where the classic solver
is slow. If you have it, substitute `mamba` for `conda` above.

### From a clone (development, or the pinned environment)

The repository ships an `environment.yml` describing the full stack, so a
clone-based install does not need the separate `conda install` step:

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git
cd FastMDXplora
mamba env create -f environment.yml || conda env create -f environment.yml
conda activate fastmdxplora
pip install -e ".[test]"
```

The `-e` makes the install editable, so your changes take effect without
reinstalling, and `[test]` adds pytest. Use `pip install .` if you only want to
run it.

### conda-forge

A single-command `conda install -c conda-forge fastmdxplora`, pulling every
dependency including the chemistry stack, is planned once the recipe clears
review. Until then, use one of the routes above.

### Optional extras

```bash
pip install "fastmdxplora[ligand]"   # OpenFF small-molecule parameterization
pip install "fastmdxplora[plumed]"   # PLUMED enhanced sampling
```

The `plumed` extra also needs the conda package:

```bash
conda install -c conda-forge openmm-plumed
```

There is also an `[md]` extra that attempts OpenMM and PDBFixer through pip.
conda-forge is preferred, because PDBFixer wheels are not reliable on every
platform.

### Windows

For analysis and reporting, a virtual environment works:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install fastmdxplora
```

If PowerShell blocks activation, allow local scripts for your account:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

You can also skip activation and call the environment's Python directly:

```powershell
.\.venv\Scripts\python.exe -m fastmdxplora.cli.main info
```

For the full chemistry workflow on Windows, use WSL2 and follow the Linux
instructions there.

### Verify

```bash
fastmdx --version    # installed version
fastmdx info         # version, phase availability, detected backends
```

If `fastmdx` is not on `PATH`, the module form is the portable fallback:

```bash
python -m fastmdxplora.cli.main --version
python -m fastmdxplora.cli.main info
```

Check which OpenMM platforms are available (CPU/CUDA/OpenCL):

```bash
python - <<'PY'
import openmm as mm
plats = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print("Available platforms:", plats)
print("CUDA available" if "CUDA" in plats else "CPU-only: simulations will run on CPU")
PY
```

### Troubleshooting

**`fastmdx` is not recognized.** The console-script directory is not on `PATH`.
Check which Python has the package and where it writes scripts:

```bash
python -m pip show fastmdxplora
python -c "import sysconfig; print(sysconfig.get_path('scripts'))"
python -m fastmdxplora.cli.main --version
```

If `import fastmdxplora` works but `fastmdx` does not, use
`python -m fastmdxplora.cli.main ...` as the portable fallback. Avoid mixing
Python installations in one terminal: the interpreter used for `pip install`
should be the one used to run the module.

**`fastmdx info` reports OpenMM or PDBFixer missing.** You have the pip-only
install. Add the chemistry stack:

```bash
conda install -c conda-forge openmm pdbfixer openmmforcefields
```

Until then, `setup` and `simulation` are skipped with a clear message, while
`analysis` and `report` continue to work.

### Graphical interface

```bash
fastmdx gui
```

Opens a local browser interface for building a study, saving it as a config
file, launching a run, and watching live telemetry with a 3D molecular viewer.
See [the GUI guide](docs/gui.md) for the cluster workflow, where you
design a config in the GUI and submit it elsewhere.

### Where to go next

- **Full walkthrough, GPU notes, and troubleshooting:** [`docs/installation.md`](docs/installation.md)
- **Already installed?** Skip to [Examples](#examples).
- **Want to contribute?** See [CONTRIBUTING.md](CONTRIBUTING.md).

## Examples

### Command line

**Run the full pipeline** (setup → simulate → analyze → report):
```bash
fastmdx explore --system protein.pdb
```
**Fetch a structure from the PDB by ID** (auto-detected, fetched from RCSB):
```bash
fastmdx explore --system 1L2Y
```
**Tune per-phase options** (flags are namespaced by phase):
```bash
fastmdx explore -s protein.pdb --setup-ph 7.4 --simulate-duration-ns 100 --simulate-platform CUDA
```
For a short local smoke test before a longer production run, use the gentle preset:
```bash
fastmdx explore -s protein.pdb --include setup simulation --simulate-preset gentle --simulate-platform CPU
```
**Run only specific phases**:
```bash
fastmdx explore -s protein.pdb --include setup simulation
```
**Run a single phase** (bare flags, no phase prefix):
```bash
fastmdx setup -s protein.pdb --ph 6.5
fastmdx simulate -s protein.pdb --output run_001 --duration-ns 50 --platform CUDA
fastmdx analyze --output run_001 --analyses rmsd rmsf rg
fastmdx report --output run_001 --no-slides
```
**Drive a whole study from a config file** (`-c` and `-config` also work):
```bash
fastmdx explore --config study.yml
```
**Generate a commented config template to edit**:
```bash
fastmdx init-config -o study.yml
```

The `-s`, `-system`, and `--system` forms are equivalent; `xplore` is an alias of `explore`.

### Python API

**Run the full pipeline**:
```python
from fastmdxplora import FastMDXplora

fmdx = FastMDXplora(system="protein.pdb")
fmdx.explore()
```
**Specify options and select phases**:
```python
fmdx = FastMDXplora(system="1L2Y")          # PDB ID, fetched from RCSB
results = fmdx.explore(
    include=["setup", "simulation", "analysis"],
    options={
        "simulation": {"duration_ns": 100, "temperature_K": 310, "platform": "CUDA"},
        "analysis":   {"include": ["rmsd", "rg", "cluster"]},
    },
)
# explore() always returns a list of runs (a single study is a list of one)
for run in results:
    print(run.run_id, run.status)
    for phase in run.phases:
        print("  ", phase.name, phase.status)
```
**Run a config file** (one system, many systems, or a parameter sweep, all the same way):
```python
fmdx = FastMDXplora(config="study.yml")
fmdx.explore()
```
**Preview a run without executing** (CLI `--dry-run`, or `dry_run=True`):
```python
FastMDXplora(config="campaign.yml").explore(dry_run=True)
```

> Recommended alias: `import fastmdxplora as fastmdx`.

See [Configuration files](#configuration-files) and [Many systems and parameter sweeps](#many-systems-and-parameter-sweeps) for the YAML format, batches, sweeps, and parallel execution.

## Configuration files

For anything beyond a quick run, capture the whole study in a single YAML file instead of a long flag list. The same file drives both the CLI and the Python API. Input is always given as a `systems:` list (even for a single system), so the file looks the same whether you study one protein or a dozen.

Generate a commented template to start from:

```bash
fastmdx init-config                    # writes fastmdxplora.yml (comprehensive)
fastmdx init-config --minimal -o study.yml   # short starter
```

A `study.yml` looks like:

```yaml
systems:
  - id: protein1
    system: protein.pdb        # PDB/CIF path, 4-char PDB ID, or sequence

output: ./my_study
include: [setup, simulation, analysis, report]

setup:
  ph: 7.4
  ion_concentration_M: 0.15

simulation:
  duration_ns: 100.0         # production length (equilibration is separate)
  temperature_K: 310.0
  platform: CUDA

analysis:
  include: [rmsd, rmsf, rg, cluster]
  selection: "name CA"
  options:
    cluster:
      methods: [kmeans, hierarchical]
      n_clusters: 5

report:
  title: "My MD Study"
```

Run it from the CLI or the API:

```bash
fastmdx explore --config study.yml     # also: -c, -config
```

```python
from fastmdxplora import FastMDXplora
FastMDXplora(config="study.yml").explore()
```

With a single system and no sweep, the output uses the familiar flat layout (`my_study/setup/`, `my_study/simulation/`, …) with the usual `manifest.json` and `resolved_config.yml`. Three things make this robust:

- **Flags override the file.** `fastmdx explore --config study.yml --simulate-duration-ns 50` keeps everything in the file but runs 50 ns. Precedence is: command-line flags / API kwargs > config file > built-in defaults.
- **Strict validation.** A typo like `pH:` (wrong case) or `simulaton:` is rejected with a did-you-mean suggestion, so a misspelled key never silently runs with the default.
- **Reproducibility.** Every run writes `resolved_config.yml`, the fully-merged configuration that actually ran (defaults + file + overrides). Feed it straight back to `--config` to reproduce the study exactly.

For a quick command-line one-off, `-s/--system` is shorthand that builds a one-element `systems` list for you:

```bash
fastmdx explore -s protein.pdb --simulate-duration-ns 50
```

## Many systems and parameter sweeps

Because input is always a `systems:` list, studying several systems is just adding entries. Add a `sweep:` block to vary parameters, and FastMDXplora runs the full cross-product, each as a complete, self-contained study.

```yaml
output: ./trpcage_campaign
include: [setup, simulation, analysis, report]

systems:
  - id: trpcage1
    system: trpcage.pdb
  - id: trpcage2
    system: trpcage.pdb
    setup: { ph: 6.5 }                 # optional per-system overrides

sweep:
  simulation.temperature_K: [300, 310, 320]   # dotted phase.option → values
  simulation.pressure_bar: [1.0, 1.2]          # multiple axes → cross-product
```

That config produces 2 systems × 3 temperatures × 2 pressures = **12 runs**. When there is more than one run, each goes in its own `runs/<id>/` subdirectory, indexed by a top-level `batch_manifest.json`, with a cross-run `comparison/` report:

```
trpcage_campaign/
  batch_manifest.json
  comparison/                                        (cross-run report)
  runs/
    trpcage1__temperature_K-300__pressure_bar-1.0/   (a full study)
    trpcage1__temperature_K-300__pressure_bar-1.2/
    ...
```

Run it exactly as any other config:

```bash
fastmdx explore --config campaign.yml
```

```python
from fastmdxplora import FastMDXplora
FastMDXplora(config="campaign.yml").explore()
```

Each run is identical in structure to a single study (its own `manifest.json`, `resolved_config.yml`, and phase directories), so existing analysis tooling works per-run unchanged. Option precedence within a run is base config < per-system overrides < swept value. Typo'd sweep axes are rejected with the valid-option list, and a failed run is recorded while the others continue.

### Cross-run comparison report

After a multi-run study, FastMDXplora automatically builds a `comparison/` report at the batch root that turns a directory of runs into a single analysis:

- **Overlays:** every run's per-frame trace (RMSD, Rg, Q-value, total SASA) drawn on one set of axes, labelled by its swept value, so divergence across the sweep is visible at a glance.
- **Trends:** each run reduced to a summary scalar (e.g. mean RMSD over the trajectory) and plotted against the swept parameter, giving a structure-property relationship.
- **`comparison_summary.csv`:** one row per run with the summary scalars, ready for further analysis.
- **`comparison_report.md`:** a written report tying the figures together, with a one-line quantitative takeaway per property (e.g. *"across temperature_K 300 → 320, mean RMSD increases 0.21 → 0.23 nm"*).

It degrades gracefully (errored runs and missing analyses are skipped) and can be turned off with `report: { comparison: false }`.

### Parallel execution

By default runs execute sequentially. An optional `execution:` block runs several at once:

```yaml
execution:
  mode: parallel          # sequential (default) | parallel
  workers: 2              # how many runs at once
  devices: [0, 1]         # GPU indices: one run pinned per device
  continue_on_error: true
```

Parallelism is process-based (each run is a subprocess, required because OpenMM contexts and the GIL don't share across threads). On GPU, the safe pattern is **one run per GPU**: list your `devices` and each worker is pinned to a distinct index round-robin. Oversubscribing a single GPU is slower than running sequentially, so `workers` should not exceed the number of devices on GPU. When `workers` is unset it defaults to one per device (GPU) or the CPU count capped at the run count (CPU).

## Outputs by phase

Each phase writes to a dedicated subdirectory under the project output root, with a structured parameters manifest so every artifact is traceable to the options that produced it.

| Phase | Key outputs |
|---|---|
| `setup` | `prepared.pdb`, `solvated.pdb`, `setup_parameters.json` |
| `simulation` | `production.dcd`, `topology.pdb`, `simulation_parameters.json` |
| `analysis` | `<analysis>/*.dat`, `<analysis>/*.png`, `analysis_manifest.json` |
| `report` | `report.md`, `dashboard.html`, `slides.pptx`, `project_bundle.zip` |

## Documentation

Documentation is available at [fastmdxplora.readthedocs.io](https://fastmdxplora.readthedocs.io) and is actively expanding.

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
