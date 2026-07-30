# FastMDXplora

> **F**ully **A**utomated **Sy**s**T**em for **M**olecular **D**ynamics e**Xplora**tion

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fjcc.70350-blue)](https://doi.org/10.1002/jcc.70350)
[![PyPI version](https://img.shields.io/pypi/v/fastmdxplora)](https://pypi.org/project/fastmdxplora/)
[![Python versions](https://img.shields.io/pypi/pyversions/fastmdxplora)](https://pypi.org/project/fastmdxplora/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml/badge.svg)](https://github.com/aai-research-lab/FastMDXplora/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/aai-research-lab/FastMDXplora/branch/main/graph/badge.svg)](https://codecov.io/gh/aai-research-lab/FastMDXplora)

---

**FastMDXplora** explores a protein's behavior end to end from a single command. Given a structure (or just a PDB ID) it performs molecular dynamics exploration all the way through setup, simulation, analysis, and reporting, then hands back publication-ready results. Start with the [Beginner's guide](https://fastmdxplora.readthedocs.io/en/latest/getting_started.html) for the complete installation, CLI, dashboard, Python, GPU, and analysis workflow. The [CLI reference](https://fastmdxplora.readthedocs.io/en/latest/cli_reference.html) lists the current flags, including the `fastmdanalysis` compatibility profile:

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

FastMDXplora runs analysis/reporting on Linux, macOS, and Windows. For the full OpenMM/PDBFixer workflow, use Linux, macOS, or Linux inside WSL2. Windows users should use native PowerShell for analysis/reporting and WSL2 for the full chemistry workflow; the [Read the Docs installation page](https://fastmdxplora.readthedocs.io/en/latest/installation.html) explains the decision and gives exact commands.

### Quick install (Linux/macOS/WSL2 full workflow)

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git   # 1
cd FastMDXplora                                                  # 2
python fastmdx install                                           # 3
```

The third command:

- detects whether conda/mamba is already installed; if not, downloads and installs **Miniforge** for your platform automatically
- creates a `fastmdxplora` conda environment with Python 3.10
- installs **OpenMM** and **PDBFixer** (the only heavy chemistry dependencies)
- installs FastMDXplora itself
- runs `fastmdx info` as a smoke test to confirm everything works

> Why the `python fastmdx …` prefix? The repo includes a tiny `fastmdx` shim at its root so the CLI can run before any `pip install` step. See [`docs/installation.md § Why python fastmdx`](docs/installation.md#why-python-fastmdx) for the full behavior and the canonical `python -m fastmdxplora.cli.main …` fallback.

Then activate and run your first simulation:

```bash
conda activate fastmdxplora
fastmdx explore --system 1L2Y
```

`1L2Y` is a small trp-cage PDB that exercises every phase on a fast turnaround. Replace it with any 4-character PDB ID or with the path to a local `.pdb` / `.cif` file.

### Three scenarios, two install paths

| Starting point | What you run |
|---|---|
| **Fresh Linux/macOS/WSL2 machine** | The 3 commands below. Miniforge is downloaded and installed automatically; Python 3.9–3.12 is needed to start the bootstrap. |
| **You already have conda or mamba** | The same 3 commands. Miniforge download is skipped; the `fastmdxplora` environment is created directly. |
| **You only want analysis + reports, no MD** | Skip the `install` step. Just `pip install fastmdxplora` then `fastmdx explore --system 1L2Y --include analyze report` (2 commands, no conda env required). To upgrade later to the full chemistry stack, see [Install from PyPI](#install-from-pypi-no-git-clone) below. |

### Install from PyPI (no git clone)

If you don't want to clone the repository — for example, you only need FastMDXplora as a library to call from a script — install the published version straight from PyPI. FastMDXplora is published under two names that resolve to the same package: `fastmdxplora` (canonical) and `fastmdx` (a one-line alias). Either command installs the same software:

```bash
pip install fastmdxplora       # canonical name
pip install fastmdx            # shorter alias; same install underneath
```

The CLI is exposed by either install as a real platform-native `fastmdx` console script on `PATH` (declared by `[project.scripts]` in `pyproject.toml`), so it behaves identically on Linux / macOS / Windows — no App Execution Aliases trap, no per-OS wrapper script.

After install:

- **Analysis + reports** (no MD): all four pip deps (MDTraj, matplotlib, scikit-learn, python-pptx) are present out of the box. Run `fastmdx explore --system 1L2Y --include analyze report`.
- **Setup + simulation** (full chemistry stack): OpenMM and PDBFixer are **not** installed by plain pip because PDBFixer wheels aren't reliable across all platforms. Run `fastmdx install` after the pip install to drop them into the `fastmdxplora` conda env (auto-installs Miniforge if needed, same flow as the git-clone path).

To upgrade a pip-only install to the full chemistry stack on any OS:

```bash
pip install fastmdxplora                                      # 1
fastmdx install                                               # 2 — creates the conda env, auto-installs Miniforge if needed
conda activate fastmdxplora                                   # 3
fastmdx explore --system 1L2Y                                 # 4 — first full sim
```

Miniforge is auto-installed on a fresh machine, OpenMM + PDBFixer drop into the `fastmdxplora` conda env, and `fastmdx info` is used as a smoke test — same flow as the git-clone path.

On Windows PowerShell, use the Python launcher and the virtual environment's
activation script:

```powershell
cd C:\Users\User\OneDrive\Documents\GitHub\FastMDXplora
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e ".[test]"
python -m fastmdxplora.cli.main --version
python -m fastmdxplora.cli.main info
```

If PowerShell blocks activation, allow local scripts for your user account:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

You can also skip activation and call the environment's Python directly:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m fastmdxplora.cli.main --version
```

### Verify

```bash
fastmdx --version    # reports the installed FastMDXplora version
fastmdx info         # version + detected backends (OpenMM, PDBFixer, ...)
```

If `fastmdx` is not on PATH, these module commands are the safest fallback:

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

### Where to go next

- **Need the full walkthrough, GPU notes, or troubleshooting?** See [`docs/installation.md`](docs/installation.md).
- **Already installed and ready to run?** Skip to [Examples](#examples).
- **Want to contribute?** See [CONTRIBUTING.md](CONTRIBUTING.md).

### Troubleshooting: `fastmdx` is not recognized

On Windows, especially with Microsoft Store Python or PowerShell, you may see:

```text
fastmdx : The term 'fastmdx' is not recognized as the name of a cmdlet, function, script file, or operable program.
```

First check which Python installed FastMDXplora and where console scripts are
written:

```powershell
python -m pip show fastmdxplora
python -c "import sys; print(sys.executable)"
python -c "import sysconfig; print(sysconfig.get_path('scripts'))"
python -m fastmdxplora.cli.main --version
```

If `import fastmdxplora` works but `fastmdx` is missing, the console-script
directory is not on PATH. Use `python -m fastmdxplora.cli.main ...` as the
portable fallback, or reinstall from the same Python to recreate the script:

```powershell
python -m pip install -e .
python -m fastmdxplora.cli.main info
```

Avoid mixing multiple Python installs in the same terminal. The Python used for
`python -m pip install ...` should be the same one used for
`python -m fastmdxplora.cli.main ...`.

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
