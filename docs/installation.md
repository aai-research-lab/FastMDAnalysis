# Installation

FastMDXplora has two dependency footprints:

- **Analysis and report** need only pip-installable packages (MDTraj, NumPy,
  SciPy, scikit-learn, matplotlib, pandas, python-pptx).
- **Setup and simulation** additionally need **OpenMM** and **PDBFixer**, which
  are distributed through **conda-forge**.

Choose the route that matches what you need. Every route below is a standard
Python installation; FastMDXplora ships no bespoke installer.

## Platform support

| Platform | Analysis + report | Full four-phase workflow |
|---|---|---|
| Linux x86_64 / aarch64 | yes | yes |
| macOS Intel / Apple Silicon | yes | yes |
| Windows (native PowerShell) | yes | use WSL2 |
| Windows via WSL2 | yes | yes |

The chemistry stack is most reliable on Linux, macOS, and WSL2, which is where
conda-forge's OpenMM builds get the most use.

## Route 1: pip, analysis and reporting

```bash
pip install fastmdxplora
```

FastMDXplora is published under two names that resolve to the same package:
`fastmdxplora` (canonical) and `fastmdx` (an alias). Either command installs
the same software.

This gives you a working analysis and report pipeline, slide deck included,
with no conda required:

```bash
fastmdx explore --system 1L2Y --include analysis report
```

If you invoke `setup` or `simulation` without the chemistry stack,
FastMDXplora detects the missing backends before running, prints the exact
install command, and exits with an error status. It does not report a
simulation as completed.

## Route 2: pip plus the chemistry stack, all four phases

Create a conda environment, add the chemistry packages from conda-forge, then
pip install FastMDXplora into that environment:

```bash
conda create -n fastmdxplora "python>=3.9,<3.14"
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer openmmforcefields
pip install fastmdxplora
```

Verify and run:

```bash
fastmdx info
fastmdx explore --system 1L2Y
```

If you already have a conda environment you want to use, only the
`conda install` and `pip install` lines are needed.

### mamba

`mamba` is a drop-in, faster replacement for the conda solver. Solving the
OpenMM stack is exactly where the classic solver is slow, so it is worth
having. Substitute `mamba` for `conda` in the commands above.

If you have no conda at all, **Miniforge** is the easiest source, since it
ships conda and mamba preconfigured for conda-forge. Download the installer
for your platform from <https://conda-forge.org/miniforge/> and follow the
instructions there.

## Route 3: from a clone

Use this if you want to modify FastMDXplora, or if you want the exact pinned
environment the project develops against. The repository ships an
`environment.yml` describing the full stack, so no separate `conda install`
step is needed:

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git
cd FastMDXplora
mamba env create -f environment.yml || conda env create -f environment.yml
conda activate fastmdxplora
pip install -e ".[test]"
```

`-e` installs in editable mode, so edits under `src/fastmdxplora/` take effect
on the next `fastmdx` invocation without reinstalling. `[test]` adds pytest and
pytest-cov. Use `pip install .` if you only want to run the software.

To validate changes, run the suite inside the activated environment:

```bash
pytest                       # full test suite
ruff check src tests         # lint with project conventions
```

For contributor conventions see
[CONTRIBUTING.md](https://github.com/aai-research-lab/FastMDXplora/blob/main/CONTRIBUTING.md).

## Route 4: conda-forge

A single-command install that pulls every dependency, including the chemistry
stack, is planned:

```bash
conda install -c conda-forge fastmdxplora
```

The recipe lives in `recipes/fastmdxplora/meta.yaml` and has not yet cleared
conda-forge review. Use Route 2 or 3 until it does.

## Windows

For analysis and reporting, a virtual environment works natively:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install fastmdxplora
python -m fastmdxplora.cli.main info
```

If PowerShell blocks activation, allow local scripts for your account:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

You can also skip activation and call the environment's Python directly:

```powershell
.\.venv\Scripts\python.exe -m fastmdxplora.cli.main info
```

For the full four-phase workflow, install
[WSL2](https://learn.microsoft.com/windows/wsl/install) with an Ubuntu
distribution and follow Route 2 or 3 inside the Ubuntu terminal. Keep PDB
files and output directories on the WSL filesystem for better I/O performance.

## Optional extras

```bash
pip install "fastmdxplora[ligand]"   # OpenFF small-molecule parameterization
pip install "fastmdxplora[plumed]"   # PLUMED enhanced sampling
pip install "fastmdxplora[test]"     # pytest, pytest-cov
pip install "fastmdxplora[docs]"     # sphinx, myst-parser, rtd theme
```

The `plumed` extra also requires the conda package:

```bash
conda install -c conda-forge openmm-plumed
```

There is also an `[md]` extra that attempts OpenMM and PDBFixer through pip.
conda-forge is preferred, because PDBFixer wheels are not available on every
platform.

## Prerequisites

- **Python 3.9 to 3.13.** See [Why this range](#why-python-39-313) below.
- **~1 GB of free disk** for the full environment (OpenMM, MDTraj, matplotlib
  and their dependencies).
- **Internet access** for package downloads and, if you use PDB IDs, for
  fetching structures from RCSB.

## Verify the install

```bash
fastmdx --version    # e.g. fastmdx 2.1.0 (FastMDXplora)
fastmdx info         # version, phase availability, backend status, citation
```

`fastmdx info` reports which backends are present. If both PDBFixer and OpenMM
say `installed`, all four phases will work end to end.

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

## If something goes wrong

| Symptom | Likely cause | Fix |
|---|---|---|
| `fastmdx info` says OpenMM or PDBFixer is missing | You have the pip-only install | `conda install -c conda-forge openmm pdbfixer openmmforcefields` in the same environment that runs `fastmdx` |
| `fastmdx` is not recognized | The console-script directory is not on `PATH` | Use `python -m fastmdxplora.cli.main ...`, or reinstall from the interpreter you are actually using |
| `conda` or `mamba` not on `PATH` after installing Miniforge | The new shell did not pick it up | Linux/macOS: `source ~/miniforge3/etc/profile.d/conda.sh`. Windows: open a fresh terminal |
| Conda cannot solve the environment | Channel priority, or a Python version outside the supported range | Add `-c conda-forge`, and use a Python between 3.9 and 3.13 |
| `git` not found | Git is not installed | Install from <https://git-scm.com/downloads>, then retry the clone |
| PDB will not download | No internet route to RCSB, or the ID is not a valid 4-character code | Use a local `.pdb` / `.cif` path instead, or check the ID |
| Setup or simulation reports missing backends | Expected on a pip-only install | Either install the chemistry stack, or run with `--include analysis report` |

## Why Python 3.9-3.13?

The chemistry phases depend on OpenMM and PDBFixer, distributed through
conda-forge. conda-forge publishes OpenMM builds for Python 3.13 on every
supported platform (linux-64, win-64, osx-64, osx-arm64), and PDBFixer is a
`noarch` package, so the whole stack works across this range.

The range is declared in one place, `fastmdxplora.MIN_PYTHON` and
`MAX_PYTHON`, and mirrored by `pyproject.toml` and `environment.yml`. Python
3.14 and newer are outside the range until the chemistry stack publishes
builds for them.

## Where to go next

- [Beginner's guide](getting_started.md) for a full first run
- [CLI reference](cli_reference.md) for the complete flag list
- [Configuration files](configuration.md) for YAML-driven studies
- [Live dashboard](gui.md) to watch a run in progress
