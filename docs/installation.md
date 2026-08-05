# Installation

FastMDXplora is packaged for conda-forge with every backend it can use, so a
single command gets the whole pipeline. The routes below exist for the cases
where that is not what you want.

The dependency footprint has three tiers:

- **Analysis and report** need only pip-installable packages (MDTraj, NumPy,
  SciPy, scikit-learn, matplotlib, pandas, python-pptx).
- **Setup and simulation** additionally need **OpenMM** and **PDBFixer**.
- **Protein-ligand preparation** additionally needs **OpenFF**,
  **openmmforcefields**, **RDKit**, and **PROPKA** — to retrieve a ligand's
  chemistry, parameterize it, and settle its protonation in the binding site.

All of those are conda-forge packages. They are on PyPI too, with varying
build quality; conda-forge is where they get the most use.

## Platform support

| Platform | Analysis + report | Full four-phase workflow |
|---|---|---|
| Linux x86_64 / aarch64 | yes | yes |
| macOS Intel / Apple Silicon | yes | yes |
| Windows (native PowerShell) | yes | use WSL2 |
| Windows via WSL2 | yes | yes |

The chemistry stack is most reliable on Linux, macOS, and WSL2, which is where
conda-forge's OpenMM builds get the most use.

## Route 1: conda-forge, everything

```bash
conda create -n fastmdxplora -c conda-forge fastmdxplora
conda activate fastmdxplora
fastmdx info
```

This installs FastMDXplora with OpenMM, PDBFixer, openmmforcefields, OpenFF,
RDKit, and PROPKA, which is every backend any phase uses. Nothing further is
needed for a protein-ligand study from a PDB identifier.

Into an environment you already have:

```bash
conda install -c conda-forge fastmdxplora
```

## Route 2: pip, analysis and reporting

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

## Route 3: pip plus the chemistry stack

For an environment where FastMDXplora itself must come from PyPI — an editable
checkout, or a version not yet on conda-forge — install the backends from
conda-forge and the package with pip:

```bash
conda create -n fastmdxplora "python>=3.9,<3.14"
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer openmmforcefields \
    openff-toolkit rdkit propka
pip install fastmdxplora
```

Drop the second line of conda packages for setup and simulation without the
ligand path. The pip extras declare the same dependencies, so
`pip install "fastmdxplora[md,ligand]"` works where wheels are available for
your platform.

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

## Route 4: from a clone

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
distribution and follow Route 1 inside the Ubuntu terminal. Keep PDB
files and output directories on the WSL filesystem for better I/O performance.

## Optional extras

```bash
pip install "fastmdxplora[ligand]"   # OpenFF small-molecule parameterization
pip install "fastmdxplora[pdf]"      # the report as a PDF as well as Markdown
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
| `fastmdx info` says OpenFF, RDKit, or PROPKA is missing | The ligand stack is absent, so setup can prepare apo systems only | `conda install -c conda-forge openff-toolkit openmmforcefields rdkit propka`, or reinstall from Route 1 |
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

## The PDF report

`fastmdxplora[pdf]` installs WeasyPrint and a Markdown converter, which turn
the written report into `report.pdf` alongside `report.md`.

WeasyPrint needs Pango, Cairo and GDK-PixBuf as system libraries. The
conda-forge package brings its own, so `conda install -c conda-forge
fastmdxplora` gives the PDF with nothing further to do. On PyPI they have to be
present already — on Debian and Ubuntu that is:

```bash
sudo apt install libpango-1.0-0 libpangoft2-1.0-0 libcairo2 libgdk-pixbuf-2.0-0
```

Where they are missing, the report phase says which library it could not find
and writes the other formats. A run is not failed for want of a fifth format.
