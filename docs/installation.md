# Installation

FastMDXplora has two installation paths:

- **Full setup + simulation:** Linux, macOS, or Linux inside WSL2. This path
  installs OpenMM and PDBFixer through conda-forge.
- **Analysis + reporting only:** Linux, macOS, or Windows with a normal Python
  environment. This path does not need OpenMM or PDBFixer.

Native Windows is suitable for analysis/reporting, dashboard use, and
development. The full chemistry workflow is documented for WSL2/Linux because
OpenMM, PDBFixer, AmberTools, and OpenFF are distributed and exercised most
consistently there. The bootstrap includes a Windows x86_64 Miniforge path, but
that native full-workflow path is not the recommended release path. Windows ARM
does not have a bundled Miniforge installer; use WSL2 or another supported
x86_64 environment.

This page covers every supported starting point and gives troubleshooting tips. If you just want to get going, the **Quick install: full workflow** section is enough.

## Quick install: full workflow

Use a Linux or macOS terminal. On Windows, open an Ubuntu WSL2 terminal first.

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git   # 1
cd FastMDXplora                                                  # 2
python fastmdx install                                           # 3
```

The third command (`install`) does everything else. It requires a working
Python 3.9–3.12 interpreter to start the bootstrap script; it can install
Miniforge and the remaining environment after that:

1. Detects whether conda or mamba is already on your `PATH`.
2. If not, downloads and installs **Miniforge** for your platform (Linux x86_64 / aarch64, macOS x86_64 / arm64, or Windows x86_64) into the standard Miniforge location.
3. Creates a `fastmdxplora` conda environment with **Python 3.10**.
4. Installs **OpenMM** and **PDBFixer** (the only heavy chemistry dependencies).
5. Installs **FastMDXplora** itself into that environment.
6. Runs `fastmdx info` as a smoke test.

> Why the `python fastmdx …` prefix? The repo includes a tiny `fastmdx` shim at its root that lets the CLI run before the package is installed. See the **Why `python fastmdx`?** section below for details.

Then activate the environment and run a first simulation:

```bash
conda activate fastmdxplora
fastmdx explore --system 1L2Y
```

`1L2Y` is a small trp-cage PDB that exercises every phase on a fast turnaround. Replace it with any 4-character PDB ID or with the path to a local `.pdb` / `.cif` file.

For a complete beginner walkthrough, including a safe gentle smoke test, the
live dashboard, the Python API, and GPU/offline-machine rules, continue to the
[Beginner's guide](getting_started.md).

### Create the environment from `environment.yml`

The repository also includes a pinned starting environment description. Use it
when you want to create the environment explicitly instead of using the
bootstrap command:

```bash
mamba env create -f environment.yml || conda env create -f environment.yml
conda activate fastmdxplora
python -m pip install -e .
fastmdx info
```

The environment file includes the core analysis/report dependencies and the
OpenMM/PDBFixer chemistry stack. Add the optional ligand and PLUMED packages
only when those workflows are needed; see [Production and GPU runs](production.md)
for the real-backend checks required before a long run.

### Windows analysis-only or local development install

On Windows PowerShell, the most reliable local development path is to use the
Python launcher and call pip through Python:

```powershell
cd C:\Users\User\OneDrive\Documents\GitHub\FastMDXplora
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e ".[test]"
python -m fastmdxplora.cli.main --version
python -m fastmdxplora.cli.main info
```

If activation is blocked by PowerShell's execution policy, allow local scripts
for your user account:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

This environment is appropriate for package development, analysis, reports,
and documentation tests. For setup/simulation on Windows, use WSL2 and repeat
the full-workflow installation there.

Activation is optional. You can also run the virtual environment's Python
directly:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m fastmdxplora.cli.main --version
```

### Optional extras

For ligand parameterization and PLUMED enhanced sampling in an editable
development checkout:

```bash
pip install -e ".[ligand]"   # OpenFF small-molecule parameterization
pip install -e ".[plumed]"   # PLUMED enhanced sampling
```

The `plumed` extra also requires the `openmm-plumed` conda package:

```bash
conda install -c conda-forge openmm-plumed
```

> Need to install for development instead? Use `python fastmdx install-e` — same flow but the local checkout is installed in editable mode.

## Prerequisites

- A **shell** (bash / zsh / PowerShell) with **internet access**.
- `git` on `PATH` (preinstalled on modern macOS and Windows 10+; on bare Linux you may need to install via your package manager).
- A terminal that supports **UTF-8** output (the CLI renders a box-drawing banner).
- **~1.5 GB of free disk** for the full install (Miniforge downloads ~150 MB; the `fastmdxplora` conda environment adds another ~800 MB of OpenMM / MDTraj / matplotlib / etc.).
- Python 3.9–3.12 is required to start the repository bootstrap command. The
  bootstrap then installs the dedicated Python 3.10 conda environment used by
  the full workflow.

## Scenarios your new user might be in

### Scenario A — Fresh Linux/macOS/WSL2 machine

You have a brand-new Linux/macOS machine (or a fresh WSL2 distro) that has no
conda or mamba yet. Python must still be available long enough to start
`python fastmdx install`.

Just run the three commands above. The third command detects the missing conda and downloads Miniforge for your OS from the pinned GitHub release used by the installer. Miniforge is then installed into the standard user installation location for that platform.

Time: roughly 5–10 minutes for Miniforge + ~5–10 minutes for environment
resolution. Disk cost: about 1 GB.

### Scenario B — You already have conda or mamba installed

Skip the auto-install. The same `python fastmdx install` command detects your existing conda/mamba, skips the Miniforge download, and creates the `fastmdxplora` environment directly.

If you don't yet have conda/mamba, Miniforge is the easiest source (it's conda + mamba + conda-forge preconfigured). The `install` command installs it for you, so you don't need to grab it manually.

### Scenario C — You only need analysis + reporting, not MD

```bash
pip install fastmdxplora
fastmdx explore --system 1L2Y --include analyze report
```

This installs FastMDXplora from PyPI directly into your system Python (no conda env required). The `analyze` and `report` phases only need pip-installable dependencies (MDTraj, matplotlib, scikit-learn, python-pptx), all of which are bundled.

The `setup` and `simulation` phases need OpenMM + PDBFixer. If they are
missing, FastMDXplora prints the exact install command and exits with an error
status instead of claiming that a simulation completed.

### If `fastmdx info` says OpenMM or PDBFixer is missing

This is the most common first-run issue. Install the chemistry backends in the
same conda environment that runs `fastmdx`:

```bash
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer
fastmdx info
```

Run `fastmdx info` again and confirm both backends say `installed`, then retry
the simulation. The dashboard and Python API use the same environment check
and show the same command when a launch is blocked. A missing backend is a
failed setup/simulation prerequisite, not a completed simulation.

If you only need to analyze an existing trajectory and write reports, do not
install the chemistry stack:

```bash
pip install fastmdxplora
fastmdx explore --system 1L2Y --include analyze report
```

### Scenario D — Windows full workflow

Install [WSL2](https://learn.microsoft.com/windows/wsl/install) with an Ubuntu
distribution, then run the Linux instructions inside the Ubuntu terminal. Keep
your PDB files and output directories in the WSL filesystem when possible for
better I/O performance.

### Scenario E — conda-forge (one command, when published)

> Coming soon. A single-command install is in progress via a `recipes/fastmdxplora/meta.yaml` recipe, which would give:
>
> ```bash
> conda install -c conda-forge fastmdxplora
> fastmdx explore --system 1L2Y
> ```
>
> Use Scenario A or B until the recipe clears review.

### Scenario F — Editable install (contributors hacking on FastMDXplora)

This is for users who want to **modify FastMDXplora's source** — adding a new analysis, fixing a bug, or contributing back upstream. The flow mirrors Scenarios A and B (clone the repo, `cd` into it, run `install`), but use `install-e` instead of `install` for editable mode. The local checkout is then installed in **editable mode** (`pip install -e .`) so any change you make under `src/fastmdxplora/` shows up the next time you run `fastmdx`.

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git   # 1
cd FastMDXplora                                                  # 2
python fastmdx install-e                                          # 3
```

What `install-e` does differently from `install`:

- Miniforge auto-install (if needed), conda env creation with Python 3.10, OpenMM + PDBFixer drop-in, and the `fastmdx info` smoke test are unchanged.
- The last step uses `pip install -e .` (editable) on the local repository checkout instead of pulling `fastmdxplora` from PyPI — your edits to `src/fastmdxplora/` immediately affect the next `fastmdx` invocation.

Then activate and run as usual:

```bash
conda activate fastmdxplora
fastmdx explore --system 1L2Y
```

To validate changes locally, install the test extras and run the suite. **Run these inside the `fastmdxplora` conda env from the previous step** (`conda activate fastmdxplora`) so the editable `src/` and `pytest` are on `PATH`:

```bash
pip install -e ".[test]"     # adds pytest, pytest-cov, ruff
pytest                       # full test suite
ruff check src tests         # lint with project conventions
```

For full contributor conventions (test requirements, coding style, PR workflow) see [CONTRIBUTING.md](https://github.com/aai-research-lab/FastMDXplora/blob/main/CONTRIBUTING.md).

## Verify the install

```bash
fastmdx --version    # e.g. fastmdx 2.0.1 (FastMDXplora)
fastmdx info         # version, detected phases, OpenMM/PDBFixer status, citation
```

`fastmdx info` reports which backends are present. If `PDBFixer: installed` and `OpenMM: installed` both say yes, all four phases will work end-to-end.

If the package imports but the `fastmdx` command is not recognized, the
console-script directory is probably not on PATH. This is common on Windows
with Microsoft Store Python or a mismatched PowerShell environment. Use the
module entrypoint as a robust fallback:

```powershell
python -m pip show fastmdxplora
python -c "import sys; print(sys.executable)"
python -c "import sysconfig; print(sysconfig.get_path('scripts'))"
python -m fastmdxplora.cli.main --version
python -m fastmdxplora.cli.main info
```

Reinstalling with the same Python can recreate the console script:

```powershell
python -m pip install -e .
```

Avoid mixing multiple Python installs in one terminal. The Python used for
`python -m pip install ...` should be the same Python used for
`python -m fastmdxplora.cli.main ...`.

To check whether a GPU-capable OpenMM platform is available:

```python
import openmm as mm
plats = [mm.Platform.getPlatform(i).getName()
         for i in range(mm.Platform.getNumPlatforms())]
print("Available platforms:", plats)
print("CUDA available" if "CUDA" in plats else "CPU-only: simulations will run on CPU")
```

## If something goes wrong

| Symptom | Likely cause | Fix |
|---|---|---|
| `git` not found | `git` is not installed (rare on modern Windows and macOS, but possible on stripped-down Linux installs) | Install Git from https://git-scm.com/downloads, then re-run the first command |
| `python` isn't recognized | `python` is not on `PATH`, so the `install` command can't even start | Install Miniforge manually from https://conda-forge.org/miniforge/ — it ships Python and conda together — then re-run the install command |
| `[FAILED] Python X.Y.Z is too new` from `python health.py` | Python ≥ 3.13 (the OpenMM / PDBFixer chemistry stack caps out at 3.12) | Use Python 3.10 or 3.11 — `install` defaults to 3.10, which is the smoothest supported version |
| `conda` / `mamba` not on PATH after Miniforge install | New shell didn't pick it up | On Linux/macOS: `source ~/miniforge3/etc/profile.d/conda.sh`. On Windows: open a fresh `cmd` or PowerShell |
| Simulation phase missing OpenMM | You used `pip install fastmdxplora` without installing the chemistry stack | `conda install -c conda-forge openmm pdbfixer` (recommended) or `pip install "fastmdxplora[md]"` (best-effort) |
| Self-heal prints a friendly install hint at exit 2 | A `setup` / `simulate` / `explore` command needs OpenMM and it's missing | Follow the install command in the hint, or use `--include analyze report` to skip chemistry phases |
| PDB won't download | No internet to RCSB, or your input wasn't a valid 4-character ID | Use a local `.pdb` / `.cif` path instead, or check the ID |

## Diagnostic entry points

| Command | What it does |
|---|---|
| `python fastmdx health` | Runs the repository doctor (verifies repo layout, deps, imports, runs a smoke test). Add `--no-fix` for diagnose-only mode. |
| `python fastmdx info` | Prints FastMDXplora version + detected backends. |
| `python fastmdx --version` | Version only. Available before `pip install` (uses the pure-Python shim in the repo root). |

`health` from inside a fresh clone is what's caught the highest-friction install bugs historically, so run it if anything seems off.

## Why `python fastmdx`?

The third step in the **Quick install: full workflow** is `python fastmdx install`, not bare `fastmdx install`. The reason is order-of-operations:

- When you first clone the repo, FastMDXplora is **not yet installed** on your system, so there's no `fastmdx` console script on `PATH` yet.
- The repo ships a tiny Python file called **`fastmdx`** at its root. It has no shebang, so it runs identically on Linux, macOS, and Windows. Running it as `python fastmdx <subcommand>` invokes a one-shot script that puts `src/` on `sys.path`, sets `PYTHONPATH`, and forwards to the CLI module — i.e. `python -m fastmdxplora.cli.main <subcommand>`.
- After you run `python fastmdx install` (or do `pip install fastmdxplora` first), a real `fastmdx` console script from `[project.scripts]` in `pyproject.toml` lands on `PATH`. From that point on, plain `fastmdx install` (and any other subcommand) works directly, with no `python fastmdx` prefix.
So the right command for your situation is:

| Your situation | Run this |
|---|---|
| Fresh clone, haven't run any install yet | `python fastmdx install` |
| Already ran `install` once (or `pip install fastmdxplora`) | `fastmdx install` (plain) |
| Want the canonical modulepath form (always works) | `python -m fastmdxplora.cli.main install` |

## Why Python 3.9–3.12 (and not 3.13)?

The chemistry phases depend on **OpenMM** and **PDBFixer**, which are primarily distributed through **conda-forge**. Their current wheels target Python 3.9–3.12. The `health.py` doctor and the `install` command both enforce this range from a single source of truth (`fastmdxplora.MIN_PYTHON = (3, 9)` and `MAX_PYTHON = (3, 13)`). Python 3.13 and newer are detected as out-of-range.

If your environment is too new, install Python 3.10 or 3.11 in a dedicated conda env, then re-run `install` (it already defaults to Python 3.10, the smoothest supported version).

## Where to go next

- **Ready to run?** Try the [Usage examples](usage_examples.md).
- **Need a specific output config?** See [Configuration files](configuration.md).
- **Want to write your own analyses or extend FastMDXplora?** See [Phases](phases.md) and the [API reference](api.md).
- **Want to contribute?** See [CONTRIBUTING.md](https://github.com/aai-research-lab/FastMDXplora/blob/main/CONTRIBUTING.md).
