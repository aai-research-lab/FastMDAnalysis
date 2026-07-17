# Installation

FastMDXplora's four phases have different dependency footprints. The analysis
and report phases work from pip alone; the setup and simulation phases need
PDBFixer and OpenMM, which are distributed primarily through conda-forge. The
full chemistry workflow is most reliable on Linux/macOS or in a Linux
environment such as WSL2. Pick the route that matches what you need.

## Full install (all four phases)

The setup and simulation chemistry stack (OpenMM, PDBFixer) installs most
reliably from conda-forge, so the full install uses the bundled
`environment.yml`. `mamba` is recommended (a faster conda solver); plain
`conda` works too. On Windows, use the WSL2 workflow below for this install.

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git
cd FastMDXplora
mamba env create -f environment.yml || conda env create -f environment.yml
conda activate fastmdxplora
pip install -e .
```

### Windows users

OpenMM itself supports Windows, but the full molecular simulation setup stack
used by FastMDXplora is most reliable on Linux/macOS. On Windows, we recommend
using WSL2 with Ubuntu, or Docker, especially for workflows involving
AmberTools, OpenFF, PDBFixer, or ligand parameterization. WSL2 lets you run a
Linux environment inside Windows; you do not need a separate Linux computer.

Windows users can still run protein-only OpenMM simulations with force fields
such as `Amber14` or `CHARMM36`. AmberTools and OpenFF are mainly needed for
ligand/small-molecule parameterization and other more complex setup workflows.

Recommended WSL2 workflow:

1. In PowerShell **as Administrator**, install Ubuntu for WSL2:

   ```powershell
   wsl --install -d Ubuntu
   wsl --update
   ```

   Restart if Windows asks you to, then open the Ubuntu application.

2. In the Ubuntu terminal, install conda (Miniforge includes conda and mamba):

   ```bash
   sudo apt update
   sudo apt install -y curl git
   curl -L -o "$HOME/Miniforge3.sh" \
     "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
   bash "$HOME/Miniforge3.sh" -b -p "$HOME/miniforge3"
   source "$HOME/miniforge3/etc/profile.d/conda.sh"
   conda init bash
   source ~/.bashrc
   ```

3. Clone FastMDXplora and create its environment inside WSL2:

   ```bash
   mkdir -p ~/src && cd ~/src
   git clone https://github.com/aai-research-lab/FastMDXplora.git
   cd FastMDXplora
   conda env create -f environment.yml
   conda activate fastmdxplora
   python -m pip install -e .
   ```

4. Run FastMDXplora commands from the Ubuntu/WSL2 terminal, for example:

   ```bash
   python -m fastmdxplora.cli.main --version
   fastmdx info
   ```

Docker is another suitable way to run the Linux-based workflow if you already
use Docker. If AmberTools or OpenFF dependency resolution fails on native
Windows, switch to WSL2/Linux (or a Linux-based Docker environment) instead of
trying to force a broken native installation.

### Native Windows development (limited scope)

For analysis/report-only work, documentation work, or limited OpenMM usage,
native Windows PowerShell can use a Python virtual environment. This is not the
recommended path for the complete setup, ligand, or parameterization workflow:

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

Activation is optional. You can also run the virtual environment's Python
directly:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m fastmdxplora.cli.main --version
```

### Optional extras

```bash
pip install -e ".[ligand]"   # OpenFF small-molecule parameterization
pip install -e ".[plumed]"   # PLUMED enhanced sampling
```

The `plumed` extra also requires the `openmm-plumed` conda package:

```bash
conda install -c conda-forge openmm-plumed
```

## Analysis and report only (from PyPI)

If you only need to analyze existing trajectories and build reports (no
setup or simulation), plain pip is enough, with no conda required:

```bash
pip install fastmdxplora
```

The `fastmdx` command and `import fastmdxplora` are available either way. The
short alias package `fastmdx` installs the same software:

```bash
pip install fastmdx
```

## Verifying the install

```bash
fastmdx --version
python -c "import fastmdxplora; print(fastmdxplora.__version__)"
```

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
import openmm
plats = [openmm.Platform.getPlatform(i).getName()
         for i in range(openmm.Platform.getNumPlatforms())]
print("CUDA available" if "CUDA" in plats else "CPU-only; simulations will run on CPU")
```
