# Beginner's guide: run your first simulation

This page is the shortest complete path from installation to a running
FastMDXplora study. It shows the same workflow three ways:

1. the command-line interface (CLI),
2. the live dashboard, and
3. the Python API.

FastMDXplora runs four phases in order:

```text
setup -> simulation -> analysis -> report
```

A run creates one output directory containing the prepared system, simulation
files, analysis results, and reports. The CLI, dashboard, and Python API all
use the same engine and write the same artifact layout.

## Before you start: choose where to run

- **Local workstation:** use it for installation, analysis/reporting, and a
  short gentle smoke test. Use a suitable remote or dedicated GPU host for a
  long production MD run when local resources are insufficient.
- **Connected HPC or cloud GPU:** verify the GPU driver, CUDA/OpenMM
  compatibility, and the environment on the actual host before production. Use
  a short real backend smoke test first; mocked PLUMED tests are not enough.
- **Offline GPU host:** manually stage a sanitized source/package and all input
  files, then run the commands using the host's existing environment. Do not
  assume SSH, a scheduler, internet access, or RCSB downloads.
  See [Production and GPU runs](production.md).

### If you are using Windows

Choose the path that matches your goal:

| Goal | Recommended environment |
|---|---|
| Analysis, reports, Python API, dashboard development, or package development | Native Windows PowerShell or Command Prompt with a Python virtual environment. |
| Full setup and molecular-dynamics simulation | Ubuntu or another supported Linux distribution inside [WSL2](https://learn.microsoft.com/windows/wsl/install). |
| GPU simulation from Windows | Run the Linux workflow inside WSL2 when the GPU and Windows driver expose CUDA to WSL2; otherwise use a separate Linux/HPC GPU host. |

Native Windows works well for the pip-installable parts of FastMDXplora and for
the test suite. The full OpenMM/PDBFixer/OpenFF chemistry stack is most
reliable on Linux, macOS, and WSL2.

For a full workflow on Windows:

1. Install WSL2 and an Ubuntu distribution using Microsoft's instructions.
2. Open the Ubuntu terminal; do not run the Linux commands in PowerShell.
3. Keep the repository and your data inside the WSL filesystem, for example
   under `~/src`, rather than under `/mnt/c`, when possible.
4. Run the Linux commands in the next section.

For analysis-only work in native PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install fastmdxplora
fastmdx info
```

If PowerShell blocks activation, either run the commands through
`.\.venv\Scripts\python.exe` directly or follow the execution-policy
guidance in [Installation](installation.md).

## 1. Install

FastMDXplora splits into two dependency footprints. Analysis and reporting are
pure pip. Setup and simulation additionally need OpenMM and PDBFixer from
conda-forge. Use Linux, macOS, or WSL2 for the full workflow.

### Analysis and reporting only

```bash
python -m venv .venv
source .venv/bin/activate                 # Windows: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install fastmdxplora
```

That is enough to analyze existing trajectories and build reports.

### Full workflow, all four phases

Create a conda environment with the chemistry stack, then install
FastMDXplora into it:

```bash
conda create -n fastmdxplora "python>=3.9,<3.14"
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer openmmforcefields
pip install fastmdxplora
```

`mamba` is a faster drop-in replacement for the conda solver and is worth
using here. If you have no conda at all, Miniforge is the easiest source:
download the installer for your platform from
<https://conda-forge.org/miniforge/>.

### From a clone

If you want to modify FastMDXplora, or want the exact environment the project
develops against, the repository ships an `environment.yml` with the whole
stack:

```bash
git clone https://github.com/aai-research-lab/FastMDXplora.git
cd FastMDXplora
mamba env create -f environment.yml || conda env create -f environment.yml
conda activate fastmdxplora
pip install -e ".[test]"
```

Optional tooling for specific workflows, such as OpenFF ligand
parameterization or an OpenMM-PLUMED build, is documented in
[Configuration](configuration.md) and the
[production guide](production.md).

## 2. Verify the environment

Always run these before a simulation:

```bash
fastmdx --version
fastmdx info
fastmdx health --no-fix
```

If the console script is not on `PATH`, use the module form with the same
arguments:

```bash
python -m fastmdxplora.cli.main --version
python -m fastmdxplora.cli.main info
```

For setup/simulation, `fastmdx info` must show OpenMM and PDBFixer as
installed. For a GPU run, also check the platforms in the exact environment
that will run the job:

```bash
python - <<'PY'
import openmm as mm
print("OpenMM platforms:", [
    mm.Platform.getPlatform(i).getName()
    for i in range(mm.Platform.getNumPlatforms())
])
PY
```

A CUDA platform appearing in this list is necessary but not sufficient: a
real CUDA/OpenMM context must initialize successfully on the target host.

## 3. Run a short CLI smoke test

Use a local structure when working offline. `1L2Y` is convenient when RCSB
access is available.

```bash
fastmdx explore \
  --system 1L2Y \
  --output runs/first_smoke \
  --include setup simulation \
  --simulate-preset gentle \
  --simulate-platform CPU
```

Or with a local file:

```bash
fastmdx explore \
  --system input/protein.pdb \
  --output runs/first_smoke \
  --include setup simulation \
  --simulate-preset gentle \
  --simulate-platform CPU
```

The `gentle` preset is a bounded smoke test. It is not a production
trajectory and should not be used as scientific evidence. Inspect the output
before moving to a longer run:

```text
runs/first_smoke/
├── manifest.json
├── resolved_config.yml
├── setup/
└── simulation/
    ├── production.dcd
    ├── topology.pdb
    ├── energy.csv
    ├── simulation.log
    ├── state_final.xml
    └── checkpoint.chk
```

## 4. Run the complete CLI workflow

Once the smoke test works, run all four phases. Set an explicit output path so
you can find the run again:

```bash
fastmdx explore \
  --system input/protein.pdb \
  --output runs/protein_001 \
  --simulate-platform CUDA \
  --simulate-device-index 0 \
  --simulate-duration-ns 100 \
  --simulate-checkpoint-interval-steps 10000
```

Important details:

- `--simulate-duration-ns` controls **production** length. NVT and NPT
  equilibration are separate settings.
- `--simulate-platform CUDA` requires a working CUDA/OpenMM backend. Do not
  silently replace it with CPU for a production run.
- Use one device per independent GPU run. On a two-GPU machine, use device
  indices `0` and `1` with separate output directories.
- Checkpoints help recovery but do not prove that a run completed.
- For long GPU runs, use the procedures in [Production and GPU runs](production.md).

Run only selected phases when continuing an existing workflow:

```bash
# Prepare and simulate now; analyze later
fastmdx explore --system input/protein.pdb \
  --output runs/protein_001 \
  --include setup simulation

# Analyze and report an existing output directory
fastmdx analyze --output runs/protein_001
fastmdx report --output runs/protein_001
```

## 5. Use the live dashboard with the CLI

The easiest dashboard workflow is to start it together with `explore`:

```bash
fastmdx explore \
  --system input/protein.pdb \
  --output runs/protein_dashboard \
  --include setup simulation analysis report \
  --simulate-preset gentle \
  --dashboard
```

The command prints a URL such as `http://127.0.0.1:8765`. Open that URL in a
browser. The dashboard shows setup/simulation status, telemetry, energy and
temperature, recent events, a molecular viewer, files, analyses, and reports.

The default bind address is loopback (`127.0.0.1`). Keep it that way unless
you understand the network and access-control implications.

Useful options:

```bash
fastmdx explore --system input/protein.pdb \
  --output runs/protein_dashboard \
  --dashboard \
  --dashboard-port 8770 \
  --dashboard-refresh-seconds 3 \
  --dashboard-open-browser
```

If a run already exists, serve it without starting a new simulation:

```bash
fastmdx gui \
  --output runs/protein_dashboard \
  --open-browser
```

The static report dashboard is different:

```text
runs/protein_dashboard/report/dashboard.html
```

Open that file after the report phase for a shareable completed-results view.
The live server is useful while a run is active; the static HTML dashboard is
for completed artifacts.

### Dashboard-first mode

Running the CLI with no subcommand opens the local simulation builder:

```bash
fastmdx
```

Use **Validate** before **Start Exploration**. The builder launches the same
canonical `explore` workflow as the CLI; it is not a separate simulation
engine. It can only launch a workflow when the server is loopback-bound.

Set `FASTMDX_NO_BROWSER=1` when working headlessly:

```bash
FASTMDX_NO_BROWSER=1 fastmdx
```

## 6. Run the same workflow from Python

A single study uses `FastMDXplora(system=...)`:

```python
from fastmdxplora import FastMDXplora

study = FastMDXplora(
    system="input/protein.pdb",
    output_dir="runs/python_study",
    options={
        "simulation": {
            "preset": "gentle",
            "platform": "CPU",
        },
    },
)

results = study.explore()
run = results[0]
print(run.status, run.output_dir)
for phase in run.phases:
    print(phase.name, phase.status, phase.artifacts)
```

For a production-sized GPU run, change the simulation options deliberately:

```python
study = FastMDXplora(
    system="input/protein.pdb",
    output_dir="runs/python_gpu",
    options={
        "simulation": {
            "platform": "CUDA",
            "device_index": "0",
            "duration_ns": 100.0,
            "checkpoint_interval_steps": 10000,
        },
    },
)
results = study.explore()
```

The Python API returns a list even for one system. Check `run.status` and each
phase status before using the artifacts.

### Python phase-by-phase control

Use the same output directory when you intentionally want separate phase
calls:

```python
from fastmdxplora import FastMDXplora

study = FastMDXplora(
    system="input/protein.pdb",
    output_dir="runs/python_phases",
)

setup = study.setup(ph=7.0, forcefield="charmm36")
print("setup:", setup.status)

simulation = study.simulate(
    platform="CPU",
    preset="gentle",
)
print("simulation:", simulation.status)

analysis = study.analyze(include=["rmsd", "rmsf", "rg"])
print("analysis:", analysis.status)

report = study.report(slides=False)
print("report:", report.status)
```

For ordinary work, `study.explore()` is safer because it preserves the phase
order and writes the complete run manifest automatically.

### Python API plus dashboard monitoring

The Python API does not replace the dashboard server. Enable telemetry in the
simulation options, run the study, and serve its output from a terminal:

```python
from fastmdxplora import FastMDXplora

study = FastMDXplora(
    system="input/protein.pdb",
    output_dir="runs/python_dashboard",
    options={"simulation": {
        "preset": "gentle",
        "platform": "CPU",
        "live_telemetry": True,
    }},
)
study.explore()
```

In another terminal, while the run is active or afterward:

```bash
fastmdx gui --output runs/python_dashboard
```

## 7. Analyze with the version 1 compatibility profile

The profile is named `v1`. It reproduces the analysis settings of
FastMDXplora version 1 so published results stay reproducible.

For a direct analysis command, use `--compat v1`:

```bash
fastmdx analyze \
  --output runs/bpti_reference \
  --trajectory trajectory/production.dcd \
  --topology trajectory/topology.pdb \
  --compat v1
```

When using `explore`, the same option is phase-prefixed:

```bash
fastmdx explore \
  --system input/bpti.pdb \
  --output runs/bpti_reference \
  --include analysis report \
  --analyze-trajectory trajectory/production.dcd \
  --analyze-topology trajectory/topology.pdb \
  --analyze-compat v1
```

The profile applies the published version 1 defaults: protein
scope, protein selection, every-second-frame loading, RMSD alignment against
frame 0, PCA, and hierarchical clustering with six clusters, plus the
compatibility defaults for the standard analyses.

It is a comparison aid, not proof of numerical reproduction. Compare the
`.dat`/`.csv` arrays, topology/selection rules, frame spacing, and analysis
metadata before comparing plots.

Other analysis controls are explicit:

```bash
fastmdx analyze --output runs/protein_001 \
  --analyses rmsd rmsf rg cluster \
  --scope protein \
  --stride 2 \
  --cluster-methods hierarchical \
  --cluster-n-clusters 6
```

## 8. Protein-ligand and PLUMED runs

For a ligand, use the ligand-capable force field and provide an SDF or MOL2:

```bash
fastmdx explore \
  --system input/protein_ligand.pdb \
  --output runs/ligand_001 \
  --setup-forcefield amber-openff \
  --setup-ligand input/ligand.sdf \
  --setup-ligand-name LIG \
  --simulate-platform CUDA \
  --simulate-duration-ns 10
```

For PLUMED, provide a script with valid atom indices and a fresh output
folder:

```bash
fastmdx explore \
  --system input/protein.pdb \
  --output runs/plumed_smoke \
  --include setup simulation \
  --simulate-platform CUDA \
  --simulate-preset gentle \
  --simulate-plumed-script input/plumed.dat
```

PLUMED is applied to production, not equilibration. Verify that `plumed.dat`,
`COLVAR`, and any `HILLS` file are written. A source-level or mocked PLUMED
test is not a real CUDA/OpenMM/`PlumedForce` smoke test.

## 9. What to inspect after a run

At minimum, check:

```bash
python - <<'PY'
import json
from pathlib import Path

root = Path("runs/first_smoke")
print("manifest:", json.loads((root / "manifest.json").read_text()).get("status"))
for path in [
    root / "resolved_config.yml",
    root / "simulation" / "production.dcd",
    root / "simulation" / "energy.csv",
    root / "simulation" / "simulation.log",
    root / "simulation" / "checkpoint.chk",
]:
    print(path, "OK" if path.exists() else "MISSING")
PY
```

A successful command is not enough: inspect the manifest, logs, telemetry,
trajectory, and analysis/report artifacts. Preserve the output directory as
the reproducibility record.

## 10. Common mistakes

- Running `--simulate-duration-ns 100` as a first test. Use `--simulate-preset gentle` first.
- Assuming `--simulate-platform CUDA` works because the package imports. Test a real CUDA context.
- Guessing at the profile flag. The valid forms are `--compat v1` and `--analyze-compat v1`.
- Reusing an output directory and overwriting valid trajectory or PLUMED history. Use a new directory for each attempt.
- Running a long simulation on an underpowered local workstation instead of a
  suitable production host.
- Using an offline GPU host without manually staging all dependencies and input files.
- Treating the static `report/dashboard.html` as a live monitor. Use `fastmdx gui` for live telemetry.

## Next pages

- [Installation and troubleshooting](installation.md)
- [CLI reference](cli_reference.md)
- [Live dashboard](gui.md)
- [Production and GPU runs](production.md)
- [Configuration files and campaigns](configuration.md)
- [The four phases](phases.md)
- [Python API reference](api.md)
