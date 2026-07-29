# FastMDXplora usage examples

Worked examples for both the command-line interface (`fastmdx`) and the
Python API. Everything here runs the same four-phase pipeline — setup →
simulation → analysis → report — whether you drive it from one flag or a
hundred-run sweep.

A note on input: a system can be a **PDB/CIF file path**, a **4-character
PDB ID** (fetched from RCSB, e.g. `1L2Y`), or a **one-letter sequence**.
The form is auto-detected, so there is no separate `--pdb-id` flag.

---

## Command-line interface

The examples use the installed `fastmdx` console script. If Windows
PowerShell says `fastmdx` is not recognized, use the module entrypoint with
the same arguments:

```powershell
python -m fastmdxplora.cli.main info
python -m fastmdxplora.cli.main explore --system protein.pdb
```

That usually means FastMDXplora is importable but the Python console-script
directory is not on PATH. See the installation troubleshooting section for
PATH checks.

### The simplest run

Run the whole pipeline on a structure file:

```bash
fastmdx explore --system protein.pdb
```

Or fetch a structure from the PDB by ID:

```bash
fastmdx explore --system 1L2Y
```

`-s`, `-system`, and `--system` are all equivalent — the single-dash
`-system` form matches the GROMACS/AMBER/NAMD convention. The `xplore`
alias works anywhere `explore` does, for the X-branding:

```bash
fastmdx xplore -s 1L2Y
```

Output lands in `./fastmdxplora_output_<timestamp>/` unless you set
`--output`.

### Tuning the run with flags

Per-phase options are namespaced by phase (`--setup-…`, `--simulate-…`,
`--analyze-…`, `--report-…`):

```bash
fastmdx explore --system protein.pdb \
    --output ./trpcage_study \
    --setup-ph 7.4 \
    --setup-ion-concentration-M 0.15 \
    --simulate-duration-ns 100.0 \
    --simulate-temperature-K 310.0 \
    --simulate-platform CUDA \
    --analyze-analyses rmsd rmsf rg cluster \
    --report-title "Trp-cage at 310 K"
```

`--simulate-duration-ns` is **production** length; equilibration (NVT/NPT)
is independent and has its own defaults.

### Choosing phases

Run only part of the pipeline with `--include` (allowlist) or
`--exclude` (denylist) — they are mutually exclusive:

```bash
# Only prepare and simulate; analyze later
fastmdx explore -s protein.pdb --include setup simulation

# Everything except the report
fastmdx explore -s protein.pdb --exclude report

# Convenience flag for the common case
fastmdx explore -s protein.pdb --no-report
```

### Running a single phase

Each phase is also its own subcommand. Here the per-phase flags are bare
(no `--simulate-` prefix), since the phase is already chosen:

```bash
fastmdx setup    --system protein.pdb --ph 6.5 --box-shape octahedron
fastmdx simulate --system protein.pdb --output ./trpcage_study --duration-ns 50.0 --platform CUDA
fastmdx analyze  --output ./trpcage_study --analyses rmsd rg --selection "name CA"
fastmdx report   --output ./trpcage_study --no-slides
```

To reproduce the published FastMDAnalysis BPTI case-study settings with the
FastMDXplora CLI, use the explicit compatibility profile. It selects protein
scope, every second frame, PCA, and six-state hierarchical clustering while
leaving FastMDXplora's normal defaults unchanged:

```bash
fastmdx analyze \
  --system bpti.pdb \
  --output ./bpti_fastmdxplora_paper \
  --trajectory bpti.dcd \
  --topology bpti.pdb \
  --compat fastmdanalysis
```

The profile runs RMSD, RMSF, radius of gyration, hydrogen bonds, SASA,
secondary structure, PCA, and hierarchical clustering with `n_clusters=6`.
It uses the full `protein` atom selection, aligned RMSD against frame 0, and
per-atom RMSF to match FastMDAnalysis when its BPTI options file omits an atom
selection. Individual method flags such as `--cluster-n-clusters` can override the
profile when testing a deliberate variant.
profile when testing a deliberate variant.

Pointing later phases at the same `--output` lets them pick up the
artifacts the earlier phases wrote. `analyze` and `report` can infer the
system from an existing run manifest; `setup` and `simulate` still need
`--system` or `--config`.

The report phase writes `report/report.md`, `report/dashboard.html`,
`report/slides.pptx`, and `report/project_bundle.zip`. Open
`<output>/report/dashboard.html` directly in a browser for the static
dashboard view.

For local live monitoring, add `--dashboard` to the workflow command. This
starts the localhost dashboard before the workflow begins and automatically
enables live telemetry when simulation runs:

```bash
fastmdx explore --system protein.pdb --output ./trpcage_study --include setup simulation --dashboard
```

PowerShell fallback when `fastmdx` is not on PATH:

```powershell
python -m fastmdxplora.cli.main explore --system protein.pdb --output .\trpcage_study --include setup simulation --dashboard
```

`fastmdx dashboard serve --output ./trpcage_study` remains available when you
want to reopen an existing run folder manually.

### MD engine controls

Integrator, pressure (in bar **or** atm), GPU device, and checkpointing:

```bash
fastmdx explore -s protein.pdb \
    --simulate-integrator langevin_middle \
    --simulate-timestep-fs 2.0 \
    --simulate-pressure-atm 1.0 \
    --simulate-device-index 0 \
    --simulate-checkpoint-interval-steps 5000
```

Supported integrators: `langevin_middle` (default), `langevin`,
`brownian`, `verlet`, `variable_langevin`, `variable_verlet`. Pressure can
be given as `--simulate-pressure-bar` or `--simulate-pressure-atm`; atm is
converted to OpenMM's native bar (1 atm = 1.01325 bar).

### Gentle simulation smoke test

Freshly repaired or very small PDB systems can be numerically fragile at full
room-temperature MD settings. For a conservative first simulation, use the
gentle preset:

```bash
fastmdx explore --system protein.pdb --output run_gentle --include setup simulation --simulate-preset gentle --simulate-platform CPU
```

The preset uses a 0.5 fs timestep, 100 K, 5/ps Langevin friction, no NPT, and
short NVT/production stages. The equivalent explicit command is:

```bash
fastmdx explore --system protein.pdb --output run_gentle --include setup simulation --simulate-duration-ns 0.001 --simulate-nvt-steps 1000 --simulate-npt-steps 0 --simulate-production-steps 1000 --simulate-timestep-fs 0.5 --simulate-temperature-K 100 --simulate-friction-per-ps 5.0 --simulate-platform CPU --simulate-precision double
```

### Skipping PDBFixer

If you already have a prepared structure, skip the fixer:

```bash
fastmdx explore -s raw.pdb --setup-fixed-pdb prepared.pdb
```

### Config files

For anything beyond a quick run, put it all in a YAML file. Generate a
fully-commented template to edit:

```bash
fastmdx init-config                          # writes fastmdxplora.yml
fastmdx init-config --minimal -o study.yml   # short starter
fastmdx init-config -o study.yml --force     # overwrite an existing file
```

Then run it:

```bash
fastmdx explore --config study.yml           # -c and -config also work
```

Flags still override the file, so you can reuse one config and tweak a
value per invocation:

```bash
fastmdx explore --config study.yml --simulate-duration-ns 50
```

**Preview without running** — print the plan (runs, systems, swept
values, output directories, phases) and exit:

```bash
fastmdx explore --config campaign.yml --dry-run
```

Every run writes `resolved_config.yml` — the exact merged configuration
that ran — so you can reproduce it later with
`fastmdx explore --config resolved_config.yml`.

### Other commands

```bash
fastmdx info        # versions, detected backends (OpenMM/PDBFixer), citation
fastmdx --cite      # just the citation
fastmdx --version
```

---

## Config file format

Input is always a `systems:` list — even for one system — so the file
looks the same shape whether you study one protein or a dozen.

### A single study

```yaml
# study.yml
systems:
  - id: trpcage
    system: trpcage.pdb        # path, PDB ID, or sequence

output: ./trpcage_study
include: [setup, simulation, analysis, report]

setup:
  ph: 7.4
  ion_concentration_M: 0.15

simulation:
  duration_ns: 100.0
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
  title: "Trp-cage at 310 K"
```

With one system and no sweep, output uses the familiar flat layout
(`trpcage_study/setup/`, `trpcage_study/simulation/`, …).

### Several systems

Add entries to the list. Each can carry its own per-phase overrides:

```yaml
# compare.yml
systems:
  - id: wildtype
    system: wt.pdb
  - id: mutant
    system: mutant.pdb
    setup: { ph: 6.5 }         # this system only

output: ./comparison
simulation:
  duration_ns: 50.0            # shared by all systems
```

### A parameter sweep

A `sweep:` block varies parameters across runs. Each axis is a dotted
`phase.option` key mapped to a list of values; multiple axes form the
full cross-product:

```yaml
# campaign.yml
systems:
  - id: trpcage1
    system: trpcage.pdb
  - id: trpcage2
    system: trpcage.pdb
    setup: { ph: 6.5 }

output: ./trpcage_campaign

sweep:
  simulation.temperature_K: [300, 310, 320]
  simulation.pressure_bar: [1.0, 1.2]
```

This is 2 systems × 3 temperatures × 2 pressures = **12 runs**. With more
than one run, each goes in `runs/<id>/` and a top-level
`batch_manifest.json` indexes them all:

```
trpcage_campaign/
  batch_manifest.json
  runs/
    trpcage1__temperature_K-300__pressure_bar-1.0/
    trpcage1__temperature_K-300__pressure_bar-1.2/
    ...
```

Within each run, option precedence is: base config < per-system override
< swept value.

### Parallel execution

By default runs go one at a time. An `execution:` block runs several at
once:

```yaml
execution:
  mode: parallel        # sequential (default) | parallel
  workers: 2            # how many runs at once
  devices: [0, 1]       # GPU indices — one run pinned per device
  continue_on_error: true
```

On GPU the safe pattern is one run per GPU: list your `devices` and each
worker is pinned to a distinct index. Don't set `workers` higher than the
number of devices on GPU — oversubscribing one GPU is slower than running
sequentially. When `workers` is unset it defaults to one per device (GPU)
or the CPU count capped at the run count (CPU).

---

## Python API

### A single study

```python
from fastmdxplora import FastMDXplora

fmdx = FastMDXplora(system="protein.pdb")
results = fmdx.explore()

for r in results:
    print(r.run_id, r.status)
    for phase in r.phases:
        print("  ", phase.name, phase.status)
```

The recommended import alias mirrors the CLI name:

```python
import fastmdxplora as fastmdx
fastmdx.FastMDXplora(system="1L2Y").explore()
```

### With options and phase selection

`options` is keyed by phase; `explore()` takes `include`/`exclude` and a
`report` convenience flag:

```python
from fastmdxplora import FastMDXplora

fmdx = FastMDXplora(
    system="1L2Y",                       # fetched from RCSB
    output_dir="./trpcage_study",
    options={
        "setup":      {"ph": 7.4, "ion_concentration_M": 0.15},
        "simulation": {"duration_ns": 100.0, "temperature_K": 310.0,
                       "platform": "CUDA", "integrator": "langevin_middle"},
        "analysis":   {"include": ["rmsd", "rmsf", "rg", "cluster"]},
    },
)

results = fmdx.explore(include=["setup", "simulation", "analysis"])
run = results[0]                          # one study -> a list of one
print("run status:", run.status)
for phase in run.phases:
    print(" ", phase.name, phase.status)
```

`include`/`exclude`/`options` can be set on the constructor or passed to
`explore()`; arguments to `explore()` take precedence.

`explore()` **always returns a list of `RunResult`** — a single study is a
list of one, a sweep is a list of many. Each `RunResult` carries
`run_id`, `system`, `status` (`"ok"`/`"error"`/`"skipped"`), `output_dir`,
`sweep_values`, and `phases` (the list of `PhaseResult` for that run).
The iteration idiom is the same no matter how many runs there are:

```python
for run in results:
    print(run.run_id, run.status)
    for phase in run.phases:
        print("  ", phase.name, phase.status)
```

### Running a single phase

Each phase has a method that returns a `PhaseResult`:

```python
from fastmdxplora import FastMDXplora

fmdx = FastMDXplora(system="protein.pdb", output_dir="./study")

setup_result = fmdx.setup(ph=6.5, box_shape="octahedron")
print(setup_result.status, setup_result.artifacts)

sim_result = fmdx.simulate(duration_ns=50.0, platform="CUDA")
fmdx.analyze(include=["rmsd", "rg"], selection="name CA")
fmdx.report(slides=False)
```

A `PhaseResult` carries `name`, `status` (`"ok"`, `"skipped"`, or
`"error"`), `output_dir`, `artifacts`, and a `message`.

### Driving from a config file (one system or many)

A config file — single system, several systems, or a sweep — runs through
the same `FastMDXplora(config=...).explore()` interface. A single-system
config writes the flat layout; many runs go in `runs/<id>/`.

```python
from fastmdxplora import FastMDXplora

# One study from a file
FastMDXplora(config="study.yml").explore()

# A whole campaign (systems × sweep) from a file — same interface
results = FastMDXplora(config="campaign.yml").explore()

for run in results:
    print(run.run_id, run.status, run.sweep_values)
```

`explore()` returns the same `list[RunResult]` here as for a single
study — one element per run. Each carries `run_id`, `system`, `status`,
`output_dir`, `sweep_values`, and its `phases`.

### Building a config in code

You don't need a file on disk — pass a config dict directly with
`config_data`:

```python
from fastmdxplora import FastMDXplora

config = {
    "output": "./scan",
    "include": ["setup", "simulation", "analysis"],
    "systems": [
        {"id": "trpcage", "system": "trpcage.pdb"},
    ],
    "sweep": {
        "simulation.temperature_K": [290, 300, 310, 320],
    },
    "execution": {"mode": "parallel", "workers": 2, "devices": [0, 1]},
}

results = FastMDXplora(config_data=config).explore()
n_ok = sum(r.status == "ok" for r in results)
print(f"{n_ok}/{len(results)} runs succeeded")
```

### Previewing a run with `--dry-run`

To see exactly what a config will do — every run, its system, swept
values, output directory, and the phases that will execute — without
running anything, use a dry run. On the CLI:

```bash
fastmdx explore --config campaign.yml --dry-run
```

In Python, pass `dry_run=True`:

```python
from fastmdxplora import FastMDXplora

planned = FastMDXplora(config="campaign.yml").explore(dry_run=True)
for run in planned:
    print(run.run_id, run.sweep_values, "->", run.output_dir)
    # run.status == "planned"; nothing was executed
```

A dry run prints the plan and returns a `list[RunResult]` with status
`"planned"` and no populated phases. Nothing is written to disk.

---

## Cross-run comparison report

When a study has more than one run, FastMDXplora automatically builds a
`comparison/` report at the batch root that aggregates the runs:

```
my_campaign/
  batch_manifest.json
  comparison/
    overlay_rmsd.png          # all runs' RMSD traces on one axes
    overlay_rg.png
    trend_rmsd.png            # mean RMSD vs the swept parameter
    trend_rg.png
    comparison_summary.csv    # one row per run, summary scalars
    comparison_report.md      # the written report
  runs/
    ...
```

Nothing extra is required — running a sweep produces it:

```bash
fastmdx explore --config campaign.yml
```

For per-frame analyses (RMSD, Rg, Q-value, total SASA) it draws an
**overlay** of every run's trace, and — when the sweep axis is numeric — a
**trend** of each run's summary scalar against that axis. The
`comparison_summary.csv` is convenient for your own plotting:

```python
import pandas as pd
df = pd.read_csv("my_campaign/comparison/comparison_summary.csv")
print(df[["temperature_K", "rmsd_mean", "rg_mean"]])
```

To turn the report off, set it in the config's report block:

```yaml
report:
  comparison: false
```

You can also (re)build it — for instance after re-running some of the
runs, or for a batch that finished earlier — with `compare()`:

```python
from fastmdxplora import FastMDXplora

# Right after a run, compare() operates on the study just produced:
fmdx = FastMDXplora(config="campaign.yml")
fmdx.explore()
fmdx.compare()

# Or rebuild for an existing batch directory:
FastMDXplora(config="campaign.yml").compare(output_dir="my_campaign")
```

---

## Reproducibility

Every run writes a `resolved_config.yml` capturing the fully-merged
configuration that actually executed (defaults + file + overrides). It is
itself a valid config, so feeding it back reproduces the run exactly:

```bash
fastmdx explore --config some_run/resolved_config.yml
```

For a batch, `batch_manifest.json` at the output root records every run,
its swept values, status, and output directory — the index for the whole
campaign.
