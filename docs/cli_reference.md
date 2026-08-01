# CLI commands

If this is your first run, start with the [Beginner's guide](getting_started.md).
This page is the compact reference after you understand the workflow.

This page is the shortest practical reference for the `fastmdx` command. If
the console script is not on your `PATH`, replace `fastmdx` with
`python -m fastmdxplora.cli.main`.

## The commands you will use most

```bash
fastmdx info
fastmdx explore --system 1L2Y
fastmdx explore --system protein.pdb --dashboard
fastmdx gui --output ./fastmdxplora_output_20260730
```

`explore` runs `setup → simulation → analysis → report`. The `xplore` alias
does the same thing. `gui --output DIR` reopens an existing output directory;
it does not run a new simulation.

## Input and output

`--system` accepts a local `.pdb`/`.cif` path, a four-character PDB ID such as
`1L2Y`, or a supported sequence input. Choose an output directory explicitly
when you want a stable location:

```bash
fastmdx explore --system protein.pdb --output ./runs/protein_001
```

The output directory contains `setup/`, `simulation/`, `analysis/`, and
`report/`. The top-level `manifest.json` and `resolved_config.yml` record what
actually ran.

## Run only selected phases

```bash
fastmdx explore --system protein.pdb --include setup simulation
fastmdx explore --system protein.pdb --exclude report
fastmdx explore --system protein.pdb --no-report
```

`--include` and `--exclude` are mutually exclusive. The phase names are
`setup`, `simulation`, `analysis`, and `report`.

The `v1` compatibility profile is opt-in. It is not called
`compact`:

```bash
# Direct analyze command
fastmdx analyze --output ./run --compat v1

# Prefixed form under explore
fastmdx explore --system protein.pdb --include analysis report \
  --analyze-compat v1
```

## Configure each phase

On `explore`, options are prefixed by phase:

```bash
fastmdx explore --system protein.pdb \
  --setup-ph 7.4 \
  --simulate-duration-ns 10 \
  --simulate-platform CUDA \
  --analyze-analyses rmsd rmsf rg cluster \
  --report-title "My protein study"
```

When calling a phase directly, remove the prefix:

```bash
fastmdx setup --system protein.pdb --ph 7.4
fastmdx simulate --system protein.pdb --output ./run --duration-ns 10
fastmdx analyze --output ./run --analyses rmsd rmsf rg
fastmdx report --output ./run --no-slides
```

`--simulate-duration-ns` means production time. NVT and NPT equilibration are
separate and can be set with `--simulate-nvt-duration-ns`,
`--simulate-npt-duration-ns`, or explicit step counts.

## Safe first run

Use the gentle preset before a long GPU run, especially with a new structure:

```bash
fastmdx explore --system protein.pdb \
  --include setup simulation \
  --simulate-preset gentle \
  --simulate-platform CPU
```

This is a smoke test, not production science. Inspect the telemetry and
artifacts, then choose production settings deliberately.

## Live dashboard

Start monitoring before the workflow begins:

```bash
fastmdx explore --system protein.pdb \
  --output ./run \
  --dashboard
```

Useful dashboard options:

| Option | Meaning | Default |
| --- | --- | --- |
| `--dashboard-host` | Local bind address | `127.0.0.1` |
| `--dashboard-port` | Starting port | `8765` |
| `--dashboard-refresh-seconds` | Browser polling interval | `3` |
| `--dashboard-frame-interval` | Live molecular-frame cadence | telemetry cadence |
| `--dashboard-max-playback-frames` | Maximum browser playback frames | `200` |
| `--dashboard-stop-on-complete` | Stop the server when the workflow ends | off |
| `--dashboard-open-browser` | Open the URL automatically | off |

Keep the host at `127.0.0.1` unless you have a specific secured network
setup. Remote dashboard binding is for read-only monitoring; workflow launch
and local-folder opening are disabled on non-loopback addresses.

Running `fastmdx` with no subcommand opens the dashboard home screen, where
you can configure and launch one workflow at a time:

```bash
fastmdx
```

The simulation builder uses the normal CLI and writes normal FastMDXplora
artifacts. It is not a separate simulation engine.

## Config files and campaigns

Generate a starter file:

```bash
fastmdx init-config --minimal -o study.yml
fastmdx explore --config study.yml --dry-run
fastmdx explore --config study.yml
```

Use a `systems:` list for one or many systems. Add a `sweep:` block to run a
parameter campaign and configure `execution:` for sequential or parallel
execution. The [configuration guide](configuration.md) has the complete YAML
shape.

## Diagnostics

```bash
fastmdx --version
fastmdx info
fastmdx health --no-fix
```

If chemistry dependencies are missing, use the command printed by the CLI or
follow [Installation](installation.md). Analysis and report-only workflows
can run without OpenMM/PDBFixer:

```bash
pip install fastmdxplora
fastmdx explore --system 1L2Y --include analyze report
```

For setup or simulation, missing OpenMM/PDBFixer is an error and the command
returns a non-zero exit status. Install the backends in the same environment:

```bash
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer
fastmdx info
```

The Python API returns a `RunResult` with `status="error"` for the same
condition. The dashboard validates the environment before launch and displays
the installation command instead of starting a run that cannot produce MD.

For long GPU runs, read [Production and GPU runs](production.md) for driver
checks, scheduler/remote launch, offline staging, PLUMED smoke tests,
monitoring, and checkpoint recovery.
