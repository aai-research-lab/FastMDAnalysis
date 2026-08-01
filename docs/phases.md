# The four phases

A FastMDXplora study runs as four phases in sequence. Each writes to its own
subdirectory under the project output root and records a structured manifest
of the parameters, software versions, and artifacts it produced, so every
result is traceable to the options that generated it.

```
  setup  ->  simulation  ->  analysis  ->  report
```

You can run the whole pipeline at once, or restrict it to specific phases with
`--include` / `--exclude` (CLI) or the equivalent API arguments.

## setup

Prepares a raw structure into a simulation-ready system. PDBFixer repairs the
input (missing atoms and residues, protonation at a chosen pH), then the
system is solvated in a water box and neutralized with ions. A named force
field is selected (CHARMM36 by default; AMBER variants and an OpenFF
small-molecule path are available). When a ligand is supplied, it is
parameterized with OpenFF and its bound pose is clash-checked before the run
proceeds.

Key outputs: `prepared.pdb`, `solvated.pdb`, `topology.pdb`, `system.xml`,
`setup_parameters.json`.

For a ligand, use `--setup-forcefield amber-openff`, provide one or more
`--setup-ligand` files, and set `--setup-ligand-name` when the residue name is
not `LIG`. Keep `--setup-keep-heterogens` when cofactors or other nonstandard
residues must remain in the prepared structure.

## simulation

Runs molecular dynamics with OpenMM: energy minimization, NVT and NPT
equilibration, then production. Integrator, thermostat/barostat, step counts,
and reporter intervals are all configurable. Optional PLUMED enhanced sampling
(metadynamics, umbrella sampling, steered MD) is applied to the production
stage, leaving equilibration unbiased.

Key outputs: `production.dcd`, `state_final.xml`, `simulation_parameters.json`,
`energy.csv`, and (when live telemetry is enabled) `live_status.json`,
`live_metrics.csv`, `live_events.log`, `live_frame.pdb`, and a bounded
`live_frames/` history. PLUMED adds `COLVAR`, `HILLS`, and `plumed.dat`.

## analysis

Computes structural and dynamic metrics from the trajectory and renders a
figure for each. The standard suite covers RMSD, RMSF, radius of gyration,
hydrogen bonds, secondary structure, clustering, dimensionality reduction,
Q-value, and dihedrals. When a ligand is present, protein-ligand analyses run
automatically: ligand pose RMSD, protein-ligand contacts with a binding-site
fingerprint, protein-ligand hydrogen bonds, and ligand RMSF. An analysis
scope (`solute`, `protein`, `ligand`, `all`) controls which atoms each metric
considers.

Key outputs: `<analysis>/*.dat`, `<analysis>/*.png`, matching publication SVG
files, per-analysis option manifests, and `analysis_manifest.json`.

The built-in analysis names are `rmsd`, `rmsf`, `rg`, `hbonds`, `ss`, `sasa`,
`qvalue`, `dihedrals`, `cluster`, and `dimred`. The compatibility profile is
named `v1` and is opt-in:

```bash
# Direct analyze command
fastmdx analyze --output run --compat v1

# Phase-prefixed option when using explore
fastmdx explore --system protein.pdb --include analysis report \
  --analyze-compat v1
```

The profile preserves the version 1 defaults where that comparison is
required.

## report

Assembles the results into shareable deliverables: a structured Markdown
report, a static browser dashboard, a slide deck, and a self-contained
project bundle.

Key outputs: `report.md`, `dashboard.html`, `slides.pptx`,
`project_bundle.zip`, dashboard assets, and report metadata. The static
dashboard can be opened directly, while `fastmdx gui --output
<run>` provides the live local view for an existing run.

## Output layout

```text
run/
├── manifest.json
├── resolved_config.yml
├── setup/
├── simulation/
├── analysis/
└── report/
```

Do not compare plots alone when validating a run. Compare the numerical
`.dat`/`.csv` outputs and the manifests, then inspect the figures and reports.
