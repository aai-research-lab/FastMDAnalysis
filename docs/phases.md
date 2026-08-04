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
input (missing atoms and residues, protonation at a chosen pH, 7.4 by
default), then the system is solvated in a water box and neutralized with
ions. The force field is chosen for you unless you name one: `auto` resolves
to `amber-openff`, which is the stack that can take a ligand.

A crystal structure carries the ligand you care about alongside the buffer
that kept the protein soluble, and a PDB record distinguishes them only by
name. Setup decides per component: a bound ligand is parameterized, a
cryoprotectant discarded, a coordinated metal kept and prepared in place. It
retrieves the bond orders, formal charges, and hydrogens a PDB record omits,
and settles the ligand's protonation from the pKa it has in the pocket rather
than in solution.

Where the structure does not determine what should be simulated, setup stops
and says what must be decided — a covalent adduct, a cofactor needing
parameters no small-molecule force field supplies, a free sugar that could be
substrate or cryoprotectant, a protonation the pH does not settle. Producing a
plausible trajectory from a guess is worse than producing none. Pass
`--setup-heterogens drop` to discard every non-standard residue instead, which
is what earlier versions did.

Key outputs: `prepared.pdb`, `solvated.pdb`, `topology.pdb`, `system.xml`,
`setup_parameters.json`, and `ligands/*.sdf` for anything parameterized.

To supply a ligand yourself rather than take it from the structure, give one or
more `--setup-ligand` files and set `--setup-ligand-name` when the residue name
is not `LIG`. A ligand supplied this way is used as given, so its protonation
and net charge are yours to state.

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
`qvalue`, `dihedrals`, `cluster`, and `dimred`. Each accepts its own settings
under `options`, which `fastmdx analyze --help` lists.

To compare against another tool, state the settings and run both. Where the
numbers agree that is a reproduction; where they do not, the difference is a
finding about the two methods rather than something to configure away.

Two analyses are worth knowing about before comparing. Clustering and
dimensionality reduction both depend on what the frames are compared in, and
FastMDXplora superposes before either: without that, the leading difference
between frames is where the molecule drifted and how it turned, not how its
shape changed. A tool that does not superpose will not agree, and the
disagreement is about method rather than settings. `--cluster-features` states
the choice for clustering.

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
