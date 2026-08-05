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
fingerprint, protein-ligand hydrogen bonds, ligand RMSF, and the typed
interactions that say what holds the ligand rather than how much it touches.
An analysis scope (`solute`, `protein`, `ligand`, `all`) controls which atoms
each metric considers -- for the analyses it applies to. A protein-ligand
measure works out both sides from the ligand's residue name and has nothing to
apply a general selection to, so it does not offer one.

Key outputs: `<analysis>/*.dat`, `<analysis>/*.png`, matching publication SVG
files, per-analysis option manifests, and `analysis_manifest.json`.

The built-in analysis names are `rmsd`, `rmsf`, `rg`, `hbonds`, `ss`, `sasa`,
`qvalue`, `dihedrals`, `cluster`, and `dimred`, plus `ligand_rmsd`,
`ligand_rmsf`, `pl_contacts`, `pl_hbonds`, and `pl_interactions` for a
protein-ligand complex. Each accepts its own settings
under `options`, which `fastmdx analyze --help` lists. A setting an analysis
does not have stops the run rather than being ignored: a measurement you asked
for and did not get is worse than a stopped run, because it looks like an
answer to the question you asked.

### What each measure is, precisely

The definition matters more than the name, and several of these are defined
more than one way in the literature. What FastMDXplora computes:

- **Radius of gyration** is mass-weighted, measured from the centre of mass,
  as GROMACS's `gyrate` and cpptraj's `radgyr` report it. `mass_weighted=false`
  gives the unweighted quantity, where every atom counts alike.
- **Q-value** follows Best, Hummer & Eaton (PNAS 2013) including its switching
  function: each native contact is judged against the distance it had natively
  and stops counting gradually rather than at a step. `beta` and
  `lambda_factor` are the paper's, and exposed. Q at the reference frame is
  therefore slightly under 1, which is inherent to a smooth measure and left
  unnormalised as published.
- **Hydrogen bonds** are counted in every frame they occur in. Baker-Hubbard
  proposes bonds above an occupancy threshold; `candidate_freq` defaults to 0
  so every bond that occurs is evaluated, because a bond present in five per
  cent of frames is present in those frames. `freq` reports how many bonds are
  persistent rather than deciding which are seen at all.
- **Clustering** compares frames by pairwise RMSD, each pair superposed
  optimally, so the distance does not depend on where the molecule sits.
  `features="coordinates"` superposes onto the first frame and compares
  directly instead -- the cheaper approximation, scaled so a distance is still
  an RMSD in nm.
- **Secondary structure** covers the protein only. DSSP returns a column for
  every residue and marks the others `NA`; a ligand has no secondary structure
  and does not appear.
- **Protein-ligand interactions** are typed rather than counted: hydrophobic
  contacts, hydrogen bonds, salt bridges, pi-stacking face-to-face and
  edge-to-face, pi-cation, halogen bonds, metal coordination and water
  bridges. Each is the published criterion, cited in its own docstring, and
  each threshold is a setting because the published values disagree.

  Three things about it are worth knowing before reading its output.

  *The chemistry is resolved before the geometry is measured.* Whether a
  nitrogen donates or a ring is aromatic is chemistry, not coordinates. Where
  the run's own setup phase settled it, that is used; otherwise a supplied
  SDF, then the Chemical Component Dictionary, then inference from the
  coordinates. Which route succeeded is recorded in the analysis options,
  because an interaction computed from inferred bond orders is a weaker claim
  than one computed from chemistry that was resolved.

  *Some interactions are refused rather than guessed.* Salt bridges and
  pi-cation interactions are claims about charge, and a ligand's charge
  inferred from coordinates is ambiguous more often than not -- guanidinium is
  +1, and -1 also balances. Where the charge was not determined those two are
  not reported, and the reason is recorded. The rest of the analysis
  continues: a ligand whose charge is unknown still has hydrogen bonds.

  *Occupancy is reported with the observation behind it.* A contact present in
  three frames of five hundred and one present in four hundred and fifty are
  both "present", and the fraction alone does not say which is which. Two
  contacts can even share a fraction and differ entirely: one present in 450
  consecutive frames formed once and stayed; one present in 450 alternating
  frames formed and broke 450 times. The number of separate appearances is
  reported alongside the fraction, the standard error counts those rather than
  frames, and thinly observed contacts are drawn hollow on the figure.

  The same applies to binding modes. Transitions between them are counted
  always, but the probabilities are given only where enough changes were seen
  for a rate to mean anything -- fewer than ten carries an uncertainty larger
  than itself, and would read as a measurement of kinetics the trajectory
  cannot support.

- **Secondary structure** uses MDTraj's DSSP implementation only. Version 1
  could also shell out to an external ``mkdssp`` binary; that route is not
  offered here, because it would make a system package a dependency of an
  analysis that already works without one. Where the two disagree, the
  disagreement is a finding about DSSP implementations and worth reporting
  rather than configuring around.

- **Contacts and hydrogen bonds** measure across the periodic boundary when
  the trajectory carries a unit cell. A solvated trajectory is not always
  imaged, and a molecule split across the boundary looks far from everything
  it is touching. `periodic=false` measures plain distances regardless.

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

Assembles the results into shareable deliverables: a written report in
Markdown and PDF, a static browser dashboard, a slide deck, and a
self-contained project bundle.

The report opens with a **methods paragraph** rather than a list of settings.
A methods section for a molecular dynamics study has to state a particular
list of things -- where the coordinates came from, how protonation was
decided, the force field and its version, the water model, the box, the ions,
the ensembles, the integrator and timestep, the thermostat and barostat with
their coupling, the cutoffs, the constraints, and how long each phase ran --
and that list is published, in the Journal of Chemical Information and
Modeling's reporting guidelines and in Communications Biology's
reproducibility checklist. Every value is already recorded by the run, so the
paragraph is assembled from them: nothing is invented, and anything the run
did not record is named as missing rather than filled in with what is usual.
The complete list of settings follows it, for repeating a run exactly rather
than describing it.

It then reports **convergence**, which is a statement about how much
independent information the trajectory holds. A frame is not an observation:
consecutive frames are nearly the same structure, so the number of independent
observations is set by how quickly a measure forgets where it was, not by how
often frames were written. On a five-thousand-frame trajectory of Trp-cage the
RMSD holds about twenty independent observations, so its uncertainty is
fourteen times what counting frames would give. Where a run is too short to
say -- and a short one usually is -- the section says what it cannot support
rather than printing a number for everything.

Key outputs: `report.md`, `report.pdf`, `dashboard.html`, `slides.pptx`,
`project_bundle.zip`, dashboard assets, and report metadata. The PDF needs the
`pdf` extra; where it is absent the run says so and writes the rest. The static
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
