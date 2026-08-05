# The four phases

```
  setup  →  simulation  →  analysis  →  report
```

Each phase writes into its own directory, records what it did, and can be run
on its own. `fastmdx explore` runs all four; `fastmdx setup`, `simulate`,
`analyze` and `report` run one.

---

## setup

Turns a structure into something that can be simulated.

It fetches by PDB identifier or reads a file, then repairs what is missing —
absent atoms, absent hydrogens, chain breaks — using PDBFixer at the pH you
give. It decides what each non-standard residue is for: a bound ligand is
parameterised, a crystallisation additive discarded, a coordinated metal kept.
Then it solvates, adds ions, and writes a serialised system.

**Where the structure does not say enough, it stops.** An unknown residue, a
ligand whose charge cannot be settled, a clash that means the pose is wrong —
each of these produces a refusal naming what could not be decided and what
would settle it, rather than a guess that runs.

For a ligand, the chemistry is resolved before anything else: from a supplied
SDF, from the Chemical Component Dictionary, or inferred from coordinates —
and which route succeeded is recorded, because an interaction computed from
inferred bond orders is a weaker claim than one from chemistry that was known.

**Writes** `prepared.pdb`, `solvated.pdb`, `system.xml`, `state.xml`,
`setup_parameters.json`, and an SDF per ligand.

---

## simulation

Runs the dynamics: energy minimisation, NVT equilibration, NPT equilibration,
production.

Beyond that it can hold parts of the system still while the solvent settles,
embed a protein in a lipid bilayer, and bias the run along a coordinate you
name. Those are covered in [Beyond a box of water](simulations.md), along with
what the phase says when a run fails.

**Writes** `production.dcd`, `energy.csv`, `checkpoint.chk`,
`simulation_parameters.json`.

---

## analysis

Measures the trajectory. Fifteen analyses, each writing its data, its figure,
and the settings it used.

**Shape and size**

| | |
|---|---|
| `rmsd` | how far the structure has moved from a reference frame |
| `rg` | radius of gyration — how compact it is |
| `sasa` | solvent-accessible surface area, total, per residue, or each residue's mean |
| `ss` | secondary structure per residue per frame, by DSSP |

**Flexibility**

| | |
|---|---|
| `rmsf` | per-atom or per-residue fluctuation about the mean |
| `dihedrals` | backbone phi, psi and omega, with the Ramachandran plot |

**Conformations**

| | |
|---|---|
| `cluster` | k-means, hierarchical and DBSCAN over the trajectory |
| `dimred` | PCA, MDS, t-SNE and UMAP projections |

**Folding**

| | |
|---|---|
| `hbonds` | hydrogen bonds per frame, and each bond's occupancy |
| `qvalue` | fraction of native contacts retained |

**The ligand**

| | |
|---|---|
| `ligand_rmsd` | how far the ligand has moved, after aligning on the protein |
| `ligand_rmsf` | which parts of the ligand move |

**Protein and ligand together**

| | |
|---|---|
| `pl_contacts` | how much of the protein the ligand touches, with a per-residue fingerprint |
| `pl_hbonds` | hydrogen bonds between them |
| `pl_interactions` | what holds the ligand: eight interaction types, each against a published criterion |

The protein-ligand four run automatically when a ligand is present. See
[Protein-ligand interactions](interactions.md) for what `pl_interactions`
measures and why some of it is refused.

### What the measures actually compute

Details that change the number, and that are worth knowing before comparing
against another tool:

- **RMSD** superposes each frame on the reference before measuring, unless
  `align: false`.
- **Radius of gyration** is mass-weighted by default, which is the physical
  definition; `mass_weighted: false` gives the geometric one.
- **Hydrogen bonds** use Baker-Hubbard at 0.25 nm and 120° by default. Both
  are settings, because published criteria disagree — PLIP allows 4.1 Å and
  100°. Bonds are counted in every frame they exist, including transient ones.
- **Secondary structure** uses MDTraj's DSSP and excludes anything DSSP cannot
  assign, so a ligand does not appear as coil.
- **Q-value** uses a switching function rather than a hard cutoff, with β and
  λ exposed.
- **Clustering** seeds k-means at 42 by default. A clustering that survives a
  change of seed is a finding; one that does not is an artefact of where the
  algorithm started, and the seed is a setting so that can be tested.
- **Contacts and hydrogen bonds** measure across the periodic boundary where
  the trajectory carries a unit cell.
- **Secondary structure** uses MDTraj's DSSP only. Version 1 could also shell
  out to an external `mkdssp`; that is not offered here, because a system
  package is a poor dependency for something that already works — and where
  the two disagree, that is a finding about DSSP worth reporting rather than
  configuring around.

**Writes** `analysis/<name>/` per analysis, each with `<name>.dat`,
`<name>.png`, `<name>.svg` and `options.json`.

---

## report

Assembles everything into things you can send.

The report opens with a **methods paragraph** rather than a list of settings.
A methods section for a molecular dynamics study has to state a particular
list of things — coordinates and their source, protonation, force field
version, water model, box, ions, ensembles, integrator, thermostat, barostat,
cutoffs, constraints, durations — and that list is published, in JCIM's
reporting guidelines and Communications Biology's reproducibility checklist.
Every value is already recorded, so the paragraph is assembled from them.
Nothing is invented, and anything the run did not record is named as missing
rather than filled in with what is usual.

It then reports **convergence**, which is a statement about how much
independent information the trajectory holds. A frame is not an observation:
consecutive frames are nearly the same structure, so the number of independent
observations depends on how quickly a measure forgets where it was, not on how
often frames were written. On a five-thousand-frame trajectory of Trp-cage the
RMSD holds about twenty independent observations, so its uncertainty is
fourteen times what counting frames would give. Where a run is too short to
say — and a short one usually is — the section says what it cannot support.

**Writes** `report.md`, `report.pdf`, `dashboard.html`, `slides.pptx`,
`project_bundle.zip`. The PDF needs the `pdf` extra; where it is absent the
run says so and writes the rest.

---

## What a run leaves behind

```
runs/study/
├── setup/          prepared and solvated structures, the system, what was decided
├── simulation/     the trajectory, the energy log, the settings used
├── analysis/       one directory per measure
├── report/         the written report, slides, dashboard, bundle
├── config.yml      the configuration this run used
└── manifest.json   every phase, artifact and parameter
```

`manifest.json` is the record of the whole run. `config.yml` is what you would
give `fastmdx explore --config` to do it again.
