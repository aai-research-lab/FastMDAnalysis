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

**Water**

| | |
|---|---|
| `water_sites` | positions a water holds through the run, and whether one molecule stays or many pass through |

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
| `pmf` | the free energy along an umbrella study's coordinate, drawn from the windows it stitched. Runs where such a study produced one, and reads its result rather than recomputing it |

`water_sites` finds the waters that are part of a binding site rather than
passing through — a water wedged between a ligand and a backbone carbonyl,
bridging a hydrogen bond neither could make alone. It reports how often a
position is occupied *and by how many distinct molecules*, because a site held
by one water throughout and a position a hundred waters pass through can share
an occupancy and mean entirely different things: the first is a molecule to
displace, the second is geometry the protein favours. It needs an explicitly
solvated run.

**Point it at a site, not at a whole protein.** Clustering links neighbours
through neighbours, so a whole-protein scope chains the entire first hydration
shell into one object — on ubiquitin that is tens of thousands of positions
and hundreds of distinct waters, which is a surface rather than a site.
FastMDXplora rejects such a cluster and says so, but the fix is a narrower
`site_selection`: a ligand, a pocket, or a handful of residues. That is a
limit of the method rather than a threshold to tune.

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

### Highlighting regions you care about

A per-residue figure with two hundred residues on the x-axis says little about
the eight that matter. Name them and they are shaded on the RMSF trace and
coloured on a cartoon of the structure:

```yaml
report:
  region_highlights:
    - label: "binding loop"
      start: 84
      end: 92
      color: "#4E79A7"
    - label: "catalytic helix"
      start: 118
      end: 131
```

`label` and `color` are optional — an unlabelled region becomes "Region 1",
an uncoloured one takes the next of six palette colours.

This produces `analysis/rmsf/rmsf_region_highlights.png` and, where PyMOL is
installed (`conda install -c conda-forge pymol-open-source`),
`report/structure_region_highlights.png` with the `.pml` script beside it so
the rendering can be adjusted by hand. Without PyMOL the RMSF figure is still
written and the report records that the structure rendering was skipped and
why.

Regions attach to **RMSF** and nothing else, because RMSF is indexed by
residue — RMSD is indexed by frame, so a residue range has no meaning on it.
RMSF analysis therefore has to have run.

A range outside the residues RMSF measured is refused, with both ranges named:
the one you asked for and the one that exists. That is usually an off-by-one
between a paper's numbering and the structure's, and seeing both makes it
obvious.

The labels are yours. FastMDXplora does not work out that residues 84 to 92
are a binding loop; it draws what you tell it to and calls it what you call
it.

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

Where FastMDXplora was run from a source checkout rather than an install, the
manifest also records the **commit**, and whether the working tree had
uncommitted changes. The version string alone is not enough there: it is
written when the package is installed, so an editable install carries whatever
it was at that moment. A study of ours came back stamped `2.3.0` for a run that
used a feature `2.3.0` did not have. The commit says what the version cannot,
and the dirty flag says when the commit does not describe the code either.

An installed copy records nothing under `source` — there is no checkout to
ask, and the version is the whole answer because the distribution was built
from a tag.
