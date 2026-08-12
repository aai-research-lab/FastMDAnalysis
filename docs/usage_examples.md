# Worked examples

Recipes for studies people actually run. Each is complete — copy it, change
the structure, and it works.

For every flag see the [CLI reference](cli_reference.md); for every setting,
[Configuration](configuration.md).

---

## A protein on its own

```bash
fastmdx explore --system 1UBQ --output runs/ubiquitin
```

Fetches ubiquitin from the PDB, prepares it, simulates it, analyses it, writes
the report. Everything else on this page is a variation.

To make it a real study rather than a smoke test, say how long:

```bash
fastmdx explore --system 1UBQ --output runs/ubiquitin \
  --simulate-duration-ns 100
```

---

## A protein with a ligand

```bash
fastmdx explore --system 181L \
  --setup-forcefield amber-openff \
  --output runs/lysozyme
```

`amber-openff` is what parameterises the ligand. The setup phase finds it,
looks its chemistry up, settles its protonation in the binding site, and
discards crystallisation additives. The analysis phase adds the
protein-ligand measures automatically.

That route needs the network: the chemistry is looked up by identifier, and a
local file carries no entry to look it up in. Where the lookup is not
available -- an offline machine, or a ligand the dictionary does not have --
supply the chemistry as a file.

**Which leaves two situations, and they want opposite things from the same
pair of files.**

### The ligand is in the structure

A complex has the ligand at its crystallographic coordinates and no chemistry.
The supplied file provides the bond orders, formal charges and aromaticity a
PDB cannot express; its coordinates are usually idealised and mean nothing
here.

**The structure wins.** The ligand is placed where the complex has it, and the
run says so. The ligand name has to match the residue name in the structure,
because that is what identifies which component to take the coordinates from.

Getting this backwards is not a small error. An ideal benzene sat seventeen
Angstroms from the cavity it belonged in: setup succeeded, the clash check
passed -- seventeen Angstroms is not a clash -- and the run was of a benzene
floating in solvent. Nothing looked wrong.

### The ligand is not in the structure

An apo protein, and a pose from docking or built by hand. **The file wins**,
because it is the only thing that knows where the ligand goes. Its coordinates
have to be a bound pose: nothing here can establish that, and a pose that is
wrong will simulate perfectly well.

Nothing has to be declared to choose between the two. If the structure holds a
residue of that name with a matching count of heavy atoms, its coordinates are
used; if it does not, the file's are. The files already say which situation it
is.

In either case, supply the chemistry like this:

```bash
fastmdx explore --system complex.pdb \
  --setup-forcefield amber-openff \
  --setup-ligand ligand.sdf \
  --output runs/complex
```

See [Protein-ligand interactions](interactions.md) for what comes out.

---

## A membrane protein

```bash
fastmdx explore --system 1AFO \
  --setup-forcefield amber14 \
  --setup-membrane POPC --setup-membrane-orient \
  --simulate-restrain "protein and not element H" \
  --output runs/glycophorin
```

`--setup-membrane-orient` turns the structure so its longest axis lies along
the membrane normal, which a PDB entry usually does not. It is checked both
ways — before, that there is a longest axis worth using; after, that the
result has the hydrophobic belt a bilayer-spanning protein has.

The restraint holds the protein while the lipids pack around it. See
[Beyond a box of water](simulations.md).

---

## Holding a structure still while it settles

```bash
fastmdx explore --system 4LYT --output runs/lysozyme \
  --simulate-restrain "protein and not element H" \
  --simulate-nvt-steps 25000 --simulate-npt-steps 25000
```

The restraint is released in stages through equilibration and is off for
production. It matters most where a pose has to survive heating: a ligand, a
membrane, a loop that was modelled in.

---

## Metadynamics

Biasing a torsion, which is a bounded coordinate and needs nothing further:

```yaml
# torsion.yml
output: runs/torsion
systems:
  - system: 1UBQ
simulation:
  duration_ns: 50
  metadynamics:
    collective_variable: torsion
    selection: "resid 25 and name N CA CB CG"
    sigma: 0.35
```

Biasing a ligand's distance from its site, which needs bounding or the ligand
leaves and never comes back:

```yaml
# unbinding.yml
output: runs/unbinding
systems:
  - system: 181L
setup:
  forcefield: amber-openff
simulation:
  duration_ns: 200
  metadynamics:
    collective_variable: ligand_distance
    site_selection: "resid 84 to 121 and name CA"
    sigma: 0.05
    funnel:
      axis_selection: "resid 130 to 134 and name CA"
```

```bash
fastmdx explore --config unbinding.yml
```

---

## Umbrella sampling

A study, not a run: the range is divided into windows, each held at its own
value of the coordinate, and the free energy is stitched from the overlap
between them. It expands into one run per window and goes through the batch
machinery, so the windows can run in parallel.

```yaml
# psi.yml
output: runs/psi
systems:
  - system: peptide.pdb
simulation:
  duration_ns: 5
  umbrella:
    collective_variable: torsion
    selection: "resid 1 and name N CA C or (resid 2 and name N)"
    from: -3.1416          # radians
    to: 2.6180
    n_windows: 12
    force_constant: 50.0
execution:
  mode: parallel
  workers: 4
```

Two things decide whether this produces anything. **The windows have to
overlap**: the width each one samples is `sqrt(kT/k)`, and neighbours further
apart than about twice that share too little for the stitching to join them.
The study reports every overlap and refuses below a threshold rather than
returning a curve nobody should read.

**A torsion is a circle.** Windows covering part of a turn measure one of the
two paths between the states and say nothing about the other — which may be
the higher. The example above tiles the full turn, from -pi to +2.618 radians
at 30-degree spacing, so the last window's neighbour is the first.

```bash
fastmdx explore --config psi.yml
```

The free energy lands in `pmf.json` with the barrier, the minima, the range
the windows covered, and — on a closed turn — how far the profile misses
meeting itself, which is what that study's statistics are worth.

---

## Steered molecular dynamics

Pulling along a coordinate at constant velocity and recording the work. Used
to force a transition that will not happen on its own, and to estimate a free
energy from the work through Jarzynski's equality across repeats.

```yaml
# pull.yml
output: runs/pull
systems:
  - system: 1UBQ
simulation:
  duration_ns: 2
  steered:
    collective_variable: radius_of_gyration
    from: 0.71             # where the coordinate starts
    to: 1.20               # where to pull it
    force_constant: 1000.0
```

`from` has to be the value the coordinate actually has at the start. It is
refused rather than guessed, because a pull that begins somewhere the molecule
is not spends its first half hauling the system to the anchor and reports work
that means nothing.

Measure it before writing the config:

```bash
fastmdx analyze --trajectory equilibrated.dcd --topology system.pdb --include rg
```

```bash
fastmdx explore --config pull.yml
```

The work profile lands in the analysis directory. On a correct pull it starts
at zero and increases; work that goes negative early is the anchor being in
the wrong place.

---

## Analysing a trajectory you already have

No simulation, from any engine that writes a format MDTraj reads:

```bash
fastmdx analyze \
  --trajectory production.dcd \
  --topology system.pdb \
  --output runs/analysis
```

Or as a config, which is easier to keep:

```yaml
output: runs/analysis
include: [analysis, report]
systems:
  - system: system.pdb
analysis:
  trajectory: production.dcd
  topology: system.pdb
  include: [rmsd, rmsf, rg, ss, sasa]
```

---

## Choosing what gets measured

```bash
fastmdx explore --system 1UBQ --output runs/focused \
  --analyze-analyses rmsd rmsf ss
```

Or the other way round — everything except the slow ones:

```bash
fastmdx explore --system 1UBQ --output runs/most \
  --analyze-exclude-analyses sasa dimred
```

The common per-analysis settings have flags of their own:

```bash
fastmdx explore --system 1UBQ --output runs/tuned \
  --analyze-cluster-n-clusters 8 \
  --analyze-cluster-methods kmeans
```

Anything else goes in a config file, where an analysis's settings are a block
under its name:

```yaml
analysis:
  include: [hbonds, cluster]
  options:
    hbonds:
      distance_cutoff: 0.30
      angle_cutoff: 150.0
    cluster:
      n_clusters: 8
```

---

## Several systems at once

```yaml
# campaign.yml
output: runs/campaign
systems:
  - system: 1UBQ
    id: wild_type
  - system: mutant.pdb
    id: L50A
simulation:
  duration_ns: 100
execution:
  mode: parallel
  workers: 2
  devices: [0, 1]        # one run pinned per GPU
  continue_on_error: true
```

```bash
fastmdx explore --config campaign.yml
```

`continue_on_error` keeps the campaign going when one system fails, which is
what you want overnight — the failures are in the manifest.

---

## Comparing runs

A study with more than one system writes a cross-run comparison as well as the
individual reports — the same measures overlaid, so a difference between
systems is visible rather than inferred from reading two reports side by side.

It is on by default for a multi-system study; `report.comparison: false` turns
it off. The campaign above produces one under `comparison/`.

---

## Reproducing a run exactly

Every run writes the config it used. To repeat one:

```bash
fastmdx explore --config runs/original/config.yml --output runs/repeat
```

For bit-for-bit reproducibility, fix the seed — without one, velocities differ
between runs and so does the trajectory:

```bash
fastmdx explore --system 1UBQ --output runs/seeded --simulate-random-seed 42
```

The report's methods section says whether a seed was fixed, because its
absence is what makes a run irreproducible.

---

## From Python

```python
import fastmdxplora as fastmdx

runs = fastmdx.FastMDXplora(
    system="181L",
    output_dir="runs/lysozyme",
    setup={"forcefield": "amber-openff"},
    simulation={"duration_ns": 100},
).explore()

print(runs[0].output_dir)
```

A config file works here too:

```python
runs = fastmdx.FastMDXplora(config="study.yml").explore()
```

See the [API reference](api.md).
