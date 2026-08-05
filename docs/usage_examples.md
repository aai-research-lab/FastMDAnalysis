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

If the structure's ligand is not in the Chemical Component Dictionary, supply
its chemistry:

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
