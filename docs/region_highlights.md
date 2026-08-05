# Highlighting regions

A per-residue figure with two hundred residues on the x-axis says little about
the eight you care about. Region highlights mark them — on the RMSF plot, and
on the structure itself.

```yaml
include: [analysis, report]

analysis:
  include: [rmsf]

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

`label` and `color` are optional: an unlabelled region becomes "Region 1", and
an uncoloured one takes the next colour from a palette chosen so adjacent
regions stay distinguishable.

---

## What it produces

**`analysis/rmsf/rmsf_region_highlights.png`** — the RMSF trace with each
region shaded and labelled.

**`report/structure_region_highlights.png`** — the regions drawn on the
structure as a cartoon, each in its colour. This one needs PyMOL:

```bash
conda install -c conda-forge pymol-open-source
```

Without it, the RMSF figure is still produced and the report records that the
structure rendering was skipped and why, rather than failing the phase. The
PyMOL script is written alongside as
**`report/structure_region_highlights.pml`**, so the rendering can be
reproduced or adjusted by hand.

---

## Why RMSF

Region highlights attach to **RMSF** and nothing else, because RMSF is indexed
by residue. RMSD is indexed by frame — it says how far the whole structure has
moved at each point in time — so a residue range has no meaning on it.

RMSF analysis therefore has to have run. Without `analysis/rmsf/rmsf.dat` the
report phase says so and names what would produce it, rather than drawing an
empty figure.

---

## What is checked

Regions are validated against the residues RMSF actually measured, so a
mistake is caught before anything is drawn:

- `start` must be at least 1
- `end` must be at least `start`
- the range must lie inside the residues present in the RMSF output

A range outside the structure is refused with both numbers named — the range
you asked for and the range that exists. That is usually an off-by-one between
a paper's numbering and the structure's, and seeing both makes it obvious.

---

## The labels are yours

FastMDXplora does not work out that residues 84 to 92 are a binding loop. It
draws what you tell it to draw and labels it what you call it, which is the
point: you know which residues matter, and the software does not.
