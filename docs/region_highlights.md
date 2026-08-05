# Highlighting regions in a figure

A per-residue figure with forty residues on the x-axis says little about the
three you care about. Region highlights mark them: a loop, a helix, the
residues around an active site.

```yaml
report:
  region_highlights:
    - label: "binding loop"
      start: 84
      end: 92
      color: "#4E79A7"
```

They apply to **RMSF**, which is indexed by residue. RMSD is indexed by frame,
so a residue range means nothing on it.

The labels are yours. FastMDXplora does not work out that a range is a loop or
a binding site — it draws what you tell it to draw, which is the point: you
know which residues matter and the software does not.

## YAML example

```yaml
include: [analysis, report]

analysis:
  include: [rmsf]
  trajectory: simulation/production.dcd
  topology: simulation/topology.pdb

report:
  title: "Trp-cage RMSF region highlights"
  region_highlights:
    - label: "example region 1"
      start: 3
      end: 7
      color: "#4E79A7"
    - label: "example region 2"
      start: 12
      end: 16
      color: "#F28E2B"
```

The labels are user-provided annotations. FastMDXplora does not infer that a
range is biologically a loop, helix, catalytic motif, or binding site.

## Outputs

When `report.region_highlights` is configured and RMSF analysis output exists,
the report phase writes:

- `analysis/rmsf/rmsf_region_highlights.png`
- `report/structure_region_highlights.png` when PyMOL rendering succeeds
- `report/structure_region_highlights.pml` when PyMOL rendering succeeds
- `report/region_highlight_summary.png`
- `report/region_highlight_manifest.json`

If PyMOL and a topology/PDB file are available, the summary includes a
PyMOL-rendered cartoon/ribbon structure panel with the same regions
highlighted. If PyMOL is unavailable, FastMDXplora still writes the RMSF
highlight plot and records the skipped structure-rendering reason in
`report/region_highlight_manifest.json`.

Install PyMOL into the environment that runs FastMDXplora when you want the
cartoon structure panel:

```bash
conda install -c conda-forge pymol-open-source
```

With micromamba:

```bash
micromamba install -c conda-forge pymol-open-source
```

## Run examples

Bash:

```bash
fastmdx explore --config region_highlights.yml
```

Windows PowerShell fallback when `fastmdx` is not on PATH:

```powershell
python -m fastmdxplora.cli.main explore --config region_highlights.yml
```

For an analysis-only/report-only workflow, provide existing trajectory and
topology paths in the `analysis:` block and set `include: [analysis, report]`.
Setup and simulation will not be run.
