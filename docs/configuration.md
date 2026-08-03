# Configuration

For anything beyond a quick run, capture the whole study in a single YAML file
instead of a long list of flags. The same file drives both the CLI
(`fastmdx explore --config study.yml`) and the Python API
(`FastMDXplora(config="study.yml")`).

Input is always given as a `systems:` list, even for a single system, so the
file looks the same whether you study one protein or a dozen.

## Structure

A config has a top level plus one block per phase (`setup`, `simulation`,
`analysis`, `report`). Every key is optional; omitted keys fall back to
sensible defaults.

```yaml
systems:
  - id: trpcage
    system: 1L2Y             # PDB ID, file path, or one-letter sequence

setup:
  ph: 7.4                    # physiological; the default
  forcefield: auto           # auto | charmm36 | amber14 | amber-fb15 | amber-openff
  heterogens: auto           # auto | drop | keep
  solvent_padding_nm: 1.0
  ion_concentration_M: 0.15

simulation:
  preset: gentle             # optional conservative smoke-test preset
  duration_ns: 10            # or set nvt/npt/production steps explicitly
  integrator: langevin_middle
  temperature_K: 300.0
  pressure_bar: 1.0
  plumed:                    # optional enhanced sampling
    enabled: false
    script: bias.dat

analysis:
  scope: solute              # solute | protein | ligand | all
  include: [rmsd, rmsf, rg]  # omit to run the full suite
  # For the published version 1 comparison profile, use the CLI flag
  # --analyze-compat v1 or direct `analyze --compat ...`.

report:
  title: "My study"
  slides: true
  region_highlights:
    - label: "example region 1"
      start: 3
      end: 7
      color: "#4E79A7"
```

## Protein-ligand studies

Supply a ligand in the `setup` block. When a ligand is present, the
ligand-aware analyses run automatically.

```yaml
setup:
  forcefield: amber-openff
  ligand: ligand.sdf
  ligand_name: LIG
  check_ligand_clashes: true
```

## Reproducing a study

Every run writes `resolved_config.yml`, the fully merged configuration that
actually ran (defaults plus your file plus any command-line overrides). Feed
it straight back to `--config` to reproduce the study exactly.

## Full field reference

The authoritative list of every option, its type, and its default lives in the
schema module, `fastmdxplora.config.schema`. The blocks and their most common
keys:

- **setup**: `ph`, `forcefield`, `ligand`, `ligand_name`,
  `solvent_padding_nm`, `box_shape`, `ion_concentration_M`,
  `nonbonded_cutoff_nm`, `constraints`, `temperature_K`.
- **simulation**: `duration_ns` (or `nvt_steps` / `npt_steps` /
  `production_steps`), `integrator`, `timestep_fs`, `temperature_K`,
  `pressure_bar`, `platform`, `precision`, `plumed`.
- **analysis**: `scope`, `selection`, `include`, `exclude`, `stride`,
  `first`, `last`, `options`. The `v1` compatibility profile is a
  CLI comparison preset rather than a YAML field; use `--compat` for direct
  `analyze` or `--analyze-compat` under `explore`.
- **report**: `title`, `author`, `document`, `slides`, `bundle`,
  `include_methods`, `include_reproducibility`, `region_highlights`,
  `comparison`.

See [Region highlight figures](region_highlights.md) for RMSF-specific
residue-region highlight examples and output paths.

`fastmdx report` also writes a static dashboard at
`<output>/report/dashboard.html` using only metrics and plots present in
the run outputs.
