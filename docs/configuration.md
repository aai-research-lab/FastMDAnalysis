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
  ph: 7.0
  forcefield: charmm36       # charmm36 | amber14 | amber-fb15 | amber-openff
  solvent_padding_nm: 1.0
  ion_concentration_M: 0.15

simulation:
  preset: gentle             # optional conservative smoke-test preset
  duration_ns: 10            # or set nvt/npt/production steps explicitly
  integrator: langevin_middle
  temperature_K: 300.0
  pressure_bar: 1.0
  checkpoint_interval_steps: auto  # every 20% of production; 0 disables
  first_aid_max_retries: 3         # retry failed MD stages from checkpoint
  restraint_type: position         # optional; currently position only
  restraint_selection: backbone    # protein | backbone | heavy | non-water | all
  restraint_k: 1000.0              # kJ mol^-1 nm^-2
  restraints_in_production: false  # remove after NVT/NPT by default
  plumed:                    # optional enhanced sampling
    enabled: false
    script: bias.dat

analysis:
  scope: solute              # solute | protein | ligand | all
  include: [rmsd, rmsf, rg]  # omit to run the full suite

report:
  title: "My study"
  slides: true
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
  `pressure_bar`, `platform`, `precision`, `checkpoint_interval_steps`,
  `first_aid_max_retries`, `restraint_type`, `restraint_selection`,
  `restraint_k`, `restraints_in_production`, `plumed`.
- **analysis**: `scope`, `selection`, `include`, `exclude`, `stride`,
  `first`, `last`, `options`.
- **report**: `title`, `author`, `document`, `slides`, `bundle`,
  `include_methods`, `include_reproducibility`, `comparison`.
