# Configuration

Everything a run does can go in one YAML file, and the same file drives the
command line, the Python API and the [GUI](gui.md).

```bash
fastmdx init-config study.yml     # a commented template
fastmdx explore --config study.yml
```

---

## The shape of it

```yaml
output: runs/study                # where everything goes

systems:                          # one or more
  - system: 1UBQ
    label: wild_type

include: [setup, simulation, analysis, report]   # which phases

setup:
  ph: 7.4
  forcefield: amber-openff

simulation:
  duration_ns: 100
  temperature_K: 310

analysis:
  include: [rmsd, rmsf, rg, ss]

report:
  author: A. Researcher
```

Every block is optional. What you leave out takes its default, and the run
records what that default was.

---

## Where the settings come from

There is no list of options to keep in step with the software, because the
options **are** the software's declaration. The same file that gives the
command line its flags gives the config file its keys and the GUI its form
fields.

Three ways to see the current set, all of them generated:

```bash
fastmdx init-config study.yml     # a template with every setting, commented
fastmdx explore --help            # the same settings as flags
fastmdx gui                       # the same settings as a form, with explanations
```

The template is the most useful of the three when writing a config: it carries
every key at its default with a sentence about what it does.

---

## Checking a config without running it

```bash
fastmdx explore --config study.yml --dry-run
```

Validates the syntax and every setting, and stops. The GUI does the same
thing from the *A config I already have* option, which will also open the file
for editing without rewriting it.

---

## Per-analysis settings

An analysis's own settings go in a block under its name:

```yaml
analysis:
  include: [hbonds, cluster, sasa]
  options:
    hbonds:
      distance_cutoff: 0.30
      sidechain_only: true
    cluster:
      methods: [kmeans, hierarchical]
      n_clusters: 8
      random_state: 7
    sasa:
      mode: average_residue
```

What each accepts is in [The four phases](phases.md), and the GUI shows it
beside each analysis.

---

## Several systems

```yaml
output: runs/campaign

systems:
  - system: 1UBQ
    label: wild_type
  - system: mutant.pdb
    label: L50A

simulation:
  duration_ns: 100

execution:
  mode: parallel          # or sequential, the default
  workers: 2
  devices: [0, 1]         # one run pinned per GPU
  continue_on_error: true
```

Settings at the top apply to every system; a `setup:` or `simulation:` block
inside a system overrides them for that one.

---

## Reproducing a run

Every run writes the config it used to `config.yml` in its output directory:

```bash
fastmdx explore --config runs/original/config.yml --output runs/repeat
```

For a trajectory that repeats exactly, fix the seed — without one, velocities
differ between runs:

```yaml
simulation:
  random_seed: 42
```

The report's methods section states whether a seed was fixed, because its
absence is what makes a run irreproducible.
