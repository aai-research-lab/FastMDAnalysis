# The Config

A FastMDXplora Config is the whole description of a study: the system, how it
is prepared, how it is simulated, what is measured, and how it is written up.
Capture that and the four phases run themselves.

It is written as YAML, which is the format rather than the thing. The same
Config drives the [GUI](gui.md), the [command line](cli_reference.md) and the
[Python API](api.md), and each of them can produce one.

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
    id: wild_type                 # optional; s1, s2... if omitted

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
    id: wild_type
  - system: mutant.pdb
    id: L50A
    simulation:                   # this system only
      temperature_K: 320

simulation:
  duration_ns: 100

execution:
  mode: parallel          # or sequential, the default
  workers: 2
  devices: [0, 1]         # one run pinned per GPU
  continue_on_error: true
```

Each entry needs a `system`. `id` is optional and names the run — it becomes
the run's directory and its label in the comparison report; without one the
systems are `s1`, `s2` and so on, which is harder to read six months later.
Two systems cannot share an `id`.

A phase block inside an entry overrides the top-level one for that system
only, which is how you vary a single setting across a campaign.

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
