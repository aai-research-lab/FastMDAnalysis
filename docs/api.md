# Python API

Everything the command line does is available from Python, through one class.

```python
import fastmdxplora as fastmdx

runs = fastmdx.FastMDXplora(system="1UBQ", output_dir="runs/ubiquitin").explore()
print(runs[0].output_dir)
```

`explore()` runs all four phases and returns one result per system.

---

## Constructing it

```python
fastmdx.FastMDXplora(
    system="1UBQ",              # a PDB identifier or a path
    output_dir="runs/study",
    options={                   # settings, by phase
        "setup": {"ph": 7.0, "forcefield": "amber-openff"},
        "simulation": {"duration_ns": 100},
        "analysis": {"include": ["rmsd", "rmsf", "ss"]},
    },
)
```

Or from a config file, which is the same file the CLI and the GUI use:

```python
fastmdx.FastMDXplora(config="study.yml")
```

The settings under `options` are exactly the config file's blocks, so anything
[Configuration](configuration.md) describes works here.

---

## Running part of it

```python
study = fastmdx.FastMDXplora(system="1UBQ", output_dir="runs/study")

study.explore(include=["setup", "simulation"])   # stop after the trajectory
study.explore(exclude=["report"])                # everything but the write-up
```

---

## What comes back

`explore()` returns a `RunResult` per system, carrying where the output went,
which phases ran, and what each produced. The same information is in
`manifest.json` in the output directory, which is what to read if the process
has ended.

```python
for run in runs:
    print(run.output_dir)
    for phase in run.phases:
        print(" ", phase.name, phase.status)
```

---

## Watching a run from Python

There is no separate Python dashboard. Serve the output directory with the
GUI from another terminal, while the Python process runs:

```bash
fastmdx gui --output runs/study
```

That gives telemetry, the molecule in 3D, and the results as they appear —
the same as for a run started from the command line.

---

## See also

- [Worked examples](usage_examples.md) — recipes, several with Python
- [Configuration](configuration.md) — every setting `options` accepts
