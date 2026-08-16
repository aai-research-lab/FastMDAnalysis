# Reading the results

A finished run is a directory, not a number. This page is the map: what is
in it, what to open first, and how to tell a measurement from something the
run could not support.

## What a run leaves behind

```
runs/study/
├── setup/          prepared and solvated structures, the system, what was decided
├── simulation/     the trajectory, the energy log, the settings used
├── analysis/       one directory per measure
├── report/         the written report, slides, dashboard, bundle
├── resolved_config.yml   every setting filled in, defaults included
└── manifest.json         every phase, artifact and parameter
```

Three files answer three different questions, and they are worth telling
apart:

| | |
|---|---|
| `resolved_config.yml` | **What was asked for**, with every default filled in. Give it back to `fastmdx explore --config` to do the run again. |
| `manifest.json` | **What happened**: every phase, artifact and parameter, plus the input structure's SHA-256, the PDB identifier and deposition date read from the file's own header, and when it was fetched. A report saying "PDB 3PTB" is checkable against the exact bytes that entered the run. |
| `analysis/<n>/options.json` | **What one measure did**: its selection, every option, the findings, and the format of the file beside it. |

**The trajectory holds the solute, not the box.**
`simulation.save_selection` defaults to `not water`, because water is nine
tenths of a solvated system and nine tenths of the file: a 20 ns run of a
small protein is about 740 MB with it and under 80 without. What is saved
comes with `trajectory_topology.pdb`, the topology that matches it, and the
whole system stays in `setup/`. Set it to `all` where the solvent is the
subject — the water-site analysis needs it and says so rather than
reporting an empty result.

Where FastMDXplora ran from a source checkout rather than an install, the
manifest also records the **commit** and whether the working tree was dirty.
The version string alone is not enough there: it is written when the package
is installed, so an editable install carries whatever it was at that moment.
A study of ours came back stamped `2.3.0` for a run that used a feature
`2.3.0` did not have. An installed copy records nothing under `source`,
because the distribution was built from a tag and the version is the whole
answer.

## Open these first

**`report/report.pdf`** is the study written up: a methods paragraph you can
paste into a manuscript, the figures, and a section saying what the run does
and does not support.

**`analysis/rmsd/rmsd.png`** answers the first question about any trajectory:
has the structure settled, or is it still moving?

**`setup/setup_parameters.json`** says what was decided about your structure
before any dynamics ran — every non-standard residue, every protonation
call. If a result surprises you, the explanation is usually here.

Or open all of it at once:

```bash
fastmdx gui --output runs/study
```

## How a number tells you what it is worth

Every measure reports through the same four registers, and reading them is
most of reading a result.

**A number, plainly.** A settled mean arrives with a standard error and the
number of independent observations behind it, not the number of frames.
Saving frames more often makes a file larger without making a measurement
better, so the effective sample count is the figure to read.

**A number, with a statement of what it does not support.** A finding
carrying `not_a_measurement` is still reported, because the value is often
the best available, and the note says what is wrong with it. The commonest
case: a run too short against its own correlation time cannot measure how
correlated it is, so the effective-sample count is an upper bound and the
true figure is smaller.

**A range instead of a number**, where the run recorded enough to bound a
quantity but not to pin it.

**No number at all**, with the reason. A free energy from a bias that never
settled, a potential of mean force across windows that do not overlap, a
binding free energy from a run that never reached bulk: each is refused by
name rather than drawn. A refusal is the most informative thing this
software produces, and it always says what would settle the question.

Averages taken on a biased run are corrected back to equilibrium where the
bias allows, and labelled as biased where it does not
([Averages on a biased run](phases.md#averages-on-a-biased-run)).

## Where each kind of result is explained

The detail that changes a number lives with the measure rather than here, so
that there is one place to correct when it changes:

- **What each analysis computes**, and the choices that move it — cutoffs,
  switching functions, alignment sets, atom typing:
  [the analysis phase](phases.md#analysis).
- **Protein--ligand interactions**, which carry the most criteria and the
  most ways to disagree with another tool:
  [Protein-ligand interactions](interactions.md).
- **Free energy** from umbrella sampling, metadynamics or steered pulling,
  and the gates each output has to pass:
  [Simulations beyond a box of water](simulations.md).
- **The file formats themselves.** Each analysis writes its data as `.dat`
  and records in its `options.json` whether that is whitespace with no
  header or comma-separated with one, together with the one-liner that
  reads it. Read the record rather than guessing from the extension.

## Comparing runs

A campaign writes `batch_manifest.json` beside the members, a `comparison/`
report overlaying their series, and `members.json`, which collects what each
member concluded.

`members.json` distinguishes two things that look identical on disk. Members
differing only by random seed are repeats of one measurement, so the spread
of their means is the **error** on it, and it is set against the error each
run estimated for itself; where the replicas spread much wider, the
single-run estimate was reading a correlation time off a trajectory too
short to contain it. Members differing by system, mutation or parameter are
different measurements, so the spread between them is the **result**. The
file says which it decided and why.
