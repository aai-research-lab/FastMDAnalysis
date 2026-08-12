# Production runs and GPUs

Everything on the other pages works on a laptop. This one is about the runs
that do not.

---

## Where a production run belongs

A real trajectory is hours to days of continuous computation. Two things
follow.

**Run it on a GPU.** OpenMM on a modern GPU is one to two orders of magnitude
faster than on a CPU. `--simulate-platform CUDA` asks for one; the default is
`auto`, which takes the fastest available and says which it chose.

**Run it somewhere it will not be interrupted.** A closed laptop lid or a
disconnected SSH session ends the run. Use a scheduler, or `nohup`, or `tmux`.

---

## Before committing hours to it

```bash
fastmdx info
```

Confirms the simulation backends are present.

```bash
fastmdx explore --system your.pdb --output runs/smoke \
  --simulate-nvt-steps 500 \
  --simulate-npt-steps 500 \
  --simulate-production-steps 5000
```

Ten picoseconds through every phase. It proves the structure prepares, the
system builds, the platform works and the report writes — which is everything
except the length. A setup problem found after six hours of production is the
same problem found in two minutes.

Then check `setup/setup_parameters.json` and read what the setup phase decided
about your structure. That is where a wrong answer is cheapest to catch.

---

## Launching one

```bash
nohup fastmdx explore \
  --config study.yml \
  --output /scratch/$USER/runs/my_run \
  --simulate-platform CUDA \
  > /scratch/$USER/runs/my_run.log 2>&1 &
```

Under a scheduler, the same command goes in the job script. Put the output on
fast local storage rather than a network filesystem: the trajectory writer
touches the file constantly, and NFS makes that the bottleneck.

---

## Watching it

The log says what is happening. For anything more, point the GUI at the output
directory:

```bash
fastmdx gui --output /scratch/$USER/runs/my_run
```

That works while the run is in progress — telemetry, energies, and the
molecule moving as frames arrive.

From your own machine against a cluster, keep the directory in sync:

```bash
rsync -az --delete user@cluster:/scratch/$USER/runs/my_run/ ~/runs/my_run/
fastmdx gui --output ~/runs/my_run
```

Re-run the `rsync` as often as you like.

---

## Choosing a length

The honest answer is that it depends on what you are asking, and the report
will tell you whether you got there. Its convergence section reports how many
**independent** observations the trajectory holds — which is far fewer than
the frame count, because consecutive frames are nearly the same structure.

A useful pattern: run something, read the convergence section, and extend if
it says the measures have not settled or rest on too little. That is more
reliable than choosing a number in advance.

---

## When it stops early

Every run writes checkpoints. If it dies, the checkpoint is in the simulation
directory and the manifest records how far it got.

If it became unstable rather than being killed, the message says which atoms
went wrong and what that points at — a ligand alone means its parameters, the
whole system at once means the integration, and the remedies differ. Nothing
is retried automatically: a run that exploded because its ligand is wrong will
explode again more slowly at half the timestep.

**The log keeps what happened, including the attempts that failed.**
`fastmdxplora.log` in the run directory is appended to rather than replaced,
so a second attempt does not erase the first -- which is usually the one
saying why the second was needed. Each invocation is separated by a banner
naming the time and the version that made it, so a directory re-run over weeks
reads as a sequence rather than a wall.

That is worth knowing when a directory has been re-run with
`--force-overwrite`: the other artefacts were replaced and the log was not, so
it describes runs whose outputs are no longer there. The banners are what tell
you which part of it is about the files you are looking at.

---

## Several runs at once

```yaml
execution:
  mode: parallel
  workers: 2
  devices: [0, 1]         # one run pinned per GPU
  continue_on_error: true
```

`workers` should not exceed the number of GPUs — two runs sharing one device
are slower than the same two in sequence.

`continue_on_error` is what you want overnight: one system failing does not
end the campaign, and the failures are in the manifest.

---

## See also

- [Worked examples](usage_examples.md) — campaigns and comparisons
- [Beyond a box of water](simulations.md) — restraints, membranes, metadynamics
