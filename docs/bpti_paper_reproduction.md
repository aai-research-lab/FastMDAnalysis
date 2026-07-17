# Reproducing the BPTI paper case study

FastMDXplora provides a `bpti-paper` project preset and the fully expanded
example at `examples/bpti_paper_6pti.yaml`. Both target the BPTI case study
in *FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics
Trajectories* (Aina & Kwan, JCC 2026).

## Use 6PTI, not 1K1N

Use PDB **6PTI**, the BPTI protein structure used for the paper case study.
Do **not** substitute **1K1N**: it is a trypsin-inhibitor complex, so it has
the enzyme as well as inhibitor/complex context. That changes the analyzed
atom population and produces different radius of gyration, SASA, hydrogen
bond, and RMSF values. The profile removes non-water heteroatoms and
crystallographic waters before protonation, leaving the intended BPTI system.

## Prerequisite: a C36m-capable OpenMM installation

The profile maps `forcefield: charmm36m` only to OpenMM's
`charmm36_2024.xml` plus `charmm36_2024/water.xml`. That XML set contains
CHARMM36m protein parameters and the compatible CHARMM-modified TIP3P water
and ions. It is not an alias for legacy `charmm36.xml`.

```bash
mamba env create -f environment.yml
conda activate fastmdxplora
pip install -e .
```

If those exact XML files are unavailable, setup stops with a message naming
the missing files and does not substitute CHARMM36. Upgrade OpenMM (8.3 or
newer is recommended) instead of continuing with a different force field.

## HPC / CUDA run

The example requests CUDA mixed precision and keeps live telemetry off.

```bash
fastmdx explore --config examples/bpti_paper_6pti.yaml --output ./bpti_6pti
```

The built-in profile is equivalent for a fresh run:

```bash
fastmdx explore --preset bpti-paper --output ./bpti_6pti
```

For separate, restart-friendly phase invocations, retain the same output
directory and config file:

```bash
fastmdx explore --config examples/bpti_paper_6pti.yaml --output ./bpti_6pti --include setup
fastmdx explore --config examples/bpti_paper_6pti.yaml --output ./bpti_6pti --include simulation
fastmdx explore --config examples/bpti_paper_6pti.yaml --output ./bpti_6pti --include analysis
fastmdx explore --config examples/bpti_paper_6pti.yaml --output ./bpti_6pti --include report
```

Production is 50,000,000 steps at 2 fs (100 ns); frames are written every
5,000 steps (10 ps), and analysis loads every second frame. The generated
`resolved_config.yml` and phase manifests record the effective settings.

## Profile settings

- Setup: pH 7.0; no retained heteroatoms or waters; CHARMM36m with
  CHARMM-compatible TIP3P; cubic 1.0 nm padding; 0.15 M NaCl; PME; HBonds;
  1.5 amu hydrogen-mass repartitioning.
- Simulation: 300 K Langevin middle, friction 1/ps, 250,000 NVT steps,
  1,000,000 NPT steps, barostat interval 50, CUDA mixed precision, and a
  checkpoint every 100,000 steps.
- Analysis: protein scope, stride 2, aligned backbone RMSD to frame 0,
  unaligned all-protein atomic RMSF, total SASA with a 0.14 nm probe,
  hierarchical clustering with `k=6`, and PCA only.

`rmsf.align` is `true` by default globally. It is explicitly `false` here,
matching the paper-era direct fluctuation-around-the-mean calculation.

## Nonbonded switching limitation

The paper reports a 1.0 nm PME real-space electrostatic cutoff and
Lennard-Jones switching from 1.0 to 1.2 nm. The OpenMM `ForceField`
`createSystem` path used by FastMDXplora exposes one shared nonbonded cutoff;
it cannot represent that mixed electrostatic/LJ cutoff scheme exactly. A
switch distance equal to the shared 1.0 nm cutoff would be a zero-width,
invalid switching region, so the runnable profile preserves the 1.0 nm
nonbonded cutoff and disables switching rather than claiming the paper's
1.0-1.2 nm treatment. Use a custom OpenMM system builder if exact split-cutoff
reproduction is required.
