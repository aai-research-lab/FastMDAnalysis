# Running somewhere else

A laptop is where a study gets designed. It is rarely where it should run: a
solvated protein of twenty thousand atoms manages perhaps ten nanoseconds a
day on a CPU and several hundred on a current GPU, which is the difference
between a result tomorrow and a result before lunch.

What travels between the two is the config file. A study designed in the GUI
runs unchanged from the command line on a cluster, because the config is the
whole description of it.

---

## The shape of it

```bash
# design here
fastmdx gui                       # or write the YAML by hand

# run there
scp study.yml protein.pdb cluster:~/work/
ssh cluster
fastmdx explore --config study.yml --output runs/study

# read here
rsync -az cluster:~/work/runs/study/ ~/runs/study/
fastmdx gui --output ~/runs/study
```

The GUI reads a finished run as happily as it drives a live one, so the
results come back to the machine where they are comfortable to look at.

---

## Installing where there is no network

Many clusters have no route to PyPI or conda-forge. FastMDXplora is pure
Python; what is awkward to install is the layer underneath it — OpenMM,
PLUMED, MDTraj, PDBFixer — and that layer is often already there, because it
is what everybody else on the machine needs too.

**Look before building anything.**

```bash
module avail 2>&1 | grep -iE "openmm|conda|cuda"
ls /opt/conda_envs /shared/envs 2>/dev/null
```

A site environment with OpenMM and PLUMED in it turns this from a container
build into a five-minute install.

**Then find what is missing.**

```bash
python -c "
for m in ('openmm','mdtraj','pdbfixer','sklearn','pandas','matplotlib',
          'scipy','yaml','pptx','PIL','numpy'):
    try:
        __import__(m); print(f'  {m:12} present')
    except ImportError:
        print(f'  {m:12} MISSING')"
python --version
```

**Fetch only those**, on a machine with a network, built for the target — not
for your laptop, which is very likely a different architecture:

```bash
pip download fastmdxplora python-pptx --no-deps \
  --platform manylinux2014_x86_64 --python-version 39 \
  --only-binary=:all: -d wheels
tar czf wheels.tgz wheels
```

`--no-deps` matters. Without it pip re-resolves the whole tree and fails on
packages that have no wheel — MDTraj among them — even where the requirement
is already satisfied on the target.

**Install on top of the site environment**, so OpenMM and PLUMED are taken
from it rather than replaced:

```bash
tar xzf wheels.tgz
pip install --user --no-index --find-links=wheels fastmdxplora python-pptx
```

Where `--user` is not allowed, a virtual environment layered over the site
one does the same job:

```bash
python -m venv --system-site-packages ~/fastmdx_env
source ~/fastmdx_env/bin/activate
pip install --no-index --find-links=wheels fastmdxplora python-pptx
```

`--system-site-packages` is what keeps the expensive half.

### Where wheels are not enough

The route above installs FastMDXplora and its pure-Python dependencies. It
does not reach the ligand path: **the OpenFF toolkit and openmmforcefields are
distributed through conda-forge only**, are not on PyPI, and no arrangement of
wheels will fetch them. A machine without them prepares proteins perfectly
well and refuses the moment a ligand needs parameterising.

There are two ways round it.

**Install from conda-forge into the environment you have.** Where the machine
can reach conda-forge, or where somebody with the rights to the site
environment can, this is the whole fix and it is one command:

```bash
conda install -c conda-forge fastmdxplora
```

FastMDXplora declares what it needs, so this brings the stack with it rather
than leaving you to assemble it. Into an existing environment that already has
OpenMM, the solver reconciles the two; where it cannot, it says so rather than
breaking what is there.

Where the machine cannot reach conda-forge but you can reach it somewhere
else, the packages can be carried:

```bash
# where there is a network, for the target's platform
conda create -p /tmp/for-transfer --download-only \
    -c conda-forge --platform linux-64 fastmdxplora
tar czf conda-packages.tgz -C /tmp/for-transfer .

# on the target
conda install --offline -c file:///path/to/unpacked fastmdxplora
```

The solve happens on the machine with the network, so it has to be told the
target's platform and Python version. That is the brittle part, and it is why
the second way exists.

**Take an image.** One attached to each release, built from
`container/fastmdx.def` and carrying the whole stack -- OpenMM with CUDA and
PLUMED, MDTraj, PDBFixer, the OpenFF toolkit, and FastMDXplora itself, all
resolved by one solver so the versions are the ones the packaging chose.
Download it where there is a network, copy the one file across, and:

```bash
apptainer exec --nv fastmdx.sif fastmdx explore --config study.yml
```

`--nv` passes the GPU through. Without it the run falls back to the CPU
platform -- which it will say, and which is worth reading rather than
discovering from the timing.

The image carries its own OpenMM, so `OPENMM_PLUGIN_DIR` is set inside it and
whatever the host has does not apply.

Building it yourself needs a machine that is x86-64 Linux **and** has a
network, which a laptop usually is not and a cluster usually does not. The
image installs FastMDXplora from conda-forge, so the release has to have
reached the feedstock -- an hour or so after a tag, sometimes longer:

```bash
apptainer build --fakeroot fastmdx.sif container/fastmdx.def
```

The definition ends by asserting that the CUDA platform and the ligand path
are both present, so a build that would have produced an image quietly
lacking them fails instead.

---

## Structures, offline

A run given a four-character identifier fetches it from RCSB, which needs a
network. Fetch them where there is one and pass the file instead:

```bash
curl -O https://files.rcsb.org/download/1UBQ.pdb
```

```yaml
systems:
  - system: 1UBQ.pdb        # rather than 1UBQ
```

The path is read from where `fastmdx` was run, not from the config that names
it.

---

## Making sure the GPU is being used

The banner names the platform, and it is worth reading:

```
Platform         CUDA (mixed precision)
```

`CPU` there on a machine with a GPU means something is wrong, and so does
`Reference` anywhere — that is OpenMM's correctness implementation, perhaps a
hundred times slower than the others. A run on it finishes and is correct, so
nothing else will tell you.

Ask for the platform rather than accepting `auto` for a first run:

```bash
fastmdx explore --config study.yml --output runs/study --simulate-platform CUDA
```

An explicit request fails loudly where the platform is unavailable. `auto`
falls back, which is right for ordinary use and wrong when you are checking
whether the hardware is reachable at all.

The run lists the platforms OpenMM found, and where one is missing says what
is known about why. Two answers come up often:

- **A library that would not load**, quoted from OpenMM — commonly
  `libcuda.so.1: cannot open shared object file`, which means the NVIDIA
  driver is not present. Correct on a head node with no GPU, and the reason
  to be on a compute node.
- **Nothing found and nothing failed**, which means nothing was tried. OpenMM
  loads its platforms from a directory fixed when it was built, and a
  relocated or shared installation moves it. Set `OPENMM_PLUGIN_DIR` to the
  `lib/plugins` beside the OpenMM library:

  ```bash
  export OPENMM_PLUGIN_DIR=/path/to/env/lib/plugins
  ```

  Worth putting in the job script rather than typing: without it every run is
  silently on `Reference`.

---

## Batch schedulers

A study is one command, so a job script is short:

```bash
#!/bin/bash
#SBATCH --job-name=fastmdx
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

export OPENMM_PLUGIN_DIR=/path/to/env/lib/plugins
fastmdx explore --config study.yml --output $SLURM_SUBMIT_DIR/runs/study \
    --simulate-platform CUDA
```

Write to a filesystem the compute nodes can see. Scratch is usually faster
than a home directory and usually not backed up, which matters for a
trajectory you cannot cheaply reproduce.

**A study of many runs schedules itself.** Umbrella sampling expands into one
run per window, and a sweep across systems into one run each; both go through
the same batch machinery:

```yaml
execution:
  mode: parallel
  workers: 4
```

Each worker takes a share of the cores rather than all of them, so the runs
divide the machine instead of fighting over it. With GPUs, list them and each
run is pinned to one:

```yaml
execution:
  mode: parallel
  devices: [0, 1, 2, 3]
```

---

## Watching from somewhere else

The GUI serves over HTTP, so an SSH tunnel is enough:

```bash
ssh -L 8000:localhost:8000 cluster
# on the cluster
fastmdx gui --output runs/study --port 8000
```

Then open `http://localhost:8000` at home. Check your site's policy first;
some clusters forbid listening sockets on compute nodes and expect this from
a login node instead.

For a long run, syncing the directory is simpler and needs no tunnel:

```bash
rsync -az --delete cluster:~/work/runs/study/ ~/runs/study/
fastmdx gui --output ~/runs/study
```

---

## What a run leaves behind

Everything needed to repeat it, in the output directory: `resolved_config.yml`
with every setting filled in, the prepared and solvated structures, the
serialised system and state, the trajectory, and the analyses and report.

`project_bundle.zip` gathers the parts worth sending to somebody else. It is
the thing to bring home if the trajectory is too large to move.
