# Production and GPU runs

This page covers the checks that should happen before a long simulation. The
commands are intentionally explicit, but paths and environment names are
placeholders until the target machine is verified.

## Where production runs belong

Use a suitable compute host for a long production MD simulation. A local
workstation is appropriate for installation, analysis/reporting, and bounded
smoke tests; a dedicated GPU host is usually preferable for long trajectories.

Common deployment models include:

1. **Local GPU:** verify the installed driver and OpenMM CUDA platform, then
   run locally with an explicit device index.
2. **Connected HPC or cloud GPU:** inspect the host and launch through the
   provider's documented SSH, scheduler, container, or session workflow.
3. **Offline or air-gapped GPU:** manually stage a sanitized package and all
   dependencies/input files, then execute the commands on the host. Do not
   assume SSH, a scheduler, remote paths, or internet access.

## Checks on any GPU host

Run these read-only checks in the same environment that will run FastMDXplora:

```bash
nvidia-smi
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
micromamba env list
micromamba run -n <verified-env> python -V
micromamba run -n <verified-env> python -m fastmdxplora.cli.main --version
micromamba run -n <verified-env> python -m fastmdxplora.cli.main info
```

Check OpenMM's actual platform list and optional packages:

```bash
micromamba run -n <verified-env> python - <<'PY'
import openmm as mm
print("OpenMM platforms:", [
    mm.Platform.getPlatform(i).getName()
    for i in range(mm.Platform.getNumPlatforms())
])
for module_name in ("openff.toolkit", "openmmforcefields", "openmmplumed"):
    try:
        __import__(module_name)
        print(module_name, "OK")
    except Exception as exc:
        print(module_name, "FAIL", type(exc).__name__, str(exc).splitlines()[0])
PY
```

A listed CUDA platform is not proof that a production context will initialize.
The actual host driver must support the CUDA/OpenMM build in the environment.
A historical failure mode was `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`; do not
retry a long run unchanged after that error. Re-check the driver, CUDA/OpenMM
build, and selected environment.

## Connected HPC: safe launch pattern

Use the existing alias or site documentation rather than placing credentials,
private keys, or sensitive host details in a script:

```bash
ssh <site-approved-alias>
cd /path/to/FastMDXplora
micromamba run -n <verified-env> fastmdx explore \
  --config /path/to/resolved_config.yml \
  --output /path/to/new_production_run
```

If a scheduler is required, use the site's documented submission command and
record the job ID. Do not guess a partition, account, module name, or scratch
path. For an interactive/background launch, use the site's approved `tmux` or
session workflow.

Before launch, use a separate fresh output directory and set checkpointing:

```bash
micromamba run -n <verified-env> fastmdx explore \
  --system /path/to/prepared_or_input.pdb \
  --output /path/to/production_run_001 \
  --simulate-platform CUDA \
  --simulate-device-index 0 \
  --simulate-duration-ns 100 \
  --simulate-checkpoint-interval-steps 10000 \
  --dashboard
```

`--dashboard` is useful for a manually monitored run, but it does not replace
scheduler/process monitoring. Keep the dashboard on loopback unless the host's
secure forwarding/access setup is understood.

## Offline or air-gapped GPU hosts

An offline host must be operated manually by its administrator or job
operator. The operator must:

1. create a sanitized staging directory;
2. upload the source/package and all input files manually;
3. run the environment/GPU checks below;
4. run a bounded smoke test;
5. launch production only after reviewing the smoke-test output; and
6. manually transfer results back when needed.

Do not upload `.env` files, credentials, private keys, tokens, or unrelated
home-directory contents. Because the machine is offline, stage all packages,
PDB/CIF files, ligand SDF/MOL2 files, topology/state files, config files, and
PLUMED assets before starting.

### Manual staging template

On a connected staging machine, create an archive from a sanitized directory
only:

```bash
mkdir -p offline_stage/fastmdxplora_source
# Copy only the reviewed source/package and required inputs into offline_stage.
# Do not copy .env files, credentials, private keys, caches, or large unrelated data.
tar -czf fastmdxplora_offline_bundle.tgz -C offline_stage .
shasum -a 256 fastmdxplora_offline_bundle.tgz
```

Transfer `fastmdxplora_offline_bundle.tgz` using the site's approved offline
transfer method. On the offline GPU host:

```bash
mkdir -p ~/fastmdxplora_jobs/job_001
cd ~/fastmdxplora_jobs/job_001
tar -xzf /path/to/fastmdxplora_offline_bundle.tgz
nvidia-smi
micromamba env list
micromamba run -n <verified-env> python -m pip install -e ./fastmdxplora_source --no-deps
micromamba run -n <verified-env> python -m fastmdxplora.cli.main info
```

Run one smoke test on one GPU device:

```bash
micromamba run -n <verified-env> fastmdx explore \
  --system inputs/protein.pdb \
  --output smoke_gpu0 \
  --include setup simulation \
  --simulate-preset gentle \
  --simulate-platform CUDA \
  --simulate-device-index 0 \
  --simulate-checkpoint-interval-steps 100
```

For an independent second smoke test, use a new output directory and device 1
when the host has another available GPU:

```bash
micromamba run -n <verified-env> fastmdx explore \
  --system inputs/protein.pdb \
  --output smoke_gpu1 \
  --include setup simulation \
  --simulate-preset gentle \
  --simulate-platform CUDA \
  --simulate-device-index 1 \
  --simulate-checkpoint-interval-steps 100
```

Do not run two jobs on the same device. Do not start production until the
selected environment and the real CUDA/OpenMM smoke test are understood.

## Real PLUMED smoke test

PLUMED code wiring and mocked tests are separate from a real backend test. On
the target GPU machine, use a minimal PLUMED script with atom indices valid for
the staged topology:

```bash
cat > inputs/plumed_smoke.dat <<'PLUMED'
d: DISTANCE ATOMS=1,2
PRINT ARG=d STRIDE=100 FILE=COLVAR
PLUMED
```

If the offline host does not have the required PLUMED environment, create the
file before transfer rather than trying to install packages online. Run a short
separate output:

```bash
micromamba run -n <verified-env> fastmdx explore \
  --system inputs/protein.pdb \
  --output smoke_plumed_gpu0 \
  --include setup simulation \
  --simulate-preset gentle \
  --simulate-platform CUDA \
  --simulate-device-index 0 \
  --simulate-plumed-script inputs/plumed_smoke.dat \
  --simulate-checkpoint-interval-steps 100
```

Confirm all of the following before production:

- OpenMM created a CUDA context and completed the bounded run.
- Energies and temperatures are finite.
- `simulation/production.dcd`, `energy.csv`, `simulation.log`, and
  `checkpoint.chk` exist.
- `simulation/plumed.dat` and `COLVAR` exist when PLUMED is enabled.
- The run manifest and telemetry report the expected status.

`HILLS` is written only when the PLUMED script contains a `HILLS` action.

## Monitoring

For an active CLI/dashboard run:

```bash
tail -f /path/to/run/simulation/simulation.log
watch -n 30 nvidia-smi
jq . /path/to/run/simulation/live_status.json
tail -n 20 /path/to/run/simulation/live_events.log
```

For a scheduler job, use the site's documented queue/accounting commands. For
a manually launched session, monitor the process, GPU, log, telemetry
timestamp, checkpoint modification time, and output-file growth together.

## Completion and recovery

Treat a run as complete only when the process or scheduler reports success,
the manifest/telemetry reports completion, and the expected files are readable.
A DCD or checkpoint file alone is not proof of completion.

For a failure or interruption:

1. preserve the original output directory;
2. record the last known step and checkpoint timestamp;
3. confirm the checkpoint belongs to the same system, topology, integrator,
   precision, device-compatible environment, and PLUMED history policy;
4. resume into a new directory or explicitly documented continuation directory;
5. avoid overwriting or blindly concatenating `production.dcd`, `COLVAR`,
   `HILLS`, `energy.csv`, and logs; and
6. validate each bounded continuation segment before proceeding.

The current runner writes `checkpoint.chk`, but the inspected CLI does not
provide a complete high-level resume command. Checkpoint recovery therefore
requires a deliberate continuation procedure rather than a blind restart.

## Related pages

- [Beginner's guide](getting_started.md)
- [Live dashboard](gui.md)
- [CLI reference](cli_reference.md)
- [Configuration](configuration.md)
