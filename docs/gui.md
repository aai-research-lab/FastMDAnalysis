# The FastMDXplora GUI

For the shortest beginner path, start with [Beginner's guide](getting_started.md).
For GPU validation, staging, and checkpoint recovery, see
[Production and GPU runs](production.md).

## Opening the GUI

```bash
fastmdx gui
```

That serves the graphical interface on <http://127.0.0.1:8765> and opens a
browser. The GUI is the umbrella for everything described on this page: a
study builder that writes a config file or launches a run, live telemetry, the
molecular viewer, analysis results, and a file browser. Pass `--output DIR` to
open an existing run, `--port` to change the port, and `--no-browser` to skip
opening a window.

Use `--output DIR` to reopen a specific run rather than the current
directory.

FastMDXplora has two dashboard views:

- **Results Dashboard**: the static `report/dashboard.html` written after
  analysis/report. It shows completed plots, generated files, report links,
  slides, and project bundles.
- **Live Simulation**: a local-only localhost view that watches lightweight
  telemetry files while a simulation is running.

The live dashboard does not replace the static report dashboard. It adds a
separate monitoring view for progress, health messages, recent events,
**true live molecular coordinates**, **completed-trajectory playback**,
**structure and ligand tooling**, and available energy/temperature samples.

## Architecture

The live dashboard is a plain-HTTP Python server with no Node or frontend
build step. Assets are vendored locally and the dashboard runs fully offline.

```
src/fastmdxplora/gui/
├── server.py                # HTTP handlers + route table
├── telemetry.py             # existing telemetry reader
├── protein_preview.py       # existing static preview generator
├── structure_info.py        # PDB chains / residues / atoms / waters / ions
├── ligand_detection.py      # ligand + cofactor detection with overrides
├── live_frames.py           # atomic live-frame file writer
├── trajectory_playback.py   # DCD -> downsampled multi-MODEL PDB for browser
├── templates/
│   └── dashboard.html       # single-page dashboard (sidebar + pages)
└── static/
    ├── dashboard.css        # black scientific visual system
    ├── dashboard.js         # navigation, polling, state, events
    ├── molecule-viewer.js   # 3Dmol viewer (live + playback)
    ├── charts.js            # live telemetry charts
    ├── aai-research-logo.svg
    ├── aai-research-mark.svg
    └── 3Dmol-min.js         # existing vendored 3Dmol bundle
```

### Page layout

The single-page app exposes:

- **Sidebar** with the AAI Research Lab logo, the FastMDXplora product name, the
  full tagline (Fully Automated SysTem for Molecular Dynamics eXploration),
  the persistent nav (New Exploration, Overview, Live Simulation, Molecular
  Viewer, Analysis, Files & Reports, Run Settings, Documentation, GitHub), and
  dynamic footer info (connection status, current run, detected platform).
- **Top bar** with system/PDB, run title, status pill, current stage, step/total,
  platform, temperature, refresh, and Pause Updates (browser-only).
- **Pages**: New Exploration, Overview, Live Simulation, Molecular Viewer,
  Analysis, Files & Reports, Run Settings.

There is one page for building a run. There used to be two -- one starting from
a protein, one from a trajectory -- and they were the same thing: both wrote a
config and ran it, differing only in which phases the config named. Built
separately, one of them offered eleven of the eighty-three settings that exist
while the other offered all of them.

### Branding

The dashboard is branded with the AAI Research Lab logo (vendored SVG, no CDN).
The tagline `Fully Automated SysTem for Molecular Dynamics eXploration` is
preserved with its original capitalization. Black / charcoal backgrounds, white
and silver typography, restrained cyan and violet accents, green for healthy,
orange for warning, red only for error. All system fonts; no remote font loads.

## Safety and network boundary

The dashboard must never terminate or interfere with an OpenMM simulation:

- All file-system writes performed by the dashboard module are wrapped in
  try/except. A failure is logged and surfaced as `unavailable` in the UI
  without stopping the surrounding workflow.
- Live frame writing runs at a separately tunable cadence
  (`--dashboard-frame-interval`, default approximate match to telemetry).
- The simulation stays authoritative. The **Pause Updates** control only
  pauses browser-side polling; it never pauses the OpenMM integrator.
- Telemetry parsing tolerates partial lines, missing files, and mid-write
  arrays.

The default server binds to `127.0.0.1`, which is the recommended and safest
mode. Do not expose the dashboard directly to the internet. When you bind to
`0.0.0.0`, `::`, or another non-loopback address, the dashboard remains useful
for monitoring but workflow launch/stop and local-folder opening are disabled.
Absolute host filesystem paths are not sent in the artifact listing.

## Recommended: start it with the workflow

For normal runs, add `--dashboard` to the command that creates or watches the
output folder. FastMDXplora starts the localhost dashboard before setup or
simulation begins, prints the URL, and automatically enables live telemetry
when simulation runs:

```bash
fastmdx explore --system local_pdbs/1L2Y.pdb \
  --output local_runs/trpcage_live_full \
  --include setup simulation analysis report \
  --simulate-preset gentle \
  --dashboard
```

If the `fastmdx` console script is not on PATH, use the module entrypoint.
This is also the recommended form for Windows PowerShell:

```powershell
python -m fastmdxplora.cli.main explore `
  --system local_pdbs\1L2Y.pdb `
  --output local_runs\trpcage_live_full `
  --include setup simulation analysis report `
  --simulate-preset gentle `
  --dashboard
```

The CLI prints the selected URL before the workflow starts:

```text
Live dashboard running at: http://127.0.0.1:8765
Watching output folder: local_runs/trpcage_live_full
Open this URL in your browser to monitor the run.
Press Ctrl+C to stop the dashboard after the workflow completes.
```

By default the server binds to `127.0.0.1:8765`, so it is local to your
machine. If the requested port is busy, FastMDXplora chooses the next
available port and prints the actual URL:

```text
Live dashboard running at: http://127.0.0.1:8766
Requested port 8765 was busy, so FastMDXplora used 8766.
Watching output folder: local_runs/my_run
```

Use `--dashboard-port` or `--dashboard-host` to customize the bind address:

```bash
fastmdx explore --system protein.pdb --output local_runs/my_run --dashboard --dashboard-port 8770
```

Binding to all interfaces can expose the dashboard on your network:

```bash
fastmdx explore --system protein.pdb --output local_runs/my_run --dashboard --dashboard-host 0.0.0.0
```

FastMDXplora prints a warning when `--dashboard-host 0.0.0.0` is used. Prefer
`--dashboard-host 127.0.0.1` for local-only access.

By default the dashboard stays open after the workflow completes so you can
inspect final Analysis Plots and Generated Files. Press Ctrl+C to stop it. Use
`--dashboard-stop-on-complete` when you want the command to exit immediately
after the workflow finishes.

## Manual fallback for existing runs

You can still reopen an existing output directory without running a workflow:

```bash
fastmdx gui --output local_runs/my_run
```

If the `fastmdx` console script is not on PATH:

```bash
python -m fastmdxplora.cli.main gui --output local_runs/my_run
```

## Saving a config file instead of launching

The **Builder** section collects the same options the CLI accepts. Two buttons
act on that selection:

- **Start Exploration** starts the run immediately on this machine.
- **Save config file** downloads a `.yml` study configuration instead.

The saved file is the canonical FastMDXplora config format, so it runs
anywhere:

```bash
fastmdx explore --config my_run.yml
```

This is the route to use when the machine that will run the simulation is not
the machine with the browser. Design the study in the GUI on your laptop, save
the config, copy it to a workstation or cluster, and submit it there. The file
is also worth keeping alongside the results: it records exactly what was run,
and command-line flags override it, so one file can drive a family of related
runs.

## Watching a run on a cluster

The GUI binds to `127.0.0.1` and serves a browser interface, so it cannot be
used directly on a compute node reached through a batch scheduler. There is no
browser there, the node is usually not addressable, and the exploration layer is
meaningless when the scheduler owns job submission.

Telemetry is written to plain files inside the output directory, so the
practical pattern is to sync the output directory to a machine that does have a
browser, and point the GUI at the copy:

```bash
# on your workstation, while the cluster job runs
rsync -az --delete user@cluster:/scratch/$USER/runs/my_run/ ~/runs/my_run/
fastmdx gui --output ~/runs/my_run
```

To keep it current, re-run the sync on a loop in a second terminal:

```bash
while true; do
  rsync -az --delete user@cluster:/scratch/$USER/runs/my_run/ ~/runs/my_run/
  sleep 60
done
```

The dashboard then updates on its normal polling interval, delayed by at most
one sync period. Telemetry, energy samples, live coordinates, and finished
analyses all transfer, because they are ordinary files.

Two practical notes. Use `--delete` so removed files do not linger in the
copy, and exclude the trajectory if it is large and you only want progress:
`--exclude 'production.dcd'` keeps the sync small, at the cost of the playback
view.

If you can open an SSH tunnel to the compute node and the job is interactive,
`ssh -L 8765:<node>:8765 user@cluster` with `fastmdx gui --output <run>` on the
node also works. The rsync route is usually less trouble.

## Live telemetry

`--dashboard` automatically implies `--simulate-live-telemetry` when simulation
runs. Existing explicit telemetry flags still work for tuning:

```bash
fastmdx explore --system protein.pdb \
  --output local_runs/my_run \
  --include setup simulation \
  --dashboard \
  --simulate-telemetry-interval 1000
```

The simulation phase writes:

- `simulation/live_status.json`
- `simulation/live_metrics.csv`
- `simulation/live_events.log`

Telemetry writing is intentionally lightweight and safe. If these files cannot
be written, the simulation continues and the live dashboard simply reports
unavailable data.

## What the Live Simulation tab shows

The live view polls the local output directory every few seconds and displays:

- current stage and status
- current step and planned step count when known
- frame count when known
- elapsed wall time and simulation time
- OpenMM platform
- latest checkpoint path when available
- energy/temperature trends when telemetry has those values
- recent events from `live_events.log`
- protein preview image when a topology/PDB is available
- optional local interactive 3D preview when a topology/PDB is available
- plain-language health explanations

Metrics are not invented. If a value has not been written yet, the dashboard
shows `not available`.

## Protein preview

The Live Simulation tab tries to show a protein image as soon as a usable
structure file exists. It checks common run artifacts such as
`simulation/topology.pdb`, `setup/topology.pdb`, `setup/solvated.pdb`, and
the original system path recorded in `manifest.json`.

If PyMOL is installed, FastMDXplora uses it for a cartoon/ribbon render and
adds extra camera padding so the whole protein is visible in the dashboard
frame. Install PyMOL from conda-forge:

```bash
conda install -c conda-forge pymol-open-source
```

or with micromamba:

```bash
micromamba install -c conda-forge pymol-open-source
```

Protein preview generation requires PyMOL so the static image is an actual
cartoon/ribbon render rather than a schematic placeholder. Preview files are
written to `report/dashboard_assets/protein_preview.png` when report assets
exist, or to `simulation/protein_preview.png` during setup/simulation-only
runs.

The preview panel uses a local vendored 3Dmol.js bundle for **Interactive 3D**
when a structure file exists. That tab shows a rotatable cartoon representation
with spectrum coloring and does not require a CDN or internet access. Use
**Spin** to call the real 3Dmol viewer spin control, **Reset view** to zoom back
to the full protein, and mouse controls to rotate or zoom.

The **PyMOL Preview** tab remains available whenever the PyMOL PNG exists. That
tab shows the PyMOL-rendered cartoon/ribbon image. If the bundled 3Dmol viewer
asset is unavailable, the dashboard can fall back to a schematic CA/backbone
trace, but that fallback is labeled as a schematic fallback, not as PyMOL. If no
usable structure file exists yet, the panel says that the preview is
unavailable. The panel polls for updates and includes a **Regenerate preview**
button.

The Results Dashboard tab also refreshes while the server is open. If analysis
or report artifacts are generated after the server starts, the Generated Files,
Analysis Plots, Report Artifacts, and embedded dashboard view update without
restarting the server. Use the **Refresh** button for an immediate rescan.

## Health messages

FastMDXplora classifies common live issues:

- **Numerical instability**: NaN/Inf positions or energies. Try a smaller
  timestep, lower temperature, stronger friction, the gentle preset, or check
  the input structure for clashes.
- **Energy spike**: energy changed sharply between samples. This can indicate
  bad contacts, an unstable timestep, or pressure/temperature coupling issues.
- **Temperature drift/spike**: temperature is far from the target. Try gentler
  heating, stronger friction, or a smaller timestep.
- **Stale telemetry**: no live update was seen recently. The simulation may be
  slow, paused, or crashed.

## Static dashboard behavior

Opening `report/dashboard.html` after a run still shows the completed Results
Dashboard. It now also includes a Live Simulation sidebar entry:

- If telemetry exists, it shows the last recorded status.
- If telemetry does not exist, it explains how to start
  `fastmdx explore ... --dashboard` or `fastmdx gui --output ...`.

## Demo or preview output

For a small demo run, use the gentle preset and telemetry:

```bash
fastmdx explore --system local_pdbs/1L2Y.pdb \
  --output local_runs/live_demo \
  --include setup simulation analysis report \
  --simulate-preset gentle \
  --dashboard
```

## New CLI options

The redesigned dashboard accepts additional options without changing existing
defaults:

| Flag | Purpose | Default |
| --- | --- | --- |
| `--dashboard-host` | Bind address for dashboard server | `127.0.0.1` |
| `--dashboard-port` | Bind port; auto-falls-forward if busy | `8765` |
| `--dashboard-refresh-seconds` | Browser polling interval | `3` |
| `--dashboard-frame-interval` | Telemetry interval at which a live frame PDB is written | matches telemetry |
| `--dashboard-max-playback-frames` | Max browser playback frames (downsampled) | `200` |
| `--dashboard-ligand-resname` | Explicit ligand residue-name override | autodetect |
| `--dashboard-binding-pocket-cutoff-A` | Pocket distance cutoff in Å | `5.0` |
| `--dashboard-open-browser` | Auto-open the dashboard URL | `False` |
| `--dashboard-stop-on-complete` | Exit server when workflow completes | `False` |

## API surface

The server exposes lightweight JSON endpoints (no database):

- `GET /api/status` — current run status snapshot
- `GET /api/metrics` — recent telemetry samples
- `GET /api/events` — recent telemetry events
- `GET /api/structure-info` — chains / residues / atoms / waters / ions / bbox
- `GET /api/ligands` — ligand instances, residue names, cofactors, explicit override
- `GET /api/live-frame-index` — current live-frame index + timestamp
- `GET /api/live-coordinates` — last-known live-frame update timestamp
- `GET /api/playback-info` — downsampled playback frames metadata
- `GET /api/files` — alias for `/api/artifacts` generated file list
- `GET /api/analyses` — alias for `/api/results` analysis summary
- `GET /api/schema` — every setting the software accepts, with its help text,
  accepted values, default, and the control to draw for it; plus each
  analysis's own settings, read from the analyses themselves
- `GET /api/browse` — folders, and files of a given `kind` (trajectory,
  structure, config), for the file picker
- `GET /api/inspect-directory` — what a folder holds and what that makes
  possible
- `POST /api/config` — the config the page describes, validated
- `POST /api/check-config` — whether a config someone already has will run
- `POST /api/load-config` — that config, in the shape the form works in
- `POST /api/run` — start what the page describes
- `POST /api/run-config` — start a config exactly as it stands
- `GET /structure/topology.pdb` — served structure (latest available)
- `GET /structure/live-frame.pdb` — latest OpenMM live frame (atomic swap)
- `GET /structure/playback.pdb` — downsampled multi-MODEL PDB for playback

All endpoints degrade gracefully: missing files, malformed JSON, or a busy
parser return a benign `{"ok": false, "reason": "..."}` document rather than
crashing.

The home-mode simulation builder also uses these POST endpoints, but only when
the server is loopback-bound:

- `POST /api/explore/validate` — validate a proposed workflow
- `POST /api/explore/start` — launch the canonical CLI workflow
- `POST /api/explore/stop` — terminate the active workflow and wait for exit

## Live molecular coordinates and trajectory playback

When the simulation phase writes frame telemetry, the runner also writes
`simulation/live_frame.pdb` atomically (tmp file + `os.replace`). The dashboard
polls the live-frame index and merges the new coordinates into the 3Dmol viewer
**without re-centering or rebuilding the UI**, so the user's camera orientation,
zoom, and visibility settings are preserved frame-to-frame.

When the run finishes, the dashboard offers a downsampled **trajectory
playback** view:

- Up to `--dashboard-max-playback-frames` frames are sampled evenly from
  the full DCD trajectory. The scientific DCD itself is never modified.
  Compilation streams the DCD in bounded chunks and caps browser requests at
  2,000 frames.
- Concurrent playback requests are serialized per run so the companion PDB
  and index cannot be replaced by competing writers.
- Controls: Play, Pause, Reverse, Previous, Next, Jump to start/end,
  frame slider, playback speed, loop, screenshot, frame download.
- Each frame is `MODEL`/`ENDMDL`-wrapped inside a single multi-MODEL PDB so
  the in-browser 3Dmol viewer can use `.mload()` for instant frame switching.

## Ligand and binding-pocket tooling

After ligand detection (or with an explicit `--dashboard-ligand-resname`):

- Center on ligand, isolate ligand, show labels, hide distant residues.
- Show binding-pocket residues within `--dashboard-binding-pocket-cutoff-A`
  (default 5.0 Å), computed per-atom (not centroid) so non-spherical ligands
  are captured correctly.
- Geometric contacts and candidate hydrogen bonds are computed by 3Dmol
  distance criteria and labeled as geometric, not as confirmed chemical
  interactions.

## Empty, waiting, and error states

Styled states are provided for:

- Waiting for structure / simulation (subtle molecular-network background)
- Protein-only run (no ligand detected)
- Telemetry stale (no recent update within refresh window)
- Simulation completed (final structure + completed trajectory)
- Simulation error (the actual known error without inventing causes)
- Viewer unavailable (falls back to the existing PyMOL preview / schematic)

All empty/waiting states carry a faint AAI Research Lab watermark.

## Manual browser testing checklist

The redesign was manually validated against:

- Protein-only simulation (e.g. trp cage) and protein-ligand simulation
  (e.g. 1L2Y/EPE)
- CPU and CUDA/OpenCL platforms
- Runs with no ligand, no telemetry yet, missing topology, completed run,
  failed run, interrupted run, and large or small trajectories
- Dashboard reload mid-run, changed dashboard port, and stale telemetry
- Narrow browser widths (sidebar collapses)
- A package install at non-checkout locations to confirm templates and
  static assets ship correctly

## Building a run

Running FastMDXplora without a subcommand starts the local dashboard before a
project exists:

```powershell
fastmdxplora
```

FastMDXplora binds the first available local port beginning at `8765`, prints
the actual URL in the terminal banner, and opens the browser when possible.
Set `FASTMDX_NO_BROWSER=1` to suppress automatic browser opening in headless
or automated environments.

The **New Exploration** page asks three questions, in the order they are
answered.

**What have you got?** A structure -- a four-character PDB identifier or a path
to a file -- or a trajectory you already have, from anywhere, or a config
written earlier. That answer sets what happens next.

**What should happen?** Setup, Simulation, Analysis and Report, pre-selected
from the first answer and adjustable. A run that only analyses is a run with
one of them ticked. Where a phase has nothing to act on it is shown but
unavailable, with the reason: a trajectory is already the result of setting up
and simulating, and there is no supported way to continue a run from one.

**Anything to change?** Every setting for the chosen phases, read from the
schema, with the help text the command line and the config template use. When
Analysis is included, each measurement appears with what it means and how to
read it. Everything has a default that works; the settings are there for when
it does not. **Reset to defaults** puts it all back.

Then the run can start here, or the config can be downloaded and run
elsewhere -- the same file either way. A laptop is a poor place to run a long
simulation; a cluster is a poor place to decide what to run.

Nothing on the page holds a list of what can be set. The controls are drawn
from what the schema reports, and the per-analysis settings from the analyses
themselves, so a new field appears without anything else being updated.

The dashboard does not implement a separate simulation engine. It writes a
config and launches the canonical `fastmdx explore --config` command; the child
process writes normal project artifacts and live telemetry, while the
already-running dashboard switches to that output directory. The config is
written beside the results, so a run can be repeated from its output directory
rather than from what somebody remembers clicking.

### A config you already have

It can be checked -- reporting YAML syntax errors, settings that do not exist,
and on success what the file will do -- then either run exactly as it stands or
opened into the form and changed. Opening never writes to it: the file may be
committed beside a paper or be the record of a run that already happened, so
anything changed is saved as a new file.

### Choosing files

Every field wanting a path has a **Browse…** button. A browser will not offer a
dialog for a path the server has to open -- the one it has uploads files to a
page -- so the walking is done by the server and drawn by the page. Only files
of the kind being asked for are listed, and folders are marked with how much
data they hold, counted by extension a few levels down.

Click **Validate** before launching. Validation checks both the form and the
Python environment behind the dashboard. If OpenMM or PDBFixer is missing,
the page shows the exact repair steps:

```bash
conda activate fastmdxplora
conda install -c conda-forge openmm pdbfixer
```

Run the command in a terminal, restart the dashboard with the same environment,
click **Validate** again, and then click **Start Exploration**. The dashboard
does not start a workflow that it already knows cannot perform chemistry setup.

Only one dashboard-launched workflow is active at a time. Scientific controls
become read-only once a workflow starts; changing browser fields does not
modify a running OpenMM context. A completed configuration can be adjusted and
launched again under a new run name.
