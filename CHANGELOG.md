# Changelog

All notable changes to FastMDXplora are documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) ·
Versioning: [SemVer 2.0.0](https://semver.org/spec/v2.0.0.html)

## [Unreleased]

Three changes here alter numbers a 2.2.0 run produced. Hydrogen-bond counts
rise, because a bond present in few frames was being left out of all of them.
The `v1` analysis profile is gone, so a command that named it now fails rather
than quietly applying a table of settings. Clustering keeps its default and
gains an alternative, which changes nothing unless asked for.

### Removed
- **`--compat v1` and `--analyze-compat v1`.** The profile applied a table of
  analysis options described as version 1's. Three of them were this
  software's own defaults, one named a parameter version 1 does not have, and
  two were neither implementation's: the hydrogen-bond entry multiplied the
  per-frame count by two, and the clustering entry asked for six clusters
  where version 1 falls back to three.

  Neither version 1 nor MDTraj counts a hydrogen bond twice — version 1 takes
  `len(baker_hubbard(frame))`, one row per donor-hydrogen-acceptor triplet —
  so the multiplier modelled nothing. And the clustering entry could not have
  reproduced version 1 at any number, because the two compare frames in
  different spaces.

  Reproducing a published result means running the same method and finding the
  same number. A flag that supplies the number removes the thing being
  checked. State the settings instead: `--scope`, `--selection`, `--stride`,
  and the per-analysis options, all of which `fastmdx analyze --help` lists.

### Fixed
- **The Q-value is measured as the paper it cites defines it.** Best, Hummer &
  Eaton give Q through a switching function: each native contact is judged
  against the distance it had natively, and stops counting gradually rather
  than at a step. A single threshold was applied to every contact instead, so
  one formed at 0.20 nm was held to the same standard as one formed at 0.44 nm.
  On a hairpin coming apart, the threshold reached zero where the published
  measure still read 0.62, and crossed a half at frame 27 where the paper's
  crossed at 89. **Q values will differ from 2.2.0, and for a folding study
  that is the result rather than a detail.** `beta` and `lambda_factor` are
  exposed, at the paper's values.

  One consequence: Q at the reference frame is now slightly under 1 rather
  than exactly 1. That is inherent to a smooth measure, and the paper leaves
  it unnormalised, so rescaling it to 1 would misstate what was measured.
- **Secondary structure leaves out what has no backbone.** DSSP returns a
  column for every residue and marks the ones without a backbone `NA`. Those
  columns were kept, and since `NA` is not in the colour map they were drawn
  as coil; worse, the residue labels were taken from the protein residues
  alone, so the two lists came out different lengths and the numbering fell
  back to counting from zero. **A protein numbered 10 to 15 was relabelled 0
  to 6, with the ligand as the last row.** Scope defaults to protein and
  ligand, so this was the ordinary case for a protein-ligand study.

  Where nothing in the selection has a backbone -- a nucleic acid, a lone
  ligand, a coarse-grained model -- the analysis now says so instead of
  drawing an empty timeline.
- **An atom column is always a number.** `Atom.serial` is supplied by the file
  a topology came from; a trajectory built in memory has none, and `None`
  became NaN in the RMSF and ligand-RMSF atom column. The saved data carried
  the NaN, and the figure cast it to the most negative integer there is and
  used that as an axis label. Where the file gives no serial, the atom's
  position is used.
- **Hydrogen bonds are counted in every frame they occur in.** Baker-Hubbard
  proposes bonds above an occupancy threshold, and only proposed bonds were
  evaluated frame by frame — so a bond present in five per cent of frames
  contributed to none of them, including the frames it was in. The series is
  the number of hydrogen bonds per frame, and it was reporting the number of
  persistent ones. **Counts will be higher than 2.2.0 reported.** Raise
  `candidate_freq` to restrict the series again.
- **`freq` reports what it names.** Its only effect had been deciding which
  bonds were proposed; it now records how many are persistent at that
  threshold, alongside how many were found, in the analysis manifest.
- **An analysis option the software cannot apply stops the run.** Every
  analysis ends its signature with `**kwargs`, which the base class stores and
  nothing reads, so a misspelled setting was accepted, ignored, and the run
  reported success. Asking for `n_clusteres` clustered at the default and said
  nothing about it.
- **`--analyze-dimred-components` had never worked.** The flag became an option
  name by dropping the analysis prefix, which gives `components`, while the
  option is `n_components`. Nothing accepted it, so anyone who set it reduced
  to two components and was told the run succeeded.
- **A default is declared once.** Three phase tables restated ninety values the
  schema also declared, and four lists of accepted values existed in the
  command line, the browser, and beside the code validating them. `DEFAULT_PH`
  had already stopped agreeing: it sat unread in the setup package at 7.0
  after the pH default became 7.4. A default must now also be one of its own
  accepted values, which cannot be stated while the two live apart.

### Added
- **`--cluster-features`**, which states what clustering compares frames in.
  `rmsd` (the default, and what 2.2.0 did) superposes every pair optimally, so
  the distance cannot depend on where the molecule sits. `coordinates`
  superposes each frame onto the first and compares directly, scaled by
  1/sqrt(n_atoms) so a distance is still an RMSD in nm — the cheaper
  approximation, and the space `ward` and k-means were defined in.

  The choice changes the answer more than any parameter here does. Comparing
  coordinates without superposing, as version 1 does, makes the leading
  difference between frames where the molecule drifted and how it turned: on a
  trajectory with a hinge motion and rigid-body drift, it recovers the two
  conformations no better than chance, where pairwise RMSD recovers them
  completely.
- **`--setup-protonation-margin`**, so a setting the software tells you to
  narrow can be reached without editing a config file.
- Every analysis now describes the settings it accepts, read from its own
  constructor and docstring rather than restated anywhere. A setting nobody
  has explained fails a check, so the description cannot fall behind.

## [2.2.0] — 2026-08-03

Protein-ligand systems from a PDB identifier alone. A crystal structure
carries the ligand you care about alongside the buffer that kept the protein
soluble, and a PDB record says which is which only by name. FastMDXplora now
decides, retrieves the chemistry the structure omits, and refuses where the
structure does not determine the answer.

### Added
- **Automatic protein-ligand preparation.** `--setup-heterogens auto` inspects
  the structure, decides what each non-standard residue means, retrieves the
  chemistry for anything worth simulating, and prepares it. A bound ligand no
  longer has to be extracted and rebuilt by hand.
- **Ligand chemistry from the Protein Data Bank.** A PDB record carries no bond
  orders, formal charges, or hydrogens, and a force field needs all three.
  These are retrieved from RCSB at the crystallographic pose, completed with
  hydrogens, and cached, so a run that worked once works again offline.
- **Protonation decided in the complex.** A ligand's pKa in a binding site is a
  property of the complex, not of the molecule: a buried acid can sit several
  units from its solution value. PROPKA is asked about the bound state, in a
  structure repaired first so the electrostatics are those of the system that
  will be simulated. States decisively away from the pH are adopted and
  reported with their environment shift; poised ones are refused.
- **Several ligands, cofactors, or copies** may be parameterized together.
  Each is placed and clash-checked against everything already present, since
  two ligands can overlap each other as readily as they can the protein.
- **`--setup-heterogens`** with `auto` (the default), `drop`, and `keep`.

### Changed
- **The default force field is now chosen for you, and it has changed.**
  `forcefield` defaults to `auto`, which resolves to `amber-openff`: the
  ff14SB protein model with TIP3P water and the OpenFF Sage small-molecule
  force field. Previously the default was CHARMM36. **Runs that relied on the
  default will produce different numbers from this release onward.** Name
  `charmm36` explicitly to keep it.

  `auto` resolves to one stack whatever the structure contains, deliberately.
  Choosing per-system would give an apo run and its holo partner different
  protein force fields, and the comparison between them is usually the point.
  CHARMM36 remains protein-only here because its native small-molecule
  partner is CGenFF; pairing it with OpenFF would be an unvalidated mixture.
- **`heterogens` now defaults to `auto`.** Previously every non-standard
  residue was discarded. **A run that prepared an apo protein from a holo
  structure will now prepare the complex, or stop and say why.** Pass
  `--setup-heterogens drop` for the old behaviour.

  What decided it is not that `auto` is convenient. Under `drop`, rhodopsin
  prepares as opsin, haemoglobin as globin, and streptavidin without its
  biotin — each completing in silence and reading like an answer to the
  question that was asked. Across thirty ordinary structures `auto` now
  prepares the ligand in twelve, matches `drop` in eight, and stops in ten,
  every stop naming what the structure left undetermined. None fails without a
  reason.
- **The default pH is physiological.** `ph` moves from 7.0 to 7.4: blood is
  7.4, and a protein studied without a stated reason otherwise is studied
  there. Cytosol sits near 7.2 and a lysosome near 4.7, so a
  compartment-specific study should say which. **Titratable groups near the
  old value may settle differently.**
- **`protonation_margin` is now a setting.** How close a ligand's pKa may come
  to the pH before setup stops rather than choose. The default of one unit is
  about the uncertainty of the pKa calculation itself, so the band marks where
  the answer is unresolved rather than merely close; narrowing it is for a
  ligand whose protonation is already known.
- **Centre-of-mass motion is removed by default**, matching OpenMM's own
  default. The previous setting let the system drift as a whole.
- **Ambiguity is refused rather than resolved.** Where a structure does not
  determine what should be simulated, setup stops and says what must be
  decided: a covalently bonded adduct, an unidentified component, a cofactor
  needing parameters a small-molecule force field cannot supply, alternate
  conformations at indistinguishable occupancy, a metal coordinated in some
  copies and not others, or a sugar that may be a glycosylation site, a
  substrate, or a cryoprotectant. Producing a plausible trajectory from a
  guess is worse than producing none.
- Removing a heterogen is announced. A bound ligand and a buffer molecule are
  indistinguishable in a PDB file, and both used to be discarded silently, so
  a run could be apo while reading as holo.

### Fixed
- **A ligand is parameterized in the charge state the pocket implies.** The pKa
  settled in the complex was computed, logged, and then discarded: the file
  handed to the small-molecule force field carried whatever protonation the
  reference chemistry held. Retinoic acid was simulated as a neutral acid and
  4-hydroxytamoxifen as a neutral amine, in both cases losing the charge that
  binds them. Benzamidine revealed it only by crashing, on a stereocentre the
  amidinium it should have been does not have.
- **An amidine is no longer read as an amine.** Its cation is delocalised and
  forms on the sp2 nitrogen; built on the sp3 one it is a different molecule.
  Benzamidine is the ligand of half the trypsin structures in the PDB.
- **Each ionizable group is settled on its own pKa.** One answer for the whole
  ligand forced an amino acid onto one side rather than building the
  zwitterion it is at pH 7.4.
- **A pyridine-type nitrogen is read.** Quinazolines, pteridines and purines
  were silently untitratable.
- **A glycan is identified by what it is bonded to.** A LINK to an asparagine,
  serine or threonine says a sugar is a glycosylation site rather than leaving
  the question open; those are dropped and the protein prepared without them.
  A sugar bonded only to other sugars, as in lysozyme's substrate, still is a
  question.
- **A nucleotide chelating a magnesium is not a covalent adduct.** Both ends of
  a LINK record were marked bonded, so ATP and GTP complexes — among the most
  common structures there are — refused as though bound to the protein.
- **A group the pKa calculation reports but the molecule does not carry no
  longer stops a run.** Folate's N10 lies between a methylene and an aromatic
  ring; read as an aliphatic secondary amine, it was given that class's model
  pKa of 10.0, when an aniline is nearer 4.6.
- **A coordinated ion is prepared, not refused.** Connectivity records name
  coordination and covalency alike, and one tying an ion to its surroundings
  made it unmatchable against a force field whose ion templates allow no
  external bonds.
- **A residue removed as a heterogen is not rebuilt.** It remained in the
  deposited sequence, so the gap it left was read as unresolved polymer and
  scheduled for reconstruction from a template that no longer existed.
- **Setup artifacts are not nested inside a second setup directory.**
- **A failed phase reports failure.** Setup swallowed any failure resolving its
  input and returned success, so a mistyped PDB identifier produced an empty
  directory and exit code 0. An absent optional backend still degrades
  gracefully, because choosing not to install OpenMM is a choice rather than a
  failure.
- Single-phase commands report why they failed; only `explore` did.
- The run summary reflects the run. Options were read under their `explore`
  spelling only, so `fastmdx setup --ph 6.5` displayed the default instead,
  for all fifteen of them.
- `scipy` and `pillow` are declared rather than arriving through
  `scikit-learn` and `python-pptx`. A test now compares what the source
  imports against what is declared.
- `netcdf4` and `umap-learn` become optional extras, matching their guarded
  imports.

### Known limitations
- A coordinated metal ion kept by `--setup-heterogens auto` can fail the
  ligand clash check: a zinc sits about 2 A from its donor atom and closer to
  that residue's hydrogens, well inside the 1.5 A threshold. Lower
  `ligand_clash_threshold_nm`, or exclude the metal.
- Hydrogens added to a retrieved ligand are placed from its own geometry,
  without reference to the surrounding protein, so a tight pocket can produce
  an apparent clash.
- No named force field loads nucleic acid parameters, so DNA and RNA cannot
  be prepared.
- Non-standard residues such as modified cysteines are treated as components
  rather than being replaced with their standard equivalents.

### Packaging
- Available from conda-forge: `conda install -c conda-forge fastmdxplora`,
  including the small-molecule stack, so protein-ligand preparation works from
  a single install command.
- The conda recipe uses the v1 format.

## [2.1.0] — 2026-08-01

A graphical interface, publication-ready figures, and Python 3.13 support.
An exploration can now be designed, launched, watched, and reviewed from a
browser, and the figures it produces are drawn for print.

### Added
- **Graphical interface.** `fastmdx gui` opens a local browser interface with
  sections for building an exploration, starting it, watching it run, viewing
  the structure and trajectory, browsing analyses, and reaching the report.
- **Config files from the GUI.** The builder can save its selection as a
  FastMDXplora `.yml` file instead of starting a run, so an exploration
  designed in the browser can be submitted anywhere, including on a cluster
  where the GUI cannot run. The GUI guide documents the rsync pattern for
  watching a cluster job from a workstation.
- **Live telemetry.** A dependency-free localhost server that watches an
  output directory while an exploration is in progress: phase progress, live
  trajectory frames, an interactive 3D molecular viewer with playback, and
  per-analysis charts. It observes and starts explorations; it does not
  reimplement any phase science.
- **Report content in the GUI.** Summary cards, per-phase progress including
  phases that did not run, trajectory statistics with means and standard
  deviations, categorised analysis sections, and quick-action links. These are
  computed once and shown identically in the browser and in the report.
- **`v1` analysis compatibility profile** (`--compat v1`): reproduces the
  scope, selection, stride, analysis set, and per-analysis options of the
  published BPTI case study from FastMDXplora version 1.
- **PDB smoke campaign** (`scripts/run_pdb_smoke_campaign.py`) for exercising
  the pipeline across many structures.
- **Report additions:** region highlights and a single-figure run summary.
- **Python 3.13 support** across the whole supported range (3.9 to 3.13).
- **Documentation:** a beginner's guide, a CLI reference, a GUI guide, a
  production-run guide, region-highlight and smoke-campaign pages, and a
  substantially expanded installation guide.

### Changed
- **Publication-ready figures.** The palette is now Okabe-Ito, which stays
  distinguishable under common colour vision deficiencies and in greyscale;
  the previous palette's red and green converged in both. Axes are closed with
  inward major and minor ticks, and tick density adapts to each panel's
  physical size, so long trajectories no longer crowd their axes. Legends get
  a translucent backing and headroom so they stay readable over dense data.
- **One interface module.** All user-interface code now lives in
  `fastmdxplora.gui`: the localhost server, the browser application and its
  assets, and the static dashboard written into a report. Both surfaces share
  one set of design tokens and display the same figures.
- **Explorer nomenclature throughout.** An exploration is started rather than a
  job launched: `Start Exploration` in the interface, `/api/explore/*`
  endpoints, `fastmdxplora.gui.exploration` in the Python API, and exploration
  wording in the terminal.
- **Each analysis has its own section.** Sections used to be merged when they
  held three figures or fewer, so an analysis appeared under its own name or a
  catch-all depending on how many plots that run produced.
- The startup banner is one left-aligned block; `info` carries the version,
  authors, DOI, phase availability, and backend status.
- CI runs on Python 3.9 through 3.13 across Linux, macOS, and Windows, with
  actions updated to Node 24 compatible versions.
- `STRUCTURE.md` rewritten to match the current tree.

### Fixed
- Simulation robustness, with clearer diagnostics when an integration step
  produces a NaN.
- Batch early-stop behaviour when an exploration fails.
- Report-only invocation and phase validation.
- Windows: path comparison, server port reuse, and platform-specific commands.
- The CLI degrades cleanly when the chemistry backends are absent rather than
  failing mid-phase.
- The startup banner printed twice, and advertised a GUI address even when no
  GUI was running.
- Figures were listed twice, once for the PNG and once for the SVG of the same
  plot.
- CLI tests asserted the exit code of a machine without OpenMM, so they passed
  in CI and failed on a working install.
- Documentation examples are covered by drift tests so they cannot go stale.

### Removed
- The bundled installer and repository doctor. Installation now follows
  standard routes: pip for analysis and reporting, pip plus conda-forge for the
  full chemistry stack, or a clone with the bundled `environment.yml`.
- The `dashboard` subcommand; `fastmdx gui --output DIR` opens the same
  interface pointed at an existing run.
- The curated chart pipeline, which redrew 17 figures in a second style that
  nothing displayed once both surfaces settled on the analysis figures.
- The live pressure metric. OpenMM's `StateDataReporter` cannot supply it, so
  the card could only ever read zero; the barostat setpoint remains in the run
  summary.
- `driftmd_workbench`, a separate package that duplicated the analysis and
  report pipeline without using it.

## [2.0.0] — 2026-05-25

**FastMDXplora** — Fully Automated SysTem for Molecular Dynamics eXploration.
A single command takes a structure (and optional bound ligand) from input to
publication-quality deliverable across four phases: setup, simulation
(including enhanced sampling), analysis (protein and protein-ligand), and
reporting.

### Packaging
- Canonical package: **`fastmdxplora`** (`import fastmdxplora`). The CLI
  command is **`fastmdx`**.
- `fastmdx` remains available on PyPI as a short alias that installs and
  re-exports `fastmdxplora`.

### Features
- End-to-end MD orchestration across four phases (setup, simulation, analysis, report).
- Named force-field selector; OpenFF ligand/cofactor parameterization with a setup-time pose clash check.
- Protein-ligand analyses: ligand pose RMSD, protein-ligand contacts + binding-site fingerprint, protein-ligand H-bonds, ligand RMSF — auto-detected for complexes.
- Analysis scope (`solute`/`protein`/`ligand`/`all`) keeping analyses off solvent.
- PLUMED enhanced sampling on the production stage (`--simulate-plumed-script`).
- Cross-platform CI (Linux/macOS/Windows); parallel batch execution and cross-run comparison.

## [0.3.0] — 2026-05-25

Enhanced sampling. FastMDXplora can now drive PLUMED collective-variable
biasing (metadynamics, umbrella sampling, steered MD, …) on the production
stage of a run, with equilibration left unbiased per standard protocol.

### Added
- **PLUMED enhanced sampling** (optional): supply a PLUMED script via `simulation.plumed` (config: `{enabled: true, script: "<inline or path to .dat>"}`) or `--simulate-plumed-script PATH` (CLI) to add collective-variable biasing — metadynamics, umbrella sampling, steered MD, etc. — to the **production** stage. Equilibration (NVT/NPT) runs unbiased, matching standard enhanced-sampling protocol; the biasing force is added just before production and the context reinitialized. PLUMED output files (COLVAR, HILLS, …) are redirected into the run's output directory, and the resolved script is saved as `plumed.dat` for reproducibility. Requires the `plumed` extra (`openmm-plumed`, installed via `conda install -c conda-forge openmm-plumed`); absent, enabling PLUMED raises a clear, actionable error.

### Changed
- Renamed the `test_md_parity.py` test module to `test_md_engine_controls.py` to match its content (MD engine controls). Trimmed the README (removed the Status and Project-family sections).

## [0.2.0] — 2026-05-25

End-to-end protein-ligand molecular dynamics. FastMDXplora can now set up,
simulate, and analyze a protein-ligand complex from a feasible bound pose:
named force fields with an OpenFF small-molecule path, a setup-time pose
sanity check, and the standard protein-ligand analysis suite, all detected
and wired automatically.

### Added
- **Protein-ligand analyses** (run automatically when a ligand is detected; `include`/`exclude` apply): in addition to ligand pose RMSD, three more commonly-reported analyses now run on protein-ligand complexes:
  - `contacts` — protein-ligand contacts, reported two ways: a per-frame count of protein residues within a cutoff (default 0.4 nm) of the ligand (`contacts.dat`), and a per-residue contact-frequency "interaction fingerprint" identifying the binding-site residues (`contacts_per_residue.csv`, also shown as the figure)
  - `pl_hbonds` — hydrogen bonds formed specifically between protein and ligand (per frame), distinct from the general intra-solute `hbonds` analysis
  - `ligand_rmsf` — per-ligand-atom fluctuation after protein alignment: the ligand's internal flexibility in the pocket
- **Ligand pose RMSD** analysis (`ligand_rmsd`): the headline protein-ligand stability metric. Each frame is rigidly aligned onto the reference using the protein (Cα by default), then RMSD is measured on the ligand atoms of the aligned coordinates — i.e. how far the ligand has moved *relative to the protein frame*, which tells you whether it holds its binding pose or drifts/unbinds. This is distinct from the standard RMSD (which aligns and measures on the same atoms). It runs automatically when a ligand is detected (from `resolved_forcefield.ligand` in the setup manifest) and is skipped for protein-only runs; `include`/`exclude` still apply. Ligand-only analyses are marked with a `requires_ligand` flag on the analysis class, and the orchestrator supplies the detected ligand residue name automatically
- **Analysis scope** (`analysis.scope` / `--analyze-scope`): a single setting controls which atoms analyses operate on — `solute` (protein + ligand, the default), `protein`, `ligand`, or `all`. It resolves to a default atom selection applied to analyses that don't set their own (the solvent-blind ones: Rg, SASA, secondary structure, Q-value, hydrogen bonds), so they no longer run on solvent/ions by accident. Analyses with a meaningful own default (the Cα-based RMSD, RMSF, clustering, dimensionality reduction) keep it. An explicit per-analysis or orchestrator-wide `selection` still overrides the scope. When a ligand is present (detected from `resolved_forcefield.ligand` in the setup manifest), `solute` and `ligand` scopes include it automatically by residue name
- **Ligand / cofactor parameterization** (protein-ligand systems): supply a small-molecule ligand as an SDF or MOL2 file via `setup.ligand` (config) or `--setup-ligand` (CLI), parameterized with an OpenFF small-molecule force field through `openmmforcefields`' `SystemGenerator`. Selected with the ligand-capable `amber-openff` named force field (AMBER ff14SB protein + TIP3P water + OpenFF Sage 2.2.1 for the ligand). Net charge is inferred from the SDF formal charges unless set explicitly via `ligand_net_charge`; the ligand residue name (`ligand_name`, default `LIG`) and small-molecule force field (`ligand_forcefield`, e.g. `openff-2.2.1` or `gaff-2.2.20`) are configurable. The supplied ligand coordinates must be a feasible bound pose (from a co-crystal structure or docking); a setup-time clash check (`check_ligand_clashes`, `ligand_clash_threshold_nm`) fails with a clear message if the pose severely overlaps the protein, rather than letting it surface as a divergent simulation later. Incoherent combinations are rejected early with clear errors (a ligand with a non-ligand-capable force field, or with a raw XML list). The resolved ligand parameterization is recorded under `resolved_forcefield.ligand` in `setup_parameters.json`. Requires the `ligand` extra (`pip install 'fastmdxplora[ligand]'`); absent, the phase degrades with an actionable install message. Ligand input is list-shaped in config for future multi-ligand support; single-ligand parameterization is implemented now
- **Named force-field selector**: pick a force field by a short, documented name via `setup.forcefield` (config) or `--setup-forcefield` (CLI) — `charmm36` (default), `amber14`, `amber-fb15`, or `amber-openff` (ligand-capable) — instead of listing raw OpenMM XML filenames. Each name resolves to the correct protein/water XML set and default water model through a single registry (`setup/forcefields.py`). The raw `force_field` XML list remains as a power-user escape hatch; specifying both a named selector and a raw list is rejected with a clear error, as is an unknown force-field name (the message lists valid choices). The resolved force field (actual XMLs + water model) is recorded under `resolved_forcefield` in `setup_parameters.json` for reproducibility, regardless of which form the user chose

### Fixed
- The ligand residue in the prepared/solvated topology is now named with the configured ligand name (default `LIG`) instead of OpenFF's default `UNK`. Previously the written `topology.pdb` labelled the ligand `UNK` while the manifest recorded `LIG`, so resname-based selection silently failed — ligand-aware analyses found no ligand atoms, and the `solute`/`ligand` analysis scopes silently excluded the ligand. The name is now set on both the ligand topology and the merged topology so it survives `Modeller.add()` across OpenMM versions
- Clustering on a trajectory with fewer frames than the requested number of clusters now fails with a clear, actionable message ("Clustering needs at least n_clusters=N frames, but the trajectory has only M...") instead of an opaque scikit-learn internals error. k-means and hierarchical clustering are guarded; DBSCAN (which doesn't take a cluster count) is unaffected
- Analyses that operate on all atoms by default (Rg, SASA, secondary structure, Q-value, hydrogen bonds) now slice the trajectory to the resolved scope/selection *before* computing, rather than processing the full solvated system. Previously several of these passed the whole trajectory straight to the underlying calculation regardless of the selection — so on a solvated complex the Q-value analysis enumerated residue pairs across ~10k water residues (tens of millions of pairs) and effectively hung, and Rg/SASA were computed over water. With the new `solute` default scope and per-analysis slicing they operate on protein (+ ligand) only — a correctness fix for any solvated run and a large speedup
- **Named force-field selector**: pick a force field by a short, documented name via `setup.forcefield` (config) or `--setup-forcefield` (CLI) — `charmm36` (default), `amber14`, `amber-fb15`, or `amber-openff` (ligand-capable) — instead of listing raw OpenMM XML filenames. Each name resolves to the correct protein/water XML set and default water model through a single registry (`setup/forcefields.py`). The raw `force_field` XML list remains as a power-user escape hatch; specifying both a named selector and a raw list is rejected with a clear error, as is an unknown force-field name (the message lists valid choices). The resolved force field (actual XMLs + water model) is recorded under `resolved_forcefield` in `setup_parameters.json` for reproducibility, regardless of which form the user chose

## [0.1.0] — 2026-05-XX

Initial claim-staking release. Establishes the project-level orchestrator
scaffolding, the four-phase API (setup, simulation, analysis, report), and
the `fastmdx` CLI.

### Added
- **Robust auto platform selection**: when `platform=auto`, the simulation runner now verifies a GPU platform (CUDA/OpenCL) can actually create a Context before committing to it, and falls back to the next candidate (ultimately CPU) if not — instead of selecting a *registered-but-unusable* platform that then fails at Context construction with a confusing error. An explicit `platform=CUDA`/`OpenCL` request is still honored as-is (the user sees the real error if their choice is broken)
- **Clear periodic-box / cutoff guard**: `prepare_system` now raises an actionable error when the nonbonded cutoff exceeds half the smallest periodic box dimension (instead of OpenMM's cryptic `NonbondedForce` message), naming the cutoff, the box, and how to fix it (increase `solvent_padding_nm` or decrease `nonbonded_cutoff_nm`)
- **`environment.yml` + git install path** — clone the repo and `mamba env create -f environment.yml || conda env create -f environment.yml` then `pip install .` to get all four phases (the OpenMM/PDBFixer chemistry stack from conda-forge) without waiting on the conda-forge package. Plain `pip install fastmdxplora` still gives the analysis + report phases on their own

### Fixed
- **Parallel execution on Windows**: spawned worker processes now reconfigure their stdout/stderr to UTF-8 (as the CLI entry point does). Previously, because workers are spawned (not forked) on Windows and bypass the CLI entry, their streams stayed on the platform codec (cp1252) and crashed with `UnicodeEncodeError` the moment the presenter printed a status glyph (✓, ▸) — so every run in `mode: parallel` failed on Windows while sequential mode succeeded
- **Headless plotting**: the analysis package now forces matplotlib's non-interactive `Agg` backend before pyplot is imported (respecting an explicit `MPLBACKEND`). Previously the backend was only forced off-Windows with no `DISPLAY`, so analyses crashed on headless machines that didn't match that gate (notably headless Windows CI) with "Can't find a usable init.tcl". FastMDXplora always writes figures to files, so a non-interactive backend is always correct
- **Cross-platform paths in reports/manifests**: figure links in the Markdown report, zip archive entry names, and the relative artifact paths recorded in manifests now use forward slashes (`as_posix()`) on every OS. Previously, on Windows these were emitted with backslashes, breaking Markdown/HTML image links and producing non-portable manifests
- **UTF-8 file/stream encoding everywhere**: all text written by FastMDXplora (reports, comparison markdown, config templates, manifests, PDB/XML artifacts) now specifies `encoding="utf-8"` explicitly, and the CLI reconfigures stdout/stderr to UTF-8 at entry. Previously, on a machine whose default locale encoding was ASCII, writing the comparison report's `→` or the config template's `—` (or printing the banner) raised `UnicodeEncodeError`
- **FastMDXplora orchestrator class** (`fastmdxplora.FastMDXplora`) — project-level coordinator following a seven-phase orchestration pattern (Aina & Kwan, JCC 2026)
- **Four phases** under `fastmdxplora.setup`, `.simulation`, `.analysis`, `.report` — each with a `run(orchestrator, output_dir, **options)` entry point and a structured parameters manifest
- **`fastmdx` CLI** with subcommands `explore` (canonical), `xplore` (X-themed alias), `setup`, `simulate`, `analyze`, `report`, `info`, plus `--version` and `--cite` flags
- **Report phase artifacts**: Markdown study report, .pptx slide deck, self-contained .zip project bundle
- **YAML configuration files**: a single config captures an entire study — input is given as a canonical `systems:` list (always a list, even for one system), plus phase selection and all per-phase options; drives both the CLI (`--config` / `-c`) and the Python API (`FastMDXplora(config=...).explore()`); `fastmdx init-config` writes a fully-commented template; strict schema validation rejects typos with did-you-mean suggestions; command-line flags override file values; every run writes a re-runnable `resolved_config.yml` for reproducibility
- **Full MD engine controls**: integrator selection (`langevin_middle`, `langevin`, `brownian`, `verlet`, `variable_langevin`, `variable_verlet`); pressure in either `pressure_bar` or `pressure_atm` (auto-converted); GPU `device_index` selection; `checkpoint_interval_steps` writing a restart-ready `.chk`; `ForceField.createSystem` pass-throughs (`nonbonded_method`, `ewald_error_tolerance`, `use_switching_function`, `switch_distance_nm`, `dispersion_correction`, `remove_cm_motion`); and `fixed_pdb` to skip PDBFixer when a prepared structure is supplied
- **Many-system & parameter-sweep mode**: the `systems:` list can hold several systems, and an optional `sweep:` of parameter axes (dotted `phase.option` keys) runs the full cross-product (systems × sweep), each as a complete self-contained study. One run writes the flat output layout; multiple runs go in `runs/<id>/` indexed by a top-level `batch_manifest.json`. Per-system option overrides and swept values merge with correct precedence (base < per-system < sweep); typo'd sweep axes are rejected with the valid-option list. An optional `execution:` block runs studies in parallel (process pool) with round-robin GPU device pinning (one run per device)
- **Cross-run comparison report**: after a multi-run study, a `comparison/` report is built automatically at the batch root — per-frame **overlays** (RMSD, Rg, Q-value, total SASA across all runs on one axes), **trend** plots of each run's summary scalar against the swept parameter, a `comparison_summary.csv`, and a written `comparison_report.md` with a quantitative takeaway per property. Degrades gracefully (errored runs / missing analyses skipped); disable with `report: { comparison: false }`; (re)build via `FastMDXplora(...).compare()` (optionally `compare(output_dir=…)` for a batch that finished earlier)
- **Dry-run / plan-only mode**: `fastmdx explore --config … --dry-run` (or `explore(dry_run=True)`) prints every run, its system, swept values, target output directory, and the phases that would execute — then exits without running anything or writing to disk
- **Uniform return shape**: `FastMDXplora.explore()` always returns a `list[RunResult]` — a single study is a list of one, a sweep is a list of many. Each `RunResult` carries `run_id`, `system`, `status`, `output_dir`, `sweep_values`, and its per-phase `PhaseResult` list in `.phases` (with a `.phase(name)` lookup helper). The single user-facing entry point is always `FastMDXplora`; the batch machinery underneath is private
- **Reproducibility manifest** (`manifest.json`) written at the project root summarizing phases executed, parameters, software versions, and DOI
- **Datasets namespace** (`fastmdxplora.datasets`) with a TrpCage placeholder
- **CI**: matrix tests on ubuntu/macos/windows × Python 3.9–3.12 (GitHub Actions)
- **PyPI**: dual-name publishing — `fastmdxplora` is the primary package, `fastmdx` is a thin alias that depends on it

### Notes
- The analysis and report phases are self-contained (no heavy runtime dependencies); the setup and simulation phases require OpenMM + PDBFixer.
