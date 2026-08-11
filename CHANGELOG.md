# Changelog

All notable changes to FastMDXplora are documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) ·
Versioning: [SemVer 2.0.0](https://semver.org/spec/v2.0.0.html)

## [Unreleased]

## [2.5.0] — 2026-08-10

**Recompute any metadynamics or umbrella result on a torsion or an angle.**
Separations on a periodic coordinate are now measured the short way round, in
the umbrella window stitching, the free energy surface, the bias each frame
felt, and the offset that corrects for it. Results on a distance or a radius
of gyration are unaffected, as are torsion studies confined to part of a turn.

Metadynamics recrossings are counted as travel between the two deepest basins
rather than past a fixed threshold, and the record names the basins.

The nonbonded cutoff comes from the force field: CHARMM36 at 1.2 nm with
switching from 1.0, AMBER hard-truncated at 1.0 with no switch. An explicit
setting still wins.

`workers` implies parallel execution and each worker takes a share of the
machine's cores. Capping groups survive setup and count as polymer.
`--force-overwrite` replaces `--force`, which still works. Analyses that
superpose are skipped where fewer than three atoms match. PLUMED's setup goes
to `simulation/plumed.log`. A study prints a progress bar.

Documentation leads with the config: one file describes a study, and the GUI,
the CLI and the Python API each build one and each run all four phases.


### Added
- **A metal in a protein site says the force field will not hold it there.**
  Standard force fields treat a metal ion as a point charge with
  Lennard-Jones terms and nothing else. That is enough for a carboxylate cage
  and often not for a metal held by histidines, so the ion drifts out of its
  site while the fold stays intact, the RMSD stays flat, and nothing else in
  the run reports it. Measured on thermolysin over 20 ns: all four calciums
  held their sites to within a tenth of an Angstrom, and the catalytic zinc
  lost His142 in the first production frame and never recovered it. Setup now
  measures which metals are actually coordinated -- rather than assuming from
  the residue name -- and names them, with the remedies. Salt stays silent.

- **The config is the study, said where somebody looks.** The README and the
  documentation index now lead with it: one file describes a study completely,
  and the GUI, the command line and the Python API each build one and each run
  all four phases. None is the primary interface and none is a subset -- they
  are generated from a single declaration, which 53 parity tests already
  enforced and nothing stated. The GUI page also understated itself: 159
  options across 15 analyses, when there are 200 across 19, and no mention of
  the 100 phase settings at all.

- **PLUMED's setup is kept beside the run.** Forty lines -- which atoms the
  collective variable is built from, the hill width, the pace, the bias
  factor -- written from C++ straight to the file descriptor, arriving in the
  middle of a progress bar. It goes to `simulation/plumed.log`: it is the only
  independent statement of what PLUMED actually did, and reading `TORSION
  between atoms 7 9 15 17` in it is how a psi selection was confirmed correct.

- **A study reads its own curve.** `pmf.json` carried a free energy and left
  the reading of it to whoever opened the file, which is how it got misread.
  The grid spans wherever the coordinate went, the windows covered a part of
  it, and a minimum taken across the whole grid references a region nothing
  visited -- a real study looked to have a 164 kJ/mol barrier that way,
  against 11 measured where the windows actually were. It now records the
  covered range and a summary: the barrier, where it sits, and every minimum
  sorted by depth. The `pmf` analysis reports them instead of only drawing a
  line.

- **A closed turn measures itself.** The two ends of a full period are the
  same place, so their free energies must agree. They arrive there from
  different windows through a chain of joins and nothing forces them to, which
  makes the difference a measurement rather than a fault: it is what the
  study's statistics are worth. Reported, not constrained -- forcing the ends
  together would make the number zero and take the information with it. A
  well-sampled synthetic profile closes to nothing on its own; a real study of
  0.2 ns windows closed to 2 kJ/mol, which is that study's uncertainty stated
  honestly. A partial turn carries `None` rather than zero, because "this was
  never a circle" and "this closed perfectly" are opposite claims.

- **The simulation phase now says why each stage exists.** `explain.py` opens
  by saying a pipeline that does all this silently "teaches nothing", and that
  held for setup and not for simulation: the entries for minimisation, NVT,
  NPT, production and the ensemble were written, and none of them could be
  reached. The runner announces through a callback carrying a message and
  nothing else, and the presenter could print an explanation from `step` but
  not from `info`, which is what stages use. Both ends existed and were not
  joined. A real run explained its protonation and its solvation, then went
  quiet for the twenty minutes in which the ensemble is chosen.

- **What NVT and NPT are for, and what choosing one does to production.**
  Both steps had an explanation and neither had a reference. They now cite
  the LiveCoMS best-practices guide and the Monte Carlo barostat OpenMM
  actually uses, and the NPT text carries the measured number rather than
  calling the error "usually a little wrong": solvation packs a box about ten
  per cent short of water, and only a barostat corrects it.

  A third explanation covers the question neither answered -- why you would
  run only one. NVT production is legitimate at a density you know is right,
  and the way you learn that density is to run NPT first and take the average
  box size from it. Going straight to NVT is a different thing entirely: it
  fixes the box at whatever solvation produced. That is the arrangement two
  runs here used, at 0.92 g/mL, and it is the one nobody intends.

- **The README says what a biased run's averages mean.** The front page
  listed the three enhanced sampling methods and what each output is and is
  not, and stopped there -- so the most distinctive thing the analyses now do,
  correcting those averages back to the ensemble somebody actually wanted, was
  not on it. The closing paragraph now names the harder case too: a run that
  does not support the thing it was for is told so in those terms.

- **The reweighting is documented.** What the analysis phase writes to
  `analysis/reweighted/`, the estimator it uses and why it carries the c(t)
  offset, which analyses are corrected and which are not, and why an umbrella
  window and a steered pull cannot be. The two details that are easy to get
  wrong are stated rather than left implicit: weighting by the final surface
  inflates early frames, and without c(t) the weights rank frames by when they
  were written instead of by where the system was. Pinned by tests, so an
  analysis that gains a correction and does not gain a sentence fails.

- **Each enhanced-sampling method reports the result it exists to produce.**
  A metadynamics run wrote its free energy surface to JSON and stopped there:
  no figure, no entry in the analysis manifest, no mention in the report.
  Sixteen analyses of the trajectory each produced a curve and a plot, and the
  one result the run was for did not. The same gap held for an umbrella
  study's potential of mean force and a steered pull's work.

  There are now three analyses -- `pmf`, `metad_surface` and `steered_work` --
  each reading what the biased run itself recorded rather than recomputing it,
  each drawing it, and each running only where such a run produced one. They
  are not analyses of the trajectory and do not need a ligand.

- **A biased trajectory says what it is.** The methods section described
  restraints well and said nothing at all about the three methods where
  biasing is the entire point. Ten analyses were reported beside the free
  energy with no distinction between them, and a reader would take a mean
  RMSD over a metadynamics run as a measurement of the system.

  The three are not alike, and one caveat covering all of them would be wrong
  about two. Metadynamics deposits a known bias, so the unbiased ensemble is
  recoverable by weighting. An umbrella window describes a system held where
  it was put, and what combines the windows is the free energy rather than an
  average across them. A steered pull is not an equilibrium ensemble at all.
  Each now gets the paragraph that is true of it.

- **The restraint ladder steps across equilibration rather than at its
  boundaries.** The strength was sampled at two points, before NVT and before
  NPT, so a four-rung ladder reached the first and third rungs and never the
  others. With `npt_steps: 0` the second sample sat inside a branch that never
  ran, and the restraint held at full strength through all of equilibration
  and dropped to zero at production -- the release all at once that a ladder
  exists to prevent, from a setting the user had written out in four steps.

- **The GUI cites the software.** The report, the slides and `fastmdx info`
  all print the citation and the GUI printed none, which made the interface
  the documentation sends a new user to first the one that never said how to
  cite the work. There is a Cite page with the reference, the DOI and the
  BibTeX entry, and the report's own dashboard carries it in its footer.

  Filled by the server from the same constants, because writing it into the
  page would have put a fourth copy beside the report's, the slides' and the
  CLI's. There were already two copies of the BibTeX entry; it is now declared
  once beside the citation it belongs to, and a test fails if it is written
  out anywhere else.

- **A progress bar during molecular dynamics.** The part of a run that takes
  the time announced how many steps it was about to take and then said nothing
  until the stage ended -- half an hour of a terminal that looked exactly like
  a hung one. Each stage now steps in chunks of a fiftieth and reports how far
  through it is, at what rate in ns/day, and how long is left, so the decision
  to wait or to stop can be made in the first few seconds rather than at the
  end.

  On a terminal it is one line, rewritten. Where the output is a file -- which
  is where a worker's output goes when several run at once -- a line of
  carriage returns is unreadable, so it prints at tenths instead.

- **A mean says how much of a run is behind it.** Ten analyses averaged over
  the whole production run without asking whether the system had settled by
  the time the averaging started, or how many *independent* observations the
  average rested on. Frames are not independent: a trajectory written every
  picosecond from a system decorrelating over a hundred has a hundred times
  fewer independent samples than frames, and an error computed as though each
  frame counted is understated -- six-fold on a test series, which is how a
  difference between two systems becomes significant on paper without being
  real.

  Both questions have one answer. The statistical inefficiency is the number
  of frames per independent sample, so it fixes where the relaxation ended
  (Chodera, J. Chem. Theory Comput. 2016, 12, 1799) and what the mean is
  worth. Below ten independent samples the mean describes the run rather than
  the system, and it is refused.

  Recorded by the base class for any analysis declaring that it produces one
  value per frame -- six of the sixteen -- so `findings` carries the mean, its
  error, the frames discarded and the independent samples behind it. Declared
  rather than inferred from the array's length: a per-atom result on a
  trajectory with as many frames as atoms would otherwise be summarised as a
  time series, and the numbers would look right.

  A run too short to measure its own correlation time is refused separately,
  because that failure flatters: the inefficiency comes back too small and the
  sample count too large. On a test series with a true inefficiency of 2000,
  four thousand frames gave 361 and an apparent eleven independent samples
  where the truth was two.

- **A metadynamics run produces a free energy surface, or refuses one.** It
  wrote PLUMED's HILLS and COLVAR and stopped: nothing read them back, so the
  module said "a run that has not converged has no free energy" while offering
  no free energy either way, and whoever wanted a surface summed the hills
  themselves and got one with nothing attached. Umbrella sampling went from a
  config block to a curve or a refusal; this went to a pair of files.

  Three conditions, each answerable from what the run already writes: the
  hills must have decayed, so the bias has flattened rather than still filling
  the landscape; the surface built from three quarters of the hills must match
  the one built from all of them; and the system must have crossed between the
  ends of the range more than once, because a barrier crossed once has been
  observed once. The evidence is written down beside the verdict either way,
  so a borderline run can be judged rather than only rejected.

- **A run records which code made it.** The manifest carried the version
  string and nothing else, and setuptools-scm writes that at install time --
  so an editable install carries whatever it was when `pip install -e .` was
  last run. A study of ours came back stamped `2.3.0` for seven windows that
  used shared setup, a feature `2.3.0` did not have: the manifest named a
  version in which the run could not have happened, and it is the number the
  report's reproducibility section prints.

  Where the package was imported from a checkout, the commit is recorded
  beside the version, and so is whether the tree had uncommitted changes. The
  checkout is found from the package rather than the working directory,
  because a run started inside another repository would otherwise record that
  repository's commit -- a precise claim about the wrong code, which is worse
  than recording none. A dirty tree is reported rather than refused: refusing
  would block the runs developers make all day, and a commit beside
  uncommitted changes does not describe what ran, so saying so is what keeps
  the commit from being decorative.

- **The banner says how to watch the run.** It computed the dashboard address
  and discarded it, on the reasoning that the GUI prints its own -- true while
  a server is running, and no help to somebody who has just started a run.
  Where none is running it prints the command that starts one, pointed at this
  run's output; where one is, the address.

### Changed
- **One page for one run.** Overview and Live Simulation showed the same run,
  and the top bar showed it a third time: three places to look for one answer,
  and neither page complete on its own. They are one page now, with what is
  happening at the top -- progress, charts, health, the structure -- and what
  has been recorded beneath it. The panels were moved rather than rebuilt, so
  every element the script reads is where it was. Two cards were both called
  progress; one is a bar for the running stage and the other a table of
  phases, and they now say so.

- **Watching a run is on by default.** The live dashboard could show nothing
  unless somebody knew to turn on a setting called telemetry -- a word that
  elsewhere means data sent to a vendor, and here means four files written
  into the run directory for the local page to read. It cost a tenth of a per
  cent, measured on a solvated system over 2,000 steps with and without, and
  the frame history is capped, so the only thing the default saved was the
  ability to watch.

- **The GUI asks to be cited where somebody will see it.** The citation page
  existed and was the last item in the sidebar, beneath Documentation and
  GitHub under a heading reading Tools -- so the one thing a scientific tool
  most needs its user to find sat after the two links that navigate away from
  it. It comes before them now, and the sidebar footer, which shows on every
  page, carries a line linking to it.

- **The interface is called the GUI.** Documentation and comments used "the
  browser" for both the interface and the thing it renders in, so the three
  interfaces read as "the command line, a config file, and the browser". A
  browser tab and the `--no-browser` flag are the literal thing and keep the
  word; everything else naming the interface says GUI.

### Removed
- **The `datasets` package.** It promised a Trp-cage trajectory that was never
  bundled: `TrpCage.traj` returned the path to a file that had never existed,
  from v0.1.0 through v2.4.0, so reading it gave a plausible string and
  passing it anywhere gave a file-not-found from inside a trajectory reader.
  Nothing replaces it. FastMDXplora fetches a deposited structure and
  simulates it, so a reference trajectory carried in the distribution would be
  megabytes of wheel for something one command produces --- and a reference
  system is a PDB identifier, not a file.

### Fixed
- **The reported free energy surface was summed the long way round.** The
  periodic separation was added to `surface_from_hills` and the only two
  places that call it were left on the default, so the number every run
  reported came from the arithmetic that had just been fixed. A converged
  10 ns metadynamics run on alanine dipeptide's psi reported 60.7 kJ/mol; the
  same hills give 25.8. Every report, dashboard, slide deck and bundle
  downstream carried the wrong figure, and the analysis phase does not
  recompute it -- it reads what the simulation phase wrote, so re-running
  analysis reproduced the error with a fresh timestamp.

- **A recrossing was counted past a line rather than between states.** The
  thresholds sat a quarter of the way in from the extremes of the hill range,
  which describes the grid and not the system. On a full turn that put them at
  -90 and +90 degrees, and a run whose minima were at -17 and 155 had one of
  them inside the dead band. The count is now travel between the two deepest
  basins, each frame assigned to the nearer one measured the short way round.
  The same run: 925 became 349, against a threshold of 4.

- **Two more analyses fitted a rotation to one point.** `cluster` and `dimred`
  superpose on `name CA` exactly as `rmsd` and `rmsf` do, and were not gated
  when those were. A capped alanine was clustered on identity rotations,
  sklearn found one distinct cluster where five were asked for, and the run
  reported `ok`. The test now asks the source which analyses superpose rather
  than trusting a list.

- **`precision` was reported as applied when the platform has no such
  setting.** OpenMM's CPU platform offers `Threads` and `DeterministicForces`
  and nothing else. The banner now says what was applied, and the diagnosis
  after a non-finite coordinate no longer suggests double precision there.

- **A path was handed back with its separators changed.** The directory the
  banner suggests watching was split with `pathlib` and rebuilt, which
  normalises separators, so a Windows caller who wrote `../runs/study` was
  given a correct path that was not the one they typed. Found by Windows CI,
  on a function three green macOS runs had passed.

- **Every force field got the same cutoff, and two of the four were wrong for
  it.** The default was 1.0 nm with switching from 0.9, applied to all of
  them. CHARMM36 is developed at 1.2 nm with switching from 1.0, so it was run
  0.2 nm short with the switch in the wrong place; the AMBER force fields are
  developed with hard truncation and were being switched, which moves a run
  away from the parameterisation rather than towards it.

  This is not a preference. A force field is fitted with a particular
  treatment of the truncation and the effects of that truncation are
  compensated in its other parameters, so the scheme belongs to the force
  field rather than beside it -- and CHARMM36 at AMBER's cutoff is a different
  force field from the one that was validated. Each registry entry now carries
  its own, an explicit `nonbonded_cutoff_nm` or `switch_distance_nm` still
  wins, and the run says which scheme it used and why.

  What OpenMM applies is the potential-based switching function; CHARMM's is
  force-based, and the toolkit does not have it. Lee et al. say so directly in
  the CHARMM-GUI Input Generator paper and prescribe this protocol for OpenMM
  regardless, having tested a range of cutoff schemes against CHARMM's own
  results. Naming it "CHARMM's switching" would claim a function that is not
  there.

- **A torsion study that covered the whole turn returned a ramp.** The bias
  on a window was `0.5 * k * (x - centre)**2`, a straight-line subtraction. On
  a circle that is wrong at the wrap: a sample at +170 degrees is ten degrees
  from a window held at -180, not three hundred and fifty. The code charged it
  933 kJ/mol instead of 0.76 -- a Boltzmann weight wrong by a factor of
  10^162 -- so that window's free energy was pushed up and every join after it
  inherited the error.

  Twelve windows tiling a full turn of alanine dipeptide's psi returned a
  monotonic slide of 180 kJ/mol with no minimum anywhere in it, and every
  check passed: twelve runs finished, no unsampled bins, every overlap above
  the threshold, no refusal. A tool that refuses when it cannot support a
  claim produced a confident wrong one. Corrected, the same windows give two
  minima at -15 and 159 degrees and two barriers between them, 10.4 kJ/mol one
  way round and 23.7 the other.

  The scope is worth stating precisely: a study confined to part of the turn
  is unaffected. The nine-window study before this one kept its windows away
  from the wrap and gave 11.0 kJ/mol before the fix and 11.5 after -- the
  difference being the periodic distance now applied throughout, where
  subtraction had been merely close. It also found the easier of the two
  barriers by luck of which arc it happened to cover.

- **Asking for workers did not ask for parallelism.** `mode` defaulted to
  sequential and `workers` was read separately, so `workers: 3` ran the runs
  one at a time and said nothing about it. Nobody sets a worker count wanting
  one at a time, and nobody lists GPUs to leave all but one idle. `mode`
  written out still wins.

- **Parallel workers competed for the whole machine.** OpenMM's CPU platform
  takes every core it can see, and a pool of workers each doing that
  oversubscribes by the worker count -- so `workers` could not be used well on
  CPU at all. Each worker is now given `cores // workers` threads, set in the
  worker process because the libraries underneath read it once at import.

- **Capping groups were stripped as heterogens.** ACE, NME and the rest
  terminate a chain and are part of the molecule. Removing them is right for
  buffer and cryoprotectant and wrong here: it left a bare alanine wearing
  atoms from its caps, and the run failed with "no template found for ALA",
  naming the residue that survived rather than the two that did not.

- **A rotation was fitted to one point.** RMSD and RMSF align on `name CA`,
  and a molecule with one alpha carbon has no superposition. MDTraj's C
  extension printed "UNCONVERGED ROTATION MATRIX. RETURNING IDENTITY" once per
  frame -- thousands of lines in a window's log -- and returned distances
  measured against no alignment at all, which look like results. Both now
  declare a minimum of three atoms and are skipped below it, the way a chain
  too short to have a fold already is.

- **`--force` did not say what it would force.** It is `--force-overwrite`
  now, with the old spelling kept as an alias so existing scripts still work.
  Three refusals still named the old one, which is worse than either name
  alone: it sends a reader to look for something the help no longer lists.
  A test scans for stragglers, because a rename is exactly the change that
  leaves messages behind.

- **A study printed a sentence a minute for an hour.** Nine windows on a
  laptop produced sixty near-identical heartbeat lines, burying the ones that
  mattered -- the windows finishing, and anything that went wrong. It draws a
  bar instead, redrawn in place where there is a terminal to redraw on and
  written out in full where the output is a file.

- **The banner said to watch a directory that would not change again.** An
  umbrella study prepares one system every window shares, and the preparation
  announced where it writes: accurate, and a second of setup finished before
  anybody could type the command. It points at the study now.

- **The NPT explanation gave one number where the effect is not one number.**
  It said solvation packs a box "roughly ten per cent short of water", from
  two runs that happened to sit at the same padding. Measured against OpenMM
  directly across box sizes, with and without a solute: pure water packs at
  0.96 to 0.99 whatever the size, a solute at 1.0 to 1.2 nm of padding brings
  it to about 0.90, and the same solvation at 2.0 nm reaches 0.96. The gap is
  the vacuum shell left around the solute, and it matters in proportion to
  how small the box is -- so it is described that way now, and the run
  reports its own number rather than a claim standing in for it.

- **The reconstructed bias was 11.1% too large on every well-tempered run.**
  PLUMED does not store the height it deposited. For a well-tempered run it
  stores that height multiplied by y/(y-1), so that summing HILLS gives the
  free energy directly -- which is the convention the free energy surface
  relies on and which stays. Reconstructing the *bias* needs it undone, and it
  was not: a run asking for 1.2 kJ/mol had 1.333 in HILLS, and every bias
  summed from that file was too large by the same ninth.

  It does not cancel. Both V and c(t) scale by the same factor, so their
  difference scales too, and the weights go as exp((V - c(t))/RT) -- a scaled
  exponent is an effective-temperature error, which sharpens the weights,
  understates the effective sample size and biases every average resting on
  them. Found by comparing against PLUMED's own record of the same quantity on
  a real 1L2Y run, which disagreed by 11.5% of the bias range while every unit
  test passed, because the test fixture wrote HILLS the convenient way rather
  than the way PLUMED writes it. The fixture now writes what PLUMED writes.

- **A hill was felt at the instant it was laid.** PLUMED prints the bias for
  a step before depositing that step's hill, and HILLS and COLVAR usually
  share a stride, so counting the hill as already felt was wrong on every
  row. On a real run it left the reconstruction out by 1.200 kJ/mol against
  PLUMED's own record -- exactly the configured hill height, which is what
  identified it once the larger height-convention error was removed. The
  frames and the c(t) checkpoints now take the same boundary, and the test
  fixture computes its reference bias the way PLUMED reports it rather than
  the way that was easier to write.

- **The collective variable was recorded as `cv`.** That is PLUMED's label for
  the column, not the name of what was biased, and it is what the provenance
  record carried. The configured name is used where it can be found, with the
  label kept beside it.

- **A run did not record its own box.** The setup record kept the atom count
  -- added because a methods section has to state it and it was only ever
  logged -- and not the periodic cell, which a methods section states for the
  same reason. Diagnosing a failed run meant reading CRYST1 out of
  `solvated.pdb` by hand to find out whether the box had ever been big enough
  for the cutoff. It now records the vectors, the perpendicular widths, the
  volume and the largest cutoff the box can carry. The widths matter rather
  than the edge lengths: for a rhombic dodecahedron the smallest width is the
  edge over root two, so reading the edge alone overstates the room by 40%.

- **The line announcing a bias said it was happening, not that it would.**
  All three methods print where their PLUMED script is written, which is
  before minimisation, in the present tense -- "Metadynamics biasing
  radius_of_gyration". Biasing is added just before production and
  equilibration runs unbiased, correctly, but the log read as though the bias
  were live from the first step. It sent the diagnosis of a real failed run
  off after the bias for a while when equilibration was the only thing that
  had run. All three now say which stage they apply to.

- **A NaN raised by the integrator never reached the diagnosis written for
  it.** The module that reads a failed state -- which atoms went non-finite,
  what residues they belong to, and what that points at -- was wired only to
  the check that runs at a stage boundary. OpenMM detects a non-finite
  coordinate during integration and throws, so that check never ran, and the
  common way a simulation dies fell through to the generic list of settings to
  try that the diagnosis exists to replace. A real 1L2Y run hit it after
  eighteen minutes and was told to lower the timestep, lower the temperature,
  raise the friction or turn off NPT, without anything knowing which applied.

  Both `step()` sites now recover the state from the context, which survives
  the exception and still holds the coordinates that went wrong. OpenMM's own
  message is kept beside the diagnosis rather than replaced by it: it says
  what the integrator noticed and links to its FAQ, while the diagnosis says
  which atoms it happened to.

- **An umbrella window and a steered pull reported their averages as though
  the run had been ordinary.** Metadynamics gained a reweighted column and an
  effective sample size beside every mean; the other two biased methods gained
  nothing, so their dashboard rows read `RMSD` exactly as a plain simulation's
  would. That is the worse case rather than the milder one: a window is held
  where it was put and a pull is not an equilibrium ensemble at all, so unlike
  metadynamics there is no corrected number to set beside the raw one.

  The analysis phase now records that a run was biased even where it cannot
  undo the bias, naming the method from the PLUMED script it wrote. The report
  leads with what the averages are not, the dashboard labels every metric as
  being of a biased ensemble, and neither invents a corrected column -- which
  would be the same numbers under a heading claiming otherwise. The wording is
  the methods section's own, so the two cannot drift apart.

- **Two config spellings with one meaning were refused.** A phase block
  present but empty -- `analysis:` with only a comment under it -- parsed as
  null and was rejected, while leaving the key out entirely was accepted:
  two spellings of the same intent behaving differently, when an option set to
  null is already read as "use the default". And `--include setup,simulation`
  reached the validator as a single string, because argparse takes a
  comma-separated list as one item. No phase or analysis name contains a
  comma, so neither reading is ambiguous. Both are now settled before
  validation, and a block of the wrong type -- a number, a list -- is still
  refused, because that is a real mistake rather than a spelling.

- **Metadynamics did not run at all.** The collective-variable plan was bound
  to the name the stage plan already used, so the next line to subscript it
  raised `'MetadynamicsPlan' object is not subscriptable`. The feature failed
  on the first line of use. Thirty-eight tests covered the module and
  twenty-three the surface; none of them ran the runner, which is where the
  two names met. A test now runs it.

- **A steered pull's work record was lost to a race.** PLUMED buffers its
  output and writes on teardown, and the record was built at the end of the
  simulation phase -- before the force was finalised. `COLVAR` existed, held
  no rows, and the record was silently not written; run by hand afterwards on
  the same file it worked. The pull now flushes as it goes.

- **A metadynamics refusal was cut mid-sentence.** Clipped at 160 characters,
  it lost the clause saying the output is still a usable snapshot of the
  filling. This is the one message whose only job is to explain why there is
  no surface, and a refusal that runs long is a refusal with something to say.

- **The surface's drift was judged where it means least.** Convergence was
  assessed over the whole grid, including regions far above the minimum that
  are visited rarely and move by several kJ/mol however long the run. Drift is
  now judged within 20 kJ/mol of the minimum, with the whole-grid figure still
  reported beside it.

- **A solvated tripeptide passed the fold gate.** The check counted every
  residue in the box -- 529 for a tripeptide in water -- so a chain far too
  short to have a fold was let through, and the Q-value analysis then sliced
  to the solute, found three residues, and failed. It counts the solute's
  residues now, and skips rather than errors.

- **Documentation that had gone stale without anything noticing.** The
  collective variables were listed as five in three places when there are
  eight, and the README said eight, so the two contradicted each other. The
  analysis count had outlived three additions. `pmf`, `metad_surface` and
  `steered_work` were filed under "Protein and ligand together", where a row
  lands if it is appended to the end of the file, implying they need a ligand.
  Each claim was true when written. They are now checked by tests, because a
  count is exactly the kind of thing that goes wrong quietly.

- **A page watching a run said every phase was "Not run".** It showed that
  run's energy, its temperature and its speed in ns/day above a table
  reporting nothing had happened. The table reads the manifest, and the
  manifest is written when a run finishes -- so "Not run", which is a claim
  that a phase did not happen, stood in for "has not finished". A phase in
  progress now names its stage, the ones before it are complete because a run
  that is simulating has finished preparing, and the ones after are pending.
  An empty directory still reads as not run, because that is a different
  thing.

- **The chart titles were drawn over the charts.** Absolutely positioned at
  the top-left of each canvas, each sat exactly where the chart draws its
  y-axis labels, and the two overlapped as soon as there was data to label.

- **The GUI re-read the whole trajectory on every mount, and said so.** Its
  terminal filled with `dcdplugin) detected standard 32-bit DCD file`, a pair
  of lines every few seconds for as long as it was open. The lines are VMD's
  DCD plugin, which MDTraj wraps, announcing every file it opens on the
  C-level file descriptor where Python's logging cannot reach it.

  The noise was the symptom. The viewer asked for the playback payload with
  `force=1` whenever the browser did not already hold one, so every mount and
  every reload bypassed the disk cache and re-streamed the trajectory to
  rebuild a file already sitting beside it. Force means the user asked for a
  rebuild. The reader is also silenced, so a legitimate read no longer writes
  to a server's terminal.

- **The Live Simulation page rendered an apparatus for reading data it did
  not have.** Opened on a run without telemetry it showed every field --
  current stage, current step, total planned steps, frames written, simulation
  time, elapsed time, checkpoint, last update -- reading "not available",
  twelve of them, above empty charts. It now says why, once, and shows nothing
  else.

  And it explained the absence by advising the reader to start the dashboard
  during a simulation, which is what somebody looking at that page has just
  done. Live telemetry is written only when a run asks for it and that is off
  by default, so the advice sent them to repeat what had already failed. The
  message names the setting and the flag that turn it on.

- **The reproducibility section implied a rerun would reproduce the study.**
  It gives an equivalent one. Solvation places water by a procedure that
  cannot be seeded -- OpenMM's `addSolvent` takes no seed -- so the same
  configuration run twice gave 37,251 atoms and then 37,763, with
  `random_seed` fixed both times: the seed fixes the dynamics, not the
  solvent. The section says so, and names what does repeat a study exactly,
  which is to simulate from the system it prepared.

- **The same observable was reported twice in one document with two different
  means.** A real study's RMSD read 0.08895 in the convergence table and
  0.09461 in its own results section, both labelled the mean, with nothing to
  say why. The table averaged the whole series; the per-analysis finding
  discarded the relaxation first, as the method requires. The table now takes
  the same settled part, so the two agree by construction, and it says how
  many frames were discarded. Drift is still measured over the whole series,
  since that is the question of whether the run was still relaxing.

- **The summary called a run settled that held six independent samples.**
  Settled and adequately sampled are different questions, and it reported only
  the first: a study whose RMSD held six independent samples and whose density
  held ten was described as having entirely settled. It now says both.

- **A study of one system printed one line and then nothing.** Each run's
  output goes to its own log so that several running at once do not interleave
  into an unreadable screen -- and that was applied to a single run too, which
  has nothing to interleave with. A config naming one system therefore said
  "Exploring 1 molecular system" and went quiet for as long as the run took.
  Output is redirected only where it would collide.

- **The reference conda recipe had fallen behind the feedstock.** It exists
  so a dependency added here reaches the package, and the traffic runs both
  ways: the feedstock learns what the conda-forge solver does and what its
  review asks for, and two such lessons had not come back. Its tests import
  the small-molecule stack, because declaring a dependency is not the same as
  it resolving and a broken solve should fail the build rather than somebody's
  first protein-ligand run. And its pip check is disabled, because
  openff-toolkit pulls in AmberTools components declaring numpy<2 against a
  numpy 2.x environment -- a reason recorded there and nowhere else, so a
  regeneration from this copy would have reintroduced a build failure already
  diagnosed. Both are back, with a test for each, and the recipe says
  graphical interface where it said browser.

- **A surface-area calculation that did not finish writing is computed
  again.** On Windows, MDTraj's `shrake_rupley` returns a final frame that was
  not fully written -- a partial value followed by zeros -- on the first call,
  and the same frame complete on a second call to the same trajectory. Frames
  other than the last come back bit-identical.

  The truncated answer identifies itself by how much of a frame reads exactly
  zero: the fault leaves a row like `[0.5354027, 0, 0, 0, 0]`, four fifths of
  it zero against none in the frames around it. Compared against the run's own
  median, because a real protein has buried residues whose surface area is
  exactly zero and which flicker between zero and a little above it as the
  structure breathes. Where a frame stands out the areas are computed again,
  up to five times, and the complete result used with the run recording that
  it had to be.

  Three things about the defect were assumed and all three were wrong. It is
  not confined to the final frame -- frames 0, 2 and 5 have been seen. A
  second call is not reliably clean. And the first attempt to recognise it
  asked whether a residue was exposed in some frames and exactly zero in
  others, which is true of every buried residue in every protein: it refused a
  solvated T4 lysozyme outright, five attempts and eighty seconds before
  giving up, on a run that was perfectly good. A test now measures the rate
  over twenty calls rather than assuming a shape, and its failure message
  carries what a fix would need.

  This is a workaround, not a rescue, and the distinction is the one this
  package draws elsewhere: a failed simulation is diagnosed and stopped
  because a salvaged trajectory looks like one that never needed salvaging.
  Here the wrong answer identifies itself and a correct one may be one call
  away.

  Without it, an unwritten frame pulls a six-frame mean down by a sixth, and
  the number reaching a report, with an error bar on it, looks like a
  measurement. Reported upstream; it affects anyone computing
  solvent-accessible surface area on Windows.

- **A SASA test was asserting a property of MDTraj.** Two routes to a
  per-residue average each called the surface-area calculation once and
  compared the results, so the test was also asserting that two identical
  calls agree. On Windows they did not: by eight per cent on one residue in
  one run, then by fifteen per cent on all five in the next, with the
  per-frame values matching every other platform and the direct route alone
  coming out low. Different Python versions failed each time -- 3.11, then 3.9
  and 3.12 -- which is nondeterminism rather than a version-specific bug.

  The two routes are now given the same areas, so the test asks about this
  package's arithmetic and nothing else, and a separate test calls the
  surface-area calculation ten times on one trajectory and requires the
  answers to be identical. If that one fails, the fault is upstream and the
  test is the reproduction to send there.

- **A mean over frames was accumulated in single precision.** `numpy.mean` on
  a float32 array sums in float32, so the average exposure of a residue lost
  digits over a long run that the per-frame route -- grouped by pandas in
  double -- did not. The two answers then differed in their last figures for
  no reason a reader could guess. Both are in double now.

- **A results section carried a figure and no number.** The mean, its
  uncertainty and the independent samples behind it are recorded for every
  per-frame analysis and were shown nowhere; each analysis now states what it
  measured, with the reason attached where the run cannot support it.

- **An analysis was described as running with default options beside a list of
  the options it ran with.** A `for ... else` runs its else when the loop
  finishes without a break, which is every time.

- **The shareable archive omitted the record of what produced it.** The bundle
  carried every output including the trajectory, and neither `manifest.json`
  nor `resolved_config.yml`: not by exclusion, but because both are written
  once every phase has finished and the bundle is built during the report
  phase, so they did not exist yet. A recipient got thirteen megabytes of
  results with no way to trace them and no file to rerun them from -- while
  the module's docstring said it contained the manifest. They are added once
  they exist.

- **The dashboard header showed the input path** where the report and the
  slides show the system's name. Its output-folder field and its
  `fastmdx gui --output ...` instruction keep their paths: that page is opened
  on the machine that produced the run, and both need one to be of any use.

- **The summary said nothing about the study.** The first section a reader
  reads was "This report was generated automatically by FastMDXplora from the
  outputs of an end-to-end molecular dynamics study" -- a statement about the
  software, in a document about their system. It now says what was simulated,
  for how long, at what temperature, and how much of the run can be
  interpreted, so the caveat arrives before the results rather than six pages
  after them.

- **The slides were a slideshow.** Twenty-one slides carrying twelve figures
  and not one number; three of them displaying filesystem paths from the
  machine that produced them; and the whole deck in 4:3, which shows on a
  modern projector with a black band down each side. The deck is 16:9, the
  path slides state how the system was built and how it was simulated, and a
  final slide gives every observable's mean with its uncertainty and the
  independent samples behind it -- from the same assessment the report uses,
  so the two cannot disagree.

- **The slide outline described a deck that no longer existed.** Three of its
  five sections read "See `setup_parameters.json`". It is built from the same
  content as the slides.

- **The methods section stated a ligand that was never there.** Run on a
  tri-alanine in water it read "The ligand LIG was parameterized with
  openff-2.2.1". Both values behind that sentence are defaults: LIG is the
  naming convention for a ligand if there is one, and the small-molecule force
  field is a property of the protein force field, present whether or not it
  was used. The sentence therefore appeared in every run. It is now written
  only where setup recorded that it prepared a ligand, which is the evidence
  rather than the setting. A methods section is the part of a report that gets
  published.

- **And stated the wrong provenance for the coordinates.** Whether they came
  from the Protein Data Bank was decided by the input being four characters
  long, so a local file named `abcd` was described as a deposited structure.
  Setup records the input form; it is read.

- **The report title carried the author's home directory.** "FastMDXplora
  Study --- /home/claude/ala3.pdb", on the first line of a document meant to
  be sent to somebody. The title names the system: `ala3`, or `181L` where the
  input was a PDB entry.

- **A requested output that could not be produced left no record.** The
  terminal warned that WeasyPrint was absent and no PDF would be written, and
  the manifest did not: read from its files afterwards, the run showed four
  formats where five were asked for, with nothing to say the fifth had been
  attempted. `not_produced.json` records what was asked for, and why it is not
  there.

- **A message advertised an input the software does not accept.** Asked to
  classify something it did not recognise, setup replied that it expected "a
  PDB file path, a 4-character PDB ID, or a one-letter amino-acid sequence" --
  and a sequence, passed, is recognised and then refused two steps later,
  because building a structure from one needs a predictor this software does
  not carry. The message now says what works, says plainly that a sequence
  does not, and the refusal names what to do instead. Found while checking a
  manuscript's claims against the code, where it had already propagated.

- **One implementation of the correlation statistic, and one verdict from
  it.** The report layer already measured independent samples from the
  autocorrelation time; a second implementation was written for the per-frame
  analyses before that was noticed. They agreed to two decimal places on every
  correlated series and disagreed on a constant one, where the report layer
  was right -- a series that never changes is one observation however many
  frames of it there are. There is now one function.

  They also asked differently whether a series can measure its own correlation
  time. The report's rule was a tenth of the run, which is the usual working
  limit and lets the flattering case through: a series with a true correlation
  of 2000 measured 361 over 4000 frames, and 361 is under a tenth of 4000. Both
  layers now halve the series and ask whether the estimate moves, so one run
  cannot be measurable in the report and unresolved in the findings.

  It lives at the top level, beside the other cross-cutting modules, rather
  than inside the analysis package: nothing registers it, it produces no
  figure, and the report asks it the same question the analyses do.

## [2.4.0] — 2026-08-05

This release is about making a simulation go where an ordinary one will not.

Most of what is interesting in a molecular system happens too rarely to see.
A ligand unbinds once a second and a simulation runs for a microsecond, so
plain molecular dynamics watches the bound state fluctuate and learns nothing
about leaving it. Enhanced sampling pays for the rare event with a bias, and
the whole difficulty is knowing what the biased run entitles you to say.

Three methods arrive together, and each is explicit about what its output is
and is not: metadynamics gives a free-energy surface if the bias converged,
steered MD gives a pathway and the work along it and *not* a free energy, and
umbrella sampling gives a potential of mean force if the windows overlap. That
last one refuses more often than it produces -- which is the point. A seven-
window study of benzene leaving the T4 lysozyme cavity was refused here
because two adjacent windows shared no sampling at all, and a tool that
returned a curve for it would have drawn a line through a region nothing
visited.

### Added

- **Umbrella sampling, from a config block to a free energy.** One block
  expands into a run per window, which is the shape the batch machinery
  already runs, so scheduling, parallelism and per-GPU pinning come from the
  code that already does them. Each window writes its own PLUMED, the COLVARs
  are read back, the approach to each window is discarded, and the sampling is
  recombined. Verified end to end from files on disk: a barrier of
  10.1 kJ/mol against a known 10.0.

  Four refusals -- windows that do not overlap, a missing window, several
  systems at once, and two ways of biasing the same coordinate. A gap is told
  apart from its commonest cause: windows that never reached their centres,
  which want the opposite remedy from windows that simply do not touch. And
  nothing is concluded from histograms a run did not fill -- below
  `minimum_samples` the overlaps are reported and no free energy is offered,
  because an overlap from tens of points is noise given a decimal place.

  The block's own settings are checked, which nothing had done: a one-letter
  misspelling of `minimum_overlap` was accepted, ignored, and the study
  stitched at the three per cent default while its author believed it had
  demanded fifteen. Every refusal above could be switched off by a typo.

  How much overlap is enough is a judgement about how much evidence a joint
  needs, so `minimum_overlap` belongs to whoever is making the claim; three
  per cent is enough to stitch and it is thin. The refusal states the
  threshold it applied, so a genuine gap can be told from a strict setting.

- **One prepared system for a set of umbrella windows.** Seven windows of one
  study came out with 37,212, 37,254, 37,436 and 37,445 atoms: four different
  systems for one measurement. Solvation does not place water the same way
  twice, so preparing each window separately makes part of the difference
  between windows the solvent rather than the restraint. The system is
  prepared once into `shared_setup/` and every window simulates from it.
  `simulation.prepared_from` is an ordinary setting, useful on its own for a
  system prepared elsewhere.

- **Steered MD.** A spring attached to a collective variable, with the anchor
  moved, dragging the system whether or not it wants to go. It gives a pathway
  and the work along it, not a free energy: the work depends on how fast the
  anchor moves, and a single fast pull spends most of it pushing water aside
  rather than on the interactions of interest, overestimating the barrier.
  Jarzynski recovers a free energy from an ensemble of pulls dominated by rare
  low-work trajectories, so it needs many repeats. The rate and the work are
  reported and no free energy is claimed.

- **Metadynamics from a named collective variable**, rather than PLUMED input
  written by hand. Eight variables, each stating what it does *not* separate
  -- the failure mode of the method is biasing something that does not
  distinguish the states that matter, after which the surface converges and
  describes a different system. Well-tempered by default, and the hill width
  is refused rather than guessed.

  Coordination counts contacts through a switching function rather than a
  step, so it does not break when a ligand rotates. Membrane depth is measured
  against the bilayer's own centre, because a membrane drifts and depth
  against a fixed plane becomes depth against nothing.

- **Walls and funnels.** An unbounded ligand run is refused: biasing a
  ligand's distance pushes it into bulk solvent, where the landscape is flat
  and the bias fills a basin that is effectively infinite. Not wrong so much
  as unfinishable.

- **Restraints, released in stages.** A structure that has just been minimised
  is not at equilibrium, and heating it lets the solute move as well as the
  solvent. Position, distance, angle and torsion restraints hold it while the
  solvent settles, stepped down through equilibration and off before
  production -- a biased production run measures the bias. On a solvated
  peptide, heavy atoms move about a sixth as far restrained as free. The
  methods section reports what was held and how it was let go.

- **Membrane systems.** A protein embedded in a POPC, POPE, DLPC, DLPE, DMPC,
  DOPC or DPPC bilayer, packed by OpenMM so no external tool is needed. The
  orientation is checked rather than assumed: `addMembrane` puts the bilayer
  in the xy plane and a PDB entry is usually in a different frame, so a
  structure can be rotated onto its longest axis, and two refusals guard that
  -- whether there is a longest axis worth rotating onto, and whether the
  result has the hydrophobic belt a bilayer-spanning protein has.

  The barostat is chosen from the topology. An ordinary one scales x, y and z
  together, squeezing a bilayer that should change thickness independently of
  its area, and area per lipid is what membrane simulations are validated
  against: the run completes and is wrong.

- **Water sites.** Most water in a simulation is bulk, but some positions are
  held throughout -- a water wedged between a ligand and a backbone carbonyl,
  bridging a hydrogen bond neither could make alone -- and displacing one
  costs entropy or gains affinity. Clustering positions and reporting
  occupancy is standard; the distinction is not. A cluster occupied in every
  frame is either one molecule that stayed, which has a residence time and is
  a molecule to displace, or a position many waters passed through, which is
  geometry the protein favours. Both are reported, and only waters near the
  solute count: bulk clusters beautifully and means nothing.

- **A diagnosis when a run fails**, read from the state it failed in. Which
  atoms went non-finite tells a wrong ligand parameter from a bilayer packing
  problem from an integration failure, and those need different remedies --
  the message used to advise a smaller timestep for all of them, which cannot
  fix a wrong parameter. Nothing is retried: a rescue that produces a
  trajectory from a broken system is worse than a failure, because the failure
  is visible.

- **An explanation of each step, while it happens.** Molecular dynamics has a
  lot of steps that are obvious once you know them and opaque before that, and
  a pipeline that does all of it silently is faster to use and teaches
  nothing. Fifteen explanations, each saying *why* rather than repeating what
  the step already said, with a citation where there is one worth following.
  On by default; `--no-explain` turns them off. A reference must carry authors
  and a year or be absent, because inventing one to look thorough would be
  worse than having none.

- **Run options in the browser.** The form was built from phase settings only,
  so a top-level setting reached the command line and the config file and not
  the browser -- which was the one interface that could not turn explanations
  on or off. A guard checks that every top-level setting either reaches the
  form or is on a list of exclusions with a reason.

### Changed

- **The settings are grouped by what they decide.** Thirty-six for setup and
  thirty-seven for simulation arrived as one flat list each, in the order they
  happened to be declared: a pH sat beside a dispersion correction, and
  finding the one you wanted meant reading all of them. Each phase now opens
  into named groups that say what they are about, and `--help` reads in the
  same sections -- one grouping declared in the schema and seen twice, rather
  than two that agree until one is edited. A setting added without being
  placed fails a test.

- **A settings block can be written in the browser.** `umbrella`, `steered`
  and `metadynamics` are blocks of several settings, not one value each. They
  reached the form -- the schema declares them -- as single-line text boxes,
  and what was typed arrived as a string no phase could read, so the browser
  was the one interface where enhanced sampling could not be set up at all.
  Each now gets a box you write the block into, one setting per line, with an
  example of the right shape showing until you type.

- **The documentation was reordered around how somebody meets it**, the
  reference pages point at what generates them, and the README leads with what
  can be studied rather than with what the software is not. Two pages
  documented a command that does not exist and an API that was never written;
  both are gone.

### Removed

- **The `gentle` preset.** It dropped the temperature to 100 K along with
  shortening the run, which is not a smoke test but different physics: water
  is ice at 100 K, and a run that survives there says nothing about whether
  the system is stable at 300 K. The smoke campaign script defaulted to it, so
  every campaign was simulating at 100 K -- the case where it mattered most,
  since a campaign exists to find out whether real structures prepare and
  simulate.

### Fixed

- **The banner described a different run from the one about to happen.** It
  reconstructed the settings from `sys.argv`, so a study driven by a config
  file printed the defaults for the step counts, the timestep and the
  temperature while using the config's values -- showing a million production
  steps while running five thousand. A banner is read at a glance and
  believed, which makes that worse than showing nothing.

- **A per-system block replaced the top-level one rather than merging**, so an
  expansion giving each window a block of only its umbrella settings discarded
  everything else the study asked for: a run requesting five thousand
  production steps ran a million, and nothing said so.

- **The umbrella pieces were committed unreachable.** The expansion was called
  by nothing, and after that the recombination was called by nothing -- a
  config with an umbrella block would have been accepted, ignored, and run
  once. The loader expands it, so every route in gets it, and the batch
  explorer recombines once the windows have finished.

- **PLUMED printed its banner three times, interleaved.** Silencing ours did
  not reach it: it writes from C++ straight to the file descriptor, where
  redirecting `sys.stdout` does not go. The descriptor moves to the run's own
  log, so nothing is lost and a worker's output sits beside its results.

- **A hydration shell clustered into a water site.** On ubiquitin the first
  shell chained into one cluster -- forty-eight thousand positions and eight
  hundred and fifty-four distinct waters, reported as a site occupied in every
  frame and described as mostly one molecule. A site now has to be compact,
  and "mostly one molecule" requires that molecule to hold most of the
  observations.

- **A run too short to show residence said it had.** Water on a protein
  surface exchanges on ten to a hundred picoseconds, so a site fully occupied
  through ten of them shows that no water left, not that one is held -- and
  every site in such a run looks the same. Below a nanosecond the finding says
  what the run cannot distinguish.

- **A water analysis failed on a system with no water.** Refusing is right in
  itself and wrong as a phase failure: an implicit-solvent run, or a
  trajectory stripped of solvent to save space, has no water sites and that is
  not an error.

- **An automatic choice was recorded as the user's.** The record said the site
  selection was "given" when the setting was `auto` -- a truthiness check
  reading a value as an absence, so `options.json` claimed a decision was
  somebody's that was not.

- **Two more counts had drifted.** The README said the trajectory is analysed
  "fifteen ways" and the `include` help said "all ten"; there are sixteen
  analyses registered. The help now says which run rather than how many -- ten
  always, water sites where there is water, five more where there is a ligand
  -- and the opening paragraph does not count at all.

- **The metadynamics help named five collective variables where there are
  eight.** It is what the browser shows beside the box, what `--help` prints,
  and what the generated config template carries, so one stale sentence was
  stale in four places. There is a test that fails when a variable is added
  and the sentence is not.

- **The Windows job failed on a test, not on the code.** `ctypes.CDLL(None)`
  asks for the running program's own symbols, which Windows does not have.
  Writing the check portably exposed a real defect: the worker restored the
  file descriptors without flushing the C runtime's buffer first, so whatever
  PLUMED had not flushed arrived in the shared terminal after the redirect was
  undone -- the leak the guard exists to prevent, arriving late.

- **`pip install "fastmdxplora[ligand]"` could not succeed.** The extra named
  `openff-toolkit`, which has no PyPI distribution at all, so pip failed to
  resolve the whole command rather than installing what it could reach. The
  toolkit is out of the extra and named where it can be had: the conda recipe,
  and the error raised at the point of use. A command that cannot work is
  worse than one that installs part of the answer and says what is left.

- **The `plumed` extra had the same defect, and worse.** `openmm-plumed` has
  no PyPI distribution either, and it was absent from the conda recipe -- so
  metadynamics would have failed at the point of use on the primary channel.
  The existing guard missed it because it sat in a list that exempts a package
  from being checked at all.

- **The reference conda recipe carried the previous version too.** It is the
  copy a dependency is added to before the feedstock, so somebody reading it
  to see what this release requires was reading a file that said it was a
  different one. It gets the same check the alias has, against the same source
  of truth.

- **The `fastmdx` alias carried the previous version.** Its version is written
  in a file where the main package takes its own from the git tag, and a
  hand-written version beside a derived one drifts. The release workflow
  refused the tag, which is a check that fires too late to be comfortable, so
  the ordinary test run now checks it against the changelog -- the thing a
  release is cut from.

- **`fastmdx info` called a phase available that could not run.** It reported
  setup and simulation as "available" three lines above OpenMM and PDBFixer as
  "missing" -- "available" meant only that the module imported and had a
  callable `run`, which is true of a phase with nothing to run it on. One
  screen contradicted itself, and the block somebody reads first was the one
  that was wrong. A phase now says what it needs, or that it is ready, or what
  it will not be able to do: setup prepares a protein without the ligand
  stack, and the report is written without a PDF, so neither is reported
  unavailable for wanting them.

- **`fastmdx info` reported two backends where the software reaches for
  eight.** A PyPI install showed OpenMM and PDBFixer present and said nothing
  about the OpenFF toolkit, which a protein-ligand setup needs and which pip
  cannot install -- so somebody reading it would conclude their install was
  complete and find out otherwise three phases later. Every backend is listed
  now, grouped by what it is for, with a command that works for each. A
  backend that is installed and will not load is told apart from one that is
  absent, because they need different remedies -- and because importing
  WeasyPrint without Pango raises from the dynamic loader rather than as an
  ImportError, which crashed the command outright.


## [2.3.0] — 2026-08-05

This release is about knowing how far to trust a number.

Ten analyses were checked against the method each cites -- a paper, a
reference implementation, or the contract of the library beneath -- and eight
defects were found. Four of them change published numbers. The pattern behind
most of them was the same: software answering where the honest response is
that the question cannot be answered from what was given, and this release
teaches it to say so instead.

Protein-ligand interactions are typed rather than counted, with the chemistry
resolved before any geometry is measured and the route recorded. An occupancy
now carries the observation behind it. The report writes a methods paragraph
against the checklists journals apply, and a convergence assessment that says
what a run cannot support. The browser interface was rebuilt around a single
page that offers every setting the schema declares, and the command line now
does the same.


### Added
- **The report as a PDF**, alongside the Markdown. A Markdown file renders
  differently in every viewer and cannot be printed with the figures where the
  text put them. Needs the `pdf` extra on PyPI, where WeasyPrint's system
  libraries have to be found separately; the conda-forge package brings its
  own. Where they are absent the run says so and writes the other formats.

- **A methods section written against the published checklists.** A methods
  paragraph for a molecular dynamics study has to state a particular list of
  things, and that list is published -- in the Journal of Chemical Information
  and Modeling's reporting guidelines and in Communications Biology's
  reproducibility checklist. Every value was already recorded; what was
  missing was the assembly. Steps are given as time, a missing random seed is
  reported as missing, and anything the run did not record is named rather
  than filled in with what is usual.

- **A convergence assessment**, which says how much independent information a
  trajectory holds. Consecutive frames are nearly the same structure, so an
  error bar computed over them is too small by the square root of the
  correlation time -- measured at four and a half times on real data. Where a
  run is too short to measure its own correlation, the section says so rather
  than reporting the independence it merely failed to rule out.

- **A flag for every setting the schema declares.** The command line's option
  table was maintained by hand: twenty-four settings had no flag, and
  fifty-three of the sixty that did carried help text that had fallen behind
  the declaration. Flags, help, types and accepted values are now derived, so
  a flag cannot offer a value the software refuses.

### Changed
- **`contacts` is `pl_contacts`**, beside `pl_hbonds` and `pl_interactions`,
  because it measures protein-ligand contacts. The old name is gone rather
  than aliased: nothing has been released under it, and the alias made the
  orchestrator run the analysis twice into the same directory.

- **What an analysis works out is recorded apart from what it was told.** Both
  reach `options.json`; only the settings reach the report, which had been
  printing two pages of raw atom-index tuples. The findings appear there as a
  sentence each.

- **The banner describes the run rather than the software.** No frame, no
  heading, no list of output formats, no feature badges. What remains is the
  system, where it is going, and the settings for each phase that will run.

- **Settings that never applied are no longer offered.** Six analyses accepted
  a general atom selection and ignored it -- a protein-ligand measure works
  out both sides from the ligand's residue name. A measurement that looks
  restricted and is not is the same defect as a count that looks complete and
  is not.

### Fixed
- The solvated atom count and the pressure the barostat held were computed and
  discarded, so a methods section had to report as unrecorded two things the
  run had printed to the terminal.
- A template failure during solvation reached the user as OpenMM's raw
  message; the explanation was wired to system creation only, which solvation
  reaches first.
- HED, the oxidised dimer of mercaptoethanol, is listed as the crystallisation
  additive it is, beside BME which it comes from.
- The analysis manifest is written after `compute` as well as before it, so
  what an analysis learns about its own run is no longer thrown away.

- **What version 1 measured and this did not.** Every analysis was compared
  against its version 1 counterpart, option by option. Most differences were
  the same setting renamed -- ``atoms`` is ``selection``, ``beta_const`` is
  ``beta``, ``linkage_method`` is ``linkage`` -- but six were real.

  Omega dihedrals are measured again, with an ``angles`` option choosing which
  torsions to compute. Omega is the peptide bond itself, near 180 degrees in
  almost every residue, and the exceptions are the finding: a cis bond, most
  often before a proline.

  MDS is offered again beside PCA, t-SNE and UMAP. It answers a different
  question -- PCA finds the directions of largest variance in the coordinates,
  MDS an arrangement preserving the distances between frames, and that
  distance is RMSD.

  Hydrogen bonds take ``distance_cutoff``, ``angle_cutoff``, ``sidechain_only``
  and ``exclude_water``. The two cutoffs were not settable at all, and were
  written in two places, so a setting reaching only one would have left the
  per-frame count disagreeing with the bonds it counted.

  Clustering takes ``random_state`` and ``n_init``. The seed was fixed at 42
  inside a function, which made every run agree with every other and hid the
  question: k-means finds a local optimum, so a clustering that survives a
  change of seed is a finding and one that does not is an artefact of where
  the algorithm started.

  Not adopted: version 1 could shell out to an external ``mkdssp`` binary for
  secondary structure. MDTraj's implementation is used instead, because a
  system package is a poor dependency for an analysis that already works.

- **Protein-ligand interactions, typed.** Counting contacts says how much of
  the protein a ligand touches; this says what is holding it -- a salt bridge
  a charge change would destroy, a hydrophobic packing that tolerates one.
  Eight interaction types, each implemented against a published criterion
  named in its own docstring.

  The chemistry is resolved before the geometry is measured, and how it was
  resolved is recorded: the run's own setup phase, a supplied SDF, the
  Chemical Component Dictionary, or inference from the coordinates. Salt
  bridges and pi-cation interactions are refused where the ligand's charge was
  not determined, because they are claims about charge and a charge inferred
  from coordinates is ambiguous more often than not.

  Occupancy carries the observation behind it. A contact present in 450
  consecutive frames and one present in 450 alternating frames are both fifty
  per cent, and only the second has an error bar worth printing. Transitions
  between binding modes are counted always and given as probabilities only
  where enough were seen for a rate to mean anything.

  Where the published criteria disagree, the disagreement is recorded rather
  than quietly resolved. PLIP allows a hydrogen bond at 4.1 A and 100 degrees
  where the literature standard is 3.5 and 120, and counts fluorine as a
  halogen bond donor where the sigma-hole literature says organic C-F does
  not. Both of PLIP's choices remain reachable as settings.

**Read this before upgrading a study in progress.** Every analysis was checked
against the method it claims to implement -- a cited paper, a reference
implementation, or the contract of the library underneath. Eight of the ten
were wrong about something, and five of those change numbers a 2.2.0 run
produced: hydrogen-bond counts rise, the Q-value becomes the switching function
its cited paper defines, the radius of gyration becomes mass-weighted as the
documentation always said it was, contacts and hydrogen bonds measure across
the periodic boundary, and secondary-structure residue numbering was shifted
whenever a ligand was present.

Three more stop rather than report something meaningless: dimensionality
reduction on a structure that does not move, secondary structure where nothing
has a backbone, and protein-ligand hydrogen bonds where the ligand's
connectivity is missing. The `v1` analysis profile is gone.

### Changed
- **One page builds a run, and it does everything.** There were two -- one
  starting from a protein, one from a trajectory -- and they were the same
  thing: both wrote a config and ran it, differing only in which phases the
  config named. Built separately, one offered eleven of the eighty-three
  settings that exist and the other offered all of them. The page that
  replaces them asks what you have, what should happen to it, and what you
  want to change, and draws every control from the schema.

  It also takes a config you already have: checked for syntax and for settings
  that do not exist, then run exactly as it stands or opened into the form and
  changed. Opening never writes to it -- a change is saved as a new file,
  because the one on disk may be committed beside a paper.

- **A trajectory you already have is enough to start.** From GROMACS, from
  AMBER, from your own script. Point the browser at the folder, choose what to
  measure, and run it. The simulation does not have to be reproduced here
  first, which is what asking for it excluded.

- **Every field wanting a path can be browsed**, and only files of the kind
  being asked for are listed.

### Fixed
- **The banner describes the run it is printing.** Built from command-line
  flags, a run started with `--config` saw none of them: it announced a
  million production steps for a run that only analysed a trajectory, and the
  default pH for a config that said otherwise.
- **The timeline shows only the stages a run can reach.** An analysis of an
  existing trajectory reaches one of seven, and six greyed out forever reads
  exactly like a run that stalled.
- **The GUI opens on something to do.** `fastmdx gui` passed the working
  directory as an active run, so it opened on the overview of a run that did
  not exist.
- **Results land where you point them.** A results folder given as a path was
  flattened into a folder name, so `/Users/someone/work` became
  `Users_someone_work` inside the launch directory.
- **Clustering declares the methods it runs.** Its docstring named one where
  the code ran two, and its default was filled in from `None` afterwards, so
  nothing reading the signature could say what it would do.

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
- **Protein-ligand hydrogen bonds need the ligand's bonds.** A donor is found
  as a nitrogen or oxygen with a hydrogen bonded to it, so a ligand whose
  connectivity is missing from the topology can only be seen to accept, never
  to donate -- and the count reported one direction under a name promising
  both. The guard that supplied missing bonds fired only when the topology had
  none at all, which a PDB carrying the protein's connectivity but no CONECT
  records for its ligand does not. The analysis now says so instead of
  counting half.
- **Distances honour the periodic box.** Contacts and hydrogen bonds measured
  plain distances whatever the trajectory carried, while the Q-value measured
  with the box on the same frames -- two analyses answering "is this near
  that" differently in one run. A solvated trajectory is not always imaged,
  and a molecule split across the boundary looks far from everything it is
  actually touching: a bound ligand sitting across the wall reported **no
  contacts at all**. Both now use the unit cell when there is one, which
  changes nothing where there is not. `periodic=False` restores the previous
  measurement.
- **The radius of gyration is mass-weighted.** The docstring said it used
  masses from the topology; `mdtraj.compute_rg` weights every atom equally
  unless told otherwise, and when given masses still measures from the
  geometric centre rather than the centre of mass. So each hydrogen counted
  for as much as each carbon, a few per cent from what GROMACS's `gyrate`,
  cpptraj's `radgyr` and a published figure report. **Rg values will differ
  from 2.2.0**, in the direction the documentation already claimed. Pass
  `mass_weighted=False` for the unweighted quantity.
- **Dimensionality reduction says when there is nothing to decompose.** A
  structure that does not move has no variance, and PCA divides each
  component's variance by the total: the ratios came out NaN and the figure
  was labelled `PC 1 (nan%)` over a scatter of coincident points, which looks
  like a result and is not one. The only sign was a numpy warning about
  dividing by zero. All three methods now refuse, since there are no
  neighbourhoods to preserve among points that are all the same point.
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
