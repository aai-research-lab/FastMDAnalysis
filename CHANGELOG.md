# Changelog

All notable changes to FastMDXplora are documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) ·
Versioning: [SemVer 2.0.0](https://semver.org/spec/v2.0.0.html)

## [Unreleased]

### Added
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

### Removed
- **The `datasets` package.** It promised a Trp-cage trajectory that was never
  bundled: `TrpCage.traj` returned the path to a file that had never existed,
  from v0.1.0 through v2.4.0, so reading it gave a plausible string and
  passing it anywhere gave a file-not-found from inside a trajectory reader.
  Nothing replaces it. FastMDXplora fetches a deposited structure and
  simulates it, so a reference trajectory carried in the distribution would be
  megabytes of wheel for something one command produces --- and a reference
  system is a PDB identifier, not a file.

### Changed
- **The interface is called the GUI.** Documentation and comments used "the
  browser" for both the interface and the thing it renders in, so the three
  interfaces read as "the command line, a config file, and the browser". A
  browser tab and the `--no-browser` flag are the literal thing and keep the
  word; everything else naming the interface says GUI.

### Fixed
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
