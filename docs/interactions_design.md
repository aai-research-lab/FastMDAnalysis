# Protein-ligand interactions: what to build, and why not to install it

## The gap

`pl_contacts` counts how many residues are near the ligand. `pl_hbonds` counts
hydrogen bonds. Neither says what kind of contact, which is what a medicinal
chemist is actually asking: is this ligand held by a salt bridge that a charge
change would destroy, or by hydrophobic packing that tolerates one?

OpenMMDL answers that through PLIP or ProLIF: hydrophobic contacts, hydrogen
bonds, salt bridges, pi-stacking, pi-cation, halogen bonds, water bridges,
metal coordination. That is the difference between counting contacts and
describing binding.

## Why not depend on either

**PLIP** is built for single structures. It takes a PDB file, so a trajectory
means writing one per frame, and it re-protonates every frame with OpenBabel --
a routine its own documentation calls non-deterministic. Across a trajectory
that would make interactions flicker for reasons unrelated to the dynamics,
and any fingerprint built on top would be measuring the protonation routine.
Measured at roughly 30-90 ms per frame on systems up to 3000 atoms.

**ProLIF** is the better tool: built for trajectories, RDKit-based, 17
interaction types, thresholds on the literature values, and about 19 ms per
frame on the same data. Its interaction rules are SMARTS, which is readable and
checkable rather than hidden behind a perception layer.

It is not taken as a dependency because it requires MDAnalysis -- `import
prolif` loads 127 MDAnalysis modules -- and this project does not depend on
MDAnalysis.

Both are worth installing in a development environment as second opinions.
Where all three agree, a rule is right; where they disagree, that is worth
finding out. That is how the doubling bug in version 1 was found.

## What we know that they have to guess

PLIP and ProLIF accept arbitrary PDB files, so they must perceive chemistry
from coordinates: which nitrogen is a donor, which ring is aromatic, what is
charged at this pH. That perception is where their disagreements come from.

This software does not have to guess:

- **The ligand's chemistry is already resolved.** Setup writes
  `setup/ligands/<resname>.sdf`, containing bond orders, formal charges and
  aromaticity from the Chemical Component Dictionary, in the protonation state
  settled against the pocket at the pH being simulated. That file is the
  ligand's chemistry, decided once, by the phase whose job it was.
- **The protein is standard residues.** Twenty of them, with known donors,
  acceptors, charges and aromatic rings. No perception required, and a
  non-standard residue is something to say so about rather than guess at.
- **Hydrogens are present and placed by the force field**, so donor geometry
  is measured rather than inferred, and nothing needs re-protonating per frame.

So the chemistry step -- the part that makes PLIP slow and its protonation
non-deterministic -- is already done before analysis begins.

## When the trajectory came from somewhere else

Most trajectories will not have been produced here. Someone arriving with a
GROMACS or AMBER run has no `setup/ligands/*.sdf`, and a PDB topology carries
connectivity but not bond orders: saving and reloading a molecule keeps the
bonds and loses the fact that one of them was double. Aromaticity and formal
charge go with it.

That is the ordinary case, not the exception, so the chemistry is resolved by
trying in order and **recording which route succeeded**, because the confidence
differs and the analysis should say how it knew:

1. **An SDF the person supplies.** They know their ligand; nothing should stop
   them saying so.
2. **The run's own `setup/ligands/<resname>.sdf`**, when analysing a run this
   software produced. Already settled, already at the right pH.
3. **The Chemical Component Dictionary, by residue name.** A ligand from a
   crystal structure usually keeps its three-letter code, and `setup/ccd.py`
   already fetches definitions by that code. The atom names have to match the
   trajectory's, which is checkable rather than assumed.
4. **Perception from coordinates**, via RDKit's `DetermineBonds`. This is what
   PLIP and ProLIF do for every structure, so it is not a poor route -- but it
   is a guess, and a guess about bond order changes which nitrogen donates and
   which ring is aromatic. Used only when the first three fail, and said so in
   the manifest.

Only where all four fail does the analysis stop, and then it says what would
let it continue -- an SDF, or a residue name the CCD knows.

The protein needs none of this while its residues are standard. A non-standard
one is reported rather than guessed at: a modified residue whose chemistry
nobody has stated is exactly the case where a silent guess becomes a silent
wrong answer.

## The rules, with their sources

Each is implemented against a citation, and the citation goes in the docstring.
Where the literature and an existing tool disagree, the disagreement is
recorded rather than silently resolved.

| interaction | criterion | source |
|---|---|---|
| hydrogen bond | D···A < 3.5 A, D-H···A > 120 deg | Baker & Hubbard 1984; McDonald & Thornton 1994 |
| hydrophobic | C···C < 4.0 A, both carbons bonded only to C or H | PLIP 4.0 A; ProLIF 4.5 A |
| salt bridge | charged centres < 4.5 A | ProLIF 4.5 A; PLIP 5.5 A |
| pi-stacking | ring centres < 5.5 A, planes near 0 or 90 deg, offset < 2.0 A | PLIP |
| pi-cation | charge to ring centre < 6.0 A, offset < 2.0 A | PLIP |
| halogen bond | X···A < 3.5 A, C-X···A near 165 deg | PLIP; Auffinger 2004 |
| metal coordination | metal to donor < 3.0 A | PLIP |
| water bridge | one water hydrogen-bonded to both | PLIP |

Note the disagreements already visible: PLIP allows a hydrogen bond at 4.1 A
and 100 degrees, where the literature standard is 3.5 A and 120. PLIP says this
is deliberate, refined against low-resolution crystal structures. For frames
with explicit hydrogens that reasoning does not apply -- the hydrogen positions
are known -- so the literature criterion is used and PLIP's is available as a
setting.

## Order of work

1. Hydrogen bonds and hydrophobic contacts. Most of the signal, and the two
   whose definitions are best established.
2. Validate against PLIP and ProLIF on the same frames. Understand every
   disagreement before continuing.
3. Salt bridges, pi-stacking, pi-cation.
4. Halogen bonds, metal coordination, water bridges.
5. Occupancy per interaction, with the count behind it. A contact present in
   3 frames of 500 and one present in 450 are both "present"; only one means
   anything.
6. Binding modes: which combinations of interactions occur together, and how
   often. Transitions between them only where the sampling supports the claim.

Step 5 is where this stops being a reimplementation. Reporting an occupancy
without saying how many observations it rests on is the same class of mistake
as reporting a per-frame hydrogen-bond count that silently excluded transient
bonds.


## What checking against the other tools found

The rules are checked against PLIP, ProLIF and MDTraj's own Baker-Hubbard on
the same frames. Where they agree, the rule is being read the same way; where
they do not, the difference is understood before continuing.

**Hydrogen bonds agree exactly.** MDTraj's Baker-Hubbard finds the same five
atom pairs on real coordinates. ProLIF finds the same four residue partners,
and the same counts once two stated differences are accounted for: it requires
a donor angle of 130 degrees where the literature standard used here is 120,
and it reports one interaction per residue pair where this reports every atom
pair. At 130 degrees counted per residue, the two agree frame for frame.

Neither difference is a defect in either. They are the reason the thresholds
are settings and the reason the count is left unreduced: which contact
represents a residue is a choice, and making it silently would hide it.

**And all three tools share one failure.** A topology without bonds cannot
answer a question about hydrogen bonds -- a donor is an atom with a hydrogen
bonded to it, and without bonds there are none. Given such a topology:

- version 1 of this software counted every bond twice, because it called
  `create_standard_bonds()` unconditionally on a topology that already had
  them;
- ProLIF returned an empty result and said nothing;
- this software, in its first version, returned no donors and said nothing.

Only the last of those is ours to fix, and it is fixed: a selection whose
hydrogens are bonded to nothing now says how many, and what would settle it.
But the pattern is worth naming, because it is the characteristic defect of
this category of software and not of any one implementation. **An answer
returned where the honest response is that the question cannot be answered
from what was given.** It is silent, it looks like a result, and nothing fails.
