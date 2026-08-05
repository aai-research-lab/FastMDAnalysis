"""Why each step happens, said while it happens.

Molecular dynamics has a lot of steps that are obvious once you know them and
opaque before that. Why is the protein put in a box of water? Why is the box
bigger than the protein? Why heat it before running it? Why is there a
separate stage where the volume can change?

A pipeline that does all of this silently is faster to use and teaches
nothing. Somebody running their first simulation ends up with a trajectory
they cannot defend, because they cannot say why any of it was done.

So FastMDXplora says what it is doing and why, and cites something worth
reading where there is one. It is on by default and `--no-explain` turns it
off, because the fourth time through it is noise.

**The explanations are keyed, not matched.** A call site names the explanation
it wants, so an explanation cannot drift onto the wrong step and a step cannot
quietly lose its explanation -- both are checked. Free-text matching would
have made the first failure silent.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Explanation", "EXPLANATIONS", "explain"]


@dataclass(frozen=True)
class Explanation:
    """Why a step happens, and where to read more."""

    why: str
    #: A citation worth following. Not every step has one, and inventing a
    #: reference to look thorough would be worse than having none.
    reference: str | None = None

    def as_text(self, width: int = 74) -> str:
        import textwrap

        body = textwrap.fill(self.why, width=width,
                             initial_indent="      ", subsequent_indent="      ")
        if self.reference:
            body += f"\n      → {self.reference}"
        return body


#: Keyed by step, so a call site asks for the explanation it means.
EXPLANATIONS: dict[str, Explanation] = {
    "heterogens": Explanation(
        why=(
            "A crystal structure contains more than the protein. Some of it "
            "is biology -- a bound ligand, a structural metal -- and some is "
            "the chemistry that made the crystal grow, like glycerol or "
            "buffer molecules. Simulating the second kind wastes computation "
            "and can hold the protein in a shape the crystal imposed. Each "
            "one is classified rather than kept or dropped wholesale."
        ),
        reference="Davis et al., The crystallographic heterogen problem, Acta Cryst D 2008",
    ),
    "protonation": Explanation(
        why=(
            "X-rays do not see hydrogens, so a crystal structure has none, "
            "and a simulation needs every one. Which histidines are "
            "protonated, and whether a glutamate is charged, depends on the "
            "local environment and on pH -- and those choices change the "
            "hydrogen bonding that holds the structure together. They are "
            "decided here rather than left to a default."
        ),
        reference="Olsson et al., PROPKA3, J Chem Theory Comput 2011",
    ),
    "ligand_chemistry": Explanation(
        why=(
            "A structure file gives a ligand's atoms and where they are, but "
            "not its bond orders, its aromaticity or its charge -- and a "
            "force field needs all three. Getting a bond order wrong moves a "
            "hydrogen, and a moved hydrogen invents or destroys a hydrogen "
            "bond. The chemistry is looked up rather than guessed wherever "
            "it can be."
        ),
        reference="Westbrook et al., The Chemical Component Dictionary, Bioinformatics 2015",
    ),
    "ligand_parameters": Explanation(
        why=(
            "The protein force field knows the twenty amino acids and "
            "nothing else, so a bound ligand has no parameters until "
            "somebody makes them. OpenFF generates them from the ligand's "
            "chemistry, which is why the chemistry had to be settled first."
        ),
        reference="Qiu et al., OpenFF 2.0.0 Sage, J Chem Theory Comput 2021",
    ),
    "solvation": Explanation(
        why=(
            "Proteins fold the way they do because of water -- the "
            "hydrophobic effect is a property of the solvent, not of the "
            "protein. Simulated in vacuum, a structure collapses in on "
            "itself. The box is padded well beyond the protein so it never "
            "interacts with its own periodic image, and ions are added to "
            "physiological concentration because charges in real cells are "
            "screened."
        ),
        reference="Jorgensen et al., Comparison of simple potential functions for water, J Chem Phys 1983",
    ),
    "membrane": Explanation(
        why=(
            "A membrane protein simulated in water is not the protein: the "
            "hydrophobic belt that normally sits in the bilayer is exposed "
            "to solvent, and the helices splay apart. The lipids are packed "
            "around it so the protein sees the environment it evolved in."
        ),
        reference="Lomize et al., OPM database and PPM server, Nucleic Acids Res 2012",
    ),
    "minimize": Explanation(
        why=(
            "The starting structure has strain in it -- atoms slightly too "
            "close, bonds slightly too long -- from the experiment, from "
            "adding hydrogens, and from dropping the protein into water. At "
            "the temperature of a simulation that strain becomes violent "
            "motion. Minimisation walks the structure downhill to a nearby "
            "arrangement with no such forces in it, before anything moves."
        ),
    ),
    "restraints": Explanation(
        why=(
            "A minimised structure is not at equilibrium, and heating it "
            "lets the solute move as well as the solvent -- side chains "
            "relax into the space crystal packing left, and a ligand drifts "
            "out of the pose that was measured. Holding the solute while the "
            "water arranges itself around it, then letting go in stages, "
            "means production starts from the structure somebody determined."
        ),
        reference="Roe & Brooks, A protocol for preparing explicitly solvated systems, J Chem Phys 2020",
    ),
    "nvt": Explanation(
        why=(
            "The system starts at zero temperature with the atoms sitting "
            "still. NVT brings it to the temperature you asked for while the "
            "volume is held fixed, which lets the water find its arrangement "
            "without the box changing size at the same time. Two things "
            "equilibrating at once is harder to diagnose when it goes wrong."
        ),
    ),
    "npt": Explanation(
        why=(
            "Now the box is allowed to change size, so the system settles to "
            "the density real water has at this temperature and pressure. "
            "Skipping this leaves the system at whatever density the "
            "solvation happened to produce, which is usually a little wrong "
            "and shows up in everything measured afterwards."
        ),
    ),
    "membrane_barostat": Explanation(
        why=(
            "A bilayer must be free to change thickness independently of its "
            "area, so pressure is coupled in the membrane plane and along "
            "the normal separately. An ordinary barostat scales all three "
            "directions together, which squeezes the membrane and gives the "
            "wrong area per lipid -- the number membrane simulations are "
            "validated against."
        ),
        reference="Chow & Ferguson, Isothermal-isobaric molecular dynamics, Comput Phys Commun 1995",
    ),
    "production": Explanation(
        why=(
            "This is the part that is analysed. Everything before it was "
            "getting the system into a state worth measuring; from here the "
            "trajectory is a sample of how the system behaves at "
            "equilibrium, and the frames written now are the ones every "
            "later number comes from."
        ),
    ),
    "metadynamics": Explanation(
        why=(
            "Some things are too slow to see by waiting -- a ligand leaving "
            "a pocket might take milliseconds, and a simulation runs for "
            "microseconds. Metadynamics fills in the energy landscape along "
            "the coordinate you named as the run proceeds, pushing the "
            "system out of wells it has already visited so it explores "
            "instead of sitting still. The free energy is recovered from the "
            "bias that was added."
        ),
        reference="Barducci et al., Well-tempered metadynamics, Phys Rev Lett 2008",
    ),
    "convergence": Explanation(
        why=(
            "Frames are not independent observations. Consecutive frames are "
            "nearly the same structure, so a thousand frames may hold only a "
            "handful of independent samples -- and an error bar computed as "
            "though they were all independent is too small by a large "
            "factor. What is reported is how much the trajectory actually "
            "supports."
        ),
        reference="Flyvbjerg & Petersen, Error estimates on averages of correlated data, J Chem Phys 1989",
    ),
    "interactions": Explanation(
        why=(
            "Counting contacts says how much of the protein a ligand "
            "touches. Typing them says what holds it there -- a salt bridge "
            "a charge change would destroy, or a hydrophobic packing that "
            "tolerates one. Those suggest different next experiments, which "
            "is why each contact is classified rather than counted."
        ),
        reference="Adasme et al., PLIP 2021, Nucleic Acids Res 2021",
    ),
}


def explain(key: str) -> Explanation | None:
    """The explanation for a step, or ``None`` where there is not one.

    Returning ``None`` rather than raising, because a missing explanation
    should not stop a run -- but the guard in the test suite fails on one, so
    it does not go unnoticed either.
    """
    return EXPLANATIONS.get(key)
