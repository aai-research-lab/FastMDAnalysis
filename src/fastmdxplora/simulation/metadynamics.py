"""Metadynamics without writing PLUMED input.

PLUMED can express almost any enhanced-sampling scheme, and the cost of that
is a language to learn before running the commonest one. Most metadynamics on
a protein-ligand system biases one of a handful of things: how far the ligand
has moved from its pose, how far apart two groups are, a torsion, or how
compact the protein is. Those do not need a language.

So this generates the input for them, and the PLUMED integration that already
exists runs it. Anything more elaborate is still written by hand and passed
through as before -- this is a shorter path to the common case, not a
replacement for the general one.

**What a collective variable commits you to.** Metadynamics fills the free
energy landscape along whatever you bias, and reports a free energy as a
function of it. If the variable does not distinguish the states that matter --
if two genuinely different arrangements have the same value -- the surface
converges and describes something that is not the system. This is the failure
mode of the method, it does not announce itself, and no amount of running
longer fixes it. Each variable below says what it separates and what it does
not.

**And a run that has not converged has no free energy.** The bias is still
growing, so the surface is still moving, and reading a barrier off it is
reading the current state of the filling rather than the landscape. What can
be checked -- whether the hills have stopped growing, whether the run has
revisited the states it left -- is reported rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

__all__ = [
    "COLLECTIVE_VARIABLES",
    "MetadynamicsPlan",
    "build_plumed_script",
    "plan_from_config",
]


#: What each variable separates, and what it does not. Written out because
#: choosing one is the decision the method turns on, and a name alone does not
#: carry it.
COLLECTIVE_VARIABLES: dict[str, str] = {
    "ligand_rmsd": (
        "how far the ligand has moved from its starting pose. Separates bound "
        "from unbound and one pose from another. Does not separate two "
        "unbound arrangements at the same distance, so the unbound basin is "
        "one broad well rather than the many states it really is."
    ),
    "ligand_distance": (
        "the distance from the ligand's centre to the binding site's. "
        "Separates depth of binding along one direction. Does not distinguish "
        "leaving by one route from leaving by another, which matters when the "
        "question is how a ligand escapes rather than how tightly it is held."
    ),
    "torsion": (
        "a single dihedral. Separates rotameric states cleanly. Does not "
        "separate anything coupled to that bond, so it suits a question about "
        "one torsion and misleads about a conformational change involving "
        "several."
    ),
    "radius_of_gyration": (
        "how compact the protein is. Separates folded from extended. Does not "
        "separate a correctly folded structure from a compact wrong one, "
        "which is the classic way a folding free energy surface comes out "
        "converged and wrong."
    ),
    "distance": (
        "the distance between two atom selections, by their centres. The "
        "general case of the ligand distance above. Does not separate two "
        "arrangements that happen to put the centres the same distance "
        "apart, which a pair of large groups often will."
    ),
}

#: Well-tempered metadynamics by default. Plain metadynamics keeps depositing
#: at full height forever, so the bias never settles and the surface never
#: converges; well-tempered shrinks the hills as a region fills, which is what
#: makes the free energy recoverable at all.
DEFAULT_BIAS_FACTOR = 10.0

#: How often to deposit, in steps. Too frequent and the hills correlate and
#: the error is underestimated; too rare and it takes forever.
DEFAULT_PACE_STEPS = 500

#: Hill height in kJ/mol. The convention is a value small compared with the
#: barriers being crossed.
DEFAULT_HEIGHT_KJMOL = 1.2


@dataclass(frozen=True)
class MetadynamicsPlan:
    """A metadynamics run, described in terms of what it biases."""

    collective_variable: str
    #: Atom selections the variable needs, already resolved to indices.
    atoms: dict[str, list[int]]
    sigma: float
    height_kjmol: float = DEFAULT_HEIGHT_KJMOL
    pace_steps: int = DEFAULT_PACE_STEPS
    bias_factor: float = DEFAULT_BIAS_FACTOR
    temperature_K: float = 300.0

    def as_record(self) -> dict[str, Any]:
        return {
            "collective_variable": self.collective_variable,
            "what_it_separates": COLLECTIVE_VARIABLES.get(
                self.collective_variable, ""),
            "sigma": self.sigma,
            "height_kjmol": self.height_kjmol,
            "pace_steps": self.pace_steps,
            "bias_factor": self.bias_factor,
            "well_tempered": self.bias_factor > 1.0,
            "n_atoms_biased": {k: len(v) for k, v in self.atoms.items()},
        }


def _plumed_list(indices: list[int]) -> str:
    """Atom indices as PLUMED writes them.

    PLUMED numbers atoms from one; everything else here numbers from zero.
    Getting this wrong biases the atom next door, which is a mistake that
    produces a plausible surface for the wrong coordinate.
    """
    return ",".join(str(int(i) + 1) for i in indices)


def plan_from_config(
    spec: dict[str, Any],
    topology: Any,
    *,
    temperature_K: float = 300.0,
    ligand_resname: str | None = None,
) -> MetadynamicsPlan:
    """Read a metadynamics block and resolve its selections to atoms."""
    import mdtraj as md

    variable = str(spec.get("collective_variable", "")).lower()
    if variable not in COLLECTIVE_VARIABLES:
        raise ValueError(
            f"Unknown collective variable {variable!r}. Available: "
            + ", ".join(sorted(COLLECTIVE_VARIABLES))
            + ". Anything else can be biased by writing PLUMED input directly "
            "and passing it as `plumed`."
        )

    mdtop = (topology if isinstance(topology, md.Topology)
             else md.Topology.from_openmm(topology))

    def select(expression: str, label: str) -> list[int]:
        found = [int(i) for i in mdtop.select(expression)]
        if not found:
            raise ValueError(
                f"The {label} selection {expression!r} matched no atoms, so "
                "there is nothing to bias."
            )
        return found

    atoms: dict[str, list[int]] = {}
    if variable in ("ligand_rmsd", "ligand_distance"):
        resname = spec.get("ligand_resname") or ligand_resname
        if not resname:
            raise ValueError(
                f"{variable} needs a ligand, and none was given or detected."
            )
        atoms["ligand"] = select(f"resname {resname}", "ligand")
        if variable == "ligand_distance":
            site = spec.get("site_selection")
            if not site:
                raise ValueError(
                    "ligand_distance measures from the ligand to a site, and "
                    "`site_selection` says where the site is -- the pocket "
                    "residues, usually. Without it there is no second point."
                )
            atoms["site"] = select(str(site), "site")
    elif variable == "distance":
        for name in ("selection_a", "selection_b"):
            expression = spec.get(name)
            if not expression:
                raise ValueError(f"`distance` needs `{name}`.")
            atoms[name] = select(str(expression), name)
    elif variable == "torsion":
        expression = spec.get("selection")
        if not expression:
            raise ValueError(
                "`torsion` needs a `selection` matching exactly four atoms, "
                "which are the dihedral."
            )
        found = select(str(expression), "torsion")
        if len(found) != 4:
            raise ValueError(
                f"A torsion is four atoms; {expression!r} matched {len(found)}."
            )
        atoms["torsion"] = found
    else:  # radius_of_gyration
        atoms["group"] = select(
            str(spec.get("selection", "protein and name CA")), "group")

    sigma = spec.get("sigma")
    if sigma is None:
        raise ValueError(
            "Metadynamics needs a `sigma`: the width of the hills, in the "
            "units of the variable being biased. It should be roughly the "
            "size of the fluctuations in that variable within a single state "
            "-- around 0.05 nm for a distance or an RMSD, around 0.35 rad for "
            "a torsion. There is no default that is right for an arbitrary "
            "coordinate, and a wrong one either smears the surface flat or "
            "takes forever to fill."
        )

    return MetadynamicsPlan(
        collective_variable=variable,
        atoms=atoms,
        sigma=float(sigma),
        height_kjmol=float(spec.get("height_kjmol", DEFAULT_HEIGHT_KJMOL)),
        pace_steps=int(spec.get("pace_steps", DEFAULT_PACE_STEPS)),
        bias_factor=float(spec.get("bias_factor", DEFAULT_BIAS_FACTOR)),
        temperature_K=float(temperature_K),
    )


def build_plumed_script(plan: MetadynamicsPlan, reference_pdb: str | None = None) -> str:
    """The PLUMED input for a plan.

    Written out rather than hidden, because it is the thing that decides what
    the run measures, and somebody checking a result should be able to read
    what was biased without reading this module.
    """
    lines = [
        "# Generated by FastMDXplora from a named collective variable.",
        f"# Biasing: {plan.collective_variable} -- "
        f"{COLLECTIVE_VARIABLES[plan.collective_variable]}",
        "",
    ]

    variable = plan.collective_variable
    if variable == "ligand_rmsd":
        if not reference_pdb:
            raise ValueError(
                "ligand_rmsd is measured against a reference structure, and "
                "none was given."
            )
        lines.append(f"cv: RMSD REFERENCE={reference_pdb} TYPE=OPTIMAL")
    elif variable == "ligand_distance":
        lines.append(f"lig: COM ATOMS={_plumed_list(plan.atoms['ligand'])}")
        lines.append(f"site: COM ATOMS={_plumed_list(plan.atoms['site'])}")
        lines.append("cv: DISTANCE ATOMS=lig,site")
    elif variable == "distance":
        lines.append(f"a: COM ATOMS={_plumed_list(plan.atoms['selection_a'])}")
        lines.append(f"b: COM ATOMS={_plumed_list(plan.atoms['selection_b'])}")
        lines.append("cv: DISTANCE ATOMS=a,b")
    elif variable == "torsion":
        lines.append(f"cv: TORSION ATOMS={_plumed_list(plan.atoms['torsion'])}")
    else:
        lines.append(f"cv: GYRATION ATOMS={_plumed_list(plan.atoms['group'])}")

    lines.append("")
    lines.append(
        "metad: METAD ARG=cv "
        f"SIGMA={plan.sigma:g} "
        f"HEIGHT={plan.height_kjmol:g} "
        f"PACE={plan.pace_steps:d} "
        f"BIASFACTOR={plan.bias_factor:g} "
        f"TEMP={plan.temperature_K:g} "
        "FILE=HILLS"
    )
    lines.append("")
    # The variable and the bias, every deposition. Without these there is no
    # way to tell afterwards whether the run converged, and a metadynamics run
    # whose convergence cannot be checked has not measured a free energy.
    lines.append(
        f"PRINT ARG=cv,metad.bias STRIDE={plan.pace_steps:d} FILE=COLVAR")
    return "\n".join(lines) + "\n"
