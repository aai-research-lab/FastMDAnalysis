"""Pulling a system along a coordinate, and what that does and does not give.

Some things do not happen on their own within reach of a simulation. A ligand
may take milliseconds to leave a pocket where a run lasts microseconds. Steered
molecular dynamics attaches a spring to a collective variable and moves the
spring's anchor, dragging the system along whether or not it wants to go.

**What this gives you is a pathway, not a free energy.** The work done pulling
depends on how fast you pull: drag a ligand out in a nanosecond and most of the
work goes into pushing water aside and straining the protein, not into
breaking the interactions you meant to measure. That dissipated work does not
cancel, and a single fast pull overestimates the barrier, sometimes by a lot.

Jarzynski's equality recovers a free energy from an ensemble of pulls, and the
average is dominated by the rare low-work trajectories -- so it needs many
repeats and converges badly when the pulling is fast. FastMDXplora does not
claim a free energy from a steered run, and reports the work so the claim can
be made deliberately by somebody who has done the repeats.

**What steered MD is genuinely good for is generating starting structures.**
Pull once, take frames along the way, and each becomes a window for umbrella
sampling -- which does give a free energy, from equilibrium sampling at each
position rather than from work done in a hurry. That is the standard use, and
it is why this exists before umbrella sampling does.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from fastmdxplora.simulation.metadynamics import (
    COLLECTIVE_VARIABLES,
    MetadynamicsPlan,
    _plumed_list,
    plan_from_config,
)

__all__ = ["SteeredPlan", "plan_steered", "build_steered_script"]

#: A spring stiff enough to drag the system without the anchor running away
#: from it, in kJ/mol/nm². Softer and the coordinate lags behind the anchor;
#: much stiffer and the pull becomes a shove.
DEFAULT_FORCE = 2000.0


@dataclass(frozen=True)
class SteeredPlan:
    """A pull: from where, to where, how hard, and over how long."""

    #: The coordinate, resolved to atoms, reusing the metadynamics layer.
    cv: MetadynamicsPlan
    to_value: float
    from_value: float | None = None
    force_constant: float = DEFAULT_FORCE
    #: Steps over which the anchor travels. Everything about how much this
    #: run means depends on this being large.
    steps: int = 500000

    def rate_per_ns(self, timestep_fs: float) -> float | None:
        """How fast the anchor moves, in units of the variable per ns.

        The number that decides whether the work means anything. Reported
        rather than checked against a threshold, because what counts as slow
        depends on the coordinate and the system -- but an order of magnitude
        is usually obvious.
        """
        if self.from_value is None or not self.steps or not timestep_fs:
            return None
        duration_ns = self.steps * timestep_fs / 1e6
        if duration_ns <= 0:
            return None
        return abs(self.to_value - self.from_value) / duration_ns

    def as_record(self) -> dict[str, Any]:
        return {
            "collective_variable": self.cv.collective_variable,
            "what_it_separates": COLLECTIVE_VARIABLES.get(
                self.cv.collective_variable, ""),
            "from": self.from_value,
            "to": self.to_value,
            "force_constant": self.force_constant,
            "steps": self.steps,
            "gives": (
                "a pathway and the work done along it, not a free energy: "
                "the work depends on how fast the anchor moved, and a single "
                "pull overestimates a barrier"
            ),
        }


def plan_steered(
    spec: dict[str, Any],
    topology: Any,
    *,
    temperature_K: float = 300.0,
    ligand_resname: str | None = None,
) -> SteeredPlan:
    """Read a steered-MD block, reusing the collective-variable machinery."""
    if "to" not in spec:
        raise ValueError(
            "Steered MD needs a `to`: the value of the collective variable "
            "to pull towards. Without it there is no destination and nothing "
            "to steer."
        )

    # The same five variables metadynamics offers, resolved the same way,
    # with the same refusals. A steered run does not need a hill width, so
    # sigma is supplied and ignored rather than demanded of the user.
    cv_spec = {k: v for k, v in spec.items()
               if k not in ("to", "from", "force_constant", "steps")}
    cv_spec.setdefault("sigma", 0.05)
    cv_spec.setdefault("unbounded", True)
    cv = plan_from_config(cv_spec, topology, temperature_K=temperature_K,
                          ligand_resname=ligand_resname)

    steps = int(spec.get("steps", 500000))
    if steps <= 0:
        raise ValueError("Steered MD needs a positive number of `steps`.")

    return SteeredPlan(
        cv=cv,
        to_value=float(spec["to"]),
        from_value=(None if spec.get("from") is None else float(spec["from"])),
        force_constant=float(spec.get("force_constant", DEFAULT_FORCE)),
        steps=steps,
    )


def build_steered_script(plan: SteeredPlan,
                         reference_pdb: str | None = None) -> str:
    """The PLUMED input for a pull.

    Written out rather than hidden, because the pulling rate is the thing
    that decides whether the result means anything and somebody checking a
    number should be able to read it.
    """
    variable = plan.cv.collective_variable
    lines = [
        "# Generated by FastMDXplora from a named collective variable.",
        f"# Pulling: {variable} -- {COLLECTIVE_VARIABLES[variable]}",
        "#",
        "# The work reported below depends on how fast the anchor moves. A",
        "# fast pull overestimates a barrier, because the work goes into",
        "# pushing solvent aside as well as into the interactions of",
        "# interest.",
        "",
    ]

    if variable == "ligand_rmsd":
        if not reference_pdb:
            raise ValueError(
                "ligand_rmsd is measured against a reference structure, and "
                "none was given."
            )
        lines.append(f"cv: RMSD REFERENCE={reference_pdb} TYPE=OPTIMAL")
    elif variable == "ligand_distance":
        lines.append(f"lig: COM ATOMS={_plumed_list(plan.cv.atoms['ligand'])}")
        lines.append(f"site: COM ATOMS={_plumed_list(plan.cv.atoms['site'])}")
        lines.append("cv: DISTANCE ATOMS=lig,site")
    elif variable == "distance":
        lines.append(f"a: COM ATOMS={_plumed_list(plan.cv.atoms['selection_a'])}")
        lines.append(f"b: COM ATOMS={_plumed_list(plan.cv.atoms['selection_b'])}")
        lines.append("cv: DISTANCE ATOMS=a,b")
    elif variable == "torsion":
        lines.append(f"cv: TORSION ATOMS={_plumed_list(plan.cv.atoms['torsion'])}")
    else:
        lines.append(f"cv: GYRATION ATOMS={_plumed_list(plan.cv.atoms['group'])}")

    lines.append("")
    # A moving restraint: the anchor starts where the system is (or where you
    # say) and arrives at the destination on the last step.
    start = ("" if plan.from_value is None
             else f"AT0={plan.from_value:g} STEP0=0 "
                  f"KAPPA0={plan.force_constant:g} ")
    lines.append(
        "pull: MOVINGRESTRAINT ARG=cv "
        + (start or f"AT0=0 STEP0=0 KAPPA0={plan.force_constant:g} ")
        + f"AT1={plan.to_value:g} STEP1={plan.steps:d} "
        f"KAPPA1={plan.force_constant:g}"
    )
    lines.append("")
    # The work is the point of recording anything: without it a steered run
    # has produced a trajectory and no number.
    lines.append(
        "PRINT ARG=cv,pull.work,pull.bias STRIDE=100 FILE=COLVAR")
    return "\n".join(lines) + "\n"
