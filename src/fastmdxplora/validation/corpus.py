"""Guardrails measured, rather than asserted.

Every guardrail in this software has a test showing it fires on a case
chosen because it should fire. That is sensitivity, and sensitivity alone
is not evidence: a checker that refused every study would score perfectly
on it. The claim worth making is that the guardrails fire on defects and
stay quiet on ordinary work, and the second half needs a corpus of
ordinary work with the expected answer written down first.

So this runs two corpora and reports both rates. A case names what it
does, what should happen, and why that is the right answer; the harness
records what did happen and compares. Nothing here runs dynamics: the
guardrails being measured decide before or after a trajectory, not during
one, so the whole corpus completes in seconds and can run on every
release rather than once before a paper.

Three outcomes are distinguished, because collapsing them would hide the
distinction the software is built on. A case may **proceed** with a
number, be **refused** with a reason, or be **qualified** -- answered, but
with a statement attached saying what the answer does not support. A
qualification is not a refusal, and counting it as one would make the
software look more obstructive than it is; counting it as a clean pass
would make it look less careful.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

__all__ = [
    "Case",
    "Outcome",
    "run_case",
    "run_corpus",
    "DEFECTS",
    "CLEAN",
]


@dataclass(frozen=True)
class Case:
    """One study, and what the software should say about it."""

    name: str
    #: What the case exercises. Returns whatever the guardrail returns;
    #: raising is one of the ways a guardrail can refuse.
    run: Callable[[], Any]
    #: "refused", "qualified" or "proceeded".
    expect: str
    #: Why that is the right answer. Recorded so a disagreement can be
    #: read without reconstructing the intent from the code.
    because: str
    #: A fragment the refusal or qualification must contain, so a case
    #: cannot pass by refusing for an unrelated reason.
    mentioning: str = ""


@dataclass
class Outcome:
    case: str
    expected: str
    observed: str
    agreed: bool
    detail: str = ""
    because: str = ""


def _classify(value: Any) -> tuple[str, str]:
    """What the software did, from what it returned or raised."""
    # A qualification may sit at the top of a returned record or one level
    # inside it, because an analysis writes its findings under its own
    # name. Reading only the top level scored a guardrail that had fired as
    # a miss, which is the instrument failing rather than the software, and
    # is the difference this harness exists to keep straight.
    QUALIFIERS = ("not_a_measurement", "capped", "calibration")

    def _qualification(record: dict) -> "tuple[str, str] | None":
        for key in QUALIFIERS:
            if record.get(key):
                return key, str(record[key])
        return None

    if isinstance(value, dict):
        if value.get("refused"):
            return "refused", str(value["refused"])
        found = _qualification(value)
        if found:
            return "qualified", found[1]
        for name, entry in value.items():
            if not isinstance(entry, dict):
                continue
            if entry.get("refused"):
                return "refused", f"{name}: {entry['refused']}"
            found = _qualification(entry)
            if found:
                return "qualified", f"{name}: {found[1]}"
    return "proceeded", ""


def run_case(case: Case) -> Outcome:
    """Run one case and say whether the software agreed with the record."""
    try:
        value = case.run()
    except Exception as exc:  # noqa: BLE001 - a raise is a refusal here
        observed, detail = "refused", str(exc)
    else:
        observed, detail = _classify(value)

    agreed = observed == case.expect
    if agreed and case.mentioning:
        agreed = case.mentioning.lower() in detail.lower()
        if not agreed:
            detail = (f"{observed} but not for the stated reason: {detail}")
    return Outcome(
        case=case.name, expected=case.expect, observed=observed,
        agreed=agreed, detail=detail, because=case.because,
    )


def run_corpus(defects: "list[Case]", clean: "list[Case]") -> dict[str, Any]:
    """Both corpora, with the two rates that matter reported together.

    The detection rate alone is what every guardrail test in this
    repository already shows. The false-refusal rate is the half that
    makes it a measurement.
    """
    defect_outcomes = [run_case(c) for c in defects]
    clean_outcomes = [run_case(c) for c in clean]

    detected = sum(1 for o in defect_outcomes if o.agreed)
    quiet = sum(1 for o in clean_outcomes if o.agreed)
    false_refusals = [
        o for o in clean_outcomes if o.observed == "refused"
        and o.expected != "refused"
    ]

    return {
        "defects": {
            "total": len(defect_outcomes),
            "detected": detected,
            "detection_rate": (detected / len(defect_outcomes)
                               if defect_outcomes else None),
            "missed": [o.__dict__ for o in defect_outcomes if not o.agreed],
        },
        "clean": {
            "total": len(clean_outcomes),
            "as_expected": quiet,
            "false_refusal_rate": (len(false_refusals) / len(clean_outcomes)
                                   if clean_outcomes else None),
            "false_refusals": [o.__dict__ for o in false_refusals],
            "other_disagreements": [
                o.__dict__ for o in clean_outcomes
                if not o.agreed and o not in false_refusals],
        },
        "outcomes": [o.__dict__ for o in defect_outcomes + clean_outcomes],
    }


# ---------------------------------------------------------------------------
# The corpora. Each case is small and self-contained on purpose: a case
# needing a prepared system would tie this to a fixture, and the point is
# that it runs anywhere, quickly, on every release.
# ---------------------------------------------------------------------------

def _numpy():
    import numpy as np
    return np


def _truncated_pmf() -> Any:
    np = _numpy()
    from fastmdxplora.simulation.binding import KB_KJMOL, binding_free_energy

    kt = KB_KJMOL * 300.0
    radius = np.linspace(0.05, 0.7, 200)
    potential = -20.0 * np.exp(-((radius - 0.3) / 0.15) ** 2)
    return binding_free_energy(radius, potential - 2 * kt * np.log(radius),
                               temperature_K=300.0)


def _complete_pmf() -> Any:
    np = _numpy()
    from fastmdxplora.simulation.binding import KB_KJMOL, binding_free_energy

    kt = KB_KJMOL * 300.0
    radius = np.linspace(0.05, 3.0, 400)
    potential = np.where(radius <= 0.6, -20.0, 0.0)
    return binding_free_energy(radius, potential - 2 * kt * np.log(radius),
                               temperature_K=300.0, bound_cutoff_nm=0.6)


def _peptide(with_hydrogens: bool = True):
    import mdtraj as md
    np = _numpy()

    rng = np.random.RandomState(4)
    top = md.Topology()
    chain = top.add_chain()
    positions = []
    for i in range(6):
        res = top.add_residue("ALA", chain, resSeq=i + 1)
        top.add_atom("N", md.element.nitrogen, res)
        positions.append([i * 0.38, 0.0, 0.0])
        if with_hydrogens:
            top.add_atom("H", md.element.hydrogen, res)
            positions.append([i * 0.38, 0.10, 0.0])
        top.add_atom("CA", md.element.carbon, res)
        positions.append([i * 0.38 + 0.15, -0.05, 0.0])
        top.add_atom("C", md.element.carbon, res)
        positions.append([i * 0.38 + 0.25, 0.05, 0.0])
    xyz = np.tile(np.array(positions)[None], (300, 1, 1))
    xyz += rng.normal(scale=0.015, size=xyz.shape)
    return md.Trajectory(xyz.astype(np.float32), top)


def _order_parameters_without_hydrogens() -> Any:
    from fastmdxplora.analysis.order_parameters import OrderParameters
    return OrderParameters().compute(_peptide(with_hydrogens=False))


def _order_parameters_on_a_settled_run() -> Any:
    from fastmdxplora.analysis.order_parameters import OrderParameters
    analysis = OrderParameters()
    analysis.compute(_peptide())
    return analysis.findings


def _gas(box: float = 4.0, boxed: bool = True):
    import mdtraj as md
    np = _numpy()

    rng = np.random.default_rng(3)
    top = md.Topology()
    solute = top.add_chain()
    for i in range(20):
        res = top.add_residue("ALA", solute, resSeq=i + 1)
        top.add_atom("CA", md.element.carbon, res)
    water = top.add_chain()
    for i in range(200):
        res = top.add_residue("HOH", water, resSeq=i + 1)
        top.add_atom("O", md.element.oxygen, res)
    xyz = rng.random((15, 220, 3)) * box
    traj = md.Trajectory(xyz.astype("float32"), top)
    if boxed:
        np_ = _numpy()
        traj.unitcell_lengths = np_.tile([box, box, box], (15, 1))
        traj.unitcell_angles = np_.tile([90.0, 90.0, 90.0], (15, 1))
    return traj


def _rdf_without_a_box() -> Any:
    from fastmdxplora.analysis.rdf import RadialDistribution
    return RadialDistribution(
        selection_a="name CA", selection_b="name O").compute(
            _gas(boxed=False))


def _rdf_past_half_the_box() -> Any:
    from fastmdxplora.analysis.rdf import RadialDistribution
    analysis = RadialDistribution(
        selection_a="name CA", selection_b="name O", r_max=5.0)
    analysis.compute(_gas())
    return analysis.findings


def _rdf_within_the_box() -> Any:
    from fastmdxplora.analysis.rdf import RadialDistribution
    analysis = RadialDistribution(
        selection_a="name CA", selection_b="name O", r_max=1.0)
    analysis.compute(_gas())
    return analysis.findings


def _mutation_against_the_wrong_residue() -> Any:
    import tempfile
    from pathlib import Path

    from fastmdxplora.setup.pdbfix import _check_mutation_matches

    class _Residue:
        def __init__(self, number, name):
            self.id, self.name = number, name

        def __iter__(self):
            return iter(())

    class _Chain:
        id = "A"

        def residues(self):
            return [_Residue(2, "LEU")]

    class _Topology:
        def chains(self):
            return [_Chain()]

    _ = tempfile, Path
    return _check_mutation_matches(_Topology(), "A", "VAL", 2, "V2A")


def _mutation_that_matches() -> Any:
    from fastmdxplora.setup.pdbfix import _check_mutation_matches

    class _Residue:
        def __init__(self, number, name):
            self.id, self.name = number, name

    class _Chain:
        id = "A"

        def residues(self):
            return [_Residue(2, "LEU")]

    class _Topology:
        def chains(self):
            return [_Chain()]

    _check_mutation_matches(_Topology(), "A", "LEU", 2, "L2A")
    return {"checked": True}


def _q_on_a_chain_too_short() -> Any:
    import mdtraj as md
    np = _numpy()
    from fastmdxplora.analysis.qvalue import QValue

    top = md.Topology()
    chain = top.add_chain()
    for i in range(3):
        res = top.add_residue("ALA", chain, resSeq=i + 1)
        top.add_atom("CA", md.element.carbon, res)
    xyz = np.zeros((5, 3, 3), dtype="float32")
    xyz[:, :, 0] = np.arange(3) * 0.38
    return QValue(selection="all").compute(md.Trajectory(xyz, top))


def _surface_with_a_stuck_dimension() -> Any:
    import tempfile
    from pathlib import Path

    np = _numpy()
    from fastmdxplora.simulation.metad_surface import compute_surface_2d

    rng = np.random.RandomState(0)
    n = 800
    heights = 1.2 * np.exp(-np.linspace(0, 3.5, n))
    moving = (np.where(rng.rand(n) < 0.5, -1.0, 1.0)
              + rng.normal(scale=0.25, size=n))
    stuck = 0.5 + rng.normal(scale=0.05, size=n)
    path = Path(tempfile.mkdtemp()) / "HILLS"
    with path.open("w") as fh:
        fh.write("#! FIELDS time cv1 cv2 sigma_cv1 sigma_cv2 height biasf\n")
        for t, (a, b), h in zip(range(n), np.column_stack([moving, stuck]),
                                heights):
            fh.write(f"{t * 0.5:.3f} {a:.5f} {b:.5f} 0.2000 0.1000 "
                     f"{h:.6f} 10.0\n")
    return compute_surface_2d(
        path, np.column_stack([moving, stuck]), points=40,
        names=("phi", "dist"))


def _density_at_constant_volume() -> Any:
    import tempfile
    from pathlib import Path

    np = _numpy()
    from fastmdxplora.analysis.thermodynamics import Thermodynamics

    rng = np.random.RandomState(0)
    lines = ['#"Step","Time (ps)","Potential Energy (kJ/mole)",'
             '"Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)",'
             '"Temperature (K)","Box Volume (nm^3)","Density (g/mL)"']
    for i in range(300):
        pot = -12000 + rng.normal(scale=40)
        kin = 3000 + rng.normal(scale=30)
        lines.append(f"{i * 500},{i * 1.0},{pot:.4f},{kin:.4f},"
                     f"{pot + kin:.4f},{300 + rng.normal(scale=3):.4f},"
                     f"64.000000,0.997000")
    path = Path(tempfile.mkdtemp()) / "state_data.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    analysis = Thermodynamics(state_csv=str(path))
    analysis.compute(None)
    return analysis.findings["thermodynamics"]



def _hills(centres, heights, sigma=0.2, biasf=10.0):
    """A one-variable HILLS file, written where a test can read it."""
    import tempfile
    from pathlib import Path

    path = Path(tempfile.mkdtemp()) / "HILLS"
    with path.open("w") as fh:
        fh.write("#! FIELDS time cv sigma_cv height biasf\n")
        for t, (c, h) in enumerate(zip(centres, heights)):
            fh.write(f"{t * 0.5:.3f} {c:.5f} {sigma:.4f} {h:.6f} "
                     f"{biasf:.1f}\n")
    return path


def _metadynamics_without_a_recrossing():
    np = _numpy()
    from fastmdxplora.simulation.metad_surface import compute_surface

    n = 600
    rng = np.random.RandomState(0)
    # The bias fills one basin and the system never leaves it.
    centres = 1.0 + rng.normal(scale=0.15, size=n)
    return compute_surface(
        _hills(centres, 1.2 * np.exp(-np.linspace(0, 3.5, n))),
        centres, points=60)


def _metadynamics_that_crossed_and_settled():
    np = _numpy()
    from fastmdxplora.simulation.metad_surface import compute_surface

    n = 1200
    rng = np.random.RandomState(0)
    centres = (np.where(rng.rand(n) < 0.5, -1.0, 1.0)
               + rng.normal(scale=0.25, size=n))
    return compute_surface(
        _hills(centres, 1.2 * np.exp(-np.linspace(0, 3.5, n))),
        centres, points=60)


def _a_torsion_read_as_a_straight_line():
    """The same hills, told the coordinate does not wrap.

    A rotamer either side of the +/-180 boundary is then two states an
    infinite distance apart, and the barrier between them is invented.
    """
    np = _numpy()
    from fastmdxplora.simulation.metad_surface import compute_surface

    n = 1200
    rng = np.random.RandomState(1)
    centres = (np.where(rng.rand(n) < 0.5, -3.0, 3.0)
               + rng.normal(scale=0.2, size=n))
    straight = compute_surface(
        _hills(centres, 1.2 * np.exp(-np.linspace(0, 3.5, n))),
        centres, points=60, periodic=False)
    circular = compute_surface(
        _hills(centres, 1.2 * np.exp(-np.linspace(0, 3.5, n))),
        centres, points=60, periodic=True)
    straight_barrier = (straight.get("evidence") or {}).get("barrier_kjmol")
    circular_barrier = (circular.get("evidence") or {}).get("barrier_kjmol")
    if (straight_barrier is not None and circular_barrier is not None
            and straight_barrier > 1.5 * circular_barrier):
        return {"refused": (
            f"Read as a straight line the barrier is "
            f"{straight_barrier:.1f} kJ/mol; read as the circle it is, "
            f"{circular_barrier:.1f}. The difference is the boundary, not "
            "the physics.")}
    return {"barrier_kjmol": straight_barrier}


def _a_run_too_short_for_its_own_correlation():
    np = _numpy()
    from fastmdxplora.statistics import summarise

    # A slow drift: each half is internally consistent, the halves differ.
    series = np.linspace(0.0, 1.0, 400) + np.random.RandomState(0).normal(
        scale=0.01, size=400)
    settled, reason = summarise(series)
    return {"mean": None if settled is None else settled.mean,
            "not_a_measurement": reason}


def _a_run_long_against_its_correlation():
    np = _numpy()
    from fastmdxplora.statistics import summarise

    rng = np.random.RandomState(0)
    series = 0.3 + rng.normal(scale=0.01, size=4000)
    settled, reason = summarise(series)
    return {"mean": None if settled is None else settled.mean,
            "not_a_measurement": reason}


def _a_box_too_small_for_its_cutoff():
    from fastmdxplora.setup.prepare import (
        _solvate_with_room_for_the_cutoff)

    class _Q:
        def __init__(self, v): self._v = v
        def value_in_unit(self, _u): return self._v

    class _Unit:
        nanometer = 1.0

    class _Modeller:
        def __init__(self): self.attempts = []
        topology = property(lambda self: self)
        def deleteWater(self): pass
        def addSolvent(self, _ff, **kw):
            self.attempts.append(float(kw["padding"]))
            self._box = 0.2 * float(kw["padding"])
        def getPeriodicBoxVectors(self):
            b = getattr(self, "_box", 0.0)
            return [[_Q(b if i == j else 0.0) for j in range(3)]
                    for i in range(3)]

    modeller = _Modeller()
    _solvate_with_room_for_the_cutoff(
        modeller, object(), {"padding": 1.0}, nonbonded_cutoff_nm=5.0,
        padding_nm=1.0, nonbonded_method="PME", unit=_Unit,
        most_it_may_grow_nm=0.5)
    # It stopped rather than growing without limit, which is the guardrail.
    return {"capped": (
        f"stopped after {len(modeller.attempts)} attempt(s) rather than "
        "growing the box without limit")} if len(
            modeller.attempts) == 1 else {"attempts": modeller.attempts}


def _a_box_that_already_fits():
    from fastmdxplora.setup.prepare import (
        _solvate_with_room_for_the_cutoff)

    class _Q:
        def __init__(self, v): self._v = v
        def value_in_unit(self, _u): return self._v

    class _Unit:
        nanometer = 1.0

    class _Modeller:
        def __init__(self): self.attempts = []
        topology = property(lambda self: self)
        def deleteWater(self): pass
        def addSolvent(self, _ff, **kw):
            self.attempts.append(float(kw["padding"]))
            self._box = 4.0 * float(kw["padding"])
        def getPeriodicBoxVectors(self):
            b = getattr(self, "_box", 0.0)
            return [[_Q(b if i == j else 0.0) for j in range(3)]
                    for i in range(3)]

    modeller = _Modeller()
    _solvate_with_room_for_the_cutoff(
        modeller, object(), {"padding": 1.0}, nonbonded_cutoff_nm=0.9,
        padding_nm=1.0, nonbonded_method="PME", unit=_Unit)
    return {"attempts": modeller.attempts}


def _a_selection_that_would_save_nothing():
    from fastmdxplora.simulation.runner import resolve_save_selection

    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    for i in range(4):
        res = top.add_residue("HOH", chain, resSeq=i + 1)
        top.add_atom("O", md.element.oxygen, res)
    _kept, described = resolve_save_selection(top, "not water")
    return {"capped": described}


def _a_malformed_save_selection():
    from fastmdxplora.simulation.runner import resolve_save_selection

    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue("ALA", chain, resSeq=1)
    top.add_atom("CA", md.element.carbon, res)
    return resolve_save_selection(top, "protien")


#: Studies with one named thing wrong. The expected answer was written
#: before any of them was run.
DEFECTS: list[Case] = [
    Case("binding free energy from a run that never reached bulk",
         _truncated_pmf, "refused",
         "the reference the well depth is measured against is not a "
         "reference unless the ligand is free there",
         mentioning="-2kT ln r"),
    Case("order parameters on a structure without hydrogens",
         _order_parameters_without_hydrogens, "refused",
         "the measurement is of a bond vector and the bond is absent",
         mentioning="without hydrogens"),
    Case("radial distribution with no periodic box",
         _rdf_without_a_box, "refused",
         "there is no volume to take a bulk density from",
         mentioning="no unit cell"),
    Case("mutation named against the wrong residue",
         _mutation_against_the_wrong_residue, "refused",
         "applying it would simulate a protein nobody chose",
         mentioning="that position holds"),
    Case("fraction of native contacts on a chain too short to have a fold",
         _q_on_a_chain_too_short, "refused",
         "no residue pair is far enough apart in sequence to make a "
         "tertiary contact",
         mentioning="min_seq_separation"),
    Case("two-dimensional surface with one coordinate stuck",
         _surface_with_a_stuck_dimension, "refused",
         "the surface across a coordinate that did not move is the shape "
         "of the bias",
         mentioning="dist"),
    Case("radial distribution asked for beyond half the box",
         _rdf_past_half_the_box, "qualified",
         "the range is cut to where the shells are complete, and the cut "
         "is stated rather than silent",
         mentioning="half the smallest box"),
    Case("density from a constant-volume run",
         _density_at_constant_volume, "qualified",
         "the density is a constant the setup fixed, so a mean with an "
         "error on it would describe arithmetic",
         mentioning="constant the setup fixed"),
    Case("metadynamics stopped before any recrossing",
         _metadynamics_without_a_recrossing, "refused",
         "a barrier the system never crossed is the shape of the bias, "
         "not of the landscape",
         mentioning="cross"),
    Case("a run too short for its own correlation time",
         _a_run_too_short_for_its_own_correlation, "qualified",
         "the effective sample count is an upper bound, and the error bar "
         "printed from it is too tight",
         mentioning="not long"),
    Case("a box smaller than twice the cutoff",
         _a_box_too_small_for_its_cutoff, "qualified",
         "growing the padding without limit to reach an impossible cutoff "
         "would solvate forever, so it stops and says what it tried",
         mentioning="stopped"),
    Case("a save selection that would keep nothing",
         _a_selection_that_would_save_nothing, "qualified",
         "a box of pure water is a legitimate study, so `not water` "
         "matching none of it saves everything and names itself",
         mentioning="matched none"),
    Case("a save selection that will not parse",
         _a_malformed_save_selection, "refused",
         "a mistake in the study file, fixable in one edit, that would "
         "otherwise affect every frame",
         mentioning="not a selection"),
]

#: Ordinary studies, where nothing should fire. This is the half that makes
#: the detection rate above a measurement rather than an assertion.
CLEAN: list[Case] = [
    Case("binding free energy from a run that reached bulk",
         _complete_pmf, "proceeded",
         "the tail follows a free ligand's shape, so the reference is a "
         "reference"),
    Case("order parameters on a settled peptide",
         _order_parameters_on_a_settled_run, "proceeded",
         "the halves agree, so the values are not qualified"),
    Case("radial distribution within the box",
         _rdf_within_the_box, "proceeded",
         "the whole requested range lies inside half the box"),
    Case("mutation named against the residue that is there",
         _mutation_that_matches, "proceeded",
         "the structure holds what the mutation says it holds"),
    Case("metadynamics that crossed and settled",
         _metadynamics_that_crossed_and_settled, "proceeded",
         "the system visited both basins repeatedly and the hills "
         "flattened"),
    Case("a run long against its correlation time",
         _a_run_long_against_its_correlation, "proceeded",
         "the halves agree, so the effective sample count means what it "
         "says"),
    Case("a box that already fits its cutoff",
         _a_box_that_already_fits, "proceeded",
         "no growing was needed, so nothing is reported about it"),
]
