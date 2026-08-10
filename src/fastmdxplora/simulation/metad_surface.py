"""A free energy surface from a metadynamics run, or a refusal saying why not.

Metadynamics wrote its hills and its trajectory of the collective variable and
stopped there. Nothing read them back, so the run produced PLUMED's output
files and no result: the module said "a run that has not converged has no free
energy" while offering no free energy either way, and somebody wanting a
surface ran ``plumed sum_hills`` themselves and got one with nothing attached.

Summing the hills is arithmetic. The part worth having is what comes with it.
A metadynamics surface is only a measurement if the bias has stopped growing
and the system has crossed the barriers more than once, and both of those are
answerable from the files the run already writes.

**The bias must have flattened.** In well-tempered metadynamics the deposited
hills shrink as a basin fills, and the height of the last hills is the direct
evidence. Hills still arriving at nearly their initial height mean the surface
is still being built, and reporting it is reporting a snapshot of a filling
process as though it were the shape of the landscape.

**The system must have come back.** A barrier crossed once has been observed
once. The height of that barrier then rests on a single event, and no amount
of sampling either side of it adds a second observation. Recrossings are what
make a barrier a measurement rather than an anecdote.

**And the surface must have stopped moving.** Built from three quarters of the
hills and from all of them, a converged surface gives the same answer. This is
the check the other two are proxies for, and the one that catches a run that
satisfies them both for the wrong reasons.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "Hills",
    "read_hills",
    "surface_from_hills",
    "recrossings",
    "compute_surface",
]

#: Below this, a barrier has been seen too few times to be a measurement
#: rather than an anecdote. A judgement, like the overlap an umbrella study
#: demands, so a study can set its own.
MINIMUM_RECROSSINGS = 4

#: The last hills, as a fraction of the first, below which the bias is taken to
#: have flattened. Well-tempered deposition decays as exp(-V/kΔT), so a tenth
#: is already deep into the tail; more than that and the basin is still filling.
SETTLED_HEIGHT_FRACTION = 0.1

#: How much the surface may still move, in kJ/mol, between three quarters of
#: the hills and all of them. Roughly thermal energy at room temperature: a
#: feature that moves by less than kT is not a feature that has moved.
SETTLED_DRIFT_KJMOL = 2.5

#: How far above its minimum a surface is still judged for drift. Beyond
#: this the free energy is estimated from too few visits to be stable, and
#: including it means reporting the reliability of the least-sampled point on
#: the grid rather than of the surface anyone reads. Twenty kJ/mol is about
#: eight times RT at 300 K: a state that rare is not one a simulation is
#: measuring, and the barrier heights above it are read off the shape rather
#: than trusted point by point.
DRIFT_CEILING_KJMOL = 20.0


@dataclass(frozen=True)
class Hills:
    """What a metadynamics run deposited, read back from its HILLS file."""

    time_ps: np.ndarray
    centre: np.ndarray
    sigma: np.ndarray
    height: np.ndarray
    #: (T + ΔT) / T. One means the run was not well-tempered, and the
    #: relationship between bias and free energy is a different one.
    bias_factor: float

    def __len__(self) -> int:
        return int(self.centre.size)


def read_hills(path: str | Path) -> Hills:
    """Read a PLUMED HILLS file for a one-dimensional collective variable.

    The columns are time, the variable, sigma, height and the bias factor.
    Comment lines carry the field names and are skipped: this reads by
    position, which is what PLUMED guarantees for a single variable.
    """
    rows: list[list[float]] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        rows.append([float(value) for value in parts[:5]])

    if not rows:
        raise ValueError(
            f"{path} holds no hills. A metadynamics run that deposited "
            "nothing has not biased anything, and there is no surface to "
            "build from it."
        )

    table = np.array(rows, dtype=float)
    factor = float(table[0, 4]) if table.shape[1] > 4 else 1.0
    return Hills(
        time_ps=table[:, 0],
        centre=table[:, 1],
        sigma=table[:, 2],
        height=table[:, 3],
        bias_factor=factor,
    )


def surface_from_hills(
    hills: Hills, grid: np.ndarray, *, upto: int | None = None,
    periodic: bool = False,
) -> np.ndarray:
    """The free energy on ``grid``, from the first ``upto`` hills.

    The bias is the sum of the deposited Gaussians. For a well-tempered run
    the bias approaches -(1 - 1/γ) F, so the free energy is the bias scaled by
    γ/(γ - 1) and negated. Without tempering the bias approaches -F directly,
    and the scaling would divide by zero, so that case is its own.

    Shifted so the lowest point is zero, because only differences mean
    anything: the absolute value of a free energy from metadynamics is set by
    where the sum happened to start.
    """
    stop = len(hills) if upto is None else int(upto)
    centres = hills.centre[:stop, None]
    sigmas = hills.sigma[:stop, None]
    heights = hills.height[:stop, None]

    # On a circle the separation is measured the short way round. A hill
    # deposited at +170 degrees is 15 degrees from a grid point at -175, not
    # 345: computed straight, it contributes nothing there, and the surface
    # grows an artificial wall exactly where the coordinate is continuous.
    separation = grid[None, :] - centres
    if periodic:
        separation = np.remainder(separation + np.pi, 2.0 * np.pi) - np.pi

    bias = np.sum(
        heights * np.exp(-(separation ** 2) / (2.0 * sigmas ** 2)),
        axis=0,
    )

    if hills.bias_factor > 1.0:
        scale = hills.bias_factor / (hills.bias_factor - 1.0)
    else:
        # Not well-tempered: the bias converges on the negative of the free
        # energy without rescaling.
        scale = 1.0

    surface = -scale * bias
    return surface - float(np.min(surface))


def basins(coordinate: np.ndarray, surface: np.ndarray,
           periodic: bool = False) -> tuple[float, float] | None:
    """The two deepest minima, which are the states a crossing goes between.

    Returns ``None`` where there is only one, because a coordinate with a
    single basin has nothing to cross to and a count of crossings would be a
    count of nothing.
    """
    finite = np.isfinite(surface)
    if finite.sum() < 3:
        return None
    where, values = coordinate[finite], surface[finite]

    if periodic and values.size > 2 and abs(
            (where[-1] - where[0]) - 2.0 * np.pi) < 1e-6:
        # A grid spanning -pi to pi inclusive names the same place twice, so
        # the wrap-around neighbour of the first point is the second to last,
        # not the last. Dropping the duplicate is simpler than special-casing
        # every comparison that follows.
        where, values = where[:-1], values[:-1]

    interior = [
        i for i in range(1, values.size - 1)
        if values[i] < values[i - 1] and values[i] < values[i + 1]
    ]
    if periodic and values.size > 2:
        # The ends are neighbours, so a minimum can sit across the join.
        if values[0] < values[1] and values[0] < values[-1]:
            interior.append(0)
        if values[-1] < values[-2] and values[-1] < values[0]:
            interior.append(values.size - 1)
    if len(interior) < 2:
        return None

    interior.sort(key=lambda i: values[i])
    return float(where[interior[0]]), float(where[interior[1]])


def transitions(values: np.ndarray, *, between: tuple[float, float],
                periodic: bool = False) -> int:
    """How many times the variable travelled from one basin to the other.

    Counted by which minimum each frame is nearer to, measured the short way
    round where the coordinate is a circle, with a frame counted as arrived
    only once it is within a quarter of the separation. A run rattling
    between them does not accumulate crossings it did not make.

    The band used to be placed a quarter of the way in from the extremes of
    the hill range, which is a statement about the grid rather than about the
    system. On a full turn that put the thresholds at -90 and +90 degrees,
    and a real study whose minima sat at -17 and 155 had one of them inside
    the dead band: the count was of excursions past an arbitrary line.
    """
    if values.size == 0:
        return 0

    def separation(a, b):
        d = np.asarray(a, dtype=float) - b
        if periodic:
            d = np.remainder(d + np.pi, 2.0 * np.pi) - np.pi
        return np.abs(d)

    first, second = between
    apart = float(separation(first, second))
    if apart <= 0:
        return 0
    close_enough = 0.25 * apart

    to_first = separation(values, first)
    to_second = separation(values, second)

    crossings = 0
    side: str | None = None
    for near_first, near_second in zip(to_first, to_second):
        if near_first <= close_enough:
            here = "first"
        elif near_second <= close_enough:
            here = "second"
        else:
            continue
        if side is not None and here != side:
            crossings += 1
        side = here
    return crossings


def recrossings(values: np.ndarray, *, low: float, high: float) -> int:
    """How many times the variable travelled from one side to the other.

    Counted with a hysteresis band: a crossing is only counted once the
    variable has reached the far region, so a run rattling around a threshold
    does not accumulate crossings it did not make. That is the difference
    between a barrier observed several times and a coordinate jittering at
    the top of one.
    """
    if values.size == 0:
        return 0

    crossings = 0
    side: str | None = None
    for value in values:
        if value <= low:
            here = "low"
        elif value >= high:
            here = "high"
        else:
            continue
        if side is not None and here != side:
            crossings += 1
        side = here
    return crossings


def compute_surface(
    hills_path: str | Path,
    colvar_values: np.ndarray | None = None,
    *,
    points: int = 200,
    minimum_recrossings: int = MINIMUM_RECROSSINGS,
    periodic: bool = False,
) -> dict[str, Any]:
    """A free energy surface, or a refusal saying what the run cannot support.

    ``colvar_values`` is the trajectory of the collective variable, which the
    run writes beside the hills. Without it the recrossing count cannot be
    made, and that is said rather than skipped -- a surface reported without
    knowing whether the system ever came back is the claim this exists to
    avoid.
    """
    hills = read_hills(hills_path)

    lowest = float(np.min(hills.centre))
    highest = float(np.max(hills.centre))
    if highest <= lowest:
        return {
            "surface": None,
            "grid": None,
            "refused": (
                "Every hill was deposited at the same value of the collective "
                "variable, so the run explored nothing. There is no surface "
                "along a coordinate that did not move."
            ),
        }

    if periodic:
        # The whole turn. A grid bounded by where hills happened to land
        # stops short of the arc they wrap around, and the coordinate does
        # not stop there.
        grid = np.linspace(-np.pi, np.pi, points)
    else:
        grid = np.linspace(lowest, highest, points)
    surface = surface_from_hills(hills, grid, periodic=periodic)
    earlier = surface_from_hills(
        hills, grid, upto=int(len(hills) * 0.75), periodic=periodic)

    # Judged where the surface means something, not everywhere on the grid.
    #
    # Both surfaces are reported against their own minimum, so they are
    # compared as shapes. The comparison used to run over every point, which
    # makes the number the worst point rather than the typical one -- and the
    # worst point is always the top of the highest barrier, estimated from a
    # handful of visits out of thousands of hills. A tripeptide's psi torsion
    # settled to 1.2 kJ/mol within 10 of its minimum and 2.3 within 20, while
    # the whole-grid figure read 5.4 from a single point at the top of a 65
    # kJ/mol barrier. On that measure no steep coordinate can ever pass,
    # however well its wells are resolved: the test could not say yes to a
    # correct answer.
    shape = surface - float(np.min(surface))
    earlier_shape = earlier - float(np.min(earlier))
    difference = np.abs(shape - earlier_shape)
    judged = shape <= DRIFT_CEILING_KJMOL
    if not judged.any():
        judged = np.ones_like(shape, dtype=bool)
    drift = float(np.max(difference[judged]))
    drift_over_the_whole_grid = float(np.max(difference))

    first = float(np.mean(hills.height[:max(1, len(hills) // 20)]))
    last = float(np.mean(hills.height[-max(1, len(hills) // 20):]))
    settled = last <= SETTLED_HEIGHT_FRACTION * first if first > 0 else False

    # Between the two deepest basins, which are the states a crossing goes
    # between. Falls back to a band placed a quarter in from the extremes of
    # the hill range where the surface has only one minimum -- there is
    # nothing to cross to, and the older measure at least says whether the
    # coordinate moved.
    span = highest - lowest
    low_edge = lowest + 0.25 * span
    high_edge = highest - 0.25 * span
    two_basins = basins(grid, surface, periodic=periodic)

    crossed: int | None = None
    if colvar_values is not None and colvar_values.size:
        sampled = np.asarray(colvar_values, dtype=float)
        if two_basins is not None:
            crossed = transitions(
                sampled, between=two_basins, periodic=periodic)
        else:
            crossed = recrossings(sampled, low=low_edge, high=high_edge)

    evidence: dict[str, Any] = {
        "basins": (None if two_basins is None
                   else [float(x) for x in two_basins]),
        "hills": len(hills),
        "bias_factor": hills.bias_factor,
        "first_hill_height_kjmol": first,
        "last_hill_height_kjmol": last,
        "drift_kjmol": drift,
        "drift_ceiling_kjmol": DRIFT_CEILING_KJMOL,
        "drift_over_the_whole_grid_kjmol": drift_over_the_whole_grid,
        "recrossings": crossed,
        # Named, because "recrossings" is a word whose definition changes the
        # number. Counting sign changes of the raw coordinate gave 97 on a run
        # this recorded as 61 -- a 60% difference, and the smaller number is
        # the honest one, since it only counts a crossing once the coordinate
        # has reached the far side. Anyone comparing their own count against
        # this one, and not knowing that, would draw the wrong conclusion
        # about their sampling.
        # It has to say which measure was used, because the two answer
        # different questions and the number alone does not distinguish them.
        "recrossings_definition": (
            (f"travel between the basins at {two_basins[0]:.3g} and "
             f"{two_basins[1]:.3g}, a frame counted as arrived once it is "
             "within a quarter of their separation, so a coordinate "
             "rattling between them does not accumulate crossings it never "
             "made")
            if two_basins is not None else
            (f"transitions between the regions below {low_edge:.3g} and "
             f"above {high_edge:.3g} -- the surface has only one minimum, so "
             "there is no second basin to travel to and this says whether "
             "the coordinate moved rather than whether it changed state")
        ),
        "barrier_kjmol": float(np.max(surface)),
    }

    reasons: list[str] = []
    if not settled:
        reasons.append(
            f"the hills are still arriving at {last:.2f} kJ/mol against "
            f"{first:.2f} at the start, so the bias is still filling the "
            "landscape rather than having flattened it -- the surface below "
            "is a snapshot of that filling"
        )
    if drift > SETTLED_DRIFT_KJMOL:
        reasons.append(
            f"the surface moved {drift:.1f} kJ/mol between three quarters of "
            f"the hills and all of them, against a tolerance of "
            f"{SETTLED_DRIFT_KJMOL:g}, so it has not stopped changing "
            f"(measured where the surface is within "
            f"{DRIFT_CEILING_KJMOL:g} kJ/mol of its minimum; over the whole "
            f"grid it moved {drift_over_the_whole_grid:.1f}, which is "
            "dominated by the top of the highest barrier and says little)"
        )
    if crossed is None:
        reasons.append(
            "the trajectory of the collective variable was not available, so "
            "there is no way to tell whether the system crossed the barriers "
            "more than once"
        )
    elif crossed < minimum_recrossings:
        reasons.append(
            f"the system crossed between the ends of the range {crossed} "
            f"time(s), against {minimum_recrossings} asked for -- a barrier "
            "crossed once has been observed once, and its height rests on "
            "that single event"
        )

    if reasons:
        return {
            # The surface, not None. The refusal says "the surface below is a
            # snapshot of that filling" and "the hills are on disk, so the
            # check can be made again" -- and then withheld it, so a run whose
            # wells had converged to half of RT left nothing to plot. "Not
            # converged" and "not available" are different claims, and only
            # the first one is true here. It is provisional, which is what the
            # refusal beside it is for.
            "surface": surface,
            "provisional": True,
            "grid": grid,
            "evidence": evidence,
            "refused": (
                "This run does not support a free energy surface: "
                + "; ".join(reasons) + ".\n\n"
                "Sampling for longer is the remedy for all three. The hills "
                "and the trajectory are on disk, so the check can be made "
                "again on a longer run without repeating this one."
            ),
        }

    return {
        "surface": surface,
        "provisional": False,
        "grid": grid,
        "evidence": evidence,
        "refused": None,
    }
