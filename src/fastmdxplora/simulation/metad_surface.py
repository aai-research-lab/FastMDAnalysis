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
    hills: Hills, grid: np.ndarray, *, upto: int | None = None
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

    bias = np.sum(
        heights * np.exp(-((grid[None, :] - centres) ** 2)
                         / (2.0 * sigmas ** 2)),
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

    grid = np.linspace(lowest, highest, points)
    surface = surface_from_hills(hills, grid)
    earlier = surface_from_hills(hills, grid, upto=int(len(hills) * 0.75))
    drift = float(np.max(np.abs(surface - earlier)))

    first = float(np.mean(hills.height[:max(1, len(hills) // 20)]))
    last = float(np.mean(hills.height[-max(1, len(hills) // 20):]))
    settled = last <= SETTLED_HEIGHT_FRACTION * first if first > 0 else False

    span = highest - lowest
    crossed: int | None = None
    if colvar_values is not None and colvar_values.size:
        crossed = recrossings(
            np.asarray(colvar_values, dtype=float),
            low=lowest + 0.25 * span, high=highest - 0.25 * span)

    evidence: dict[str, Any] = {
        "hills": len(hills),
        "bias_factor": hills.bias_factor,
        "first_hill_height_kjmol": first,
        "last_hill_height_kjmol": last,
        "drift_kjmol": drift,
        "recrossings": crossed,
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
            f"{SETTLED_DRIFT_KJMOL:g}, so it has not stopped changing"
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
            "surface": None,
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
        "grid": grid,
        "evidence": evidence,
        "refused": None,
    }
