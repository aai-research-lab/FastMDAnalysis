"""A free energy along a coordinate, from equilibrium sampling at each point.

Steered molecular dynamics drags a system along a coordinate and reports the
work, which depends on how fast it was dragged. Umbrella sampling does the
opposite: it holds the system at a series of positions, lets each equilibrate,
and recombines the sampling into a potential of mean force. Nothing is
hurried, so nothing is dissipated, and the result is a free energy rather than
an upper bound on one.

**It works only if adjacent windows sample overlapping ranges.** The
recombination stitches histograms together, and where two neighbours never
visit the same value there is nothing to stitch: the free energy on one side
cannot be placed relative to the other, and the curve through that gap is
interpolation dressed as a measurement. Overlap is therefore checked and a
gap is reported rather than bridged.

That check is the reason this module exists rather than a call to an external
WHAM program. Every implementation computes the same PMF; few of them say
when the windows could not support one.

**Where the starting structures come from matters.** Windows started from a
single structure are strained at the far end of the range, and the strain
relaxes into the sampling as drift. The usual source is a steered run: pull
once, take a frame near each window's centre, and each window begins near
where it will sit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "Window",
    "expand_umbrella",
    "UmbrellaPlan",
    "plan_windows",
    "windows_as_sweep",
    "collect_samples",
    "overlap_between",
    "compute_pmf",
]

#: Boltzmann's constant in kJ/mol/K, so a PMF comes out in kJ/mol.
KB_KJ = 0.008314462618


@dataclass(frozen=True)
class Window:
    """One position along the coordinate, held there by a spring."""

    index: int
    centre: float
    force_constant: float

    def plumed_lines(self, cv_lines: list[str]) -> str:
        return "\n".join(cv_lines + [
            "",
            f"restraint: RESTRAINT ARG=cv AT={self.centre:g} "
            f"KAPPA={self.force_constant:g}",
            "",
            "PRINT ARG=cv,restraint.bias STRIDE=100 FILE=COLVAR",
        ]) + "\n"


@dataclass(frozen=True)
class UmbrellaPlan:
    """Where the windows are, and how hard each holds."""

    windows: tuple[Window, ...]
    collective_variable: str
    #: Steps to discard at the start of each window before it counts. A
    #: window begins away from where it will settle, and counting that
    #: approach as sampling biases the histogram towards where it started.
    equilibration_steps: int = 0

    def as_record(self) -> dict[str, Any]:
        return {
            "collective_variable": self.collective_variable,
            "n_windows": len(self.windows),
            "from": self.windows[0].centre if self.windows else None,
            "to": self.windows[-1].centre if self.windows else None,
            "force_constant": (self.windows[0].force_constant
                               if self.windows else None),
            "equilibration_steps": self.equilibration_steps,
        }


def plan_windows(spec: dict[str, Any]) -> UmbrellaPlan:
    """Read an umbrella block into a set of windows.

    Positions may be given explicitly, or as a range and a count. The count
    is the decision that matters: too few and adjacent windows do not overlap,
    which no amount of sampling repairs.
    """
    variable = str(spec.get("collective_variable", "")).lower()
    if not variable:
        raise ValueError(
            "Umbrella sampling needs a `collective_variable` -- the coordinate "
            "the free energy is a function of."
        )

    force = spec.get("force_constant")
    if force is None:
        raise ValueError(
            "Umbrella sampling needs a `force_constant`: how firmly each "
            "window is held, in kJ/mol per unit of the variable squared. It "
            "sets how far a window wanders, and therefore whether neighbours "
            "overlap -- too stiff and they do not, too soft and the system "
            "escapes towards the nearest minimum. There is no value that is "
            "right for an arbitrary coordinate."
        )
    force = float(force)

    if "centres" in spec or "centers" in spec:
        centres = [float(c) for c in (spec.get("centres") or spec.get("centers"))]
        if len(centres) < 2:
            raise ValueError("Umbrella sampling needs at least two windows.")
    else:
        for key in ("from", "to", "n_windows"):
            if spec.get(key) is None:
                raise ValueError(
                    "Umbrella windows are given either as `centres`, or as "
                    "`from`, `to` and `n_windows`. "
                    f"`{key}` is missing."
                )
        count = int(spec["n_windows"])
        if count < 2:
            raise ValueError("Umbrella sampling needs at least two windows.")
        centres = list(np.linspace(float(spec["from"]), float(spec["to"]),
                                   count))

    return UmbrellaPlan(
        windows=tuple(Window(index=i, centre=c, force_constant=force)
                      for i, c in enumerate(centres)),
        collective_variable=variable,
        equilibration_steps=int(spec.get("equilibration_steps", 0)),
    )


def expand_umbrella(config: dict[str, Any]) -> dict[str, Any]:
    """Turn a config with an umbrella block into one with a run per window.

    An umbrella job is one system held at many positions, which is the shape
    the batch machinery already runs -- so the windows become `systems`
    entries differing in the position each holds, and the scheduling,
    parallelism and per-GPU pinning come from the code that already does
    those things.

    Returns the config unchanged where there is no umbrella block.
    """
    simulation = config.get("simulation") or {}
    spec = simulation.get("umbrella")
    if not spec:
        return config

    if config.get("systems") and len(config["systems"]) > 1:
        raise ValueError(
            "Umbrella sampling holds one system at many positions. This "
            f"config has {len(config['systems'])} systems, and a window set "
            "for each would be several separate free energies -- run them as "
            "separate studies so each has its own windows and its own "
            "overlap check."
        )

    plan = plan_windows(spec)
    base = dict(config.get("systems", [{}])[0]) if config.get("systems") else {}

    systems = []
    for window in plan.windows:
        entry = dict(base)
        entry["id"] = f"window_{window.index:02d}"
        # Everything the window needs to bias itself, resolved per run.
        entry["simulation"] = dict(entry.get("simulation") or {})
        entry["simulation"]["umbrella"] = dict(
            {k: v for k, v in spec.items()
             if k not in ("centres", "centers", "from", "to", "n_windows")},
            centre=window.centre,
            force_constant=window.force_constant,
            index=window.index,
        )
        systems.append(entry)

    expanded = dict(config)
    expanded["systems"] = systems
    # The block has been expanded; leaving it would have every run try to
    # expand it again.
    expanded["simulation"] = {k: v for k, v in simulation.items()
                              if k != "umbrella"}
    expanded["_umbrella_plan"] = plan
    return expanded


def windows_as_sweep(plan: UmbrellaPlan) -> list[dict[str, Any]]:
    """The windows as a sweep, which is the shape the batch machinery runs.

    An umbrella job is one system at many restraint positions, which is the
    same shape as a parameter sweep -- so the runs are expanded and scheduled
    by the machinery that already does that, rather than by a second one.
    """
    return [
        {
            "id": f"window_{w.index:02d}",
            "umbrella_centre": w.centre,
            "umbrella_force_constant": w.force_constant,
        }
        for w in plan.windows
    ]


def collect_samples(
    output_dir: Any,
    plan: UmbrellaPlan,
    *,
    equilibration_fraction: float = 0.2,
) -> dict[int, np.ndarray]:
    """Read each window's sampling back from its COLVAR file.

    The first part of each window is discarded. A window begins away from
    where it will settle, and counting the approach as sampling biases the
    histogram towards where the run started -- which is the one place the
    free energy is guaranteed not to be flat.
    """
    from pathlib import Path

    root = Path(output_dir)
    samples: dict[int, np.ndarray] = {}
    missing: list[str] = []

    for window in plan.windows:
        colvar = root / f"window_{window.index:02d}" / "simulation" / "COLVAR"
        if not colvar.is_file():
            missing.append(f"window_{window.index:02d}")
            continue
        rows = []
        for line in colvar.read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    rows.append(float(parts[1]))
                except ValueError:
                    continue
        if not rows:
            missing.append(f"window_{window.index:02d} (empty)")
            continue
        values = np.asarray(rows)
        cut = int(len(values) * equilibration_fraction)
        samples[window.index] = values[cut:]

    if missing:
        raise FileNotFoundError(
            "These windows produced no sampling: " + ", ".join(missing) + ". "
            "A free energy cannot be computed from a partial set, because "
            "the windows either side of a missing one have nothing between "
            "them to stitch through."
        )
    return samples


def overlap_between(a: np.ndarray, b: np.ndarray, bins: int = 50) -> float:
    """How much two windows' sampling shares ground, from 0 to 1.

    The overlap coefficient: the area shared by two normalised histograms.
    Zero means they never visited the same value, and no recombination can
    place one relative to the other.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return 0.0

    lo = min(a.min(), b.min())
    hi = max(a.max(), b.max())
    if hi <= lo:
        return 1.0 if abs(a.mean() - b.mean()) < 1e-12 else 0.0

    edges = np.linspace(lo, hi, bins + 1)
    pa, _ = np.histogram(a, bins=edges, density=True)
    pb, _ = np.histogram(b, bins=edges, density=True)
    width = edges[1] - edges[0]
    return float(np.minimum(pa, pb).sum() * width)


def compute_pmf(
    samples: dict[int, np.ndarray],
    plan: UmbrellaPlan,
    *,
    temperature_K: float = 300.0,
    bins: int = 60,
    minimum_overlap: float = 0.03,
) -> dict[str, Any]:
    """A potential of mean force, or a refusal saying why not.

    ``samples`` maps a window's index to the collective-variable values it
    sampled, after equilibration has been discarded.

    The overlap between neighbours is checked first. Where a pair does not
    share ground, the free energy on one side cannot be placed relative to the
    other, and a curve drawn through the gap is interpolation presented as a
    measurement.
    """
    ordered = [w for w in plan.windows if w.index in samples]
    if len(ordered) < 2:
        raise ValueError(
            "A potential of mean force needs at least two windows with "
            f"sampling in them; {len(ordered)} were given."
        )

    gaps = []
    overlaps = []
    for left, right in zip(ordered, ordered[1:]):
        shared = overlap_between(samples[left.index], samples[right.index])
        overlaps.append({
            "between": [left.index, right.index],
            "centres": [left.centre, right.centre],
            "overlap": shared,
        })
        if shared < minimum_overlap:
            gaps.append((left, right, shared))

    if gaps:
        described = "; ".join(
            f"windows {a.index} and {b.index} (at {a.centre:g} and "
            f"{b.centre:g}) share {s:.1%}"
            for a, b, s in gaps
        )
        return {
            "pmf": None,
            "overlaps": overlaps,
            "refused": (
                "Adjacent windows do not overlap, so no free energy can be "
                f"computed across the gap: {described}. Recombination "
                "stitches histograms together, and where two neighbours never "
                "visit the same value there is nothing to stitch -- a curve "
                "drawn through the gap would be interpolation presented as a "
                "measurement.\n\n"
                "More windows between them, or a softer force constant so "
                "each wanders further, will close it. Sampling for longer "
                "will not."
            ),
        }

    # WHAM, iterated to self-consistency. MBAR is better where it is
    # available, and this needs no dependency for the common case.
    kT = KB_KJ * float(temperature_K)
    everything = np.concatenate([samples[w.index] for w in ordered])
    edges = np.linspace(everything.min(), everything.max(), bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    counts = np.array([np.histogram(samples[w.index], bins=edges)[0]
                       for w in ordered], dtype=float)
    n_per_window = counts.sum(axis=1)
    bias = np.array([
        0.5 * w.force_constant * (centres - w.centre) ** 2 for w in ordered
    ])

    free_energies = np.zeros(len(ordered))
    for _ in range(2000):
        weights = np.exp((free_energies[:, None] - bias) / kT)
        denominator = (n_per_window[:, None] * weights).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            probability = counts.sum(axis=0) / denominator
        probability = np.nan_to_num(probability)
        updated = -kT * np.log(
            np.clip((probability[None, :] * np.exp(-bias / kT)).sum(axis=1),
                    1e-300, None))
        updated -= updated[0]
        if np.max(np.abs(updated - free_energies)) < 1e-6:
            free_energies = updated
            break
        free_energies = updated

    with np.errstate(divide="ignore"):
        pmf = -kT * np.log(np.clip(probability, 1e-300, None))
    pmf -= pmf[np.isfinite(pmf)].min()

    return {
        "pmf": {"coordinate": centres.tolist(), "free_energy_kjmol": pmf.tolist()},
        "overlaps": overlaps,
        "refused": None,
        "temperature_K": float(temperature_K),
        "n_windows": len(ordered),
    }
