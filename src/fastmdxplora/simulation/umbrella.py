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
    "plan_from_expanded",
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
    #: How much two neighbours must share before their free energies can be
    #: placed relative to one another. Three per cent is enough to stitch and
    #: thin: on a real study, pairs at seven per cent passed while a reader
    #: might reasonably want fifteen. It is a judgement about how much
    #: evidence a joint needs, so it belongs to whoever is making the claim.
    minimum_overlap: float = 0.03
    #: How many values a window must have recorded before its histogram means
    #: anything. An overlap is the area two histograms share, and a histogram
    #: from tens of points is mostly noise -- so a run short enough to be a
    #: smoke test will produce overlaps, and gaps, that are arithmetic rather
    #: than evidence. Like the overlap threshold, it is a judgement, so a
    #: study can set its own.
    minimum_samples: int = 200

    def as_record(self) -> dict[str, Any]:
        return {
            "collective_variable": self.collective_variable,
            "n_windows": len(self.windows),
            "from": self.windows[0].centre if self.windows else None,
            "to": self.windows[-1].centre if self.windows else None,
            "force_constant": (self.windows[0].force_constant
                               if self.windows else None),
            "equilibration_steps": self.equilibration_steps,
            "minimum_overlap": self.minimum_overlap,
            "minimum_samples": self.minimum_samples,
        }


#: Every key the umbrella block reads. Declared beside the reader so the two
#: cannot drift: a setting added below without being added here is refused,
#: which is a loud failure in a test rather than a quiet one in a study.
_UMBRELLA_OWN_KEYS: frozenset[str] = frozenset({
    "force_constant", "centres", "centers", "from", "to", "n_windows",
    "equilibration_steps", "minimum_overlap", "minimum_samples",
    # Written by the expansion onto each window, and read back from it.
    "centre", "index",
})


def _accepted_keys() -> frozenset[str]:
    """Umbrella's own settings, plus whatever names a collective variable.

    Taken from the collective-variable layer rather than copied, because the
    layer is shared -- a study biasing a ligand's depth in a bilayer names the
    bilayer the same way whichever method does the biasing, and a second copy
    of that list would go stale the next time a variable is added.
    """
    from fastmdxplora.simulation.metadynamics import COLLECTIVE_VARIABLE_KEYS

    return _UMBRELLA_OWN_KEYS | COLLECTIVE_VARIABLE_KEYS


def check_umbrella_keys(spec: dict[str, Any]) -> None:
    """Refuse a setting the block does not have.

    Nothing else checked this. `minimum_ovelap: 0.15` -- one letter -- was
    accepted, ignored, and the study stitched at the three per cent default
    while its author believed it had demanded fifteen. The run completes and
    produces a curve, and the guard that was supposed to stand behind that
    curve never existed.

    Every refusal this module makes can be switched off by a typo, so the
    typo is what has to be caught.
    """
    from fastmdxplora.config.loader import ConfigError, _suggest

    accepted = _accepted_keys()
    unknown = sorted(set(spec) - accepted)
    if not unknown:
        return
    named = ", ".join(
        f"'{key}'{_suggest(key, set(accepted))}" for key in unknown)
    raise ConfigError(
        f"Unknown umbrella setting{'s' if len(unknown) > 1 else ''}: {named}. "
        "Accepted: " + ", ".join(sorted(accepted - {"centre", "index"})) + "."
    )


def plan_windows(spec: dict[str, Any]) -> UmbrellaPlan:
    """Read an umbrella block into a set of windows.

    Positions may be given explicitly, or as a range and a count. The count
    is the decision that matters: too few and adjacent windows do not overlap,
    which no amount of sampling repairs.
    """
    check_umbrella_keys(spec)

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
        # Plain floats: numpy scalars serialise as "!!python/object" in YAML
        # and are not JSON, so a config written back out would not reload.
        centres = [float(c) for c in np.linspace(
            float(spec["from"]), float(spec["to"]), count)]

    return UmbrellaPlan(
        windows=tuple(Window(index=i, centre=c, force_constant=force)
                      for i, c in enumerate(centres)),
        collective_variable=variable,
        equilibration_steps=int(spec.get("equilibration_steps", 0)),
        minimum_overlap=float(spec.get("minimum_overlap", 0.03)),
        minimum_samples=int(spec.get("minimum_samples", 200)),
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
        # Every shared simulation setting, then this window's own. A
        # per-system block replaces the top-level one rather than merging, so
        # a block holding only the umbrella settings silently discarded the
        # step counts, the timestep and everything else the study asked for.
        merged = {k: v for k, v in simulation.items() if k != "umbrella"}
        merged.update(entry.get("simulation") or {})
        entry["simulation"] = merged
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
    # The plan is deliberately not stashed in the config. It was, and
    # validation rejected the extra key -- correctly, since a config is what
    # a user wrote and not a place to hide state. Each window carries its own
    # block, so the set can be rebuilt from the runs when they finish.
    return expanded


def plan_from_expanded(config: dict[str, Any]) -> UmbrellaPlan | None:
    """Rebuild the window set from a config that has already been expanded.

    Returns ``None`` where this is not an umbrella study.
    """
    systems = config.get("systems") or []
    windows = []
    variable = ""
    equilibration = 0
    minimum = 0.03
    fewest = 200
    for entry in systems:
        block = ((entry.get("simulation") or {}).get("umbrella")
                 if isinstance(entry, dict) else None)
        if not block or block.get("centre") is None:
            continue
        variable = str(block.get("collective_variable", variable))
        equilibration = int(block.get("equilibration_steps", equilibration))
        minimum = float(block.get("minimum_overlap", minimum))
        fewest = int(block.get("minimum_samples", fewest))
        windows.append(Window(
            index=int(block.get("index", len(windows))),
            centre=float(block["centre"]),
            force_constant=float(block["force_constant"]),
        ))
    if len(windows) < 2:
        return None
    return UmbrellaPlan(
        windows=tuple(sorted(windows, key=lambda w: w.index)),
        collective_variable=variable,
        equilibration_steps=equilibration,
        minimum_overlap=minimum,
        minimum_samples=fewest,
    )


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
    directories: dict[int, Any],
    *,
    equilibration_fraction: float = 0.2,
) -> dict[int, np.ndarray]:
    """Read each window's sampling back from its COLVAR file.

    Takes the directories rather than working them out. The first version
    guessed -- ``<output>/window_00/simulation/COLVAR`` -- and the runs are
    actually at ``<output>/runs/window-00/simulation/COLVAR``: under a
    ``runs`` directory, with the identifier slugged. Two mistakes in one
    path, neither visible until a real study finished and found nothing.
    The caller knows where it put things.

    The first part of each window is discarded. A window begins away from
    where it will settle, and counting the approach as sampling biases the
    histogram towards where the run started -- which is the one place the
    free energy is guaranteed not to be flat.
    """
    from pathlib import Path

    samples: dict[int, np.ndarray] = {}
    missing: list[str] = []

    for index, directory in sorted(directories.items()):
        colvar = Path(directory) / "simulation" / "COLVAR"
        if not colvar.is_file():
            missing.append(f"window {index} ({colvar})")
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
            missing.append(f"window {index} (empty COLVAR)")
            continue
        values = np.asarray(rows)
        cut = int(len(values) * equilibration_fraction)
        samples[index] = values[cut:]

    if missing:
        raise FileNotFoundError(
            "These windows produced no sampling: " + "; ".join(missing) + ". "
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


def windows_that_drifted(
    samples: dict[int, np.ndarray], plan: "UmbrellaPlan"
) -> list[dict[str, Any]]:
    """Windows whose sampling is not where the restraint was told to hold it.

    A window is a restraint at a position, and the histogram it produces is
    supposed to sit around that position. When it does not, the restraint has
    lost: usually because every window started from the same structure and the
    spring never dragged it out, so the far windows relax back towards the
    bound state and pile up on each other.

    That is a different failure from a genuine gap, and it wants the opposite
    remedy. A softer force constant -- the advice for windows that do not
    reach each other -- lets a window that has already slipped slip further.

    "Not where it was told" is measured against the spacing to its neighbour,
    because that is the distance the plan intends between one window and the
    next: sampling further from its own centre than half that is sampling
    somewhere another window was supposed to be.
    """
    ordered = [w for w in plan.windows if w.index in samples]
    drifted: list[dict[str, Any]] = []
    for position, window in enumerate(ordered):
        neighbours = [
            abs(other.centre - window.centre)
            for offset in (-1, 1)
            if 0 <= position + offset < len(ordered)
            for other in [ordered[position + offset]]
        ]
        if not neighbours:
            continue
        allowed = 0.5 * min(neighbours)
        sat_at = float(np.median(samples[window.index]))
        away = abs(sat_at - window.centre)
        if away > allowed:
            drifted.append({
                "window": window.index,
                "centre": window.centre,
                "sampled_at": sat_at,
                "away_by": away,
            })
    return drifted


def windows_with_too_little_sampling(
    samples: dict[int, np.ndarray], plan: "UmbrellaPlan"
) -> list[dict[str, Any]]:
    """Windows whose histograms are too sparse to mean anything.

    An overlap is the area two histograms share. From tens of points a
    histogram is mostly noise, so a run short enough to be a smoke test
    produces overlaps -- and gaps -- that are arithmetic rather than
    evidence. A free energy stitched from them would carry that noise as
    structure and look exactly like a result.
    """
    thin = [
        {"window": window.index, "samples": int(len(samples[window.index]))}
        for window in plan.windows
        if window.index in samples
        and len(samples[window.index]) < plan.minimum_samples
    ]
    return thin


def compute_pmf(
    samples: dict[int, np.ndarray],
    plan: UmbrellaPlan,
    *,
    temperature_K: float = 300.0,
    bins: int = 60,
    minimum_overlap: float | None = None,
) -> dict[str, Any]:
    """A potential of mean force, or a refusal saying why not.

    ``samples`` maps a window's index to the collective-variable values it
    sampled, after equilibration has been discarded.

    The overlap between neighbours is checked first. Where a pair does not
    share ground, the free energy on one side cannot be placed relative to the
    other, and a curve drawn through the gap is interpolation presented as a
    measurement.
    """
    # The plan carries it, so a study's own threshold applies wherever the
    # recombination happens rather than only where somebody remembered to
    # pass it.
    if minimum_overlap is None:
        minimum_overlap = plan.minimum_overlap

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

    drifted = windows_that_drifted(samples, plan)
    thin = windows_with_too_little_sampling(samples, plan)

    if thin:
        # Said before the gap and before the drift, because both are read off
        # histograms that this run did not fill. A refusal naming a specific
        # pair would be a precise claim resting on tens of points.
        fewest = min(t["samples"] for t in thin)
        reason = (
            f"{len(thin)} of {len(ordered)} windows recorded fewer than "
            f"{plan.minimum_samples} values after equilibration was discarded "
            f"-- the thinnest has {fewest}. An overlap is the area two "
            "histograms share, and a histogram from tens of points is mostly "
            "noise, so the overlaps below are arithmetic rather than "
            "evidence. A free energy stitched from them would carry that "
            "noise as structure and look exactly like a result.\n\n"
        )
        if drifted:
            where = "; ".join(
                f"window {d['window']} was held at {d['centre']:g} and "
                f"sampled around {d['sampled_at']:.3g}"
                for d in drifted
            )
            reason += (
                f"{len(drifted)} of them are also not at their centres: "
                f"{where}. A run this short explains that on its own -- a "
                "restraint needs time to pull a system to where it is held -- "
                "so whether the windows are too softly held cannot be told "
                "apart from a run that ended before they arrived. Sample for "
                "longer first, and read this again.\n\n"
            )
        reason += (
            f"`minimum_samples` ({plan.minimum_samples}) is a judgement about "
            "how much evidence a histogram needs, and a study that wants a "
            "smoke test to reach the recombination can lower it."
        )
        return {
            "pmf": None,
            "overlaps": overlaps,
            "drifted": drifted,
            "thin": thin,
            "refused": reason,
        }

    if gaps:
        described = "; ".join(
            f"windows {a.index} and {b.index} (at {a.centre:g} and "
            f"{b.centre:g}) share {s:.1%}"
            for a, b, s in gaps
        )
        reason = (
            "Adjacent windows do not overlap, so no free energy can be "
            f"computed across the gap: {described}. Recombination "
            "stitches histograms together, and where two neighbours never "
            "visit the same value there is nothing to stitch -- a curve "
            "drawn through the gap would be interpolation presented as a "
            "measurement.\n\n"
        )
        if drifted:
            # The windows are not where they were told to be, so the gap is
            # a symptom. Saying "use a softer force constant" here would make
            # it worse, which is what the generic advice would have said.
            where = "; ".join(
                f"window {d['window']} was held at {d['centre']:g} and "
                f"sampled around {d['sampled_at']:.3g}"
                for d in drifted
            )
            reason += (
                f"{len(drifted)} of {len(ordered)} windows did not sample "
                f"where they were held: {where}. That is the gap's cause "
                "rather than a shortage of windows -- a restraint that has "
                "lost its window leaves the ground between them unvisited "
                "however many more are added.\n\n"
                "Windows started from one structure are strained at the far "
                "end of the range and relax back towards the bound state. "
                "Seed each window from a steered run near its own centre, or "
                "hold them harder with a larger `force_constant`. A softer "
                "one will make this worse."
            )
        else:
            reason += (
                f"Each pair must share at least {minimum_overlap:.0%} "
                "(`minimum_overlap`). More windows between them, or a softer "
                "force constant so each wanders further, will close a gap. "
                "Sampling for longer will not."
            )
        return {
            "pmf": None,
            "overlaps": overlaps,
            "drifted": drifted,
            "thin": thin,
            "refused": reason,
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
