"""One table from a campaign's members, and what their spread means.

A campaign leaves each member in its own directory with its own analyses.
Reading them is arithmetic; knowing what the spread across them measures is
not, and it is the whole of why this exists.

**Replicas and variants are opposites, and they look identical on disk.**
Members that differ only by random seed are repeats of one measurement, so
the spread of their means is the error on it. Members that differ by
system, mutation or parameter are different measurements, so the spread
between them is the result. Averaging the second kind reports a mutant
series' biology as noise; quoting the first kind as a finding reports noise
as biology. This tells them apart from the campaign's own sweep record
rather than from the shape of the numbers, and says which it decided.

**Replicas are what makes a reported uncertainty checkable.** Every recorded
mean carries a standard error estimated from one trajectory's
autocorrelation. Ten trajectories give the same quantity a second way, as
the spread of ten means, and the two should agree. Where the single-run
estimate is much the smaller, the run was too short to see its own
correlation time and the error bars it printed were too tight. That
comparison is only available to a seed sweep, which is why it is reported
only there.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "SEED_AXES",
    "aggregate_members",
    "member_directory",
    "read_member_findings",
]

#: Sweep axes that make members repeats rather than different studies.
#: Written out rather than matched on the word "seed", because a study
#: sweeping something merely named seed-like is not replicating anything.
SEED_AXES = frozenset({
    "simulation.random_seed",
    "random_seed",
    # Prepared from different random placements of hydrogens and ions,
    # and otherwise the same study: repeats of one measurement.
    "setup.random_seed",
    # Prepared from different models of one ensemble: the same molecule
    # started from different structures, which is what independent starts
    # are.
    "setup.model",
})

#: How far the single-run error estimate may sit from the spread of the
#: replicas' means before the estimate is called optimistic. A factor of
#: two: the estimate is a variance, so agreeing to better than a factor of
#: two on ten replicas is already at the limit of what ten can resolve.
CALIBRATION_FACTOR = 2.0


def member_directory(batch_dir: str | Path, run: dict[str, Any]) -> Path:
    """Where a member's results are: relative to the batch where the
    manifest says so and the folder is there, so a campaign that was moved
    or copied reads its own members; otherwise the path recorded. A
    campaign moved without that record found every member missing, or, with
    another campaign at the old path, read that one's."""
    root = Path(batch_dir)
    relative = run.get("output_dir_relative")
    if relative and (root / str(relative)).is_dir():
        return root / str(relative)
    run_id = run.get("run_id")
    if run_id and (root / "runs" / str(run_id)).is_dir():
        return root / "runs" / str(run_id)
    return Path(str(run.get("output_dir", "")))


def read_member_findings(run_dir: str | Path) -> dict[str, dict[str, Any]]:
    """Every analysis's findings in one member, keyed by analysis name."""
    root = Path(run_dir)
    out: dict[str, dict[str, Any]] = {}
    for options in sorted(root.glob("analysis/*/options.json")):
        try:
            record = json.loads(options.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        name = record.get("analysis") or options.parent.name
        findings = record.get("findings")
        if isinstance(findings, dict):
            out[str(name)] = findings
    return out


def _members_are_replicas(manifest: dict[str, Any]) -> tuple[bool, str]:
    """Whether the campaign repeated one study or compared several."""
    sweep = manifest.get("sweep") or {}
    axes = set(sweep)
    systems = {
        (run.get("system") or run.get("structure"))
        for run in manifest.get("runs", [])
    }
    systems.discard(None)

    if len(systems) > 1:
        return False, (
            f"the members cover {len(systems)} different systems, so the "
            "spread between them is what the campaign set out to compare"
        )
    if not axes:
        return False, (
            "the campaign has no sweep axis, so there is nothing saying the "
            "members are repeats of one another"
        )
    if axes <= SEED_AXES:
        return True, (
            "the members differ only by random seed, so they are repeats of "
            "one estimate and the spread of their means is its error"
        )
    return False, (
        "the members differ by "
        + ", ".join(sorted(axes))
        + ", so the spread between them is a result rather than an error"
    )


def _finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _share_the_start(root: Path, runs: list[dict[str, Any]]) -> None:
    """Replicas' means from the start they share, as the stopping rule takes
    them (:func:`fastmdxplora.simulation.stopping.share_the_start`): each
    run's own equilibration detection keeps a little of the relaxation they
    share from their one structure, and the replicas' spread and pooled
    mean then carry it. Written into each member's record, which says so."""
    from fastmdxplora.simulation.stopping import share_the_start

    members = [member_directory(root, run) for run in runs]
    names = set.intersection(*[
        {p.parent.name for p in member.glob("analysis/*/options.json")} for member in members])
    for name in sorted(names):
        share_the_start(members, name)


def aggregate_members(batch_dir: str | Path) -> dict[str, Any]:
    """Collect a campaign's members into one comparison.

    Returns the per-member values for every analysis that reported a
    recorded mean, and, where the members are replicas, the comparison
    between the error each run estimated for itself and the spread the
    replicas actually show.
    """
    root = Path(batch_dir)
    manifest_path = root / "batch_manifest.json"
    if not manifest_path.is_file():
        return {"members": [], "refused": (
            f"No batch_manifest.json in {root}, so there is no campaign here "
            "to aggregate. A single run is not a campaign."
        )}

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    runs = [r for r in manifest.get("runs", []) if r.get("status") == "ok"]
    if len(runs) < 2:
        return {"members": [], "refused": (
            f"{len(runs)} member(s) finished, and a comparison across "
            "members needs at least two."
        )}

    replicas, why = _members_are_replicas(manifest)
    if replicas:
        _share_the_start(root, runs)

    per_analysis: dict[str, list[dict[str, Any]]] = {}
    for run in runs:
        findings = read_member_findings(member_directory(root, run))
        for analysis, record in findings.items():
            mean = record.get("mean") or {}
            if "mean" not in mean:
                continue
            per_analysis.setdefault(analysis, []).append({
                "member": run.get("run_id"),
                "sweep_values": run.get("sweep_values") or {},
                "mean": float(mean["mean"]),
                # NaN where the run could not estimate one: None, as for a
                # run that recorded none (NaN is not JSON a page can read).
                "standard_error": _finite_or_none(mean.get("standard_error")),
                "effective_samples": mean.get("effective_samples"),
                # Carried through rather than summarised away: a member that
                # said its own number was not a measurement must not be
                # averaged into one silently.
                "not_a_measurement": mean.get("not_a_measurement"),
            })

    table: dict[str, Any] = {}
    for analysis, rows in sorted(per_analysis.items()):
        if len(rows) < 2:
            continue
        means = np.array([r["mean"] for r in rows], dtype=float)
        entry: dict[str, Any] = {
            "members": rows,
            "mean_of_means": float(np.mean(means)),
            "spread_of_means": float(np.std(means, ddof=1)),
            "spread_is": (
                "the error on one estimate" if replicas
                else "the difference between the studies compared"),
        }

        # Only an error a run could estimate is a prediction: a run too short
        # to give one records NaN (read as None above), and NaN compared
        # with anything is false, so the spread was said to be "within the
        # factor" of "nan" and the error "supported by repeating it".
        errors = [r["standard_error"] for r in rows
                  if r["standard_error"] is not None and r["standard_error"] > 0]
        if replicas and not errors:
            observed = float(np.std(means, ddof=1))
            entry["predicted_standard_error"] = None
            entry["observed_standard_error"] = observed
            entry["ratio_observed_to_predicted"] = None
            entry["calibration"] = (
                f"The replicas spread {observed:.3g}. No run gave an error "
                "of its own for this observable, so there is nothing to "
                "check the spread against: the spread is the only error "
                "there is."
            )
        elif replicas and errors:
            predicted = float(np.mean(errors))
            observed = float(np.std(means, ddof=1))
            entry["predicted_standard_error"] = predicted
            entry["observed_standard_error"] = observed
            entry["ratio_observed_to_predicted"] = float(observed / predicted)
            qualified = [r["member"] for r in rows if r["not_a_measurement"]]
            if qualified:
                entry["qualified_by_their_own_runs"] = qualified
            # A prediction from some of the runs is said as such.
            some = ("" if len(errors) == len(rows) else
                    f" (from the {len(errors)} of {len(rows)} runs that "
                    "gave an error of their own)")
            if observed > CALIBRATION_FACTOR * predicted:
                entry["calibration"] = (
                    f"The replicas spread {observed:.3g} where a single run "
                    f"predicted {predicted:.3g}{some}, more than the factor of "
                    f"{CALIBRATION_FACTOR:g} that ten replicas can resolve. "
                    "The single-run estimate reads the correlation time off "
                    "a trajectory too short to contain it, so the error bars "
                    "printed for this observable are too tight and the "
                    "replicas are the honest figure."
                    + (" The runs said as much themselves."
                       if qualified else "")
                )
            else:
                entry["calibration"] = (
                    f"The replicas spread {observed:.3g} against a predicted "
                    f"{predicted:.3g}{some}, within the factor of "
                    f"{CALIBRATION_FACTOR:g} that this many replicas can "
                    "resolve: the error this observable reports for itself "
                    "is supported by repeating it."
                )
        table[analysis] = entry

    return {
        "members": [r.get("run_id") for r in runs],
        "replicas": replicas,
        "why": why,
        "analyses": table,
        "refused": None if table else (
            "No analysis recorded a mean in two or more members, so "
            "there is nothing to compare across them."
        ),
    }
