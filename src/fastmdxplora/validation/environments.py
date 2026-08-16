"""One study description, run in several places, compared field by field.

The claim is that a configuration is the study: hand it to a container on a
cluster, a conda installation on a workstation and a wheel from PyPI, and
the same thing happens. That is checkable, and it is checkable more
precisely than by comparing results, because two of the three artefacts a
run leaves are supposed to be identical regardless of where it ran.

**The resolved configuration must match exactly.** It is every default
filled in, and a default that resolves differently on another machine is
the failure this file exists to catch: the study said one thing and two
machines heard two.

**The input digests must match exactly.** The manifest records the SHA-256
of the structure that entered the run. Two runs of one study that started
from different bytes are not two runs of one study, whatever their
configurations say.

**The results must not be expected to match exactly, and saying so is part
of the comparison.** Solvation places water from a generator whose state is
not fixed by the dynamics seed, so two runs of one configuration begin from
different coordinates and diverge from there. Reporting that as a failure
would make an honest limitation look like a defect; ignoring it would make
a real difference invisible. The observables are therefore compared against
the spread across environments rather than for equality, and the reason is
carried with the result.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

__all__ = [
    "ENVIRONMENTAL_KEYS",
    "compare_environments",
    "digest_of",
]

#: Keys whose values are about where a run happened rather than what it
#: was. A difference in these is not a difference in the study, and
#: reporting one as a discrepancy would bury the discrepancies that matter
#: under a list of paths.
ENVIRONMENTAL_KEYS = frozenset({
    "output_dir", "output", "run_id", "started_at", "finished_at",
    "hostname", "platform", "device_index", "precision", "threads",
    "version", "versions", "python", "duration_s", "seed_source",
})


def digest_of(path: str | Path) -> str:
    """SHA-256 of a file, for comparing inputs byte for byte."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _flatten(value: Any, prefix: str = "") -> dict[str, Any]:
    """A nested configuration as one level of dotted keys.

    Compared flat because a nested diff reports the outermost block that
    changed, and the useful answer is which setting did.
    """
    out: dict[str, Any] = {}
    if isinstance(value, dict):
        for key, item in value.items():
            out.update(_flatten(item, f"{prefix}{key}."))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            out.update(_flatten(item, f"{prefix}{index}."))
    else:
        out[prefix.rstrip(".")] = value
    return out


def _read_config(run_dir: Path) -> dict[str, Any]:
    for name in ("resolved_config.yml", "resolved_config.yaml"):
        path = run_dir / name
        if path.is_file():
            try:
                import yaml

                return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            except ImportError:  # pragma: no cover - yaml is a dependency
                return {}
    return {}


def _read_manifest(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "manifest.json"
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _input_digest(manifest: dict[str, Any]) -> str | None:
    """The structure's digest, wherever the manifest recorded it."""
    for key in ("input_sha256", "structure_sha256", "sha256"):
        if manifest.get(key):
            return str(manifest[key])
    structure = manifest.get("structure")
    if isinstance(structure, dict):
        for key in ("sha256", "digest"):
            if structure.get(key):
                return str(structure[key])
    return None


def compare_environments(
    runs: "dict[str, str | Path]",
) -> dict[str, Any]:
    """Compare one study as it ran in several environments.

    ``runs`` maps a name for each environment to that run's directory.
    """
    if len(runs) < 2:
        return {"agreed": None, "refused": (
            f"{len(runs)} environment(s) given, and a comparison across "
            "environments needs at least two."
        )}

    configs: dict[str, dict[str, Any]] = {}
    digests: dict[str, str | None] = {}
    missing: list[str] = []
    for name, directory in runs.items():
        path = Path(directory)
        config = _read_config(path)
        if not config:
            missing.append(name)
            continue
        configs[name] = _flatten(config)
        digests[name] = _input_digest(_read_manifest(path))

    if len(configs) < 2:
        return {"agreed": None, "refused": (
            "No resolved configuration was found for "
            + ", ".join(missing)
            + ". A run writes one; a directory without it is not a run this "
            "can compare."
        )}

    keys = set().union(*(set(c) for c in configs.values()))
    study_differences: dict[str, dict[str, Any]] = {}
    environmental: dict[str, dict[str, Any]] = {}
    for key in sorted(keys):
        values = {name: config.get(key) for name, config in configs.items()}
        if len({json.dumps(v, sort_keys=True, default=str)
                for v in values.values()}) == 1:
            continue
        leaf = key.split(".")[-1]
        if leaf in ENVIRONMENTAL_KEYS:
            environmental[key] = values
        else:
            study_differences[key] = values

    unique_digests = {d for d in digests.values() if d is not None}
    digests_agree = len(unique_digests) <= 1
    unknown_digests = [n for n, d in digests.items() if d is None]

    return {
        "environments": sorted(configs),
        "agreed": not study_differences and digests_agree,
        "study_differences": study_differences,
        "environmental_differences": environmental,
        "input_digest": {
            "agree": digests_agree,
            "digests": digests,
            "not_recorded_for": unknown_digests,
        },
        "missing_configuration": missing,
        "what_is_not_compared": (
            "Results. Solvation places water from a generator the dynamics "
            "seed does not fix, so two runs of one configuration start from "
            "different coordinates and diverge from there. Their observables "
            "should agree within the spread between replicas, not exactly, "
            "and exact reproduction needs the prepared system carried over "
            "rather than only the configuration."
        ),
        "refused": None,
    }
