"""Config file loading, validation, and merging.

Loads a YAML configuration, validates it strictly against the schema
registry (:mod:`fastmdxplora.config.schema`), and produces a normalized
structure ready to drive :class:`fastmdxplora.FastMDXplora`.

Validation is strict by design: unknown keys raise
:class:`ConfigError` with a did-you-mean suggestion, and values whose
type doesn't match the schema raise with a clear message. A typo'd
config that silently runs with defaults is the worst failure mode in
science (you think you set ``ph: 7.4``, you actually ran the default,
and your results are subtly wrong with no indication why) — so we never
silently ignore.

Override precedence (highest wins):

  1. Explicit flags / kwargs supplied at call time
  2. Values in the config file
  3. Built-in phase defaults

This module implements (1) beating (2). The phase ``run()`` functions
implement (2) beating (3) via their own ``DEFAULTS`` tables.
"""

from __future__ import annotations

import difflib
from pathlib import Path
from typing import Any

from fastmdxplora.config.schema import (
    PHASE_KEYS,
    PHASE_SCHEMAS,
    TOP_LEVEL,
    TOP_LEVEL_KEYS,
    PhaseSchema,
)

VALID_PHASES = set(PHASE_KEYS)


class ConfigError(ValueError):
    """Raised for any problem loading or validating a config file."""


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_config_file(path: str | Path) -> dict[str, Any]:
    """Read and parse a YAML config file. Returns the raw dict.

    Raises
    ------
    ConfigError
        If the file is missing, unreadable, not valid YAML, or doesn't
        parse to a mapping at the top level.
    """
    import yaml

    p = Path(path)
    if not p.exists():
        raise ConfigError(f"Config file not found: {p}")
    try:
        with p.open(encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        raise ConfigError(f"Failed to parse YAML in {p}: {exc}") from exc

    if data is None:
        # Empty file — treat as empty config (valid; everything defaults)
        return {}
    if not isinstance(data, dict):
        raise ConfigError(
            f"Config file {p} must contain a YAML mapping at the top level, "
            f"got {type(data).__name__}."
        )
    return data


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _suggest(key: str, valid: set[str]) -> str:
    """Return a ' (did you mean X?)' suffix if a close match exists.

    Case-insensitive: a pure case mismatch (``pH`` vs ``ph``) is one of
    the most common config typos and short keys fall below difflib's
    default ratio when case differs, so we check case-folded matches
    first, then fall back to fuzzy matching.
    """
    # Exact case-insensitive match first (handles pH -> ph, PH -> ph, etc.)
    lower_map = {v.lower(): v for v in valid}
    if key.lower() in lower_map and lower_map[key.lower()] != key:
        return f" (did you mean '{lower_map[key.lower()]}'?)"
    # Fuzzy match on the case-folded forms
    matches = difflib.get_close_matches(
        key.lower(), list(lower_map.keys()), n=1, cutoff=0.6
    )
    if matches:
        return f" (did you mean '{lower_map[matches[0]]}'?)"
    return ""


def _type_name(t: type | tuple[type, ...]) -> str:
    if isinstance(t, tuple):
        return " or ".join(_friendly_type(x) for x in t)
    return _friendly_type(t)


def _friendly_type(t: type) -> str:
    return {
        str: "string",
        int: "integer",
        float: "number",
        bool: "true/false",
        list: "list",
        dict: "mapping",
    }.get(t, t.__name__)


def _check_type(value: Any, expected: type | tuple[type, ...]) -> bool:
    """Type check with YAML-friendly coercions.

    YAML parses ``1`` as int and ``1.0`` as float. A field declared
    ``float`` should accept an int (``temperature_K: 300`` is fine), so
    we accept int wherever float is allowed. We also reject bool where
    int/float is expected (YAML ``true`` is a Python bool, which is an
    int subclass — without this guard ``ph: true`` would pass an int
    check).
    """
    # bool is a subclass of int — handle it explicitly so numeric fields
    # don't silently accept booleans.
    allowed = expected if isinstance(expected, tuple) else (expected,)
    if isinstance(value, bool):
        return bool in allowed
    # int acceptable wherever float is expected
    if isinstance(value, int) and (float in allowed or int in allowed):
        return True
    if isinstance(value, float) and float in allowed:
        return True
    return isinstance(value, allowed)


def _validate_block(
    block: dict[str, Any],
    schema: PhaseSchema,
    *,
    context: str,
) -> None:
    """Validate one phase block (or the top-level scalars) against a schema."""
    valid = schema.field_names()
    for key, value in block.items():
        if key not in valid:
            raise ConfigError(
                f"Unknown {context} option '{key}'{_suggest(key, valid)}. "
                f"Valid options: {', '.join(sorted(valid))}."
            )
        fld = schema.get(key)
        assert fld is not None  # guaranteed by the membership check
        if value is None:
            # Explicit null is allowed — means "use the default".
            continue
        if not _check_type(value, fld.type):
            raise ConfigError(
                f"{context} option '{key}' should be {_type_name(fld.type)}, "
                f"got {type(value).__name__} ({value!r})."
            )


def _validate_phase_list(value: Any, *, field_name: str) -> None:
    """Validate top-level include/exclude phase names after type checking."""
    if value is None:
        return
    unknown = [phase for phase in value if not isinstance(phase, str) or phase not in VALID_PHASES]
    if unknown:
        valid = ", ".join(PHASE_KEYS)
        raise ConfigError(
            f"Top-level '{field_name}' contains unknown phase(s): "
            f"{', '.join(repr(p) for p in unknown)}. Valid phases: {valid}."
        )


def validate_config(data: dict[str, Any], *, require_systems: bool = False) -> None:
    """Strictly validate a parsed config dict against the schema.

    Parameters
    ----------
    data : dict
        The parsed config.
    require_systems : bool, default False
        If True, the config must define a non-empty ``systems`` list.
        The full pipeline (CLI / BatchExplorer) sets this; unit tests that
        validate fragments leave it False.

    Raises
    ------
    ConfigError
        On unknown top-level keys, unknown per-phase keys, type
        mismatches, mutually-exclusive include/exclude, a missing
        ``systems`` list (when required), or a malformed execution block.
    """
    # Top-level keys: scalar fields + phase block names
    for key in data:
        if key not in TOP_LEVEL_KEYS:
            raise ConfigError(
                f"Unknown top-level key '{key}'{_suggest(key, TOP_LEVEL_KEYS)}. "
                f"Valid: {', '.join(sorted(TOP_LEVEL_KEYS))}."
            )

    # Validate the top-level scalar fields (output, verbose, include,
    # exclude). Phase blocks and batch keys (systems/sweep/execution)
    # have their own structure and are validated separately below.
    from fastmdxplora.config.schema import BATCH_KEYS, EXECUTION

    top_scalars = {
        k: v for k, v in data.items()
        if k not in PHASE_KEYS and k not in BATCH_KEYS
    }
    _validate_block(top_scalars, TOP_LEVEL, context="top-level")

    # include / exclude mutual exclusion
    if data.get("include") and data.get("exclude"):
        raise ConfigError(
            "Config sets both 'include' and 'exclude' at the top level; "
            "they are mutually exclusive."
        )
    _validate_phase_list(data.get("include"), field_name="include")
    _validate_phase_list(data.get("exclude"), field_name="exclude")

    # `systems` is the canonical (and required) way to specify input.
    if require_systems and not data.get("systems"):
        raise ConfigError(
            "Config must define a `systems:` list (the canonical way to "
            "specify input). Even a single system goes in the list, e.g.\n"
            "  systems:\n"
            "    - {id: protein1, system: protein.pdb}"
        )

    # Validate batch keys (systems / sweep) via the batch layer.
    if data.get("systems") is not None or data.get("sweep") is not None:
        from fastmdxplora.batch.sweep import (
            SweepError,
            normalize_summary_for_validation,
        )
        try:
            normalize_summary_for_validation(data)
        except SweepError as exc:
            raise ConfigError(str(exc)) from exc

    # Validate the execution block (parallelism settings)
    if data.get("execution") is not None:
        execution = data["execution"]
        if not isinstance(execution, dict):
            raise ConfigError(
                f"The 'execution' block must be a mapping, "
                f"got {type(execution).__name__}."
            )
        _validate_block(execution, EXECUTION, context="execution")
        mode = execution.get("mode")
        if mode is not None and mode not in ("sequential", "parallel"):
            raise ConfigError(
                f"execution.mode must be 'sequential' or 'parallel', "
                f"got {mode!r}."
            )

    # Validate each per-phase block
    for phase in PHASE_KEYS:
        if phase not in data:
            continue
        block = data[phase]
        if not isinstance(block, dict):
            raise ConfigError(
                f"The '{phase}' block must be a mapping of options, "
                f"got {type(block).__name__}."
            )
        _validate_block(block, PHASE_SCHEMAS[phase], context=f"{phase}")

    simulation = data.get("simulation", {})
    if isinstance(simulation, dict):
        checkpoint_interval = simulation.get("checkpoint_interval_steps")
        if isinstance(checkpoint_interval, str):
            raw = checkpoint_interval.strip().lower()
            if raw != "auto":
                try:
                    int(raw)
                except ValueError as exc:
                    raise ConfigError(
                        "simulation option 'checkpoint_interval_steps' should be "
                        "an integer, 0, or 'auto'."
                    ) from exc
        elif isinstance(checkpoint_interval, int) and checkpoint_interval < 0:
            raise ConfigError(
                "simulation option 'checkpoint_interval_steps' must be >= 0 or 'auto'."
            )

        retries = simulation.get("first_aid_max_retries")
        if isinstance(retries, int) and retries < 0:
            raise ConfigError(
                "simulation option 'first_aid_max_retries' must be >= 0."
            )

        restraint_type = simulation.get("restraint_type")
        if restraint_type is not None and str(restraint_type).lower() not in {"none", "position"}:
            raise ConfigError(
                "simulation option 'restraint_type' must be 'position' or 'none'."
            )

        selection = simulation.get("restraint_selection")
        if selection is not None:
            allowed = {"protein", "backbone", "heavy", "non-water", "nonwater", "all"}
            if str(selection).strip().lower().replace("_", "-") not in allowed:
                raise ConfigError(
                    "simulation option 'restraint_selection' must be one of: "
                    "protein, backbone, heavy, non-water, all."
                )

        restraint_k = simulation.get("restraint_k")
        if isinstance(restraint_k, (int, float)) and restraint_k < 0:
            raise ConfigError(
                "simulation option 'restraint_k' must be >= 0."
            )

    # analysis include/exclude mutual exclusion (nested)
    analysis = data.get("analysis", {})
    if isinstance(analysis, dict) and analysis.get("include") and analysis.get("exclude"):
        raise ConfigError(
            "The 'analysis' block sets both 'include' and 'exclude'; "
            "they are mutually exclusive."
        )


# ---------------------------------------------------------------------------
# Phase-option extraction
# ---------------------------------------------------------------------------
def phase_options(data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Extract the per-phase option blocks from a validated config.

    Returns ``{phase: {option: value, ...}}`` with ``None`` values
    dropped so phase ``DEFAULTS`` apply. Used by :class:`BatchExplorer`
    to assemble the base options shared by every run.
    """
    options: dict[str, dict[str, Any]] = {}
    for phase in PHASE_KEYS:
        if phase in data and isinstance(data[phase], dict):
            block = {k: v for k, v in data[phase].items() if v is not None}
            if block:
                options[phase] = block
    return options
