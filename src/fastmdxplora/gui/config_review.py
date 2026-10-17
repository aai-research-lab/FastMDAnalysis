"""Choose changes between two Configs without changing either version.

The browser sends choices by their place in the server-derived list. Display
strings and dotted setting names are never used to edit a Config. A compound
setting (systems, sweep, an umbrella or stopping rule) is one choice so its
internal structure and its decision record travel together.
"""

from __future__ import annotations

from copy import deepcopy
import json
from typing import Any

from fastmdxplora.config.loader import (
    ConfigError,
    _ONE_THING_TWO_WORDS,
    _settle_the_words_that_mean_one_thing,
    canonical_phase_keys,
    validate_config,
)
from fastmdxplora.config.schema import PHASE_KEYS
from fastmdxplora.gui.config_builder import refused
from fastmdxplora.gui.plan import changes_between, plan_of, sourced
from fastmdxplora.gui.schema_payload import label_and_unit

_MISSING = object()
_METADATA = {"agent", "agent_model", "decisions"}
_BLOCKS = {*PHASE_KEYS, "execution"}


def _fields(before: dict, after: dict) -> list[tuple[tuple[str, ...], ...]]:
    paths = []
    for key in dict.fromkeys([*after, *before]):
        if key in _METADATA:
            continue
        old, new = before.get(key), after.get(key)
        if key in _BLOCKS and all(v is None or isinstance(v, dict) for v in (old, new)):
            paths.extend((key, name) for name in dict.fromkeys([*(new or {}), *(old or {})]))
        else:
            paths.append((key,))
    for key in dict.fromkeys([*(after.get("decisions") or {}), *(before.get("decisions") or {})]):
        if key in _METADATA:
            continue
        parts = key.split(".", 1)
        path = tuple(parts) if parts[0] in _BLOCKS else (key,)
        if path not in paths:
            paths.append(path)
    groups = []
    for path in paths:
        group = (path,)
        if len(path) == 2:
            for established, general in _ONE_THING_TWO_WORDS.get(path[0], ()):
                if path[1] in (established, general):
                    group = ((path[0], established), (path[0], general))
                    break
        if group not in groups:
            groups.append(group)
    return groups


def _value(config: dict, path: tuple[str, ...]) -> Any:
    value = config
    for part in path:
        if not isinstance(value, dict) or part not in value:
            return _MISSING
        value = value[part]
    return value


def _said(value: Any) -> str:
    if value is _MISSING:
        return "not set"
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _changes(before: dict, after: dict) -> list[tuple[tuple[tuple[str, ...], ...], dict]]:
    displayed = {row["setting"]: row for row in changes_between(before, after)}
    old_decisions, new_decisions = before.get("decisions") or {}, after.get("decisions") or {}
    rows = []
    for paths in _fields(before, after):
        path = paths[0]
        setting = ".".join(path)
        old, new = _value(before, path), _value(after, path)
        names = [".".join(p) for p in paths]
        reason_changed = any(old_decisions.get(n) != new_decisions.get(n) for n in names)
        if old == new and not reason_changed:
            continue
        row = dict(displayed.get(setting) or {
            "setting": setting, "label": label_and_unit(path[-1])[0],
            "before": _said(old), "after": _said(new),
        })
        row["id"] = str(len(rows))
        if len(paths) > 1:
            row["aliases"] = names
        if reason_changed:
            row["before_decision"] = deepcopy(old_decisions.get(setting))
            row["after_decision"] = deepcopy(new_decisions.get(setting))
            row["before_decisions"] = {n: deepcopy(old_decisions[n]) for n in names if n in old_decisions}
            row["after_decisions"] = {n: deepcopy(new_decisions[n]) for n in names if n in new_decisions}
        rows.append((paths, row))
    return rows


def review_config(payload: dict[str, Any]) -> dict[str, Any]:
    """Derive choices, or validate the selected combination. No file IO.

    With no ``accepted`` list this only lists the changes. With a list,
    every other change is rejected. The validator works on copies because
    it normalises aliases and expands umbrella windows in place.
    """
    if not isinstance(payload, dict):
        return {"ok": False, "error": "Review needs the previous and proposed Configs."}
    before, after = payload.get("before"), payload.get("after")
    if not isinstance(before, dict) or not isinstance(after, dict):
        return {"ok": False, "error": "Review needs the previous and proposed Configs."}
    try:
        validate_config(deepcopy(before), require_systems=True)
        validate_config(deepcopy(after), require_systems=True)
    except ConfigError as exc:
        return {"ok": False, "error": str(exc), **refused(exc)}
    # Settle spellings without expanding an umbrella into systems. The
    # validator's full normalisation still runs on separate copies above.
    before, after = deepcopy(before), deepcopy(after)
    for config in (before, after):
        canonical_phase_keys(config)
        _settle_the_words_that_mean_one_thing(config)
    changes = _changes(before, after)
    rows = [row for _, row in changes]
    if "accepted" not in payload:
        return {"ok": True, "changes": rows}
    accepted = payload["accepted"]
    ids = {row["id"] for row in rows}
    if (not isinstance(accepted, list) or any(not isinstance(i, str) for i in accepted)
            or len(set(accepted)) != len(accepted) or not set(accepted) <= ids):
        return {"ok": False, "error": "Choose only changes in this Config review."}
    chosen = set(accepted)
    config = deepcopy(before)
    decisions = deepcopy(before.get("decisions") or {})
    for paths, row in changes:
        if row["id"] not in chosen:
            continue
        for path in paths:
            value = _value(after, path)
            if value is not _MISSING or _value(config, path) is not _MISSING:
                parent = config
                for part in path[:-1]:
                    if not isinstance(parent.get(part), dict):
                        parent[part] = {}
                    parent = parent[part]
                if value is _MISSING:
                    parent.pop(path[-1], None)
                    if len(path) == 2 and not parent and path[0] not in after:
                        config.pop(path[0], None)
                else:
                    parent[path[-1]] = deepcopy(value)
            setting = ".".join(path)
            decision = (after.get("decisions") or {}).get(setting, _MISSING)
            if decision is _MISSING:
                decisions.pop(setting, None)
            else:
                decisions[setting] = deepcopy(decision)
    for name in _METADATA - {"decisions"}:
        if name in after:
            config[name] = deepcopy(after[name])
        else:
            config.pop(name, None)
        if name in (after.get("decisions") or {}):
            decisions[name] = deepcopy(after["decisions"][name])
        else:
            decisions.pop(name, None)
    if decisions:
        config["decisions"] = decisions
    else:
        config.pop("decisions", None)
    try:
        validate_config(deepcopy(config), require_systems=True)
    except ConfigError as exc:
        return {"ok": False, "error": str(exc), "changes": rows, **refused(exc)}
    import yaml

    return {
        "ok": True, "config": config,
        "yaml": yaml.safe_dump(config, sort_keys=False),
        "plan": sourced(plan_of(config), config),
        "changes": changes_between(after, config),
        "review": {"accepted": [r["setting"] for r in rows if r["id"] in chosen],
                   "rejected": [r["setting"] for r in rows if r["id"] not in chosen]},
    }
