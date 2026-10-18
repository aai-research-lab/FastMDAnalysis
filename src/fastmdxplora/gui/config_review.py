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
from fastmdxplora.gui.plan import _said as _said_as_shown, changes_between, plan_of, sourced
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


def _said(value: Any, unit: str = "") -> str:
    if value is _MISSING:
        return "not set"
    return _said_as_shown(value, unit)


#: Where a decision came from, as the plan's tags name it (agent-panel.js).
_SOURCE_SAID = {"person": "you asked", "fastmdx-defaults.yml": "your default",
                "agent": "Agent\u2019s choice"}


def _reason(decision: Any, unit: str) -> str:
    """A decision record as a sentence: why, where from, the alternatives."""
    if not isinstance(decision, dict):
        return "No reason recorded"
    why = str(decision.get("why") or "").strip() or "No reason recorded"
    notes = []
    source = decision.get("source")
    if source:
        notes.append(_SOURCE_SAID.get(str(source), str(source)))
    if "alternatives" in decision:
        others = decision["alternatives"]
        said = (", ".join(_said(v, unit) for v in others) or "none"
                if isinstance(others, list) else _said(others, unit))
        notes.append(f"alternatives: {said}")
    return f"{why} ({'; '.join(notes)})" if notes else why


def _shown(displayed: list[dict], names: list[str], path: tuple[str, ...],
           old: Any, new: Any) -> dict[str, Any]:
    """A choice as the version's change list says it, never as raw JSON.

    One line of `changes_between` is used as it is; a compound setting
    said in several lines keeps each as a part. Mirrored names are one
    setting, so only the first name with a line is read.
    """
    label, unit = label_and_unit(path[-1])
    for name in names:
        parts = [row for row in displayed
                 if row["setting"] == name or row["setting"].startswith(name + ".")]
        if parts:
            break
    if len(parts) == 1:
        return {"label": parts[0]["label"], "before": parts[0]["before"],
                "after": parts[0]["after"]}
    if parts:
        return {"label": label,
                "before": "; ".join(f"{p['label']} {p['before']}" for p in parts),
                "after": "; ".join(f"{p['label']} {p['after']}" for p in parts),
                "parts": [{k: p[k] for k in ("label", "before", "after")} for p in parts]}
    before, after = _said(old, unit), _said(new, unit)
    if old != new and before == after:
        before = json.dumps(old, ensure_ascii=False, sort_keys=True, default=str)
        after = json.dumps(new, ensure_ascii=False, sort_keys=True, default=str)
    shown = {"label": label, "before": before, "after": after}
    if old == new:
        shown["same_value"] = True
    return shown


def _changes(before: dict, after: dict) -> list[tuple[tuple[tuple[str, ...], ...], dict]]:
    displayed = changes_between(before, after)
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
        row = {"setting": setting, **_shown(displayed, names, path, old, new)}
        row["id"] = str(len(rows))
        if len(paths) > 1:
            row["aliases"] = names
            row["alias_labels"] = [label_and_unit(p[-1])[0] for p in paths]
        if reason_changed:
            unit = label_and_unit(path[-1])[1]
            row["before_reason"] = _reasons(old_decisions, names, unit)
            row["after_reason"] = _reasons(new_decisions, names, unit)
            row["before_decision"] = deepcopy(old_decisions.get(setting))
            row["after_decision"] = deepcopy(new_decisions.get(setting))
            row["before_decisions"] = {n: deepcopy(old_decisions[n]) for n in names if n in old_decisions}
            row["after_decisions"] = {n: deepcopy(new_decisions[n]) for n in names if n in new_decisions}
        rows.append((paths, row))
    return rows


def _reasons(decisions: dict, names: list[str], unit: str) -> str:
    """The reasons recorded under any of a setting's names, each said once."""
    said = list(dict.fromkeys(_reason(decisions[n], unit) for n in names if n in decisions))
    return "; ".join(said) or "No reason recorded"


def _settled(config: dict) -> dict:
    """A copy with phase lists and mirrored names settled, umbrellas kept."""
    config = canonical_phase_keys(deepcopy(config))
    _settle_the_words_that_mean_one_thing(config)
    return config


def _one_name_each(changes: list[dict]) -> list[dict]:
    """Changes between settled Configs, a mirrored setting under one name."""
    established = {f"{phase}.{g}": f"{phase}.{e}"
                   for phase, pairs in _ONE_THING_TWO_WORDS.items() for e, g in pairs}
    said = {row["setting"] for row in changes}
    return [row for row in changes if established.get(row["setting"]) not in said]


def _receipt(rows: list[dict], chosen: set[str]) -> str:
    """The choices as the person saw them: label and old to new value."""
    def said(row: dict) -> str:
        if row.get("same_value"):
            return f"{row['label']} (its reason)"
        return f"{row['label']} {row['before']} \u2192 {row['after']}"

    accepted = [said(r) for r in rows if r["id"] in chosen]
    rejected = [said(r) for r in rows if r["id"] not in chosen]
    return " ".join(f"{word}: {'; '.join(items)}." for word, items
                    in (("Accepted", accepted), ("Rejected", rejected)) if items)


def _involved(rows: list[dict], refusal: dict) -> list[str]:
    """The rows a refusal's remedy names, by the settings it lists."""
    named = [str(s) for s in (refusal.get("remedy") or {}).get("settings") or ()]

    def touches(row: dict) -> bool:
        return any(name == s or s.startswith(name + ".") or name.startswith(s + ".")
                   for name in row.get("aliases") or [row["setting"]] for s in named)

    return [row["id"] for row in rows if touches(row)]


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
    # Compared with spellings settled, without expanding an umbrella into
    # systems; the validator's full normalisation runs on copies above.
    # Values are taken from each version as written, so a reviewed version
    # holds the spelling its source used and nothing the review did not
    # show (round 1 review: an added `setup_from` listed as a change).
    written = {"before": canonical_phase_keys(deepcopy(before)),
               "after": canonical_phase_keys(deepcopy(after))}
    before, after = _settled(before), _settled(after)
    changes = _changes(before, after)
    for paths, row in changes:
        # Two names said only where a version was written with both.
        used = {p for p in paths for c in written.values() if _value(c, p) is not _MISSING}
        if len(used) < 2:
            row.pop("alias_labels", None)
    rows = [row for _, row in changes]
    if "accepted" not in payload:
        return {"ok": True, "changes": rows}
    accepted = payload["accepted"]
    ids = {row["id"] for row in rows}
    if (not isinstance(accepted, list) or any(not isinstance(i, str) for i in accepted)
            or len(set(accepted)) != len(accepted) or not set(accepted) <= ids):
        return {"ok": False, "error": "Choose only changes in this Config review."}
    chosen = set(accepted)
    config = deepcopy(written["before"])
    decisions = deepcopy(before.get("decisions") or {})
    for paths, row in changes:
        if row["id"] not in chosen:
            continue
        for path in paths:
            value = _value(written["after"], path)
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
        found = refused(exc)
        return {"ok": False, "error": str(exc), "changes": rows, **found,
                "involved": _involved(rows, found)}
    import yaml

    return {
        "ok": True, "config": config,
        "yaml": yaml.safe_dump(config, sort_keys=False),
        "plan": sourced(plan_of(config), config),
        "changes": _one_name_each(changes_between(after, _settled(config))),
        "review": {"accepted": [r["setting"] for r in rows if r["id"] in chosen],
                   "rejected": [r["setting"] for r in rows if r["id"] not in chosen],
                   "said": _receipt(rows, chosen)},
    }
