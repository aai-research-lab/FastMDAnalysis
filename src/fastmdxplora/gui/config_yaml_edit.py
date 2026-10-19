"""Validate an edited Config and locate refusals in its original YAML text."""
from __future__ import annotations

from copy import deepcopy
import math
from typing import Any

import yaml
from yaml.nodes import MappingNode, ScalarNode, SequenceNode

from fastmdxplora.config.loader import ALIAS_GROUPS, ConfigError, _validate_block, canonical_phase_keys, validate_config
from fastmdxplora.config.schema import EXECUTION, PHASE_SCHEMAS, TOP_LEVEL
from fastmdxplora.gui.config_builder import (
    _lists_as_lists, build_config, export_config, phases_it_runs, refused, state_from_config,
)

MAX_TEXT = 1_000_000
MAX_NODES = 20_000
MAX_DEPTH = 60
MAX_WINDOWS = 1000


class _YamlError(ConfigError):
    default_code = "config.file.unparseable"

    def __init__(self, context, context_mark, problem, problem_mark):
        super().__init__(problem)
        self.problem_mark = problem_mark


class _StrictLoader(yaml.SafeLoader):
    def construct_object(self, node, deep=False):
        try:
            return super().construct_object(node, deep=deep)
        except _YamlError:
            raise
        except (ValueError, OverflowError) as exc:
            raise _YamlError(None, None, "The YAML scalar cannot be read: " + str(exc), node.start_mark) from exc

    def construct_mapping(self, node, deep=False):
        # Explicit duplicate keys are usually an editing mistake. YAML merge
        # defaults may still be overridden by an explicit key, as in the CLI.
        seen = set()
        for key, _ in node.value:
            name = "<<" if key.tag == "tag:yaml.org,2002:merge" else self.construct_object(key, deep=deep)
            try:
                if name in seen:
                    raise _YamlError(None, None,
                        "Duplicate YAML key: " + str(name), key.start_mark)
                seen.add(name)
            except TypeError as exc:
                raise _YamlError(None, None,
                    "Config mapping keys must be scalar values.", key.start_mark) from exc
        return super().construct_mapping(node, deep=deep)


def _bounded(text: str) -> None:
    stack = []
    for count, event in enumerate(yaml.parse(text), 1):
        if count > MAX_NODES:
            raise ConfigError("The YAML has too many nodes to edit here.", code="config.file.unparseable")
        if isinstance(event, (yaml.MappingStartEvent, yaml.SequenceStartEvent)):
            stack.append(event.anchor)
            if len(stack) > MAX_DEPTH:
                raise _YamlError(None, None,
                    "The YAML is nested too deeply.", event.start_mark)
        elif isinstance(event, (yaml.MappingEndEvent, yaml.SequenceEndEvent)):
            stack.pop()
        elif isinstance(event, yaml.AliasEvent) and event.anchor in stack:
            raise _YamlError(None, None,
                "Recursive YAML aliases cannot describe a Config.", event.start_mark)


def _locations(node, prefix=(), active=None, budget=None):
    active = set() if active is None else active
    budget = [MAX_NODES] if budget is None else budget
    budget[0] -= 1
    if budget[0] < 0:
        raise ConfigError("Expanded YAML aliases exceed the editor's node limit.", code="config.file.unparseable")
    if id(node) in active:
        raise ConfigError("Recursive YAML aliases cannot describe a Config.", code="config.file.unparseable")
    active.add(id(node))
    result = {}
    if isinstance(node, MappingNode):
        for key, value in node.value:
            if not isinstance(key, ScalarNode):
                continue
            path = prefix + (key.value,)
            result[path] = {"line": key.start_mark.line + 1, "column": key.start_mark.column + 1,
                            "path": ".".join(path)}
            result.update(_locations(value, path, active, budget))
    elif isinstance(node, SequenceNode):
        for i, value in enumerate(node.value):
            path = prefix + (str(i),)
            result[path] = {"line": value.start_mark.line + 1, "column": value.start_mark.column + 1,
                            "path": ".".join(path)}
            result.update(_locations(value, path, active, budget))
    active.remove(id(node))
    return result


def _guard_values(value, path=()):
    """Bound work that the shared validator's normalisation can generate."""
    if isinstance(value, float) and not math.isfinite(value):
        raise ConfigError("Numeric Config values must be finite.", code="config.option.out_of_range",
                          option=path[-1] if path else None, context=".".join(path[:-1]) or "top-level")
    if isinstance(value, dict):
        if path[-2:] == ("simulation", "umbrella"):
            for name in ("centres", "centers", "n_windows"):
                count = value.get(name)
                if count is None:
                    continue
                if name in {"centres", "centers"} and not isinstance(count, list):
                    raise ConfigError("Umbrella centres must be a YAML list.", code="config.option.wrong_type",
                                      option=name, context=".".join(path))
                try:
                    count = len(count) if name in {"centres", "centers"} else int(count)
                except (TypeError, ValueError, OverflowError) as exc:
                    raise ConfigError("The umbrella window count must be an integer.", code="config.option.wrong_type",
                                      option=name, context=".".join(path)) from exc
                if count > MAX_WINDOWS:
                    raise ConfigError("The YAML editor checks at most 1000 umbrella windows. Use the CLI for a larger study.",
                        code="config.option.out_of_range", option=name, context=".".join(path), maximum=MAX_WINDOWS)
            # The umbrella reader accepts numeric strings inside its block.
            # Check their converted values before it can create windows.
            for name in ("centres", "centers", "centre", "from", "to", "force_constant"):
                raw = value.get(name)
                for number in raw if isinstance(raw, list) else [raw]:
                    if number is None:
                        continue
                    try:
                        finite = math.isfinite(float(number))
                    except (TypeError, ValueError, OverflowError):
                        continue  # The reader reports malformed values.
                    if not finite:
                        raise ConfigError("Umbrella numeric values must be finite.", code="config.option.out_of_range",
                                          option=name, context=".".join(path))
        for key, item in value.items():
            if (not isinstance(key, (str, int, float, bool, type(None)))
                    or isinstance(key, float) and not math.isfinite(key)):
                raise ConfigError("Config mapping keys must be finite text or scalar values.",
                                  code="config.option.wrong_type")
            _guard_values(item, path + (str(key),))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            _guard_values(item, path + (str(i),))
    elif not isinstance(value, (str, int, float, bool, type(None))):
        raise ConfigError("Config values must use JSON-compatible YAML mappings, lists and scalar values.",
                          code="config.option.wrong_type", option=path[-1] if path else None,
                          context=".".join(path[:-1]) or "top-level")


def _diagnostic(exc, locations, prefix=()):
    answer = refused(exc)
    record = answer["refusal"]
    details = record.get("details") or {}
    context = details.get("context")
    context = () if context in (None, "top-level") else tuple(str(context).split("."))
    options = details.get("options") or [details.get("option")]
    found = []
    for option in options:
        if not isinstance(option, str):
            continue
        wanted = prefix + context + (option,)
        direct = locations.get(wanted)
        if direct is None:
            # Only the loader's declared aliases may identify another key.
            aliases = next((group for group in ALIAS_GROUPS if option in group), ())
            candidates = [locations[prefix + context + (name,)] for name in aliases
                          if prefix + context + (name,) in locations]
            if len(candidates) == 1:
                direct = candidates[0]
        if direct and direct not in found:
            found.append(direct)
    mark = getattr(exc, "problem_mark", None)
    if mark is not None:
        found = [{"line": mark.line + 1, "column": mark.column + 1, "path": ""}]
    return {"ok": False, "error": str(exc), **answer,
            "diagnostics": [{"message": str(exc), "code": record["code"],
                             "locations": found}]}


def _canonical(config):
    chosen = _lists_as_lists(deepcopy(config))
    # Default filling must see the same phase lists and alias spellings as
    # the CLI, including comma-joined legacy include/exclude keys.
    validate_config(chosen, require_systems=True)
    full = export_config(chosen, full=True)
    if not full.get("ok"):
        return None
    result = yaml.safe_load(full["yaml"])
    for field in TOP_LEVEL.fields:
        result.setdefault(field.name, deepcopy(field.default))
    execution = result.get("execution")
    if not isinstance(execution, dict):
        execution = result["execution"] = {}
    for field in EXECUTION.fields:
        execution.setdefault(field.name, deepcopy(field.default))
    result["include_phase"] = phases_it_runs(result)
    result.pop("exclude_phase", None)
    return result


def _lost(expected, actual, path=()):
    if isinstance(expected, dict) and isinstance(actual, dict):
        found = []
        for key, value in expected.items():
            if key not in actual:
                # An absent optional null has the same meaning in the loader.
                if value is not None:
                    found.append(path + (str(key),))
            else:
                found.extend(_lost(value, actual[key], path + (str(key),)))
        return found
    if isinstance(expected, list) and isinstance(actual, list) and len(expected) == len(actual):
        return [p for i, (a, b) in enumerate(zip(expected, actual))
                for p in _lost(a, b, path + (str(i),))]
    return [] if expected == actual else [path]


def validate_yaml(request: dict[str, Any]) -> dict[str, Any]:
    """No files are read or written and no execution is started."""
    text = request.get("text")
    try:
        size = len(text.encode("utf-8")) if isinstance(text, str) else MAX_TEXT + 1
    except UnicodeError:
        size = MAX_TEXT + 1
    if not isinstance(text, str) or size > MAX_TEXT:
        return _diagnostic(ConfigError("Config YAML must be text within 1 MB.", code="config.file.unparseable"), {})
    locations = {}
    try:
        _bounded(text)
        tree = yaml.compose(text, Loader=yaml.SafeLoader)
        locations = _locations(tree) if tree is not None else {}
        top = {path[0]: place for path, place in locations.items() if len(path) == 1}
        for name, place in canonical_phase_keys(top).items():
            locations.setdefault((name,), place)
        original = yaml.load(text, Loader=_StrictLoader)
        if not isinstance(original, dict):
            raise ConfigError("Config YAML must be a mapping.", code="config.file.not_a_mapping")
        _guard_values(original)
        config = deepcopy(original)
        validate_config(config, require_systems=True)
        # Per-system overrides are checked with the same schema before they
        # reach a phase, without fetching structures or expanding a sweep.
        for i, system in enumerate(original.get("systems") or []):
            if not isinstance(system, dict):
                continue
            try:
                effective = deepcopy(original)
                effective["systems"] = [{"system": system["system"]}]
                for phase in PHASE_SCHEMAS:
                    if isinstance(system.get(phase), dict):
                        effective[phase] = {**(original.get(phase) or {}), **system[phase]}
                simulation = effective.get("simulation") or {}
                window = simulation.get("umbrella")
                if isinstance(window, dict) and "centre" in window and not any(
                        name in window for name in ("centres", "centers", "n_windows", "from", "to")):
                    # A single window has already been expanded. Check its
                    # declared fields, then validate the rest of the effective
                    # Config without interpreting that window as a new plan.
                    from fastmdxplora.simulation.umbrella import check_umbrella_keys

                    _validate_block(simulation, PHASE_SCHEMAS["simulation"], context="simulation")
                    check_umbrella_keys(window)
                    effective["simulation"] = {key: value for key, value in simulation.items() if key != "umbrella"}
                _guard_values(effective)
                validate_config(effective, require_systems=True)
            except ConfigError as exc:
                return _diagnostic(exc, locations, ("systems", str(i)))
        for_form = canonical_phase_keys(deepcopy(original))
        for name in ("include_phase", "exclude_phase"):
            if name in config:
                for_form[name] = deepcopy(config[name])
        loaded = state_from_config(for_form)
        if not loaded.get("ok"):
            return {**loaded, "diagnostics": [{"message": loaded.get("error"), "code": None, "locations": []}]}
        _guard_values(loaded["state"])
        answer = {"ok": True, "state": loaded["state"], "diagnostics": [], "apply_ok": False}
        if isinstance(request.get("form"), dict):
            _guard_values(request["form"])
            proposed = build_config(request["form"], full=True)
            _guard_values(proposed)
            validate_config(proposed, require_systems=True)
            expected, actual = _canonical(original), _canonical(proposed)
            if expected is None or actual is None:
                return {**answer, "mapping_reason": "The Builder mapping could not be checked. Keep the YAML before applying."}
            losses = list(dict.fromkeys(_lost(expected, actual) + _lost(actual, expected)))
            answer["apply_ok"] = not losses
            if losses:
                answer["mapping_reason"] = "The Builder cannot retain these YAML settings: " + ", ".join(".".join(p) for p in losses[:15]) + ". Keep the YAML or revise it before applying."
                answer["diagnostics"] = [{"message": answer["mapping_reason"], "code": None,
                    "locations": [locations[p] for p in losses if p in locations][:15]}]
        return answer
    except yaml.YAMLError as exc:
        mark = getattr(exc, "problem_mark", None)
        message = getattr(exc, "problem", None) or "The YAML could not be parsed."
        result = _diagnostic(ConfigError(message, code="config.file.unparseable"), {})
        if mark is not None:
            result["diagnostics"][0]["locations"] = [{"line": mark.line + 1, "column": mark.column + 1, "path": ""}]
        return result
    except ConfigError as exc:
        return _diagnostic(exc, locations)
    except (TypeError, ValueError, OverflowError, KeyError, AttributeError, RecursionError):
        return _diagnostic(ConfigError(
            "The YAML contains a key or value the Config schema cannot interpret. Use text keys and the declared value types.",
            code="config.option.wrong_type"), {})
