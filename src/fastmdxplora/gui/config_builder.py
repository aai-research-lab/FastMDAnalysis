"""Turn what the page holds into a config file, and say so if it will not work.

The point of writing a config rather than starting a run directly is that the
same file works somewhere else. A laptop is a poor place to run fifty
nanoseconds; a cluster is a poor place to decide what to run. So the browser
builds the file, and the file goes wherever the compute is::

    fastmdx explore --config exploration.yml

Nothing is generated here that the command line would reject. The config is put
through the same validator the CLI uses before it is handed back, so a
configuration that will fail on the cluster fails in the browser instead, while
the person is still looking at the form that produced it.

Values equal to their default are left out. A file that restates forty settings
the software would have chosen anyway hides the three the person actually
decided, and those three are the ones a reader of the methods section needs.
"""

from __future__ import annotations

from typing import Any

from fastmdxplora.config.loader import ConfigError, validate_config
from fastmdxplora.config.schema import PHASE_SCHEMAS

__all__ = ["build_config", "config_yaml"]


def _defaults_for(phase: str) -> dict[str, Any]:
    return {f.name: f.default for f in PHASE_SCHEMAS[phase].fields}


def _coerce(value: Any, field: Any) -> Any:
    """Bring a value in from a form, where everything arrives as text.

    A field's type may be several types at once -- a temperature is written
    either 300 or 300.0 -- so the accepted types are taken as a set rather
    than compared one at a time.
    """
    if value is None or value == "":
        return None
    accepted = field.type if isinstance(field.type, tuple) else (field.type,)

    if bool in accepted:
        if isinstance(value, bool):
            return value
        return str(value).strip().lower() in {"true", "1", "yes", "on"}
    if list in accepted and isinstance(value, str):
        return [part for part in (p.strip() for p in value.split(",")) if part]
    if float in accepted or int in accepted:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return value
        # A whole number goes in as one where the field takes an int, so the
        # file reads 300 rather than 300.0.
        if int in accepted and number.is_integer():
            return int(number)
        return number
    return value


def build_config(state: dict[str, Any]) -> dict[str, Any]:
    """A config dict from the page's state, with the defaults left out.

    ``state`` holds a block per phase, plus the top-level keys the run needs
    (``system``, ``output``, ``include``). Unknown phases and unknown fields
    are dropped rather than passed through: the page should not be able to
    write a file the command line will refuse to read.
    """
    config: dict[str, Any] = {}

    for key in ("output", "include", "exclude", "verbose"):
        value = state.get(key)
        if value not in (None, "", [], {}):
            config[key] = value

    # A run names its inputs under `systems`, one entry each, so that a
    # config describing one protein and one describing forty have the same
    # shape. The page offers a single system; the file it writes is the same
    # file a batch would use, with one entry in the list.
    system = state.get("system")
    if system:
        entry: dict[str, Any] = {"system": system}
        if state.get("system_id"):
            entry["id"] = state["system_id"]
        config["systems"] = [entry]
    elif state.get("systems"):
        config["systems"] = state["systems"]

    for phase, group in PHASE_SCHEMAS.items():
        block = state.get(phase)
        if not isinstance(block, dict):
            continue
        defaults = _defaults_for(phase)
        fields = {f.name: f for f in group.fields}
        kept: dict[str, Any] = {}
        for name, raw in block.items():
            field = fields.get(name)
            if field is None:
                continue
            value = _coerce(raw, field)
            if value is None:
                continue
            # Restating a default tells a reader nothing and buries the
            # settings that were actually chosen.
            if value == defaults.get(name):
                continue
            kept[name] = value
        if kept:
            config[phase] = kept

    return config


def config_yaml(state: dict[str, Any]) -> dict[str, Any]:
    """The config as text, validated, with what it will do stated plainly.

    Returns ``ok`` false and the validator's own message where the settings
    would be refused, so the failure arrives here rather than on the cluster an
    hour later.
    """
    import yaml

    config = build_config(state)

    try:
        validate_config(config)
    except ConfigError as exc:
        return {
            "ok": False,
            "error": str(exc),
            "yaml": None,
            "settings_changed": sum(
                len(v) for k, v in config.items() if isinstance(v, dict)
            ),
        }

    header = (
        "# Written by the FastMDXplora dashboard.\n"
        "#\n"
        "# Only settings that differ from the defaults appear here, so what\n"
        "# this file says is what was decided. Run it anywhere:\n"
        "#\n"
        "#     fastmdx explore --config this-file.yml\n"
        "#\n"
    )
    body = yaml.safe_dump(config, sort_keys=False, default_flow_style=False)
    return {
        "ok": True,
        "error": None,
        "yaml": header + body,
        "settings_changed": sum(
            len(v) for k, v in config.items() if isinstance(v, dict)
        ),
    }
