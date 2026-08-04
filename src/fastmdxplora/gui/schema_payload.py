"""The settings, described well enough for a browser to draw them.

The dashboard's form was written by hand, one control at a time, which is why
it offered eleven of the eighty-three settings that exist: every new field had
to be noticed and added, and a field nobody noticed was simply unreachable from
the browser. Whole phases were missing -- there was no way to configure an
analysis or a report at all.

So the form is not written any more. This turns the schema into something a
browser can render: a name, a type, what it means, what values it accepts, and
what it does if left alone. Adding a field to the schema puts a control in the
dashboard, and nothing has to be kept in step by hand.

The per-analysis settings are read the same way, from the analyses themselves,
via ``analysis.describe``. Those need the chemistry stack, so they are asked
for separately and their absence is reported rather than raised: a dashboard
that cannot offer clustering options is worth more than one that will not load.
"""

from __future__ import annotations

from typing import Any

from fastmdxplora.config.schema import PHASE_SCHEMAS

__all__ = ["schema_payload", "field_payload"]


#: How a value is best asked for. The schema says what a field *is*; this says
#: what a browser should put on the screen for it. Kept here rather than in the
#: schema because it is a fact about drawing, not about configuring.
_CONTROL_FOR_TYPE = {
    bool: "checkbox",
    int: "number",
    float: "number",
    str: "text",
    list: "list",
    dict: "mapping",
}


def _control(field: Any) -> str:
    if field.choices:
        return "select"
    return _CONTROL_FOR_TYPE.get(field.type, "text")


def _jsonable(value: Any) -> Any:
    """Values reach the browser as JSON, and tuples are not JSON."""
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, dict)):
        return value
    return str(value)


def field_payload(field: Any) -> dict[str, Any]:
    """One setting, as much as a browser needs to offer it."""
    return {
        "name": field.name,
        "type": getattr(field.type, "__name__", str(field.type)),
        "control": _control(field),
        "help": field.help,
        "default": _jsonable(field.default),
        "choices": _jsonable(field.choices),
        "example": _jsonable(field.example),
        "required": bool(getattr(field, "required", False)),
    }


def _analysis_options() -> dict[str, Any]:
    """What each analysis can be told, read from the analyses themselves.

    Importing them needs MDTraj. Where it is absent the dashboard should still
    load and still offer everything else, so the failure is described rather
    than raised.
    """
    try:
        import fastmdxplora.analysis  # noqa: F401  (populates the registry)
        from fastmdxplora.analysis.describe import (
            describe_all, explain_analysis,
        )
    except Exception as exc:  # noqa: BLE001
        return {
            "available": False,
            "reason": (
                "Per-analysis settings need the analysis stack, which is not "
                f"installed here ({type(exc).__name__}). Everything else on "
                "this page still works."
            ),
            "analyses": {},
            "explanations": {},
        }

    analyses: dict[str, Any] = {}
    explanations: dict[str, Any] = {}
    for name, options in describe_all().items():
        explanations[name] = explain_analysis(name)
        analyses[name] = [
            {
                "name": option.name,
                "control": (
                    "select" if option.choices
                    else _CONTROL_FOR_TYPE.get(type(option.default), "text")
                ),
                "help": option.help,
                "default": _jsonable(option.default),
                "choices": _jsonable(option.choices),
                "type_text": option.type_text,
                # Settings every analysis shares come from the base class and
                # are worth grouping apart from the ones that make an analysis
                # what it is.
                "shared": option.owner == "Analysis",
            }
            for option in options
        ]
    return {
        "available": True,
        "reason": None,
        "analyses": analyses,
        "explanations": explanations,
    }


def schema_payload() -> dict[str, Any]:
    """Every setting the software accepts, ready to be drawn."""
    phases = {
        phase: {
            "fields": [field_payload(f) for f in group.fields],
            "description": getattr(group, "description", None),
        }
        for phase, group in PHASE_SCHEMAS.items()
    }
    return {
        "phases": phases,
        "analysis_options": _analysis_options(),
    }
