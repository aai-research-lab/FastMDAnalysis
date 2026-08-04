"""What each analysis can be told, read from where it is already written.

An analysis declares its settings twice over: the constructor names them and
gives their defaults, and the class docstring says what each one means, what
type it takes, and sometimes which values are allowed. Both are already there,
already accurate, and already what a reader of the source consults.

So nothing is declared again here. This module reads those two, and everything
that needs to describe an analysis -- the config template, the browser's
controls, a help message -- asks it rather than keeping a list of its own. A
list kept elsewhere would be a second thing to update when an analysis gains an
option, and the update that gets forgotten is the one nobody is looking at.

The cost is that describing an analysis means importing it, and importing it
means MDTraj. Callers that may run without the analysis stack import this
lazily and degrade, the same way a phase does when its backend is absent.
"""

from __future__ import annotations

import inspect
import re
from dataclasses import dataclass
from typing import Any

__all__ = ["OptionDoc", "describe_analysis", "describe_all", "undocumented_options"]


#: A numpy-style parameter entry: the name, then a type after the colon, then
#: an indented body running until a line that is not indented. One entry may
#: introduce several names at once -- ``title, xlabel, ylabel : str`` -- so the
#: name is matched anywhere in the comma-separated list that opens the line.
_PARAM = (
    r"^(?P<names>{name}|[\w, ]*\b{name}\b[\w, ]*)\s*:\s*"
    r"(?P<type>[^\n]*)\n(?P<body>(?:[ \t]+\S.*\n?)+)"
)

#: A type written as a set of literals -- ``{"ward", "complete"}`` -- states
#: the accepted values. Only the braced part is read: the rest of the line
#: carries the default, which is quoted too, and taking the whole line offered
#: the default twice.
_LITERAL_SET = re.compile(r"^\{(?P<values>[^}]*)\}")


@dataclass(frozen=True)
class OptionDoc:
    """One setting an analysis accepts."""

    name: str
    default: Any
    type_text: str | None
    choices: tuple[str, ...] | None
    help: str | None
    #: The class that introduced it. Settings every analysis shares come from
    #: the base and are worth grouping apart from the ones that make an
    #: analysis what it is.
    owner: str = ""

    @property
    def documented(self) -> bool:
        return bool(self.help)


def _clean(text: str) -> str:
    """Collapse a docstring paragraph to one line, without its markup."""
    joined = " ".join(text.split())
    return joined.replace("``", "")


def _choices_from(type_text: str | None) -> tuple[str, ...] | None:
    if not type_text:
        return None
    match = _LITERAL_SET.match(type_text.strip())
    if match is None:
        return None
    values = re.findall(r'"([^"]+)"', match.group("values"))
    # Order is kept, duplicates dropped: a docstring may list a value twice.
    seen: list[str] = []
    for value in values:
        if value not in seen:
            seen.append(value)
    return tuple(seen) or None


def _docstring_for(cls: type) -> str:
    """The class docstring, which is where parameters are written.

    ``inspect.getdoc(cls.__init__)`` is not it: these constructors carry no
    docstring of their own, and asking returns ``object.__init__``'s inherited
    one, which describes nothing and matches nothing.
    """
    parts = []
    for klass in cls.__mro__:
        if klass is object:
            continue
        if klass.__doc__:
            parts.append(inspect.cleandoc(klass.__doc__))
        # The base class documents its own parameters on __init__ rather than
        # on the class, so both are read. Taken from __dict__ so an inherited
        # docstring is not counted twice.
        init = klass.__dict__.get("__init__")
        if init is not None and init.__doc__:
            parts.append(inspect.cleandoc(init.__doc__))
    return "\n\n".join(parts)


def describe_analysis(name: str) -> tuple[OptionDoc, ...]:
    """Every setting ``name`` accepts, with what its own source says about it."""
    from fastmdxplora.analysis.orchestrator import get_analysis_class

    cls = get_analysis_class(name)
    doc = _docstring_for(cls)

    options: list[OptionDoc] = []
    seen: set[str] = set()
    for klass in cls.__mro__:
        init = klass.__dict__.get("__init__")
        if init is None:
            continue
        for param in inspect.signature(init).parameters.values():
            if param.kind is not inspect.Parameter.KEYWORD_ONLY:
                continue
            if param.name in seen:
                continue
            seen.add(param.name)
            match = re.search(
                _PARAM.format(name=re.escape(param.name)), doc, re.M
            )
            type_text = match.group("type").strip() if match else None
            options.append(
                OptionDoc(
                    name=param.name,
                    default=(
                        None
                        if param.default is inspect.Parameter.empty
                        else param.default
                    ),
                    type_text=type_text,
                    choices=_choices_from(type_text),
                    help=_clean(match.group("body")) if match else None,
                    owner=klass.__name__,
                )
            )
    return tuple(sorted(options, key=lambda o: o.name))


def describe_all() -> dict[str, tuple[OptionDoc, ...]]:
    """Every analysis and every setting it accepts."""
    from fastmdxplora.analysis.orchestrator import available_analyses

    return {name: describe_analysis(name) for name in available_analyses()}


def undocumented_options() -> dict[str, tuple[str, ...]]:
    """Settings an analysis accepts but its docstring does not explain.

    A setting nobody has described cannot be offered in a form or a template
    with any useful label, so this is what stands between the browser and a
    complete one. Reported rather than raised: a missing sentence should fail a
    check, not a run.
    """
    missing = {}
    for name, options in describe_all().items():
        absent = tuple(o.name for o in options if not o.documented)
        if absent:
            missing[name] = absent
    return missing
