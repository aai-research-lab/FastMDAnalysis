"""Your defaults: the values a workspace's studies start from.

FastMDXplora has its own defaults (300 K, 2 fs, ...), written in the
schema. A lab has its usual choices too: its assays run at 310 K, it
simulates with one force field, it runs three replicas. Written once in
``fastmdx-defaults.yml`` in the folder the studies are kept in, they fill
what a new study leaves unset, from every interface alike: the command line,
the GUI's Run and Config Builder, the Agent and an AI app.

The rules, each for a reason:

- **It fills only what a study leaves unset.** A value the config gives,
  even `null` (FastMDXplora's default, chosen), is the study's. A setting
  absent from the config is the only one filled.
- **It never decides what is simulated.** `systems`, `output`, `sweep` and
  the study-level keys are refused in it: the file holds how a lab
  simulates, not what.
- **A continuation is left alone.** A study that continues another
  (`simulation.resume_from`) keeps that study's settings, and one prepared
  from another (`simulation.setup_from`) keeps its setup.
- **Each value it fills is said.** It is recorded in the study's
  `decisions` with `source: fastmdx-defaults.yml` and the file's `why`, so
  the report's **Why these settings** says which values were the lab's,
  and the resolved config records the value itself, so the study runs the
  same again without the file.

The file is found as git finds its repository: the nearest
``fastmdx-defaults.yml`` in the study's folder or a folder above it,
stopping at the home folder.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["DEFAULTS_FILE", "YourDefaults", "defaults_for", "find_defaults",
           "read_defaults", "refused_with", "with_defaults", "write_defaults"]

DEFAULTS_FILE = "fastmdx-defaults.yml"

#: What the file may set: how a study is simulated, never what.
_BLOCKS = ("setup", "simulation", "analysis", "report", "execution")

def _one_of() -> tuple[frozenset[str], ...]:
    """Settings that answer one question between them: a config that gives
    one is not given another from the file. `analysis.include` filled beside
    the config's `analysis.exclude` was refused as setting both, and the
    file's `pressure_bar` beside a study's `pressure_atm` won, as bar does,
    so the study ran at the lab's pressure and not its own; the same for a
    stage's steps beside its duration."""
    from fastmdxplora.config.loader import SAME_SETTING

    return (frozenset({"analysis.include", "analysis.exclude"}), *SAME_SETTING)


def _answers(name: str, given: set[str]) -> bool:
    """Whether settings ``given`` (dotted names) already set ``name``."""
    if name in given:
        return True
    return any(name in group and group & given for group in _one_of())


@dataclass(frozen=True)
class YourDefaults:
    """A defaults file, read and checked."""

    path: Path
    #: The settings it gives, block by block, as a config holds them.
    values: dict[str, dict[str, Any]] = field(default_factory=dict)
    #: Why each is the usual value, by its dotted name, where the file says.
    whys: dict[str, str] = field(default_factory=dict)

    def settings(self) -> list[tuple[str, Any]]:
        """Every setting it gives, as (dotted name, value), in file order."""
        return [(f"{block}.{key}", value) for block, given in self.values.items()
                for key, value in given.items()]

    def said(self) -> str:
        """For the Agent and a person: each setting, its value and why."""
        import json

        lines = []
        for name, value in self.settings():
            why = self.whys.get(name)
            shown = value if isinstance(value, str) else json.dumps(value)
            lines.append(f"- `{name}`: {shown}" + (f" ({why})" if why else ""))
        return "\n".join(lines)


def find_defaults(start: str | Path, *, home: str | Path | None = None) -> Path | None:
    """The nearest defaults file in ``start`` or a folder above it, stopping
    at the home folder (or the top of the file system, for a folder outside
    it); None where there is none. ``start`` need not exist yet: a study's
    folder is named before it is made."""
    here = Path(start).expanduser()
    if not here.is_absolute():
        here = Path.cwd() / here
    stop = Path(home).expanduser() if home is not None else Path.home()
    try:
        stop = stop.resolve()
    except OSError:  # pragma: no cover - a home folder that cannot be resolved
        pass
    for folder in [here, *here.parents]:
        try:
            resolved = folder.resolve()
        except OSError:  # pragma: no cover - an unreadable parent
            continue
        candidate = resolved / DEFAULTS_FILE
        if candidate.is_file():
            return candidate
        if resolved == stop:
            return None
    return None


def read_defaults(path: str | Path) -> YourDefaults:
    """The file at ``path``, checked: a mapping of the blocks it may set,
    each setting one the schema knows with a value it accepts, and a `why`
    only for a setting it gives. Raises ConfigError saying what is wrong
    and in which file."""
    import yaml

    from fastmdxplora.config.loader import (
        ConfigError, canonical_phase_keys, validate_config)

    where = Path(path)
    try:
        data = yaml.safe_load(where.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ConfigError(f"{where} could not be read: {exc}",
                          code="config.file.missing", path=str(where)) from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"{where} is not YAML that can be read: {exc}",
                          code="config.file.unparseable", path=str(where)) from exc
    if data is None:
        return YourDefaults(where)
    if not isinstance(data, dict):
        raise ConfigError(
            f"{where} holds your defaults as a mapping of blocks (`simulation:`, "
            f"`setup:`, ...); it holds a {type(data).__name__}.",
            code="config.file.not_a_mapping", path=str(where),
            found_type=type(data).__name__)
    canonical_phase_keys(data)
    decisions = data.pop("decisions", None) or {}
    for key in data:
        if key not in _BLOCKS:
            raise ConfigError(
                f"{where} sets `{key}`. Your defaults say how a study is simulated "
                f"({', '.join(_BLOCKS)}), never what is: `systems`, `output`, `sweep` "
                "and the study's own keys belong to each study's config.",
                code="config.option.not_permitted", option=str(key),
                context=DEFAULTS_FILE, permitted=list(_BLOCKS))
        if not isinstance(data[key], dict):
            raise ConfigError(
                f"{where}: `{key}` is a mapping of settings.",
                code="config.option.wrong_type", option=str(key), context=DEFAULTS_FILE,
                expected_type="mapping", found_type=type(data[key]).__name__)
    try:
        # Checked as part of a config: every setting known, every value one
        # the schema accepts. On a copy: the check settles some values in
        # place, and the file is filled in as written.
        import copy

        validate_config(copy.deepcopy(data))
    except ConfigError as exc:
        raise ConfigError(f"{where}: {exc}", code=exc.code,
                          **exc.refusal.details) from exc
    given = {f"{block}.{key}" for block, settings in data.items() for key in settings}
    if not isinstance(decisions, dict):
        raise ConfigError(
            f"{where}: `decisions` is a mapping from a setting's dotted name to "
            "`{why: ...}`.", code="config.option.wrong_type", option="decisions",
            context=DEFAULTS_FILE, expected_type="mapping",
            found_type=type(decisions).__name__)
    whys: dict[str, str] = {}
    for name, said in decisions.items():
        if str(name) not in given:
            raise ConfigError(
                f"{where} gives a reason for `{name}`, which it does not set. A "
                "reason here is for one of your defaults.",
                code="config.option.unknown", option=str(name), context=DEFAULTS_FILE,
                permitted=sorted(given))
        why = said.get("why") if isinstance(said, dict) else said
        if not isinstance(why, str) or not why.strip():
            raise ConfigError(
                f"{where}: the reason for `{name}` is `{{why: ...}}`, a sentence.",
                code="config.option.missing_companion", option=str(name),
                context=DEFAULTS_FILE, requires=["why"])
        whys[str(name)] = " ".join(why.split())
    return YourDefaults(where, {block: dict(settings) for block, settings in data.items()},
                        whys)


def write_defaults(path: str | Path, *, given: dict[str, Any] | None = None,
                   unset: list[str] | tuple[str, ...] | None = None,
                   whys: dict[str, str | None] | None = None) -> YourDefaults:
    """The defaults file at ``path`` written with the settings ``given`` (by
    dotted name), without those ``unset``, and with each reason in ``whys``
    (an empty one removes it), from GUI Settings' This folder.

    What is written is read back as `read_defaults` reads any defaults
    file, before it takes the file's place: a value the schema refuses or
    a setting the file may not hold is refused with the same words, and
    the file is left as it was. A reason goes with its setting: one for a
    setting the file does not give is not written.
    The file is written whole from what it says, its blocks in the order
    the schema gives them and the reasons last, so comments in it would
    be lost: Settings refuses a file with comments before it asks
    (`gui/settings_page.py`). It is created where there is none.
    """
    import os
    import stat
    import tempfile

    import yaml

    from fastmdxplora.config.loader import ConfigError

    where = Path(path)
    data: dict[str, Any] = {}
    if where.is_file():
        try:
            text = where.read_text(encoding="utf-8")
        except OSError as exc:
            raise ConfigError(f"{where} could not be read: {exc}",
                              code="config.file.missing", path=str(where)) from exc
        # Checked first: a file already wrong is said as it is, not
        # written over with what was asked.
        read_defaults(where)
        loaded = yaml.safe_load(text)
        data = dict(loaded) if isinstance(loaded, dict) else {}
    decisions = data.pop("decisions", None) or {}
    decisions = dict(decisions) if isinstance(decisions, dict) else {}

    def split(name: str) -> tuple[str, str]:
        block, _, key = str(name).partition(".")
        if block not in _BLOCKS or not key or "." in key:
            raise ConfigError(
                f"`{name}` is not a setting your defaults may hold: they say how a study "
                f"is simulated ({', '.join(_BLOCKS)}), one setting at a time.",
                code="config.option.not_permitted", option=str(name),
                context=DEFAULTS_FILE, permitted=list(_BLOCKS))
        return block, key

    for name, value in (given or {}).items():
        block, key = split(name)
        held = data.get(block)
        data[block] = dict(held) if isinstance(held, dict) else {}
        data[block][key] = value
    for name in unset or ():
        block, key = split(name)
        held = data.get(block)
        if isinstance(held, dict):
            held = dict(held)
            held.pop(key, None)
            data[block] = held
        decisions.pop(name, None)
    for name, why in (whys or {}).items():
        split(name)
        said = " ".join(str(why or "").split())
        if said:
            decisions[name] = {"why": said}
        else:
            decisions.pop(name, None)
    data = {block: data[block] for block in _BLOCKS
            if isinstance(data.get(block), dict) and data[block]}
    given_now = {f"{block}.{key}" for block, settings in data.items() for key in settings}
    # A reason goes with its setting.
    decisions = {name: said for name, said in decisions.items() if name in given_now}
    if decisions:
        data["decisions"] = decisions
    where.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=".fastmdx-defaults.", suffix=".yml",
                                         dir=where.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as out:
            yaml.safe_dump(data, out, sort_keys=False, allow_unicode=True,
                           default_flow_style=False)
        # Readable as the file was, or as a new file is here: a temporary
        # file is its owner's alone, and a lab's shared file would become so.
        if where.is_file():
            os.chmod(temporary, stat.S_IMODE(where.stat().st_mode))
        else:
            mask = os.umask(0)
            os.umask(mask)
            os.chmod(temporary, 0o666 & ~mask)
        try:
            read_defaults(temporary)
        except ConfigError as exc:
            raise ConfigError(str(exc).replace(temporary, str(where)), code=exc.code,
                              **exc.refusal.details) from exc
        os.replace(temporary, where)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return read_defaults(where)


def defaults_for(start: str | Path, *, home: str | Path | None = None) -> YourDefaults | None:
    """The defaults that apply to a study in ``start``, or None."""
    found = find_defaults(start, home=home)
    return read_defaults(found) if found is not None else None


def _swept(config: dict[str, Any]) -> set[str]:
    sweep = config.get("sweep")
    return {str(name) for name in sweep} if isinstance(sweep, dict) else set()


def with_defaults(config: dict[str, Any],
                  defaults: YourDefaults | None) -> tuple[dict[str, Any], list[str]]:
    """The config with your defaults filling what it leaves unset, and the
    dotted names of the settings filled. The config given is not changed.

    A value the config gives, `null` included, is kept; a setting the study
    sweeps is left to its sweep; a continuation is left whole, and a study
    prepared from another keeps that one's setup. Each setting filled is
    recorded in `decisions` with the file as its source, unless the config
    already says why it has its value.
    """
    from fastmdxplora.naming import default_system_id

    if defaults is None or not defaults.values:
        return config, []
    simulation = config.get("simulation") if isinstance(config.get("simulation"), dict) else {}
    if simulation.get("resume_from"):
        return config, []
    skipped_blocks = {"setup"} if simulation.get("setup_from") else set()
    swept = _swept(config)
    out = dict(config)
    decisions = dict(out.get("decisions") or {}) if isinstance(
        out.get("decisions"), dict) else {}
    # What each system entry sets of its own, by dotted name. A setting a
    # system gives is that system's, as one the study gives is the study's.
    entries = config.get("systems") if isinstance(config.get("systems"), list) else []
    systems = [dict(entry) if isinstance(entry, dict) else entry for entry in entries]
    own = [{f"{block}.{key}" for block in _BLOCKS
            if isinstance(entry, dict) and isinstance(entry.get(block), dict)
            for key in entry[block]} for entry in systems]
    filled: list[str] = []
    for block, settings in defaults.values.items():
        if block in skipped_blocks:
            continue
        held = out.get(block)
        if held is not None and not isinstance(held, dict):
            continue
        merged = dict(held or {})
        answered = {f"{block}.{key}" for key in merged}
        for key, value in settings.items():
            name = f"{block}.{key}"
            if _answers(name, swept) or _answers(name, answered):
                continue
            set_by = [i for i, names in enumerate(own) if _answers(name, names)]
            if set_by and len(set_by) == len(systems):
                # Every system gives its own: the file's value would be
                # used by no run, and recorded as the reason for values it
                # did not give.
                continue
            why = defaults.whys.get(name) or f"Your usual value, from {DEFAULTS_FILE}."
            if set_by:
                # Some systems give their own: the file's value goes to the
                # others, each in its own block, so it can neither win over a
                # system's other spelling of it nor be said of a system that
                # did not take it.
                took = []
                for i, entry in enumerate(systems):
                    if i in set_by or not isinstance(entry, dict):
                        continue
                    entry_block = entry.get(block)
                    if entry_block is not None and not isinstance(entry_block, dict):
                        continue
                    entry[block] = {**(entry_block or {}), key: value}
                    # Named as the run folders are, so each run can tell
                    # whether the reason is its own (`decisions_for_system`).
                    took.append(str(entry.get("id") or default_system_id(i)))
                if not took:
                    continue
                filled.append(name)
                if name not in decisions:
                    decisions[name] = {
                        "why": f"{why.rstrip().rstrip('.')}{_FOR}"
                               f"{', '.join(took)}{_THE_OTHERS}",
                        "source": DEFAULTS_FILE}
                continue
            merged[key] = value
            filled.append(name)
            if name not in decisions:
                decisions[name] = {"why": why, "source": DEFAULTS_FILE}
        if merged:
            out[block] = merged
    if filled and systems:
        out["systems"] = systems
    if filled:
        out["decisions"] = decisions
    return out, filled


#: How a reason filled for some systems only says which: written by
#: `with_defaults`, read back per run by `decisions_for_system`.
_FOR = ". For "
_THE_OTHERS = "; the other systems give their own."


def decisions_for_system(decisions: Any, system_id: str) -> Any:
    """A study's ``decisions`` as one of its systems' runs records them.

    A default that only some systems took is recorded study-wide as "For
    plain; the other systems give their own." A run of a system that gave
    its own value recorded that reason for it, so the report of a run at
    320 K said the lab's 310 K was why. Here the reason stays only in the
    runs of the systems that took it, said as the file says it.
    """
    if not isinstance(decisions, dict):
        return decisions
    kept: dict[str, Any] = {}
    for name, said in decisions.items():
        why = said.get("why") if isinstance(said, dict) else None
        if (isinstance(why, str) and said.get("source") == DEFAULTS_FILE
                and why.endswith(_THE_OTHERS) and _FOR in why):
            reason, _, names = why[:-len(_THE_OTHERS)].rpartition(_FOR)
            if system_id not in [part.strip() for part in names.split(",")]:
                continue
            said = {**said, "why": f"{reason}."}
        kept[name] = said
    return kept


def refused_with(config: dict[str, Any], filled: dict[str, Any],
                 defaults: YourDefaults | None, *, require_systems: bool = False) -> Any:
    """The refusal your defaults cause, naming the file, where the config is
    accepted without them and refused with them; None otherwise. A refusal
    the config has of its own is left to be said as its own."""
    import copy

    from fastmdxplora.config.loader import ConfigError, validate_config

    if defaults is None or filled is config:
        return None
    try:
        validate_config(copy.deepcopy(config), require_systems=require_systems)
    except ConfigError:
        return None
    try:
        validate_config(copy.deepcopy(filled), require_systems=require_systems)
    except ConfigError as exc:
        return ConfigError(
            f"With your defaults ({defaults.path}) filled in, the study is refused: "
            f"{exc} Set that setting in the study's config, or change {DEFAULTS_FILE}.",
            code=exc.code, **exc.refusal.details)
    return None
