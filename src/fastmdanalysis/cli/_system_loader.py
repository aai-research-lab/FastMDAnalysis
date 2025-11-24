# src/fastmdanalysis/cli/_system_loader.py
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Sequence, Set, Union


# --------------------------- parsing & normalization ---------------------------

def _read_system_file(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    name = path.name.lower()
    if name.endswith((".yaml", ".yml")):
        try:
            import yaml  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "YAML system file provided but PyYAML is not installed. "
                "Install with: pip install PyYAML"
            ) from e
        data = yaml.safe_load(text) or {}
    else:
        data = json.loads(text) or {}
    if not isinstance(data, dict):
        raise ValueError("System file must parse to a mapping/object.")
    return data


def _apply_aliases(cfg: MutableMapping[str, Any]) -> None:
    """Normalize friendly aliases into canonical keys without overwriting explicit canonicals."""
    aliases = {
        "traj": "trajectory",
        "top": "topology",
        "selection": "atoms",
        "outdir": "output",
    }
    for src, dst in aliases.items():
        if src in cfg and dst not in cfg:
            cfg[dst] = cfg[src]


def _normalize_bool_like(x: Any) -> Any:
    if isinstance(x, str):
        s = x.strip().lower()
        if s in {"true", "yes", "1"}:
            return True
        if s in {"false", "no", "0"}:
            return False
    return x


def _resolve_paths(cfg: Mapping[str, Any], base_dir: Path) -> Dict[str, Any]:
    """Resolve trajectory (str|list[str]) and topology (str) relative to base_dir."""
    out: Dict[str, Any] = dict(cfg)

    def _abs(p: str) -> str:
        pp = Path(p)
        return str(pp if pp.is_absolute() else (base_dir / pp).resolve())

    traj = cfg.get("trajectory")
    if traj is not None:
        if isinstance(traj, (list, tuple)):
            out["trajectory"] = [_abs(str(t)) for t in traj]
        else:
            out["trajectory"] = _abs(str(traj))

    top = cfg.get("topology")
    if top is not None:
        out["topology"] = _abs(str(top))

    return out


# --------------------------- CLI integration helpers ---------------------------

def _infer_cli_specified(parser: argparse.ArgumentParser, argv: List[str]) -> Set[str]:
    """
    Return set of action.dest that the user explicitly set on the CLI.
    Works for `--flag value` and `--flag=value` forms; also short options.
    """
    specified: Set[str] = set()
    tokens = list(argv)  # preserve order; we’ll scan
    for act in parser._actions:
        opts = getattr(act, "option_strings", []) or []
        for opt in opts:
            # exact token or assignment form --opt=val
            if opt in tokens or any(t.startswith(opt + "=") for t in tokens):
                specified.add(act.dest)
                break
    return specified


def _coerce_value(value: Any, act: argparse.Action) -> Any:
    """
    Best-effort coercion of YAML scalars/lists to the action's expected type.
    - Obeys nargs ('*' or '+') by wrapping scalars into a list
    - Applies the action.type callable if present
    - Normalizes booleans for common flags
    """
    # Normalize known bool-likes regardless of action.type being None
    if getattr(act, "dest", "") in {"verbose", "strict", "stop_on_error", "slides"}:
        value = _normalize_bool_like(value)

    wants_list = getattr(act, "nargs", None) in ("+", "*")
    v = value
    if wants_list and not isinstance(v, (list, tuple)):
        v = [v]

    t = getattr(act, "type", None)
    if t is not None:
        if isinstance(v, list):
            v = [t(x) for x in v]
        else:
            v = t(v)
    return v


# ------------------------------ public entry point ------------------------------

def merge_system_into_args(
    *,
    parser: argparse.ArgumentParser,
    argv: List[str],
    args: argparse.Namespace,
    system_path: Union[str, Path],
) -> argparse.Namespace:
    """
    Load YAML/JSON and populate args for any dest **not** set on CLI.
    CLI wins; YAML only fills in omissions.

    The YAML can use canonical keys with aliases:
      trajectory|traj, topology|top, selection->atoms, outdir->output.

    Paths for trajectory/topology are resolved relative to the YAML file directory.
    """
    p = Path(system_path).expanduser()
    cfg = _read_system_file(p)

    # Normalize aliases
    _apply_aliases(cfg)

    # Resolve relative paths against the config file directory
    cfg = _resolve_paths(cfg, base_dir=p.parent)

    # Detect which CLI flags the user explicitly set (supports --opt=val form)
    specified = _infer_cli_specified(parser, argv)

    # Fill args for any dest not explicitly set on CLI
    for act in parser._actions:
        dest = getattr(act, "dest", None)
        if not dest or dest == argparse.SUPPRESS or dest == "help":
            continue
        if dest in specified:
            continue  # CLI takes precedence
        if dest in cfg:
            try:
                setattr(args, dest, _coerce_value(cfg[dest], act))
            except Exception:
                setattr(args, dest, cfg[dest])

    # Convenience alias if someone used 'output_dir' in YAML
    if "output_dir" in cfg and not getattr(args, "output", None):
        setattr(args, "output", cfg["output_dir"])

    # Attach the raw dict of per-analysis options for the handler to merge
    if "options" in cfg and not getattr(args, "_system_options", None):
        setattr(args, "_system_options", cfg["options"])

    # Provenance
    setattr(args, "_system_file", str(p))
    return args
