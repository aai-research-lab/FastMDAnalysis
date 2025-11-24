#src/fastmdanalysis/cli/_common.py
from __future__ import annotations

import os
import sys
import json
import glob
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Sequence

import argparse


def make_common_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument(
        "--frames", type=str, default=None,
        help="Frame selection as 'start,stop,stride' (e.g., '0,-1,10'). Negative indices allowed.",
    )
    p.add_argument(
        "--atoms", type=str, default=None,
        help='Global atom selection string (e.g., "protein", "protein and name CA").',
    )
    p.add_argument(
        "--verbose", action="store_true",
        help="Print detailed log messages to the screen.",
    )
    return p


def add_file_args(subparser: argparse.ArgumentParser) -> None:
    """
    Add standard IO args.

    --trajectory accepts one or more values and supports:
      - space-separated file paths:   -traj a.dcd b.dcd
      - comma-separated tokens:       -traj "a.dcd,b.dcd"
      - glob patterns:                -traj "Q*.dcd"
    All tokens are expanded (env vars, ~, globs) and validated to exist.
    """
    subparser.add_argument(
        "-traj", "--trajectory",
        nargs="+", required=True,
        help="One or more trajectory files (space- or comma-separated) and/or glob patterns."
    )
    subparser.add_argument(
        "-top", "--topology",
        required=True,
        help="Path to topology file"
    )
    subparser.add_argument(
        "-o", "--output",
        default=None,
        help="Output directory name"
    )


def setup_logging(output_dir: str, verbose: bool, command: str) -> logging.Logger:
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_filename = Path(output_dir) / f"{command}.log"
    
    # Clear any existing handlers to avoid duplicates
    logger = logging.getLogger()
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()
    
    # Create new handlers
    file_handler = logging.FileHandler(log_filename, mode="w")
    stream_handler = logging.StreamHandler(sys.stdout)
    
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[file_handler, stream_handler],
    )
    
    logger.info("FastMDAnalysis command: %s", " ".join(sys.argv))
    return logger


def parse_frames(frames_str: Optional[str]) -> Optional[Tuple[int, int, int]]:
    if not frames_str:
        return None
    try:
        tup = tuple(map(int, frames_str.split(",")))
        if len(tup) != 3:
            raise ValueError
        return tup  # (start, stop, stride)
    except ValueError:
        raise SystemExit("Invalid --frames format. Expected 'start,stop,stride' (e.g., '0,-1,10').")


def _split_commas_and_ws(token: str) -> List[str]:
    """
    Split a token by commas and whitespace. Collapses repeats and strips empties.
    Examples:
      'a,b,c' -> ['a','b','c']
      'a  b, c' -> ['a','b','c']
    """
    parts = re.split(r"[,\s]+", token.strip())
    return [p for p in parts if p]


def _expand_one(piece: str) -> List[str]:
    """Expand ~, env vars, and glob patterns; return list of matches or the path if literal exists."""
    expanded = os.path.expanduser(os.path.expandvars(piece))
    
    # Always try glob first - it handles relative paths correctly
    # even for literal-looking paths without glob characters
    matches = sorted(glob.glob(expanded))
    
    if matches:
        return matches
    
    # Fallback: check if it exists as a literal path
    # Use Path.resolve() to handle relative paths properly
    try:
        resolved_path = Path(expanded).resolve()
        if resolved_path.exists():
            return [str(resolved_path)]
    except (OSError, ValueError):
        pass
    
    return []

def _dedupe_preserve_order(paths: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def expand_trajectory_args(tokens: Sequence[str]) -> List[str]:
    """
    Normalize a list of tokens into a validated list of file paths.
    """
    pieces: List[str] = []
    for tok in tokens:
        for piece in _split_commas_and_ws(tok):
            pieces.append(piece)

    expanded_paths: List[str] = []
    missing_literals: List[str] = []
    for piece in pieces:
        matches = _expand_one(piece)
        if matches:
            expanded_paths.extend(matches)
        else:
            # keep track of unresolved literals or unmatched globs
            missing_literals.append(piece)

    files = _dedupe_preserve_order(expanded_paths)

    if not files:
        msg = "No trajectory files found from input: " + ", ".join(tokens)
        if missing_literals:
            msg += f"\nUnresolved paths/patterns: {', '.join(missing_literals)}"
        raise SystemExit(msg)

    # Warn/raise on partially unresolved inputs (helps catch typos early)
    # Here we choose to raise with a clear message.
    if missing_literals:
        raise SystemExit(
            "Some trajectory paths/patterns did not resolve to existing files:\n  - "
            + "\n  - ".join(missing_literals)
        )

    return files


def normalize_topology_arg(path: str) -> str:
    p = os.path.expanduser(os.path.expandvars(path))
    if not os.path.exists(p):
        raise SystemExit(f"Topology file not found: {path}")
    return p


def build_instance(
    trajectory: Sequence[str] | str,
    topology: str,
    frames: Optional[Tuple[int, int, int]],
    atoms: Optional[str],
):
    """
    Construct FastMDAnalysis. The core class supports either a single path or a list of paths.
    Upstream code should pass the output of expand_trajectory_args(...) for 'trajectory'.
    """
    from fastmdanalysis import FastMDAnalysis
    return FastMDAnalysis(trajectory, topology, frames=frames, atoms=atoms)


def coerce_scalar(s: str) -> Any:
    """Best-effort type coercion for CLI string values."""
    sl = s.lower()
    if sl == "none":
        return None
    if sl in ("true", "yes", "on"):
        return True
    if sl in ("false", "no", "off"):
        return False
    try:
        if "." in s:
            return float(s)
        return int(s)
    except ValueError:
        if "," in s:
            return [coerce_scalar(x) for x in s.split(",")]
        return s


def parse_opt_pairs(pairs: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Parse repeated --opt ANALYSIS.PARAM=VALUE into a nested dict.
    Example:
      --opt rmsd.ref=0 --opt cluster.methods=dbscan,hierarchical --opt cluster.n_clusters=5
    -> {"rmsd":{"ref":0}, "cluster":{"methods":["dbscan","hierarchical"], "n_clusters":5}}
    """
    out: Dict[str, Dict[str, Any]] = {}
    for item in pairs or []:
        if "=" not in item:
            raise SystemExit(f"--opt expects 'analysis.param=value', got: {item}")
        left, value = item.split("=", 1)
        if "." not in left:
            raise SystemExit(f"--opt expects 'analysis.param=value', got: {item}")
        analysis, param = left.split(".", 1)
        analysis = analysis.strip().lower()
        param = param.strip()
        value = coerce_scalar(value.strip())
        out.setdefault(analysis, {})[param] = value
    return out


def validate_options_mapping(obj: Any, src: str) -> Dict[str, Dict[str, Any]]:
    """
    Validate that 'obj' is a mapping: {analysis: {param: value}}.
    Returns a normalized nested dict (analysis names lowercased).
    """
    if not isinstance(obj, dict):
        raise SystemExit(f"{src}: top-level must be a mapping of analyses to option dicts")
    out: Dict[str, Dict[str, Any]] = {}
    for k, v in obj.items():
        if not isinstance(k, str):
            raise SystemExit(f"{src}: analysis keys must be strings; got key={k!r}")
        if not isinstance(v, dict):
            raise SystemExit(f"{src}: value for analysis '{k}' must be a mapping of parameters")
        out[k.lower()] = dict(v)
    return out


def load_options_file(path: str) -> Dict[str, Dict[str, Any]]:
    """
    Load options from YAML or JSON file.
    - Expands ~ and environment variables in the path.
    - Supports .yml/.yaml via PyYAML (if installed) and .json via json stdlib.
    """
    p = Path(os.path.expanduser(os.path.expandvars(path)))
    if not p.exists():
        raise SystemExit(f"--options: file not found: {p}")

    suffix = p.suffix.lower()
    try:
        if suffix in (".yml", ".yaml"):
            try:
                import yaml  # type: ignore
            except Exception:
                raise SystemExit(
                    f"--options: YAML file provided but PyYAML is not installed. "
                    f"Install with: pip install PyYAML"
                )
            with p.open("r", encoding="utf-8") as fh:
                data = yaml.safe_load(fh)
        elif suffix == ".json":
            with p.open("r", encoding="utf-8") as fh:
                data = json.load(fh)
        else:
            raise SystemExit(f"--options: unsupported file type '{suffix}'. Use .yml/.yaml or .json")
    except SystemExit:
        raise
    except Exception as e:
        raise SystemExit(f"--options: failed to parse {p}: {e}")

    return validate_options_mapping(data, src=str(p))


def deep_merge_options(base: Dict[str, Dict[str, Any]], override: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Per-analysis shallow merge; override wins."""
    out = {k: dict(v) for k, v in base.items()}
    for analysis, params in (override or {}).items():
        dst = out.setdefault(analysis, {})
        for k, v in params.items():
            dst[k] = v
    return out


__all__ = [
    "make_common_parser",
    "add_file_args",
    "setup_logging",
    "parse_frames",
    "expand_trajectory_args",
    "normalize_topology_arg",
    "build_instance",
    "coerce_scalar",
    "parse_opt_pairs",
    "validate_options_mapping",
    "load_options_file",
    "deep_merge_options",
]