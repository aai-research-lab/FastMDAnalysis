# src/fastmdanalysis/cli/analyze.py
from __future__ import annotations

from typing import Any, Dict, List, Optional
import argparse
import logging

from ._common import add_file_args, load_options_file, parse_opt_pairs, deep_merge_options


def register(subparsers: argparse._SubParsersAction, common_parser: argparse.ArgumentParser) -> None:
    """
    Register the 'analyze' subcommand.

    Notes
    -----
    - The top-level CLI (main.py) supports constructing the object via a system file:
        FastMDAnalysis(system=<yaml/json or dict>)
      When that path is used, trajectory/topology need not be provided on the CLI.
    - This parser still exposes trajectory/topology for the non-system path.
    """
    p = subparsers.add_parser(
        "analyze",
        parents=[common_parser],
        help="Run multiple analyses (include/exclude) with optional slides",
        conflict_handler="resolve",
    )

    # System config (handled by main.py during construction)
    p.add_argument(
        "-system",
        "--system",
        metavar="FILE",
        help=(
            "YAML/JSON system file. The main entrypoint constructs FastMDAnalysis(system=FILE), "
            "so defaults like include/exclude/options/output/slides/strict/stop_on_error are carried by the instance. "
            "CLI flags here still override those defaults."
        ),
        default=None,
    )

    # Optional trajectory/topology/frames/atoms flags for the non-system path
    add_file_args(p)
    # Make traj/top NOT hard-required here; main.py enforces requirements when no -system is provided.
    for a in p._actions:
        if getattr(a, "dest", None) in {"trajectory", "topology"}:
            a.required = False  # relax; main.py decides based on presence/absence of --system

    # Selection (nargs="+"); also accept comma-separated at runtime (see _coerce_list below)
    p.add_argument(
        "--include", nargs="+",
        help='Analyses to include (default: "all"). Example: --include rmsd rmsf rg',
    )
    p.add_argument(
        "--exclude", nargs="+",
        help="Analyses to exclude from the full set.",
    )

    # Options file and inline overrides
    p.add_argument(
        "--options", type=str, default=None, metavar="FILE",
        help="Path to options file (YAML .yml/.yaml or JSON .json). Matches the API 'options' schema.",
    )
    p.add_argument(
        "--opt", action="append", default=[], metavar="ANALYSIS.PARAM=VALUE",
        help="(Optional) Override/add a specific option (repeatable). Example: --opt rmsd.ref=0",
    )

    # Behavior flags
    p.add_argument(
        "--stop-on-error", action="store_true",
        help="Abort on first analysis error (default: continue).",
    )
    p.add_argument(
        "--slides", nargs="?", const=True, metavar="OUT.pptx",
        help="Create a PowerPoint deck of figures (optionally specify output path).",
    )
    p.add_argument(
        "--strict", action="store_true",
        help="Enable strict mode: raise errors for unknown options (default: log warnings).",
    )

    # Register handler
    p.set_defaults(_handler=handle)


def _coerce_list(val: Optional[List[str]]) -> Optional[List[str]]:
    """Accept nargs='+' lists and also split any comma-separated tokens the user provides."""
    if not val:
        return None
    out: List[str] = []
    for tok in val:
        if tok is None:
            continue
        s = str(tok).strip()
        if not s:
            continue
        if "," in s:
            out.extend([t.strip() for t in s.split(",") if t.strip()])
        else:
            out.append(s)
    return out or None


def handle(args: argparse.Namespace, fastmda, logger: logging.Logger) -> None:
    """
    Handler called by cli/main.py.

    Precedence for options:
      system options (lowest) <- file options <- inline --opt (highest)

    Precedence for other knobs:
      include/exclude/slides/strict/stop_on_error/output:
        CLI value if provided, else instance defaults from system=, else orchestrator defaults.
    """
    # -------------------- Build per-analysis options with precedence --------------------
    # From instance (if constructed with system=)
    system_opts: Dict[str, Dict[str, Any]] = getattr(fastmda, "_system_options", {}) or {}

    # From file
    file_options: Dict[str, Dict[str, Any]] = {}
    if args.options:
        file_options = load_options_file(args.options)

    # From CLI --opt
    cli_overrides = parse_opt_pairs(args.opt)

    # Merge: system <- file <- CLI
    options = deep_merge_options(system_opts, file_options)
    options = deep_merge_options(options, cli_overrides)

    # -------------------- Selection & behavior --------------------
    include = _coerce_list(args.include) or getattr(fastmda, "_system_include", None)
    exclude = _coerce_list(args.exclude) or getattr(fastmda, "_system_exclude", None)

    # Slides: None means "use instance default if any"
    slides = args.slides if args.slides is not None else getattr(fastmda, "_system_slides", None)

    # Strict/stop_on_error: CLI True wins, otherwise instance default
    strict = bool(args.strict or getattr(fastmda, "_system_strict", False))
    stop_on_error = bool(args.stop_on_error or getattr(fastmda, "_system_stop_on_error", False))

    # Output: CLI override > instance default > orchestrator default
    output = args.output if getattr(args, "output", None) else getattr(fastmda, "_system_output", None)

    # -------------------- Run --------------------
    results = fastmda.analyze(
        include=include,
        exclude=exclude,
        options=options if options else None,
        stop_on_error=stop_on_error,
        verbose=True,          # keep progress prints
        slides=slides,          # bool or OUT.pptx or None
        strict=strict,          # strict mode flag
        output=output,          # output directory
    )

    # -------------------- Summary --------------------
    print("\nSummary:")
    for name, res in results.items():
        if name == "slides":
            continue
        status = "OK" if res.ok else f"FAIL ({type(res.error).__name__}: {res.error})"
        print(f"  {name:<10} {status:<50} {res.seconds:6.2f}s")

    # Slides reporting
    if slides:
        sres = results.get("slides")
        if sres and sres.ok:
            print(f"\n[fastmda] Slide deck: {sres.value}")
            logger.info("Slide deck created: %s", sres.value)
        elif sres:
            print(f"\n[fastmda] Slides failed: {sres.error}")
            logger.error("Slides failed: %s", sres.error)
