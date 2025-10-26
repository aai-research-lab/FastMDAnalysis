# FastMDAnalysis/src/fastmdanalysis/analysis/analyze.py
"""
Unified analysis orchestrator for FastMDAnalysis.

This module provides a bound method `FastMDAnalysis.analyze(...)` that:
- Runs multiple analysis routines in a single call.
- Supports include/exclude selection with a canonical default order.
- Accepts per-analysis keyword options (filtered against each method's signature).
- Optionally builds a PowerPoint slide deck of figures produced during the run
  via the top-level `slides` argument (bool or explicit output path).
- Collects all generated folders/files into a single analyze output directory.

Binding (done once, typically in fastmdanalysis/__init__.py):
----------------------------------------------------------------
from .analysis.analyze import analyze as _analyze
FastMDAnalysis.analyze = _analyze
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, List, Tuple, Union
import inspect
import warnings
import time
import shutil
import math

# Slide deck utilities (timestamped filename handled inside slideshow.py)
from ..utils.slideshow import slide_show, gather_figures


# Canonical analysis names in preferred execution order.
_DEFAULT_ORDER: Tuple[str, ...] = (
    "rmsd",
    "rmsf",
    "rg",
    "hbonds",
    "ss",
    "sasa",
    "dimred",
    "cluster",
)


@dataclass
class AnalysisResult:
    """Container for per-analysis outcomes."""
    name: str
    ok: bool
    value: Any = None
    error: Optional[BaseException] = None
    seconds: float = 0.0


def _discover_available(self) -> List[str]:
    """Return the subset of _DEFAULT_ORDER implemented on this instance."""
    available: List[str] = []
    for name in _DEFAULT_ORDER:
        meth = getattr(self, name, None)
        if callable(meth):
            available.append(name)
    return available


def _validate_options(options: Optional[Mapping[str, Mapping[str, Any]]]) -> Dict[str, Dict[str, Any]]:
    """Ensure options is a nested mapping {analysis_name: {kw: value}}."""
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise TypeError("options must be a mapping of {analysis_name: {kw: value}}")
    norm: Dict[str, Dict[str, Any]] = {}
    for analysis, kwargs in options.items():
        if not isinstance(kwargs, Mapping):
            raise TypeError(f"options['{analysis}'] must be a mapping of keyword arguments")
        norm[analysis] = dict(kwargs)
    return norm


def _final_list(
    available: Sequence[str],
    include: Optional[Sequence[str]],
    exclude: Optional[Sequence[str]],
) -> List[str]:
    """
    Resolve final ordered list of analyses to run.

    - If include is None or ['all'], start from all available in default order.
    - Else, keep only included ones (preserving _DEFAULT_ORDER ordering).
    - Then drop any in exclude.
    """
    avail_set = set(available)

    if include is None or (len(include) == 1 and str(include[0]).lower() == "all"):
        candidates = [name for name in _DEFAULT_ORDER if name in avail_set]
    else:
        want = {s.lower() for s in include}
        unknown = want - set(_DEFAULT_ORDER)
        if unknown:
            warnings.warn(
                f"Unknown analyses in include: {sorted(unknown)}; valid names: {_DEFAULT_ORDER}"
            )
        candidates = [name for name in _DEFAULT_ORDER if (name in avail_set and name in want)]

    if exclude:
        drop = {s.lower() for s in exclude}
        candidates = [name for name in candidates if name not in drop]

    if not candidates:
        raise ValueError("No analyses to run after applying include/exclude.")
    return candidates


def _filter_kwargs(callable_obj, kwargs: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Only pass keyword arguments that the callable actually accepts.
    If the callable has **kwargs, pass everything through unchanged.
    """
    if not kwargs:
        return {}
    sig = inspect.signature(callable_obj)
    if any(p.kind == p.VAR_KEYWORD for p in sig.parameters.values()):
        return dict(kwargs)
    accepted = {
        name
        for name, p in sig.parameters.items()
        if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)
    }
    return {k: v for k, v in kwargs.items() if k in accepted}


def _unique_path(dest: Path) -> Path:
    """
    If dest exists, append _1, _2, ... until a free path is found.
    Works for both files and directories.
    """
    if not dest.exists():
        return dest
    stem = dest.stem
    suffix = dest.suffix
    parent = dest.parent
    i = 1
    while True:
        candidate = parent / f"{stem}_{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


def _print_summary(results: Dict[str, AnalysisResult], analyze_outdir: Path) -> None:
    """Pretty-print a compact summary table."""
    if not results:
        return
    names = [k for k in results.keys() if k != "slides"]
    width = max(6, max((len(n) for n in names), default=6))
    print("\nSummary:")
    for n in names:
        r = results[n]
        status = "OK" if r.ok else "FAIL"
        print(f"  {n.ljust(width)}  {status:<7}  {r.seconds:>8.2f}s")
    if "slides" in results:
        s = results["slides"]
        if s.ok and s.value:
            print(f"\n[fastmda] Slide deck: {s.value}")
        elif not s.ok:
            print(f"\n[fastmda] Slide deck: FAILED ({s.error})")
    print(f"\n[fastmda] Output collected in: {analyze_outdir.resolve()}")


def _inject_cluster_defaults(self, opts: Dict[str, Dict[str, Any]], plan: Sequence[str]) -> None:
    """
    Ensure cluster runs with all methods by default and has n_clusters when needed.
    """
    if "cluster" not in plan:
        return
    ck = opts.setdefault("cluster", {})

    # methods default: all
    methods = ck.get("methods", "all")
    if isinstance(methods, str):
        # Allow comma-separated string; normalize
        methods_list = [m.strip().lower() for m in methods.split(",")]
    else:
        methods_list = [str(m).lower() for m in methods]

    if "all" in methods_list:
        methods_list = ["dbscan", "kmeans", "hierarchical"]
    ck["methods"] = methods_list

    # n_clusters default if needed
    if any(m in ("kmeans", "hierarchical") for m in methods_list) and "n_clusters" not in ck:
        try:
            n_frames = int(getattr(self.traj, "n_frames", 0))
        except Exception:
            n_frames = 0
        # Heuristic: clamp between 2 and 6; ~sqrt(n_frames/10) if available, else 3
        if n_frames > 0:
            guess = max(2, min(6, int(round(math.sqrt(max(1.0, n_frames / 10.0))))))
        else:
            guess = 3
        ck["n_clusters"] = guess


def run(
    self,
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    options: Optional[Mapping[str, Mapping[str, Any]]] = None,
    stop_on_error: bool = False,
    verbose: bool = True,
    slides: Optional[Union[bool, str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
) -> Dict[str, AnalysisResult]:
    """
    Execute multiple analyses on the current FastMDAnalysis instance and
    collect all generated outputs into a single directory.

    Parameters
    ----------
    include, exclude, options, stop_on_error, verbose, slides
        See docstring above.
    output : str | Path | None
        Parent directory into which all per-analysis output folders and the slide deck
        will be moved. Defaults to "analyze_output".

    Returns
    -------
    Dict[str, AnalysisResult]
        Per-analysis outcomes. If slides were requested, includes a "slides" entry.
    """
    available = _discover_available(self)
    plan = _final_list(available, include, exclude)
    opts = _validate_options(options)

    # Inject cluster defaults so kmeans & hierarchical run by default
    _inject_cluster_defaults(self, opts, plan)

    analyze_outdir = Path(output) if output is not None else Path("analyze_output")
    analyze_outdir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, AnalysisResult] = {}

    if verbose:
        print(f"[FastMDAnalysis] Running {len(plan)} analyses: {', '.join(plan)}")

    run_t0 = time.time()

    # Record original outdirs to move later
    per_analysis_outdirs: Dict[str, Path] = {}

    for name in plan:
        fn = getattr(self, name, None)
        if not callable(fn):
            warnings.warn(f"Skipping '{name}' (not implemented on this instance).")
            continue

        kw = _filter_kwargs(fn, opts.get(name, {}))

        if verbose and opts.get(name):
            dropped = set(opts[name].keys()) - set(kw.keys())
            if dropped:
                warnings.warn(f"Ignoring unsupported options for '{name}': {sorted(dropped)}")

        if verbose:
            print(f"  • {name}() ...", end="", flush=True)

        t0 = time.perf_counter()
        try:
            value = fn(**kw)
            ok = True
            err = None
            try:
                outdir = Path(getattr(value, "outdir"))
                if outdir.exists():
                    per_analysis_outdirs[name] = outdir
            except Exception:
                pass
        except BaseException as e:
            ok = False
            value = None
            err = e
            if verbose:
                print(" failed")
            if stop_on_error:
                raise
        finally:
            dt = time.perf_counter() - t0

        results[name] = AnalysisResult(name=name, ok=ok, value=value, error=err, seconds=dt)
        if verbose and ok:
            print(f" done ({dt:.2f}s)")

    # ---- Move per-analysis outputs under analyze_outdir/<analysis> -----------
    moved_dirs: List[Path] = []
    for name, src_dir in per_analysis_outdirs.items():
        if src_dir.exists():
            dest_dir = analyze_outdir / name
            # ensure unique (won't clobber prior runs)
            if dest_dir.exists():
                dest_dir = _unique_path(dest_dir)
            # Skip if already there
            try:
                if src_dir.resolve() == dest_dir.resolve():
                    moved_dirs.append(dest_dir)
                    continue
            except Exception:
                pass
            dest_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src_dir), str(dest_dir))
            moved_dirs.append(dest_dir)
            # Patch the analysis object's outdir if available
            if results[name].ok and results[name].value is not None:
                try:
                    setattr(results[name].value, "outdir", dest_dir)
                except Exception:
                    pass
            if results[name].value is None:
                results[name].value = {"outdir": dest_dir}

    # ---- Slides: create AFTER moving, scanning only moved_dirs (no duplicates)
    if slides:
        t0 = time.perf_counter()
        try:
            # Only scan the folders produced in THIS run
            roots: List[Union[str, Path]] = moved_dirs.copy()
            images = gather_figures(roots, since_epoch=run_t0 - 5)  # tolerant window

            if not images:
                raise FileNotFoundError("No figures found to include in slide deck.")

            # Let slideshow pick a timestamped name, then we relocate if needed
            deck_path = slide_show(
                images=images,
                outpath=None,
                title="FastMDAnalysis — Analysis Slides",
                subtitle=f"{len(images)} figure(s) — generated {time.strftime('%Y-%m-%d %H:%M:%S')}",
            )
            deck_path = Path(deck_path)

            # Move deck into analyze_outdir if not already there
            if deck_path.parent.resolve() != analyze_outdir.resolve():
                dest = analyze_outdir / deck_path.name
                if dest.exists():
                    dest = _unique_path(dest)
                shutil.move(str(deck_path), str(dest))
                deck_path = dest

            results["slides"] = AnalysisResult(
                name="slides", ok=True, value=deck_path, error=None, seconds=time.perf_counter() - t0
            )
            if verbose:
                print(f"[FastMDAnalysis] Slides created: {deck_path}")
        except BaseException as e:
            results["slides"] = AnalysisResult(
                name="slides", ok=False, value=None, error=e, seconds=time.perf_counter() - t0
            )
            if verbose:
                warnings.warn(f"Slide creation failed: {e}")

    if verbose:
        _print_summary(results, analyze_outdir)

    return results


def analyze(
    self,
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    options: Optional[Mapping[str, Mapping[str, Any]]] = None,
    stop_on_error: bool = False,
    verbose: bool = True,
    slides: Optional[Union[bool, str, Path]] = None,
    output: Optional[Union[str, Path]] = None,
) -> Dict[str, AnalysisResult]:
    """
    Public façade so callers can do: fastmda.analyze(...)
    """
    return run(
        self,
        include=include,
        exclude=exclude,
        options=options,
        stop_on_error=stop_on_error,
        verbose=verbose,
        slides=slides,
        output=output,
    )