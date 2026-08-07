"""Tests for the cross-run comparison report (``batch/compare.py``).

The sandbox can't run real MD/analysis, so these tests fabricate the
on-disk artifacts a completed batch would leave behind — per-run
``analysis/<name>/<name>.dat`` files plus a ``batch_manifest.json`` — and
exercise the comparison builder against them.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from fastmdxplora.batch.compare import build_comparison_report


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_batch(
    tmp_path: Path,
    *,
    temps=(300, 310, 320),
    analyses=("rmsd", "rg"),
    statuses=None,
    n_frames=50,
) -> Path:
    """Fabricate a batch output dir with per-run analysis data + manifest."""
    root = tmp_path / "batch"
    runs = []
    statuses = statuses or {}
    rng = np.random.default_rng(0)
    for T in temps:
        rid = f"trpcage__temperature-K-{T}"
        rout = root / "runs" / rid
        status = statuses.get(T, "ok")
        if status == "ok":
            for analysis in analyses:
                d = rout / "analysis" / analysis
                d.mkdir(parents=True, exist_ok=True)
                base = 0.2 + 0.002 * (T - 300)
                series = base + 0.01 * rng.standard_normal(n_frames).cumsum() / 10
                np.savetxt(d / f"{analysis}.dat", series, fmt="%.8e")
        runs.append({
            "run_id": rid, "system": "prot.pdb", "status": status,
            "output_dir": str(rout),
            "sweep_values": {"setup.temperature_K": T},
            "phases": [], "message": "",
        })
    root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "tool": "FastMDXplora", "n_runs": len(temps),
        "sweep": {"setup.temperature_K": list(temps)},
        "runs": runs,
    }
    (root / "batch_manifest.json").write_text(json.dumps(manifest, indent=2))
    return root


# ---------------------------------------------------------------------------
# Core behaviour
# ---------------------------------------------------------------------------
class TestComparisonReport:
    def test_builds_full_report(self, tmp_path):
        root = _make_batch(tmp_path)
        cmp_dir = build_comparison_report(root)
        assert cmp_dir is not None
        names = {p.name for p in cmp_dir.iterdir()}
        # Overlays + trends for both analyses, plus CSV + markdown
        assert "overlay_rmsd.png" in names
        assert "overlay_rg.png" in names
        assert "trend_rmsd.png" in names
        assert "trend_rg.png" in names
        assert "comparison_summary.csv" in names
        assert "comparison_report.md" in names

    def test_summary_csv_one_row_per_run(self, tmp_path):
        root = _make_batch(tmp_path)
        cmp_dir = build_comparison_report(root)
        with (cmp_dir / "comparison_summary.csv").open() as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 3
        # The swept axis column is present and numeric scalars are filled
        assert "temperature_K" in rows[0]
        assert rows[0]["rmsd_mean"]  # non-empty
        assert {r["temperature_K"] for r in rows} == {"300", "310", "320"}

    def test_markdown_mentions_analyses_and_trend(self, tmp_path):
        root = _make_batch(tmp_path)
        cmp_dir = build_comparison_report(root)
        md = (cmp_dir / "comparison_report.md").read_text(encoding="utf-8")
        assert "Cross-run comparison report" in md
        assert "3 successful runs" in md
        assert "RMSD" in md
        assert "temperature_K" in md
        # A quantitative takeaway sentence is present
        assert "increases" in md or "decreases" in md or "is flat" in md
        # Figures are referenced by relative name
        assert "overlay_rmsd.png" in md
        assert "trend_rmsd.png" in md

    def test_trend_uses_first_sweep_axis(self, tmp_path):
        root = _make_batch(tmp_path)
        cmp_dir = build_comparison_report(root)
        # Trend figure exists and is a non-trivial PNG
        trend = cmp_dir / "trend_rmsd.png"
        assert trend.is_file() and trend.stat().st_size > 1000


# ---------------------------------------------------------------------------
# Graceful degradation
# ---------------------------------------------------------------------------
class TestComparisonDegradation:
    def test_no_manifest_returns_none(self, tmp_path):
        assert build_comparison_report(tmp_path / "nonexistent") is None

    def test_fewer_than_two_ok_runs_returns_none(self, tmp_path):
        # Only one successful run -> nothing to compare
        root = _make_batch(tmp_path, temps=(300,))
        assert build_comparison_report(root) is None

    def test_errored_runs_excluded(self, tmp_path):
        # Three runs but one errored; the other two still compare
        root = _make_batch(tmp_path, statuses={320: "error"})
        cmp_dir = build_comparison_report(root)
        assert cmp_dir is not None
        with (cmp_dir / "comparison_summary.csv").open() as fh:
            rows = list(csv.DictReader(fh))
        # Only the 2 ok runs appear
        assert {r["temperature_K"] for r in rows} == {"300", "310"}

    def test_missing_analysis_data_skipped(self, tmp_path):
        # Runs exist and are ok, but have no analysis .dat files
        root = tmp_path / "batch"
        runs = []
        for T in (300, 310):
            rid = f"r-{T}"
            (root / "runs" / rid).mkdir(parents=True, exist_ok=True)
            runs.append({
                "run_id": rid, "system": "p.pdb", "status": "ok",
                "output_dir": str(root / "runs" / rid),
                "sweep_values": {"setup.temperature_K": T},
                "phases": [], "message": "",
            })
        (root / "batch_manifest.json").write_text(json.dumps({
            "sweep": {"setup.temperature_K": [300, 310]}, "runs": runs,
        }))
        # No analysis data anywhere -> nothing comparable -> None
        assert build_comparison_report(root) is None

    def test_no_sweep_still_overlays(self, tmp_path):
        """A multi-system batch with no sweep: overlays yes, trends no."""
        root = tmp_path / "batch"
        runs = []
        rng = np.random.default_rng(1)
        for sysname in ("wt", "mutant"):
            rid = sysname
            d = root / "runs" / rid / "analysis" / "rmsd"
            d.mkdir(parents=True, exist_ok=True)
            np.savetxt(d / "rmsd.dat", 0.2 + 0.01 * rng.standard_normal(40).cumsum() / 10)
            runs.append({
                "run_id": rid, "system": f"{sysname}.pdb", "status": "ok",
                "output_dir": str(root / "runs" / rid),
                "sweep_values": {}, "phases": [], "message": "",
            })
        (root / "batch_manifest.json").write_text(json.dumps({
            "sweep": {}, "runs": runs,
        }))
        cmp_dir = build_comparison_report(root)
        assert cmp_dir is not None
        names = {p.name for p in cmp_dir.iterdir()}
        # Overlay exists; no trend (no numeric sweep axis)
        assert "overlay_rmsd.png" in names
        assert "trend_rmsd.png" not in names


# ---------------------------------------------------------------------------
# Auto-trigger through BatchExplorer
# ---------------------------------------------------------------------------
class TestComparisonAutoTrigger:
    def test_batch_builds_comparison_when_data_present(self, tmp_path, monkeypatch):
        """A multi-run batch auto-builds the comparison when analysis data exists.

        We stub the per-run execution so the run dirs get fake analysis
        .dat files (the sandbox can't run real MD), then confirm
        BatchExplorer's post-run hook produces the comparison report.
        """
        import numpy as np
        from fastmdxplora.batch import explorer as exp_mod
        from fastmdxplora.orchestrator import RunResult

        def fake_execute_run(spec_dict, run_out, include, exclude, verbose,
                             device, quiet=True, force=False):
            out = Path(run_out)
            T = spec_dict["sweep_values"]["setup.temperature_K"]
            for analysis in ("rmsd", "rg"):
                d = out / "analysis" / analysis
                d.mkdir(parents=True, exist_ok=True)
                base = 0.2 + 0.002 * (T - 300)
                np.savetxt(d / f"{analysis}.dat",
                           base + 0.01 * np.arange(30) / 30)
            return RunResult(
                run_id=spec_dict["run_id"], system=spec_dict["system"],
                status="ok", output_dir=out,
                sweep_values=spec_dict["sweep_values"], phases=[],
            )

        monkeypatch.setattr(exp_mod, "_execute_run", fake_execute_run)

        cfg = tmp_path / "c.yml"
        cfg.write_text(f"""
output: {tmp_path / 'b'}
include: [setup, analysis]
systems:
  - {{id: a, system: prot.pdb}}
sweep:
  setup.temperature_K: [300, 310, 320]
""")
        from fastmdxplora import FastMDXplora
        FastMDXplora(config=str(cfg)).explore()

        cmp_dir = tmp_path / "b" / "comparison"
        assert cmp_dir.is_dir()
        assert (cmp_dir / "comparison_report.md").is_file()
        assert (cmp_dir / "overlay_rmsd.png").is_file()
        assert (cmp_dir / "trend_rmsd.png").is_file()

    def test_comparison_can_be_disabled(self, tmp_path, monkeypatch):
        import numpy as np
        from fastmdxplora.batch import explorer as exp_mod
        from fastmdxplora.orchestrator import RunResult

        def fake_execute_run(spec_dict, run_out, include, exclude, verbose,
                             device, quiet=True, force=False):
            out = Path(run_out)
            d = out / "analysis" / "rmsd"
            d.mkdir(parents=True, exist_ok=True)
            np.savetxt(d / "rmsd.dat", np.arange(30) / 100)
            return RunResult(
                run_id=spec_dict["run_id"], system=spec_dict["system"],
                status="ok", output_dir=out,
                sweep_values=spec_dict["sweep_values"], phases=[],
            )

        monkeypatch.setattr(exp_mod, "_execute_run", fake_execute_run)

        cfg = tmp_path / "c.yml"
        cfg.write_text(f"""
output: {tmp_path / 'b'}
include: [setup, analysis]
systems:
  - {{id: a, system: prot.pdb}}
sweep:
  setup.temperature_K: [300, 310]
report:
  comparison: false
""")
        from fastmdxplora import FastMDXplora
        FastMDXplora(config=str(cfg)).explore()
        # comparison disabled -> no comparison dir
        assert not (tmp_path / "b" / "comparison").exists()

    def test_compare_method_rebuilds(self, tmp_path):
        """fmdx.compare(output_dir=...) rebuilds the report standalone."""
        root = _make_batch(tmp_path)
        from fastmdxplora import FastMDXplora
        fmdx = FastMDXplora(config_data={
            "output": str(root),
            "systems": [{"id": "r", "system": "p.pdb"}],
            "sweep": {"setup.temperature_K": [300, 310, 320]},
        })
        cmp_dir = fmdx.compare(output_dir=root)
        assert cmp_dir is not None
        assert (cmp_dir / "comparison_report.md").is_file()

    def test_compare_without_dir_needs_explore_first(self, tmp_path):
        """compare() with no output_dir and no prior run errors clearly."""
        from fastmdxplora import FastMDXplora
        fmdx = FastMDXplora(config_data={
            "systems": [{"id": "r", "system": "p.pdb"}],
            "sweep": {"setup.temperature_K": [300, 310]},
        })
        # config-driven object hasn't run yet -> no output_dir known
        with pytest.raises(ValueError, match="output directory"):
            fmdx.compare()
