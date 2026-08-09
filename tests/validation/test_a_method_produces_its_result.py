"""What has to be true of an enhanced-sampling run the pipeline completed.

Run these deliberately, because each makes a real simulation:

    pytest -m network tests/validation -v

Three or four minutes for all three methods on a fourteen-atom peptide. An
ordinary `pytest` skips them.

They are marked `network` alongside the structure corpus not because they
fetch anything -- the peptide is written from a string here -- but because
they are the same kind of test: slow, real, and run on purpose rather than on
every commit.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from .methods import INAPPLICABLE, METHODS, TRIPEPTIDE, Method, names

pytestmark = pytest.mark.network


@pytest.fixture(scope="session")
def method_run(workspace):
    """Run each method once, and remember what it produced."""
    cache: dict[str, dict] = {}

    def _for(name: str) -> dict:
        if name in cache:
            return cache[name]
        method = next(m for m in METHODS if m.name == name)
        directory = workspace / f"method-{name}"
        directory.mkdir(parents=True, exist_ok=True)

        # The peptide beside the config, because the config names it
        # relatively and a study is run from wherever it sits.
        (directory / "tripeptide.pdb").write_text(TRIPEPTIDE, encoding="utf-8")
        config = directory / f"{name}.yml"
        config.write_text(yaml.safe_dump(method.config), encoding="utf-8")

        output = directory / "run"
        completed = subprocess.run(
            [sys.executable, "-m", "fastmdxplora.cli.main", "explore",
             "--config", str(config), "--output", str(output), "--force"],
            capture_output=True, text=True, timeout=3600, cwd=directory)
        cache[name] = {
            "method": method,
            "output": output,
            "returncode": completed.returncode,
            "log": completed.stdout + completed.stderr,
        }
        return cache[name]

    return _for


def _find(result: dict, relative: str) -> Path | None:
    """Locate a file the run produced, by name rather than by path.

    Written as a path at first, and twice wrong: PLUMED writes `HILLS`
    relative to the working directory rather than into `simulation/`, and a
    single-run study keeps its run at the top where a multi-run study nests
    it. What the test means is that the file exists somewhere under the
    study, which is what it now asks.
    """
    name = Path(relative).name
    for candidate in sorted(result["output"].rglob(name)):
        if candidate.is_file():
            return candidate
    return None


class TestEveryMethodCompletes:
    @pytest.mark.parametrize("name", names())
    def test_it_runs(self, name, method_run) -> None:
        result = method_run(name)
        assert result["returncode"] == 0, (
            f"{result['method'].description} failed:\n" + result["log"][-2500:]
        )

    @pytest.mark.parametrize("name", names())
    def test_a_failure_is_explained_rather_than_traced(
        self, name, method_run
    ) -> None:
        """Three methods failed in five minutes of first real use, and each
        printed only which run had stopped. The reason existed in the
        manifest and not on the screen."""
        result = method_run(name)
        if result["returncode"] == 0:
            return
        assert "Traceback" not in result["log"], result["log"][-2500:]


class TestEveryMethodProducesItsResult:
    """Umbrella wrote a curve, metadynamics a surface, steered a work
    record -- and each stopped there. The result a method exists for reached
    a JSON file and no further: no figure, no manifest entry, no report."""

    @pytest.mark.parametrize("name", names())
    def test_the_files_it_promises_are_there(self, name, method_run) -> None:
        result = method_run(name)
        if result["returncode"] != 0:
            pytest.skip("did not run; covered by the completion test")

        method: Method = result["method"]
        for relative in method.produces:
            assert _find(result, relative) is not None, (
                f"{method.name}: nothing named {Path(relative).name} was "
                f"written anywhere under {result['output']}"
            )

    @pytest.mark.parametrize(
        "name", [m.name for m in METHODS if m.figure])
    def test_the_result_is_drawn(self, name, method_run) -> None:
        """A number in a JSON file is not a result anyone reads."""
        result = method_run(name)
        if result["returncode"] != 0:
            pytest.skip("did not run; covered by the completion test")

        figure = result["method"].figure
        assert _find(result, f"{figure}.png") is not None, (
            f"{figure}.png was not drawn")


class TestUmbrellaStitchesItsWindows:
    def test_the_windows_overlap(self, method_run) -> None:
        """Recombination joins histograms where they share ground. Windows
        too far apart, or too soft to hold, leave a gap that no amount of
        sampling closes."""
        result = method_run("umbrella")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "pmf.json").read_text(encoding="utf-8"))
        if record.get("refused"):
            pytest.skip(f"refused: {record['refused'][:120]}")

        overlaps = [o["overlap"] for o in record.get("overlaps", [])]
        assert overlaps, "no overlaps were reported"
        assert min(overlaps) > 0, "adjacent windows shared nothing"

    def test_an_unsampled_bin_is_not_a_number(self, method_run) -> None:
        """They arrived as -kT*log(1e-300) -- 1724 kJ/mol, seven hundred
        times RT -- sitting among neighbours of eleven and thirteen."""
        result = method_run("umbrella")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "pmf.json").read_text(encoding="utf-8"))
        curve = (record.get("pmf") or {}).get("free_energy_kjmol") or []
        for value in curve:
            assert value is None or value < 500, (
                "a free energy of hundreds of kJ/mol on a smooth coordinate "
                "is a clip leaking into the output, not a barrier")

    def test_the_comparison_knows_what_it_compared(self, method_run) -> None:
        """Five windows are one experiment sampled in pieces. Their mean
        radius of gyration rose monotonically with window index because each
        was held further apart, and the report tabulated that as five runs
        disagreeing."""
        result = method_run("umbrella")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        report = result["output"] / "comparison" / "comparison_report.md"
        if not report.is_file():
            pytest.skip("no comparison was written")
        text = report.read_text(encoding="utf-8")
        assert "one experiment sampled in pieces" in text
        assert text.index("Free energy") < text.index("## RMSD")


class TestMetadynamicsJudgesItsOwnConvergence:
    def test_the_hills_shrink(self, method_run) -> None:
        """Well-tempered deposition lowers the hills as a basin fills, and
        that decay is the convergence signal."""
        import numpy as np

        result = method_run("metadynamics")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        hills = np.loadtxt(_find(result, "HILLS"), comments="#")
        assert len(hills) > 10, "too few hills to say anything"
        heights = hills[:, 3]
        assert heights[-5:].mean() < heights[:5].mean(), (
            "the hills are not shrinking, so the bias is not well-tempered")

    def test_a_provisional_surface_is_still_drawn(self, method_run) -> None:
        """A run whose bias has not settled still describes the landscape it
        has filled. Withholding the picture leaves a reader with a sentence
        where they could have had both."""
        result = method_run("metadynamics")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "metadynamics_surface.json")
            .read_text(encoding="utf-8"))
        assert record.get("free_energy_kjmol"), (
            "no surface was returned, refused or not")
        if record.get("refused"):
            assert record.get("provisional") is True

    def test_the_drift_says_where_it_was_judged(self, method_run) -> None:
        """Measured over the whole grid, the number is the worst point --
        always the top of the highest barrier, estimated from a handful of
        visits."""
        result = method_run("metadynamics")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "metadynamics_surface.json")
            .read_text(encoding="utf-8"))
        evidence = record.get("evidence") or {}
        assert "drift_ceiling_kjmol" in evidence
        assert "recrossings_definition" in evidence


class TestASteeredPullReportsItsWork:
    def test_the_work_and_the_curve_are_both_there(self, method_run) -> None:
        """A pull that accumulated work smoothly met resistance all the way;
        one that accumulated it in a step snapped past something. The total
        is the same either way."""
        result = method_run("steered")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "steered_work.json").read_text(encoding="utf-8"))
        assert record["work_kjmol"] > 0
        assert len(record["trajectory"]["work_kjmol"]) > 10

    def test_it_records_where_the_pull_reached(self, method_run) -> None:
        """The anchor travels the whole way and the system lags behind; that
        lag is the dissipation."""
        result = method_run("steered")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "steered_work.json").read_text(encoding="utf-8"))
        assert "requested_to" in record and "to" in record

    def test_it_says_the_work_is_not_a_free_energy(self, method_run) -> None:
        result = method_run("steered")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        record = json.loads(
            _find(result, "steered_work.json").read_text(encoding="utf-8"))
        assert "Jarzynski" in record["work_is_not_a_free_energy"]


class TestAnAnalysisThatDoesNotApplyIsNotAFailure:
    """The half of the gating a corpus of real proteins cannot test, because
    on a real protein everything applies."""

    @pytest.mark.parametrize("analysis", sorted(INAPPLICABLE))
    def test_it_is_not_run_at_all(self, analysis, method_run) -> None:
        result = method_run("steered")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        manifest = json.loads(
            _find(result, "analysis_manifest.json").read_text(encoding="utf-8"))
        assert analysis not in manifest["plan"], (
            f"{analysis} was planned on a tripeptide, where it "
            f"{INAPPLICABLE[analysis]}")

    def test_nothing_else_failed(self, method_run) -> None:
        """An analysis that did not apply was reported as one that broke, and
        the distinction is the whole point."""
        result = method_run("steered")
        if result["returncode"] != 0:
            pytest.skip("did not run")

        manifest = json.loads(
            _find(result, "analysis_manifest.json").read_text(encoding="utf-8"))
        failed = {name: r.get("message", "")
                  for name, r in manifest["results"].items()
                  if r.get("status") != "ok"}
        assert not failed, failed
