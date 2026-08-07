"""What a run in progress says about itself.

Three defects of one shape: a page reading a file that does not hold the
answer, and reporting the gap as a value.
"""

from __future__ import annotations

import json
from pathlib import Path
from urllib.request import urlopen

from fastmdxplora.gui import protein_preview
from fastmdxplora.gui.protein_preview import (
    STRUCTURE_CANDIDATES,
    SYSTEM_CANDIDATES,
    find_structure,
    find_system,
)
from fastmdxplora.gui.structure_info import (
    _MAX_PDB_BYTES_FOR_INLINE_SCAN,
    _MAX_PDB_BYTES_FOR_SYSTEM_SCAN,
    count_structure,
)
from fastmdxplora.gui.telemetry import TelemetryWriter, read_status

from tests.test_gui import DashboardConfig, start_test_server


def _protein_pdb() -> str:
    return "\n".join(
        f"ATOM  {index:>5}  CA  ALA A{index:>4}"
        f"{index:>12.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00 10.00           C"
        for index in range(1, 5)
    )


def _solvated_pdb() -> str:
    """Protein, ligand, water and ions: what setup actually builds."""
    return "\n".join([
        _protein_pdb(),
        "HETATM  101  C1  BNZ B   1       5.000   0.000   0.000  1.00 10.00           C",
        "HETATM  102  C2  BNZ B   1       6.000   0.000   0.000  1.00 10.00           C",
        "HETATM  201  O   HOH C   1      10.000   0.000   0.000  1.00 10.00           O",
        "HETATM  202  O   HOH C   2      11.000   0.000   0.000  1.00 10.00           O",
        "HETATM  301 NA    NA D   1      20.000   0.000   0.000  1.00 10.00          NA",
        "END",
    ])


class TestThePanelsReadTheSystemThatWasSimulated:
    """`setup/prepared.pdb` is written before solvation and before the ligand
    is put back, so counting from it describes a different structure than the
    one that ran."""

    def test_the_system_is_looked_for_before_the_stages_toward_it(self) -> None:
        assert SYSTEM_CANDIDATES.index("setup/solvated.pdb") < SYSTEM_CANDIDATES.index(
            "setup/prepared.pdb"
        )
        assert SYSTEM_CANDIDATES.index("simulation/topology.pdb") < SYSTEM_CANDIDATES.index(
            "setup/prepared.pdb"
        )

    def test_what_is_drawn_is_still_the_solute(self, tmp_path: Path) -> None:
        """The viewer's preference is unchanged: a solvated topology sent whole
        to a browser is mostly water."""
        run = tmp_path / "run"
        for rel in ("setup/prepared.pdb", "simulation/topology.pdb"):
            path = run / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(_solvated_pdb(), encoding="utf-8")

        assert find_structure(run) == run / "setup/prepared.pdb"
        assert STRUCTURE_CANDIDATES[0] == "setup/prepared.pdb"

    def test_but_what_is_counted_is_the_solvated_system(self, tmp_path: Path) -> None:
        run = tmp_path / "run"
        prepared = run / "setup" / "prepared.pdb"
        prepared.parent.mkdir(parents=True)
        prepared.write_text(_protein_pdb() + "\nEND\n", encoding="utf-8")
        topology = run / "simulation" / "topology.pdb"
        topology.parent.mkdir(parents=True)
        topology.write_text(_solvated_pdb(), encoding="utf-8")

        assert find_system(run) == topology

    def test_prepared_is_counted_when_it_is_all_there_is(self, tmp_path: Path) -> None:
        """Setup still going: describe what exists rather than nothing."""
        run = tmp_path / "run"
        prepared = run / "setup" / "prepared.pdb"
        prepared.parent.mkdir(parents=True)
        prepared.write_text(_protein_pdb() + "\nEND\n", encoding="utf-8")

        assert find_system(run) == prepared

    def test_the_endpoint_reports_the_ligand_water_and_ions(self, tmp_path: Path) -> None:
        """The whole complaint, end to end: BNZ in explicit solvent read as
        `Ligands none, Water 0, Ions 0`."""
        run = tmp_path / "run"
        prepared = run / "setup" / "prepared.pdb"
        prepared.parent.mkdir(parents=True)
        prepared.write_text(_protein_pdb() + "\nEND\n", encoding="utf-8")
        topology = run / "simulation" / "topology.pdb"
        topology.parent.mkdir(parents=True)
        topology.write_text(_solvated_pdb(), encoding="utf-8")

        server, base_url = start_test_server(run, config=DashboardConfig())
        try:
            payload = json.loads(urlopen(f"{base_url}/api/structure-info", timeout=5).read())
        finally:
            server.shutdown()
            server.server_close()

        assert payload["valid"] is True
        assert "BNZ" in payload["ligand_resnames"]
        assert payload["water_residues"] == 2
        assert payload["ions"] == 1

    def test_a_water_box_is_not_refused_as_too_large(self, tmp_path: Path) -> None:
        """Roughly 100k atoms is 8 MB of PDB, so the viewer's ceiling would
        have answered the counting question with `too-large`."""
        big = tmp_path / "topology.pdb"
        line = (
            "HETATM  201  O   HOH C   1      10.000   0.000   0.000  1.00 10.00           O\n"
        )
        with big.open("w", encoding="utf-8") as handle:
            handle.write(_protein_pdb() + "\n")
            handle.write(line * 130_000)
            handle.write("END\n")
        assert big.stat().st_size > _MAX_PDB_BYTES_FOR_INLINE_SCAN

        assert count_structure(big).get("reason") == "too-large"
        counted = count_structure(big, max_bytes=_MAX_PDB_BYTES_FOR_SYSTEM_SCAN)
        assert counted["valid"] is True
        assert counted["water_residues"] == 1

    def test_the_counting_endpoints_ask_the_counting_question(self) -> None:
        source = Path(protein_preview.__file__).with_name("server.py").read_text(
            encoding="utf-8"
        )
        for marker in ("def _structure_info_payload", "def _ligands_payload"):
            body = source.split(marker, 1)[1].split("\ndef ", 1)[0]
            assert "find_system(root)" in body, marker
            assert "max_bytes=_MAX_PDB_BYTES_FOR_SYSTEM_SCAN" in body, marker


class TestNothingDescribesItselfAsNotAvailable:
    """A placeholder that travels to the page is read as a value."""

    def test_the_writer_does_not_invent_a_stage(self, tmp_path: Path) -> None:
        writer = TelemetryWriter(tmp_path / "run" / "simulation", enabled=True)
        writer.write_status(status="running")

        assert read_status(tmp_path / "run")["stage"] is None

    def test_a_stage_that_was_reported_is_kept(self, tmp_path: Path) -> None:
        writer = TelemetryWriter(tmp_path / "run" / "simulation", enabled=True)
        writer.mark_stage("nvt", "current", status="running")

        assert read_status(tmp_path / "run")["stage"] == "nvt"

    def test_the_page_calls_a_run_without_one_starting(self) -> None:
        source = (
            Path(protein_preview.__file__).with_name("static") / "dashboard.js"
        ).read_text(encoding="utf-8")
        assert 'status.stage || "starting"' in source
        assert 'status.stage || "not available"' not in source


class TestThePrecisionIsRecordedWhileTheRunIsGoing:
    """The simulation manifest carries it, and is written when the run ends."""

    def test_the_writer_records_both(self, tmp_path: Path) -> None:
        writer = TelemetryWriter(
            tmp_path / "run" / "simulation",
            enabled=True,
            platform="CUDA",
            precision="mixed",
            precision_applied=True,
        )
        writer.write_status(status="running")

        status = read_status(tmp_path / "run")
        assert status["precision"] == "mixed"
        assert status["precision_applied"] is True

    def test_the_runner_passes_what_it_selected(self) -> None:
        from fastmdxplora.simulation import runner

        source = Path(runner.__file__).read_text(encoding="utf-8")
        assert "precision=precision," in source
        # `Precision` is a GPU-platform property, so what was requested and
        # what is in force are not the same claim.
        assert 'precision_applied=bool(platform_props.get("Precision"))' in source

    def test_a_none_on_disk_does_not_erase_a_live_value(self, tmp_path: Path) -> None:
        """The orchestrator opens the file to mark the setup phase before the
        runner exists. Merging its defaults back over the runner's own status
        is how platform came to read empty during a run."""
        run = tmp_path / "run"
        orchestrator = TelemetryWriter(run / "simulation", enabled=True)
        orchestrator.mark_stage("setup", "current", status="running")

        runner_writer = TelemetryWriter(
            run / "simulation",
            enabled=True,
            platform="CUDA",
            precision="mixed",
            precision_applied=True,
            total_steps=500,
        )
        runner_writer.write_status(stage="minimization", status="running", current_step=0)

        status = read_status(run)
        assert status["platform"] == "CUDA"
        assert status["precision"] == "mixed"
        assert status["total_planned_steps"] == 500


class TestTheViewerDrawsWhatWasSimulated:
    """The ligand is put back when the system is built, so the prepared solute
    is the one file that never has it."""

    def test_the_served_structure_has_the_ligand(self, tmp_path: Path) -> None:
        run = tmp_path / "run"
        prepared = run / "setup" / "prepared.pdb"
        prepared.parent.mkdir(parents=True)
        prepared.write_text(_protein_pdb() + "\nEND\n", encoding="utf-8")
        topology = run / "simulation" / "topology.pdb"
        topology.parent.mkdir(parents=True)
        topology.write_text(_solvated_pdb(), encoding="utf-8")

        server, base_url = start_test_server(run)
        try:
            body = urlopen(f"{base_url}/structure/topology.pdb", timeout=5).read().decode()
        finally:
            server.shutdown()
            server.server_close()

        assert "BNZ" in body
        assert "ALA" in body
        # ...and without the water box, which is what kept the viewer on the
        # prepared solute in the first place.
        assert "HOH" not in body
        assert " NA D" not in body

    def test_a_run_with_only_a_prepared_solute_still_draws(self, tmp_path: Path) -> None:
        run = tmp_path / "run"
        prepared = run / "setup" / "prepared.pdb"
        prepared.parent.mkdir(parents=True)
        prepared.write_text(_protein_pdb() + "\nEND\n", encoding="utf-8")

        server, base_url = start_test_server(run)
        try:
            body = urlopen(f"{base_url}/structure/topology.pdb", timeout=5).read().decode()
        finally:
            server.shutdown()
            server.server_close()

        assert "ALA" in body

    def test_the_version_stamp_follows_the_served_file(self, tmp_path: Path) -> None:
        """A stamp taken from a file nobody serves leaves a stale copy in the
        browser when the served one changes."""
        run = tmp_path / "run"
        prepared = run / "setup" / "prepared.pdb"
        prepared.parent.mkdir(parents=True)
        prepared.write_text(_protein_pdb() + "\nEND\n", encoding="utf-8")
        topology = run / "simulation" / "topology.pdb"
        topology.parent.mkdir(parents=True)
        topology.write_text(_solvated_pdb(), encoding="utf-8")

        source = Path(protein_preview.__file__).read_text(encoding="utf-8")
        body = source.split("def _structure_info(", 1)[1].split("\ndef ", 1)[0]
        assert "find_system(root)" in body


class TestTheConnectionRowReportsTheConnection:
    """It was filled with the run's own status, so a page whose server had
    gone away went on reporting the last thing a live run said about itself."""

    def _dashboard_js(self) -> str:
        return (
            Path(protein_preview.__file__).with_name("static") / "dashboard.js"
        ).read_text(encoding="utf-8")

    def test_the_run_status_no_longer_fills_it(self) -> None:
        source = self._dashboard_js()
        body = source.split("function applyStatus(", 1)[1].split("\n  function ", 1)[0]
        assert "sidebar-connection-state" not in body
        assert "sidebar-status-dot" not in body

    def test_the_connection_does(self) -> None:
        source = self._dashboard_js()
        body = source.split("function updateConnectionState(", 1)[1].split(
            "\n  function ", 1
        )[0]
        assert "sidebar-connection-state" in body
        for word in ("live", "reconnecting", "lost"):
            assert f'"{word}"' in body, word

    def test_the_stale_mark_is_removed_again(self) -> None:
        """Added and left, a dot reports a fault that has passed."""
        source = self._dashboard_js()
        body = source.split("function updateConnectionState(", 1)[1].split(
            "\n  function ", 1
        )[0]
        assert "classList.toggle" in body
        assert "classList.add" not in body


class TestTheStatusCardsDoNotContradictThePhaseTable:
    """`manifest.json` is written when the run ends, so during a run these
    cards read "Unknown / No run manifest found" and "0 / not available"
    directly above a table correctly saying Setup ok, Simulation running."""

    def _cards(self, root: Path) -> dict[str, tuple[str, str]]:
        from fastmdxplora.gui.report_dashboard import _summary_cards

        built = _summary_cards(
            project_root=root, manifest={}, analysis_manifest={}, sim_manifest={}
        )
        return {card.label: (card.value, card.detail) for card in built}

    def _reached(self, root: Path, *, stage: str) -> None:
        order = ("setup", "minimization", "nvt", "npt", "production")
        writer = TelemetryWriter(root / "simulation", enabled=True)
        for name in order[: order.index(stage)]:
            writer.mark_stage(name, "completed", status="running")
        writer.mark_stage(stage, "current", status="running")

    def test_a_run_in_production_has_finished_setup(self, tmp_path: Path) -> None:
        run = tmp_path / "run"
        self._reached(run, stage="production")

        cards = self._cards(run)
        assert cards["Project status"][0] == "Running"
        assert cards["Phases completed"] == ("1", "setup")

    def test_a_run_still_preparing_has_finished_nothing(self, tmp_path: Path) -> None:
        run = tmp_path / "run"
        self._reached(run, stage="setup")

        cards = self._cards(run)
        assert cards["Project status"][0] == "Running"
        assert cards["Phases completed"][0] == "0"
        # ...and says so as "not yet", not as a verdict about the run.
        assert cards["Phases completed"][1] == "none finished yet"

    def test_an_empty_folder_is_a_different_claim(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty"
        empty.mkdir()

        cards = self._cards(empty)
        assert cards["Project status"][0] == "Not run"
        assert cards["Phases completed"][1] == "nothing has run"

    def test_no_card_reports_absence_as_a_value(self, tmp_path: Path) -> None:
        run = tmp_path / "run"
        self._reached(run, stage="npt")

        for value, detail in self._cards(run).values():
            assert "not available" not in f"{value} {detail}".lower()

    def test_a_finished_run_still_reads_from_its_manifest(self, tmp_path: Path) -> None:
        """The live path is a fallback, not a replacement."""
        from fastmdxplora.gui.report_dashboard import _summary_cards

        run = tmp_path / "run"
        self._reached(run, stage="production")
        manifest = {
            "phases": [
                {"name": "setup", "status": "ok"},
                {"name": "simulation", "status": "ok"},
            ]
        }
        built = _summary_cards(
            project_root=run, manifest=manifest, analysis_manifest={}, sim_manifest={}
        )
        cards = {card.label: (card.value, card.detail) for card in built}
        assert cards["Project status"][1] == "Recorded from manifest"
        assert cards["Phases completed"] == ("2", "setup, simulation")


class TestTheChartsDrawWhatIsRecorded:
    def _text(self, name: str) -> str:
        base = Path(protein_preview.__file__).parent
        folder = "templates" if name.endswith(".html") else "static"
        return (base / folder / name).read_text(encoding="utf-8")

    def test_total_energy_is_recorded(self) -> None:
        from fastmdxplora.gui.telemetry import METRIC_FIELDS, _normalise_energy_header

        assert "total_energy" in METRIC_FIELDS
        assert _normalise_energy_header("Total Energy (kJ/mole)") == "total_energy"

    def test_and_therefore_drawn(self) -> None:
        """It was written to disk from the first run and plotted in none of
        them: potential energy alone does not show whether the integration
        is holding."""
        assert '"total_energy"' in self._text("charts.js")
        assert 'data-chart="total_energy"' in self._text("dashboard.html")

    def test_every_chart_has_a_series_behind_it(self) -> None:
        import re

        drawn = set(re.findall(r'data-chart="(\w+)"', self._text("dashboard.html")))
        configured = set(re.findall(r'\{key: "(\w+)"', self._text("charts.js")))
        assert drawn == configured, drawn ^ configured

    def test_the_title_sits_above_its_canvas(self) -> None:
        """A title drawn over the canvas lands on the y-axis labels."""
        import re

        html = self._text("dashboard.html")
        rows = re.findall(r'<div class="chart-row">(.*?)<canvas', html, re.S)
        assert len(rows) == 5, len(rows)
        for row in rows:
            assert "chart-title" in row

    def test_the_chart_head_is_styled(self) -> None:
        """Markup added without a rule for it stacked as blocks against the
        left margin, and the canvas overflowed a clipping row by the height
        of the title."""
        import re

        html = self._text("dashboard.html")
        css = self._text("dashboard.css")

        used = set(re.findall(r'class="(chart-[\w-]+)"', html))
        used |= set(re.findall(r'class="(chart-[\w-]+) ', html))
        for name in used:
            assert f".{name}" in css, name
        assert "flex-direction: column" in css.split(".chart-row {")[-1][:200]

    def test_each_chart_carries_its_current_value(self) -> None:
        """The number and its trend, in one place. A separate strip repeated
        four of them without their history and was a third register a metric
        had to be added to."""
        import re

        html = self._text("dashboard.html")
        plotted = set(re.findall(r'data-chart="(\w+)"', html))
        valued = set(re.findall(r'data-chart-value="(\w+)"', html))
        assert plotted == valued, plotted ^ valued
        assert 'data-metric="potential_energy"' not in html
        assert "renderMetricCards" not in self._text("dashboard.js")


class TestOneHealthPanel:
    def test_the_duplicate_card_is_gone(self) -> None:
        base = Path(protein_preview.__file__).parent
        html = (base / "templates" / "dashboard.html").read_text(encoding="utf-8")
        js = (base / "static" / "dashboard.js").read_text(encoding="utf-8")

        assert html.count("Simulation health") == 1
        assert "live-health-message" not in html
        # ...and nothing is still writing to it.
        assert "live-health" not in js


class TestThePayloadsDoNotManufactureAVerdict:
    """`_display_value` existed to turn absence into the words "not available",
    which then reached the page as the system name in the top bar."""

    def _running(self, root: Path) -> None:
        writer = TelemetryWriter(root / "simulation", enabled=True, platform="CPU")
        writer.mark_stage("production", "current", status="running")

    def test_the_top_bar_gets_no_name_rather_than_a_wrong_one(
        self, tmp_path: Path
    ) -> None:
        from fastmdxplora.gui.server import _system_info

        run = tmp_path / "run"
        self._running(run)
        assert _system_info(run, {}, {}, {})["system"] == ""

    def test_what_is_known_is_still_reported(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _system_info

        run = tmp_path / "run"
        self._running(run)
        assert _system_info(run, {}, {}, {})["platform"] == "CPU"

    def test_absence_is_a_dash(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _system_info

        run = tmp_path / "run"
        self._running(run)
        info = _system_info(run, {}, {}, {})
        assert info["atoms"] == "—"
        assert "not available" not in " ".join(info.values())

    def test_a_run_without_a_manifest_is_in_progress_not_unknown(
        self, tmp_path: Path
    ) -> None:
        from fastmdxplora.gui.server import _summary_records

        run = tmp_path / "run"
        self._running(run)
        rows = {row["label"]: row["value"] for row in _summary_records(run, {}, {}, {})}
        assert rows["Project status"] == "in progress"
        assert "not available" not in " ".join(rows.values())

    def test_and_an_empty_folder_is_neither(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _summary_records

        empty = tmp_path / "nothing"
        empty.mkdir()
        rows = {row["label"]: row["value"] for row in _summary_records(empty, {}, {}, {})}
        assert rows["Project status"] == "not run"


class TestTheFrameCountIsCountedNotPlanned:
    def test_frames_written_over_steps_actually_run(self) -> None:
        from fastmdxplora.simulation.runner import _frames_written

        assert _frames_written(50_000, 250) == 200
        # A run stopped a fifth of the way in wrote a fifth of the frames.
        assert _frames_written(10_000, 250) == 40
        assert _frames_written(0, 250) == 0
        assert _frames_written(50_000, None) == 0

    def test_the_sampler_reports_it_while_the_run_goes(self) -> None:
        """It used to arrive once, when production ended, so a page watching a
        run showed a dash where the frame count goes for the whole of it."""
        import inspect

        from fastmdxplora.simulation import runner

        signature = inspect.signature(runner._run_md_stage_with_live_metrics)
        assert "trajectory_interval_steps" in signature.parameters

        body = inspect.getsource(runner._run_md_stage_with_live_metrics)
        assert "frame_count=frames_written" in body

    def test_equilibration_reports_no_count_rather_than_zero(self) -> None:
        """No trajectory is being written there, so a zero would be a claim."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        passes = source.count("trajectory_interval_steps=trajectory_interval_steps")
        assert passes == 1, passes

    def test_the_recorded_total_is_not_the_plan_restated(self) -> None:
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        # Once, and only as the denominator: what was planned is the right
        # source for `planned_frames`, and the wrong one for what was written.
        assert source.count('plan["production_steps"] // trajectory_interval_steps') == 1
        planned = source.split("planned_frames = ", 1)[1].split(")", 1)[0]
        assert 'plan["production_steps"] // trajectory_interval_steps' in planned
        assert "production_start_step" in source


class TestARunThatFinishedSaysSo:
    """The project-level writer marks setup, marks analysis and report, and
    records the run ending. It was created only when a config named
    `live_telemetry` outright, so a run leaving it at its default got
    telemetry from the runner and none of those."""

    def _writer_for(self, tmp_path: Path, options: dict):
        from fastmdxplora.orchestrator import FastMDXplora

        class Stub:
            output_dir = tmp_path
            options: dict = {}

        return FastMDXplora._dashboard_writer(Stub(), options)

    def test_the_default_is_enough(self, tmp_path: Path) -> None:
        assert self._writer_for(tmp_path, {"simulation": {}}) is not None

    def test_and_the_schema_is_where_the_default_lives(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        assert PHASE_SCHEMAS["simulation"].get("live_telemetry").default is True

    def test_turning_it_off_still_turns_it_off(self, tmp_path: Path) -> None:
        options = {"simulation": {"live_telemetry": False}}
        assert self._writer_for(tmp_path, options) is None

    def test_a_completed_run_is_not_stale_telemetry(self, tmp_path: Path) -> None:
        """`analyze_health` gates staleness on the run still running, which is
        right -- it fired here only because nothing recorded the end."""
        from fastmdxplora.gui.telemetry import analyze_health

        old = "2020-01-01T00:00:00+00:00"
        running = analyze_health({"status": "running", "last_update_timestamp": old}, [])
        finished = analyze_health({"status": "completed", "last_update_timestamp": old}, [])

        assert running["state"] == "warning"
        assert finished["state"] != "warning"

    def test_a_run_that_does_not_simulate_writes_no_telemetry(
        self, tmp_path: Path
    ) -> None:
        """Telemetry lives under `simulation/`, so an analysis-only run would
        otherwise create that folder and claim to have recorded something."""
        from fastmdxplora.orchestrator import FastMDXplora

        class Stub:
            output_dir = tmp_path
            options: dict = {}

        writer = FastMDXplora._dashboard_writer(
            Stub(), {"simulation": {}}, ["analysis", "report"]
        )
        assert writer is None
        assert not (tmp_path / "simulation").exists()


class TestTheTopBarCarriesProgressNotTrivia:
    def _text(self, name: str) -> str:
        base = Path(protein_preview.__file__).parent
        folder = "templates" if name.endswith(".html") else "static"
        return (base / folder / name).read_text(encoding="utf-8")

    def test_the_platform_is_not_reported_twice(self) -> None:
        """It is static, and the sidebar footer already carries it."""
        html = self._text("dashboard.html")
        assert 'id="sidebar-platform"' in html
        assert 'id="topbar-platform"' not in html

    def test_one_observable_is_not_promoted_above_the_others(self) -> None:
        """Temperature sat in the top bar while density, energy and speed did
        not, and it appeared twice more on the same screen."""
        assert 'id="topbar-temperature"' not in self._text("dashboard.html")

    def test_what_replaced_them_is_progress(self) -> None:
        html = self._text("dashboard.html")
        js = self._text("dashboard.js")
        assert 'id="topbar-progress"' in html
        assert 'id="topbar-eta"' in html
        assert 'setText("topbar-eta", computeETA(status))' in js

    def test_the_output_folder_is_where_the_run_is_identified(self) -> None:
        """It is navigation, not a measurement, so it left the grid of things
        measured about the run."""
        from fastmdxplora.gui.report_dashboard import _summary_cards

        html = self._text("dashboard.html")
        assert 'id="sidebar-output-folder"' in html

        labels = [
            card.label
            for card in _summary_cards(
                project_root=Path("/tmp"),
                manifest={"phases": [{"name": "setup", "status": "ok"}]},
                analysis_manifest={},
                sim_manifest={},
            )
        ]
        assert "Output folder" not in labels


class TestInapplicableIsNotUnavailable:
    """A count has one value, not a sample. The standard-deviation column
    said "not available", which reads as a measurement that was attempted."""

    def test_a_count_has_no_standard_deviation(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.report_dashboard import _metric_rows

        rows = {
            row.metric: row
            for row in _metric_rows(tmp_path, {"n_frames": 200, "n_atoms": 36_843})
        }
        assert rows["Frame count"].average == "200"
        assert rows["Frame count"].stddev == "—"
        assert rows["Atom count"].stddev == "—"

    def test_an_empty_table_says_what_is_missing(self) -> None:
        from fastmdxplora.gui import report_dashboard

        source = Path(report_dashboard.__file__).read_text(encoding="utf-8")
        assert "No analysis outputs to summarise yet." in source
        assert '<td colspan="4" class="table-empty">not available' not in source

    def test_the_report_header_does_not_name_an_absent_system(self) -> None:
        from fastmdxplora.gui import report_dashboard

        source = Path(report_dashboard.__file__).read_text(encoding="utf-8")
        assert 'if system else "not available"' not in source


class TestTheSuiteDoesNotReportABusyMachineAsADefect:
    def test_one_ceiling_for_every_request(self) -> None:
        from fastmdxplora.gui import protein_preview

        tests_dir = Path(protein_preview.__file__).parents[3] / "tests"
        source = (tests_dir / "test_gui.py").read_text(encoding="utf-8")
        assert "HTTP_TIMEOUT = " in source
        # Not one raised number among a dozen that were left alone.
        assert "timeout=5" not in source


class TestTheBannerDescribesTheRunItIsStarting:
    """The frame interval was read from the command line alone. A run driven
    by a config file typed none of it, so the banner computed a default and
    announced a frame every 100 steps for a run writing one every 250 --
    500 frames promised, 200 written."""

    @staticmethod
    def _banner(argv, **fields) -> str:
        import contextlib
        import io
        import sys

        from fastmdxplora.utils.presenter import SessionPresenter

        original = sys.argv
        sys.argv = ["fastmdx"] + argv
        buffer = io.StringIO()
        try:
            with contextlib.redirect_stdout(buffer):
                SessionPresenter(stream=buffer).banner(**fields)
        finally:
            sys.argv = original
        return buffer.getvalue()

    def test_the_config_is_read(self, tmp_path: Path) -> None:
        config = tmp_path / "study.yml"
        config.write_text(
            "simulation:\n"
            "  production_steps: 50000\n"
            "  trajectory_interval_steps: 250\n",
            encoding="utf-8",
        )

        text = self._banner(["explore", "--config", str(config)])

        assert "every 250 production steps" in text
        # 100 is what `trajectory_interval_for(50000)` computes, and what the
        # banner used to print regardless of what the run was told to do.
        assert "every 100 production steps" not in text

    def test_the_command_line_still_wins(self, tmp_path: Path) -> None:
        config = tmp_path / "study.yml"
        config.write_text(
            "simulation:\n  trajectory_interval_steps: 250\n", encoding="utf-8"
        )

        text = self._banner(
            ["explore", "--config", str(config),
             "--simulate-trajectory-interval-steps", "500"]
        )

        assert "every 500 production steps" in text


class TestTheTopBarNamesTheRunFromItsFirstMinute:
    """Reading only the manifest -- written when a run ends -- the top bar
    showed a placeholder for the whole run. Setup records the system before
    simulation starts."""

    def _setup_written(self, root: Path, system: str) -> None:
        import json

        (root / "setup").mkdir(parents=True)
        (root / "setup" / "setup_parameters.json").write_text(
            json.dumps({"phase": "setup", "input": {"system": system}}),
            encoding="utf-8",
        )

    def test_the_name_comes_from_setup_during_a_run(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _system_info

        run = tmp_path / "run"
        self._setup_written(run, "181L")
        assert _system_info(run, {}, {}, {})["system"] == "181L"

    def test_the_manifest_still_wins_once_there_is_one(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _system_info

        run = tmp_path / "run"
        self._setup_written(run, "181L")
        assert _system_info(run, {"system": "1UBQ"}, {}, {})["system"] == "1UBQ"

    def test_with_no_name_anywhere_the_browser_chooses(self, tmp_path: Path) -> None:
        """Empty, not a dash: a dash would override the page's own label."""
        from fastmdxplora.gui.server import _system_info

        run = tmp_path / "run"
        run.mkdir()
        assert _system_info(run, {}, {}, {})["system"] == ""


class TestTheTimelineShowsTheStagesThatRan:
    """A phase block in a config is not a statement about what runs.
    `write_resolved_config` writes only phases with non-empty options, and it
    writes them into the output directory when the run ends -- so a run whose
    analysis and report used defaults showed seven stages while running and
    five the moment it finished."""

    def _config(self, root: Path, text: str) -> None:
        root.mkdir(parents=True, exist_ok=True)
        (root / "resolved_config.yml").write_text(text, encoding="utf-8")

    def _manifest(self, root: Path, *names: str) -> None:
        import json

        root.mkdir(parents=True, exist_ok=True)
        (root / "manifest.json").write_text(
            json.dumps({"phases": [{"name": n, "status": "ok"} for n in names]}),
            encoding="utf-8",
        )

    def test_a_finished_run_shows_what_it_ran(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.telemetry import run_phases, run_stages

        run = tmp_path / "run"
        self._config(run, "setup: {ph: 7.4}\nsimulation: {nvt_steps: 5000}\n")
        self._manifest(run, "setup", "simulation", "analysis", "report")

        assert run_phases(run) == ["setup", "simulation", "analysis", "report"]
        assert len(run_stages(run)) == 7

    def test_a_configured_phase_is_not_an_included_one(self, tmp_path: Path) -> None:
        """Without a manifest, a config naming two phase blocks says nothing
        about the other two -- so assume all rather than hide them."""
        from fastmdxplora.gui.telemetry import run_phases, run_stages

        run = tmp_path / "run"
        self._config(run, "setup: {ph: 7.4}\nsimulation: {nvt_steps: 5000}\n")

        assert run_phases(run) == []
        assert len(run_stages(run)) == 7

    def test_include_is_still_honoured(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.telemetry import run_phases

        run = tmp_path / "run"
        self._config(run, "include: [analysis, report]\n")
        assert run_phases(run) == ["analysis", "report"]

    def test_exclude_is_honoured_too(self, tmp_path: Path) -> None:
        """Stated outright, unlike a phase block, so it can be trusted."""
        from fastmdxplora.gui.telemetry import run_phases

        run = tmp_path / "run"
        self._config(run, "exclude: [report]\n")
        assert "report" not in run_phases(run)
        assert "setup" in run_phases(run)


class TestEveryTabIsReachable:
    """The Ligand tab was clipped off the right edge of the information
    panel: four flex items that refuse to shrink below their labels, on one
    unwrapped line, inside a container with `overflow: hidden`."""

    def _text(self, name: str) -> str:
        base = Path(protein_preview.__file__).parent
        folder = "templates" if name.endswith(".html") else "static"
        return (base / folder / name).read_text(encoding="utf-8")

    def test_the_strip_wraps(self) -> None:
        css = self._text("dashboard.css")
        rule = css.split(".info-tabs {", 1)[1].split("}", 1)[0]
        assert "flex-wrap: wrap" in rule

    def test_a_tab_keeps_its_label(self) -> None:
        """`flex: 1` is `1 1 0%`; a basis of zero with `min-width: auto` is
        what made the line overflow rather than the tabs share it."""
        css = self._text("dashboard.css")
        rule = css.split("\n.info-tab {", 1)[1].split("}", 1)[0]
        assert "flex: 1 1 auto" in rule

    def test_all_four_tabs_are_declared(self) -> None:
        import re

        html = self._text("dashboard.html")
        tabs = re.findall(r'class="info-tab[^"]*"[^>]*data-tab="(\w+)"', html)
        panes = re.findall(r'class="info-pane[^"]*" data-tab="(\w+)"', html)
        assert set(tabs) == {"structure", "simulation", "selection", "ligand"}
        assert set(tabs) == set(panes)


class TestAFactAppearsOnce:
    """Progress, ETA and step were on the page three times each -- the top
    bar, the Simulation progress card, and the Exploration status hero -- so
    a reader had to check three places to see whether they agreed."""

    def _text(self, name: str) -> str:
        base = Path(protein_preview.__file__).parent
        folder = "templates" if name.endswith(".html") else "static"
        return (base / folder / name).read_text(encoding="utf-8")

    def test_the_hero_keeps_only_what_it_says_best(self) -> None:
        html = self._text("dashboard.html")
        assert 'id="hero-status-text"' in html
        assert 'id="hero-stage"' in html
        for gone in ("hero-progress-fill", "hero-progress-pct", "hero-sim-time",
                     "hero-elapsed", "hero-eta", "hero-step"):
            assert gone not in html, gone

    def test_nothing_still_writes_to_them(self) -> None:
        js = self._text("dashboard.js")
        for gone in ("hero-progress-fill", "hero-sim-time", "hero-elapsed",
                     "hero-eta", "hero-step"):
            assert gone not in js, gone

    def test_the_numbers_still_have_a_home(self) -> None:
        """Removed from the hero, not from the page."""
        html = self._text("dashboard.html")
        for kept in ("live-progress-fill", "live-sim-time", "live-elapsed-cell",
                     "live-eta-cell", "live-step-cell", "topbar-progress", "topbar-eta"):
            assert kept in html, kept

    def test_no_cell_opens_with_a_verdict(self) -> None:
        """Placeholders the first poll overwrites. Until then a dash, not a
        claim that the value could not be obtained."""
        import re

        html = self._text("dashboard.html")
        filled = re.findall(r'id="(?:live|hero|topbar|sidebar|health)[\w-]*"[^>]*>([^<]*)<', html)
        assert "not available" not in filled


class TestAFinishedRunIsNotProgressing:
    """Every check in `analyze_health` reads the last sample. When nothing was
    wrong it fell through to "Normal progress / Simulation is progressing
    normally" -- present tense, on a page opened hours after the run ended."""

    def test_a_completed_run_says_so(self) -> None:
        from fastmdxplora.gui.telemetry import analyze_health

        health = analyze_health({"status": "completed"}, [])
        assert health["state"] == "ok"
        assert health["message"] == "Completed"
        assert "progressing" not in health["explanation"]

    def test_a_running_one_still_reads_the_same(self) -> None:
        from fastmdxplora.gui.telemetry import analyze_health

        health = analyze_health({"status": "running"}, [])
        assert health["message"] == "Normal progress"

    def test_a_finished_run_that_went_wrong_still_says_so(self) -> None:
        """The terminal case is about tense, not about suppressing faults."""
        from fastmdxplora.gui.telemetry import analyze_health

        health = analyze_health(
            {"status": "completed"}, [{"temperature": float("nan")}]
        )
        assert health["state"] == "failed"

    def test_the_advice_matches_the_default(self) -> None:
        """The no-telemetry text said the feature was off by default and told
        you to switch it on. It has been on by default since the schema
        changed, so that was a fix for a problem that no longer exists."""
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        from fastmdxplora.gui.telemetry import analyze_health

        assert PHASE_SCHEMAS["simulation"].get("live_telemetry").default is True
        explanation = analyze_health({}, [])["explanation"]
        assert "off by default" not in explanation
        assert "on by default" in explanation


class TestThePanelSaysWhatHeldTheLigand:
    """Three rows read "Requires analysis output" unconditionally -- on a run
    that had produced `analysis/pl_interactions/`, they were asking for
    something already on disk. And they could not distinguish a contact type
    looked for and not found from one never looked for."""

    def _table(self, root: Path, *rows: str) -> None:
        folder = root / "analysis" / "pl_interactions"
        folder.mkdir(parents=True)
        header = (
            "kind,ligand_atom,protein_atom,residue,frames_present,"
            "frames_total,occupancy,episodes,standard_error,well_sampled\n"
        )
        (folder / "pl_interactions.dat").write_text(
            header + "".join(rows), encoding="utf-8"
        )

    def test_an_unanalysed_run_says_so(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _ligand_interactions

        assert _ligand_interactions(tmp_path)["analysed"] is False

    def test_the_kinds_are_counted(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _ligand_interactions

        self._table(
            tmp_path,
            "hydrophobic,1,10,LEU99,180,200,0.90,3,0.02,True\n",
            "hydrophobic,2,14,VAL111,40,200,0.20,7,0.03,False\n",
            "hydrogen_bond,3,22,SER117,120,200,0.60,2,0.02,True\n",
        )
        kinds = _ligand_interactions(tmp_path)["kinds"]

        # Two rows, two residues -- and where the same residue contributes
        # several atom pairs it is still one residue.
        assert kinds["hydrophobic"]["residues"] == 2
        assert kinds["hydrophobic"]["pairs"] == 2
        assert kinds["hydrophobic"]["best_occupancy"] == 0.90
        assert kinds["hydrophobic"]["best_residue"] == "LEU99"
        assert kinds["hbonds"]["residues"] == 1

    def test_looked_for_and_not_found_is_its_own_answer(self, tmp_path: Path) -> None:
        """Different from never looked for, which is what the row used to say
        in both cases."""
        from fastmdxplora.gui.server import _ligand_interactions

        self._table(tmp_path, "hydrophobic,1,10,LEU99,180,200,0.90,3,0.02,True\n")
        measured = _ligand_interactions(tmp_path)

        assert measured["analysed"] is True
        assert measured["kinds"]["salt_bridges"]["residues"] == 0

    def test_a_thinly_observed_contact_is_counted_as_such(self, tmp_path: Path) -> None:
        """A contact seen in a handful of frames is not the observation a bare
        count implies."""
        from fastmdxplora.gui.server import _ligand_interactions

        self._table(
            tmp_path,
            "hydrophobic,1,10,LEU99,180,200,0.90,3,0.02,True\n",
            "hydrophobic,2,14,VAL111,40,200,0.20,7,0.03,False\n",
        )
        assert _ligand_interactions(tmp_path)["kinds"]["hydrophobic"]["thinly_sampled"] == 1

    def test_the_rows_are_no_longer_hardcoded(self) -> None:
        base = Path(protein_preview.__file__).parent
        js = (base / "static" / "dashboard.js").read_text(encoding="utf-8")
        assert "<td>Requires analysis output</td>" not in js
        assert 'interactionCell(info, "salt_bridges")' in js

    def test_atom_pairs_are_not_reported_as_contacts(self, tmp_path: Path) -> None:
        """`pl_interactions.dat` has one row per ligand-atom / protein-atom
        pair. Benzene against three residues produced thirteen rows, and
        reporting that read as though the ligand were held by thirteen
        separate things."""
        from fastmdxplora.gui.server import _ligand_interactions

        rows = [
            f"hydrophobic,{i},{20 + i},LEU99,180,200,0.84,3,0.02,True\n"
            for i in range(1, 8)
        ] + [
            f"hydrophobic,{i},{40 + i},VAL111,30,200,0.15,9,0.03,False\n"
            for i in range(1, 6)
        ] + ["hydrophobic,9,60,ILE78,150,200,0.75,2,0.02,True\n"]
        self._table(tmp_path, *rows)

        hydrophobic = _ligand_interactions(tmp_path)["kinds"]["hydrophobic"]
        assert hydrophobic["pairs"] == 13
        assert hydrophobic["residues"] == 3
        assert hydrophobic["best_residue"] == "LEU99"
        # VAL111's every contact is thinly observed; the other two are not.
        assert hydrophobic["thinly_sampled"] == 1

    def test_the_numbering_names_the_structure_it_belongs_to(self) -> None:
        """The ligand's chain and residue id are OpenMM's, from the solvated
        system. The crystal numbering the setup log reports is not recorded in
        any artifact the page can read, so the rows say which one this is
        rather than leaving it to be mistaken for the other."""
        base = Path(protein_preview.__file__).parent
        js = (base / "static" / "dashboard.js").read_text(encoding="utf-8")
        assert "Chain (simulated)" in js
        assert "<th>Chain</th>" not in js


class TestTheViewerRespondsToClicksAndToggles:
    def _text(self, name: str) -> str:
        base = Path(protein_preview.__file__).parent
        return (base / "static" / name).read_text(encoding="utf-8")

    def test_atoms_are_made_clickable_after_they_exist(self) -> None:
        """3Dmol sets the property on the atoms currently selected. Called at
        viewer creation, before any model is added, it set it on nothing --
        so every atom arrived unclickable and the selection panel waited for
        an event that could not be raised."""
        js = self._text("molecule-viewer.js")
        install = js.split("function installPdb(", 1)[1].split("\n  function ", 1)[0]
        assert "setClickable" in install
        assert "setHoverable" in install

        creation = js.split("function ensureMainViewer(", 1)[1].split("\n  function ", 1)[0]
        assert "setClickable" not in creation

    def test_water_and_ions_are_fetched_rather_than_restyled(self) -> None:
        """They are stripped from the copy the browser gets, so restyling
        could never reveal them."""
        js = self._text("molecule-viewer.js")
        assert "solvent=1" in js
        handler = js.split('data-vis]")', 1)[1].split("\n    });", 1)[0]
        assert 'which === "water"' in handler
        assert "onStructureUpdated" in handler

    def test_the_route_serves_the_solvated_system_on_request(
        self, tmp_path: Path
    ) -> None:
        from urllib.request import urlopen

        run = tmp_path / "run"
        sim = run / "simulation"
        sim.mkdir(parents=True)
        (sim / "topology.pdb").write_text(_solvated_pdb(), encoding="utf-8")

        server, base_url = start_test_server(run)
        try:
            plain = urlopen(f"{base_url}/structure/topology.pdb", timeout=30).read().decode()
            full = urlopen(
                f"{base_url}/structure/topology.pdb?solvent=1", timeout=30
            ).read().decode()
        finally:
            server.shutdown()
            server.server_close()

        assert "HOH" not in plain and "BNZ" in plain
        assert "HOH" in full and "BNZ" in full


class TestTheFileListCanBeReadWithoutGuessing:
    def _text(self, name: str) -> str:
        base = Path(protein_preview.__file__).parent
        folder = "templates" if name.endswith(".html") else "static"
        return (base / folder / name).read_text(encoding="utf-8")

    def test_a_row_names_the_file_within_its_group(self) -> None:
        """Sixteen analyses each write an `options.json`. Listed by bare
        filename they were sixteen identical rows, told apart only by their
        byte counts -- and the payload already carried the path that
        distinguishes them."""
        js = self._text("dashboard.js")
        row = js.split("function fileRowHtml(file) {", 1)[1].split("\n  }", 1)[0]
        assert 'String(file.path || "").split("/").slice(1)' in row
        assert row.index("within") < row.index("file.name")

    def test_the_payload_carries_the_path_that_distinguishes_them(
        self, tmp_path: Path
    ) -> None:
        from fastmdxplora.gui.server import _artifact_records

        for analysis in ("rmsd", "rmsf"):
            target = tmp_path / "analysis" / analysis
            target.mkdir(parents=True)
            (target / "options.json").write_text("{}", encoding="utf-8")

        paths = {record["path"] for record in _artifact_records(tmp_path)}
        assert paths == {"analysis/rmsd/options.json", "analysis/rmsf/options.json"}
        # The names alone do not tell them apart, which is what the row used
        # to show.
        assert {record["name"] for record in _artifact_records(tmp_path)} == {
            "options.json"
        }

    def test_the_viewer_scratch_is_folded_away(self) -> None:
        """The live page writes a rolling PDB history as it goes. Sixty-four
        of them listed flat buried the trajectory, the checkpoint and the
        final state. They now land in the run record, which folds as a whole
        -- the same treatment for the same reason, applied to manifests and
        options files too."""
        from fastmdxplora.gui.server import _artifact_label

        label, group = _artifact_label(
            "simulation/live_frames/frame_000001_nvt_000000001000.pdb"
        )
        assert group == "record"
        assert "scratch" in label.lower()

        js = self._text("dashboard.js")
        assert 'key === "record"' in js
        assert "file-fold" in js
        assert ".file-fold" in self._text("dashboard.css")

    def test_a_row_says_what_the_file_is(self) -> None:
        """Deciding whether to download 85 MB should not require knowing that
        `production.dcd` is the trajectory."""
        from fastmdxplora.gui.server import _artifact_label

        assert _artifact_label("simulation/production.dcd") == (
            "Production trajectory", "simulation",
        )
        assert _artifact_label("report/report.pdf")[0] == "Report (PDF)"
        assert _artifact_label("analysis/rmsd/options.json")[0] == "rmsd: options used"

    def test_the_groups_are_purposes_not_directories(self) -> None:
        """`production.dcd` and the dashboard's own scratch shared a group
        because they share a folder."""
        from fastmdxplora.gui.server import _artifact_label

        assert _artifact_label("simulation/production.dcd")[1] == "simulation"
        assert _artifact_label("simulation/live_status.json")[1] == "record"
        assert _artifact_label("analysis/rmsd/rmsd.png")[1] == "figures"
        assert _artifact_label("analysis/rmsd/rmsd.dat")[1] == "analysis"

    def test_the_record_is_kept_rather_than_hidden(self, tmp_path: Path) -> None:
        from fastmdxplora.gui.server import _artifact_records

        scratch = tmp_path / "simulation" / "live_frames"
        scratch.mkdir(parents=True)
        (scratch / "frame_000001_nvt_000000001000.pdb").write_text("x", encoding="utf-8")

        records = _artifact_records(tmp_path)
        assert len(records) == 1
        assert records[0]["group"] == "record"
        assert records[0]["href"].startswith("/artifacts/")

    def test_the_analysis_tab_does_not_offer_the_report(self) -> None:
        """The markdown report, slides, bundle and analysis manifest are the
        report phase's deliverables, and the Report tab lists all of them with
        size and date. Four of them repeated here made this the place people
        looked."""
        js = self._text("dashboard.js")
        section = js.split("function renderAnalysisSections(payload) {", 1)[1].split(
            "\n  function ", 1
        )[0]
        assert "quick_actions" not in section
        assert "actionHost.hidden = true" in section
