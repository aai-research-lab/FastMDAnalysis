"""New study shows every step from the start, and what the study costs here.

The Config Builder showed one step, Start, until a starting point was
chosen; the other four appeared after, so the form's shape was learned by
clicking. Its Run on this machine sat in the accent colour though it could
not run. And what a study costs was one line, the time, under the picture
of the system.

Every step is shown now, as a study from a structure has them, there to
read and not to change until a starting point is chosen. Run is quiet, as
Send is with nothing written, until there is a study to run. And beside the
summary, as a compiler's diagnostics, what the study costs here: the speed
a run makes on this machine's own timing (ns/day), the time every run
takes, the disk its trajectories take against what is free where it is
written, and the GPU memory a run needs against what the GPUs have free,
as a Run here chooses its GPU. A cost that does not fit is an error, one
that would leave the disk all but full a warning.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from collections import namedtuple

import pytest

from fastmdxplora.gui import preview
from fastmdxplora.gui.preview import DCD_BYTES_PER_ATOM, DCD_BYTES_PER_FRAME, system_preview
from tests.test_the_builder_says_what_setup_will_build import _peptide


@pytest.fixture(autouse=True)
def _own_config(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.setattr(preview, "_platform", lambda requested, precision: (
        "CPU" if requested in ("auto", "CPU") else requested, precision))
    monkeypatch.setattr(preview, "_gpus_now", lambda: None)


def _study(tmp_path, **simulation) -> dict:
    path = tmp_path / "peptide.pdb"
    path.write_text(_peptide(6) + "END\n", encoding="utf-8")
    return {"system": str(path), "output": str(tmp_path / "out"),
            "include_phase": ["setup", "simulation", "analysis"],
            "simulation": {"duration_ns": 1, **simulation}}


def _lines(answer) -> dict:
    assert answer["ok"], answer
    return {line["kind"]: line for line in answer["cost"]}


def test_speed_time_and_disk_are_worked_out_from_this_machine(tmp_path):
    from fastmdxplora.cost import calibrate

    calibrate(particles=30_000, steps=5_000, seconds=42.0, platform_name="CPU", precision="mixed")
    answer = system_preview(_study(tmp_path, platform="CPU"))
    lines = _lines(answer)
    particles = answer["estimate"]["particles"]
    per_step = 42.0 / (30_000 * 5_000) * particles
    # 2 fs a step, as many steps a day as the timing allows.
    per_day = 2e-6 * 86400 / per_step
    assert lines["speed"]["said"].startswith(f"about {per_day:,.3g} ns/day on CPU, for about "
                                            f"{particles:,} particles")
    assert lines["time"]["said"].startswith("about ") and "1,250,000 steps in all" in lines["time"]["said"]
    # 500,000 production steps, a frame every 250 of them; the water left out.
    estimate = answer["estimate"]
    atoms = particles - 3 * estimate["waters"]
    size = 2000 * (DCD_BYTES_PER_ATOM * atoms + DCD_BYTES_PER_FRAME)
    disk = lines["disk"]
    assert disk["level"] == "note"
    assert disk["said"].startswith(f"about {preview._size(size)} of trajectory: 2,000 frames of "
                                   f"{atoms:,} atoms (the atoms but water)")
    assert "free where it is written" in disk["said"]
    assert lines["gpu"]["said"] == "runs on the CPU, so no GPU memory is needed"
    assert all(line["level"] == "note" for line in answer["cost"])


def test_a_disk_too_small_is_an_error_and_one_nearly_full_a_warning(tmp_path, monkeypatch):
    Usage = namedtuple("Usage", "total used free")
    study = _study(tmp_path, platform="CPU", save_selection="all")
    answer = system_preview(study)
    trajectory = 2000 * (DCD_BYTES_PER_ATOM * answer["estimate"]["particles"] + DCD_BYTES_PER_FRAME)
    assert "(every atom)" in _lines(answer)["disk"]["said"]
    monkeypatch.setattr(shutil, "disk_usage", lambda path: Usage(0, 0, trajectory // 2))
    disk = _lines(system_preview(study))["disk"]
    assert disk["level"] == "error" and disk["said"].endswith(": it does not fit")
    monkeypatch.setattr(shutil, "disk_usage", lambda path: Usage(0, 0, int(trajectory * 1.1)))
    disk = _lines(system_preview(study))["disk"]
    assert disk["level"] == "warning" and disk["said"].endswith("all but full")
    # A study run until it is determined is priced at its ceiling.
    study["simulation"]["stop_when"] = {"max_duration_ns": 2, "observables": ["rmsd"]}
    assert " at most: " in _lines(system_preview(study))["disk"]["said"]


def test_the_gpu_memory_a_run_needs_against_what_is_free(tmp_path, monkeypatch):
    from fastmdxplora.remote import gpu_room
    from fastmdxplora.remote.gpu_room import HERE, Gpu, Room

    gpu_room.learn(HERE, job="small", particles=1_000, peak_mb=500, gpu="GPU-0000aaaa-1111")
    gpu_room.learn(HERE, job="large", particles=100_000, peak_mb=2_500, gpu="GPU-0000aaaa-1111")
    room = Room(gpus=(Gpu(0, "GPU-0000aaaa-1111", "RTX 4090", 24_564, 1_000, 23_564, 0),))
    monkeypatch.setattr(preview, "_gpus_now", lambda: room)
    gpu = _lines(system_preview(_study(tmp_path, platform="CUDA")))["gpu"]
    assert gpu["level"] == "note"
    assert gpu["said"].startswith("about ") and " GB of GPU memory a run (from 2 runs in mixed " \
        "precision measured on this computer" in gpu["said"]
    assert gpu["said"].endswith("on GPU 0 (RTX 4090), 23.0 GB free of 24.0 GB")
    full = Room(gpus=(Gpu(0, "GPU-0000aaaa-1111", "RTX 4090", 24_564, 24_464, 100, 0),))
    monkeypatch.setattr(preview, "_gpus_now", lambda: full)
    gpu = _lines(system_preview(_study(tmp_path, platform="CUDA")))["gpu"]
    assert gpu["level"] == "error"
    # No GPU read here: said, not guessed.
    monkeypatch.setattr(preview, "_gpus_now", lambda: None)
    gpu = _lines(system_preview(_study(tmp_path, platform="CUDA")))["gpu"]
    assert gpu == {"kind": "gpu", "level": "note",
                   "said": "not checked: no GPU here that nvidia-smi reads"}


def test_a_machine_never_timed_and_a_study_that_simulates_nothing(tmp_path):
    lines = _lines(system_preview(_study(tmp_path, platform="CPU")))
    assert lines["speed"]["said"].startswith("ns/day not known: this machine has not been timed")
    assert "time" not in lines
    study = _study(tmp_path)
    study["include_phase"] = ["setup"]
    assert system_preview(study)["cost"] == [{
        "kind": "time", "level": "note",
        "said": "Nothing is simulated, so no time, trajectory or GPU memory is counted."}]


LOOK = """() => ({
    links: [...document.querySelectorAll('#run-steps .builder-step-link')]
        .filter((l) => !l.hidden).map((l) => l.dataset.step),
    shown: [...document.querySelectorAll('section.builder-step')].filter((s) => !s.hidden)
        .map((s) => s.dataset.step),
    inert: [...document.querySelectorAll('section.builder-step')].filter((s) => s.inert)
        .map((s) => s.dataset.step),
    waiting: document.querySelectorAll('.builder-step-waiting').length,
    fields: document.querySelectorAll('#run-protocol-body .builder-field').length,
    run: (() => {
        const b = document.getElementById('run-start-button');
        return [b.disabled, getComputedStyle(b).backgroundColor];
    })(),
    cost: [...document.querySelectorAll('#run-cost-lines .builder-cost-line')]
        .map((l) => [l.dataset.level || '', l.dataset.kind || '', l.textContent]),
    costShown: !document.getElementById('run-cost').hidden,
})"""


def test_every_step_is_shown_from_the_start_and_the_cost_beside_the_summary(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.cost import calibrate
    from fastmdxplora.gui.server import start_dashboard_session

    calibrate(particles=30_000, steps=5_000, seconds=42.0, platform_name="CPU", precision="mixed")
    path = tmp_path / "peptide.pdb"
    path.write_text(_peptide(6) + "END\n", encoding="utf-8")
    session = start_dashboard_session(output=str(tmp_path), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#run", wait_until="domcontentloaded")
            page.wait_for_selector("#run-protocol-body .builder-field", state="attached")
            fresh = page.evaluate(LOOK)
            # A click in a step not yet open takes no focus.
            page.locator("#run-protocol-body input").first.click(force=True, timeout=5000)
            focused = page.evaluate("() => !!document.activeElement.closest('#run-step-protocol')")
            page.click('#run-start .builder-start-card >> nth=0')
            page.wait_for_function("() => !document.querySelector('.builder-step-waiting')")
            chosen = page.evaluate(LOOK)
            page.fill("#run-systems .builder-system-structure", str(path))
            page.wait_for_function("() => document.querySelectorAll("
                                   "'#run-cost-lines .builder-cost-line[data-kind]').length === 4")
            costed = page.evaluate(LOOK)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    steps = ["start", "system", "protocol", "analysis", "review"]
    assert fresh["links"] == steps and fresh["shown"] == steps
    assert fresh["inert"] == ["system", "protocol", "analysis", "review"] and fresh["waiting"] == 4
    assert fresh["fields"] > 5 and focused is False
    # Quiet: no accent behind it while it cannot run.
    assert fresh["run"] == [True, "rgba(0, 0, 0, 0)"]
    assert not fresh["costShown"]
    assert chosen["inert"] == [] and chosen["waiting"] == 0 and chosen["costShown"]
    assert chosen["cost"][0][2].startswith("Name a structure")
    kinds = [kind for _level, kind, _said in costed["cost"]]
    assert kinds == ["speed", "time", "disk", "gpu"]
    assert all(level == "note" for level, _kind, _said in costed["cost"])
    assert "ns/day on CPU" in costed["cost"][0][2]
    assert costed["run"][0] is False and costed["run"][1] != "rgba(0, 0, 0, 0)"


def test_a_fast_run_s_speed_is_said_in_figures_not_an_exponent(tmp_path):
    """Three figures were asked of a format that writes 1.23e+03 for a run
    making a thousand ns a day or more: a small system on a GPU."""
    from fastmdxplora.cost import calibrate

    calibrate(particles=30_000, steps=5_000, seconds=0.01, platform_name="CPU", precision="mixed")
    answer = system_preview(_study(tmp_path, platform="CPU"))
    said = _lines(answer)["speed"]["said"]
    particles = answer["estimate"]["particles"]
    per_day = 2e-6 * 86400 / (0.01 / (30_000 * 5_000) * particles)
    assert per_day >= 1000
    import math

    places = 2 - math.floor(math.log10(per_day))
    assert "e+" not in said
    assert said.startswith(f"about {round(per_day, places):,.0f} ns/day on CPU")
    assert preview._three_figures(56.78) == "56.8" and preview._three_figures(0.01234) == "0.0123"


def _disk(config, timed=False):
    from types import SimpleNamespace

    estimate = SimpleNamespace(particles=30_000, waters=9_000, water_model="tip3p")
    time = ({"ok": True, "seconds": 1000.0, "runs": 1, "steps": 1, "platform": "CPU"} if timed
            else {"ok": False, "reason": "This machine has not been timed yet."})
    return next(line["said"] for line in preview.study_cost(config, estimate, time)
                if line["kind"] == "disk")


def test_a_study_run_until_determined_is_counted_at_the_runner_s_own_interval(tmp_path):
    """The interval was worked out from the ceiling's length, but the runner
    picks it from the first piece's production and every continuation keeps
    it: a 10 ns first piece with a 100 ns ceiling writes up to 20,000 frames,
    not 2,000."""
    from fastmdxplora.simulation.lengths import steps_in
    from fastmdxplora.simulation.runner import trajectory_interval_for

    simulation = {"duration_ns": 10, "timestep_fs": 2,
                  "stop_when": {"max_duration_ns": 100, "measure": "rmsd"}}
    said = _disk({"output": str(tmp_path), "simulation": simulation})
    frames = steps_in(100, 2.0) // trajectory_interval_for(steps_in(10, 2.0))
    assert frames == 20_000
    assert f" at most: {frames:,} frames of 3,000 atoms" in said, said
    # A set interval is the runner's too.
    simulation["trajectory_interval_steps"] = 5_000
    assert " at most: 10,000 frames of 3,000 atoms" in _disk({"output": str(tmp_path),
                                                              "simulation": simulation})


def test_every_run_of_a_study_is_counted_on_a_machine_not_timed(tmp_path):
    """The runs were counted from the time line, so on a machine never timed
    a study of five seeds was costed as one run's disk."""
    config = {"output": str(tmp_path), "systems": [{"id": "s", "system": "1UBQ"}],
              "sweep": {"simulation.random_seed": [1, 2, 3, 4, 5]},
              "simulation": {"duration_ns": 10, "timestep_fs": 2}}
    for timed in (False, True):
        assert "in each of 5 runs" in _disk(config, timed=timed)


def test_the_gpu_line_counts_the_gpus_runs_going_here_hold(tmp_path, monkeypatch):
    """A Run here leaves out the GPUs its workspace's runs hold; the cost's
    GPU line counted none held, so a GPU a run just took read as free."""
    from types import SimpleNamespace

    from fastmdxplora import gpu_here, runs_here

    going = [{"folder": str(tmp_path / "other"), "by": "the GUI", "pid": 4242,
              "gpu": {"uuids": ["GPU-1"], "need_mb": 6000}}]
    asked: list = []
    monkeypatch.setattr(runs_here, "going_in", lambda folder, walk=True: going
                        if Path(folder) == tmp_path else [])
    monkeypatch.setattr(preview, "_gpus_now", lambda: SimpleNamespace(gpus=[object()]))
    monkeypatch.setattr("fastmdxplora.remote.send.simulates_on_gpu", lambda config, cpus: True)
    monkeypatch.setattr(gpu_here, "choice_here",
                        lambda config, folder, held, room=None: (asked.append(held), (None, []))[1])
    preview._gpu_line({"output": str(tmp_path / "new"), "simulation": {}}, None)
    assert [[(h.uuid, h.need_mb) for h in held] for held in asked] == [[("GPU-1", 6000)]]


def test_each_run_of_a_sweep_is_counted_from_its_own_settings(tmp_path):
    """Every run was counted as the study's base: a sweep over the
    production length (0.04 and 1 ns) was said as two runs of 200 frames,
    where the runner writes 200 and 2,000."""
    from fastmdxplora.simulation.lengths import steps_in
    from fastmdxplora.simulation.runner import trajectory_interval_for

    config = {"output": str(tmp_path), "systems": [{"id": "s", "system": "1UBQ"}],
              "sweep": {"simulation.duration_ns": [0.04, 1]},
              "simulation": {"duration_ns": 0.04, "timestep_fs": 2}}
    frames = [steps_in(ns, 2.0) // trajectory_interval_for(steps_in(ns, 2.0)) for ns in (0.04, 1)]
    assert frames == [200, 2000]
    said = _disk(config)
    assert "over 2 runs: 200 to 2,000 frames a run, of 3,000 atoms (the atoms but water)" in said, said
    size = sum(f * (DCD_BYTES_PER_ATOM * 3000 + DCD_BYTES_PER_FRAME) for f in frames)
    assert said.startswith(f"about {preview._size(size)} of trajectory over 2 runs"), said
