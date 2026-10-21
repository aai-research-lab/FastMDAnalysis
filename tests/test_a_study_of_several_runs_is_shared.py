"""A study of several runs is shared as one file and opened as the study it was.

Every study of a paper is one of several runs: seeds of one preparation,
arms and mutants, umbrella windows. `fastmdx report --share` refused them, so
each run went alone and what joined them (the record of runs, the table of
members, the comparison) did not go at all. Packed now as its folder lays it
out, each run a dataset of its own files and phases, with a prepared system a
run simulated from carried once: kept outside the study, in `prepared/`, its
paths said as that folder, so a run opened elsewhere finds it as it did here.
"""

from __future__ import annotations

import json
import socket
import zipfile
from pathlib import Path, PurePosixPath

import pytest

from fastmdxplora.sharing import PROFILE, STUDY, ShareRefused
from fastmdxplora.sharing.crate import METADATA, check
from fastmdxplora.sharing.opening import open_shared
from fastmdxplora.sharing.pack import share_study

AUTHORS = [("Adekunle Aina", "0000-0002-8215-7452")]
COMPUTER = "labbox07"


def _trajectory(folder: Path) -> None:
    md = pytest.importorskip("mdtraj")
    import numpy as np

    topology = md.Topology()
    residue = topology.add_residue("ALA", topology.add_chain())
    for name in ("N", "CA", "C"):
        topology.add_atom(name, md.element.carbon, residue)
    xyz = np.zeros((4, 3, 3), dtype=np.float32)
    xyz[:, 1, 1] = 0.15
    xyz[:, 2, 2] = 0.15
    trajectory = md.Trajectory(xyz, topology)
    trajectory[0].save_pdb(str(folder / "trajectory_topology.pdb"))
    trajectory.save_dcd(str(folder / "production.dcd"))


def _prepared(folder: Path) -> Path:
    """A prepared system as setup leaves one."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "system.xml").write_text("<System atoms='3'/>", encoding="utf-8")
    (folder / "state.xml").write_text("<State/>", encoding="utf-8")
    (folder / "topology.pdb").write_text("ATOM      1  N   ALA A   1\nEND\n",
                                         encoding="utf-8")
    (folder / "setup_parameters.json").write_text(json.dumps({
        "forcefield": "amber14", "output_dir": str(folder),
        "ligands": [{"name": "BEN", "charge": 1}]}), encoding="utf-8")
    (folder / "ligands").mkdir(exist_ok=True)
    (folder / "ligands" / "BEN.sdf").write_text(f"BEN\n  {folder}\n\nM  END\n",
                                                encoding="utf-8")
    return folder


def _run(root: Path, run_id: str, prepared: Path, given: str, *,
         status: str = "ok") -> Path:
    """A run of the study, simulated from ``prepared``, as the batch lays it out."""
    from fastmdxplora.simulation.pipeline import prepared_system_reference

    run = root / "runs" / run_id
    for folder in ("simulation", "analysis/rmsd", "report", "agent/conversations"):
        (run / folder).mkdir(parents=True, exist_ok=True)
    _trajectory(run / "simulation")
    phases = [{"name": "setup", "status": "skipped", "taken_from": str(prepared)}] + [
        {"name": name, "status": status if name == "analysis" else "ok",
         "output_dir": f"{run}/{name}", "started_at": "2026-10-07T01:00:00+00:00",
         "finished_at": f"2026-10-07T0{2 + len(run_id)}:00:00+00:00",
         "produced_by": {"version": "2.6.0", "host": COMPUTER}}
        for name in ("simulation", "analysis", "report")]
    (run / "manifest.json").write_text(json.dumps({
        "tool": "FastMDXplora", "version": "2.6.0", "system": "1ABC",
        "output_dir": str(run), "phases": phases}), encoding="utf-8")
    (run / "resolved_config.yml").write_text(
        f"systems:\n- system: 1ABC\noutput: {run}\nsimulation:\n  setup_from: {given}\n",
        encoding="utf-8")
    (run / "simulation" / "simulation_parameters.json").write_text(json.dumps({
        "parameters": {"timestep_fs": 2.0, "integrator": "langevin_middle",
                       "setup_from": given},
        "resolved": {"trajectory_interval_steps": 500, "production_steps": 2000},
        "prepared_system": prepared_system_reference(prepared, run, given)}),
        encoding="utf-8")
    (run / "simulation" / "live_status.json").write_text(json.dumps({
        "status": "completed", "host": COMPUTER}), encoding="utf-8")
    (run / "analysis" / "analysis_manifest.json").write_text(json.dumps({
        "trajectory_input": f"{run}/simulation/production.dcd"}), encoding="utf-8")
    (run / "analysis" / "rmsd" / "rmsd.dat").write_text("0 0.1\n1 0.2\n", encoding="utf-8")
    (run / "report" / "report.md").write_text(f"# Run {run_id}\n\nRan on {COMPUTER}.\n",
                                              encoding="utf-8")
    # Never shared: a run's conversations, its process, its tags.
    (run / "agent" / "conversations" / "c1.json").write_text("{}", encoding="utf-8")
    (run / ".fastmdxplora_run.json").write_text('{"pid": 1}', encoding="utf-8")
    (run / "study_tags.json").write_text('{"tags": ["mine"]}', encoding="utf-8")
    return run


def _campaign(tmp_path: Path, *, b_status: str = "ok") -> tuple[Path, Path]:
    """Three seeds of one system: `a` and `c` simulated from a preparation
    kept outside the study (`a` naming its study's folder, `c` the system's
    own), `b` from the study's own shared one. Returns the study and the
    preparation's study."""
    root = tmp_path / "campaign"
    elsewhere = _prepared(tmp_path / "preparation" / "setup").parent
    (elsewhere / "manifest.json").write_text(json.dumps({
        "phases": [{"name": "setup", "status": "ok",
                    "produced_by": {"host": "prepbox3"}}]}), encoding="utf-8")
    shared = _prepared(root / "shared_setup" / "setup")
    _run(root, "a", elsewhere / "setup", str(elsewhere))
    _run(root, "b", shared, str(shared), status=b_status)
    _run(root, "c", elsewhere / "setup", str(elsewhere / "setup"))
    (root / "comparison").mkdir()
    (root / "comparison" / "comparison_report.md").write_text(
        f"# Cross-run comparison report\n\nRuns in {root}/runs.\n", encoding="utf-8")
    (root / "members.json").write_text(json.dumps({
        "members": [{"run": "a", "dir": f"{root}/runs/a"},
                    {"run": "b", "dir": f"{root}/runs/b"}]}), encoding="utf-8")
    (root / "resolved_config.yml").write_text(
        f"systems:\n- system: 1ABC\nsweep:\n  simulation.random_seed: [1, 2, 3]\n"
        f"output: {root}\n",
        encoding="utf-8")
    runs = [{"run_id": run, "status": "ok", "output_dir": f"{root}/runs/{run}",
             "output_dir_relative": f"runs/{run}"} for run in ("a", "b", "c")]
    (root / "batch_manifest.json").write_text(json.dumps({
        "tool": "FastMDXplora", "kind": "batch", "version": "2.6.0",
        "config": "/home/someone/configs/seeds.yml", "output_dir": str(root),
        "n_runs": 3, "systems": [{"id": "1ABC", "system": "1ABC"}],
        "planned": [{"run_id": "a"}, {"run_id": "b"}, {"run_id": "c"}], "runs": runs}),
        encoding="utf-8")
    (root / "agent" / "conversations").mkdir(parents=True)
    (root / "agent" / "conversations" / "c0.json").write_text("{}", encoding="utf-8")
    return root, elsewhere


def _shared(tmp_path: Path, **options) -> tuple[Path, dict, Path, Path]:
    root, elsewhere = _campaign(tmp_path)
    out = tmp_path / "out" / "campaign.zip"
    return out, share_study(root, out, authors=AUTHORS, **options), root, elsewhere


def _by_id(crate: dict) -> dict:
    return {e["@id"]: e for e in crate["@graph"]}


def test_a_study_of_several_runs_is_packed_whole(tmp_path):
    out, said, _, _ = _shared(tmp_path)
    names = set(zipfile.ZipFile(out).namelist())
    assert {"batch_manifest.json", "resolved_config.yml", "members.json",
            "comparison/comparison_report.md", METADATA,
            "runs/a/manifest.json", "runs/a/resolved_config.yml",
            "runs/a/simulation/production.dcd", "runs/a/simulation/trajectory_topology.pdb",
            "runs/a/analysis/rmsd/rmsd.dat", "runs/b/report/report.md"} <= names
    # The preparation kept outside the study is carried once, whole enough
    # to be found and started from; the study's own shared one alike.
    assert {"prepared/preparation/setup/system.xml", "prepared/preparation/setup/state.xml",
            "prepared/preparation/setup/topology.pdb",
            "prepared/preparation/setup/setup_parameters.json",
            "prepared/preparation/setup/ligands/BEN.sdf",
            "shared_setup/setup/system.xml", "shared_setup/setup/state.xml"} <= names
    for never in ("runs/a/agent/conversations/c1.json", "agent/conversations/c0.json",
                  "runs/a/.fastmdxplora_run.json", "runs/a/study_tags.json",
                  "prepared/preparation/manifest.json"):
        assert never not in names
    crate = said["crate"]
    root = _by_id(crate)["./"]
    assert {"@id": PROFILE} in root["conformsTo"] and PROFILE.endswith("/1.1")
    assert root["fmx:runs"] == [{"@id": "runs/a/"}, {"@id": "runs/b/"}, {"@id": "runs/c/"}]
    assert "in 3 runs" in root["name"]
    # Two runs of one preparation, named two ways: carried once.
    assert not any(name.startswith("prepared/") and not name.startswith(
        "prepared/preparation/") for name in names)
    listed = check(crate)
    assert set(listed) == names - {METADATA}


def test_each_run_is_a_dataset_of_its_files_and_phases(tmp_path):
    _, said, _, _ = _shared(tmp_path)
    by_id = _by_id(said["crate"])
    member = by_id["runs/a/"]
    assert member["@type"] == "Dataset" and member["fmx:system"] == "1ABC"
    assert member["mainEntity"] == {"@id": "runs/a/resolved_config.yml"}
    assert {"@id": "runs/a/simulation/production.dcd"} in member["hasPart"]
    assert all(part["@id"].startswith("runs/a/") for part in member["hasPart"])
    assert by_id["runs/a/simulation/production.dcd"]["fmx:run"] == "a"
    assert by_id["runs/a/manifest.json"]["fmx:role"] == "manifest"
    assert by_id["batch_manifest.json"]["fmx:role"] == "manifest"
    simulation = by_id["#a/simulation"]
    assert simulation["@type"] == "CreateAction" and simulation["name"] == \
        "Simulation of the run a"
    # Each run's simulation read the prepared system it simulated from.
    assert {"@id": "prepared/preparation/setup/topology.pdb"} in simulation["object"]
    assert {"@id": "shared_setup/setup/topology.pdb"} in by_id["#b/simulation"]["object"]
    assert {"@id": "runs/a/simulation/production.dcd"} in simulation["result"]
    assert {"@id": "#a/simulation"} in member["mentions"]
    carried = by_id["prepared/preparation/"]
    assert carried["fmx:usedBy"] == [{"@id": "runs/a/"}, {"@id": "runs/c/"}]


def test_no_folder_or_computer_of_its_own_is_left(tmp_path):
    out, said, root, elsewhere = _shared(tmp_path)
    with zipfile.ZipFile(out) as archive:
        texts = {name: archive.read(name).decode("utf-8") for name in archive.namelist()
                 if name.endswith((".json", ".yml", ".md", ".sdf"))}
    for name, text in texts.items():
        for left in (str(tmp_path), COMPUTER, "prepbox3", socket.gethostname()):
            assert left not in text, (name, left)
    reference = json.loads(texts["runs/a/simulation/simulation_parameters.json"])
    assert reference["prepared_system"]["resolved"] == f"{STUDY}/prepared/preparation/setup"
    assert reference["parameters"]["setup_from"] == f"{STUDY}/prepared/preparation"
    named_alone = json.loads(texts["runs/c/simulation/simulation_parameters.json"])
    assert named_alone["parameters"]["setup_from"] == f"{STUDY}/prepared/preparation/setup"
    assert "__FASTMDX_HOME__/configs/seeds.yml" in texts["batch_manifest.json"]


def test_it_opens_as_the_study_it_was(tmp_path):
    from fastmdxplora.gui.workspace import study_state
    from fastmdxplora.simulation.pipeline import setup_records_of

    out, said, root, elsewhere = _shared(tmp_path)
    opened = open_shared(str(out), into=tmp_path / "elsewhere")
    batch = json.loads((opened / "batch_manifest.json").read_text(encoding="utf-8"))
    assert batch["output_dir"] == str(opened.resolve())
    # Each run finds the prepared system it simulated from inside the study
    # it is now, the one carried in and the study's own, by its content.
    assert setup_records_of(opened / "runs" / "a") == \
        (opened / "prepared" / "preparation" / "setup").resolve()
    assert setup_records_of(opened / "runs" / "b") == \
        (opened / "shared_setup" / "setup").resolve()
    assert setup_records_of(opened / "runs" / "c") == \
        (opened / "prepared" / "preparation" / "setup").resolve()
    config = (opened / "runs" / "a" / "resolved_config.yml").read_text(encoding="utf-8")
    assert f"setup_from: {opened.resolve()}/prepared/preparation" in config
    assert study_state(opened) == "completed"
    assert json.loads((opened / "shared_from.json").read_text())["files_checked"] == \
        said["files"]


def test_it_packs_to_the_same_bytes_twice(tmp_path):
    root, _ = _campaign(tmp_path)
    first = share_study(root, tmp_path / "a.zip", authors=AUTHORS)
    second = share_study(root, tmp_path / "b.zip", authors=AUTHORS)
    assert first["sha256"] == second["sha256"]


def test_a_run_that_did_not_finish_keeps_the_study_back(tmp_path):
    root, _ = _campaign(tmp_path, b_status="error")
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "shared.zip", authors=AUTHORS)
    assert refused.value.refusal.code == "environment.share.not_finished"
    assert str(refused.value).startswith(f"The run b of the study in {root} did not finish "
                                         "(analysis); share the study once every run has, "
                                         f"after `fastmdx resume {root}`.")
    assert not (tmp_path / "shared.zip").exists()


@pytest.mark.parametrize("planned", [[], [{"run_id": "../x"}]])
def test_a_record_of_runs_that_names_no_run_here_is_refused(tmp_path, planned):
    root, _ = _campaign(tmp_path)
    batch = json.loads((root / "batch_manifest.json").read_text())
    batch["planned"], batch["runs"] = planned, []
    (root / "batch_manifest.json").write_text(json.dumps(batch))
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "shared.zip", authors=AUTHORS)
    assert refused.value.refusal.code == "environment.share.not_a_study"


def test_a_run_s_conversation_in_an_archive_is_never_opened(tmp_path):
    """Packing never writes one; an archive holding one is not opened."""
    out, _, _, _ = _shared(tmp_path)
    crafted = tmp_path / "crafted.zip"
    with zipfile.ZipFile(out) as source, zipfile.ZipFile(crafted, "w") as sink:
        for info in source.infolist():
            sink.writestr(info, source.read(info))
        sink.writestr("runs/a/agent/conversations/c9.json", "{}")
    with pytest.raises(ShareRefused) as refused:
        open_shared(str(crafted), into=tmp_path / "elsewhere")
    assert refused.value.refusal.code == "environment.share.unsafe"
    assert not any((tmp_path / "elsewhere").iterdir())


def _of_one_run(root: Path, run: str = "a") -> None:
    """``root`` made a study of the one run ``run``."""
    (root / "resolved_config.yml").write_text(f"output: {root}\n", encoding="utf-8")
    (root / "batch_manifest.json").write_text(json.dumps({
        "version": "2.6.0", "output_dir": str(root), "planned": [{"run_id": run}],
        "runs": [{"run_id": run, "status": "ok"}]}), encoding="utf-8")


def test_a_preparation_moved_with_the_study_is_found_as_a_run_finds_it(tmp_path):
    """Review of round 1: a study fetched back or copied with its
    preparation has a stale recorded path, and the run finds its system
    relative to itself; packing looked only where the path said, and carried
    nothing."""
    import shutil

    from fastmdxplora.simulation.pipeline import setup_records_of

    root, elsewhere = _campaign(tmp_path)
    moved = tmp_path / "moved"
    moved.mkdir()
    shutil.move(str(root), moved / "campaign")
    shutil.move(str(elsewhere), moved / "preparation")
    root = moved / "campaign"
    assert setup_records_of(root / "runs" / "a") == (moved / "preparation" / "setup").resolve()
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    names = set(zipfile.ZipFile(out).namelist())
    assert "prepared/preparation/setup/system.xml" in names, sorted(
        n for n in names if "prepared" in n)
    # The study's own, found where it is now, not where it ran.
    assert {"shared_setup/setup/system.xml", "shared_setup/setup/state.xml"} <= names


def test_a_run_whose_preparation_is_not_found_keeps_the_study_back(tmp_path):
    """Review of round 1: shared without it, the run was opened with no setup
    record, so its report said setup never ran."""
    import shutil

    root, elsewhere = _campaign(tmp_path)
    shutil.rmtree(elsewhere)
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert refused.value.refusal.code == "environment.share.not_finished"
    assert "The run a " in str(refused.value) and "not found here" in str(refused.value)
    assert not (tmp_path / "o.zip").exists()


def test_a_window_seeded_from_a_pull_outside_carries_the_pull_s_preparation(tmp_path):
    """Review of round 1: a window's seed is in the study and keeps no setup
    record; the preparation it was taken from was outside and not carried."""
    import shutil

    from fastmdxplora.simulation.pipeline import setup_records_of
    from fastmdxplora.simulation.seeding import record_the_preparation

    pull = _prepared(tmp_path / "pull" / "setup")
    root = tmp_path / "umbrella"
    seed = root / "seeds" / "window-00"
    seed.mkdir(parents=True)
    for name in ("system.xml", "state.xml", "topology.pdb"):
        shutil.copy2(pull / name, seed / name)
    record_the_preparation(seed, pull, frame=3, centre=0.5)
    _run(root, "w0", seed, str(seed))
    _of_one_run(root, "w0")
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    names = set(zipfile.ZipFile(out).namelist())
    assert {"seeds/window-00/system.xml", "seeds/window-00/seeded_from.json",
            "prepared/pull/setup/setup_parameters.json"} <= names
    seeded = json.loads(zipfile.ZipFile(out).read("seeds/window-00/seeded_from.json"))
    assert seeded["relative_to_run"] == "../../prepared/pull/setup"
    shutil.rmtree(tmp_path / "pull")
    opened = open_shared(str(out), into=tmp_path / "x")
    assert setup_records_of(opened / "runs" / "w0") == \
        (opened / "prepared" / "pull" / "setup").resolve()


def test_the_way_from_a_run_to_its_preparation_names_no_home_folder(tmp_path, monkeypatch):
    """Review of round 1: a study in scratch simulated from a preparation in
    the home folder recorded `../../../home/<name>/prep/setup`, which kept
    the person's name and no check found."""
    home = tmp_path / "home" / "rosalind"
    monkeypatch.setenv("HOME", str(home))
    prep = _prepared(home / "prep" / "setup").parent
    root = tmp_path / "scratch" / "job" / "campaign"
    root.mkdir(parents=True)
    _run(root, "a", prep / "setup", str(prep) + "/")
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    with zipfile.ZipFile(out) as archive:
        texts = {n: archive.read(n).decode("utf-8") for n in archive.namelist()
                 if n.endswith((".json", ".yml", ".sdf"))}
    assert not any("rosalind" in text for text in texts.values())
    reference = json.loads(texts["runs/a/simulation/simulation_parameters.json"])
    assert reference["prepared_system"]["relative_to_run"] == "../../prepared/prep/setup"
    # Named with a slash at its end, as typed.
    assert reference["parameters"]["setup_from"] == f"{STUDY}/prepared/prep/"


def test_a_prepared_system_is_packed_byte_for_byte(tmp_path):
    """Review of round 1: `system.xml` was read for paths and names like any
    text, so a word of it that was the computer's name changed its bytes,
    and the run opened elsewhere no longer found its system by them."""
    from fastmdxplora.simulation.pipeline import prepared_system_reference, setup_records_of

    root, elsewhere = _campaign(tmp_path)
    system = elsewhere / "setup" / "system.xml"
    system.write_text(f"<System name='{COMPUTER}'/>", encoding="utf-8")
    for name in ("a", "c"):
        run = root / "runs" / name
        record = json.loads((run / "simulation" / "simulation_parameters.json").read_text())
        record["prepared_system"] = prepared_system_reference(
            elsewhere / "setup", run, record["prepared_system"]["given"])
        (run / "simulation" / "simulation_parameters.json").write_text(json.dumps(record))
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    assert zipfile.ZipFile(out).read("prepared/preparation/setup/system.xml") == \
        system.read_bytes()
    opened = open_shared(str(out), into=tmp_path / "x")
    assert setup_records_of(opened / "runs" / "a") == \
        (opened / "prepared" / "preparation" / "setup").resolve()


def test_only_the_runs_of_its_record_go_and_its_own_folders_are_kept(tmp_path):
    """Review of round 1: a run's folder the record of runs does not name (an
    earlier attempt) went in with no dataset of its own; a run named with a
    backslash was taken for a folder; and a study's own `prepared/` folder
    could be written over by a preparation carried in."""
    root, _ = _campaign(tmp_path)
    (root / "runs" / "old").mkdir()
    (root / "runs" / "old" / "notes.md").write_text("an earlier attempt\n")
    (root / "prepared" / "preparation").mkdir(parents=True)
    (root / "prepared" / "preparation" / "mine.txt").write_text("the study's own\n")
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    names = set(zipfile.ZipFile(out).namelist())
    assert not any(n.startswith("runs/old/") for n in names)
    assert "prepared/preparation/mine.txt" in names
    assert "prepared/preparation-2/setup/system.xml" in names
    batch = json.loads((root / "batch_manifest.json").read_text())
    batch["planned"] = [{"run_id": "a\\..\\..\\x"}]
    (root / "batch_manifest.json").write_text(json.dumps(batch))
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "p.zip", authors=AUTHORS)
    assert refused.value.refusal.code == "environment.share.not_a_study"


def test_a_run_still_going_is_named(tmp_path):
    """Review of round 1: the refusal said runs were still going, not which."""
    root, _ = _campaign(tmp_path)
    (root / "runs" / "b" / "simulation" / "live_status.json").write_text(
        json.dumps({"status": "running"}), encoding="utf-8")
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert str(refused.value) == (f"The run b of the study in {root} is still going; share "
                                  "the study once every run has finished.")


def _texts_of(out: Path) -> dict[str, str]:
    with zipfile.ZipFile(out) as archive:
        return {n: archive.read(n).decode("utf-8") for n in archive.namelist()
                if n.endswith((".json", ".yml", ".md", ".sdf")) and n != METADATA}


def test_a_study_moved_with_its_preparation_keeps_no_path_of_where_it_ran(tmp_path):
    """Review of round 2: moved together after it ran, its records still
    named the folders of the computer it ran on, and its Config named a
    folder that was no longer there."""
    import shutil

    ran = tmp_path / "cluster" / "scratch" / "kunle" / "job7"
    ran.mkdir(parents=True)
    root, elsewhere = _campaign(ran)
    here = tmp_path / "fetched"
    shutil.move(str(ran), here)
    out = tmp_path / "o.zip"
    share_study(here / "campaign", out, authors=AUTHORS)
    texts = _texts_of(out)
    for name, text in texts.items():
        assert str(ran) not in text, name
    opened = open_shared(str(out), into=tmp_path / "x")
    config = (opened / "runs" / "a" / "resolved_config.yml").read_text(encoding="utf-8")
    assert f"setup_from: {opened.resolve()}/prepared/preparation" in config


@pytest.mark.parametrize("where", ["a run's folder not in the record", "a dot folder"])
def test_a_system_in_a_folder_the_archive_leaves_out_is_carried(tmp_path, where):
    """Review of round 2: a run's system inside the study, in a folder the
    archive leaves out, was dropped and the run opened with none; refused
    then (round 3), the refusal's remedy led to another refusal."""
    import shutil

    from fastmdxplora.simulation.pipeline import setup_records_of

    root = tmp_path / "campaign"
    folder = {"a run's folder not in the record": root / "runs" / "old" / "setup",
              "a dot folder": root / ".prep" / "setup"}[where]
    _prepared(folder)
    _run(root, "a", folder, str(folder.parent))
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    names = set(zipfile.ZipFile(out).namelist())
    carried = {"a run's folder not in the record": "prepared/old/setup/system.xml",
               "a dot folder": "prepared/prep/setup/system.xml"}[where]
    assert carried in names
    assert not any(n.startswith(("runs/old/", ".prep/")) for n in names)
    shutil.rmtree(folder.parent)
    opened = open_shared(str(out), into=tmp_path / "x")
    assert setup_records_of(opened / "runs" / "a") == \
        (opened / PurePosixPath(carried).parent).resolve()


def test_a_system_inside_the_study_kept_as_links_is_refused(tmp_path):
    """Review of round 2: a run's system kept as links was dropped."""
    root = tmp_path / "campaign"
    real = _prepared(tmp_path / "kept")
    folder = root / "linked"
    folder.mkdir(parents=True)
    for name in ("system.xml", "state.xml", "topology.pdb", "setup_parameters.json"):
        (folder / name).symlink_to(real / name)
    _run(root, "a", folder, str(folder))
    _of_one_run(root)
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert refused.value.refusal.code == "environment.share.not_finished"
    # Every link at once, not one a try.
    assert ("its setup_parameters.json, state.xml, system.xml, topology.pdb are links"
            in str(refused.value))
    assert not (tmp_path / "o.zip").exists()


def test_each_run_s_own_recorded_folder_is_said_as_its_place(tmp_path):
    """Review of round 3: a run whose own record gives another folder than
    the record of runs (carried on elsewhere) kept that folder's path."""
    root, _ = _campaign(tmp_path)
    manifest = root / "runs" / "a" / "manifest.json"
    found = json.loads(manifest.read_text())
    found["output_dir"] = "/somewhere/else/ws/clean-own"
    manifest.write_text(json.dumps(found))
    out = tmp_path / "o.zip"
    said = share_study(root, out, authors=AUTHORS)
    texts = _texts_of(out)
    assert "/somewhere/else" not in texts["runs/a/manifest.json"]
    assert json.loads(texts["runs/a/manifest.json"])["output_dir"] == f"{STUDY}/runs/a"
    opened = open_shared(str(out), into=tmp_path / "x")
    assert opened.name == "1abc-in-3-runs"
    assert _by_id(said["crate"])["runs/a/"]["fmx:system"] == "1ABC"


def test_a_home_reached_through_a_link_is_said_as_home(tmp_path, monkeypatch):
    """Review of round 2: only the home folder's path as named was said as
    home, not where a link to it leads, which records written from there
    give."""
    real_home = tmp_path / "gpfs" / "u7x3q"
    real_home.mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (tmp_path / "home" / "rosalind").symlink_to(real_home)
    monkeypatch.setenv("HOME", str(tmp_path / "home" / "rosalind"))
    prep = _prepared(real_home / "prep" / "setup").parent
    (prep / "setup" / "setup_parameters.json").write_text(json.dumps({
        "structure": f"{real_home}/structures/x.pdb"}), encoding="utf-8")
    root = tmp_path / "scratch" / "campaign"
    root.mkdir(parents=True)
    _run(root, "a", prep / "setup", str(prep))
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    texts = _texts_of(out)
    assert not any("u7x3q" in text for text in texts.values())
    assert "__FASTMDX_HOME__/structures/x.pdb" in \
        texts["prepared/prep/setup/setup_parameters.json"]


def test_a_preparation_named_through_a_link_is_said_as_carried(tmp_path):
    """Review of round 2: a path under a preparation named through a link
    was left as written."""
    real = _prepared(tmp_path / "realprep" / "setup").parent
    (tmp_path / "linkprep").symlink_to(real)
    root = tmp_path / "campaign"
    run = _run(root, "a", real / "setup", str(tmp_path / "linkprep"))
    with (run / "resolved_config.yml").open("a", encoding="utf-8") as config:
        config.write(f"report:\n  structure: {tmp_path}/linkprep/setup/topology.pdb\n")
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    config = _texts_of(out)["runs/a/resolved_config.yml"]
    assert str(tmp_path) not in config
    assert f"structure: {STUDY}/prepared/realprep/setup/topology.pdb" in config


def test_one_study_s_two_systems_are_carried_in_one_folder(tmp_path):
    """Review of round 2: a study holding two systems the runs used was
    carried twice, and its name said as one of them only."""
    from fastmdxplora.simulation.pipeline import setup_records_of

    other = tmp_path / "other"
    _prepared(other / "setup")
    _prepared(other / "shared_setup" / "setup")
    root = tmp_path / "campaign"
    _run(root, "a", other / "shared_setup" / "setup", str(other))
    _run(root, "c", other / "setup", str(other / "setup"))
    _of_one_run(root)
    batch = json.loads((root / "batch_manifest.json").read_text())
    batch["planned"] = [{"run_id": "a"}, {"run_id": "c"}]
    (root / "batch_manifest.json").write_text(json.dumps(batch))
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    names = set(zipfile.ZipFile(out).namelist())
    assert {"prepared/other/setup/system.xml",
            "prepared/other/shared_setup/setup/system.xml"} <= names
    assert not any(n.startswith("prepared/other-2") for n in names)
    texts = _texts_of(out)
    assert f"setup_from: {STUDY}/prepared/other\n" in texts["runs/a/resolved_config.yml"]
    opened = open_shared(str(out), into=tmp_path / "x")
    assert setup_records_of(opened / "runs" / "a") == \
        (opened / "prepared" / "other" / "shared_setup" / "setup").resolve()
    assert setup_records_of(opened / "runs" / "c") == \
        (opened / "prepared" / "other" / "setup").resolve()


def test_every_computer_s_name_its_records_give_is_said_as_one(tmp_path):
    """Review of round 2: an umbrella's pull run on another node kept its
    name, and a name too short to be told from a word (`vm`) kept in a
    record's own field for it."""
    root, _ = _campaign(tmp_path)
    (root / "seed_pull").mkdir()
    (root / "seed_pull" / "manifest.json").write_text(json.dumps({"phases": [
        {"name": "simulation", "status": "ok", "produced_by": {"host": "pullnode9"}}]}))
    (root / "seed_pull" / "pull.log").write_text("pulled on pullnode9\n")
    (root / "runs" / "a" / "simulation" / "live_status.json").write_text(
        json.dumps({"status": "completed", "host": "vm"}))
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    with zipfile.ZipFile(out) as archive:
        assert b"pullnode9" not in archive.read("seed_pull/manifest.json")
    live = json.loads(_texts_of(out)["runs/a/simulation/live_status.json"])
    assert live["host"] == "__FASTMDX_HOST__"


def test_a_run_never_started_or_failed_is_said_with_the_command(tmp_path):
    """Review of round 2: the remedy gave `fastmdx resume` with no folder."""
    import shutil

    root, _ = _campaign(tmp_path)
    shutil.rmtree(root / "runs" / "c")
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert str(refused.value) == (f"The run c of the study in {root} never ran; run it with "
                                  f"`fastmdx resume {root}`, then share the study.")
    batch = json.loads((root / "batch_manifest.json").read_text())
    batch["planned"] = batch["planned"][:2]
    batch["runs"][1]["status"] = "error"
    batch["runs"] = batch["runs"]
    (root / "batch_manifest.json").write_text(json.dumps(batch))
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert str(refused.value) == (f"The run b of the study in {root} ended with an error; "
                                  "share the study once every run has finished, after "
                                  f"`fastmdx resume {root}`.")
    batch["runs"] = 5
    (root / "batch_manifest.json").write_text(json.dumps(batch))
    share_study(root, tmp_path / "p.zip", authors=AUTHORS)  # read without results


def test_a_study_of_one_run_is_said_so(tmp_path):
    root = tmp_path / "campaign"
    _run(root, "a", _prepared(root / "shared_setup" / "setup"),
         str(root / "shared_setup" / "setup"))
    _of_one_run(root)
    said = share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert _by_id(said["crate"])["./"]["name"].endswith("in one run")


def test_a_preparation_outside_holding_a_link_is_refused(tmp_path):
    """Review of round 3: a system outside the study whose files are links
    was carried without them, and its runs opened with none."""
    kept = _prepared(tmp_path / "kept")
    prep = tmp_path / "prep" / "setup"
    prep.mkdir(parents=True)
    for name in ("system.xml", "state.xml", "topology.pdb", "setup_parameters.json"):
        (prep / name).symlink_to(kept / name)
    root = tmp_path / "campaign"
    _run(root, "a", prep, str(prep.parent))
    _of_one_run(root)
    with pytest.raises(ShareRefused) as refused:
        share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    assert refused.value.refusal.code == "environment.share.not_finished"
    assert "holds 4 links" in str(refused.value)
    for name in ("system.xml", "state.xml", "topology.pdb", "setup_parameters.json"):
        assert str(prep / name) in str(refused.value)


def test_a_folder_carried_from_is_one_dataset(tmp_path):
    """Review of round 3: a study's two systems gave two entities of one
    `@id`, each with half the runs that used it."""
    other = tmp_path / "other"
    _prepared(other / "setup")
    _prepared(other / "shared_setup" / "setup")
    root = tmp_path / "campaign"
    _run(root, "a", other / "shared_setup" / "setup", str(other))
    _run(root, "c", other / "setup", str(other / "setup"))
    _of_one_run(root)
    batch = json.loads((root / "batch_manifest.json").read_text())
    batch["planned"] = [{"run_id": "a"}, {"run_id": "c"}]
    (root / "batch_manifest.json").write_text(json.dumps(batch))
    said = share_study(root, tmp_path / "o.zip", authors=AUTHORS)
    ids = [e["@id"] for e in said["crate"]["@graph"]]
    assert len(ids) == len(set(ids))
    carried = _by_id(said["crate"])["prepared/other/"]
    assert sorted(u["@id"] for u in carried["fmx:usedBy"]) == ["runs/a/", "runs/c/"]


def test_a_name_alone_keeps_what_follows_it(tmp_path):
    """Review of round 3: the slash after a study folder's name was taken
    with it before a character that ends a path (`/x/prep/*`)."""
    from fastmdxplora.sharing.pack import _Scrubber

    scrub = _Scrubber(["/data/camp"], [], {"/data/prep": ("prepared/prep", True)})
    assert scrub("copied /data/prep/* here")[0] == f"copied {STUDY}/prepared/prep/* here"
    assert scrub("setup_from: /data/prep\n")[0] == f"setup_from: {STUDY}/prepared/prep\n"
    assert scrub("/data/prep/manifest.json")[0] == "/data/prep/manifest.json"


def test_a_preparation_moved_into_the_study_keeps_no_path_of_where_it_was(tmp_path):
    """Review of round 3: a study fetched back with its preparation copied
    into it kept the path the preparation had where it ran."""
    import shutil

    ran = tmp_path / "cluster" / "job7"
    ran.mkdir(parents=True)
    _campaign(ran)
    here = tmp_path / "fetched"
    shutil.move(str(ran), here)
    root = here / "campaign"
    shutil.move(str(here / "preparation"), root / "preparation")
    for run in ("a", "c"):
        record = root / "runs" / run / "simulation" / "simulation_parameters.json"
        found = json.loads(record.read_text())
        found["prepared_system"]["relative_to_run"] = "../../preparation/setup"
        record.write_text(json.dumps(found))
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    for name, text in _texts_of(out).items():
        assert str(ran) not in text, name
    assert not any(n.startswith("prepared/") for n in zipfile.ZipFile(out).namelist())


def test_a_folder_named_for_the_person_is_not_carried_by_that_name(tmp_path, monkeypatch):
    """Review of round 3: `/scratch/<login>/setup` named by itself was carried
    as `prepared/<login>/`, the person's name in every path of it."""
    import getpass

    monkeypatch.setattr(getpass, "getuser", lambda: "u7x3q")
    prep = _prepared(tmp_path / "scratch" / "u7x3q" / "setup")
    root = tmp_path / "campaign"
    _run(root, "a", prep, str(prep))
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    names = zipfile.ZipFile(out).namelist()
    assert not any("u7x3q" in n for n in names)
    assert "prepared/preparation/setup/system.xml" in names


def test_the_recorded_seeds_study_is_shared_and_opened(tmp_path):
    """A study of three replicas as FastMDXplora wrote one (the records of
    tests/data/report_traced/seeds), named for its system, not its file."""
    import shutil

    from fastmdxplora.gui.exploration import runs_of_a_study
    from fastmdxplora.gui.workspace import study_state

    root = tmp_path / "seeds"
    shutil.copytree(Path(__file__).parent / "data" / "report_traced" / "seeds", root)
    said = share_study(root, tmp_path / "seeds.zip", authors=AUTHORS)
    assert _by_id(said["crate"])["./"]["name"] == "A FastMDXplora study of tri-ala, in 3 runs"
    opened = open_shared(str(tmp_path / "seeds.zip"), into=tmp_path / "x")
    assert opened.name == "tri-ala-in-3-runs"
    assert study_state(opened) == "completed"
    assert [r["state"] for r in runs_of_a_study(opened)] == ["completed"] * 3


def test_a_run_s_recorded_folder_never_says_another_run_s_paths(tmp_path):
    """Review of series 2: two runs carried on in one job folder, or a run
    recording another's folder, had one run's paths said as the other's."""
    root, _ = _campaign(tmp_path)
    for run in ("a", "b"):
        manifest = root / "runs" / run / "manifest.json"
        found = json.loads(manifest.read_text())
        found["output_dir"] = "/scratch/job"
        manifest.write_text(json.dumps(found))
        (root / "runs" / run / "analysis" / "analysis_manifest.json").write_text(
            json.dumps({"trajectory_input": "/scratch/job/simulation/production.dcd"}))
    manifest = root / "runs" / "c" / "manifest.json"
    found = json.loads(manifest.read_text())
    found["output_dir"] = str(root / "runs" / "b")
    manifest.write_text(json.dumps(found))
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    texts = _texts_of(out)
    for run in ("a", "b"):
        read = json.loads(texts[f"runs/{run}/analysis/analysis_manifest.json"])
        assert not read["trajectory_input"].startswith(f"{STUDY}/runs/"), run
    assert f"{STUDY}/runs/c" not in texts["runs/b/manifest.json"]


def test_a_preparation_s_name_is_said_as_carried_before_a_run_s_folder(tmp_path):
    """Review of series 2: a run that ran in the folder a preparation is
    named by had that name said as the run's folder."""
    root = tmp_path / "campaign"
    job = tmp_path / "x" / "job"
    prep = _prepared(job / "setup")
    _run(root, "a", prep, str(job))
    _of_one_run(root)
    manifest = root / "runs" / "a" / "manifest.json"
    found = json.loads(manifest.read_text())
    found["output_dir"] = str(job)
    manifest.write_text(json.dumps(found))
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    record = json.loads(_texts_of(out)["runs/a/simulation/simulation_parameters.json"])
    assert record["parameters"]["setup_from"] == f"{STUDY}/prepared/job"


def test_a_system_in_a_folder_set_aside_is_carried(tmp_path):
    """Review of series 2: a system in `previous/` was refused as a link."""
    from fastmdxplora.simulation.pipeline import setup_records_of

    root = tmp_path / "campaign"
    folder = root / "previous" / "prep" / "setup"
    _prepared(folder)
    _run(root, "a", folder, str(folder.parent))
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    assert "prepared/prep/setup/system.xml" in zipfile.ZipFile(out).namelist()
    opened = open_shared(str(out), into=tmp_path / "x")
    assert setup_records_of(opened / "runs" / "a") == \
        (opened / "prepared" / "prep" / "setup").resolve()


def test_a_preparation_named_from_home_is_said_as_carried(tmp_path, monkeypatch):
    """Review of series 2: `setup_from: ~/prep` stayed as typed, so the
    Config run again elsewhere read the reader's own `~/prep`."""
    home = tmp_path / "home" / "rosa"
    home.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    prep = _prepared(home / "prep" / "setup").parent
    root = tmp_path / "scratch" / "campaign"
    root.mkdir(parents=True)
    _run(root, "a", prep / "setup", "~/prep")
    _of_one_run(root)
    out = tmp_path / "o.zip"
    share_study(root, out, authors=AUTHORS)
    config = _texts_of(out)["runs/a/resolved_config.yml"]
    assert "~/prep" not in config and f"setup_from: {STUDY}/prepared/prep" in config


def test_a_name_another_computer_s_home_holds_is_said_as_a_refusal_or_found(tmp_path):
    """Review of series 2: `~otheruser/...`, read here, raised an error that
    was no refusal."""
    import shutil

    from fastmdxplora.simulation.pipeline import setup_records_of

    ran = tmp_path / "cluster"
    root = ran / "campaign"
    prep = _prepared(root / "shared_setup" / "setup")
    _run(root, "a", prep, "~nosuchuser42/campaign/shared_setup/setup")
    _of_one_run(root)
    here = tmp_path / "here"
    shutil.move(str(ran), here)
    out = tmp_path / "o.zip"
    share_study(here / "campaign", out, authors=AUTHORS)
    assert not any("nosuchuser42" in text for text in _texts_of(out).values())
    opened = open_shared(str(out), into=tmp_path / "x")
    assert setup_records_of(opened / "runs" / "a") == \
        (opened / "shared_setup" / "setup").resolve()


def test_a_copied_preparation_keeps_no_path_of_where_it_was_made(tmp_path):
    """Review of series 2: a preparation made in one folder and copied to
    another kept, in its own setup record and ligands, the folder it was
    made in."""
    import shutil

    made = tmp_path / "cluster" / "scratch" / "kaina" / "prep7"
    _prepared(made / "setup")
    copy = tmp_path / "preps" / "prep7"
    shutil.copytree(made, copy)
    root = tmp_path / "campaign"
    _run(root, "a", copy / "setup", str(copy))
    _of_one_run(root)
    out = tmp_path / "o.zip"
    printed: list[str] = []
    said = share_study(root, out, authors=AUTHORS, said=printed.append)
    for name, text in _texts_of(out).items():
        assert str(made) not in text, name
    assert json.loads(_texts_of(out)["prepared/prep7/setup/setup_parameters.json"])[
        "output_dir"] == f"{STUDY}/prepared/prep7/setup"
    assert printed[-1].endswith(f"kB), SHA-256 {said['sha256']}.")
