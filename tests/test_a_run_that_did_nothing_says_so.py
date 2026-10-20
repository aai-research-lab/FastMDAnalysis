"""A study that could not do what it was asked says so, in its exit code and
before it runs where that can be known.

Seen by the weekly test of 10-20: a study whose RMSD selection MDTraj could
not parse said "rmsd x error ... Cannot use literals as truth", "analysis
complete" and exited 0; a sequence as the system, which the minimal template
offered, was recorded as a setup that went "ok", with nothing prepared, and
exited 0; `--dry-run` passed a structure file that is not there, a system
that is no PDB identifier, and a selection that does not parse, each found
only when setup or the analysis reached it; a 4 fs timestep with nothing
constrained was worth knowing in the GUI's Check it and said by no dry run;
and `fastmdx config -f` into a folder that does not exist was a traceback.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from test_report_wiring import _make_traj_files  # noqa: E402

from fastmdxplora.cli.main import main  # noqa: E402

BAD = "protein and nam CA"


def _config(path: Path, **data) -> Path:
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


def _analysis_study(tmp_path: Path, analysis: dict, n_frames: int = 30) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    root = tmp_path / "study"
    _make_traj_files(root, n_residues=4, n_frames=n_frames)
    return _config(tmp_path / "study.yml",
                   systems=[{"id": "s1", "system": str(root / "simulation" / "topology.pdb")}],
                   output=str(root), include_phase=["analysis"], analysis=analysis)


NONE = "resname ZZZ"


def test_an_analysis_asked_for_that_fails_is_a_failed_run(tmp_path, capsys):
    # A selection that parses and matches nothing is found only as the
    # analysis reads the trajectory.
    config = _analysis_study(tmp_path / "bad", {"include": ["rmsd", "rg"],
                                                "options": {"rmsd": {"selection": NONE}}})
    assert main(["explore", "--config", str(config)]) == 1
    said = capsys.readouterr().err
    assert "the study asked for failed" in said and "rmsd" in said
    assert "matched zero atoms" in said
    # The same study with a selection that matches succeeds.
    config = _analysis_study(tmp_path / "good", {"include": ["rmsd", "rg"]})
    assert main(["explore", "--config", str(config)]) == 0
    assert "asked for failed" not in capsys.readouterr().err
    # An analysis of the default plan that frames this few cannot support
    # is said in its own row, and the run is not failed for it.
    config = _analysis_study(tmp_path / "short", {}, n_frames=3)
    assert main(["explore", "--config", str(config)]) == 0
    record = json.loads((tmp_path / "short" / "study" / "analysis"
                         / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert any(r["status"] == "error" for r in record["results"].values())


def test_only_an_analysis_named_is_asked_for(tmp_path, capsys):
    """A selection for the analyses (`select_atoms`) asks for none of them:
    a short run's clustering, of the default plan, failing on too few frames
    failed the run once `select_atoms` was set."""
    config = _analysis_study(tmp_path / "selected", {"select_atoms": "name CA"}, n_frames=3)
    assert main(["explore", "--config", str(config)]) == 0
    assert "asked for failed" not in capsys.readouterr().err
    record = json.loads((tmp_path / "selected" / "study" / "analysis"
                         / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert record["results"]["cluster"]["status"] == "error"


def test_too_few_frames_never_fails_the_run(tmp_path, capsys):
    """Named, and the run too short for it: the run did what it was asked."""
    config = _analysis_study(tmp_path / "named", {"include": ["rmsd", "cluster"]}, n_frames=3)
    assert main(["explore", "--config", str(config)]) == 0
    assert "asked for failed" not in capsys.readouterr().err


def test_an_analysis_a_system_names_is_asked_for(tmp_path, capsys):
    from fastmdxplora.cli.main import _analyses_asked_for_that_failed

    class Run:
        run_id = "s1"

        def __init__(self, folder):
            self.output_dir = folder

    folder = tmp_path / "run"
    (folder / "analysis").mkdir(parents=True)
    (folder / "analysis" / "analysis_manifest.json").write_text(json.dumps(
        {"results": {"rmsd": {"status": "error", "message": "rmsd: matched zero atoms"}}}))
    config = {"systems": [{"id": "s1", "system": "x.pdb",
                           "analysis": {"include": ["rmsd"]}}]}
    assert _analyses_asked_for_that_failed(config, [Run(folder)]) == [
        "rmsd (matched zero atoms)"]
    assert _analyses_asked_for_that_failed({"systems": [{"system": "x.pdb"}]},
                                           [Run(folder)]) == []


def test_a_bad_selection_is_refused_before_a_run_writes_anything(tmp_path, capsys):
    config = _analysis_study(tmp_path, {"select_atoms": BAD})
    before = sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*"))
    assert main(["explore", "--config", str(config)]) == 2
    assert "is not an MDTraj selection" in capsys.readouterr().err
    assert sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*")) == before


def test_a_sequence_is_a_setup_that_failed(tmp_path, capsys):
    from fastmdxplora import FastMDXplora

    out = tmp_path / "seq"
    results = FastMDXplora(system="MKTAYIAKQRQISFVKSHFSRQ", output_dir=out).explore(
        include_phase=["setup"])
    assert results[0].status != "ok"
    assert "structure predictor" in (results[0].message or "") + " ".join(
        str(getattr(p, "message", "") or "") for p in results[0].phases)
    assert (out / "setup" / "input.sequence").is_file()
    # Given to the command line, it is refused before anything is written.
    rc = main(["explore", "--system", "MKTAYIAKQRQISFVKSHFSRQ", "--output",
               str(tmp_path / "cli"), "--include-phase", "setup"])
    assert rc == 2
    assert "structure predictor" in capsys.readouterr().err
    assert not (tmp_path / "cli").exists()


@pytest.mark.parametrize("system, says", [
    ("not_here.pdb", "No structure at 'not_here.pdb'"),
    ("1L2YA", "Could not classify system input '1L2YA'"),
    ("MKTAYIAKQRQISFVKSHFSRQ", "reads as an amino-acid sequence"),
])
def test_a_dry_run_refuses_only_a_system_setup_cannot_read(tmp_path, monkeypatch, capsys,
                                                           system, says):
    monkeypatch.chdir(tmp_path)
    _make_traj_files(tmp_path / "made")
    here = tmp_path / "made" / "simulation" / "topology.pdb"
    # A PDB identifier and a file that is there pass, and so does any system
    # where setup does not run.
    fine = _config(tmp_path / "fine.yml", systems=[{"id": "a", "system": "1L2Y"},
                                                   {"id": "b", "system": str(here)}])
    assert main(["explore", "--config", str(fine), "--dry-run"]) == 0
    unset = _config(tmp_path / "unset.yml", systems=[{"id": "a", "system": system}],
                    include_phase=["report"])
    assert main(["explore", "--config", str(unset), "--dry-run"]) == 0
    capsys.readouterr()

    config = _config(tmp_path / "s.yml", systems=[{"id": "a", "system": system}])
    assert main(["explore", "--config", str(config), "--dry-run"]) == 2
    assert says in capsys.readouterr().err
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "fine.yml", "made", "s.yml", "unset.yml"]


@pytest.mark.parametrize("block", [
    {"analysis": {"select_atoms": BAD}},
    {"analysis": {"options": {"rmsd": {"selection": BAD}}}},
    {"simulation": {"save_selection": BAD}},
])
def test_a_dry_run_refuses_a_selection_that_does_not_parse(tmp_path, monkeypatch, capsys,
                                                           block):
    monkeypatch.chdir(tmp_path)
    config = _config(tmp_path / "s.yml", systems=[{"id": "a", "system": "1L2Y"}], **block)
    assert main(["explore", "--config", str(config), "--dry-run"]) == 2
    said = capsys.readouterr().err
    assert repr(BAD) in said and "is not an MDTraj selection" in said


def test_a_dry_run_says_what_check_it_says_is_worth_knowing(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    config = _config(tmp_path / "s.yml", systems=[{"id": "a", "system": "1L2Y"}],
                     setup={"constraints": "None"}, simulation={"timestep_fs": 4})
    assert main(["explore", "--config", str(config), "--dry-run"]) == 0
    said = capsys.readouterr().out
    assert "Worth knowing before it runs:" in said
    assert "(timestep_fs) A 4 fs timestep is longer than the" in said


def test_a_template_into_a_folder_that_is_not_there_is_one_line(tmp_path, capsys):
    rc = main(["config", "-f", str(tmp_path / "nodir" / "x.yml")])
    assert rc == 2
    said = capsys.readouterr().err
    assert said.count("\n") == 1 and "cannot write" in said
    assert "Traceback" not in said
