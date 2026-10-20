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


def test_an_analysis_asked_for_that_fails_is_a_failed_run(tmp_path, capsys):
    config = _analysis_study(tmp_path / "bad", {"include": ["rmsd", "rg"],
                                                "select_atoms": BAD})
    assert main(["explore", "--config", str(config)]) == 1
    said = capsys.readouterr().err
    assert "the study asked for failed" in said and "rmsd" in said
    assert "Cannot use literals as truth" in said
    # The same study with a selection that parses succeeds.
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


def test_a_selection_the_study_set_counts_without_include(tmp_path, capsys):
    config = _analysis_study(tmp_path, {"select_atoms": BAD})
    assert main(["explore", "--config", str(config)]) == 1
    assert "rmsd" in capsys.readouterr().err


def test_a_sequence_is_a_setup_that_failed(tmp_path, capsys):
    out = tmp_path / "seq"
    rc = main(["explore", "--system", "MKTAYIAKQRQISFVKSHFSRQ", "--output", str(out),
               "--include-phase", "setup"])
    assert rc == 1
    said = capsys.readouterr()
    assert "structure predictor" in said.out + said.err
    assert (out / "setup" / "input.sequence").is_file()


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
