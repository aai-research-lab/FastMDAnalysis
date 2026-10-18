"""A study of several runs given no output is named by its first system.

`BatchExplorer` named the folder from `self.config`, an attribute it never
sets, so a study of runs started from Python with no `output` stopped with
an AttributeError before it ran.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from fastmdxplora.batch.explorer import BatchExplorer


@pytest.fixture
def structure(tmp_path: Path) -> Path:
    path = tmp_path / "protein.pdb"
    path.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "END\n", encoding="utf-8")
    return path


def test_no_output_names_the_folder_by_the_first_system(tmp_path, structure,
                                                         monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    explorer = BatchExplorer(config_data={
        "systems": [{"id": "chignolin", "system": str(structure)},
                    {"id": "other", "system": str(structure)}],
        "sweep": {"simulation.random_seed": [1, 2]},
    })
    assert re.fullmatch(r"fastmdxplora_chignolin_study_\d{14}",
                        str(explorer.output_dir))


def test_an_output_given_is_kept(tmp_path, structure) -> None:
    explorer = BatchExplorer(config_data={
        "output": str(tmp_path / "out"),
        "systems": [{"id": "chignolin", "system": str(structure)}],
        "sweep": {"simulation.random_seed": [1, 2]},
    })
    assert explorer.output_dir == tmp_path / "out"
