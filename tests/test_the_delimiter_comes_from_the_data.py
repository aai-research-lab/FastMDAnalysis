"""The delimiter is a property of the data, not of the note above it.

`_numeric_column` sniffed the delimiter from the file's first line. That
held while a whitespace `.dat` began with a number, and stopped holding the
moment `save_data` started writing a `#` line naming the layout: the note
reads "whitespace-delimited, no column header", and the comma in it made a
whitespace file parse as comma-separated. `names` then came back wider than
the array, and the column lookup raised `IndexError`.

It surfaced for `sasa` and not for `ligand_rmsd`, because that one's
preferred column name happened to match the first comma-split token, so the
lookup landed on index 0 and the values came back correct. A test written
against `ligand_rmsd` alone would have passed.

The wording could have been changed instead. A message that must never
contain a comma is a trap for whoever edits it next, so the reader takes
its delimiter from a data row.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.validation.cross_tool import _numeric_column


class _Array(Analysis):
    name = "sasa"

    def compute(self, traj):  # pragma: no cover - not run
        raise NotImplementedError

    def plot(self, result, ax):  # pragma: no cover - not run
        raise NotImplementedError


class _Frame(Analysis):
    name = "pl_interactions"

    def compute(self, traj):  # pragma: no cover - not run
        raise NotImplementedError

    def plot(self, result, ax):  # pragma: no cover - not run
        raise NotImplementedError


class TestEveryShapeOnDisk:
    def test_whitespace_under_a_hash_note(self, tmp_path):
        """What this release writes. The note names the analysis and
        contains a comma; the data has one column and no commas at all."""
        path = _Array(output_dir=tmp_path).save_data(
            np.array([1.0, 2.0, 3.0]), tmp_path / "sasa.dat")

        assert np.allclose(_numeric_column(path, prefer="area"), [1.0, 2.0, 3.0])

    def test_the_note_is_not_mistaken_for_column_names(self, tmp_path):
        path = _Array(output_dir=tmp_path).save_data(
            np.array([[0.0, 1.0], [1.0, 2.0]]), tmp_path / "two.dat")

        assert np.allclose(_numeric_column(path, prefer="nothing"), [1.0, 2.0])

    def test_comma_with_a_header_row(self, tmp_path):
        path = _Frame(output_dir=tmp_path).save_data(
            pd.DataFrame({"frame": [0, 1], "distance_nm": [0.4, 0.5]}),
            tmp_path / "pl.dat")

        assert np.allclose(
            _numeric_column(path, prefer="distance"), [0.4, 0.5])

    def test_a_2_5_4_file_with_no_note(self, tmp_path):
        """Every run already on disk. These must keep reading."""
        path = tmp_path / "legacy.dat"
        np.savetxt(path, np.array([[0.0, 0.1], [1.0, 0.2]]), fmt="%.8e")

        assert np.allclose(_numeric_column(path, prefer="nothing"), [0.1, 0.2])

    def test_a_plumed_fields_header(self, tmp_path):
        path = tmp_path / "colvar.dat"
        path.write_text("#! FIELDS time d1\n0.0 0.35\n1.0 0.41\n")

        assert np.allclose(_numeric_column(path, prefer="d1"), [0.35, 0.41])

    def test_a_file_of_nothing_but_notes_says_so(self, tmp_path):
        path = tmp_path / "empty.dat"
        path.write_text("# a note and no rows\n")

        with pytest.raises(ValueError, match="no data rows"):
            _numeric_column(path, prefer="anything")
