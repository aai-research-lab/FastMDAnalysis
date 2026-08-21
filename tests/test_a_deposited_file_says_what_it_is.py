"""Two `.dat` conventions under one extension, and now each says which.

`pl_interactions.dat` is comma-separated with a header row; `rdf.dat` is
whitespace with none. Same extension, same directory tree, and a reader who
guesses wrong gets a `ValueError` at best and a column of `NaN` at worst.

`options.json` has recorded the answer since 2.5.4, which serves a reader
who knows to look and still has the run directory. A deposited file
outlives its directory -- it goes into a supplementary archive, an email, a
student's home folder -- so the whitespace form now carries its layout on a
`#` line. `np.loadtxt` skips those by default, which is what keeps this a
note to the reader rather than a change to the format.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fastmdxplora.analysis.base import Analysis


class _Array(Analysis):
    name = "arrayish"

    def compute(self, traj):  # pragma: no cover - not run
        raise NotImplementedError

    def plot(self, result, ax):  # pragma: no cover - not run
        raise NotImplementedError


class _Frame(Analysis):
    name = "framish"

    def compute(self, traj):  # pragma: no cover - not run
        raise NotImplementedError

    def plot(self, result, ax):  # pragma: no cover - not run
        raise NotImplementedError


@pytest.fixture
def values():
    return np.array([0.0, 1.5, 2.25, 3.125])


class TestTheFileNamesItsOwnLayout:
    def test_the_whitespace_form_says_how_to_read_it(self, values, tmp_path):
        path = _Array(output_dir=tmp_path).save_data(values, tmp_path / "a.dat")
        first = path.read_text().splitlines()[0]

        assert first.startswith("#")
        assert "arrayish" in first
        assert "np.loadtxt" in first

    def test_saying_so_does_not_change_what_it_holds(self, values, tmp_path):
        """The point of a `#` line: `np.loadtxt` skips it untold."""
        path = _Array(output_dir=tmp_path).save_data(values, tmp_path / "a.dat")

        assert np.allclose(np.loadtxt(path), values)

    def test_the_comma_form_is_left_alone(self, tmp_path):
        """Its header row already names its columns; a `#` line above one
        would break `pd.read_csv` for a reader who did nothing wrong."""
        frame = pd.DataFrame({"residue": ["LEU118"], "occupancy": [0.42]})
        path = _Frame(output_dir=tmp_path).save_data(frame, tmp_path / "b.dat")

        assert not path.read_text().startswith("#")
        assert pd.read_csv(path)["residue"][0] == "LEU118"

    def test_the_record_and_the_file_agree(self, values, tmp_path):
        """`options.json` and the `#` line must not drift apart."""
        analysis = _Array(output_dir=tmp_path)
        path = analysis.save_data(values, tmp_path / "a.dat")

        recorded = analysis._data_format
        assert "whitespace" in recorded["layout"]
        assert "np.loadtxt" in recorded["read_with"]
        assert "whitespace" in path.read_text().splitlines()[0]
