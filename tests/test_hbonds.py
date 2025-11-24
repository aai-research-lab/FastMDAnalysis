# tests/test_hbonds.py
from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest

import fastmdanalysis.analysis.hbonds as hbonds_mod
from fastmdanalysis.analysis.hbonds import HBondsAnalysis, AnalysisError


# ---- Minimal fakes for topology/trajectory ----------------------------------

class _FakeTopology:
    """
    Minimal topology stub to drive the hbonds logic.

    Parameters control behavior of:
      - n_bonds (can raise, or be 0/nonzero)
      - select(query) (mapping + optional raise for 'protein')
      - create_standard_bonds() (may raise, may set bonds to nonzero)
    """
    def __init__(
        self,
        n_atoms: int,
        n_bonds: int = 1,
        select_map: dict[str, list[int] | None] | None = None,
        raise_on_n_bonds: bool = False,
        raise_on_create_bonds: bool = False,
        raise_on_protein_select: bool = False,
        bonds_after_create: int | None = 1,
    ):
        self._n_atoms = n_atoms
        self._n_bonds = n_bonds
        self._select_map = select_map or {}
        self._raise_on_n_bonds = raise_on_n_bonds
        self._raise_on_create_bonds = raise_on_create_bonds
        self._raise_on_protein_select = raise_on_protein_select
        self._bonds_after_create = bonds_after_create

    @property
    def n_bonds(self) -> int:
        if self._raise_on_n_bonds:
            raise RuntimeError("n_bonds access failed")
        return self._n_bonds

    def select(self, query: str):
        if self._raise_on_protein_select and query == "protein":
            raise RuntimeError("protein selection failed")
        # Return explicit mapping if present, else "all atoms"
        return self._select_map.get(query, list(range(self._n_atoms)))

    def create_standard_bonds(self):
        if self._raise_on_create_bonds:
            raise RuntimeError("bond creation failed")
        if self._bonds_after_create is not None:
            self._n_bonds = self._bonds_after_create


class _FakeFrame:
    """Single-frame stub; baker_hubbard() ignores its internals for these tests."""
    pass


class _FakeTrajectory:
    """
    Minimal trajectory stub:
      - Provides .topology, .n_frames, .n_atoms
      - __getitem__ returns a one-frame-like object
      - .atom_slice returns a pre-plumbed sliced trajectory (or self)
    """
    def __init__(self, n_frames: int, topology: _FakeTopology, atom_slice_return: "_FakeTrajectory" | None = None):
        self.n_frames = n_frames
        self.topology = topology
        self.n_atoms = topology._n_atoms
        self._atom_slice_return = atom_slice_return

    def __getitem__(self, idx):
        return _FakeFrame()

    def atom_slice(self, idx):
        return self._atom_slice_return or self


# ---- Helpers ----------------------------------------------------------------

def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ok")
    return path


@pytest.fixture
def baker_hubbard_counts(monkeypatch):
    """
    Stub md.baker_hubbard to yield 0, 1, 2 H-bonds on successive calls.
    This lets us verify per-frame counting deterministically.
    """
    calls = {"n": 0}
    empt = np.empty((0, 3), dtype=int)
    one = np.array([[0, 1, 2]], dtype=int)
    two = np.array([[0, 1, 2], [0, 1, 2]], dtype=int)

    def fake_baker(_frame, periodic=False):
        i = calls["n"]
        calls["n"] += 1
        return [empt, one, two][min(i, 2)]

    monkeypatch.setattr(hbonds_mod.md, "baker_hubbard", fake_baker)
    return calls


# ---- Tests ------------------------------------------------------------------

def test_hbonds_run_basic_no_fallback(tmp_path, monkeypatch, baker_hubbard_counts):
    # Work trajectory has bonds -> no fallback, plot is auto-called
    topo = _FakeTopology(n_atoms=5, n_bonds=1)
    traj = _FakeTrajectory(n_frames=3, topology=topo)

    h = HBondsAnalysis(traj)
    h.outdir = tmp_path

    # Avoid disk I/O complexity for data, and spy on auto-plot
    saved = {}
    plot_calls = {"n": 0}
    monkeypatch.setattr(HBondsAnalysis, "_save_data", lambda self, arr, stem: saved.setdefault("data", arr))
    def fake_plot(self, *args, **kwargs):
        plot_calls["n"] += 1
        return _touch(tmp_path / "hbonds.png")
    monkeypatch.setattr(HBondsAnalysis, "plot", fake_plot)

    res = h.run()
    assert res["fallback"] is False
    assert res["selection_used"] == "all atoms"
    assert h.data.shape == (3, 1)
    assert h.data.flatten().tolist() == [0, 1, 2]
    assert np.array_equal(saved["data"], h.data)
    assert plot_calls["n"] == 1
    assert (tmp_path / "hbonds.png").exists()


def test_hbonds_atoms_selection_empty_raises():
    # Atoms string maps to empty selection -> immediate AnalysisError
    topo = _FakeTopology(n_atoms=10, select_map={"name CA": []}, n_bonds=1)
    traj = _FakeTrajectory(n_frames=1, topology=topo)
    h = HBondsAnalysis(traj, atoms="name CA")
    with pytest.raises(AnalysisError):
        h._prepare_work_trajectory()


def test_hbonds_fallback_to_protein_and_note_file(tmp_path, monkeypatch, baker_hubbard_counts):
    # Initial work has zero bonds (even after create), but protein slice has bonds -> fallback used
    work_topo = _FakeTopology(n_atoms=10, n_bonds=0, bonds_after_create=0)  # remains unbonded
    protein_topo = _FakeTopology(n_atoms=6, n_bonds=1)
    protein_traj = _FakeTrajectory(n_frames=3, topology=protein_topo)
    # Selecting "protein" returns non-empty; slicing returns protein_traj
    work_topo._select_map = {"protein": list(range(6))}
    traj = _FakeTrajectory(n_frames=3, topology=work_topo, atom_slice_return=protein_traj)

    h = HBondsAnalysis(traj)
    h.outdir = tmp_path

    # Spy auto-plot to avoid real mpl; still exercise run() path
    plot_calls = {"n": 0}
    monkeypatch.setattr(HBondsAnalysis, "_save_data", lambda self, arr, stem: None)
    def fake_plot(self, *args, **kwargs):
        plot_calls["n"] += 1
        return _touch(tmp_path / "hbonds.png")
    monkeypatch.setattr(HBondsAnalysis, "plot", fake_plot)

    res = h.run()
    assert res["fallback"] is True
    assert res["selection_used"] == "protein (fallback)"
    assert plot_calls["n"] == 1
    assert (tmp_path / "hbonds_NOTE.txt").exists()


def test_hbonds_fallback_protein_empty_uses_full_traj_label_all_atoms(tmp_path, monkeypatch, baker_hubbard_counts):
    """
    Make the initial 'work' trajectory an atoms subset (Cα-only) with no bonds,
    then make 'protein' selection empty → fallback to FULL trajectory where
    bonds can be created. This ensures used_fallback=True and label "all atoms (fallback)".
    """
    # Full trajectory topology: will gain bonds when create_standard_bonds() is called
    full_topo = _FakeTopology(n_atoms=8, n_bonds=0, bonds_after_create=1, select_map={"protein": []})
    full_traj = _FakeTrajectory(n_frames=2, topology=full_topo)

    # Sliced trajectory for atoms selection (no bonds even after create)
    ca_topo = _FakeTopology(n_atoms=3, n_bonds=0, bonds_after_create=0)
    ca_traj = _FakeTrajectory(n_frames=2, topology=ca_topo)

    # Atom slicing returns the CA-only traj
    full_traj._atom_slice_return = ca_traj
    # Atoms selection must return non-empty indices to build the slice
    full_topo._select_map["name CA"] = [0, 2, 4]

    h = HBondsAnalysis(full_traj, atoms="name CA")
    h.outdir = tmp_path

    plot_calls = {"n": 0}
    monkeypatch.setattr(HBondsAnalysis, "_save_data", lambda self, arr, stem: None)
    def fake_plot(self, *args, **kwargs):
        plot_calls["n"] += 1
        return _touch(tmp_path / "hbonds.png")
    monkeypatch.setattr(HBondsAnalysis, "plot", fake_plot)

    res = h.run()
    assert res["fallback"] is True
    assert res["selection_used"] == "all atoms (fallback)"
    assert plot_calls["n"] == 1
    assert (tmp_path / "hbonds_NOTE.txt").exists()


def test_hbonds_fallback_still_no_bonds_raises(monkeypatch):
    # Work has no bonds; protein slice also ends up with no bonds -> final AnalysisError
    work_topo = _FakeTopology(n_atoms=10, n_bonds=0, bonds_after_create=0, select_map={"protein": list(range(4))})
    protein_topo = _FakeTopology(n_atoms=4, n_bonds=0, bonds_after_create=0)
    protein_traj = _FakeTrajectory(n_frames=1, topology=protein_topo)
    traj = _FakeTrajectory(n_frames=1, topology=work_topo, atom_slice_return=protein_traj)

    h = HBondsAnalysis(traj)
    with pytest.raises(AnalysisError):
        h.run()


def test__has_bonds_handles_exception():
    # n_bonds property raising should be treated as "no bonds"
    topo = _FakeTopology(n_atoms=3, n_bonds=0, raise_on_n_bonds=True)
    traj = _FakeTrajectory(n_frames=1, topology=topo)
    assert HBondsAnalysis._has_bonds(traj) is False


def test_plot_without_data_raises():
    # plot() with no data set should raise AnalysisError
    topo = _FakeTopology(n_atoms=3, n_bonds=1)
    traj = _FakeTrajectory(n_frames=1, topology=topo)
    h = HBondsAnalysis(traj)
    with pytest.raises(AnalysisError):
        h.plot()


def test_plot_custom_kwargs_and_color_branch(tmp_path, monkeypatch):
    """
    Exercise plot() without invoking real Matplotlib internals that can fail on
    some versions (e.g., _NoValueType sentinel). We stub pyplot functions enough
    to verify branches, input coercion, and _save_plot invocation.
    """
    # Stub pyplot surface
    class _DummyLabel:
        def set_fontsize(self, *a, **k):
            pass

    class _DummyAxis:
        def __init__(self):
            self.label = _DummyLabel()

        def set_major_locator(self, *a, **k):
            pass

    class _DummyTitle:
        def __init__(self):
            self._size = 20

        def set_fontsize(self, size, *a, **k):
            self._size = size

        def get_fontsize(self):
            return self._size

    class _DummyAxes:
        def __init__(self):
            self.xaxis = _DummyAxis()
            self.yaxis = _DummyAxis()
            self.title = _DummyTitle()

        def plot(self, *a, **k):
            pass

        def set_title(self, *a, **k):
            pass

        def set_xlabel(self, *a, **k):
            pass

        def set_ylabel(self, *a, **k):
            pass

        def grid(self, *a, **k):
            pass

        def tick_params(self, *a, **k):
            pass

        def get_xticklabels(self):
            return []

    class _DummyFig:
        def subplots(self, *a, **k):
            return _DummyAxes()

    monkeypatch.setattr(hbonds_mod.plt, "figure", lambda *a, **k: _DummyFig())
    monkeypatch.setattr(hbonds_mod.plt, "plot", lambda *a, **k: None)
    monkeypatch.setattr(hbonds_mod.plt, "title", lambda *a, **k: None)
    monkeypatch.setattr(hbonds_mod.plt, "xlabel", lambda *a, **k: None)
    monkeypatch.setattr(hbonds_mod.plt, "ylabel", lambda *a, **k: None)
    monkeypatch.setattr(hbonds_mod.plt, "grid", lambda *a, **k: None)
    monkeypatch.setattr(hbonds_mod.plt, "gca", lambda *a, **k: _DummyAxes())
    monkeypatch.setattr(hbonds_mod.plt, "close", lambda *a, **k: None)

    topo = _FakeTopology(n_atoms=3, n_bonds=1)
    traj = _FakeTrajectory(n_frames=1, topology=topo)
    h = HBondsAnalysis(traj)
    h.outdir = tmp_path
    h.data = np.array([[0.0], [2.0], [1.0]], dtype=float)

    monkeypatch.setattr(HBondsAnalysis, "_save_plot", lambda self, fig, stem: _touch(tmp_path / f"{stem}.png"))

    out = h.plot(
        color="black",
        linestyle="--",
        marker="o",
        title="HB vs Frame",
        xlabel="f",
        ylabel="count",
    )
    assert Path(out).exists()


def test_run_wraps_unknown_exception(monkeypatch):
    # Force an unexpected exception in run() to hit the generic wrapper
    topo = _FakeTopology(n_atoms=3, n_bonds=1)
    traj = _FakeTrajectory(n_frames=1, topology=topo)
    h = HBondsAnalysis(traj)

    def boom(*args, **kwargs):
        raise ValueError("unexpected")

    monkeypatch.setattr(HBondsAnalysis, "_prepare_work_trajectory", boom)
    with pytest.raises(AnalysisError) as ei:
        h.run()
    assert "Hydrogen bonds analysis failed: unexpected" in str(ei.value)


def test_hbonds_work_create_bonds_raises_but_has_bonds(tmp_path, monkeypatch, baker_hubbard_counts):
    """
    Work traj already has bonds, but create_standard_bonds raises.
    We should ignore that and proceed without fallback.
    """
    topo = _FakeTopology(n_atoms=5, n_bonds=1, raise_on_create_bonds=True)
    traj = _FakeTrajectory(n_frames=2, topology=topo)
    h = HBondsAnalysis(traj)
    h.outdir = tmp_path

    monkeypatch.setattr(HBondsAnalysis, "_save_data", lambda self, arr, stem: None)
    plot_calls = {"n": 0}
    def fake_plot(self, *args, **kwargs):
        plot_calls["n"] += 1
        return _touch(tmp_path / "hbonds.png")
    monkeypatch.setattr(HBondsAnalysis, "plot", fake_plot)

    res = h.run()
    assert res["fallback"] is False
    assert plot_calls["n"] == 1

def test_hbonds_fallback_protein_select_raises_uses_full_traj(tmp_path, monkeypatch, baker_hubbard_counts):
    """
    Selecting 'protein' raises → fallback should use FULL trajectory.
    We simulate that .atom_slice(full_range) returns a distinct fallback trajectory
    whose bonds can be created, while the initial work traj remains unbonded.
    The label should be 'protein (fallback)' per the implementation.
    """
    class _SmartSliceTrajectory(_FakeTrajectory):
        def __init__(self, n_frames, topology):
            super().__init__(n_frames, topology)
            self._fallback_traj = None

        def set_fallback(self, fb_traj):
            self._fallback_traj = fb_traj

        def atom_slice(self, idx):
            # If idx is the full-range [0..n_atoms-1], return the dedicated fallback traj.
            try:
                arr = np.asarray(idx)
                if arr.ndim == 1 and len(arr) == self.n_atoms and np.all(arr == np.arange(self.n_atoms)):
                    return self._fallback_traj or self
            except Exception:
                pass
            return super().atom_slice(idx)

    # Work traj: no bonds even after create; selecting 'protein' raises
    work_topo = _FakeTopology(n_atoms=6, n_bonds=0, bonds_after_create=0, raise_on_protein_select=True)
    work_traj = _SmartSliceTrajectory(n_frames=2, topology=work_topo)

    # Fallback traj: bonds will be created successfully
    fb_topo = _FakeTopology(n_atoms=6, n_bonds=0, bonds_after_create=1)
    fb_traj = _FakeTrajectory(n_frames=2, topology=fb_topo)
    work_traj.set_fallback(fb_traj)

    h = HBondsAnalysis(work_traj)
    h.outdir = tmp_path

    # Keep I/O light and assert the auto-plot call
    monkeypatch.setattr(HBondsAnalysis, "_save_data", lambda self, arr, stem: None)
    plot_calls = {"n": 0}
    def fake_plot(self, *args, **kwargs):
        plot_calls["n"] += 1
        return _touch(tmp_path / "hbonds.png")
    monkeypatch.setattr(HBondsAnalysis, "plot", fake_plot)

    res = h.run()
    assert res["fallback"] is True
    assert res["selection_used"] == "protein (fallback)"
    assert plot_calls["n"] == 1





def test_hbonds_atoms_selection_label_and_no_note(tmp_path, monkeypatch, baker_hubbard_counts):
    """
    Atoms selection is valid → label equals selection string and no fallback.
    Confirm no NOTE file is written.
    """
    base = _FakeTopology(n_atoms=10, n_bonds=1, select_map={"name CA": [0, 2, 4]})
    sliced_topo = _FakeTopology(n_atoms=3, n_bonds=1)
    sliced_traj = _FakeTrajectory(n_frames=3, topology=sliced_topo)
    traj = _FakeTrajectory(n_frames=3, topology=base, atom_slice_return=sliced_traj)

    h = HBondsAnalysis(traj, atoms="name CA")
    h.outdir = tmp_path

    monkeypatch.setattr(HBondsAnalysis, "_save_data", lambda self, arr, stem: None)
    plot_calls = {"n": 0}
    def fake_plot(self, *args, **kwargs):
        plot_calls["n"] += 1
        return _touch(tmp_path / "hbonds.png")
    monkeypatch.setattr(HBondsAnalysis, "plot", fake_plot)

    res = h.run()
    assert res["fallback"] is False
    assert res["selection_used"] == "name CA"
    assert plot_calls["n"] == 1
    assert not (tmp_path / "hbonds_NOTE.txt").exists()

