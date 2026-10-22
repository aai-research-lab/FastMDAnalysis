"""Every frame's solute is checked against its own periodic images, as
``gmx mindist -pi`` checks it, and every analysis of those frames reads what
the check found.

A chignolin chain came to 0.261 nm of its own copy in a study whose check of
the chain's two ends passed it.
"""

from __future__ import annotations

import itertools
import json

import mdtraj as md
import numpy as np
import pytest
from scipy.spatial.distance import cdist

import fastmdxplora.analysis  # noqa: F401 -- registers the analyses
from fastmdxplora.analysis.image_check import (
    CLOSE_NM,
    CUTOFF_NM,
    SoluteImage,
    solute_image_record,
)


def _chain(path: np.ndarray, *, bonded: bool = True) -> md.Topology:
    """One residue per point of ``path``, a CB beside each CA and a hydrogen
    on it, the heavy atoms bonded along the chain."""
    top = md.Topology()
    chain = top.add_chain()
    previous = None
    for i in range(len(path)):
        residue = top.add_residue("ALA", chain, resSeq=i + 1)
        ca = top.add_atom("CA", md.element.carbon, residue)
        cb = top.add_atom("CB", md.element.carbon, residue)
        top.add_atom("HA", md.element.hydrogen, residue)
        if bonded:
            top.add_bond(ca, cb)
            if previous is not None:
                top.add_bond(previous, ca)
        previous = ca
    return top


def _positions(path: np.ndarray) -> np.ndarray:
    """CA on the path, CB 0.15 nm off it in z, HA 0.1 nm off it in -z."""
    xyz = []
    for point in path:
        xyz += [point, point + [0.0, 0.0, 0.15], point - [0.0, 0.0, 0.1]]
    return np.asarray(xyz, dtype=np.float64)


def _dodecahedron(a: float) -> np.ndarray:
    """OpenMM's reduced rhombic dodecahedron, as setup builds it."""
    return np.array([[a, 0, 0], [0, a, 0], [a / 2, a / 2, a * np.sqrt(2) / 2]])


def _octahedron(a: float) -> np.ndarray:
    """OpenMM's reduced truncated octahedron."""
    return np.array([[a, 0, 0], [a / 3, 2 * np.sqrt(2) * a / 3, 0],
                     [-a / 3, np.sqrt(2) * a / 3, np.sqrt(6) * a / 3]])


def _brute(points: np.ndarray, box: np.ndarray, cells: int = 3) -> float:
    """Every pair against every image within three cells, no shortcuts."""
    best = np.inf
    for shift in itertools.product(range(-cells, cells + 1), repeat=3):
        if shift == (0, 0, 0):
            continue
        image = points + np.asarray(shift, dtype=np.float64) @ box
        best = min(best, float(cdist(points, image).min()))
    return best


def _wandering(n_frames: int, box: np.ndarray, seed: int,
               bonded: bool = True) -> tuple[md.Trajectory, np.ndarray]:
    """A 14-residue chain, a random walk of 0.38 nm steps from compact to
    stretched along x, each frame turned at random, moved, and then wrapped
    atom by atom into the cell, as an engine may write it. With it, each
    frame's whole chain as the trajectory stores it: the wrapped coordinates
    put back by the whole cells each atom was moved, in the trajectory's own
    cell, in float64."""
    rng = np.random.default_rng(seed)
    frames = []
    for f in range(n_frames):
        steps = rng.normal(size=(14, 3))
        steps[:, 0] += 4.5 * f / max(n_frames - 1, 1)
        steps /= np.linalg.norm(steps, axis=1)[:, None]
        path = np.cumsum(0.38 * steps, axis=0)
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        frames.append(_positions(path) @ q.T + rng.uniform(-2, 2, 3))
    xyz = np.asarray(frames)
    boxes = np.repeat(box[None], n_frames, axis=0)
    cells = np.floor(np.einsum("fai,fij->faj", xyz, np.linalg.inv(boxes)))
    wrapped = xyz - np.einsum("faj,fji->fai", cells, boxes)
    traj = md.Trajectory(wrapped.astype(np.float32),
                         _chain(np.zeros((14, 3)), bonded=bonded))
    traj.unitcell_vectors = boxes.astype(np.float32)
    traj.time = np.arange(n_frames) * 25.0
    stored = np.asarray(traj.unitcell_vectors, dtype=np.float64)
    whole = traj.xyz.astype(np.float64) + np.einsum("faj,fji->fai", cells, stored)
    # Split across the cell's faces as stored, so making it whole is tested.
    assert np.linalg.norm(np.diff(traj.xyz[:, ::3], axis=1), axis=2).max() > 1.0
    return traj, whole


class TestTheDistanceIsTheBruteForceOne:
    @pytest.mark.parametrize("box", [
        np.diag([3.2, 3.2, 3.2]),
        _dodecahedron(3.4),
        _octahedron(3.6),
    ], ids=["cube", "dodecahedron", "octahedron"])
    def test_to_a_millionth_of_a_nanometre(self, box):
        traj, whole = _wandering(6, box, seed=1)
        heavy = traj.topology.select("not element H")
        distances, record = solute_image_record(traj, np.arange(traj.n_atoms))
        cells = np.asarray(traj.unitcell_vectors, dtype=np.float64)
        for f in range(traj.n_frames):
            expected = _brute(whole[f][heavy], cells[f])
            assert distances[f] == pytest.approx(expected, abs=1e-6), (f, expected)
        assert record["frames"] == 6
        assert record["made_whole"] == "along its bonds"
        assert record["nearest_nm"] == pytest.approx(distances.min())

    def test_hydrogens_are_not_counted(self):
        traj, whole = _wandering(3, np.diag([3.0, 3.0, 3.0]), seed=2)
        distances, record = solute_image_record(traj, np.arange(traj.n_atoms))
        heavy = traj.topology.select("not element H")
        assert record["heavy_atoms"] == heavy.size == 28
        cells = np.asarray(traj.unitcell_vectors, dtype=np.float64)
        for f in range(3):
            assert distances[f] == pytest.approx(_brute(whole[f][heavy], cells[f]), abs=1e-6)

    def test_without_bonds_the_chain_is_walked_in_order(self):
        traj, whole = _wandering(4, _dodecahedron(3.4), seed=3, bonded=False)
        distances, record = solute_image_record(traj, np.arange(traj.n_atoms))
        heavy = traj.topology.select("not element H")
        assert "atom order" in record["made_whole"]
        cells = np.asarray(traj.unitcell_vectors, dtype=np.float64)
        for f in range(4):
            assert distances[f] == pytest.approx(_brute(whole[f][heavy], cells[f]), abs=1e-6)


def _stretched(n_frames: int = 40, box: float = 4.0, reaching: int = 10,
               touching: int = 3) -> md.Trajectory:
    """A 10-residue chain, CA 0.38 nm apart along x in a cube: folded back on
    itself in most frames, stretched straight (3.42 nm, the CB at its ends
    0.58 nm from the image of the other end) in ``reaching`` of them, and in
    ``touching`` of those pushed to 0.4 nm."""
    straight = np.c_[0.38 * np.arange(10), np.zeros(10), np.zeros(10)]
    hairpin = np.r_[straight[:5], straight[:5][::-1] + [0.0, 0.5, 0.0]]
    frames = []
    for f in range(n_frames):
        if f < reaching:
            path = straight.copy()
            if f < touching:
                path[:, 0] *= (box - 0.4) / (9 * 0.38)
        else:
            path = hairpin
        frames.append(_positions(path) + 0.3)
    traj = md.Trajectory(np.asarray(frames, dtype=np.float32), _chain(straight))
    traj.unitcell_lengths = np.full((n_frames, 3), box, dtype=np.float32)
    traj.unitcell_angles = np.full((n_frames, 3), 90.0, dtype=np.float32)
    traj.time = np.arange(n_frames) * 10.0
    return traj


class TestAChainAcrossASmallBox:
    def test_the_frames_under_the_cutoff_are_counted(self):
        traj = _stretched()
        distances, record = solute_image_record(traj, np.arange(traj.n_atoms))
        assert record["reaches_its_image"] is True
        assert record["frames_under_cutoff"] == 10
        assert record["frames_under_close"] == 3
        assert record["fraction_under_cutoff"] == pytest.approx(0.25)
        assert record["nearest_nm"] == pytest.approx(0.4, abs=1e-5)
        assert record["nearest_frame"] in (0, 1, 2)
        assert record["nearest_time_ns"] == pytest.approx(record["nearest_frame"] / 100)
        assert {record["nearest_atoms"][0].split()[0],
                record["nearest_atoms"][1].split()[0]} == {"ALA1", "ALA10"}
        assert "its own copy" in record["said"] and "10 of 40" in record["said"]
        assert np.all(distances[10:] > CUTOFF_NM)

    def test_a_chain_clear_of_its_image_says_so(self):
        traj = _stretched(reaching=0, touching=0)
        _, record = solute_image_record(traj, np.arange(traj.n_atoms))
        assert record["reaches_its_image"] is False
        assert "no closer than" in record["said"]

    def test_it_is_a_series_with_its_figure(self, tmp_path):
        result = SoluteImage(output_dir=tmp_path).run(_stretched())
        assert result.status == "ok", result.message
        assert np.loadtxt(result.data_path).shape == (40,)
        assert result.figure_path.is_file()
        found = json.loads((tmp_path / "image_check" / "options.json").read_text())["findings"]
        assert found["solute_image"]["frames_under_cutoff"] == 10
        # The check's own mean carries what it found, where the Overview,
        # the Agent and the report read a mean's reason.
        assert "its own copy" in found["mean"]["not_a_measurement"]

    def test_the_close_threshold_is_half_a_nanometre(self):
        assert CLOSE_NM == 0.5 and CUTOFF_NM == 1.0


def _hairpin_reaching_its_image() -> md.Trajectory:
    """A hairpin whose two ends sit side by side (their images 2.6 nm away)
    and whose turn, 2.4 nm up, comes 0.6 nm from the image of its base: the
    ends pass, the chain does not."""
    up = np.c_[np.zeros(8), 0.343 * np.arange(8), np.zeros(8)]
    down = (up + [0.4, 0.0, 0.0])[::-1]
    path = np.r_[up, down]
    top = md.Topology()
    chain = top.add_chain()
    previous = None
    for i in range(len(path)):
        residue = top.add_residue("ALA", chain, resSeq=i + 1)
        atom = top.add_atom("CA", md.element.carbon, residue)
        if previous is not None:
            top.add_bond(previous, atom)
        previous = atom
    rng = np.random.default_rng(0)
    xyz = path[None] + 0.2 + rng.normal(0.0, 0.005, (200, len(path), 3))
    traj = md.Trajectory(xyz.astype(np.float32), top)
    traj.unitcell_lengths = np.full((200, 3), 3.0, dtype=np.float32)
    traj.unitcell_angles = np.full((200, 3), 90.0, dtype=np.float32)
    traj.time = np.arange(200) * 10.0
    return traj


class TestTheEndToEndDistanceReadsIt:
    def test_a_chain_whose_ends_pass_is_marked_by_the_whole_chain(self, tmp_path):
        from fastmdxplora.analysis.end_to_end import EndToEndDistance

        analysis = EndToEndDistance(output_dir=tmp_path)
        assert analysis.run(_hairpin_reaching_its_image()).status == "ok"
        found = analysis.findings
        assert found["periodic"]["nearest_self_image_nm"] > 2.0     # the ends
        assert found["periodic"]["nearest_self_image_whole_chain_nm"] < 1.0
        assert "its own copy" in found["not_a_measurement"]
        assert not np.isfinite(found["mean"]["standard_error"])

    def test_handed_the_finding_it_reads_it_rather_than_its_own(self):
        from fastmdxplora.analysis.end_to_end import EndToEndDistance

        analysis = EndToEndDistance()
        analysis.solute_image = {"reaches_its_image": True, "said": "Handed: its own copy.",
                                 "nearest_nm": 0.3}
        traj = _stretched(reaching=0, touching=0)
        analysis.compute(traj)
        assert "Handed: its own copy." in analysis.findings["not_a_measurement"]


@pytest.fixture(scope="module")
def analysed(tmp_path_factory):
    from fastmdxplora.analysis import AnalysisOrchestrator

    folder = tmp_path_factory.mktemp("study")
    traj = _hairpin_reaching_its_image()
    traj[0].save_pdb(str(folder / "top.pdb"))
    traj.save_dcd(str(folder / "traj.dcd"))
    orchestrator = AnalysisOrchestrator(str(folder / "traj.dcd"), str(folder / "top.pdb"),
                                        output_dir=folder / "analysis")
    results = orchestrator.run(include=["rg", "end_to_end", "cluster", "image_check"])
    return folder / "analysis", results


class TestEveryAnalysisOfTheFramesReadsIt:
    def _found(self, folder, name):
        return json.loads((folder / name / "options.json").read_text())["findings"]

    def test_the_check_runs_first(self, analysed):
        _, results = analysed
        assert list(results)[0] == "image_check"
        assert all(r.status == "ok" for r in results.values()), {
            n: r.message for n, r in results.items()}

    def test_each_mean_carries_the_one_reason_once(self, analysed):
        folder, _ = analysed
        said = self._found(folder, "image_check")["solute_image"]["said"]
        for name in ("rg", "end_to_end", "image_check"):
            mean = self._found(folder, name)["mean"]
            assert mean["not_a_measurement"].count(said) == 1, name
            assert not np.isfinite(mean["standard_error"]), name
            assert mean["error_withheld_because"], name

    def test_an_analysis_with_no_mean_says_it_too(self, analysed):
        folder, _ = analysed
        found = self._found(folder, "cluster")
        assert found["solute_image"]["reaches_its_image"] is True
        assert "its own copy" in found["solute_image"]["said"]

    def test_the_report_says_it(self, analysed):
        from fastmdxplora.report.document import _findings_notes

        folder, _ = analysed
        for name in ("rg", "cluster", "image_check"):
            notes = _findings_notes(self._found(folder, name))
            assert sum("its own copy" in note for note in notes) == 1, (name, notes)

    def test_a_study_without_a_box_does_not_plan_it(self, tmp_path):
        from fastmdxplora.analysis import AnalysisOrchestrator

        traj = _hairpin_reaching_its_image()
        traj.unitcell_vectors = None
        traj[0].save_pdb(str(tmp_path / "top.pdb"))
        traj.save_dcd(str(tmp_path / "traj.dcd"))
        orchestrator = AnalysisOrchestrator(str(tmp_path / "traj.dcd"), str(tmp_path / "top.pdb"),
                                            output_dir=tmp_path / "analysis")
        assert "image_check" not in orchestrator._build_plan(None, None)

    def test_a_study_with_a_box_plans_it_first(self, analysed):
        from fastmdxplora.analysis import AnalysisOrchestrator

        folder, _ = analysed
        orchestrator = AnalysisOrchestrator(str(folder.parent / "traj.dcd"),
                                            str(folder.parent / "top.pdb"),
                                            output_dir=folder.parent / "again")
        assert orchestrator._build_plan(None, None)[0] == "image_check"
