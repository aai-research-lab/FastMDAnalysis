"""A clustering of a long run writes no frame-by-frame matrix by default, and
one too large for the machine is refused with the stride that fits.

A chignolin study extended to 128,000 frames was killed by the kernel while
it clustered them: the float64 matrix (131 GB) and its float32 copy (66 GB),
and every study before it saved that matrix twice beside its results.
"""

from __future__ import annotations

import json
import logging
import math
import tracemalloc

import mdtraj as md
import numpy as np
import pytest

import fastmdxplora.analysis  # noqa: F401 -- registers the analyses
from fastmdxplora.analysis import cluster as clustering
from fastmdxplora.analysis.cluster import (
    MEMORY_SHARE,
    MOST_MAP_CELLS,
    Cluster,
    _pairwise_rmsd,
    memory_needed,
)


def _chain(n_frames: int, n_atoms: int = 4, seed: int = 0) -> md.Trajectory:
    """Many frames of a few alpha carbons, wandering about a shape, with a
    clock of 10 ps a frame."""
    rng = np.random.default_rng(seed)
    top = md.Topology()
    chain = top.add_chain()
    for i in range(n_atoms):
        residue = top.add_residue("ALA", chain, resSeq=i + 1)
        top.add_atom("CA", md.element.carbon, residue)
    shape = rng.normal(0.0, 0.4, (n_atoms, 3))
    drift = np.cumsum(rng.normal(0.0, 0.01, (n_frames, n_atoms, 3)), axis=0)
    traj = md.Trajectory((shape[None] + drift).astype(np.float32), top)
    traj.time = np.arange(n_frames) * 10.0
    return traj


#: Past the bound, so the map is of blocks of two frames.
LONG = MOST_MAP_CELLS + 100


@pytest.fixture(scope="module")
def clustered(tmp_path_factory):
    out = tmp_path_factory.mktemp("long")
    traj = _chain(LONG)
    result = Cluster(output_dir=out).run(traj)
    return traj, result


class TestALongRunWritesNoMatrixOverTheBound:

    def test_nothing_of_the_frames_square_is_written(self, clustered):
        traj, result = clustered
        assert result.status == "ok", result.message
        folder = result.output_dir
        assert not (folder / "hierarchical_distance_matrix.npy").exists()
        assert not (folder / "cluster_rmsd_matrix_full.npy").exists()
        assert (folder / "hierarchical_linkage.npy").is_file()
        saved = np.load(folder / "cluster_rmsd_matrix.npz")
        cells = math.ceil(traj.n_frames / 2)
        assert int(saved["block"]) == 2
        assert saved["rmsd_nm"].shape == (cells, cells)
        assert max(saved["rmsd_nm"].shape) <= MOST_MAP_CELLS
        # Every clustered frame's index and time, for the axis.
        assert saved["frames"].shape == (traj.n_frames,)
        assert np.allclose(saved["time_ns"], traj.time / 1000.0)

    def test_the_map_is_the_block_means_of_the_matrix(self, clustered):
        traj, result = clustered
        saved = np.load(result.output_dir / "cluster_rmsd_matrix.npz")
        atoms = traj.topology.select("name CA")
        full = _pairwise_rmsd(traj, atoms)
        starts = np.arange(0, traj.n_frames, 2)
        sizes = np.diff(np.append(starts, traj.n_frames))
        expected = (np.add.reduceat(np.add.reduceat(full, starts, axis=0), starts, axis=1)
                    / np.outer(sizes, sizes))
        assert np.allclose(saved["rmsd_nm"], expected, atol=1e-6)

    def test_the_block_is_recorded_and_said_with_the_figure(self, clustered):
        _, result = clustered
        found = json.loads((result.output_dir / "options.json").read_text())["findings"]
        assert found["matrix"]["block"] == 2
        assert "blocks of 2" in found["matrix"]["said"]
        assert found["memory"]["bytes_needed"] == memory_needed(
            LONG, ["kmeans", "hierarchical"], "rmsd", "average")

    def test_the_report_s_caption_says_the_block(self, clustered):
        from fastmdxplora.report.document import _caption_of_figure

        _, result = clustered
        caption = _caption_of_figure(result.output_dir, "cluster", "cluster_rmsd_matrix")
        assert "blocks of 2 frames" in caption

    def test_the_gui_reads_the_same_file(self, clustered, tmp_path):
        from fastmdxplora.gui.figure_data import MOST_CELLS, _rmsd_matrix

        _, result = clustered
        root = tmp_path / "study"
        folder = root / "analysis" / "cluster"
        folder.mkdir(parents=True)
        (folder / "cluster_rmsd_matrix.npz").write_bytes(
            (result.output_dir / "cluster_rmsd_matrix.npz").read_bytes())
        payload = _rmsd_matrix(root)
        saved = np.load(folder / "cluster_rmsd_matrix.npz")
        cells = saved["rmsd_nm"].shape[0]
        further = math.ceil(cells / MOST_CELLS)
        # Blocks of blocks: each cell is 2 x `further` frames.
        block = 2 * further
        assert payload["block"] == block
        assert len(payload["values"]) == math.ceil(cells / further)
        # Means over those frames, the last block's fewer frames weighed as
        # few: the same numbers as averaging the whole matrix so.
        traj, _ = clustered
        full = _pairwise_rmsd(traj, traj.topology.select("name CA"))
        starts = np.arange(0, traj.n_frames, block)
        sizes = np.diff(np.append(starts, traj.n_frames))
        expected = (np.add.reduceat(np.add.reduceat(full, starts, axis=0), starts, axis=1)
                    / np.outer(sizes, sizes))
        assert np.allclose(np.array(payload["values"], dtype=float), expected, rtol=1e-4)


class TestTheFullMatrixOnlyWhenAsked:
    def test_it_is_written_with_its_size_said_first(self, tmp_path, caplog):
        traj = _chain(120)
        with caplog.at_level(logging.INFO, logger="fastmdx.analysis.cluster"):
            result = Cluster(output_dir=tmp_path, save_full_matrix=True).run(traj)
        assert result.status == "ok", result.message
        path = result.output_dir / "cluster_rmsd_matrix_full.npy"
        assert path in result.artifacts
        written = np.load(path)
        assert written.dtype == np.float32 and written.shape == (120, 120)
        assert np.allclose(written, _pairwise_rmsd(traj, traj.topology.select("name CA")),
                           atol=1e-6)
        said = [r.getMessage() for r in caplog.records if "whole" in r.getMessage()]
        assert said and "120 by 120" in said[0] and "0.1 MB" in said[0]

    def test_the_setting_is_offered(self):
        from fastmdxplora.analysis.describe import describe_analysis

        names = {option.name for option in describe_analysis("cluster")}
        assert {"save_full_matrix", "stride"} <= names


class TestTheMemoryCheck:
    N = 600

    def _refused(self, tmp_path, monkeypatch, available: int):
        monkeypatch.setattr(clustering, "memory_available", lambda: (available, "mocked"))
        return Cluster(output_dir=tmp_path).run(_chain(self.N))

    def test_past_the_memory_it_is_refused_with_the_stride_that_fits(self, tmp_path, monkeypatch):
        from fastmdxplora.refusals import CODES

        methods = ["kmeans", "hierarchical"]
        available = memory_needed(self.N, methods, "rmsd", "average") // 5
        result = self._refused(tmp_path, monkeypatch, available)

        assert result.status == "error"
        found = json.loads((result.output_dir / "options.json").read_text())["findings"]
        assert found["refused"]["code"] == "analysis.memory.matrix_too_large"
        assert "analysis.memory.matrix_too_large" in {code.id for code in CODES}
        stride = found["memory"]["stride_that_fits"]
        fits = math.ceil(self.N / stride)
        budget = MEMORY_SHARE * available
        assert memory_needed(fits, methods, "rmsd", "average") <= budget
        assert memory_needed(math.ceil(self.N / (stride - 1)), methods, "rmsd",
                             "average") > budget
        # Said with the frames, the bytes and the setting.
        assert f"{self.N:,} frames" in result.message
        assert f"stride: {stride}" in result.message
        assert "available" in result.message
        assert found["memory"]["bytes_available"] == available

        # And the stride it names does fit.
        again = Cluster(output_dir=tmp_path / "again", stride=stride).run(_chain(self.N))
        assert again.status == "ok", again.message
        labels = np.loadtxt(again.output_dir / "cluster_kmeans.dat", delimiter=",",
                            skiprows=1)
        assert labels.shape[0] == fits
        assert list(labels[:3, 0]) == [0, stride, 2 * stride]

    def test_the_refusal_carries_its_numbers(self, monkeypatch):
        from fastmdxplora.refusals import StudyError, refusal_of

        monkeypatch.setattr(clustering, "memory_available", lambda: (10**6, "mocked"))
        with pytest.raises(StudyError) as raised:
            Cluster().compute(_chain(self.N))
        refusal = refusal_of(raised.value)
        assert refusal.code == "analysis.memory.matrix_too_large"
        assert refusal.details["frames"] == self.N
        assert refusal.details["bytes_available"] == 10**6
        assert refusal.details["stride"] > 1

    def test_within_the_memory_it_clusters(self, tmp_path, monkeypatch):
        result = self._refused(tmp_path, monkeypatch, 10**12)
        assert result.status == "ok", result.message

    def test_the_other_analyses_of_the_study_still_run(self, tmp_path, monkeypatch):
        from fastmdxplora.analysis import AnalysisOrchestrator

        traj = _chain(self.N)
        traj[0].save_pdb(str(tmp_path / "top.pdb"))
        traj.save_dcd(str(tmp_path / "traj.dcd"))
        monkeypatch.setattr(clustering, "memory_available", lambda: (10**6, "mocked"))
        results = AnalysisOrchestrator(str(tmp_path / "traj.dcd"), str(tmp_path / "top.pdb"),
                                       output_dir=tmp_path / "analysis").run(
            include=["rg", "cluster"])
        assert results["cluster"].status == "error"
        assert "stride:" in results["cluster"].message
        assert results["rg"].status == "ok"

    def test_where_memory_cannot_be_read_it_says_so_and_clusters(self, tmp_path, monkeypatch):
        monkeypatch.setattr(clustering, "memory_available", lambda: (None, "not read"))
        result = Cluster(output_dir=tmp_path).run(_chain(80))
        assert result.status == "ok"
        found = json.loads((result.output_dir / "options.json").read_text())["findings"]
        assert found["memory"]["bytes_available"] is None

    def test_memory_is_read_from_the_machine(self):
        available, where = clustering.memory_available()
        assert available is None or available > 0
        assert where


class TestTheEstimateIsWhatTheMethodsMake:
    """The bytes the check counts against the peak the methods traced: never
    under it by more than terms in n alone, and not far over it."""

    @pytest.mark.parametrize("options", [
        {},
        {"methods": ["kmeans"]},
        {"methods": ["hierarchical"], "linkage": "ward"},
        {"methods": ["dbscan"], "eps": 5.0},
        {"features": "coordinates", "methods": ["hierarchical"]},
    ])
    def test_the_peak_is_within_the_estimate(self, options):
        n = 700
        traj = _chain(n)
        analysis = Cluster(**options)
        tracemalloc.start()
        try:
            before = tracemalloc.get_traced_memory()[0]
            analysis.compute(traj)
            peak = tracemalloc.get_traced_memory()[1] - before
        finally:
            tracemalloc.stop()
        estimate = memory_needed(n, analysis.methods, analysis.features, analysis.linkage)
        assert peak <= estimate + 2000 * n, (peak / n**2, estimate / n**2)
        assert estimate <= 1.35 * peak, (peak / n**2, estimate / n**2)

    def test_k_means_on_coordinates_makes_no_matrix(self):
        analysis = Cluster(features="coordinates", methods=["kmeans", "dbscan"],
                           linkage="ward")
        analysis.compute(_chain(300))
        assert analysis._distances is None
        assert memory_needed(10**6, ["kmeans"], "coordinates", "average") == 0
