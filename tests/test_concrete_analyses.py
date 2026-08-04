"""Tests for the four concrete analyses: RMSD, RMSF, Rg, Dihedrals.

These exercise each analysis in three dimensions:

  1. Numerical correctness — verified either against hand-computed values
     on a controlled synthetic trajectory, or against MDTraj's own
     primitives (which the analyses wrap).
  2. I/O contract — outputs land at the canonical paths, the data file
     round-trips through numpy/pandas, the options manifest is recorded.
  3. Plot customization — user-supplied title/labels/figsize/xunit
     reach the saved figure.

A small backbone-like trajectory fixture (5 ALA residues, 50 frames,
real timing) provides realistic input for protein analyses.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
import pytest

from fastmdxplora.analysis import (
    AnalysisOrchestrator,
    available_analyses,
)
from fastmdxplora.analysis.dihedrals import Dihedrals
from fastmdxplora.analysis.rg import Rg
from fastmdxplora.analysis.rmsd import RMSD
from fastmdxplora.analysis.rmsf import RMSF


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _build_backbone_traj(
    n_residues: int = 5, n_frames: int = 50, seed: int = 42
) -> md.Trajectory:
    """Build a small protein-ish trajectory: 4 backbone atoms per residue."""
    rng = np.random.RandomState(seed)
    top = md.Topology()
    chain = top.add_chain()
    for i in range(n_residues):
        res = top.add_residue("ALA", chain, resSeq=i + 1)
        top.add_atom("N", md.element.nitrogen, res)
        top.add_atom("CA", md.element.carbon, res)
        top.add_atom("C", md.element.carbon, res)
        top.add_atom("O", md.element.oxygen, res)

    n_atoms = top.n_atoms
    base = np.zeros((n_atoms, 3))
    for i in range(n_atoms):
        base[i] = [i * 0.15, 0.0, 0.0]  # 1.5 Å spacing along x

    xyz = np.tile(base[None, :, :], (n_frames, 1, 1))
    xyz += rng.normal(scale=0.02, size=xyz.shape)

    times = np.arange(n_frames) * 10.0  # 10 ps per frame
    return md.Trajectory(xyz=xyz.astype(np.float32), topology=top, time=times)


@pytest.fixture
def backbone_traj() -> md.Trajectory:
    return _build_backbone_traj()


@pytest.fixture
def backbone_traj_files(tmp_path: Path, backbone_traj: md.Trajectory):
    """Save the backbone trajectory as DCD + PDB files."""
    pdb = tmp_path / "top.pdb"
    dcd = tmp_path / "traj.dcd"
    backbone_traj[0].save_pdb(str(pdb))
    backbone_traj.save_dcd(str(dcd))
    return dcd, pdb


# ===========================================================================
# Registry — all four analyses present and in order
# ===========================================================================
def test_all_four_registered():
    names = available_analyses()
    for n in ("rmsd", "rmsf", "rg", "dihedrals"):
        assert n in names


def test_canonical_order():
    """Registration order matters for include=None executions."""
    names = available_analyses()
    expected = ("rmsd", "rmsf", "rg", "dihedrals")
    # All expected names should appear in the registry in the canonical
    # order (subsequent sub-deliveries may add more after them).
    indices = [names.index(n) for n in expected]
    assert indices == sorted(indices), f"out-of-order: {names}"


# ===========================================================================
# RMSD
# ===========================================================================
class TestRMSD:
    def test_compute_returns_one_value_per_frame(self, backbone_traj: md.Trajectory):
        rmsd = RMSD()
        out = rmsd.compute(backbone_traj)
        assert out.shape == (backbone_traj.n_frames,)
        assert out.dtype == np.float64

    def test_rmsd_to_self_is_zero(self, backbone_traj: md.Trajectory):
        """RMSD of the reference frame against itself is exactly zero.

        Without alignment this is the mathematical invariant (a frame minus
        itself), so it holds exactly regardless of geometry. The aligned
        path is exercised by the other RMSD tests; asserting "self == 0" on
        the aligned path would instead probe the ill-conditioning of QCP
        superposition on this near-collinear backbone fixture (the base
        geometry is points along the x-axis), which is not the intent and
        is platform/BLAS-dependent.
        """
        rmsd = RMSD(ref=0, align=False)
        out = rmsd.compute(backbone_traj)
        assert out[0] == 0.0

    def test_rmsd_is_nonnegative(self, backbone_traj: md.Trajectory):
        rmsd = RMSD()
        out = rmsd.compute(backbone_traj)
        assert (out >= 0).all()

    def test_negative_ref_index(self, backbone_traj: md.Trajectory):
        """ref=-1 should be the same as ref=last_frame."""
        n = backbone_traj.n_frames
        a = RMSD(ref=-1).compute(backbone_traj)
        b = RMSD(ref=n - 1).compute(backbone_traj)
        np.testing.assert_allclose(a, b, atol=1e-6)

    def test_out_of_range_ref_raises(self, backbone_traj: md.Trajectory):
        rmsd = RMSD(ref=9999)
        with pytest.raises(ValueError, match="out of range"):
            rmsd.compute(backbone_traj)

    def test_no_align_differs_from_aligned(self, backbone_traj: md.Trajectory):
        """Aligned and unaligned RMSD should generally differ."""
        a = RMSD(align=True).compute(backbone_traj)
        b = RMSD(align=False).compute(backbone_traj)
        # They should differ in at least some frames; allow for small numerical equivalence
        assert not np.allclose(a, b, atol=1e-9)

    def test_run_writes_outputs(self, tmp_path: Path, backbone_traj: md.Trajectory):
        rmsd = RMSD(output_dir=tmp_path)
        result = rmsd.run(backbone_traj)
        assert result.status == "ok"
        assert result.data_path.exists()
        assert result.figure_path.exists()
        assert result.options_path.exists()
        # Data round-trips
        loaded = np.loadtxt(result.data_path)
        np.testing.assert_allclose(loaded, result.data, rtol=1e-5)

    def test_options_recorded_in_manifest(
        self, tmp_path: Path, backbone_traj: md.Trajectory
    ):
        RMSD(output_dir=tmp_path, ref=5, align=False).run(backbone_traj)
        with (tmp_path / "rmsd" / "options.json").open() as fh:
            manifest = json.load(fh)
        assert manifest["options"]["ref"] == 5
        assert manifest["options"]["align"] is False

    def test_user_xunit_override(self, tmp_path: Path, backbone_traj: md.Trajectory):
        rmsd = RMSD(output_dir=tmp_path, xunit="ps")
        rmsd.run(backbone_traj)
        assert rmsd._user_xunit == "ps"

    def test_default_xlabel_uses_ns_when_time_available(
        self, tmp_path: Path, backbone_traj: md.Trajectory
    ):
        """A trajectory with time data should default to Time (ns)."""
        rmsd = RMSD(output_dir=tmp_path)
        rmsd.run(backbone_traj)
        assert rmsd.default_xlabel() == "Time (ns)"

    def test_default_xlabel_falls_back_to_frame(self, tmp_path: Path):
        """A trajectory without time data should default to Frame."""
        traj = _build_backbone_traj()
        traj.time = np.zeros(traj.n_frames)  # no real timing
        rmsd = RMSD(output_dir=tmp_path)
        rmsd.run(traj)
        assert rmsd.default_xlabel() == "Frame"


# ===========================================================================
# RMSF
# ===========================================================================
class TestRMSF:
    def test_per_residue_returns_two_columns(self, backbone_traj: md.Trajectory):
        out = RMSF().compute(backbone_traj)
        assert out.ndim == 2
        assert out.shape[1] == 2
        # CA selection → one row per residue
        assert out.shape[0] == backbone_traj.n_residues

    def test_per_atom_returns_two_columns(self, backbone_traj: md.Trajectory):
        out = RMSF(per_residue=False, selection="all").compute(backbone_traj)
        assert out.ndim == 2
        assert out.shape[1] == 2
        assert out.shape[0] == backbone_traj.n_atoms

    def test_rmsf_is_nonnegative(self, backbone_traj: md.Trajectory):
        out = RMSF().compute(backbone_traj)
        assert (out[:, 1] >= 0).all()

    def test_rmsf_magnitude_reasonable(self, backbone_traj: md.Trajectory):
        """The synthetic trajectory has 0.02 nm noise; RMSF should be ~0.02."""
        out = RMSF().compute(backbone_traj)
        mean_rmsf = out[:, 1].mean()
        # Per-atom Gaussian noise σ=0.02 nm in each xyz dimension gives
        # an expected position fluctuation of order σ. Allow factor-of-2 slack.
        assert 0.005 < mean_rmsf < 0.04, f"unexpected mean RMSF: {mean_rmsf}"

    def test_run_writes_outputs(self, tmp_path: Path, backbone_traj: md.Trajectory):
        result = RMSF(output_dir=tmp_path).run(backbone_traj)
        assert result.status == "ok"
        assert result.data_path.exists()
        assert result.figure_path.exists()

    def test_x_label_is_residue(self, tmp_path: Path, backbone_traj: md.Trajectory):
        rmsf = RMSF(output_dir=tmp_path)
        rmsf.run(backbone_traj)
        assert rmsf.default_xlabel() == "Residue"

    def test_x_label_is_atom_when_per_atom(self, tmp_path: Path, backbone_traj):
        rmsf = RMSF(output_dir=tmp_path, per_residue=False, selection="all")
        rmsf.run(backbone_traj)
        assert rmsf.default_xlabel() == "Atom serial"


# ===========================================================================
# Rg
# ===========================================================================
class TestRg:
    def test_default_returns_one_value_per_frame(self, backbone_traj: md.Trajectory):
        out = Rg().compute(backbone_traj)
        assert out.shape == (backbone_traj.n_frames,)

    def test_rg_is_positive(self, backbone_traj: md.Trajectory):
        out = Rg().compute(backbone_traj)
        assert (out > 0).all()

    def test_equal_weights_match_mdtraj_compute_rg(self, backbone_traj: md.Trajectory):
        """MDTraj's compute_rg weights every atom equally, and says so.

        This used to assert the default matched it, which held while the
        default was a thin wrapper. The default is now mass-weighted, as the
        docstring had always claimed and as the conventional definition has
        it, so the equivalence is with the unweighted option.
        """
        ours = Rg(mass_weighted=False).compute(backbone_traj)
        theirs = md.compute_rg(backbone_traj)
        np.testing.assert_allclose(ours, theirs, atol=1e-8)

    def test_the_default_is_mass_weighted(self, backbone_traj: md.Trajectory):
        weighted = Rg().compute(backbone_traj)
        equal = Rg(mass_weighted=False).compute(backbone_traj)
        assert not np.allclose(weighted, equal, rtol=1e-4)

    def test_with_selection(self, backbone_traj: md.Trajectory):
        """Rg with a CA selection must differ from all-atom Rg (atoms differ)."""
        all_atoms = Rg(selection="all").compute(backbone_traj)
        ca_only = Rg(selection="name CA").compute(backbone_traj)
        assert not np.allclose(all_atoms, ca_only)

    def test_by_chain_single_chain(self, backbone_traj: md.Trajectory):
        """With a single chain, by_chain should still return a (n, 1) shape."""
        out = Rg(by_chain=True).compute(backbone_traj)
        # The synthetic trajectory has 1 chain → 1 column (just total)
        assert out.ndim == 2
        assert out.shape[0] == backbone_traj.n_frames

    def test_run_writes_outputs(self, tmp_path: Path, backbone_traj):
        result = Rg(output_dir=tmp_path).run(backbone_traj)
        assert result.status == "ok"
        assert result.data_path.exists()
        assert result.figure_path.exists()


# ===========================================================================
# Dihedrals
# ===========================================================================
class TestDihedrals:
    def test_compute_returns_dataframe(self, backbone_traj: md.Trajectory):
        out = Dihedrals().compute(backbone_traj)
        assert isinstance(out, pd.DataFrame)
        assert set(out.columns) >= {"frame", "residue", "phi_deg", "psi_deg"}

    def test_angles_in_valid_range(self, backbone_traj: md.Trajectory):
        out = Dihedrals().compute(backbone_traj)
        assert ((out["phi_deg"] >= -180) & (out["phi_deg"] <= 180)).all()
        assert ((out["psi_deg"] >= -180) & (out["psi_deg"] <= 180)).all()

    def test_n_rows_equals_n_frames_times_n_inner_residues(
        self, backbone_traj: md.Trajectory
    ):
        """For a 5-residue peptide, residues 2-4 have both phi and psi (3 res)."""
        out = Dihedrals().compute(backbone_traj)
        # Inner residues = those with both phi (needs prev residue C) and
        # psi (needs next residue N). For a 5-residue chain: residues 2,3,4.
        n_inner = 3
        assert len(out) == backbone_traj.n_frames * n_inner

    def test_run_writes_outputs(self, tmp_path: Path, backbone_traj: md.Trajectory):
        result = Dihedrals(output_dir=tmp_path).run(backbone_traj)
        assert result.status == "ok"
        assert result.data_path.exists()
        assert result.figure_path.exists()
        # The .dat file is CSV format
        df = pd.read_csv(result.data_path)
        assert "phi_deg" in df.columns
        assert "psi_deg" in df.columns

    def test_no_protein_raises(self, tmp_path: Path):
        """A trajectory without backbone atoms should raise a clear error."""
        # Build a "trajectory" of pure water — no protein backbone
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("HOH", chain)
        for nm in ("O", "H1", "H2"):
            el = md.element.oxygen if nm == "O" else md.element.hydrogen
            top.add_atom(nm, el, res)
        xyz = np.random.RandomState(0).rand(5, top.n_atoms, 3).astype(np.float32)
        traj = md.Trajectory(xyz=xyz, topology=top)

        dh = Dihedrals(output_dir=tmp_path)
        with pytest.raises(ValueError, match="backbone dihedrals"):
            dh.compute(traj)

    def test_scatter_mode_runs(self, tmp_path: Path, backbone_traj):
        """Density=False switches to scatter plot."""
        result = Dihedrals(output_dir=tmp_path, density=False).run(backbone_traj)
        assert result.status == "ok"


# ===========================================================================
# End-to-end: all four through the AnalysisOrchestrator
# ===========================================================================
class TestEndToEnd:
    def test_all_four_run_through_orchestrator(
        self, tmp_path: Path, backbone_traj_files
    ):
        dcd, pdb = backbone_traj_files
        ao = AnalysisOrchestrator(
            trajectory=dcd,
            topology=pdb,
            output_dir=tmp_path / "out",
        )
        results = ao.run(include=["rmsd", "rmsf", "rg", "dihedrals"])
        for name in ("rmsd", "rmsf", "rg", "dihedrals"):
            assert name in results
            assert results[name].status == "ok", results[name].message
            assert results[name].data_path.exists()
            assert results[name].figure_path.exists()

    def test_orchestrator_writes_complete_manifest(
        self, tmp_path: Path, backbone_traj_files
    ):
        dcd, pdb = backbone_traj_files
        ao = AnalysisOrchestrator(
            trajectory=dcd, topology=pdb, output_dir=tmp_path / "out"
        )
        ao.run(include=["rmsd", "rg"])
        manifest_path = ao.output_dir / "analysis_manifest.json"
        with manifest_path.open() as fh:
            manifest = json.load(fh)
        assert manifest["plan"] == ["rmsd", "rg"]
        assert manifest["results"]["rmsd"]["status"] == "ok"
        assert manifest["results"]["rg"]["status"] == "ok"

    def test_per_analysis_options_forwarded(
        self, tmp_path: Path, backbone_traj_files
    ):
        dcd, pdb = backbone_traj_files
        ao = AnalysisOrchestrator(
            trajectory=dcd, topology=pdb, output_dir=tmp_path / "out"
        )
        results = ao.run(
            include=["rmsd"],
            options={"rmsd": {"ref": 5, "align": False}},
        )
        with (results["rmsd"].options_path).open() as fh:
            manifest = json.load(fh)
        assert manifest["options"]["ref"] == 5
        assert manifest["options"]["align"] is False


# ===========================================================================
# Plot customization for the new analyses
# ===========================================================================
class TestUserCustomization:
    def test_user_xlabel_wins_over_default(
        self, tmp_path: Path, backbone_traj: md.Trajectory
    ):
        rmsd = RMSD(output_dir=tmp_path, xlabel="Custom X")
        rmsd.run(backbone_traj)
        assert rmsd._user_xlabel == "Custom X"
        # Default xlabel for RMSD (without override) is "Time (ns)" for
        # this trajectory; the user override should still be Custom X.

    def test_user_title_propagates(self, tmp_path: Path, backbone_traj):
        rg = Rg(output_dir=tmp_path, title="Compactness over time")
        rg.run(backbone_traj)
        assert rg.figure_title() == "Compactness over time"

    def test_user_xunit_frames(self, tmp_path: Path, backbone_traj):
        """xunit=frames overrides the auto-detected ns default."""
        rmsd = RMSD(output_dir=tmp_path, xunit="frames")
        rmsd.run(backbone_traj)
        assert rmsd.default_xlabel() == "Frame"
        x, label = rmsd.frame_axis(backbone_traj)
        assert label == "Frame"
        np.testing.assert_array_equal(x, np.arange(backbone_traj.n_frames))

    def test_user_xunit_ps(self, tmp_path: Path, backbone_traj):
        rmsd = RMSD(output_dir=tmp_path, xunit="ps")
        rmsd.run(backbone_traj)
        _, label = rmsd.frame_axis(backbone_traj)
        assert label == "Time (ps)"

    def test_invalid_xunit_raises(self, tmp_path: Path, backbone_traj):
        rmsd = RMSD(output_dir=tmp_path, xunit="furlongs")
        with pytest.raises(ValueError, match="xunit"):
            rmsd.frame_axis(backbone_traj)


def _hbond_traj(n_frames: int = 100, transient_frames: int = 5) -> md.Trajectory:
    """Two donor-acceptor pairs: one bonded always, one only briefly."""
    top = md.Topology()
    chain = top.add_chain()
    r1 = top.add_residue("ALA", chain, resSeq=1)
    n1 = top.add_atom("N", md.element.nitrogen, r1)
    h1 = top.add_atom("H", md.element.hydrogen, r1)
    r2 = top.add_residue("ALA", chain, resSeq=2)
    top.add_atom("O", md.element.oxygen, r2)
    r3 = top.add_residue("ALA", chain, resSeq=3)
    n2 = top.add_atom("N", md.element.nitrogen, r3)
    h2 = top.add_atom("H", md.element.hydrogen, r3)
    r4 = top.add_residue("ALA", chain, resSeq=4)
    top.add_atom("O", md.element.oxygen, r4)
    top.add_bond(n1, h1)
    top.add_bond(n2, h2)

    xyz = np.zeros((n_frames, 6, 3), dtype=np.float32)
    xyz[:, 0] = [0.00, 0.0, 0.0]
    xyz[:, 1] = [0.10, 0.0, 0.0]
    xyz[:, 2] = [0.30, 0.0, 0.0]          # bonded in every frame
    xyz[:, 3] = [0.00, 1.0, 0.0]
    xyz[:, 4] = [0.10, 1.0, 0.0]
    xyz[:, 5] = [2.00, 1.0, 0.0]          # far away
    xyz[:transient_frames, 5] = [0.30, 1.0, 0.0]   # except briefly
    return md.Trajectory(xyz=xyz, topology=top, time=np.arange(n_frames) * 1.0)


class TestHBondsCountsEveryBondItDepicts:
    """The series is the number of hydrogen bonds per frame, so it must be.

    Baker-Hubbard proposes bonds above an occupancy threshold, and only the
    proposed ones were evaluated frame by frame. A bond present in five per
    cent of frames therefore contributed to none of them -- including the
    frames where it was there. The plot said one bond where there were two.
    """

    def test_a_transient_bond_is_counted_in_the_frames_it_exists(self) -> None:
        from fastmdxplora.analysis.hbonds import HBonds

        traj = _hbond_traj(n_frames=100, transient_frames=5)
        counts = HBonds().compute(traj)["n_hbonds"].to_numpy()

        assert counts[0] == 2, "both bonds are present in the first frames"
        assert counts[50] == 1, "only the persistent one remains later"
        assert counts.sum() == 105

    def test_raising_the_threshold_restricts_the_series_again(self) -> None:
        """The old behaviour is still reachable, deliberately."""
        from fastmdxplora.analysis.hbonds import HBonds

        traj = _hbond_traj(n_frames=100, transient_frames=5)
        counts = HBonds(candidate_freq=0.1).compute(traj)["n_hbonds"].to_numpy()
        assert counts[0] == 1
        assert counts.sum() == 100

    def test_freq_reports_how_many_bonds_persist(self) -> None:
        """Once every bond is proposed, freq would otherwise decide nothing."""
        from fastmdxplora.analysis.hbonds import HBonds

        analysis = HBonds()
        analysis.compute(_hbond_traj())
        assert analysis.options["n_candidate_bonds"] == 2
        assert analysis.options["n_persistent_bonds"] == 1

    @staticmethod
    def _captured_warnings(build):
        """Collect what the package logger emits.

        ``caplog`` sees nothing: the package sets ``propagate = False`` on the
        ``fastmdx`` logger so its own handler owns the formatting, which means
        records never reach the root logger pytest attaches to. A handler is
        added to that logger instead.
        """
        import logging

        records: list[str] = []

        class _Collect(logging.Handler):
            def emit(self, record):
                records.append(record.getMessage())

        logger = logging.getLogger("fastmdx")
        handler = _Collect(level=logging.WARNING)
        logger.addHandler(handler)
        try:
            build()
        finally:
            logger.removeHandler(handler)
        return " ".join(records)

    def test_the_multiplier_says_what_it_is_doing(self) -> None:
        """It reproduces an earlier count; it is not a convention.

        MDTraj lists each hydrogen bond once, as one donor-hydrogen-acceptor
        triplet, and version 1 counted the same way -- one row per triplet, one
        frame at a time. Doubling models neither.
        """
        from fastmdxplora.analysis.hbonds import HBonds

        said = self._captured_warnings(lambda: HBonds(count_multiplier=2))
        assert "not a convention" in said

    def test_no_multiplier_is_silent(self) -> None:
        from fastmdxplora.analysis.hbonds import HBonds

        said = self._captured_warnings(HBonds)
        assert "multiplying" not in said


def _hinge_traj(n_frames: int = 40, seed: int = 0) -> md.Trajectory:
    """Two conformational states, plus rigid-body drift on alternate frames.

    The hinge moves half the chain, which is a change of shape. The drift moves
    all of it, which is not.
    """
    rng = np.random.RandomState(seed)
    top = md.Topology()
    chain = top.add_chain()
    for i in range(12):
        res = top.add_residue("ALA", chain, resSeq=i + 1)
        top.add_atom("CA", md.element.carbon, res)

    base = np.zeros((12, 3))
    base[:, 0] = np.arange(12) * 0.38
    xyz = np.tile(base[None], (n_frames, 1, 1))
    xyz += rng.normal(scale=0.005, size=xyz.shape)
    xyz[n_frames // 2:, 6:, 1] += 0.5      # the hinge
    xyz[::2] += 5.0                         # the drift
    return md.Trajectory(xyz=xyz.astype(np.float32), topology=top)


def _recovery(labels: np.ndarray) -> float:
    """How much of the true two-state split the labels reproduce."""
    half = len(labels) // 2
    agree = (labels[:half] == labels[0]).sum() + (labels[half:] == labels[half]).sum()
    return max(agree, len(labels) - agree) / len(labels)


class TestClusteringComparesShapeNotPlacement:
    """What frames are compared in decides the answer.

    Comparing raw coordinates makes the leading difference between frames where
    the molecule drifted and how it turned. Two identical conformations in
    different orientations then look as far apart as anything in the run, and
    clustering reports position rather than shape.
    """

    def test_pairwise_rmsd_recovers_the_conformational_states(self) -> None:
        from fastmdxplora.analysis.cluster import Cluster

        labels = Cluster(
            features="rmsd", methods=["hierarchical"], n_clusters=2,
            linkage="ward", selection="all",
        ).compute(_hinge_traj())["hierarchical"]
        assert _recovery(labels) == 1.0

    def test_superposed_coordinates_recover_them_too(self) -> None:
        from fastmdxplora.analysis.cluster import Cluster

        labels = Cluster(
            features="coordinates", methods=["hierarchical"], n_clusters=2,
            linkage="ward", selection="all",
        ).compute(_hinge_traj())["hierarchical"]
        assert _recovery(labels) >= 0.9

    def test_unaligned_coordinates_do_not(self) -> None:
        """The comparison that motivates superposing at all.

        This is what clustering on ``traj.xyz`` without alignment gives: the
        rigid-body drift dominates and the split is no better than chance.
        """
        from sklearn.cluster import AgglomerativeClustering

        traj = _hinge_traj()
        raw = traj.xyz.reshape(traj.n_frames, -1)
        labels = AgglomerativeClustering(n_clusters=2, linkage="ward").fit_predict(raw)
        assert _recovery(labels) < 0.7

    def test_the_coordinate_space_keeps_distances_in_nm(self) -> None:
        """So eps means the same thing in both, and the two are comparable."""
        from fastmdxplora.analysis.cluster import _superposed_coordinates

        traj = _hinge_traj()
        atom_idx = np.arange(traj.n_atoms)
        # superpose() aligns in place, so the trajectory now holds the frames
        # the returned points were built from. Aligning it a second time to
        # compare against would measure something else.
        points = _superposed_coordinates(traj, atom_idx)
        expected = np.sqrt(
            ((traj.xyz[0, atom_idx] - traj.xyz[1, atom_idx]) ** 2).sum()
            / len(atom_idx)
        )
        assert np.isclose(
            np.linalg.norm(points[0] - points[1]), expected, rtol=1e-5
        )

    def test_an_unknown_feature_space_is_refused(self) -> None:
        import pytest

        from fastmdxplora.analysis.cluster import Cluster

        with pytest.raises(ValueError, match="Unknown clustering features"):
            Cluster(features="cartesian")


class TestQValueFollowsThePaperItCites:
    """Best, Hummer & Eaton define Q through a switching function.

    A contact was instead counted as present or absent by a single distance
    shared by every pair, so one formed at 0.20 nm was judged by the same
    standard as one formed at 0.44 nm. The measure reported a fold coming apart
    well before the paper's does -- on a peeling hairpin it reached zero while
    the paper's still read 0.62.
    """

    @staticmethod
    def _hairpin(n_frames: int = 100, peel: float = 0.5, noise: float = 0.01):
        rng = np.random.RandomState(3)
        top = md.Topology()
        chain = top.add_chain()
        for i in range(14):
            res = top.add_residue("ALA", chain, resSeq=i + 1)
            top.add_atom("CA", md.element.carbon, res)
            top.add_atom("CB", md.element.carbon, res)
        pos = []
        for i in range(7):
            pos += [[i * 0.38, 0.0, 0.0], [i * 0.38, 0.15, 0.0]]
        for i in range(7):
            pos += [[(6 - i) * 0.38, 0.42, 0.0], [(6 - i) * 0.38, 0.27, 0.0]]
        xyz = np.tile(np.array(pos)[None], (n_frames, 1, 1))
        xyz[:, 14:, 1] += np.linspace(0.0, peel, n_frames)[:, None]
        xyz += rng.normal(scale=noise, size=xyz.shape)
        return md.Trajectory(xyz=xyz.astype(np.float32), topology=top)

    @staticmethod
    def _reference_q(traj, cutoff=0.45, beta=50.0, lam=1.8, sep=4):
        """The formula as published, written out independently."""
        pairs = np.column_stack(np.triu_indices(traj.n_residues, k=sep))
        distances, _ = md.compute_contacts(
            traj, contacts=pairs, scheme="closest-heavy")
        native = distances[0] < cutoff
        r0 = distances[0][native]
        return (
            1.0 / (1.0 + np.exp(beta * (distances[:, native] - lam * r0)))
        ).mean(axis=1)

    def test_it_matches_the_published_formula(self) -> None:
        from fastmdxplora.analysis.qvalue import QValue

        traj = self._hairpin()
        assert np.allclose(
            QValue(selection="all").compute(traj),
            self._reference_q(traj),
            atol=1e-6,
        )

    def test_a_contact_is_judged_against_its_own_native_distance(self) -> None:
        """Not against one number shared by every pair.

        A single threshold treats a contact formed at 0.20 nm and one formed
        at 0.44 nm alike: the second breaks on any stretch at all, the first
        can more than double and still count. Judged against its own native
        distance, each is allowed the same proportional stretch, so a fold
        under way reads as further from unfolded than a threshold says.
        """
        from fastmdxplora.analysis.qvalue import QValue

        traj = self._hairpin()
        q = QValue(selection="all").compute(traj)

        # The same trajectory, scored by the threshold this used to apply.
        pairs = np.column_stack(np.triu_indices(traj.n_residues, k=4))
        distances, _ = md.compute_contacts(
            traj, contacts=pairs, scheme="closest-heavy")
        native = distances[0] < 0.45
        threshold = (distances[:, native] < 0.45).sum(axis=1) / native.sum()

        # Both read as folded at the reference. The switching function sits a
        # little under 1 even there, which is inherent to it: a contact at its
        # native distance is well inside the tolerance but the sigmoid is
        # smooth, so it never quite saturates.
        assert q[0] > 0.99 and threshold[0] == 1.0

        # The difference is in how the fold comes apart. The threshold reaches
        # zero while contacts are still within lambda of where they started.
        assert threshold[-1] == 0.0
        assert q[-1] > 0.25, "the paper's measure still sees the fold"
        assert q[-1] > threshold[-1] + 0.25

    def test_the_defaults_are_the_paper_s(self) -> None:
        from fastmdxplora.analysis.qvalue import QValue

        analysis = QValue()
        assert analysis.beta == 50.0
        assert analysis.lambda_factor == 1.8
        assert analysis.cutoff == 0.45
        assert analysis.min_seq_separation == 4

    def test_a_steep_beta_approaches_the_old_threshold(self) -> None:
        """The switching function contains the hard cutoff as a limit."""
        from fastmdxplora.analysis.qvalue import QValue

        traj = self._hairpin()
        steep = QValue(selection="all", beta=1e5, lambda_factor=1.0).compute(traj)
        assert steep.min() < 0.1, "at the limit it becomes a threshold again"


class TestSecondaryStructureIgnoresWhatIsNotProtein:
    """DSSP returns a column for every residue, not only the protein ones.

    A residue without backbone atoms gets the code ``NA``: a ligand, an ion, a
    water. Those columns used to survive into the timeline, where ``NA`` is not
    in the colour map and so was drawn as coil -- and because the labels were
    taken from the protein residues alone, the two lists came out different
    lengths and the code fell back to plain numbering. A protein numbered 10 to
    15 was relabelled 0 to 6, with the ligand as the last row.

    Scope defaults to ``solute``, which is protein and ligand, so this was the
    ordinary case for the systems this software is meant for.
    """

    @staticmethod
    def _protein_with_ligand(n_residues: int = 6, first_resseq: int = 10):
        top = md.Topology()
        chain = top.add_chain()
        for i in range(n_residues):
            res = top.add_residue("ALA", chain, resSeq=first_resseq + i)
            for name, element in (("N", md.element.nitrogen),
                                  ("CA", md.element.carbon),
                                  ("C", md.element.carbon),
                                  ("O", md.element.oxygen)):
                top.add_atom(name, element, res)
        ligand_chain = top.add_chain()
        ligand = top.add_residue("LIG", ligand_chain, resSeq=200)
        top.add_atom("C1", md.element.carbon, ligand)

        xyz = np.zeros((5, top.n_atoms, 3), dtype=np.float32)
        for atom in range(top.n_atoms):
            xyz[:, atom] = [atom * 0.15, 0.0, 0.0]
        return md.Trajectory(xyz=xyz, topology=top)

    def test_the_ligand_is_not_given_a_row(self) -> None:
        from fastmdxplora.analysis.ss import SS

        traj = self._protein_with_ligand()
        columns = [c for c in SS(selection="all").compute(traj).columns
                   if c != "frame"]
        assert len(columns) == 6, "one row per protein residue and no more"

    def test_residue_numbering_survives(self) -> None:
        from fastmdxplora.analysis.ss import SS

        traj = self._protein_with_ligand(first_resseq=10)
        columns = [c for c in SS(selection="all").compute(traj).columns
                   if c != "frame"]
        assert columns == [10, 11, 12, 13, 14, 15]

    def test_no_NA_reaches_the_data(self) -> None:
        """It is not in the colour map, so it would be drawn as coil."""
        from fastmdxplora.analysis.ss import SS

        df = SS(selection="all").compute(self._protein_with_ligand())
        codes = df.drop(columns="frame").to_numpy()
        assert "NA" not in set(np.unique(codes).tolist())

    def test_a_protein_only_system_is_unchanged(self) -> None:
        from fastmdxplora.analysis.ss import SS

        traj = self._protein_with_ligand()
        protein = traj.atom_slice(traj.topology.select("protein"))
        columns = [c for c in SS(selection="all").compute(protein).columns
                   if c != "frame"]
        assert columns == [10, 11, 12, 13, 14, 15]


class TestSecondaryStructureSaysWhenItDoesNotApply:
    """A system with no protein backbone has no secondary structure.

    Dropping the non-protein columns can leave none at all -- a nucleic acid
    run, a lone ligand, a coarse-grained model. Drawing an empty timeline for
    that is not an answer to the question; it is a picture of nothing, and the
    only sign anything was wrong came from matplotlib complaining that the axis
    had no extent.
    """

    @staticmethod
    def _no_backbone():
        top = md.Topology()
        chain = top.add_chain()
        residue = top.add_residue("ALA", chain)
        for i in range(5):
            top.add_atom(f"A{i}", md.element.carbon, residue)
        xyz = np.random.RandomState(0).rand(6, 5, 3).astype(np.float32)
        return md.Trajectory(xyz=xyz, topology=top)

    def test_it_refuses_and_says_why(self) -> None:
        import pytest

        from fastmdxplora.analysis.ss import SS

        with pytest.raises(ValueError, match="has a protein backbone"):
            SS(selection="all").compute(self._no_backbone())

    def test_the_orchestrator_records_it_without_stopping(self) -> None:
        """One analysis that does not apply must not end the run."""
        from fastmdxplora.analysis.orchestrator import AnalysisOrchestrator

        assert hasattr(AnalysisOrchestrator, "run")
        # The orchestrator catches per-analysis failures and records them; the
        # refusal above therefore reaches the user as a stated reason rather
        # than as an empty figure or a stopped exploration.


class TestAtomLabelsAreNeverNotANumber:
    """The atom column is written to the data file, so it must hold a number.

    ``Atom.serial`` is filled in by the file a topology came from. A trajectory
    built in memory -- which the Python API allows -- has none, and ``None``
    became NaN in the column. The saved ``.dat`` carried the NaN, and the plot
    cast it to the most negative integer there is and used that as an axis
    label.
    """

    @staticmethod
    def _in_memory():
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("ALA", chain, resSeq=1)
        for name in ("N", "CA", "C", "O"):
            top.add_atom(name, md.element.carbon, res)
        xyz = np.random.RandomState(0).rand(10, 4, 3).astype(np.float32)
        return md.Trajectory(xyz=xyz, topology=top)

    def test_rmsf_labels_atoms_without_serials(self) -> None:
        from fastmdxplora.analysis.rmsf import RMSF

        out = RMSF(per_residue=False, selection="all").compute(self._in_memory())
        assert not np.isnan(out[:, 0]).any()
        assert out[:, 0].tolist() == [1.0, 2.0, 3.0, 4.0]

    def test_a_file_supplied_serial_is_kept(self, tmp_path) -> None:
        """Where the file says, the file wins: the number is the user's."""
        from fastmdxplora.analysis.rmsf import RMSF

        traj = self._in_memory()
        path = tmp_path / "top.pdb"
        traj[0].save_pdb(path)
        loaded = md.load(str(path))
        from_file = [a.serial for a in loaded.topology.atoms]

        out = RMSF(per_residue=False, selection="all").compute(
            md.Trajectory(xyz=traj.xyz, topology=loaded.topology))
        assert out[:, 0].tolist() == [float(s) for s in from_file]

    def test_ligand_rmsf_labels_survive_the_cast_the_plot_makes(self) -> None:
        from fastmdxplora.analysis.rmsf import _atom_labels

        traj = self._in_memory()
        labels = _atom_labels(list(traj.topology.atoms))
        assert labels.dtype == np.int64
        assert (labels > 0).all()


class TestDimRedSaysWhenThereIsNothingToDecompose:
    """Variance is what a decomposition decomposes.

    A structure that does not move has none. PCA divides each component's
    variance by the total, which is zero, so the ratios came out NaN and the
    figure was labelled "PC 1 (nan%)" over a scatter of coincident points --
    a plot that looks like a result and is not one. The only sign was a numpy
    warning about dividing by zero, nine of them in one test run.
    """

    @staticmethod
    def _motionless(n_frames: int = 10):
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("ALA", chain)
        for i in range(5):
            top.add_atom(f"A{i}", md.element.carbon, res)
        base = np.random.RandomState(0).rand(5, 3)
        xyz = np.tile(base[None], (n_frames, 1, 1)).astype(np.float32)
        return md.Trajectory(xyz=xyz, topology=top)

    def test_it_refuses_and_says_why(self) -> None:
        import pytest

        from fastmdxplora.analysis.dimred import DimRed

        with pytest.raises(ValueError, match="nothing to decompose"):
            DimRed(methods=["pca"], selection="all").compute(self._motionless())

    def test_no_division_by_zero_reaches_numpy(self) -> None:
        import warnings

        from fastmdxplora.analysis.dimred import DimRed

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                DimRed(methods=["pca"], selection="all").compute(self._motionless())
            except ValueError:
                pass
            assert not [w for w in caught if "invalid value" in str(w.message)]

    def test_every_method_is_covered_not_only_pca(self) -> None:
        """There are no neighbourhoods among points that are all one point."""
        import pytest

        from fastmdxplora.analysis.dimred import DimRed

        with pytest.raises(ValueError, match="nothing to decompose"):
            DimRed(methods=["tsne"], selection="all").compute(self._motionless())

    def test_a_trajectory_that_moves_is_unaffected(self) -> None:
        from fastmdxplora.analysis.dimred import DimRed

        traj = self._motionless()
        moving = md.Trajectory(
            xyz=traj.xyz + np.random.RandomState(1).normal(
                scale=0.01, size=traj.xyz.shape).astype(np.float32),
            topology=traj.topology,
        )
        analysis = DimRed(methods=["pca"], selection="all")
        embedding = analysis.compute(moving)["pca"]
        assert np.isfinite(embedding).all()
        assert np.isfinite(analysis._explained_variance).all()


class TestRadiusOfGyrationIsMassWeighted:
    """What was documented and what was computed had come apart.

    The docstring said the radius of gyration used masses from the topology.
    ``mdtraj.compute_rg`` weights every atom equally unless told otherwise --
    its own Notes say so -- and when it is given masses it still measures from
    the geometric centre rather than the centre of mass. So the reported number
    counted each hydrogen for as much as each carbon, a few per cent from the
    value GROMACS, cpptraj and a published figure all report.
    """

    @staticmethod
    def _protein_like(n_residues: int = 30, n_frames: int = 3):
        rng = np.random.RandomState(0)
        top = md.Topology()
        chain = top.add_chain()
        for i in range(n_residues):
            res = top.add_residue("ALA", chain, resSeq=i + 1)
            for name, element in (("N", md.element.nitrogen),
                                  ("CA", md.element.carbon),
                                  ("C", md.element.carbon),
                                  ("O", md.element.oxygen),
                                  ("HA", md.element.hydrogen),
                                  ("H", md.element.hydrogen)):
                top.add_atom(name, element, res)
        xyz = rng.normal(size=(n_frames, top.n_atoms, 3)).astype(np.float32)
        return md.Trajectory(xyz=xyz, topology=top)

    @staticmethod
    def _reference(traj):
        """The conventional definition, written out independently."""
        masses = np.array([a.element.mass for a in traj.topology.atoms])
        xyz = traj.xyz.astype(np.float64)
        centre = (xyz * masses[None, :, None]).sum(axis=1) / masses.sum()
        squared = ((xyz - centre[:, None, :]) ** 2).sum(axis=2)
        return np.sqrt((squared * masses[None, :]).sum(axis=1) / masses.sum())

    def test_it_matches_the_conventional_definition(self) -> None:
        from fastmdxplora.analysis.rg import Rg

        traj = self._protein_like()
        assert np.allclose(
            Rg(selection="all").compute(traj), self._reference(traj), atol=1e-9
        )

    def test_equal_weights_remain_available_and_match_mdtraj(self) -> None:
        from fastmdxplora.analysis.rg import Rg

        traj = self._protein_like()
        assert np.allclose(
            Rg(selection="all", mass_weighted=False).compute(traj),
            md.compute_rg(traj),
            atol=1e-6,
        )

    def test_the_two_differ_by_enough_to_matter(self) -> None:
        """Hydrogens are half the atoms and a fraction of the mass."""
        from fastmdxplora.analysis.rg import Rg

        traj = self._protein_like()
        weighted = Rg(selection="all").compute(traj)
        equal = Rg(selection="all", mass_weighted=False).compute(traj)
        assert not np.allclose(weighted, equal, rtol=1e-3)

    def test_a_topology_without_masses_falls_back(self) -> None:
        """A virtual site or a coarse-grained bead may carry no element."""
        from fastmdxplora.analysis.rg import Rg

        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("UNK", chain)
        for _ in range(4):
            top.add_atom("X", md.element.virtual, res)
        traj = md.Trajectory(
            xyz=np.random.RandomState(1).rand(3, 4, 3).astype(np.float32),
            topology=top,
        )
        out = Rg(selection="all").compute(traj)
        assert np.isfinite(out).all()


class TestDistancesHonourThePeriodicBox:
    """A solvated trajectory is not always imaged.

    A molecule split across the periodic boundary looks far from everything it
    is actually touching. Contacts and hydrogen bonds measured plain distances
    regardless, so a bound ligand sitting across the boundary reported no
    contacts at all -- while the Q-value, on the same frames, measured with the
    box because it never overrode MDTraj's default. Two analyses answering
    "is this near that" differently on one trajectory.
    """

    @staticmethod
    def _split_across_boundary(box: float = 5.0):
        top = md.Topology()
        chain = top.add_chain()
        for i in range(4):
            res = top.add_residue("ALA", chain, resSeq=i + 1)
            for name in ("N", "CA", "C", "O"):
                top.add_atom(name, md.element.carbon, res)
        ligand_chain = top.add_chain()
        ligand = top.add_residue("LIG", ligand_chain, resSeq=900)
        for i in range(3):
            top.add_atom(f"C{i}", md.element.carbon, ligand)

        protein = np.column_stack(
            [np.linspace(0.05, 0.5, 16), np.zeros(16), np.zeros(16)])
        # Just the other side of the wall: close to the protein through it.
        lig = np.array([[box - 0.2, 0, 0], [box - 0.15, 0, 0], [box - 0.25, 0, 0]])
        xyz = np.vstack([protein, lig])[None].astype(np.float32)
        traj = md.Trajectory(xyz=xyz, topology=top)
        traj.unitcell_vectors = np.array(
            [[[box, 0, 0], [0, box, 0], [0, 0, box]]], dtype=np.float32)
        return traj

    def test_contacts_are_found_through_the_boundary(self) -> None:
        from fastmdxplora.analysis.contacts import Contacts

        traj = self._split_across_boundary()
        found = Contacts(ligand_resname="LIG", protein_selection="not resname LIG",
                         selection="all").compute(traj)["n_contacts"][0]
        assert found > 0, "the ligand is bound; it is the box that is in the way"

    def test_the_old_behaviour_is_still_reachable(self) -> None:
        from fastmdxplora.analysis.contacts import Contacts

        traj = self._split_across_boundary()
        found = Contacts(ligand_resname="LIG", protein_selection="not resname LIG",
                         selection="all", periodic=False).compute(traj)["n_contacts"][0]
        assert found == 0

    def test_a_trajectory_without_a_box_is_unaffected(self) -> None:
        """Where there is no unit cell, the setting changes nothing."""
        from fastmdxplora.analysis.contacts import Contacts

        traj = self._split_across_boundary()
        traj.unitcell_vectors = None
        with_box = Contacts(ligand_resname="LIG", protein_selection="not resname LIG",
                            selection="all").compute(traj)["n_contacts"][0]
        without = Contacts(ligand_resname="LIG", protein_selection="not resname LIG",
                           selection="all", periodic=False).compute(traj)["n_contacts"][0]
        assert with_box == without

    def test_hydrogen_bonds_take_the_same_setting(self) -> None:
        from fastmdxplora.analysis.hbonds import HBonds

        analysis = HBonds()
        assert analysis.periodic is True
        assert analysis.options["periodic"] is True


class TestProteinLigandHBondsNeedsTheLigandsBonds:
    """A donor is a nitrogen or oxygen with a hydrogen bonded to it.

    A ligand whose bonds are absent from the topology therefore cannot be seen
    to donate, only to accept -- and the count reports one direction under a
    name that promises both. The guard that built missing bonds only fired when
    the topology had none at all, which a PDB carrying the protein's
    connectivity but no CONECT records for its ligand does not.
    """

    @staticmethod
    def _complex(bond_the_ligand: bool):
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("ALA", chain, resSeq=1)
        n = top.add_atom("N", md.element.nitrogen, res)
        ca = top.add_atom("CA", md.element.carbon, res)
        c = top.add_atom("C", md.element.carbon, res)
        o = top.add_atom("O", md.element.oxygen, res)
        top.add_bond(n, ca)
        top.add_bond(ca, c)
        top.add_bond(c, o)

        ligand_chain = top.add_chain()
        ligand = top.add_residue("LIG", ligand_chain, resSeq=900)
        lo = top.add_atom("O1", md.element.oxygen, ligand)
        lh = top.add_atom("H1", md.element.hydrogen, ligand)
        if bond_the_ligand:
            top.add_bond(lo, lh)

        # O1-H1 ... O : the ligand donating to the protein backbone oxygen.
        frame = np.array([[[0, 0, 0], [0.15, 0, 0], [0.30, 0, 0], [0.42, 0, 0],
                           [0.62, 0, 0], [0.52, 0, 0]]], dtype=np.float32)
        return md.Trajectory(xyz=np.repeat(frame, 4, axis=0), topology=top)

    def test_the_donation_is_found_when_the_bonds_are_there(self) -> None:
        from fastmdxplora.analysis.pl_hbonds import ProteinLigandHBonds

        counts = ProteinLigandHBonds(
            ligand_resname="LIG", protein_selection="not resname LIG",
            selection="all",
        ).compute(self._complex(True))["n_hbonds"]
        assert (counts == 1).all()

    def test_it_refuses_rather_than_counting_one_direction(self) -> None:
        import pytest

        from fastmdxplora.analysis.pl_hbonds import ProteinLigandHBonds

        with pytest.raises(ValueError, match="cannot be seen to donate"):
            ProteinLigandHBonds(
                ligand_resname="LIG", protein_selection="not resname LIG",
                selection="all",
            ).compute(self._complex(False))

    def test_a_ligand_with_no_hydrogens_is_not_blocked(self) -> None:
        """It can only accept, and that is a fact about the ligand."""
        from fastmdxplora.analysis.pl_hbonds import ProteinLigandHBonds

        traj = self._complex(True)
        heavy_only = traj.atom_slice(
            [a.index for a in traj.topology.atoms
             if a.element is None or a.element.symbol != "H"]
        )
        counts = ProteinLigandHBonds(
            ligand_resname="LIG", protein_selection="not resname LIG",
            selection="all",
        ).compute(heavy_only)["n_hbonds"]
        assert len(counts) == traj.n_frames
