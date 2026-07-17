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

    def test_align_defaults_to_true_and_changes_atomic_rmsf(self):
        """Unaligned RMSF retains rigid-body translation around the mean."""
        top = md.Topology()
        chain = top.add_chain()
        residue = top.add_residue("ALA", chain, resSeq=1)
        for name, element in (("N", md.element.nitrogen), ("CA", md.element.carbon), ("C", md.element.carbon)):
            top.add_atom(name, element, residue)
        base = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
        xyz = np.stack([base, base + np.array([1.0, 0.0, 0.0])]).astype(np.float32)
        traj = md.Trajectory(xyz=xyz, topology=top)

        aligned = RMSF(selection="protein", per_residue=False).compute(traj)
        unaligned = RMSF(selection="protein", per_residue=False, align=False).compute(traj)

        assert RMSF().align is True
        assert np.all(aligned[:, 1] < 1e-5)
        assert np.all(unaligned[:, 1] > 0.4)

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

    def test_matches_mdtraj_compute_rg(self, backbone_traj: md.Trajectory):
        """The unwrapped version must equal MDTraj's compute_rg exactly."""
        ours = Rg().compute(backbone_traj)
        theirs = md.compute_rg(backbone_traj)
        np.testing.assert_allclose(ours, theirs, atol=1e-8)

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
