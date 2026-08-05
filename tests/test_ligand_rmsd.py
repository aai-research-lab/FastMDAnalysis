"""Tests for ligand pose RMSD and the ligand-aware orchestrator routing.

LigandRMSD aligns each frame on the protein, then measures the RMSD of the
ligand atoms — the headline "does the ligand stay in its pose" metric. The
orchestrator runs ligand-only analyses automatically when a ligand is present
(detected from the setup manifest) and skips them otherwise; include/exclude
still work.
"""

from __future__ import annotations

import mdtraj as md
import numpy as np
import pytest

from fastmdxplora.analysis.ligand_rmsd import LigandRMSD
from fastmdxplora.analysis.orchestrator import AnalysisOrchestrator, _REGISTRY


@pytest.fixture
def protein_ligand_traj():
    """4-frame trajectory: 6-residue protein + a 3-atom LIG, with the ligand
    shifted by a known amount each frame so RMSD is predictable."""
    lines = []
    ai = 1
    offsets = {"N": (0, 0, 0), "CA": (1.5, 0, 0), "C": (1.5, 1.5, 0), "O": (0, 1.5, 1.0)}
    for ri in range(1, 7):
        bz = ri * 3.8
        for nm, el in [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]:
            dx, dy, dz = offsets[nm]
            lines.append(
                f"ATOM  {ai:5d}  {nm:<3s} ALA A{ri:4d}    "
                f"{ri*4.0+dx:8.3f}{dy:8.3f}{bz+dz:8.3f}  1.00  0.00           {el}"
            )
            ai += 1
    for k, (x, y, z) in enumerate([(10.0, 2.0, 5.0), (11.0, 2.0, 5.0), (10.5, 3.0, 5.0)]):
        lines.append(
            f"HETATM{ai:5d}  C{k+1:<2d} LIG L   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
        ai += 1
    lines.append("END")
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".pdb")
    os.write(fd, ("\n".join(lines) + "\n").encode())
    os.close(fd)
    base = md.load(path)
    lig = base.topology.select("resname LIG")
    frames = [base.xyz[0].copy() for _ in range(4)]
    for f in range(1, 4):
        frames[f][lig] += 0.1 * f
    return md.Trajectory(np.stack(frames), base.topology)


class TestLigandRMSD:
    def test_computes_and_reference_is_zero(self, protein_ligand_traj, tmp_path):
        lr = LigandRMSD(ligand_resname="LIG", output_dir=tmp_path)
        result = lr.compute(protein_ligand_traj)
        assert result.shape == (4,)
        assert result[0] == pytest.approx(0.0, abs=1e-6)

    def test_increases_with_displacement(self, protein_ligand_traj, tmp_path):
        lr = LigandRMSD(ligand_resname="LIG", output_dir=tmp_path)
        result = lr.compute(protein_ligand_traj)
        # Ligand shifted by 0.1*f in x,y,z each frame -> sqrt(3)*0.1*f.
        assert result[1] == pytest.approx(np.sqrt(3) * 0.1, abs=1e-3)
        assert result[3] > result[2] > result[1] > result[0]

    def test_requires_ligand_resname(self, tmp_path):
        with pytest.raises(ValueError, match="requires `ligand_resname`"):
            LigandRMSD(output_dir=tmp_path)

    def test_unknown_resname_raises_on_compute(self, protein_ligand_traj, tmp_path):
        lr = LigandRMSD(ligand_resname="ZZZ", output_dir=tmp_path)
        result = lr.run(protein_ligand_traj)
        assert result.status == "error"
        assert "ZZZ" in result.message

    def test_runs_end_to_end(self, protein_ligand_traj, tmp_path):
        lr = LigandRMSD(ligand_resname="LIG", output_dir=tmp_path)
        result = lr.run(protein_ligand_traj)
        assert result.status == "ok"
        assert (result.output_dir / "ligand_rmsd.dat").exists()
        assert (result.output_dir / "ligand_rmsd.png").exists()


class TestLigandAwareRouting:
    def test_registered_and_marked(self):
        assert "ligand_rmsd" in _REGISTRY
        assert _REGISTRY["ligand_rmsd"].requires_ligand is True

    def _orch(self, tmp_path, ligand_resname):
        # A minimal 2-atom protein PDB is enough to drive plan logic.
        pdb = tmp_path / "po.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        return AnalysisOrchestrator(
            str(pdb), output_dir=tmp_path / "out", ligand_resname=ligand_resname
        )

    def test_skipped_without_ligand(self, tmp_path):
        ao = self._orch(tmp_path, None)
        assert "ligand_rmsd" not in ao._build_plan(None, None)

    def test_included_with_ligand(self, tmp_path):
        ao = self._orch(tmp_path, "LIG")
        assert "ligand_rmsd" in ao._build_plan(None, None)

    def test_explicit_include_honored_without_ligand(self, tmp_path):
        ao = self._orch(tmp_path, None)
        assert ao._build_plan(["ligand_rmsd"], None) == ["ligand_rmsd"]

    def test_exclude_works_with_ligand(self, tmp_path):
        ao = self._orch(tmp_path, "LIG")
        assert "ligand_rmsd" not in ao._build_plan(None, ["ligand_rmsd"])


# ---------------------------------------------------------------------------
# The other ligand-aware analyses: contacts, ligand_rmsf, pl_hbonds
# ---------------------------------------------------------------------------
@pytest.fixture
def complex_in_contact():
    """Compact 3D protein (8 ALA) + a 3-atom LIG within ~0.35 nm of a couple
    residues, 4 frames. Non-degenerate geometry so superpose is well-behaved
    and contacts actually form."""
    top = md.Topology()
    chain = top.add_chain()
    coords = []
    for ri in range(8):
        res = top.add_residue("ALA", chain, resSeq=ri + 1)
        cx, cy, cz = (ri % 3) * 0.4, (ri // 3) * 0.4, (ri % 2) * 0.3
        for nm, el, off in [
            ("N", md.element.nitrogen, (0, 0, 0)),
            ("CA", md.element.carbon, (0.15, 0.05, 0)),
            ("C", md.element.carbon, (0.25, 0.15, 0.05)),
            ("O", md.element.oxygen, (0.2, 0.25, 0.1)),
            ("CB", md.element.carbon, (0.1, 0.1, 0.2)),
        ]:
            top.add_atom(nm, el, res)
            coords.append((cx + off[0], cy + off[1], cz + off[2]))
    ligres = top.add_residue("LIG", chain, resSeq=100)
    for k, off in enumerate([(0.2, 0.1, 0.15), (0.3, 0.15, 0.15), (0.25, 0.05, 0.25)]):
        top.add_atom(f"C{k+1}", md.element.carbon, ligres)
        coords.append(off)
    base = np.array(coords, dtype=np.float32)
    lig = top.select("resname LIG")
    frames = [base.copy() for _ in range(4)]
    for f in range(1, 4):
        frames[f][lig] += 0.02 * f
    return md.Trajectory(np.stack(frames), top)


class TestContacts:
    def test_detects_contacts_and_fingerprint(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.contacts import Contacts

        c = Contacts(ligand_resname="LIG", cutoff=0.4, output_dir=tmp_path)
        df = c.compute(complex_in_contact)
        assert list(df.columns) == ["frame", "n_contacts"]
        assert df["n_contacts"].sum() > 0
        # Per-residue fingerprint has frequencies in (0, 1].
        assert not c._per_residue.empty
        assert (c._per_residue["contact_frequency"] <= 1.0).all()
        assert (c._per_residue["contact_frequency"] > 0.0).all()

    def test_writes_both_files(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.contacts import Contacts

        c = Contacts(ligand_resname="LIG", cutoff=0.4, output_dir=tmp_path)
        result = c.run(complex_in_contact)
        assert result.status == "ok"
        assert (result.output_dir / "pl_contacts.dat").exists()
        assert (result.output_dir / "pl_contacts_per_residue.csv").exists()

    def test_requires_ligand_resname(self, tmp_path):
        from fastmdxplora.analysis.contacts import Contacts

        with pytest.raises(ValueError, match="requires `ligand_resname`"):
            Contacts(output_dir=tmp_path)

    def test_unknown_resname_errors(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.contacts import Contacts

        c = Contacts(ligand_resname="ZZZ", output_dir=tmp_path)
        result = c.run(complex_in_contact)
        assert result.status == "error"
        assert "ZZZ" in result.message


class TestLigandRMSF:
    def test_computes_per_atom(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.ligand_rmsf import LigandRMSF

        lr = LigandRMSF(ligand_resname="LIG", output_dir=tmp_path)
        result = lr.compute(complex_in_contact)
        assert result.shape == (3, 2)  # 3 ligand atoms, (serial, rmsf)
        assert (result[:, 1] >= 0).all()

    def test_runs_end_to_end(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.ligand_rmsf import LigandRMSF

        result = LigandRMSF(ligand_resname="LIG", output_dir=tmp_path).run(
            complex_in_contact
        )
        assert result.status == "ok"
        assert (result.output_dir / "ligand_rmsf.dat").exists()

    def test_requires_ligand_resname(self, tmp_path):
        from fastmdxplora.analysis.ligand_rmsf import LigandRMSF

        with pytest.raises(ValueError, match="requires `ligand_resname`"):
            LigandRMSF(output_dir=tmp_path)


class TestProteinLigandHBonds:
    def test_computes_per_frame(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.pl_hbonds import ProteinLigandHBonds

        plh = ProteinLigandHBonds(ligand_resname="LIG", output_dir=tmp_path)
        df = plh.compute(complex_in_contact)
        assert list(df.columns) == ["frame", "n_hbonds"]
        assert len(df) == 4
        assert (df["n_hbonds"] >= 0).all()  # count is non-negative

    def test_runs_end_to_end(self, complex_in_contact, tmp_path):
        from fastmdxplora.analysis.pl_hbonds import ProteinLigandHBonds

        result = ProteinLigandHBonds(
            ligand_resname="LIG", output_dir=tmp_path
        ).run(complex_in_contact)
        assert result.status == "ok"
        assert (result.output_dir / "pl_hbonds.dat").exists()

    def test_requires_ligand_resname(self, tmp_path):
        from fastmdxplora.analysis.pl_hbonds import ProteinLigandHBonds

        with pytest.raises(ValueError, match="requires `ligand_resname`"):
            ProteinLigandHBonds(output_dir=tmp_path)


class TestAllLigandAnalysesAutoRoute:
    def test_all_four_skipped_without_ligand(self, tmp_path):
        pdb = tmp_path / "po.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        ao = AnalysisOrchestrator(str(pdb), output_dir=tmp_path / "o", ligand_resname=None)
        plan = ao._build_plan(None, None)
        for n in ("ligand_rmsd", "ligand_rmsf", "pl_contacts", "pl_hbonds"):
            assert n not in plan

    def test_all_four_run_with_ligand(self, tmp_path):
        pdb = tmp_path / "po.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        ao = AnalysisOrchestrator(str(pdb), output_dir=tmp_path / "o", ligand_resname="LIG")
        plan = ao._build_plan(None, None)
        for n in ("ligand_rmsd", "ligand_rmsf", "pl_contacts", "pl_hbonds"):
            assert n in plan
