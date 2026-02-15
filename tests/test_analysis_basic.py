import numpy as np

def test_rmsd_basic(fastmda):
    a = fastmda.rmsd(ref=0)
    assert hasattr(a, "data")
    assert isinstance(a.data, np.ndarray)
    assert a.data.size > 0
    assert np.isfinite(a.data).all()

def test_rmsf_basic(fastmda):
    a = fastmda.rmsf()
    assert hasattr(a, "data")
    assert a.data.ndim in (1, 2)
    assert np.isfinite(a.data).all()

def test_rg_basic(fastmda):
    a = fastmda.rg()
    assert hasattr(a, "data")
    assert np.isfinite(a.data).all()

def test_hbonds_basic(fastmda):
    a = fastmda.hbonds()
    assert hasattr(a, "data")

def test_ss_basic(fastmda):
    a = fastmda.ss()
    assert hasattr(a, "data")

def test_sasa_basic(fastmda):
    a = fastmda.sasa(probe_radius=0.14)
    assert hasattr(a, "data")

def test_phi_basic(fastmda):
    a = fastmda.phi()
    assert hasattr(a, "data")
    assert isinstance(a.data, np.ndarray)
    assert a.data.size > 0
    assert np.isfinite(a.data).all()

def test_psi_basic(fastmda):
    a = fastmda.psi()
    assert hasattr(a, "data")
    assert isinstance(a.data, np.ndarray)
    assert a.data.size > 0
    assert np.isfinite(a.data).all()

def test_omega_basic(fastmda):
    a = fastmda.omega()
    assert hasattr(a, "data")
    assert isinstance(a.data, np.ndarray)
    assert a.data.size > 0
    assert np.isfinite(a.data).all()

def test_dihedrals_basic(fastmda):
    a = fastmda.dihedrals()
    assert hasattr(a, "data")
    assert isinstance(a.data, dict)
    assert "phi_avg" in a.data
    assert "psi_avg" in a.data
    assert "omega_avg" in a.data


def test_compute_stat_outputs(fastmda, tmp_path):
    rmsd = fastmda.rmsd(ref=0, compute_stat=True, output=str(tmp_path / "rmsd"))
    assert "rmsd_stats" in rmsd.results
    assert np.isfinite(rmsd.results["rmsd_stats"]["mean"])
    assert np.isfinite(rmsd.results["rmsd_stats"]["std"])

    rmsf = fastmda.rmsf(compute_stat=True, output=str(tmp_path / "rmsf"))
    assert "rmsf_stats" in rmsf.results
    assert np.isfinite(rmsf.results["rmsf_stats"]["mean"])
    assert np.isfinite(rmsf.results["rmsf_stats"]["std"])

    rg = fastmda.rg(compute_stat=True, output=str(tmp_path / "rg"))
    assert "rg_stats" in rg.results
    assert np.isfinite(rg.results["rg_stats"]["mean"])
    assert np.isfinite(rg.results["rg_stats"]["std"])

    sasa = fastmda.sasa(probe_radius=0.14, compute_stat=True, output=str(tmp_path / "sasa"))
    assert "total_sasa_stats" in sasa.results
    assert np.isfinite(sasa.results["total_sasa_stats"]["mean"])
    assert np.isfinite(sasa.results["total_sasa_stats"]["std"])

    qvalue = fastmda.qvalue(compute_stat=True, output=str(tmp_path / "qvalue"))
    assert "qvalue_stats" in qvalue.results
    assert np.isfinite(qvalue.results["qvalue_stats"]["mean"])
    assert np.isfinite(qvalue.results["qvalue_stats"]["std"])

