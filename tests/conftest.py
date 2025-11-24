# tests/conftest.py
# Enable Matplotlib's built-in pytest fixtures (e.g., `matplotlib`)
pytest_plugins = ("matplotlib.testing.conftest",)

import warnings
from pathlib import Path
import pytest
import mdtraj as md

from fastmdanalysis.utils import load_trajectory
from fastmdanalysis import FastMDAnalysis

# Silence the benign MDTraj CRYST1 warning in this dataset
warnings.filterwarnings(
    "ignore",
    message="Unlikely unit cell vectors detected in PDB file likely resulting from a dummy CRYST1 record",
    category=UserWarning,
    module="mdtraj",
)

@pytest.fixture(scope="session")
def dataset_paths():
    """Return (traj_path, top_path) for the small TrpCage dataset."""
    from fastmdanalysis import datasets
    if hasattr(datasets, "TrpCage"):
        ds = datasets.TrpCage
        return ds.traj, ds.top
    # fallback for older naming
    ds = getattr(datasets, "trp_cage")
    return ds.traj, ds.top

@pytest.fixture(scope="session")
def traj(dataset_paths):
    traj_path, top_path = dataset_paths
    # Subsample frames to keep CI fast
    t = load_trajectory(traj_path, top_path, frames=(0, None, 10))
    t.topology.create_standard_bonds()
    return t

@pytest.fixture(scope="function")
def outdir(tmp_path):
    return tmp_path

@pytest.fixture(scope="function")
def fastmda(dataset_paths):
    traj_path, top_path = dataset_paths
    # Same frame stride for consistency with 'traj' fixture
    return FastMDAnalysis(traj_path, top_path, frames=(0, None, 10), atoms="protein")

# --- Compatibility fixture: behave like pytest-matplotlib's `matplotlib` ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

@pytest.fixture
def matplotlib(mpl_test_settings):
    """Provide a minimal 'matplotlib' fixture without requiring pytest-matplotlib."""
    try:
        yield matplotlib
    finally:
        plt.close("all")


@pytest.fixture
def minimal_md_files():
    """Provide minimal MD trajectory and topology files for testing"""
    data_dir = Path(__file__).parent.parent / "src" / "fastmdanalysis" / "data"
    traj_file = data_dir / "trp_cage.dcd"
    top_file = data_dir / "trp_cage.pdb"
    return traj_file, top_file
