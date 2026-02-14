"""Tests for Q-value (Fraction of Native Contacts) analysis."""

import json
import numpy as np
import pytest
from pathlib import Path


def test_qvalue_basic(fastmda):
    """Test basic Q-value computation."""
    a = fastmda.qvalue()
    assert hasattr(a, "data")
    assert isinstance(a.data, np.ndarray)
    assert a.data.size > 0
    assert a.data.ndim == 2
    assert a.data.shape[1] == 1  # Should be (N_frames, 1)
    assert np.isfinite(a.data).all()
    assert np.all((a.data >= 0) & (a.data <= 1))  # Q should be in [0, 1]


def test_qvalue_reference_frame(fastmda):
    """Test Q-value with different reference frame."""
    # Test with first frame (default)
    a1 = fastmda.qvalue(reference_frame=0)
    assert a1.data is not None

    # Test with different frame (if available)
    if fastmda.traj.n_frames > 1:
        a2 = fastmda.qvalue(reference_frame=1)
        assert a2.data is not None
        # Q-values might differ with different reference frame
        # but both should be valid
        assert np.all((a2.data >= 0) & (a2.data <= 1))


def test_qvalue_metadata(fastmda, tmp_path):
    """Test that metadata file is generated correctly."""
    output_dir = tmp_path / "q_output"
    a = fastmda.qvalue(output=str(output_dir))

    # Check metadata file exists
    metadata_path = output_dir / "qvalue_metadata.json"
    assert metadata_path.exists()

    # Load and validate metadata
    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    assert "native_contacts_count" in metadata
    assert "reference_frame" in metadata
    assert "beta_const_nm_inv" in metadata
    assert "lambda_const" in metadata
    assert "native_cutoff_nm" in metadata
    assert "n_frames" in metadata
    assert "n_atoms" in metadata

    # Validate values
    assert metadata["native_contacts_count"] > 0
    assert metadata["reference_frame"] == 0
    assert metadata["beta_const_nm_inv"] == 50.0
    assert metadata["lambda_const"] == 1.8
    assert metadata["native_cutoff_nm"] == 0.45


def test_qvalue_custom_parameters(fastmda):
    """Test Q-value with custom parameters."""
    a = fastmda.qvalue(
        reference_frame=0,
        beta_const=60.0,
        lambda_const=1.6,
        native_cutoff=0.50
    )

    assert a.data is not None
    assert np.all((a.data >= 0) & (a.data <= 1))

    # Check metadata has custom values
    assert a.metadata["beta_const_nm_inv"] == 60.0
    assert a.metadata["lambda_const"] == 1.6
    assert a.metadata["native_cutoff_nm"] == 0.50


def test_qvalue_output_files(fastmda, tmp_path):
    """Test that output files are generated correctly."""
    output_dir = tmp_path / "q_output"
    a = fastmda.qvalue(output=str(output_dir))

    # Check data file
    data_path = output_dir / "qvalue.dat"
    assert data_path.exists()

    # Check plot file
    plot_path = output_dir / "qvalue.png"
    assert plot_path.exists()

    # Check metadata file
    metadata_path = output_dir / "qvalue_metadata.json"
    assert metadata_path.exists()

    # Verify data file format
    data = np.loadtxt(data_path, skiprows=1)  # Skip header
    assert data.ndim == 1 or data.shape[1] == 1
    assert len(data) == fastmda.traj.n_frames


def test_qvalue_invalid_reference_frame(fastmda):
    """Test Q-value with invalid reference frame."""
    with pytest.raises(Exception):  # Should raise AnalysisError
        fastmda.qvalue(reference_frame=fastmda.traj.n_frames + 10)


def test_qvalue_results_dict(fastmda):
    """Test that Q-value returns proper results dictionary."""
    a = fastmda.qvalue()

    assert isinstance(a.results, dict)
    assert "qvalue" in a.results
    assert np.array_equal(a.results["qvalue"], a.data)


def test_qvalue_plot_generated(fastmda, tmp_path):
    """Test that plot is generated with proper annotations."""
    output_dir = tmp_path / "q_output"
    a = fastmda.qvalue(output=str(output_dir))

    # Check that plot file was created
    plot_path = output_dir / "qvalue.png"
    assert plot_path.exists()
    assert plot_path.stat().st_size > 0  # File has content
