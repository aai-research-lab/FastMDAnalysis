"""
Tests for the validation script validate_fastmda.py
"""
import pytest
import sys
import subprocess
from pathlib import Path


def test_validation_script_exists():
    """Test that the validation script exists and is executable."""
    script_path = Path(__file__).parent.parent / "validate_fastmda.py"
    assert script_path.exists(), "validate_fastmda.py should exist"
    assert script_path.is_file(), "validate_fastmda.py should be a file"


def test_validation_help():
    """Test that the validation script can display help."""
    script_path = Path(__file__).parent.parent / "validate_fastmda.py"
    result = subprocess.run(
        [sys.executable, str(script_path), "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, "Help command should succeed"
    assert "FastMDAnalysis" in result.stdout, "Help should mention FastMDAnalysis"
    assert "--frames" in result.stdout, "Help should document --frames"
    assert "--atoms" in result.stdout, "Help should document --atoms"


@pytest.mark.slow
def test_validation_runs_with_small_dataset(tmp_path):
    """Test that the validation script runs successfully with a small dataset."""
    script_path = Path(__file__).parent.parent / "validate_fastmda.py"
    
    # Run with a very small frame selection for speed
    result = subprocess.run(
        [
            sys.executable, str(script_path),
            "--frames", "0:10:5",
            "--atoms", "protein",
            "--output-dir", str(tmp_path / "validation_test")
        ],
        capture_output=True,
        text=True,
        timeout=120
    )
    
    # Check that it ran without crashing
    assert result.returncode == 0, f"Validation should succeed. stderr: {result.stderr}"
    
    # Check that output files were created
    output_dir = tmp_path / "validation_test"
    assert (output_dir / "validation_report.json").exists(), "JSON report should be created"
    assert (output_dir / "validation_summary.csv").exists(), "CSV summary should be created"
    
    # Check that the output mentions key analyses
    assert "RMSD" in result.stdout or "RMSD" in result.stderr
    assert "RMSF" in result.stdout or "RMSF" in result.stderr
    assert "SASA" in result.stdout or "SASA" in result.stderr


def test_validation_csv_format(tmp_path):
    """Test that the CSV output has the correct format."""
    script_path = Path(__file__).parent.parent / "validate_fastmda.py"
    
    # Run validation
    result = subprocess.run(
        [
            sys.executable, str(script_path),
            "--frames", "0:10:5",
            "--atoms", "protein",
            "--output-dir", str(tmp_path / "validation_csv_test")
        ],
        capture_output=True,
        text=True,
        timeout=120
    )
    
    assert result.returncode == 0, "Validation should succeed"
    
    # Read CSV and check columns
    csv_path = tmp_path / "validation_csv_test" / "validation_summary.csv"
    assert csv_path.exists(), "CSV should exist"
    
    with open(csv_path, 'r') as f:
        header = f.readline().strip()
        required_columns = [
            'analysis_name', 'backend', 'metric', 'status', 'shape_match',
            'max_abs_diff', 'mean_abs_diff', 'rmse', 'mismatch_count', 'detail',
            'fastmda_min', 'fastmda_max', 'fastmda_mean', 'fastmda_std',
            'ref_min', 'ref_max', 'ref_mean', 'ref_std',
            'fastmda_shape', 'ref_shape'
        ]
        for col in required_columns:
            assert col in header, f"CSV should contain column '{col}'"


def test_validation_json_format(tmp_path):
    """Test that the JSON output is valid and has expected structure."""
    import json
    
    script_path = Path(__file__).parent.parent / "validate_fastmda.py"
    
    # Run validation
    result = subprocess.run(
        [
            sys.executable, str(script_path),
            "--frames", "0:10:5",
            "--atoms", "protein",
            "--output-dir", str(tmp_path / "validation_json_test")
        ],
        capture_output=True,
        text=True,
        timeout=120
    )
    
    assert result.returncode == 0, "Validation should succeed"
    
    # Read and parse JSON
    json_path = tmp_path / "validation_json_test" / "validation_report.json"
    assert json_path.exists(), "JSON should exist"
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Check that it's a list
    assert isinstance(data, list), "JSON should contain a list of results"
    assert len(data) > 0, "JSON should contain at least one result"
    
    # Check that each result has expected keys
    for result_item in data:
        assert 'name' in result_item, "Each result should have a 'name'"
        assert 'status' in result_item, "Each result should have a 'status'"
