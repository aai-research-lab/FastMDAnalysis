# tests/test_init_coverage.py
"""
Tests to improve coverage of src/fastmdanalysis/__init__.py
"""
import pytest
import tempfile
import json
import yaml
from pathlib import Path
from fastmdanalysis import FastMDAnalysis, _normalize_frames, _load_system_config


def test_normalize_frames_string_comma_separated():
    """Test frame normalization with comma-separated string"""
    result = _normalize_frames("0,100,2")
    assert result == (0, 100, 2)


def test_normalize_frames_string_colon_separated():
    """Test frame normalization with colon-separated string"""
    result = _normalize_frames("0:100:2")
    assert result == (0, 100, 2)


def test_normalize_frames_with_none():
    """Test frame normalization with None values"""
    result = _normalize_frames("None,100,2")
    assert result == (None, 100, 2)
    
    result = _normalize_frames("0,None,2") 
    assert result == (0, None, 2)


def test_normalize_frames_tuple():
    """Test frame normalization with tuple input"""
    result = _normalize_frames((0, 100, 2))
    assert result == (0, 100, 2)


def test_normalize_frames_none():
    """Test frame normalization with None input"""
    result = _normalize_frames(None)
    assert result is None


def test_normalize_frames_invalid():
    """Test frame normalization with invalid input"""
    with pytest.raises(TypeError):
        _normalize_frames("invalid")
    
    with pytest.raises(TypeError):
        _normalize_frames((0, 100))  # Wrong length


def test_load_system_config_json(tmp_path):
    """Test loading system config from JSON file"""
    config_data = {
        "trajectory": "test.dcd",
        "topology": "test.pdb", 
        "frames": "0,100,2",
        "atoms": "protein",
        "include": ["rmsd", "rmsf"]
    }
    
    config_file = tmp_path / "test_config.json"
    with open(config_file, 'w') as f:
        json.dump(config_data, f)
    
    result = _load_system_config(str(config_file))
    assert result["trajectory"] == "test.dcd"
    assert result["topology"] == "test.pdb"
    assert result["include"] == ["rmsd", "rmsf"]


def test_load_system_config_yaml(tmp_path):
    """Test loading system config from YAML file"""
    config_data = {
        "trajectory": "test.dcd",
        "topology": "test.pdb",
        "atoms": "protein and name CA"
    }
    
    config_file = tmp_path / "test_config.yaml"
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f)
    
    result = _load_system_config(str(config_file))
    assert result["trajectory"] == "test.dcd"
    assert result["topology"] == "test.pdb"


def test_load_system_config_dict():
    """Test loading system config from dict"""
    config_data = {
        "traj": "test.dcd",  # alias
        "top": "test.pdb",   # alias  
        "atoms": "protein"
    }
    
    result = _load_system_config(config_data)
    assert result["trajectory"] == "test.dcd"
    assert result["topology"] == "test.pdb"


def test_load_system_config_invalid_type():
    """Test loading system config with invalid type"""
    with pytest.raises(TypeError):
        _load_system_config(123)  # Not a path or mapping


def test_fastmda_constructor_with_system_config(tmp_path, minimal_md_files):
    """Test FastMDAnalysis constructor with system config file"""
    traj_path, top_path = minimal_md_files
    
    config_data = {
        "trajectory": str(traj_path),
        "topology": str(top_path),
        "frames": "0,10,2",
        "atoms": "protein"
    }
    
    config_file = tmp_path / "system_config.yaml"
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f)
    
    # Test with system config file
    fa = FastMDAnalysis(system=str(config_file))
    assert fa.traj.n_frames == 5  # 10 frames with stride 2: 0,2,4,6,8
    assert fa.default_atoms == "protein"


def test_fastmda_constructor_with_system_dict(minimal_md_files):
    """Test FastMDAnalysis constructor with system config dict"""
    traj_path, top_path = minimal_md_files
    
    config_dict = {
        "trajectory": str(traj_path),
        "topology": str(top_path),
        "atoms": "all"
    }
    
    fa = FastMDAnalysis(system=config_dict)
    assert fa.traj.n_frames > 0
    assert fa.default_atoms == "all"


def test_fastmda_constructor_keyword_aliases(minimal_md_files):
    """Test FastMDAnalysis constructor with keyword aliases"""
    traj_path, top_path = minimal_md_files
    
    # Test 'trajectory' alias
    fa1 = FastMDAnalysis(trajectory=str(traj_path), topology=str(top_path))
    assert fa1.traj.n_frames > 0
    
    # Test 'traj' alias  
    fa2 = FastMDAnalysis(traj=str(traj_path), top=str(top_path))
    assert fa2.traj.n_frames > 0


def test_fastmda_constructor_multiple_trajectories(minimal_md_files):
    """Test FastMDAnalysis constructor with multiple trajectories"""
    traj_path, top_path = minimal_md_files
    
    # Pass same trajectory twice to simulate multiple files
    fa = FastMDAnalysis(traj_file=[str(traj_path), str(traj_path)], top_file=str(top_path))
    # Should concatenate both trajectories
    assert fa.full_traj.n_frames > 0


def test_fastmda_missing_required_args():
    """Test FastMDAnalysis constructor with missing required arguments"""
    with pytest.raises(ValueError):
        FastMDAnalysis()  # No trajectory or topology
    
    with pytest.raises(ValueError):
        FastMDAnalysis(traj_file="test.dcd")  # Missing topology


def test_fastmda_analysis_methods_with_atoms(minimal_md_files):
    """Test analysis methods with custom atom selections"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path), atoms="protein")
    
    # Test each analysis method with custom atoms
    analysis = fa.rmsd(atoms="name CA")
    assert analysis is not None
    
    analysis = fa.rmsf(atoms="backbone") 
    assert analysis is not None
    
    analysis = fa.rg(atoms="all")
    assert analysis is not None


def test_fastmda_cluster_with_methods(minimal_md_files):
    """Test cluster analysis with different method parameters"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    # Test with specific methods
    analysis = fa.cluster(methods=["dbscan", "kmeans"], atoms="protein")
    assert analysis is not None


def test_fastmda_analyze_with_system_defaults(minimal_md_files):
    """Test analyze method using system config defaults"""
    traj_path, top_path = minimal_md_files
    
    config_dict = {
        "trajectory": str(traj_path),
        "topology": str(top_path), 
        "include": ["rmsd", "rmsf"],
        "output": "test_output"
    }
    
    fa = FastMDAnalysis(system=config_dict)
    
    # analyze should use system defaults
    result = fa.analyze()
    assert result is not None


def test_fastmda_properties(minimal_md_files):
    """Test FastMDAnalysis properties"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    # Test default properties
    assert hasattr(fa, 'figdir')
    assert hasattr(fa, 'outdir')
    assert fa.figdir == "figures"
    assert fa.outdir == "results"

