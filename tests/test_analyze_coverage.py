"""
Tests to improve coverage of src/fastmdanalysis/analysis/analyze.py
"""
import pytest
import tempfile
import warnings
from pathlib import Path
from unittest.mock import patch, MagicMock
from fastmdanalysis import FastMDAnalysis


class MockAnalysis:
    """Mock analysis class for testing"""
    def __init__(self, outdir=None):
        self.outdir = outdir if outdir else Path("mock_output")


def test_discover_available(minimal_md_files):
    """Test _discover_available function"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    # Import the function directly
    from fastmdanalysis.analysis.analyze import _discover_available
    
    available = _discover_available(fa)
    expected = ['rmsd', 'rmsf', 'rg', 'hbonds', 'ss', 'sasa', 'dimred', 'cluster']
    assert set(available) == set(expected)


def test_validate_options():
    """Test _validate_options function"""
    from fastmdanalysis.analysis.analyze import _validate_options
    
    # Test None input
    result = _validate_options(None)
    assert result == {}
    
    # Test valid options
    options = {
        'rmsd': {'ref': 0},
        'rmsf': {'atoms': 'protein'}
    }
    result = _validate_options(options)
    assert result['rmsd']['ref'] == 0
    assert result['rmsf']['atoms'] == 'protein'
    
    # Test invalid options type
    with pytest.raises(TypeError):
        _validate_options("invalid")
    
    # Test invalid nested type
    with pytest.raises(TypeError):
        _validate_options({'rmsd': 'invalid'})


def test_final_list_all():
    """Test _final_list with 'all' include"""
    from fastmdanalysis.analysis.analyze import _final_list
    
    available = ['rmsd', 'rmsf', 'rg']
    
    # Test include=None (should get all available)
    result = _final_list(available, include=None, exclude=None)
    assert result == ['rmsd', 'rmsf', 'rg']
    
    # Test include=['all']
    result = _final_list(available, include=['all'], exclude=None)
    assert result == ['rmsd', 'rmsf', 'rg']


def test_final_list_include_exclude():
    """Test _final_list with specific include and exclude"""
    from fastmdanalysis.analysis.analyze import _final_list
    
    available = ['rmsd', 'rmsf', 'rg', 'hbonds']
    
    # Test specific include
    result = _final_list(available, include=['rmsd', 'rg'], exclude=None)
    assert result == ['rmsd', 'rg']
    
    # Test exclude
    result = _final_list(available, include=None, exclude=['rmsf'])
    assert result == ['rmsd', 'rg', 'hbonds']
    
    # Test both include and exclude
    result = _final_list(available, include=['rmsd', 'rmsf', 'rg'], exclude=['rmsf'])
    assert result == ['rmsd', 'rg']


def test_final_list_unknown_analysis_warning():
    """Test _final_list warns about unknown analyses"""
    from fastmdanalysis.analysis.analyze import _final_list
    
    available = ['rmsd', 'rmsf']
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = _final_list(available, include=['rmsd', 'unknown'], exclude=None)
        
        # Should warn about unknown analysis but still return valid ones
        assert len(w) == 1
        assert "unknown" in str(w[0].message).lower()
        assert result == ['rmsd']


def test_final_list_no_analyses_error():
    """Test _final_list raises error when no analyses to run"""
    from fastmdanalysis.analysis.analyze import _final_list
    
    available = ['rmsd', 'rmsf']
    
    with pytest.raises(ValueError, match="No analyses to run"):
        _final_list(available, include=['hbonds'], exclude=None)  # hbonds not in available
    
    with pytest.raises(ValueError, match="No analyses to run"):
        _final_list(available, include=['rmsd'], exclude=['rmsd'])  # excluded the only one


def test_filter_kwargs():
    """Test _filter_kwargs function"""
    from fastmdanalysis.analysis.analyze import _filter_kwargs
    
    def test_func(a, b=2, *, c=3):
        return a + b + c
    
    kwargs = {'a': 1, 'b': 2, 'c': 3, 'd': 4}  # 'd' should be filtered out
    
    result = _filter_kwargs(test_func, kwargs)
    assert 'a' in result
    assert 'b' in result  
    assert 'c' in result
    assert 'd' not in result  # Unknown parameter filtered out


def test_filter_kwargs_var_keyword():
    """Test _filter_kwargs with **kwargs function"""
    from fastmdanalysis.analysis.analyze import _filter_kwargs
    
    def test_func_var(a, **kwargs):
        return a, kwargs
    
    kwargs = {'a': 1, 'b': 2, 'c': 3}  # All should pass through
    
    result = _filter_kwargs(test_func_var, kwargs)
    assert result == kwargs  # All kwargs passed through due to **kwargs


def test_unique_path(tmp_path):
    """Test _unique_path function"""
    from fastmdanalysis.analysis.analyze import _unique_path
    
    base_file = tmp_path / "test.txt"
    
    # First call should return original path
    result1 = _unique_path(base_file)
    assert result1 == base_file
    
    # Create the file
    base_file.write_text("test")
    
    # Second call should return unique path
    result2 = _unique_path(base_file)
    assert result2 == tmp_path / "test_1.txt"
    
    # Create that one too
    result2.write_text("test")
    
    # Third call should return next unique path
    result3 = _unique_path(base_file)
    assert result3 == tmp_path / "test_2.txt"


def test_dedupe_paths(tmp_path):
    """Test _dedupe_paths function"""
    from fastmdanalysis.analysis.analyze import _dedupe_paths
    
    # Create some test files
    file1 = tmp_path / "file1.txt"
    file2 = tmp_path / "file2.txt"
    file1.write_text("test1")
    file2.write_text("test2")
    
    # Test with duplicates
    paths = [file1, file2, file1, file2]  # Duplicates
    result = _dedupe_paths(paths)
    assert len(result) == 2
    assert file1 in result
    assert file2 in result


def test_inject_cluster_defaults(minimal_md_files):
    """Test _inject_cluster_defaults function"""
    from fastmdanalysis.analysis.analyze import _inject_cluster_defaults
    
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    opts = {}
    plan = ['cluster']
    
    # Import and call the function
    from fastmdanalysis.analysis.analyze import _inject_cluster_defaults
    _inject_cluster_defaults(fa, opts, plan)
    
    # Should inject cluster defaults
    assert 'cluster' in opts
    assert 'methods' in opts['cluster']
    assert 'n_clusters' in opts['cluster']
    assert opts['cluster']['methods'] == ['dbscan', 'kmeans', 'hierarchical']


def test_inject_cluster_defaults_custom_methods():
    """Test _inject_cluster_defaults with custom methods"""
    from fastmdanalysis.analysis.analyze import _inject_cluster_defaults
    
    class MockSelf:
        traj = MagicMock()
        traj.n_frames = 100
    
    opts = {'cluster': {'methods': 'dbscan,kmeans'}}
    plan = ['cluster']
    
    _inject_cluster_defaults(MockSelf(), opts, plan)
    
    assert opts['cluster']['methods'] == ['dbscan', 'kmeans']
    # Should not inject n_clusters for dbscan+kmeans (only kmeans needs it)


def test_inject_cluster_defaults_not_in_plan():
    """Test _inject_cluster_defaults when cluster not in plan"""
    from fastmdanalysis.analysis.analyze import _inject_cluster_defaults
    
    opts = {'cluster': {'methods': 'all'}}
    plan = ['rmsd', 'rmsf']  # No cluster
    
    _inject_cluster_defaults(MagicMock(), opts, plan)
    
    # Should not modify opts since cluster not in plan
    assert opts['cluster']['methods'] == 'all'


def test_run_with_strict_mode(minimal_md_files):
    """Test run function with strict mode"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    # Mock the analyses to avoid long runtimes
    with patch.object(fa, 'rmsd') as mock_rmsd:
        mock_rmsd.return_value = MockAnalysis()
        
        results = fa.analyze(
            include=['rmsd'],
            strict=True,
            verbose=False
        )
        
        assert 'rmsd' in results
        assert results['rmsd'].ok


def test_run_with_stop_on_error(minimal_md_files):
    """Test run function with stop_on_error"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    # Mock rmsd to raise an error
    with patch.object(fa, 'rmsd') as mock_rmsd:
        mock_rmsd.side_effect = ValueError("Test error")
        
        # With stop_on_error=False, should continue
        results = fa.analyze(
            include=['rmsd'],
            stop_on_error=False,
            verbose=False
        )
        
        assert 'rmsd' in results
        assert not results['rmsd'].ok
        assert isinstance(results['rmsd'].error, ValueError)


def test_run_with_custom_output_dir(minimal_md_files):
    """Test run function with custom output directory"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "custom_output"
        
        with patch.object(fa, 'rmsd') as mock_rmsd:
            mock_rmsd.return_value = MockAnalysis()
            
            results = fa.analyze(
                include=['rmsd'],
                output=str(output_dir),
                verbose=False
            )
            
            assert output_dir.exists()
            assert 'rmsd' in results


def test_run_analysis_with_outdir(minimal_md_files):
    """Test run function with analysis that has outdir attribute"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))
    
    # Create a temporary directory for the mock analysis
    with tempfile.TemporaryDirectory() as tmpdir:
        analysis_outdir = Path(tmpdir) / "analysis_out"
        analysis_outdir.mkdir()
        
        mock_analysis = MockAnalysis(outdir=analysis_outdir)
        
        with patch.object(fa, 'rmsd') as mock_rmsd:
            mock_rmsd.return_value = mock_analysis
            
            results = fa.analyze(
                include=['rmsd'],
                verbose=False
            )
            
            assert 'rmsd' in results
            assert results['rmsd'].ok


def test_run_with_unknown_options_warning(minimal_md_files):
    """Test run function warns about unknown options"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))

    with patch.object(fa, 'rmsd') as mock_rmsd:
        mock_rmsd.return_value = MockAnalysis()
        
        # Mock the _filter_kwargs to simulate dropping unknown params
        with patch('fastmdanalysis.analysis.analyze._filter_kwargs') as mock_filter:
            mock_filter.return_value = {}  # Simulate all params being filtered out
            
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                
                results = fa.analyze(
                    include=['rmsd'],
                    options={'rmsd': {'unknown_param': 'value'}},
                    verbose=True
                )
                
                # Should warn about unknown parameter being dropped
                assert len(w) >= 1
                warning_found = any('unknown_param' in str(warning.message) for warning in w)
                assert warning_found


def test_analyze_alias(minimal_md_files):
    """Test that analyze is an alias for run"""
    traj_path, top_path = minimal_md_files
    fa = FastMDAnalysis(str(traj_path), str(top_path))

    with patch.object(fa, 'rmsd') as mock_rmsd:
        mock_rmsd.return_value = MockAnalysis()

        # Both should work the same - import run directly
        from fastmdanalysis.analysis.analyze import run
        
        results1 = fa.analyze(include=['rmsd'], verbose=False)
        results2 = run(fa, include=['rmsd'], verbose=False)  # FIXED: call run as function

        assert 'rmsd' in results1
        assert 'rmsd' in results2


def test_print_summary():
    """Test _print_summary function"""
    from fastmdanalysis.analysis.analyze import _print_summary, AnalysisResult
    from io import StringIO
    import sys

    results = {
        'rmsd': AnalysisResult('rmsd', True, seconds=1.5),
        'rmsf': AnalysisResult('rmsf', False, error=ValueError("test"), seconds=2.0)
    }

    output_dir = Path("/test/output")

    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()

    try:
        _print_summary(results, output_dir)
        output = captured_output.getvalue()

        assert "Summary:" in output
        assert "rmsd" in output
        assert "rmsf" in output
        assert "OK" in output
        assert "FAIL" in output
        expected_path = str(output_dir.resolve())
        assert expected_path in output
    finally:
        sys.stdout = old_stdout  



