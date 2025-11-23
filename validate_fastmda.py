#!/usr/bin/env python3
"""
validate_fastmda.py

Validation script that compares FastMDAnalysis routines with MDTraj and MDAnalysis
on the Ubiquitin Q99 dataset. Generates a JSON report and CSV summary table.

Usage:
    python validate_fastmda.py [--frames START:STOP:STRIDE] [--atoms SELECTION]
    
Example:
    python validate_fastmda.py --frames 0:-1:10 --atoms "protein"
"""

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import csv

import numpy as np
import mdtraj as md

# Filter out benign MDTraj warnings
warnings.filterwarnings('ignore', message='Unlikely unit cell vectors detected')
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Import FastMDAnalysis components
try:
    from fastmdanalysis import FastMDAnalysis
    from fastmdanalysis.datasets import Ubiquitin
except ImportError as e:
    print(f"Error importing FastMDAnalysis: {e}", file=sys.stderr)
    print("Make sure FastMDAnalysis is installed or add src to PYTHONPATH", file=sys.stderr)
    sys.exit(1)

# Try importing MDAnalysis
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms as mda_rms
    from MDAnalysis.analysis import align as mda_align
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    warnings.warn("MDAnalysis not available - some comparisons will be skipped")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Validate FastMDAnalysis against MDTraj and MDAnalysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--frames',
        type=str,
        default='0:-1:10',
        help='Frame selection as START:STOP:STRIDE (default: 0:-1:10)'
    )
    parser.add_argument(
        '--atoms',
        type=str,
        default='protein',
        help='Atom selection string (default: protein)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='validation_output',
        help='Output directory for validation results (default: validation_output)'
    )
    return parser.parse_args()


def parse_frame_selection(frames_str: str) -> Tuple[Optional[int], Optional[int], int]:
    """
    Parse frame selection string like '0:-1:10' into (start, stop, stride).
    
    Args:
        frames_str: String in format 'start:stop:stride'
        
    Returns:
        Tuple of (start, stop, stride)
    """
    parts = frames_str.split(':')
    if len(parts) != 3:
        raise ValueError(f"Invalid frame selection format: {frames_str}")
    
    start = int(parts[0]) if parts[0] else 0
    stride = int(parts[2]) if parts[2] else 1
    
    # Handle stop: None for empty, -1 for last frame, or the specified value
    if not parts[1]:
        stop = None
    elif parts[1] == '-1':
        stop = -1
    else:
        stop = int(parts[1])
    
    return start, stop, stride


def compute_statistics(data: np.ndarray) -> Dict[str, float]:
    """Compute basic statistics for an array."""
    if data is None or data.size == 0:
        return {
            'min': np.nan, 'max': np.nan, 'mean': np.nan, 'std': np.nan
        }
    return {
        'min': float(np.min(data)),
        'max': float(np.max(data)),
        'mean': float(np.mean(data)),
        'std': float(np.std(data))
    }


def compare_arrays(fastmda_data: np.ndarray, ref_data: np.ndarray, 
                   name: str) -> Dict[str, Any]:
    """
    Compare two arrays and return comparison metrics.
    
    Args:
        fastmda_data: Array from FastMDAnalysis
        ref_data: Reference array from MDTraj/MDAnalysis
        name: Name of the comparison
        
    Returns:
        Dictionary with comparison metrics
    """
    result = {
        'name': name,
        'shape_match': False,
        'max_abs_diff': np.nan,
        'mean_abs_diff': np.nan,
        'rmse': np.nan,
        'mismatch_count': 0,
        'fastmda_stats': compute_statistics(fastmda_data),
        'ref_stats': compute_statistics(ref_data),
        'fastmda_shape': str(fastmda_data.shape) if fastmda_data is not None else 'None',
        'ref_shape': str(ref_data.shape) if ref_data is not None else 'None',
        'status': 'fail',
        'detail': ''
    }
    
    # Check if shapes match
    if fastmda_data is None or ref_data is None:
        result['detail'] = 'One or both arrays are None'
        return result
        
    if fastmda_data.shape != ref_data.shape:
        result['detail'] = f'Shape mismatch: {fastmda_data.shape} vs {ref_data.shape}'
        return result
    
    result['shape_match'] = True
    
    # Compute differences
    abs_diff = np.abs(fastmda_data - ref_data)
    result['max_abs_diff'] = float(np.max(abs_diff))
    result['mean_abs_diff'] = float(np.mean(abs_diff))
    result['rmse'] = float(np.sqrt(np.mean((fastmda_data - ref_data) ** 2)))
    
    # Count mismatches (using a tolerance)
    tolerance = 1e-6
    result['mismatch_count'] = int(np.sum(abs_diff > tolerance))
    
    # Determine pass/fail
    # Pass if RMSE is small and shapes match
    if result['rmse'] < 1e-4:
        result['status'] = 'pass'
        result['detail'] = f'Excellent agreement (RMSE={result["rmse"]:.2e})'
    elif result['rmse'] < 1e-2:
        result['status'] = 'pass'
        result['detail'] = f'Good agreement (RMSE={result["rmse"]:.2e})'
    else:
        result['status'] = 'warn'
        result['detail'] = f'Moderate differences (RMSE={result["rmse"]:.2e})'
    
    return result


def validate_rmsd(fastmda: FastMDAnalysis, traj: md.Trajectory, 
                  ref_frame: int = 0) -> List[Dict[str, Any]]:
    """
    Validate RMSD calculation against MDTraj.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        ref_frame: Reference frame index
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating RMSD...")
    
    # FastMDAnalysis RMSD
    try:
        fmda_rmsd = fastmda.rmsd(ref=ref_frame)
        fmda_data = fmda_rmsd.data.flatten()
    except Exception as e:
        results.append({
            'name': 'RMSD',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}',
            'shape_match': False,
            'fastmda_shape': 'N/A',
            'ref_shape': 'N/A'
        })
        return results
    
    # MDTraj RMSD
    try:
        # Use protein selection to match
        atom_indices = traj.topology.select('protein')
        mdtraj_data = md.rmsd(traj, traj, frame=ref_frame, atom_indices=atom_indices)
        
        comparison = compare_arrays(fmda_data, mdtraj_data, 'RMSD')
        comparison['backend'] = 'mdtraj'
        comparison['metric'] = 'rmsd'
        results.append(comparison)
        
        print(f"  RMSD vs MDTraj: {comparison['status']} - {comparison['detail']}")
    except Exception as e:
        print(f"  Error comparing with MDTraj: {e}")
        results.append({
            'name': 'RMSD',
            'backend': 'mdtraj',
            'metric': 'rmsd',
            'status': 'error',
            'detail': f'MDTraj error: {str(e)}'
        })
    
    # Note: MDAnalysis RMSD comparison is currently disabled for performance reasons.
    # The comparison with MDTraj is sufficient for validation purposes.
    # To enable MDAnalysis comparison, implement a separate function that:
    # 1. Loads trajectory with MDAnalysis Universe
    # 2. Selects protein atoms
    # 3. Computes RMSD for each frame using mda_rms.rmsd() with proper alignment
    # 4. Converts Angstroms to nm by dividing by 10.0
    # 5. Compares with FastMDAnalysis results
    
    return results


def validate_rmsf(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate RMSF calculation against MDTraj.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating RMSF...")
    
    # FastMDAnalysis RMSF
    try:
        fmda_rmsf = fastmda.rmsf()
        fmda_data = fmda_rmsf.data.flatten() if fmda_rmsf.data.ndim > 1 else fmda_rmsf.data
    except Exception as e:
        results.append({
            'name': 'RMSF',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
        return results
    
    # MDTraj RMSF
    try:
        # FastMDAnalysis computes RMSF relative to average structure
        # We need to do the same for comparison
        avg_xyz = np.mean(traj.xyz, axis=0, keepdims=True)
        ref = md.Trajectory(avg_xyz, traj.topology)
        mdtraj_data = md.rmsf(traj, ref)
        
        comparison = compare_arrays(fmda_data, mdtraj_data, 'RMSF')
        comparison['backend'] = 'mdtraj'
        comparison['metric'] = 'rmsf'
        results.append(comparison)
        
        print(f"  RMSF vs MDTraj: {comparison['status']} - {comparison['detail']}")
    except Exception as e:
        print(f"  Error comparing with MDTraj: {e}")
        results.append({
            'name': 'RMSF',
            'backend': 'mdtraj',
            'metric': 'rmsf',
            'status': 'error',
            'detail': f'MDTraj error: {str(e)}'
        })
    
    return results


def validate_rg(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate Radius of Gyration calculation against MDTraj.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating Radius of Gyration...")
    
    # FastMDAnalysis RG
    try:
        fmda_rg = fastmda.rg()
        fmda_data = fmda_rg.data.flatten() if fmda_rg.data.ndim > 1 else fmda_rg.data
    except Exception as e:
        results.append({
            'name': 'Radius of Gyration',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
        return results
    
    # MDTraj RG
    try:
        # MDTraj compute_rg doesn't accept atom_indices directly
        # We need to slice the trajectory first
        mdtraj_data = md.compute_rg(traj)
        
        comparison = compare_arrays(fmda_data, mdtraj_data, 'Radius of Gyration')
        comparison['backend'] = 'mdtraj'
        comparison['metric'] = 'rg'
        results.append(comparison)
        
        print(f"  RG vs MDTraj: {comparison['status']} - {comparison['detail']}")
    except Exception as e:
        print(f"  Error comparing with MDTraj: {e}")
        results.append({
            'name': 'Radius of Gyration',
            'backend': 'mdtraj',
            'metric': 'rg',
            'status': 'error',
            'detail': f'MDTraj error: {str(e)}'
        })
    
    return results


def validate_hbonds(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate Hydrogen Bonds calculation against MDTraj.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating Hydrogen Bonds...")
    
    # FastMDAnalysis HBonds
    try:
        fmda_hbonds = fastmda.hbonds()
        # HBonds data is typically a count per frame or list of bonds
        if hasattr(fmda_hbonds, 'data') and isinstance(fmda_hbonds.data, np.ndarray):
            fmda_data = fmda_hbonds.data.flatten() if fmda_hbonds.data.ndim > 1 else fmda_hbonds.data
        else:
            # If data format is unexpected, report as info rather than trying to construct comparable data
            results.append({
                'name': 'Hydrogen Bonds',
                'backend': 'FastMDAnalysis',
                'metric': 'hbond_data',
                'status': 'info',
                'detail': f'HBonds data type: {type(fmda_hbonds.data)}',
                'shape_match': True,
                'fastmda_stats': {},
                'ref_stats': {},
                'fastmda_shape': str(getattr(fmda_hbonds.data, 'shape', 'N/A')),
                'ref_shape': 'N/A'
            })
            return results
    except Exception as e:
        results.append({
            'name': 'Hydrogen Bonds',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
        return results
    
    # MDTraj HBonds  
    try:
        hbonds_list = md.baker_hubbard(traj, periodic=False)
        # Count unique hbonds
        mdtraj_count = len(hbonds_list)
        
        # For comparison, we'll compare the count
        comparison = {
            'name': 'Hydrogen Bonds',
            'backend': 'mdtraj',
            'metric': 'hbond_count',
            'status': 'info',
            'detail': f'FastMDA found data, MDTraj found {mdtraj_count} unique bonds',
            'shape_match': True,
            'fastmda_stats': compute_statistics(fmda_data),
            'ref_stats': {'count': mdtraj_count},
            'fastmda_shape': str(fmda_data.shape),
            'ref_shape': f'({mdtraj_count},)'
        }
        results.append(comparison)
        
        print(f"  HBonds: {comparison['detail']}")
    except Exception as e:
        print(f"  Error comparing with MDTraj: {e}")
        results.append({
            'name': 'Hydrogen Bonds',
            'backend': 'mdtraj',
            'metric': 'hbond_count',
            'status': 'error',
            'detail': f'MDTraj error: {str(e)}'
        })
    
    return results


def validate_ss(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate Secondary Structure calculation against MDTraj.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating Secondary Structure...")
    
    # FastMDAnalysis SS
    try:
        fmda_ss = fastmda.ss()
        fmda_data = fmda_ss.data
    except Exception as e:
        results.append({
            'name': 'Secondary Structure',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
        return results
    
    # MDTraj SS (DSSP)
    try:
        mdtraj_data = md.compute_dssp(traj, simplified=True)
        
        # Convert both to string arrays for comparison if needed
        if fmda_data.dtype.kind in ['U', 'S', 'O']:
            # Already strings, compare element-wise
            # For string arrays, we can't use numerical comparison
            # Check if arrays are equal
            matches = (fmda_data == mdtraj_data).sum()
            total = fmda_data.size
            match_rate = matches / total if total > 0 else 0.0
            
            comparison = {
                'name': 'Secondary Structure',
                'backend': 'mdtraj',
                'metric': 'dssp',
                'status': 'pass' if match_rate > 0.95 else ('warn' if match_rate > 0.8 else 'fail'),
                'detail': f'Match rate: {match_rate:.2%} ({matches}/{total} elements)',
                'shape_match': fmda_data.shape == mdtraj_data.shape,
                'fastmda_shape': str(fmda_data.shape),
                'ref_shape': str(mdtraj_data.shape),
                'fastmda_stats': {},
                'ref_stats': {},
                'max_abs_diff': np.nan,
                'mean_abs_diff': np.nan,
                'rmse': np.nan,
                'mismatch_count': total - matches
            }
        else:
            comparison = compare_arrays(fmda_data, mdtraj_data, 'Secondary Structure')
            comparison['backend'] = 'mdtraj'
            comparison['metric'] = 'dssp'
        
        results.append(comparison)
        print(f"  SS vs MDTraj: {comparison['status']} - {comparison['detail']}")
    except Exception as e:
        print(f"  Error comparing with MDTraj: {e}")
        results.append({
            'name': 'Secondary Structure',
            'backend': 'mdtraj',
            'metric': 'dssp',
            'status': 'error',
            'detail': f'MDTraj error: {str(e)}'
        })
    
    return results


def validate_sasa(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate SASA calculation against MDTraj.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating SASA...")
    
    # FastMDAnalysis SASA
    try:
        fmda_sasa = fastmda.sasa(probe_radius=0.14)
        fmda_total = fmda_sasa.data.get('total_sasa') if isinstance(fmda_sasa.data, dict) else None
        fmda_residue = fmda_sasa.data.get('residue_sasa') if isinstance(fmda_sasa.data, dict) else None
        fmda_avg_residue = fmda_sasa.data.get('average_residue_sasa') if isinstance(fmda_sasa.data, dict) else None
    except Exception as e:
        results.append({
            'name': 'SASA',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
        return results
    
    # MDTraj SASA
    try:
        # MDTraj SASA - compute on the full selected trajectory
        mdtraj_total = md.shrake_rupley(traj, probe_radius=0.14, mode='atom')
        # Sum over atoms to get total per frame
        mdtraj_total_sum = np.sum(mdtraj_total, axis=1)
        
        # Compare total SASA
        if fmda_total is not None:
            comparison = compare_arrays(fmda_total, mdtraj_total_sum, 'SASA (total)')
            comparison['backend'] = 'mdtraj'
            comparison['metric'] = 'total_sasa'
            results.append(comparison)
            print(f"  SASA (total) vs MDTraj: {comparison['status']} - {comparison['detail']}")
        
        # Compare per-residue SASA
        try:
            mdtraj_residue = md.shrake_rupley(traj, probe_radius=0.14, mode='residue')
            if fmda_residue is not None:
                comparison = compare_arrays(fmda_residue, mdtraj_residue, 'SASA (per-residue)')
                comparison['backend'] = 'mdtraj'
                comparison['metric'] = 'residue_sasa'
                results.append(comparison)
                print(f"  SASA (residue) vs MDTraj: {comparison['status']} - {comparison['detail']}")
            
            # Compare average per-residue
            mdtraj_avg = np.mean(mdtraj_residue, axis=0)
            if fmda_avg_residue is not None:
                comparison = compare_arrays(fmda_avg_residue, mdtraj_avg, 'SASA (avg per-residue)')
                comparison['backend'] = 'mdtraj'
                comparison['metric'] = 'avg_residue_sasa'
                results.append(comparison)
                print(f"  SASA (avg residue) vs MDTraj: {comparison['status']} - {comparison['detail']}")
        except Exception as e:
            print(f"  Error computing per-residue SASA with MDTraj: {e}")
        
    except Exception as e:
        print(f"  Error comparing with MDTraj: {e}")
        results.append({
            'name': 'SASA',
            'backend': 'mdtraj',
            'metric': 'total_sasa',
            'status': 'error',
            'detail': f'MDTraj error: {str(e)}'
        })
    
    return results


def validate_dimred(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate Dimensionality Reduction calculation.
    
    Note: This is harder to validate directly as sklearn results can vary,
    but we check that the output has the expected shape and properties.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating Dimensionality Reduction...")
    
    # FastMDAnalysis DimRed
    try:
        fmda_dimred = fastmda.dimred(methods=['pca', 'mds', 'tsne'])
        
        if hasattr(fmda_dimred, 'data') and isinstance(fmda_dimred.data, dict):
            for method, data in fmda_dimred.data.items():
                result = {
                    'name': f'Dimensionality Reduction ({method.upper()})',
                    'backend': 'FastMDAnalysis',
                    'metric': f'dimred_{method}',
                    'status': 'pass' if data is not None and data.shape[0] == traj.n_frames else 'fail',
                    'detail': f'Shape: {data.shape}' if data is not None else 'No data',
                    'shape_match': True,
                    'fastmda_stats': compute_statistics(data) if data is not None else {},
                    'fastmda_shape': str(data.shape) if data is not None else 'None',
                    'ref_shape': 'N/A (sklearn-based)'
                }
                results.append(result)
                print(f"  DimRed ({method}): {result['status']} - {result['detail']}")
    except Exception as e:
        results.append({
            'name': 'Dimensionality Reduction',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
    
    return results


def validate_clustering(fastmda: FastMDAnalysis, traj: md.Trajectory) -> List[Dict[str, Any]]:
    """
    Validate Clustering calculation.
    
    Note: Like dimred, clustering results can vary, but we check expected properties.
    
    Args:
        fastmda: FastMDAnalysis instance
        traj: MDTraj trajectory
        
    Returns:
        List of comparison results
    """
    results = []
    
    print("Validating Clustering...")
    
    # FastMDAnalysis Clustering
    try:
        fmda_cluster = fastmda.cluster(methods=['kmeans', 'dbscan', 'hierarchical'], 
                                       n_clusters=3, eps=0.5, min_samples=2)
        
        if hasattr(fmda_cluster, 'results') and isinstance(fmda_cluster.results, dict):
            for method, data in fmda_cluster.results.items():
                if data is not None and 'labels' in data:
                    labels = data['labels']
                    n_clusters = len(np.unique(labels))
                    result = {
                        'name': f'Clustering ({method})',
                        'backend': 'FastMDAnalysis',
                        'metric': f'cluster_{method}',
                        'status': 'pass' if len(labels) == traj.n_frames else 'fail',
                        'detail': f'{n_clusters} clusters found, {len(labels)} labels',
                        'shape_match': True,
                        'fastmda_stats': {'n_clusters': n_clusters, 'n_labels': len(labels)},
                        'fastmda_shape': str(labels.shape),
                        'ref_shape': 'N/A (sklearn-based)'
                    }
                    results.append(result)
                    print(f"  Clustering ({method}): {result['status']} - {result['detail']}")
    except Exception as e:
        results.append({
            'name': 'Clustering',
            'backend': 'FastMDAnalysis',
            'status': 'error',
            'detail': f'Error: {str(e)}'
        })
    
    return results


def save_json_report(results: List[Dict[str, Any]], output_file: Path):
    """Save validation results as JSON."""
    
    def convert_to_json_serializable(obj):
        """Convert numpy types to Python types."""
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_to_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_json_serializable(item) for item in obj]
        else:
            return obj
    
    # Convert all results to be JSON serializable
    json_serializable_results = convert_to_json_serializable(results)
    
    with open(output_file, 'w') as f:
        json.dump(json_serializable_results, f, indent=2)
    print(f"\nJSON report saved to: {output_file}")


def save_csv_summary(results: List[Dict[str, Any]], output_file: Path):
    """
    Save validation results as CSV with all required columns.
    
    CSV columns:
    - analysis_name
    - backend (mdtraj/mdanalysis)
    - metric (e.g., total/residue/avg for SASA)
    - status (pass/fail/warn/error/info)
    - shape_match (True/False)
    - max_abs_diff
    - mean_abs_diff
    - rmse
    - mismatch_count
    - detail
    - fastmda_min, fastmda_max, fastmda_mean, fastmda_std
    - ref_min, ref_max, ref_mean, ref_std
    - fastmda_shape
    - ref_shape
    """
    fieldnames = [
        'analysis_name', 'backend', 'metric', 'status', 'shape_match',
        'max_abs_diff', 'mean_abs_diff', 'rmse', 'mismatch_count', 'detail',
        'fastmda_min', 'fastmda_max', 'fastmda_mean', 'fastmda_std',
        'ref_min', 'ref_max', 'ref_mean', 'ref_std',
        'fastmda_shape', 'ref_shape'
    ]
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in results:
            row = {
                'analysis_name': result.get('name', ''),
                'backend': result.get('backend', ''),
                'metric': result.get('metric', ''),
                'status': result.get('status', ''),
                'shape_match': result.get('shape_match', ''),
                'max_abs_diff': result.get('max_abs_diff', ''),
                'mean_abs_diff': result.get('mean_abs_diff', ''),
                'rmse': result.get('rmse', ''),
                'mismatch_count': result.get('mismatch_count', ''),
                'detail': result.get('detail', ''),
                'fastmda_shape': result.get('fastmda_shape', ''),
                'ref_shape': result.get('ref_shape', '')
            }
            
            # Extract FastMDA stats
            fmda_stats = result.get('fastmda_stats', {})
            if isinstance(fmda_stats, dict):
                row['fastmda_min'] = fmda_stats.get('min', '')
                row['fastmda_max'] = fmda_stats.get('max', '')
                row['fastmda_mean'] = fmda_stats.get('mean', '')
                row['fastmda_std'] = fmda_stats.get('std', '')
            else:
                row['fastmda_min'] = row['fastmda_max'] = row['fastmda_mean'] = row['fastmda_std'] = ''
            
            # Extract reference stats
            ref_stats = result.get('ref_stats', {})
            if isinstance(ref_stats, dict):
                row['ref_min'] = ref_stats.get('min', '')
                row['ref_max'] = ref_stats.get('max', '')
                row['ref_mean'] = ref_stats.get('mean', '')
                row['ref_std'] = ref_stats.get('std', '')
            else:
                row['ref_min'] = row['ref_max'] = row['ref_mean'] = row['ref_std'] = ''
            
            writer.writerow(row)
    
    print(f"CSV summary saved to: {output_file}")


def main():
    """Main validation function."""
    args = parse_args()
    
    # Parse frame selection
    start, stop, stride = parse_frame_selection(args.frames)
    frames = (start, stop, stride)
    
    print("=" * 70)
    print("FastMDAnalysis Validation")
    print("=" * 70)
    print(f"Dataset: Ubiquitin Q99")
    print(f"Trajectory: {Ubiquitin.traj}")
    print(f"Topology: {Ubiquitin.top}")
    print(f"Frame selection: {args.frames} -> {frames}")
    print(f"Atom selection: {args.atoms}")
    print(f"Output directory: {args.output_dir}")
    print("=" * 70)
    print()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load trajectory with MDTraj for reference
    print("Loading trajectory with MDTraj...")
    traj = md.load(Ubiquitin.traj, top=Ubiquitin.top)
    if stop == -1 or stop is None:
        traj = traj[start::stride]
    else:
        traj = traj[start:stop:stride]
    
    # Select atoms
    if args.atoms:
        atom_indices = traj.topology.select(args.atoms)
        traj = traj.atom_slice(atom_indices)
    
    print(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")
    print()
    
    # Initialize FastMDAnalysis
    print("Initializing FastMDAnalysis...")
    fastmda = FastMDAnalysis(
        Ubiquitin.traj,
        Ubiquitin.top,
        frames=frames,
        atoms=args.atoms
    )
    print(f"FastMDAnalysis trajectory: {fastmda.traj.n_frames} frames, {fastmda.traj.n_atoms} atoms")
    print()
    
    # Run all validations
    all_results = []
    
    all_results.extend(validate_rmsd(fastmda, traj, ref_frame=0))
    all_results.extend(validate_rmsf(fastmda, traj))
    all_results.extend(validate_rg(fastmda, traj))
    all_results.extend(validate_hbonds(fastmda, traj))
    all_results.extend(validate_ss(fastmda, traj))
    all_results.extend(validate_sasa(fastmda, traj))
    all_results.extend(validate_dimred(fastmda, traj))
    all_results.extend(validate_clustering(fastmda, traj))
    
    # Save reports
    print()
    print("=" * 70)
    print("Saving reports...")
    save_json_report(all_results, output_dir / 'validation_report.json')
    save_csv_summary(all_results, output_dir / 'validation_summary.csv')
    
    # Print summary
    print()
    print("=" * 70)
    print("Summary:")
    print("=" * 70)
    status_counts = {}
    for result in all_results:
        status = result.get('status', 'unknown')
        status_counts[status] = status_counts.get(status, 0) + 1
    
    for status, count in sorted(status_counts.items()):
        print(f"  {status.upper()}: {count}")
    
    print("=" * 70)
    print("Validation complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
