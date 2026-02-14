#!/usr/bin/env python
"""Quick validation test for Q-value analysis."""

from fastmdanalysis import FastMDAnalysis, datasets
import tempfile
import json
from pathlib import Path

# Load test data
ds = datasets.TrpCage
fastmda = FastMDAnalysis(ds.traj, ds.top, frames=(0, None, 10), atoms='protein')

# Test Q-value analysis
with tempfile.TemporaryDirectory() as tmpdir:
    result = fastmda.qvalue(reference_frame=0, output=tmpdir)

    # Check outputs
    outdir = Path(tmpdir)
    dat_file = outdir / 'qvalue.dat'
    png_file = outdir / 'qvalue.png'
    meta_file = outdir / 'qvalue_metadata.json'

    print('Testing Q-value analysis output files:')
    print(f'  Data file exists: {dat_file.exists()}')
    print(f'  Plot file exists: {png_file.exists()}')
    print(f'  Metadata file exists: {meta_file.exists()}')

    # Load and display metadata
    with open(meta_file) as f:
        metadata = json.load(f)

    print('\nMetadata:')
    nc = metadata['native_contacts_count']
    ref = metadata['reference_frame']
    beta = metadata['beta_const_nm_inv']
    lam = metadata['lambda_const']
    cutoff = metadata['native_cutoff_nm']
    nf = metadata['n_frames']
    print(f'  Native contacts: {nc}')
    print(f'  Reference frame: {ref}')
    print(f'  Beta constant: {beta} nm^-1')
    print(f'  Lambda constant: {lam}')
    print(f'  Native cutoff: {cutoff} nm')
    print(f'  Frames analyzed: {nf}')

    # Check Q-value range
    q_min = result.data.min()
    q_max = result.data.max()
    q_mean = result.data.mean()
    print('\nQ-value statistics:')
    print(f'  Min: {q_min:.4f}')
    print(f'  Max: {q_max:.4f}')
    print(f'  Mean: {q_mean:.4f}')
    print('\nAnalysis completed successfully!')
