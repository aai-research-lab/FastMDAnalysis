# FastMDXplora Test Run

_Generated: 2026-06-27 (UTC)_  
_Tool: FastMDXplora_

## Summary

This report was generated automatically by FastMDXplora from the outputs of an end-to-end molecular dynamics study.

## Methods

### System preparation

The input system was prepared using FastMDXplora's automated setup pipeline with the following parameters:

- **ph**: `7.4`
- **keep\_heterogens**: `False`
- **keep\_water**: `False`
- **fixed\_pdb**: `None`
- **forcefield**: `charmm36`
- **force\_field**: `None`
- **water\_model**: `None`
- **ligand**: `None`
- **ligand\_forcefield**: `None`
- **ligand\_name**: `LIG`
- **ligand\_net\_charge**: `None`
- **check\_ligand\_clashes**: `True`
- **ligand\_clash\_threshold\_nm**: `0.15`
- **solvent\_padding\_nm**: `1.0`
- **box\_shape**: `cube`
- **ion\_positive**: `Na+`
- **ion\_negative**: `Cl-`
- **ion\_concentration\_M**: `0.15`
- **neutralize**: `True`
- **nonbonded\_method**: `PME`
- **nonbonded\_cutoff\_nm**: `1.0`
- **ewald\_error\_tolerance**: `0.0005`
- **use\_switching\_function**: `True`
- **switch\_distance\_nm**: `None`
- **dispersion\_correction**: `True`
- **remove\_cm\_motion**: `False`
- **constraints**: `HBonds`
- **rigid\_water**: `True`
- **hydrogen\_mass\_amu**: `None`
- **temperature\_K**: `300.0`

### Molecular dynamics simulation

Production MD was performed with the following simulation parameters:

- **preset**: `None`
- **minimize**: `True`
- **minimize\_tolerance\_kjmol\_per\_nm**: `10.0`
- **minimize\_max\_iterations**: `0`
- **nvt\_steps**: `100`
- **npt\_steps**: `0`
- **production\_steps**: `100`
- **duration\_ns**: `None`
- **nvt\_duration\_ns**: `None`
- **npt\_duration\_ns**: `None`
- **integrator**: `langevin_middle`
- **integrator\_error\_tolerance**: `0.001`
- **timestep\_fs**: `0.5`
- **temperature\_K**: `100.0`
- **friction\_per\_ps**: `5.0`
- **pressure\_bar**: `None`
- **pressure\_atm**: `None`
- **barostat\_frequency**: `25`
- **random\_seed**: `None`
- **platform**: `CPU`
- **precision**: `double`
- **device\_index**: `None`
- **trajectory\_interval\_steps**: `10`
- **state\_interval\_steps**: `1000`
- **checkpoint\_interval\_steps**: `0`
- **plumed**: `None`

## Results

Analysis was performed on a trajectory of 10 frames and 1618 residues.
Analyses performed: rmsd, rg.

### RMSD

**Parameters:**
- `selection`: `name CA`
- `ref`: `0`
- `align`: `True`

![rmsd — rmsd](analysis/rmsd/rmsd.png)

**Figure:** RMSD tracks structural deviation from the reference conformation\. A plateau suggests the system reached a relatively stable structural state, while large drift may indicate ongoing conformational change\.

### RG

**Parameters:**
- `selection`: `protein`
- `by_chain`: `False`

![rg — rg](analysis/rg/rg.png)

**Figure:** Radius of gyration tracks global molecular compactness over time\. Stable values suggest the system maintains a consistent overall size, while large shifts may indicate compaction or expansion\.


### Simulation Energy Diagnostics

FastMDXplora generated diagnostic plots from `simulation/energy.csv` to summarize simulation stability.

![Energy Trace](simulation/energy_trace.png)

**Figure:** Potential, kinetic, and total energy over time provide a quick diagnostic of simulation stability\. Large discontinuities, runaway values, or sudden jumps may indicate unstable integration settings or problematic system preparation\.

![Simulation Health](simulation/simulation_health.png)

**Figure:** Temperature, density, and volume traces summarize ensemble health\. Stable temperature and gradually stabilizing density or volume are useful indicators that equilibration and production settings behaved as expected\.

Energy summary table: [`simulation/energy_summary.csv`](simulation/energy_summary.csv)


## Discussion

_This section is intended for the user to complete. FastMDXplora provides the analytical scaffolding; scientific interpretation remains the researcher's responsibility._

## Citation

If you use FastMDXplora in your work, please cite:

> Aina, A.; Kwan, D. FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics Trajectories. J. Comput. Chem. 2026, 47, e70350. DOI: 10.1002/jcc.70350

BibTeX:

```bibtex
@article{aina2026fastmd,
  author  = {Aina, Adekunle and Kwan, Derrick},
  title   = {FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics Trajectories},
  journal = {Journal of Computational Chemistry},
  volume  = {47},
  number  = {8},
  pages   = {e70350},
  year    = {2026},
  doi     = {10.1002/jcc.70350},
}
```

## Reproducibility

- **FastMDXplora version**: `2.0.1.dev21+g0c273923c.d20260627`
- **Python**: `3.13.13`
- **Platform**: `Windows-11-10.0.26200-SP0`
- **System input**: `1L2Y`
- **Output directory**: `test_run`

Per-phase parameter manifests for phases in this workflow are preserved at `setup/setup_parameters.json`, `simulation/simulation_parameters.json`, `analysis/analysis_manifest.json`. The complete session manifest is at `manifest.json` at the project root.
