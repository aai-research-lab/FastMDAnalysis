# FastMDXplora Study: tri\-ala

_Generated: 2026-10-18 (UTC)_
_Software: FastMDXplora_
_Dashboard: [dashboard.html](dashboard.html)_

## Summary

The prepared system contained 1,409 atoms and was simulated for 4 ps at 300 K. Of 6 observables assessed, 2 had equilibrated, 1 had not, 3 could not be judged from a run this length; see [Convergence](#convergence) before using any average from this run.

## Methods

**System preparation.** Starting coordinates were taken from `/tmp/D/real/tri-ala.pdb`. The file itself is recorded in the setup manifest by its SHA-256 digest (e1c61caf90ab), which identifies it independently of the path it was read from. Missing heavy atoms and hydrogens were added with PDBFixer at pH 7.0, which sets the protonation of ionizable side chains. Heterogens: non-standard residues were decided per component: any ligand that could be parameterised was kept and prepared, crystallographic water and buffer were removed, and setup stopped rather than guess where the structure did not determine what to simulate. The protein was described by the amber14-all.xml, amber14/tip3p.xml force field with the tip3p water model. The complex was solvated in a dodecahedron box sized to leave at least 1.40 nm between the solute and its nearest periodic image (grown from the 1.20 nm asked for, so that the box is at least twice the cutoff across), with Na+/Cl- ions at 0.15 M. The solvated system contained 1,409 atoms. Hydrogens and ions were placed with random seed 859956268 (`setup.random_seed`), which reproduces the solvated system. Long-range electrostatics were treated with PME and a real-space cutoff of 0.90 nm. Bonds were constrained with `HBonds`, and water was held rigid.

**Simulation protocol.** The system was energy-minimized before dynamics. Equilibration comprised 0.4 ps in the NVT ensemble followed by 0.4 ps in the NPT ensemble. Production dynamics were run for 4 ps in the NPT ensemble. Equations of motion were integrated with the Langevin middle-scheme integrator and a 2 fs timestep, at 300.0 K with a friction coefficient of 1.0 ps⁻¹. Pressure was maintained at 1.0 bar, with the barostat applied every 25 steps. Coordinates were written every 0.04 ps. Simulations used OpenMM on the CPU platform. No random seed was fixed, so velocities differ between runs.

**Analysis.** Each per-frame quantity (`end_to_end`, `hbonds`, `rg`, `rmsd` and `sasa`) was averaged over the frames after an equilibration period detected automatically, as the latest start that keeps at least 90% of the most statistically independent samples any start keeps (after Chodera, J. Chem. Theory Comput. 2016, 12, 1799). Standard errors account for the correlation between frames through the statistical inefficiency, the number of frames per independent sample, estimated from the autocorrelation function summed over Geyer's initial positive sequence (Stat. Sci. 1992, 7, 473) and corrected for the bias of an autocorrelation taken about the sample mean. Where frames were discarded, the error is that of the whole run scaled to the frames kept, unless the discard at least doubled the independent samples. An error is given only where the frames averaged hold at least 10 independent samples and span at least 25 times their own statistical inefficiency. Frames whose statistical inefficiency, before the correction for the sample mean, is under 2 are exempt from the second condition; their error is still computed with the corrected inefficiency. The equilibration periods discarded ran from 60 to 66 of the 100 frames analysed. No error is given for `end_to_end`, `hbonds` and `rg`; the report says why for each.

**Software.** System setup, simulation and analysis were performed with FastMDXplora 0.1.0.

**Tools.** FastMDXplora uses matplotlib 3.11.2, MDTraj 1.11.1.post2, NumPy 2.4.6, OpenMM 8.6.1, pandas 3.0.6, PDBFixer 1.12.0, scikit-learn 1.9.1, SciPy 1.17.1.

### Every setting used

The paragraphs above say what a methods section says; this is the complete record, for anyone repeating the run exactly.

### System preparation

The input system was prepared using FastMDXplora's automated setup pipeline with the following parameters:

- **agent**: `None`
- **ph**: `7.0`
- **heterogens**: `auto`
- **mutations**: `None`
- **mutation\_chain**: `None`
- **protonation\_margin**: `1.0`
- **residue\_states**: `None`
- **model**: `None`
- **replace\_nonstandard\_residues**: `True`
- **chains**: `None`
- **build\_missing\_termini**: `False`
- **keep\_heterogens**: `False`
- **keep\_water**: `False`
- **fixed\_pdb**: `None`
- **preparation\_audit**: `False`
- **forcefield**: `auto`
- **force\_field**: `None`
- **water\_model**: `None`
- **ligand**: `None`
- **ligand\_forcefield**: `None`
- **ligand\_name**: `LIG`
- **ligand\_resname**: `None`
- **ligand\_net\_charge**: `None`
- **ligand\_pose**: `auto`
- **check\_ligand\_clashes**: `True`
- **ligand\_clash\_threshold\_nm**: `0.15`
- **membrane**: `None`
- **membrane\_orient**: `False`
- **membrane\_orientation\_checked**: `False`
- **membrane\_center\_z\_nm**: `None`
- **solvent\_padding\_nm**: `1.2`
- **box\_shape**: `dodecahedron`
- **ion\_positive**: `Na+`
- **ion\_negative**: `Cl-`
- **ion\_concentration\_M**: `0.15`
- **neutralize**: `True`
- **random\_seed**: `859956268`, drawn for this study
- **nonbonded\_method**: `PME`
- **nonbonded\_cutoff\_nm**: `0.9`
- **ewald\_error\_tolerance**: `0.0005`
- **use\_switching\_function**: `None`
- **switch\_distance\_nm**: `None`
- **dispersion\_correction**: `True`
- **remove\_cm\_motion**: `True`
- **constraints**: `HBonds`
- **rigid\_water**: `True`
- **hydrogen\_mass\_amu**: `None`
- **temperature\_K**: `300.0`

### Molecular dynamics simulation

Production MD was performed with the following simulation parameters:

- **agent**: `None`
- **duration\_ns**: `0.004`
- **nvt\_duration\_ns**: `0.0004`
- **npt\_duration\_ns**: `0.0004`
- **nvt\_steps**: `200`
- **ensemble**: `npt`
- **npt\_steps**: `200`
- **production\_steps**: `2000`
- **setup\_from**: `None`
- **prepared\_from**: `None`
- **extra\_ns**: `None`
- **stop\_when**: `None`
- **resume\_unsealed**: `False`
- **resume\_from**: `None`
- **minimize**: `True`
- **integrator**: `langevin_middle`
- **integrator\_error\_tolerance**: `0.001`
- **minimize\_tolerance\_kjmol\_per\_nm**: `10.0`
- **minimize\_max\_iterations**: `0`
- **timestep\_fs**: `2`
- **temperature\_K**: `300.0`
- **pressure\_bar**: `1.0`
- **pressure\_atm**: `0.98692`
- **friction\_per\_ps**: `1.0`
- **barostat\_frequency**: `25`
- **random\_seed**: `None`
- **platform**: `CPU`
- **precision**: `mixed`
- **device\_index**: `None`
- **trajectory\_interval\_steps**: `20`
- **state\_interval\_steps**: `1000`
- **save\_selection**: `not water`
- **checkpoint\_interval\_steps**: `10000`
- **live\_telemetry**: `True`
- **telemetry\_interval**: `1000`
- **dashboard\_ligand\_resname**: `None`
- **dashboard\_binding\_pocket\_cutoff\_A**: `5.0`
- **dashboard\_max\_playback\_frames**: `2000`
- **restrain**: `None`
- **restraint\_release**: `None`
- **restrain\_production**: `False`
- **umbrella**: `None`
- **steered**: `None`
- **metadynamics**: `None`
- **plumed**: `None`

## Results

Analysis was performed on a trajectory of 100 frames: 3 protein residues and 2 other residues (ions and any cofactors kept by the save selection).

Analyses performed: rmsd, rmsf, rg, hbonds, ss, sasa, dihedrals, order\_parameters, thermodynamics, end\_to\_end, moments\_of\_inertia, cluster, dimred.

### Analysis Summary Figure

![Analysis summary](analysis_summary.png)

_Panel inclusion and skipped optional source figures are recorded in `analysis_summary_manifest.json`._

### RMSD

Root-Mean-Square Deviation (RMSD).

**Parameters:**

- `selection`: `protein and name CA`
- `ref`: `0`
- `align`: `True`

**Findings:**

- Mean over the equilibrated part of the run: 0.0078 ± 0.0010 nm, from 16 independent samples after discarding 66 frames.

![rmsd](../analysis/rmsd/rmsd.png)

### RMSF

Root-Mean-Square Fluctuation (RMSF).

**Parameters:**

- `selection`: `protein and name CA`
- `ref`: `0`
- `per_residue`: `True`
- `start`: `equilibrated`
- `per_residue_average`: `sqrt(sum_i m_i MSF_i / sum_i m_i) over the residue's selected atoms, as gmx rmsf -res`

![rmsf](../analysis/rmsf/rmsf.png)

### RG

Radius of Gyration (Rg).

**Parameters:**

- `selection`: `protein`
- `by_chain`: `False`
- `mass_weighted`: `True`

**Findings:**

- Mean over the equilibrated part of the run: 0.3355 nm, from fewer than one independent sample after discarding 66 frames.
- This run is not long against its own correlation time: 34 frames hold 0.5 independent samples by their own estimate, and a correlated run shorter than 25 of them reads its correlation time low (only frames whose statistical inefficiency, before the correction for the sample mean, is under 2 are exempt). That count is an upper bound, so an error computed from it would be a lower bound -- and none is reported here rather than one that is wrong in a knowable direction. On ten replicas of one system differing only by seed, errors computed from one run were about ten times smaller than the spread of the ten means. The remedy is a longer run, or better, replicas.

![rg](../analysis/rg/rg.png)

### Hbonds

Hydrogen bond analysis.

**Parameters:**

- `selection`: `protein`
- `method`: `baker_hubbard`
- `criterion`: `Baker-Hubbard: H...A < 2.5 A, D-H...A > 120 deg`
- `periodic`: `True`
- `freq`: `0.1`
- `candidate_freq`: `0.0`
- `count_multiplier`: `1`
- `n_persistent_bonds`: `0`
- `n_candidate_bonds`: `0`

**Findings:**

- Mean over the equilibrated part of the run: 0, from 1 independent sample after discarding 66 frames.
- This run is not long against its own correlation time: 34 frames hold 1.0 independent samples by their own estimate, and a correlated run shorter than 25 of them reads its correlation time low (only frames whose statistical inefficiency, before the correction for the sample mean, is under 2 are exempt). That count is an upper bound, so an error computed from it would be a lower bound -- and none is reported here rather than one that is wrong in a knowable direction. On ten replicas of one system differing only by seed, errors computed from one run were about ten times smaller than the spread of the ten means. The remedy is a longer run, or better, replicas.

![hbonds](../analysis/hbonds/hbonds.png)

### SS

Secondary structure assignment.

**Parameters:**

- `selection`: `protein`
- `simplified`: `True`

![ss](../analysis/ss/ss.png)

### SASA

Solvent-Accessible Surface Area (SASA).

**Parameters:**

- `selection`: `protein`
- `mode`: `total`
- `probe_radius`: `0.14`
- `n_sphere_points`: `960`
- `with_ligand`: `False`

**Findings:**

- Mean over the equilibrated part of the run: 4.7592 ± 0.0075 nm², from 24 independent samples after discarding 60 frames.

![sasa](../analysis/sasa/sasa.png)

### Dihedrals

Backbone dihedrals and the Ramachandran plot.

**Parameters:**

- `density`: `True`
- `angles`: `['phi', 'psi', 'omega']`
- `bins`: `72`

![dihedrals](../analysis/dihedrals/dihedrals.png)

### Order\_Parameters

Backbone N--H generalised order parameters, S^2.

**Parameters:**

- `selection`: `protein`
- `align_selection`: `protein and name CA`
- `ref`: `0`
- `definition`: `Lipari-Szabo generalised order parameter, from the closed form of the internal second-rank correlation plateau after superposition`

![order\_parameters](../analysis/order_parameters/order_parameters.png)

### Thermodynamics

Density, energy and temperature, from the record the run already kept.


![thermodynamics](../analysis/thermodynamics/thermodynamics.png)

### End\_To\_End

The distance between the two ends of a chain, per frame.

**Parameters:**

- `selection`: `protein`
- `atom`: `CA`
- `by_chain`: `True`

**Findings:**

- Mean over the equilibrated part of the run: 0.7004 nm, from fewer than one independent sample after discarding 66 frames.
- This run is not long against its own correlation time: 34 frames hold 0.7 independent samples by their own estimate, and a correlated run shorter than 25 of them reads its correlation time low (only frames whose statistical inefficiency, before the correction for the sample mean, is under 2 are exempt). That count is an upper bound, so an error computed from it would be a lower bound -- and none is reported here rather than one that is wrong in a knowable direction. On ten replicas of one system differing only by seed, errors computed from one run were about ten times smaller than the spread of the ten means. The remedy is a longer run, or better, replicas.

![end\_to\_end](../analysis/end_to_end/end_to_end.png)

### Moments\_Of\_Inertia

The three principal moments of inertia, per frame.

**Parameters:**

- `selection`: `protein`

![moments\_of\_inertia](../analysis/moments_of_inertia/moments_of_inertia.png)

### Cluster

Conformational clustering.

**Parameters:**

- `selection`: `protein and name CA`
- `methods`: `['kmeans', 'hierarchical']`
- `features`: `rmsd`
- `n_clusters`: `5`
- `eps`: `0.2`
- `min_samples`: `5`
- `linkage`: `average`
- `random_state`: `42`
- `n_init`: `10`

![cluster: hierarchical](../analysis/cluster/cluster_hierarchical.png)

![cluster: hierarchical counts](../analysis/cluster/cluster_hierarchical_counts.png)

![cluster: hierarchical dendrogram](../analysis/cluster/cluster_hierarchical_dendrogram.png)

![cluster: kmeans](../analysis/cluster/cluster_kmeans.png)

![cluster: kmeans counts](../analysis/cluster/cluster_kmeans_counts.png)

![cluster: rmsd matrix](../analysis/cluster/cluster_rmsd_matrix.png)

### Dimred

Dimensionality reduction.

**Parameters:**

- `selection`: `protein and name CA`
- `methods`: `['pca']`
- `n_components`: `2`
- `perplexity`: `30.0`
- `n_neighbors`: `15`
- `min_dist`: `0.1`
- `random_state`: `42`
- `landscape_bins`: `40`

![dimred: pca](../analysis/dimred/dimred_pca.png)

![dimred: pca landscape](../analysis/dimred/dimred_pca_landscape.png)


## Convergence

A frame is not an observation. Consecutive frames of a trajectory are nearly the same structure, so the number of independent observations is set by how quickly each observable forgets where it was, not by how often frames were written. The uncertainties below count the former.

| observable | frames | discarded | independent | mean | uncertainty | equilibrated |
|---|---|---|---|---|---|---|
| potential energy | 2 | 0 | 2.0 | -20,008 | not determined | too short to say |
| temperature | 2 | 0 | 2.0 | 289.3 | not determined | too short to say |
| density | 2 | 0 | 2.0 | 0.9977 | not determined | too short to say |
| RMSD | 100 | 66 | 16.1 | 0.0078 | 0.0010 | yes |
| Rg | 100 | 66 | 0.5 | 0.3355 | not determined | no |
| SASA | 100 | 60 | 24.2 | 4.7591 | 0.0075 | yes |

### The checks this run was held to

| check | result | on what |
|---|---|---|
| Each observable equilibrates before it is averaged | **failed** | still moving: Rg |
| Each observable's correlation time is resolved by the run | **failed** | not resolved: Rg |
| Each mean rests on at least 10 independent samples | **failed** | too few: density, potential energy, temperature |
| The mean temperature is within 5 K of the target | **failed** | 289.3 K against 300.0 K |
| The potential energy's trend stays under 1 kJ/mol per ns per atom | not judged | 2 energy records, too few to judge |

### What this run cannot support

- The mean temperature was 289.3 K against a target of 300.0 K. A thermostat that does not hold its target means the ensemble is not the one the settings describe.
- Still moving in one direction: Rg. The run has not finished equilibrating, so averages over it describe the approach rather than the state.
- Too short to say whether it has equilibrated: density, potential energy, temperature. Drift is judged by comparing the start of the run against the end, and a series this short has no start and end to compare.
- Too short to tell how correlated it is: Rg. The correlation time is estimated by summing until the series forgets itself, and a run this length cannot see a memory longer than half of it. Any independence reported here is the independence the estimate failed to rule out, not independence that was found.
- Too few independent samples to average: density, potential energy, temperature. Consecutive frames are nearly the same structure, so the number of independent observations is set by how fast the observable forgets, not by how often it was written out.

**What would support it.** Relative shape anisotropy, asphericity, radius of gyration, end-to-end distance, acylindricity, polar SASA and hydrophobic SASA withheld their means for want of sampling: at least 78 ps more production should give each 25 independent samples (a floor: the run is too short to resolve its own correlation time, so run that much and estimate again); at this run's own speed on CPU, about 3 min. This extends the study in place, joins the segments and reruns the analyses:

```yaml
simulation:
  resume_from: /tmp/D/real/clean
  extra_ns: 0.078
```


## Discussion

_Everything above is a record: what was built, what was run, what was computed from it, and how much of it the run supports. What none of it says is what the system was doing, or whether that answers the question the study was for. That judgement is yours, and this section is where it goes._

_Software can report that a mean rests on eight independent samples. It cannot tell you whether eight is enough for the claim you intend to make._

## Citation

If you use FastMDXplora in your work, please cite:

> Aina, A.; Kwan, D. FastMDAnalysis: Software for Automated Analysis of Molecular Dynamics Trajectories. J. Comput. Chem. 2026, 47, e70350. DOI: 10.1002/jcc.70350

BibTeX:

```bibtex
@article{aina2026fastmd,
  author  = {Aina, Adekunle and Kwan, Derrick},
  title   = {FastMDAnalysis: Software for Automated Analysis of Molecular
             Dynamics Trajectories},
  journal = {Journal of Computational Chemistry},
  volume  = {47},
  number  = {8},
  pages   = {e70350},
  year    = {2026},
  doi     = {10.1002/jcc.70350},
}
```

## Reproducibility

- **FastMDXplora version**: `0.1.0`
- **Source commit**: `be6ea867cb0b`
- **Python**: `3.11.17`
- **Platform**: `Linux-6.18.44-fc-v139-x86_64-with-glibc2.39`
- **System input**: `/tmp/D/real/tri-ala.pdb`
- **Output directory**: `/tmp/D/real/clean`

Per-phase parameter manifests for phases in this workflow are preserved at `setup/setup_parameters.json`, `simulation/simulation_parameters.json`, `analysis/analysis_manifest.json`. The complete session manifest is at `manifest.json` at the project root.

**What rerunning the configuration repeats.** Running `resolved_config.yml` again prepares the same system: it records the random seed setup placed hydrogens and ions with (859956268, drawn for this study), and setup minimises on one CPU thread, so the same structure and settings give the same atoms in the same places. No `simulation.random_seed` was given, so a rerun starts its dynamics from velocities drawn afresh: a new trajectory of the same system. To simulate this very system without preparing it again, point `simulation.setup_from` at its `setup/` directory.
