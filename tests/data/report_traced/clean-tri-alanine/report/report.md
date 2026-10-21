# FastMDXplora Study: tri\-ala

_Generated: 2026-10-18 (UTC)_
_Software: FastMDXplora_
_Dashboard: [dashboard.html](dashboard.html)_

## Summary

The prepared system contained 1,412 atoms and was simulated for 40 ps at 300 K. All 4 observables assessed had equilibrated, and each mean is determined.

## Methods

**System preparation.** Starting coordinates were taken from `/tmp/D/real/tri-ala.pdb`. The file itself is recorded in the setup manifest by its SHA-256 digest (e1c61caf90ab), which identifies it independently of the path it was read from. Missing heavy atoms and hydrogens were added with PDBFixer at pH 7.0, which sets the protonation of ionizable side chains. Heterogens: non-standard residues were decided per component: any ligand that could be parameterised was kept and prepared, crystallographic water and buffer were removed, and setup stopped rather than guess where the structure did not determine what to simulate. The protein was described by the amber14-all.xml, amber14/tip3p.xml force field with the tip3p water model. The complex was solvated in a dodecahedron box sized to leave at least 1.40 nm between the solute and its nearest periodic image (grown from the 1.20 nm asked for, so that the box is at least twice the cutoff across), with Na+/Cl- ions at 0.15 M. The solvated system contained 1,412 atoms. Hydrogens and ions were placed with random seed 1651082175 (`setup.random_seed`), which reproduces the solvated system. Long-range electrostatics were treated with PME and a real-space cutoff of 0.90 nm. Bonds were constrained with `HBonds`, and water was held rigid.

**Simulation protocol.** The system was energy-minimized before dynamics. Equilibration comprised 0.4 ps in the NVT ensemble followed by 0.4 ps in the NPT ensemble. Production dynamics were run for 40 ps in the NPT ensemble. Equations of motion were integrated with the Langevin middle-scheme integrator and a 2 fs timestep, at 300.0 K with a friction coefficient of 1.0 ps⁻¹. Pressure was maintained at 1.0 bar, with the barostat applied every 25 steps. Coordinates were written every 0.2 ps. Simulations used OpenMM on the CPU platform. No random seed was fixed, so velocities differ between runs.

**Analysis.** The per-frame quantity `rmsd` was averaged over the frames after an equilibration period detected automatically, as the latest start that keeps at least 90% of the most statistically independent samples any start keeps (after Chodera, J. Chem. Theory Comput. 2016, 12, 1799). Standard errors account for the correlation between frames through the statistical inefficiency, the number of frames per independent sample, estimated from the autocorrelation function summed over Geyer's initial positive sequence (Stat. Sci. 1992, 7, 473) and corrected for the bias of an autocorrelation taken about the sample mean. Where frames were discarded, the error is that of the whole run scaled to the frames kept, unless the discard at least doubled the independent samples. An error is given only where the frames averaged hold at least 10 independent samples and span at least 25 times their own statistical inefficiency. Frames whose statistical inefficiency, before the correction for the sample mean, is under 2 are exempt from the second condition; their error is still computed with the corrected inefficiency. The first 10 of the 200 frames analysed were discarded as equilibration.

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
- **random\_seed**: `1651082175`, drawn for this study
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
- **duration\_ns**: `0.04`
- **nvt\_duration\_ns**: `0.0004`
- **npt\_duration\_ns**: `0.0004`
- **nvt\_steps**: `200`
- **ensemble**: `npt`
- **npt\_steps**: `200`
- **production\_steps**: `20000`
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
- **trajectory\_interval\_steps**: `100`
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

Analysis was performed on a trajectory of 200 frames: 3 protein residues and 2 other residues (ions and any cofactors kept by the save selection).

Analyses performed: rmsd, thermodynamics.

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

- Mean over the equilibrated part of the run: 0.01037 ± 0.00085 nm, from 82 independent samples after discarding 10 frames.

![rmsd](../analysis/rmsd/rmsd.png)

### Thermodynamics

Density, energy and temperature, from the record the run already kept.


![thermodynamics](../analysis/thermodynamics/thermodynamics.png)


## Convergence

A frame is not an observation. Consecutive frames of a trajectory are nearly the same structure, so the number of independent observations is set by how quickly each observable forgets where it was, not by how often frames were written. The uncertainties below count the former.

| observable | frames | discarded | independent | mean | uncertainty | equilibrated |
|---|---|---|---|---|---|---|
| potential energy | 20 | 1 | 12.8 | -19,798 | 43 | yes |
| temperature | 20 | 2 | 13.4 | 298.0 | 2.4 | yes |
| density | 20 | 2 | 16.0 | 0.9989 | 0.0025 | yes |
| RMSD | 200 | 10 | 82.3 | 0.01037 | 0.00085 | yes |

### The checks this run was held to

| check | result | on what |
|---|---|---|
| Each observable equilibrates before it is averaged | passed | on density, potential energy, RMSD, temperature |
| Each observable's correlation time is resolved by the run | passed | on density, potential energy, RMSD, temperature |
| Each mean rests on at least 10 independent samples | passed | fewest: potential energy, 13 |
| The mean temperature is within 5 K of the target | passed | 298.0 K against 300.0 K |
| The potential energy's trend stays under 1 kJ/mol per ns per atom | passed | 0.59 kJ/mol per ns per atom |

**What would support it.** Total energy withheld its mean for want of sampling: 2 ps more production should give it 10 independent samples; at this run's own speed on CPU, about 1 min. This extends the study in place, joins the segments and reruns the analyses:

```yaml
simulation:
  resume_from: /tmp/D/ws/clean-tri-alanine
  extra_ns: 0.002
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
- **Output directory**: `/tmp/D/ws/clean-tri-alanine`

Per-phase parameter manifests for phases in this workflow are preserved at `setup/setup_parameters.json`, `simulation/simulation_parameters.json`, `analysis/analysis_manifest.json`. The complete session manifest is at `manifest.json` at the project root.

**What rerunning the configuration repeats.** Running `resolved_config.yml` again prepares the same system: it records the random seed setup placed hydrogens and ions with (1651082175, drawn for this study), and setup minimises on one CPU thread, so the same structure and settings give the same atoms in the same places. No `simulation.random_seed` was given, so a rerun starts its dynamics from velocities drawn afresh: a new trajectory of the same system. To simulate this very system without preparing it again, point `simulation.setup_from` at its `setup/` directory.
