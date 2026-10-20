# Cross-run comparison report

This study compared **3 successful runs** over `random_seed` (3 values).

## RMSD

![RMSD overlay](overlay_rmsd.png)

## Radius of gyration

![Radius of gyration overlay](overlay_rg.png)

## Total SASA

## Summary

| Run | mean RMSD (nm) | mean Radius of gyration (nm) | mean Total SASA (nm²) |
|---|---|---|---|
| s1__random-seed-1 | 0.0070 ± 0.0010 | 0.31396 ± 0.00096 | 4.6651 ± 0.0065 |
| s1__random-seed-2 | 0.0080 ± 0.0011 |  |  |
| s1__random-seed-3 | 0.0090 ± 0.0012 | 0.3355 \* | 4.7592 ± 0.0075 |

Each mean is over the frames after equilibration, with its standard error, as its run's analysis recorded it.

\* Not determined: the run's analysis found the series too short, or with too few independent samples, for its mean to be determined.

Full table: `comparison_summary.csv`.

## Across the replicas

These runs are replicas: the members differ only by random seed, so they are repeats of one estimate and the spread of their means is its error. That spread is set here against the error each run estimated for itself.

| Analysis | mean of the means | spread of the means | error each run estimated |
|---|---|---|---|
| end_to_end | 0.6827 | 0.025 | 0.0024 |
| hbonds | 0.002778 | 0.0039 | 0.0056 |
| rg | 0.3248 | 0.015 | 0.00096 |
| rmsd | 0.008 | 0.001 | 0.0011 |
| sasa | 4.712 | 0.067 | 0.007 |

**end_to_end.** The replicas spread 0.0251 where a single run predicted 0.00238 (from the 1 of 2 runs that gave an error of their own), more than the factor of 2 that ten replicas can resolve. The single-run estimate reads the correlation time off a trajectory too short to contain it, so the error bars printed for this observable are too tight and the replicas are the honest figure. The runs said as much themselves.

**hbonds.** The replicas spread 0.00393 against a predicted 0.00559 (from the 1 of 2 runs that gave an error of their own), within the factor of 2 that this many replicas can resolve: the error this observable reports for itself is supported by repeating it.

**rg.** The replicas spread 0.0153 where a single run predicted 0.000961 (from the 1 of 2 runs that gave an error of their own), more than the factor of 2 that ten replicas can resolve. The single-run estimate reads the correlation time off a trajectory too short to contain it, so the error bars printed for this observable are too tight and the replicas are the honest figure. The runs said as much themselves.

**rmsd.** The replicas spread 0.001 against a predicted 0.0011, within the factor of 2 that this many replicas can resolve: the error this observable reports for itself is supported by repeating it.

**sasa.** The replicas spread 0.0665 where a single run predicted 0.007, more than the factor of 2 that ten replicas can resolve. The single-run estimate reads the correlation time off a trajectory too short to contain it, so the error bars printed for this observable are too tight and the replicas are the honest figure.
