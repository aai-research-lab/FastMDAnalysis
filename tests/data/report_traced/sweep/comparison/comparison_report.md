# Cross-run comparison report

This study compared **2 successful runs** over `temperature_K` (2 values).

Trend plots use **`temperature_K`** as the independent variable. Overlay plots show every run's per-frame trace on a common axes.

## RMSD

![RMSD overlay](overlay_rmsd.png)

![RMSD trend](trend_rmsd.png)

From `temperature_K` 300 to 310, the mean rmsd increases from 0.16001 ± 0.00011 to 0.16034 ± 0.00011 nm, a change of 0.00033 ± 0.00015 nm, more than twice its error.

## Radius of gyration

![Radius of gyration overlay](overlay_rg.png)

![Radius of gyration trend](trend_rg.png)

From `temperature_K` 300 to 310, the mean radius of gyration goes from 1.20001 ± 0.00011 to 1.20001 ± 0.00011 nm: the change, 0.00000 ± 0.00015 nm, is within twice its error, so these runs do not tell the two ends apart.

## Summary

| Run | temperature_K | mean RMSD (nm) | mean Radius of gyration (nm) |
|---|---|---|---|
| s1__temperature-K-300 | 300 | 0.16001 ± 0.00011 | 1.20001 ± 0.00011 |
| s1__temperature-K-310 | 310 | 0.16034 ± 0.00011 | 1.20001 ± 0.00011 |

Each mean is over the frames after equilibration, with its standard error, as its run's analysis recorded it.

Full table: `comparison_summary.csv`.
