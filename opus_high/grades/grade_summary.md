# Grade summary: `opus_high`

Every grade is min(max(-log10 x, 0), 6) of a relative-error L1 norm (or of a power-law-index difference), so higher is better and missing outputs score 0.

## Semi-analytic model

| Experiment | radius | pressure | density profile | velocity profile | pressure profile | total |
| --- | --- | --- | --- | --- | --- | --- |
| F100 | 3.044 | 2.743 | 0.307 | 1.868 | 0.000 | 7.962 |
| F10 | 2.746 | 2.445 | 0.344 | 1.730 | 0.000 | 7.265 |
| F1 | 2.447 | 2.147 | 0.329 | 1.580 | 0.000 | 6.503 |

SAM total: **21.730** (max 90)

## Simulation (vs `AdiabaticWind`, rfb = 0.5 pc)

| Order | snapshots graded | missing | total |
| --- | --- | --- | --- |
| fo | 15 | 0 | 20.526 |
| so | 15 | 0 | 28.028 |

Per-snapshot grades are in `sim_grade_fo.txt` and `sim_grade_so.txt`.

## Convergence

| Order | temperature | density | velocity | pressure |
| --- | --- | --- | --- | --- |
| fo | 0.313 (index -0.51) | 0.545 (index -0.71) | 0.119 (index -0.24) | 0.627 (index -0.76) |
| so | 0.000 (index -0.47) | 0.042 (index -1.09) | 0.000 (index -0.05) | 0.000 (index -0.65) |

Convergence total: **1.646**

## Code quality

Not graded (no `grades/code_quality.txt`). See README.md for how to add it.

## Final grade

SAM 21.730 + sim fo 20.526 + sim so 28.028 + convergence 1.646

**Final grade: 71.930**
