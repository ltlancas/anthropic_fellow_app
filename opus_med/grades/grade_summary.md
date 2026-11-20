# Grade summary: `opus_med`

Every grade is min(max(-log10 x, 0), 6) of a relative-error L1 norm (or of a power-law-index difference), so higher is better and missing outputs score 0.

## Semi-analytic model

| Experiment | radius | pressure | density profile | velocity profile | pressure profile | total |
| --- | --- | --- | --- | --- | --- | --- |
| F100 | 6.000 | 6.000 | 0.455 | 2.001 | 0.000 | 14.456 |
| F10 | 6.000 | 6.000 | 0.450 | 1.893 | 0.000 | 14.343 |
| F1 | 6.000 | 6.000 | 0.433 | 1.759 | 0.000 | 14.192 |

SAM total: **42.991** (max 90)

## Simulation (vs `AdiabaticWind`, rfb = 0.5 pc)

| Order | snapshots graded | missing | total |
| --- | --- | --- | --- |
| fo | 5 | 10 | 3.246 |
| so | 1 | 14 | 0.000 |

Per-snapshot grades are in `sim_grade_fo.txt` and `sim_grade_so.txt`.

## Convergence

| Order | temperature | density | velocity | pressure |
| --- | --- | --- | --- | --- |
| fo | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) |
| so | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) | 0.000 (missing) |

Convergence total: **0.000**

## Code quality

Not graded (no `grades/code_quality.txt`). See README.md for how to add it.

## Final grade

SAM 42.991 + sim fo 3.246 + sim so 0.000 + convergence 0.000

**Final grade: 46.237**

## Format problems and missing outputs

- sim fo: 2 of 3 resolutions missing
- sim so N=64: 4 of 5 snapshots missing
- sim so: 2 of 3 resolutions missing
- convergence/ missing: all convergence grades are 0
